"""Portable fixed-orbit comparison and checkpointed Scheme B experiments.

All outputs go to --output; existing research caches are never output targets.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import binary_single as bs
import binary_single_montecarlo as mc
import haloconcentration as hc
import p_a_j as prior


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2,
                                   allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def source_hashes():
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(ROOT.glob('*.py'))}


def provenance():
    return dict(python=sys.version, platform=platform.platform(),
                packages={n: importlib.metadata.version(n) for n in
                          ('numpy', 'scipy', 'matplotlib', 'astropy', 'hmf', 'camb')},
                sources=source_hashes())


def fixed_orbit_run(output, reference=None):
    """Real simulation outputs for the same saved inputs on both computers."""
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    config = bs.BinarySingleConfig()
    distribution = prior.AppendixDOrbitalDistribution()
    sample = distribution.sample(128, np.random.default_rng(20260925))
    regenerated_a, regenerated_e = sample.semi_major_axis_pc, sample.eccentricity
    if reference is None:
        a, e = regenerated_a, regenerated_e
    else:
        with np.load(reference / 'fixed_orbits.npz', allow_pickle=False) as saved:
            a, e = saved['initial_a'].copy(), saved['initial_e'].copy()
    arrays = dict(initial_a=a, initial_e=e, sampler_a=regenerated_a,
                  sampler_e=regenerated_e)
    times, redshifts = mc.global_time_grid(12, 200, hc.cosmology_for_model(config.halo_model))
    arrays.update(time_edges=times, redshift_edges=redshifts)
    rows = []
    for z in (12, 0):
        shell = bs.build_shell_environment(1.15e12, z, config)
        for field in ('shell_mass_msun', 'environment_density_msun_pc3',
                      'velocity_dispersion_km_s', 'hard_semimajor_axis_pc'):
            arrays[f'z{z}_{field}'] = getattr(shell, field)
        for i in (0, 4):
            for mode in ('full', 'gw_only'):
                environment = (shell if mode == 'full' else replace(shell,
                    environment_density_msun_pc3=np.zeros_like(shell.shell_mass_msun)))
                event = bs.evolve_binary_batch(a, e, environment, i, config)
                prefix = f'z{z}_R{i+1}_{mode}_'
                for field in ('hard_mask', 'merged_mask', 'invalid_mask',
                              'k_above_calibration_mask', 'merger_time_myr'):
                    arrays[prefix+field] = getattr(event, field)
                arrays[prefix+'final_a'] = event.final_state.semi_major_axis_pc
                arrays[prefix+'final_e'] = event.final_state.eccentricity
                rows.append(dict(z=z, shell=i+1, mode=mode, sampled=event.sample_count,
                    hard=event.hard_count, merged=event.merger_count, invalid=event.invalid_count))
    np.savez_compressed(output / 'fixed_orbits.npz', **arrays)
    metadata = dict(provenance=provenance(), config=asdict(config),
                    distribution=asdict(distribution), seed=20260925,
                    elapsed_seconds=time.perf_counter()-started, rows=rows)
    if sys.platform != 'win32':
        import resource
        metadata['peak_rss_mib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    write_json(output / 'fixed_orbits.json', metadata)
    if reference is not None:
        baseline = json.loads((reference / 'fixed_orbits.json').read_text(encoding='utf-8'))
        comparisons = []
        with np.load(reference / 'fixed_orbits.npz', allow_pickle=False) as expected:
            for name, actual in arrays.items():
                target = expected[name]
                if actual.dtype.kind in 'biu' or name.startswith('initial_'):
                    ok = np.array_equal(actual, target)
                    comparisons.append(dict(array=name, passed=bool(ok), tolerance='exact'))
                else:
                    atol = 1e-14 if name.endswith(('_a', '_pc')) else 1e-8
                    ok = np.allclose(actual, target, rtol=2e-5, atol=atol, equal_nan=True)
                    finite = np.isfinite(actual) & np.isfinite(target)
                    error = float(np.max(np.abs(actual[finite]-target[finite]))) if finite.any() else 0.0
                    comparisons.append(dict(array=name, passed=bool(ok),
                        rtol=2e-5, atol=atol, max_absolute_error=error))
        same_sources = baseline['provenance']['sources'] == metadata['provenance']['sources']
        report = dict(passed=bool(same_sources and all(r['passed'] for r in comparisons)),
                      same_sources=same_sources, comparisons=comparisons,
                      baseline_platform=baseline['provenance']['platform'],
                      runtime_platform=metadata['provenance']['platform'])
        write_json(output / 'comparison.json', report)
        if not report['passed']:
            raise RuntimeError('本地/云端对照存在差异；请查看 comparison.json，暂不扩大样本。')
    print(json.dumps(metadata, ensure_ascii=False, indent=2), flush=True)


def cohort_job(output, seed, mode, samples, accretion, timestep, entry):
    config = bs.BinarySingleConfig(global_timestep_myr=timestep,
                                   ksi=0.0 if mode == 'gw_only' else 1.0)
    population = mc.CohortPopulationConfig(initial_samples_per_shell=samples,
                    accretion_samples_per_shell=accretion, entry_orbit_model=entry)
    identity = dict(mass_msun=1.5e6, config=asdict(config), population=asdict(population),
                    seed=seed, sources=source_hashes(), distribution=asdict(prior.AppendixDOrbitalDistribution()))
    job_hash = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:12]
    job_dir = output / f'B_{entry}_{mode}_seed{seed}_{job_hash}'
    job_dir.mkdir(parents=True, exist_ok=True)
    if (job_dir / 'summary.json').exists():
        print('Completed job retained:', job_dir, flush=True)
        return json.loads((job_dir / 'summary.json').read_text(encoding='utf-8'))
    write_json(job_dir / 'config.json', identity)
    started = time.perf_counter()

    def progress(completed, total):
        write_json(job_dir / 'progress.json', dict(completed=completed, total=total,
                    session_elapsed_seconds=time.perf_counter()-started, resumable=True))
        print(f'{job_dir.name}: {completed}/{total}', flush=True)

    history = mc.run_cohort_shell_monte_carlo(1.5e6, config=config,
        population_config=population, random_seed=seed, progress_callback=progress,
        checkpoint_path=job_dir / 'checkpoint.npz', resume=True)
    arrays = {name: value for name, value in vars(history).items() if isinstance(value, np.ndarray)}
    for i, pop in enumerate(history.final_populations):
        for name, value in vars(pop).items():
            arrays[f'R{i+1}_{name}'] = value
    with (job_dir / 'history.npz.partial').open('wb') as stream:
        np.savez_compressed(stream, **arrays)
    (job_dir / 'history.npz.partial').replace(job_dir / 'history.npz')
    zmid = mc.redshift_from_cosmic_age_myr(
        0.5*(history.time_edges_myr[:-1]+history.time_edges_myr[1:]), hc.cosmology_for_model(config.halo_model))
    dt = np.diff(history.time_edges_myr)
    rate = history.rate_shell_per_year.sum(axis=1)
    np.savetxt(job_dir / 'rates.csv', np.column_stack((zmid, dt, rate,
               history.merged_weight_per_step.sum(axis=1))), delimiter=',',
               header='z_mid,duration_myr,rate_per_year,merged_weight', comments='')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(zmid, rate, '.-', lw=0.8)
    ax.set(xlabel='z (time-bin midpoint)', ylabel='Raw bin rate [yr^-1 halo^-1]',
           title=f'Scheme B {entry}, {mode}; seed={seed}')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(job_dir / 'rates.png', dpi=160)
    plt.close(fig)
    summary = dict(job=job_dir.name, seed=seed, mode=mode, entry=entry,
                   samples_drawn=int(history.drawn_per_boundary.sum()),
                   inside_merged_weight=float(history.merged_weight_per_step.sum()),
                   outside_merged_weight=float(history.outside_merged_weight_per_boundary.sum()),
                   final_alive_weight=float(history.surviving_weight_at_boundary[-1].sum()),
                   balance_max_abs=float(np.max(np.abs(history.population_balance_residual))),
                   session_elapsed_seconds=time.perf_counter()-started,
                   provenance=provenance(), note='Small-sample execution pilot; no quantitative convergence claim.')
    write_json(job_dir / 'progress.json', dict(completed=len(dt), total=len(dt), status='complete'))
    write_json(job_dir / 'summary.json', summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', choices=['reference', 'compare', 'cohort'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path, default=ROOT / 'cloud' / 'reference')
    parser.add_argument('--seed', type=int, default=20260923)
    parser.add_argument('--samples', type=int, default=128)
    parser.add_argument('--accretion', type=int, default=8)
    parser.add_argument('--timestep', type=float, default=200)
    parser.add_argument('--mode', choices=['full', 'gw_only'], default='full')
    parser.add_argument('--entry', choices=['gw_aged', 'pristine'], default='gw_aged')
    args = parser.parse_args()
    if args.task in ('reference', 'compare'):
        fixed_orbit_run(args.output, args.reference if args.task == 'compare' else None)
    else:
        cohort_job(args.output, args.seed, args.mode, args.samples, args.accretion, args.timestep, args.entry)
