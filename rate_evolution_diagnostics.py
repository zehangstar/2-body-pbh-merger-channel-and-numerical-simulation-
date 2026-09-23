"""Reproducible rate-trend diagnostics; no fit to the target redshift slope."""

import argparse
import csv
from dataclasses import asdict, replace
import json
import hashlib
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import binary_single as bs
import binary_single_montecarlo as mc
import haloconcentration as hc
import p_a_j as prior

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'pic_and_data/rate_evolution_diagnosis'
BASE = ROOT / 'pic_and_data/r_bs_perhalo'
META = json.loads((BASE / 'run_config.json').read_text(encoding='utf-8'))
CFG = bs.BinarySingleConfig(**META['config'])
M0 = META['halo_M0_msun']


def save_rows(name, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_json(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')


def paired_history():
    """Replay exactly the baseline initial conditions; preserve event pairing."""
    data = dict(np.load(BASE / 'raw_shell_history.npz'))
    distribution = prior.AppendixDOrbitalDistribution(**META['prior'])
    rng = np.random.default_rng(META['seed'])
    rows = []
    start = time.perf_counter()
    for k, z in enumerate(data['redshift'][:-1]):
        duration = float(data['time'][k+1] - data['time'][k])
        cfg = replace(CFG, global_timestep_myr=duration)
        shell = bs.build_shell_environment(M0, float(z), cfg)
        noenv = replace(shell, environment_density_msun_pc3=np.zeros(len(shell.shell_mass_msun)))
        for i, mass in enumerate(shell.shell_mass_msun):
            # Current baseline uses one complete chunk per shell.
            sample = distribution.sample(CFG.sample_count_per_shell_and_step, rng)
            full = bs.evolve_binary_batch(sample.semi_major_axis_pc, sample.eccentricity, shell, i, cfg)
            gw = bs.evolve_binary_batch(sample.semi_major_axis_pc, sample.eccentricity, noenv, i, cfg)
            plus = int(np.count_nonzero(full.merged_mask & ~gw.merged_mask))
            minus = int(np.count_nonzero(gw.merged_mask & ~full.merged_mask))
            n = full.sample_count
            population = cfg.pbh_fraction_total * cfg.fraction_binary * mass / (cfg.primary_mass_msun + cfg.secondary_mass_msun)
            factor = population / n / (duration * 1e6)
            # D = I(full) - I(GW), with matched initial conditions.
            paired_count_variance = (plus + minus - (plus-minus)**2/n) * n/(n-1)
            rows.append(dict(step=k, shell=i+1, z_start=float(z), duration_myr=duration,
                mass=mass, halo_mass=shell.halo_mass_msun, rho=shell.environment_density_msun_pc3[i],
                velocity=shell.velocity_dispersion_km_s[i], ah=shell.hard_semimajor_axis_pc[i],
                sampled=n, hard=full.hard_count, full=full.merger_count, gw=gw.merger_count,
                plus=plus, minus=minus, invalid=full.invalid_count+gw.invalid_count,
                baseline_difference=full.merger_count-int(data['merged'][k,i]),
                full_rate=full.merger_count*factor, gw_rate=gw.merger_count*factor,
                excess_rate=(plus-minus)*factor, excess_se=np.sqrt(paired_count_variance)*factor))
        if (k+1)%10 == 0 or k+1 == len(data['merged']):
            print(f'paired {k+1}/{len(data["merged"])}: {time.perf_counter()-start:.1f}s', flush=True)
    save_rows('paired_shells.csv', rows)
    print('baseline mismatches',sum(r['baseline_difference'] != 0 for r in rows),'invalid',sum(r['invalid'] for r in rows),flush=True)


class Figure19Prior:
    def sample(self, n, rng):
        return prior.sample_figure19_orbital_parameters(n, rng)


def variants():
    rows = []
    for name, model, distribution in [
        ('paper_order_prada', CFG.halo_model, Figure19Prior()),
        ('strict_ludlow', 'ludlow16', prior.AppendixDOrbitalDistribution(**META['prior'])),
    ]:
        cfg = replace(CFG, halo_model=model)
        result = mc.run_raw_shell_monte_carlo(M0, cfg, distribution, META['seed'], META['chunk_size'],
            lambda k,n: print(f'{name}: {k}/{n}',flush=True) if k%20 == 0 or k == n else None)
        np.savez_compressed(OUT / f'{name}.npz', time=result.elapsed_time_edges_myr,
            redshift=result.redshift_edges, merged=result.mergers_per_step,
            sampled=result.sampled_per_step, invalid=result.invalid_per_step)
        for k,z in enumerate(result.redshift_edges[:-1]):
            shell = bs.build_shell_environment(M0,float(z),cfg)
            dt = result.time_edges_myr[k+1]-result.time_edges_myr[k]
            pop = cfg.pbh_fraction_total*cfg.fraction_binary*shell.shell_mass_msun/(cfg.primary_mass_msun+cfg.secondary_mass_msun)
            rate = np.sum(pop*result.mergers_per_step[k]/result.sampled_per_step[k])/(dt*1e6)
            zmid = float(mc.redshift_from_cosmic_age_myr((result.time_edges_myr[k+1]+result.time_edges_myr[k])/2,hc.cosmology_for_model(model)))
            rows.append(dict(variant=name,step=k,z_mid=zmid,rate=rate,invalid=int(result.invalid_per_step[k].sum())))
    save_rows('variants.csv',rows)


def duration_controls():
    """Freeze the same z=0 halo and initial sample; vary only the observation window."""
    shell = bs.build_shell_environment(M0,0.0,CFG)
    noenv = replace(shell,environment_density_msun_pc3=np.zeros(len(shell.shell_mass_msun)))
    distribution = prior.AppendixDOrbitalDistribution(**META['prior'])
    sample = distribution.sample(20000,np.random.default_rng(2319))
    rows=[]
    for i in [0,4]:
        for dt in [84.03621897459016,100.0,200.0,400.0]:
            cfg=replace(CFG,global_timestep_myr=dt)
            a,e=sample.semi_major_axis_pc,sample.eccentricity
            full=bs.evolve_binary_batch(a,e,shell,i,cfg)
            gw=bs.evolve_binary_batch(a,e,noenv,i,cfg)
            rows.append(dict(shell=i+1,duration_myr=dt,sampled=len(a),full=full.merger_count,
                gw=gw.merger_count,full_probability_per_myr=full.merger_count/len(a)/dt,
                gw_probability_per_myr=gw.merger_count/len(a)/dt,
                invalid=full.invalid_count+gw.invalid_count))
    save_rows('duration_controls.csv',rows)


def paper_reference():
    import pymupdf
    source=ROOT / 'tmp/arxiv_2408.06515v2_source/plot_1526417.97_X.pdf'
    paths=pymupdf.open(source)[0].get_drawings()
    xticks=[paths[i]['items'][0][1].x for i in [2,4,6,8,10]]
    yticks=[paths[i]['items'][0][1].y for i in [44,46]]
    xcal=np.polyfit(xticks,[-3.,-2.,-1.,0.,1.],1)
    ycal=np.polyfit(yticks,[-10.,-9.],1)
    rows=[]
    for name,index in [('simulation',64),('polynomial',65)]:
        path=paths[index]
        points=[path['items'][0][1]]+[item[2] for item in path['items']]
        for point in points:
            rows.append(dict(curve=name,z=10**np.polyval(xcal,point.x),rate=10**np.polyval(ycal,point.y)))
    save_rows('paper_fig14_vector.csv',rows)
    smooth=[r for r in rows if r['curve']=='polynomial'][::-1]
    anchors=[]
    for z in [.003,.01,.1,1.,5.,8.]:
        rate=10**np.interp(np.log10(z),np.log10([r['z'] for r in smooth]),np.log10([r['rate'] for r in smooth]))
        anchors.append(dict(z=z,rate=rate))
    save_rows('paper_fig14_anchors.csv',anchors)
    save_json('paper_reference_metadata.json',dict(source=str(source.relative_to(ROOT)),
        sha256=hashlib.sha256(source.read_bytes()).hexdigest(),x_log10_calibration=xcal.tolist(),
        y_log10_calibration=ycal.tolist(),simulation_path=64,polynomial_path=65,
        mass_filename_msun=1526417.97,caption_mass_msun=1.5e6))


def numerical_check():
    from scipy.integrate import solve_ivp, quad
    rows=[]
    for z in [0.,12.]:
        shell=bs.build_shell_environment(M0,z,CFG)
        sample=prior.AppendixDOrbitalDistribution(**META['prior']).sample(160,np.random.default_rng(123))
        for i in [0,4]:
            batch=bs.evolve_binary_batch(sample.semi_major_axis_pc,sample.eccentricity,shell,i,CFG)
            # Reference one well-spaced subset; same ODE, independent DOP853 algorithm.
            for j in np.flatnonzero(batch.hard_mask)[::8]:
                a=sample.semi_major_axis_pc[j]; e=sample.eccentricity[j]
                cutoff=CFG.merger_radius_gm_c2*bs.G_SI*60*bs.SOLAR_MASS_KG/bs.C_SI**2/bs.PARSEC_M
                rhs=lambda u,y: bs.log_a_orbital_rhs(u,y,a,shell.environment_density_msun_pc3[i],shell.velocity_dispersion_km_s[i],shell.hard_semimajor_axis_pc[i],CFG)
                def time_event(u,y):
                    return y[0]-CFG.global_timestep_myr
                time_event.terminal=True
                time_event.direction=1
                sol=solve_ivp(rhs,(0,np.log(a/cutoff)),[0.,np.log((1-e)*(1+e))],method='DOP853',rtol=1e-10,atol=[1e-10,1e-12],events=time_event,max_step=.1)
                merged=sol.success and len(sol.t_events[0])==0
                rows.append(dict(z=z,shell=i+1,sample=int(j),success=bool(sol.success),
                    reference_merged=merged,batch_merged=bool(batch.merged_mask[j]),
                    time_error_myr=float(abs(sol.y[0,-1]-batch.merger_time_myr[j])) if merged and batch.merged_mask[j] else 0.))
    save_rows('dop853.csv',rows)
    cos=hc.cosmology_for_model(CFG.halo_model)
    h0=100*cos['h']*1000/(1e6*bs.PARSEC_M)*bs.MYR_S
    numerical_time=quad(lambda z: 1/((1+z)*h0*np.sqrt(cos['omega_m0']*(1+z)**3+cos['omega_lambda0'])),0,12,epsabs=1e-8)[0]
    analytic_time=float(mc.cosmic_age_myr(0,cos)-mc.cosmic_age_myr(12,cos))
    save_json('numerical_check.json',dict(orbit_count=len(rows),invalid=sum(not r['success'] for r in rows),
        classification_differences=sum(r['reference_merged']!=r['batch_merged'] for r in rows),
        max_time_error_myr=max(r['time_error_myr'] for r in rows),
        numerical_elapsed_myr=numerical_time,analytic_elapsed_myr=analytic_time,
        elapsed_difference_myr=numerical_time-analytic_time))


def summarize():
    paired=np.genfromtxt(OUT/'paired_shells.csv',delimiter=',',names=True)
    baseline=np.genfromtxt(BASE/'rate_per_bin.csv',delimiter=',',names=True)
    raw=dict(np.load(BASE/'raw_shell_history.npz'))
    steps=len(baseline)
    nshell=raw['merged'].shape[1]
    observed=paired['full'].reshape(steps,nshell)
    if not np.array_equal(observed,raw['merged']):
        raise ValueError('Paired replay does not match the original baseline.')
    rows=[]
    for k in range(steps):
        p=paired[paired['step']==k]
        full=p['full_rate'].sum(); gw=p['gw_rate'].sum(); excess=p['excess_rate'].sum()
        nbin=CFG.pbh_fraction_total*CFG.fraction_binary*baseline['halo_mass_msun'][k]/60
        rows.append(dict(step=k,z_mid=baseline['z_mid'][k],duration_myr=p['duration_myr'][0],
            halo_mass=baseline['halo_mass_msun'][k],binary_population=nbin,
            full_rate=full,gw_rate=gw,excess_rate=excess,excess_se=np.sqrt(np.sum(p['excess_se']**2)),
            gw_fraction=gw/full,full_rate_per_binary=full/nbin,
            weighted_p=full/nbin*p['duration_myr'][0]*1e6,
            positive_events=int(p['plus'].sum()),negative_events=int(p['minus'].sum())))
    save_rows('paired_rates.csv',rows)
    summary=dict(baseline_sha256=hashlib.sha256((BASE/'raw_shell_history.npz').read_bytes()).hexdigest(),
        baseline_metadata=META,replayed_sample_count=int(paired['sampled'].sum()),
        replay_mismatches=int(np.count_nonzero(observed-raw['merged'])),
        invalid_count=int(paired['invalid'].sum()),first_bin=rows[0],last_full_bin=rows[-2],short_last_bin=rows[-1],
        gw_weighted_cumulative=float(sum(r['gw_rate']*r['duration_myr']*1e6 for r in rows)),
        full_weighted_cumulative=float(sum(r['full_rate']*r['duration_myr']*1e6 for r in rows)),
        excess_weighted_cumulative=float(sum(r['excess_rate']*r['duration_myr']*1e6 for r in rows)))
    save_json('summary.json',summary)
    a=np.genfromtxt(OUT/'paired_rates.csv',delimiter=',',names=True)
    ref=np.genfromtxt(OUT/'paper_fig14_vector.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    smooth=ref[ref['curve']=='polynomial']
    variants=np.genfromtxt(OUT/'variants.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    fullbin=np.isclose(a['duration_myr'],CFG.global_timestep_myr)
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    ax=axes[0]
    for name,color,label in [('full_rate','black','A: full evolution'),('gw_rate','tab:orange','A: matched GW-only'),('excess_rate','tab:blue','A: full - GW (diagnostic)')]:
        ax.loglog(a['z_mid'][fullbin],a[name][fullbin],color=color,label=label)
        ax.plot(a['z_mid'][~fullbin],a[name][~fullbin],'D',mfc='none',color=color)
    ax.loglog(smooth['z'],smooth['rate'],'--',color='tab:green',label='Paper Fig.14 polynomial')
    ax.set(xlabel='Redshift z',ylabel='Rate [yr^-1 halo^-1]',title='Paired initial conditions; diamonds = short final bin')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    ax=axes[1]
    for field,label in [('binary_population','Binary population weight'),('full_rate_per_binary','Full probability / time per binary'),('full_rate','Full halo rate')]:
        ax.loglog(a['z_mid'][fullbin],a[field][fullbin]/a[field][-2],label=label)
    ax.set(xlabel='Redshift z',ylabel='Relative to last complete 200 Myr bin',title='Rate = population x specific probability/time')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(OUT/'rate_decomposition.png',dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5))
    ax.loglog(a['z_mid'][fullbin],a['full_rate'][fullbin],label='Strict / Prada12-HMF')
    for name in np.unique(variants['variant']):
        v=variants[variants['variant']==name]
        ax.loglog(v['z_mid'][:-1],v['rate'][:-1],label=name)
    ax.loglog(smooth['z'],smooth['rate'],'k--',label='Paper Fig.14 polynomial')
    ax.set(xlabel='Redshift z',ylabel='Rate [yr^-1 halo^-1]',title='One model change at a time; complete bins')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(OUT/'model_controls.png',dpi=180); plt.close(fig)
    comparisons=[]
    for z in [.01,.1,1.,5.,8.]:
        row=dict(z=z,paper=10**np.interp(np.log10(z),np.log10(smooth['z'][::-1]),np.log10(smooth['rate'][::-1])))
        for name in ['full_rate','gw_rate','excess_rate']:
            row[name]=10**np.interp(np.log10(z),np.log10(a['z_mid'][::-1]),np.log10(a[name][::-1]))
        comparisons.append(row)
    save_rows('matched_redshift_comparison.csv',comparisons)


def cohort_gw_clock_check():
    """Reproduce B/GW terminal deaths by one uninterrupted isolated integration."""
    metadata=json.loads((BASE/'cohort_aged_gw_config.json').read_text(encoding='utf-8'))
    data=dict(np.load(BASE/'cohort_aged_gw.npz'))
    cfg=bs.BinarySingleConfig(**metadata['config'])
    pop=mc.CohortPopulationConfig(**metadata['population_config'])
    rng=np.random.default_rng(metadata['seed'])
    dist=prior.AppendixDOrbitalDistribution(**metadata['prior'])
    age=float(mc.cosmic_age_myr(0,hc.cosmology_for_model(cfg.halo_model)))-pop.formation_age_myr
    cfg=replace(cfg,global_timestep_myr=age,eccentricity_growth_model='zero')
    state=bs.build_shell_environment(M0,0,cfg)
    direct=np.zeros(data['supplied'].shape[1]); invalid=0
    for k in range(len(data['redshift'])):
        for i in range(len(direct)):
            if data['supplied'][k,i] <= 0:
                continue
            n=pop.initial_samples_per_shell if k==0 else pop.accretion_samples_per_shell
            samples=dist.sample(n,rng)
            isolated=replace(state,environment_density_msun_pc3=np.zeros(len(direct)),
                hard_semimajor_axis_pc=np.full(len(direct),samples.semi_major_axis_pc.max()))
            events=bs.evolve_binary_batch(samples.semi_major_axis_pc,samples.eccentricity,isolated,i,cfg)
            direct[i]+=events.merger_count*data['supplied'][k,i]/n
            invalid+=events.invalid_count
    split=data['outside_merged'].sum(axis=0)+data['merged_weight'].sum(axis=0)
    save_json('cohort_gw_clock_check.json',dict(direct_merged_by_shell=direct.tolist(),
        inherited_merged_by_shell=split.tolist(),difference_by_shell=(direct-split).tolist(),
        invalid=invalid,maximum_weight_difference=float(np.abs(direct-split).max())))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('task',choices=['paired','variants','duration','reference','numerical','summarize','cohort_clock'])
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    {'paired':paired_history,'variants':variants,'duration':duration_controls,
     'reference':paper_reference,'numerical':numerical_check,'summarize':summarize,
     'cohort_clock':cohort_gw_clock_check}[args.task]()
