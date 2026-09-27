"""Fig.12 可重复诊断与绘图入口；输出至 pic_and_data/fig12_diagnosis。

python fig12_diagnostics.py --task integration
python fig12_diagnostics.py --task runs --variants strict,strict_tight,gw_only,paper_order
python fig12_diagnostics.py --task sensitivity
python fig12_diagnostics.py --task plot

paper_order 是已有 Fig.19 有效参数诊断口径，不是作者未公开的生成器。
本脚本不进行人口重加权、HMF 积分或 cohort 继承。
"""

import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

import binary_single as bs
import binary_single_montecarlo as mc
import haloconcentration as hc
import halostructure as hs
import p_a_j as prior

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "pic_and_data" / "fig12_diagnosis"
M0 = 1.15e12
SEED = 20260922
COLORS = ["blue", "darkgreen", "red", "darkturquoise", "magenta", "olive", "black", "royalblue", "forestgreen", "firebrick"]


def save_rows(name, rows):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_json(name, data):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


class Figure19DiagnosticPrior:
    def sample(self, n, rng):
        return prior.sample_figure19_orbital_parameters(n, rng)


def extract_reference():
    """从 PDF 线段与刻度定位读取，不用 PNG 像素估计。"""
    import pymupdf
    source = ROOT / "tmp/arxiv_2408.06515v2_source/Merger3body_N1.1_c_v1.pdf"
    page = pymupdf.open(source)[0]
    paths = page.get_drawings()
    x_ticks = [paths[i]["items"][0][1].x for i in [2, 4, 6, 8, 10, 12]]
    y_ticks = [paths[i]["items"][0][1].y for i in [14, 16, 18, 20]]
    x_to_time = np.polyfit(x_ticks, np.arange(6) * 2500.0, 1)
    y_to_log_count = np.polyfit(y_ticks, [2., 3., 4., 5.], 1)
    curves = [p for p in paths if len(p["items"]) > 50]
    if len(curves) != 10:
        raise ValueError("原始 PDF 路径结构发生变化，需要重新识别曲线。")
    table = np.empty((67, 10))
    rows = []
    for shell, path in enumerate(curves):
        segments = []
        for item in path["items"]:
            if item[0] == "l" and abs(item[1].y - item[2].y) < 1e-7 and item[2].x > item[1].x:
                t_left = float(np.polyval(x_to_time, item[1].x))
                t_right = float(np.polyval(x_to_time, item[2].x))
                count = float(10**np.polyval(y_to_log_count, item[1].y))
                segments.append((t_left, count))
                rows.append(dict(shell=shell + 1, t_label_left_myr=t_left, t_label_right_myr=t_right, cumulative_from_vector=count))
        segments = np.array(segments)
        index = np.searchsorted(segments[:, 0], np.arange(67) * 200.0 + 0.001, side="right") - 1
        table[:, shell] = segments[index, 1]
    save_rows("paper_vector_segments.csv", rows)
    np.savetxt(OUTPUT / "paper_reference.csv", np.column_stack((np.arange(67)*200.0, np.arange(1,68)*200.0,table)), delimiter=",", header="t_label_myr,t_completed_if_first_bin_at_zero_myr," + ",".join(f"R{i}" for i in range(1,11)), comments="")
    save_json("paper_reference_metadata.json", dict(source=str(source.relative_to(ROOT)),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),x_to_time_coefficients=x_to_time.tolist(),y_to_log_count_coefficients=y_to_log_count.tolist(),first_label_myr=0,last_label_myr=13200,last_horizontal_end_myr=13400,number_of_bins=67,endpoint_counts=table[-1].tolist(),explanation="原图在 t=0 已有第一片计数；t_completed=t_label+200 仅用于时间片对齐，保存原标签不覆盖。"))
    return table


def integration_comparison():
    cfg = bs.BinarySingleConfig(integration_method="adaptive_log_a")
    shell = bs.build_shell_environment(M0, 12, cfg)
    samples = prior.AppendixDOrbitalDistribution().sample(500, np.random.default_rng(12))
    a, e = samples.semi_major_axis_pc, samples.eccentricity
    reference = bs.evolve_binary_batch(a, e, shell, 0, replace(cfg,integration_rtol=1e-10,integration_log_j2_atol=1e-12,integration_time_atol_myr=1e-10))
    rows = []
    for method, value in [("euler",2.0),("euler",0.2),("euler",0.02),("adaptive_log_a",1e-4),("adaptive_log_a",1e-6),("adaptive_log_a",1e-8)]:
        conf = replace(cfg,integration_method=method,**({"local_timestep_myr":value} if method=="euler" else {"integration_rtol":value}))
        start = time.perf_counter()
        event = bs.evolve_binary_batch(a,e,shell,0,conf)
        both = event.merged_mask & reference.merged_mask
        rows.append(dict(method=method,step_or_rtol=value,hard=event.hard_count,merged=event.merger_count,invalid=event.invalid_count,classification_difference=int(np.count_nonzero(event.merged_mask != reference.merged_mask)),max_merger_time_error_myr=float(np.max(np.abs(event.merger_time_myr[both]-reference.merger_time_myr[both]))),seconds=time.perf_counter()-start))
        print(rows[-1],flush=True)
    save_rows("integration_convergence.csv",rows)
    # Independent integration algorithm, with the same transformed ODE.
    dop_rows = []
    cutoff = cfg.merger_radius_gm_c2 * bs.G_SI * 60 * bs.SOLAR_MASS_KG / bs.C_SI**2 / bs.PARSEC_M
    stop_event = lambda u,y: y[0] - cfg.global_timestep_myr
    stop_event.terminal = True
    stop_event.direction = 1
    nominal = bs.evolve_binary_batch(a,e,shell,0,cfg)
    for i in np.flatnonzero(reference.hard_mask)[::5]:
        rhs = lambda u,y: bs.log_a_orbital_rhs(u,y,a[i],shell.environment_density_msun_pc3[0],shell.velocity_dispersion_km_s[0],shell.hard_semimajor_axis_pc[0],cfg)
        sol = solve_ivp(rhs,(0,np.log(a[i]/cutoff)),[0,np.log((1-e[i])*(1+e[i]))],method="DOP853",rtol=1e-10,atol=[1e-10,1e-12],events=stop_event,max_step=.1)
        is_merged = sol.success and len(sol.t_events[0])==0
        actual_a = 0.0 if is_merged else a[i]*np.exp(-sol.t[-1])
        actual_e = np.sqrt(-np.expm1(sol.y[1,-1]))
        dop_rows.append(dict(sample_index=int(i),success=sol.success,reference_merged=bool(is_merged),batch_merged=bool(nominal.merged_mask[i]),reference_time_myr=float(sol.y[0,-1]),batch_time_myr=float(nominal.merger_time_myr[i]) if nominal.merged_mask[i] else 200.0,a_relative_error=0.0 if is_merged else float(abs(nominal.final_state.semi_major_axis_pc[i]/actual_a-1)),e_absolute_error=float(abs(nominal.final_state.eccentricity[i]-actual_e))))
    save_rows("dop853_comparison.csv",dop_rows)
    # Chain rule against the existing SI-unit Eq.19/25 implementation.
    use = reference.hard_mask
    physical = bs.orbital_derivatives(a[use],e[use],shell.environment_density_msun_pc3[0],shell.velocity_dispersion_km_s[0],shell.hard_semimajor_axis_pc[0])
    du_dt = -physical.da_pc_myr/a[use]
    transformed = bs.log_a_orbital_rhs(np.zeros(use.sum()),np.column_stack((np.zeros(use.sum()),np.log((1-e[use])*(1+e[use])))),a[use],shell.environment_density_msun_pc3[0],shell.velocity_dispersion_km_s[0],shell.hard_semimajor_axis_pc[0],cfg)
    direct = np.column_stack((1/du_dt,-2*e[use]*physical.de_per_myr/((1-e[use])*(1+e[use]))/du_dt))
    # Circular GW solution has a^4=a0^4-4 beta*t.
    beta = 64/5*bs.G_SI**3*(60*bs.SOLAR_MASS_KG)*(30*bs.SOLAR_MASS_KG)**2/bs.C_SI**5*bs.MYR_S/bs.PARSEC_M**4
    circular_a = (4*beta*np.array([100.,500.])+cutoff**4)**.25
    zero_shell = replace(shell,environment_density_msun_pc3=np.zeros(10))
    circular = bs.evolve_binary_batch(circular_a,np.zeros(2),zero_shell,0,cfg)
    cutoff_rows = []
    for radius in [3.,6.,60.]:
        event=bs.evolve_binary_batch(a,e,shell,0,replace(cfg,merger_radius_gm_c2=radius))
        both=event.merged_mask & nominal.merged_mask
        cutoff_rows.append(dict(radius_gm_c2=radius,merged=event.merger_count,invalid=event.invalid_count,max_delta_time_myr=float(np.max(abs(event.merger_time_myr[both]-nominal.merger_time_myr[both])))))
    save_rows("termination_radius.csv",cutoff_rows)
    save_json("integration_summary.json",dict(chain_rule_max_relative_error=float(np.max(abs(transformed-direct)/np.maximum(abs(direct),1e-100))),dop853_orbits=len(dop_rows),dop853_success=all(r['success'] for r in dop_rows),dop853_classification_differences=sum(r['reference_merged']!=r['batch_merged'] for r in dop_rows),dop853_max_time_error_myr=max(abs(r['reference_time_myr']-r['batch_time_myr']) for r in dop_rows),dop853_max_survivor_a_relative_error=max(r['a_relative_error'] for r in dop_rows),circular_merger_time_myr=float(circular.merger_time_myr[0]),circular_expected_merger_time_myr=100.,circular_survivor_a_relative_error=float(abs(circular.final_state.semi_major_axis_pc[1]/(circular_a[1]**4-4*beta*200)**.25-1)),circular_invalid=circular.invalid_count))


def run_history(variant,n,seed):
    cfg = bs.BinarySingleConfig(sample_count_per_shell_and_step=n,integration_method="adaptive_log_a")
    distribution = prior.AppendixDOrbitalDistribution()
    if variant == "strict_tight":
        cfg = replace(cfg,integration_rtol=1e-8,integration_log_j2_atol=1e-11,integration_time_atol_myr=1e-10)
    elif variant == "gw_only":
        cfg = replace(cfg,ksi=0.0)
    elif variant == "paper_order":
        distribution = Figure19DiagnosticPrior()
    elif variant == "taper":
        cfg = replace(cfg,eccentricity_growth_model="sesana_taper")
    elif variant != "strict":
        raise ValueError(variant)
    start=time.perf_counter()
    def progress(i,total):
        if i%10==0 or i==total:
            print(f"{variant}: {i}/{total}, {time.perf_counter()-start:.1f}s",flush=True)
    result=mc.run_raw_shell_monte_carlo(M0,cfg,distribution,seed,5000,progress)
    arrays=dict(time=result.elapsed_time_edges_myr,redshift=result.redshift_edges,sampled=result.sampled_per_step,hard=result.hard_per_step,merged=result.mergers_per_step,invalid=result.invalid_per_step,k_above=result.k_above_calibration_per_step,cumulative=result.cumulative_mergers)
    OUTPUT.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(OUTPUT/f"{variant}.npz",**arrays)
    save_json(f"{variant}_config.json",dict(config=asdict(cfg),samples=n,seed=seed,chunk_size=5000,halo_M0_msun=M0,scheme="A independent redraw",prior=asdict(distribution) if variant!="paper_order" else dict(name="Figure19DiagnosticPrior",j_scale=prior.FIGURE19_EFFECTIVE_ANGULAR_MOMENTUM_SCALE,a_scale_pc=prior.FIGURE19_EFFECTIVE_SEMIMAJOR_AXIS_SCALE_PC),elapsed_seconds=time.perf_counter()-start,invalid=int(result.invalid_per_step.sum())))
    columns=["elapsed_time_myr","redshift"]+[f"R{i}_cumulative" for i in range(1,11)]
    np.savetxt(OUTPUT/f"{variant}_cumulative.csv",np.column_stack((arrays['time'],arrays['redshift'],arrays['cumulative'])),delimiter=",",header=",".join(columns),comments="")
    rows=[]
    for t in range(result.time_step_count):
        for i in range(10):
            rows.append(dict(step=t+1,shell=i+1,t_end_myr=arrays['time'][t+1],redshift_start=arrays['redshift'][t],**{key:int(arrays[key][t,i]) for key in ['sampled','hard','merged','invalid','k_above']}))
    save_rows(f"{variant}_counts.csv",rows)
    print(variant,"endpoint",arrays['cumulative'][-1],"invalid",arrays['invalid'].sum(),flush=True)


def frozen_sensitivity(n=5000):
    cfg=bs.BinarySingleConfig(integration_method="adaptive_log_a")
    samples=prior.AppendixDOrbitalDistribution().sample(n,np.random.default_rng(SEED))
    rows=[]
    for z in [12.,3.,0.]:
        base=bs.build_shell_environment(M0,z,cfg)
        for name in ["baseline","K_taper","K_zero","velocity_outer_boundary","mu2","ludlow16","prada_rho_eq"]:
            conf=replace(cfg,eccentricity_growth_model={"K_taper":"sesana_taper","K_zero":"zero"}.get(name,"sesana_endpoint"))
            shell=base
            draw=samples
            if name=="mu2":
                conf=replace(conf,velocity_mu=2.)
                shell=bs.build_shell_environment(M0,z,conf)
            elif name=="ludlow16":
                conf=replace(conf,halo_model="ludlow16")
                shell=bs.build_shell_environment(M0,z,conf)
            elif name=="velocity_outer_boundary":
                conf=replace(conf,velocity_radius_strategy="outer_boundary")
                shell=bs.build_shell_environment(M0,z,conf)
            elif name=="prada_rho_eq":
                cosmology=hc.cosmology_for_model(cfg.halo_model)
                density=prior.matter_density_at_equality_msun_pc3(omega_m0=cosmology['omega_m0'],h=cosmology['h'])
                draw=replace(prior.AppendixDOrbitalDistribution(),rho_eq_msun_pc3=density).sample(n,np.random.default_rng(SEED))
            for i in [0,1,4,9]:
                event=bs.evolve_binary_batch(draw.semi_major_axis_pc,draw.eccentricity,shell,i,conf)
                rows.append(dict(redshift=z,shell=i+1,variant=name,samples=n,rho=shell.environment_density_msun_pc3[i],velocity=shell.velocity_dispersion_km_s[i],a_h=shell.hard_semimajor_axis_pc[i],hard=event.hard_count,merged=event.merger_count,invalid=event.invalid_count))
            print("sensitivity",z,name,flush=True)
    save_rows("frozen_sensitivity.csv",rows)


def prior_and_paired_diagnostics():
    """记录初态尾部与同一轨道的环境/GW 事件交集。"""
    n = 1_000_000
    rows = []
    prior_summary = {}
    for name, distribution in [("strict", prior.AppendixDOrbitalDistribution()), ("paper_order", Figure19DiagnosticPrior())]:
        draw = distribution.sample(n, np.random.default_rng(22))
        a, e = draw.semi_major_axis_pc, draw.eccentricity
        prior_summary[name] = dict(samples=n, seed=22, mean_a_pc=float(a.mean()), median_a_pc=float(np.median(a)), mean_e=float(e.mean()), median_e=float(np.median(e)), fraction_e_above_09=float(np.mean(e > .9)))
        for ah in [1.07e-2, 1.79e-5, 3.17e-6, 1e-4]:
            mask = a <= ah
            rows.append(dict(prior=name, a_h_pc=ah, n_sample=n, n_hard=int(mask.sum()), hard_percent=100*mask.mean(), hard_e_median=float(np.median(e[mask]))))
    save_rows('prior_hard_tail.csv', rows)
    save_json('prior_summary.json', prior_summary)
    config = bs.BinarySingleConfig()
    samples = prior.AppendixDOrbitalDistribution().sample(5000, np.random.default_rng(SEED))
    rows = []
    for z in [12., 3., 0.]:
        shell = bs.build_shell_environment(M0, z, config)
        no_environment = replace(shell, environment_density_msun_pc3=np.zeros(10))
        for i in range(10):
            full = bs.evolve_binary_batch(samples.semi_major_axis_pc, samples.eccentricity, shell, i, config)
            gw = bs.evolve_binary_batch(samples.semi_major_axis_pc, samples.eccentricity, no_environment, i, config)
            rows.append(dict(redshift=z, shell=i+1, samples=5000, hard=full.hard_count, full_only=int(np.count_nonzero(full.merged_mask & ~gw.merged_mask)), gw_only=int(np.count_nonzero(gw.merged_mask & ~full.merged_mask)), both=int(np.count_nonzero(full.merged_mask & gw.merged_mask)), neither_hard=int(np.count_nonzero(full.hard_mask & ~full.merged_mask & ~gw.merged_mask)), invalid=full.invalid_count+gw.invalid_count))
    save_rows('paired_event_outcomes.csv', rows)


def plot_results():
    paper=extract_reference()
    names=[n for n in ['strict','strict_tight','gw_only','paper_order','taper'] if (OUTPUT/f'{n}.npz').exists()]
    runs={n:dict(np.load(OUTPUT/f'{n}.npz')) for n in names}
    configs={n:json.loads((OUTPUT/f'{n}_config.json').read_text(encoding='utf-8')) for n in names}
    rows=[]
    old=np.genfromtxt(ROOT/'pic_and_data/fig12_schemeA_strict_prada12_hmf_counts.csv',delimiter=',',names=True)
    for i in range(10):
        row=dict(shell=i+1,paper_endpoint=float(paper[-1,i]),old_euler_raw=int(np.sum(old[f'R{i+1}_merged'])),old_euler_invalid=int(np.sum(old[f'R{i+1}_invalid'])))
        for n in names:
            data=runs[n]
            sample=configs[n]['samples']
            counts=data['merged'][:,i]
            raw=int(counts.sum())
            scaled=raw*2e6/sample
            # Independent Bernoulli trials in each time bin; scaled Monte Carlo SE.
            se=np.sqrt(np.sum(counts*(1-counts/sample)*sample/(sample-1)))*2e6/sample
            row.update({f'{n}_raw':raw,f'{n}_scaled':scaled,f'{n}_MC_se':se,f'{n}_zero_count_upper95_scaled':float(-np.log(.05)*2e6/sample) if raw==0 else '',f'{n}_ratio_to_paper':scaled/paper[-1,i],f'{n}_hard':int(data['hard'][:,i].sum()),f'{n}_invalid':int(data['invalid'][:,i].sum()),f'{n}_first67_scaled':float(data['merged'][:67,i].sum()*2e6/sample)})
        rows.append(row)
    save_rows('shell_summary.csv',rows)
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    panels=[('strict','Adaptive; strict joint'),('gw_only','GW only; same hard cut'),('paper_order','Fig.19 effective-prior diagnostic')]
    for ax,(name,title) in zip([axes[0,0],axes[0,1],axes[1,0]],panels):
        if name not in runs:
            continue
        data=runs[name]
        scale=2e6/configs[name]['samples']
        for i,color in enumerate(COLORS):
            ax.step(data['time'],np.where(data['cumulative'][:,i]>0,data['cumulative'][:,i]*scale,np.nan),where='post',color=color,label=f'R{i+1}')
            ax.step(np.arange(1,68)*200,paper[:,i],where='post',color=color,ls=':',alpha=.6,lw=1)
        ax.set(yscale='log',xlabel='Elapsed time [Myr]',ylabel='Cumulative mergers (scaled estimate)',title=title)
        ax.grid(alpha=.2)
    axes[0,0].legend(ncol=2,fontsize=8)
    axes[1,0].text(.03,.03,'Dotted: paper vectors, first bin aligned at 200 Myr',transform=axes[1,0].transAxes,fontsize=8)
    ax=axes[1,1]
    shell=np.arange(1,11)
    ax.plot(shell,[r['old_euler_raw']*400/r['paper_endpoint'] for r in rows],'x--',label='Old Euler strict')
    for name in names:
        if name=='strict_tight': continue
        ax.plot(shell,[r[f'{name}_ratio_to_paper'] for r in rows],'o-',label=name)
    ax.axhline(1,color='black',lw=1)
    ax.set(yscale='log',xlabel='Shell',ylabel='Estimate / paper endpoint',title='Amplitude mismatch after numerical repair',xticks=shell)
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(OUTPUT/'fig12_diagnostic_comparison.png',dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.5))
    data=runs['strict']
    for i,color in enumerate(COLORS):
        ax.step(data['time'],np.where(data['cumulative'][:,i]>0,data['cumulative'][:,i],np.nan),where='post',color=color,label=f'R{i+1}')
    ax.set(yscale='log',xlabel='Elapsed time [Myr]',ylabel='Raw cumulative mergers',title=f"Stable integration | strict joint | {configs['strict']['samples']:,} draws/shell/bin")
    ax.legend(ncol=2)
    ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(OUTPUT/'fig12_corrected_raw.png',dpi=180)
    plt.close(fig)
    if 'strict_tight' in runs:
        save_json('full_history_convergence.json',dict(max_bin_count_difference=int(np.max(abs(runs['strict']['merged']-runs['strict_tight']['merged']))),number_of_changed_bins=int(np.count_nonzero(runs['strict']['merged']!=runs['strict_tight']['merged'])),final_count_difference=(runs['strict']['cumulative'][-1]-runs['strict_tight']['cumulative'][-1]).tolist()))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',choices=['integration','runs','sensitivity','prior','plot'],required=True)
    parser.add_argument('--samples',type=int,default=5000)
    parser.add_argument('--seed',type=int,default=SEED)
    parser.add_argument('--variants',default='strict,strict_tight,gw_only,paper_order')
    args=parser.parse_args()
    if args.task=='integration': integration_comparison()
    elif args.task=='runs':
        for variant in args.variants.split(','): run_history(variant,args.samples,args.seed)
    elif args.task=='sensitivity': frozen_sensitivity(args.samples)
    elif args.task=='prior': prior_and_paired_diagnostics()
    else: plot_results()
