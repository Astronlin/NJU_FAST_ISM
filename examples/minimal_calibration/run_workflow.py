#!/usr/bin/env python3
"""Run the real Sun17 example through reference, calibration, baseline and table.

Uses the normal FISH commands; omits only spatial gridding. Outputs go to a fresh
workspace and include native fit records and a Doppler-aligned spectral table.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from prepare_viewer import prepare
from fastfish.commands.shell import FishShell
from fastfish.processing.preview import PreviewDataset, reconstruct_calibrated
from fastfish.processing.stages import stage_path
from fastfish.baseline.components import reconstruct_components
from fastfish.rfi.detection import frequency_mask


def summarize(session):
    state=json.loads(Path(session).read_text());root=Path(state['root'])
    ds=PreviewDataset(state['config'],state['selection'],memory_mib=128)
    params=state['params'];cfg=dict(ds.cfg,**params)
    folders={name:stage_path(root,name,ds,params) for name in ['reference','calibration','baseline','resample']}
    summaries={name:json.loads((folder/'summary.json').read_text()) for name,folder in folders.items()}
    record=json.loads((folders['resample']/'processed.json').read_text())[0]
    product=record['products'][0]
    table=np.load(product['data'],mmap_mode='r');rfi=np.load(product['rfi'],mmap_mode='r')
    assert np.isfinite(table).all() and np.isfinite(rfi).all()
    assert np.all((rfi>=0)&(rfi<=1))
    table_manifest=json.loads((root/'table_manifest.json').read_text())
    output_f=np.asarray(table_manifest['frequency_ghz'])
    stem=Path(ds.files[0]['entry']['path']).stem
    calfolder=folders['calibration']/'Sun17_A_01/M01'
    basefolder=folders['baseline']/'Sun17_A_01/M01'
    with np.load(calfolder/(stem+'_calibration.npz')) as p:
        coeff=p['coefficients'][:1].copy();edges=p['reference_edges_ghz'].copy()
    with np.load(basefolder/(stem+'_baseline.npz')) as p:
        components=json.loads(str(p['component_records_json'][0]))
    reports=json.loads((basefolder/(stem+'_fit_reports.json')).read_text())
    masks=json.loads((calfolder/(stem+'_intervals.json')).read_text())
    f=ds.frequency(0);raw=ds.raw(0,0)[:1];scales=ds.scales(0)[:1]
    pol=reconstruct_calibrated(raw,f,coeff,edges)
    residual,models=reconstruct_components(pol,f,scales,cfg,components,rows=slice(0,1),return_components=True)
    v=(1-f*scales[0]/cfg['rest_freq_ghz'])*299792.458
    clean=(abs(v)>cfg['protect_kms'])&(abs(v)<cfg['local_fit_kms'])&~frequency_mask(f,masks[0]['intervals'])
    stats=[]
    for p in range(2):
        x=residual[0,clean,p];centre=np.median(x)
        stats.append(dict(polarization=['AA','BB'][p],first_spectrum_clean_median_K=float(centre),
            first_spectrum_clean_robust_scatter_K=float(1.4826*np.median(abs(x-centre)))))
    summary=dict(dataset='Sun17_A_01 M01 real extract',selected_spectra=len(ds.ra),
        stages={name:dict(spectra=value['spectra'],scan_beams=value['scan_beams']) for name,value in summaries.items()},
        native_channels=len(f),output_shape=list(table.shape),finite_fraction=float(np.isfinite(table).mean()),
        output_velocity_kms=((1-output_f[[0,-1]]/cfg['rest_freq_ghz'])*299792.458).tolist(),
        RFI_intervals_per_chunk=[x['intervals'] for x in masks],
        baseline_mode='ordered components, AA/BB fitted separately',
        baseline_period_at_limit_chunks=sum(bool(x.get('at_search_boundary')) for x in reports),
        baseline_first_spectrum=stats,gridding_performed=False,
        checks_passed=True,note='Finite tables and clean-channel residual statistics do not establish baseline truth inside protected H I or absolute calibration accuracy.')
    (root/'workflow_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Reuse SAVED fits, never fit a second model for the figure.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    display=32
    def mean(a):return a[:len(a)//display*display].reshape(-1,display).mean(axis=1)
    vm=mean(v)
    for p,label in enumerate(['AA','BB']):
        ax=axs[0,p];obs=pol[0,:,p];total=sum((m[0,:,p] for m in models),np.zeros(len(f)))
        ax.plot(vm,mean(obs),color='.6',lw=.7,label='Calibrated')
        if models:ax.plot(vm,mean(models[0][0,:,p]),lw=.9,label='First component')
        ax.plot(vm,mean(total),color='black',lw=1,label='Total baseline')
        centre=float(np.median(obs[abs(v)>700]));ax.set_ylim(centre-3,centre+3)
        ax.set(title=label+' | baseline components; H I peak clipped',ylabel='Temperature (K)',xlabel='LSRK velocity (km/s)')
        ax=axs[1,p]
        if models:
            ax.plot(vm,mean(pol[0,:,p]-models[0][0,:,p]),color='.6',lw=.7,label='After first component')
            local=sum((m[0,:,p] for m in models[1:]),np.zeros(len(f)))
            ax.plot(vm,mean(local),color='#247c43',lw=1,label='Remaining local model')
        ax.plot(vm,mean(residual[0,:,p]),color=['#1764ab','#c24f17'][p],lw=.8,label='Final residual')
        ax.set(xlim=(-600,600),ylim=(-2,2),title=label+' | residual; H I peak clipped',ylabel='Residual (K)',xlabel='LSRK velocity (km/s)')
        for a in [axs[0,p],axs[1,p]]:
            a.axvspan(-cfg['protect_kms'],cfg['protect_kms'],color='.8',alpha=.3,label='Protected H I')
            for lo,hi in masks[0]['intervals']:
                bounds=(1-np.array([lo,hi])/1000*scales[0]/cfg['rest_freq_ghz'])*299792.458
                a.axvspan(min(bounds),max(bounds),color='red',alpha=.18)
            a.grid(alpha=.2);a.legend(fontsize=8)
    fig.savefig(root/'baseline.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    velocity=(1-output_f/cfg['rest_freq_ghz'])*299792.458
    ax.step(velocity,table.mean(axis=0),where='mid',color='#1764ab',lw=.8,label='Mean of 32 spectra; n_b=2')
    ax.axhline(0,color='.5',ls='--',lw=.7);ax.set(xlabel='LSRK velocity (km/s)',ylabel='Baseline-subtracted temperature (K)',title='Doppler-aligned output table; no spatial gridding')
    ax.legend();ax.grid(alpha=.2);fig.savefig(root/'table-spectrum.png',dpi=150);plt.close(fig)
    from calibrate import plot_diagnostics
    pack=ds.reference(0,0,params)
    plot_diagnostics(raw,f,ds.files[0]['coords']['mjd'][:1],pack,dict(pol=pol),root/'calibration.png')
    print(json.dumps(summary,indent=2));print('Table manifest:',root/'table_manifest.json')
    return summary


def run(output,n_t_cal=10,n_t_apply=10,baseline_components=None):
    session=prepare(output)
    shell=FishShell(viewer='none');shell.path=session
    settings=dict(n_t_cal=n_t_cal,n_t_apply=n_t_apply)
    if baseline_components is not None:settings['baseline_components']=baseline_components
    for key,value in settings.items():
        shell.onecmd('set '+key+' '+json.dumps(value))
        if shell.last_error:raise RuntimeError(shell.last_error)
    for command in ['tsyscal run','wait','applycal run','wait','baseline run','wait','table run','wait']:
        print('fish> '+command,flush=True);shell.onecmd(command)
        if shell.last_error:raise RuntimeError(shell.last_error)
    return summarize(session)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--n-t-cal',type=int,default=10)
    parser.add_argument('--n-t-apply',type=int,default=10)
    parser.add_argument('--components-file',help='Optional JSON list of baseline component specifications')
    args=parser.parse_args()
    components=json.loads(Path(args.components_file).read_text()) if args.components_file else None
    run(args.output,args.n_t_cal,args.n_t_apply,components)
