"""Server-backed IC443 preview windows. Only Run-all actions write products."""
from pathlib import Path
import argparse
import fcntl
import html
import json
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from fastfish.processing.preview import PreviewDataset
from fastfish.baseline.windows import parse_windows,format_windows
from fastfish.viewer.baseline_workspace import attach_workspace
from fastfish.processing.stages import stage_path
from fastfish.processing.pipeline import dump,read,source_hash
from fastfish.rfi.detection import C


def display_data(x,columns,start=None,end=None,limit=2200):
    """Display-only means plus min/max envelopes; native samples on zoom."""
    x=np.asarray(x);ids=np.arange(len(x))
    if start is not None and end is not None:
        ids=ids[(x>=min(start,end))&(x<=max(start,end))]
    science=[k for k in columns if k in ('AA','BB','average','observed','residual','residual_AA','residual_BB','mean')]
    if len(ids)>2*limit:
        sizes=np.full(limit,len(ids)//limit,dtype=int);sizes[:len(ids)%limit]+=1
        starts=np.r_[0,np.cumsum(sizes)[:-1]]
        result=dict(x=np.add.reduceat(x[ids],starts)/sizes)
        for k,value in columns.items():
            y=np.asarray(value)[ids];valid=~np.isnan(y)
            count=np.add.reduceat(valid.astype(int),starts)
            total=np.add.reduceat(np.where(valid,y,0.),starts)
            result[k]=np.divide(total,count,out=np.full(limit,np.nan),where=count>0)
            if k in science:
                result[k+'_low']=np.fmin.reduceat(y,starts)
                result[k+'_high']=np.fmax.reduceat(y,starts)
    else:
        result=dict(x=x[ids],**{k:np.asarray(v)[ids] for k,v in columns.items()})
        for k in science:result[k+'_low']=result[k];result[k+'_high']=result[k]
    return result


def build_document(doc,config,selection,root,initial_params=None,initial_index=None,on_settings=None,allow_batches=True):
    from bokeh.events import Tap,SelectionGeometry
    from bokeh.layouts import column,row
    from bokeh.models import (ColumnDataSource,Range1d,LinearAxis,Button,Div,Spinner,
        Select,CheckboxGroup,Panel,Tabs,BoxAnnotation,CDSView,BooleanFilter,TextAreaInput)
    from bokeh.plotting import figure
    from bokeh.transform import linear_cmap
    from bokeh.palettes import Category10
    ds=PreviewDataset(config,selection);root=Path(root)
    if initial_params:ds.cfg.update(initial_params)
    arguments=doc.session_context.request.arguments if doc.session_context is not None else {}
    initial_view=arguments.get('view',[b'baseline' if on_settings else b'reference'])[0]
    if isinstance(initial_view,bytes):initial_view=initial_view.decode()
    if initial_view=='where':
        footprint=figure(width=1000,height=680,title='FAST · '+str(ds.cfg.get('name','observation'))+' · pointing footprint',x_axis_label='RA (deg)',y_axis_label='Dec (deg)',tools='pan,wheel_zoom,box_zoom,reset',active_scroll='wheel_zoom',output_backend='webgl')
        footprint.circle(ds.ra,ds.dec,size=2,alpha=.3)
        doc.add_root(column(Div(text='<h2>Pointing footprint</h2><p>{} spectra; coordinate metadata only. No SCI spectra read or calibrated.</p>'.format(len(ds.ra))),footprint));doc.title='FISH · footprint';return
    executor=ThreadPoolExecutor(max_workers=1);state=dict(generation=0,index=ds.closest(ds.cfg['header']['CRVAL1'],ds.cfg['header']['CRVAL2']),point=None,region=None,region_result=None,closed=False,job=None,chain=None)
    if initial_index is not None:state['index']=max(0,min(len(ds.ra)-1,int(initial_index)))
    title=Div(text='<h1 style="margin:0">FAST · '+html.escape(str(ds.cfg.get('name','observation')))+' </h1><p>Full native band retained for fitting</p>',sizing_mode='stretch_width')
    status=Div(text='Loading the selected chunk…',sizing_mode='stretch_width',name='preview_status')
    note=Div(text='Wide views use display-only means and min/max envelopes; zoom reveals native channels. Preview is read-only. Click a pointing centre; wheel/box zoom and pan are available. RFI values are retained and marked. Each broad period uses the same ≤32-spectrum context as Run all.',sizing_mode='stretch_width')
    scan=Select(title='Scan',value='All',options=['All']+ds.cfg['scans'],width=230,name='scan_filter')
    beams=sorted(set(x['beam'] for x in ds.files));beam=Select(title='Beam',value='All',options=['All']+beams,width=120,name='beam_filter')
    index=Spinner(title='Spectrum index',low=0,high=len(ds.ra)-1,step=1,value=state['index'],width=160,name='spectrum_index')
    selection_info=Div(text=str(len(ds.ra))+' selected spectra · catalogue '+ds.cfg['rfi_profile'],sizing_mode='stretch_width')
    arguments=doc.session_context.request.arguments if doc.session_context is not None else {}
    compact=arguments.get('compact',[b'0'])[0] in (b'1','1')
    layout_choice=Select(title='Window layout',value='Side by side',options=['Side by side','Stacked'],width=180,name='window_layout')
    controls=column(row(scan,beam,index),row(layout_choice),selection_info,sizing_mode='stretch_width')
    splits=[];maps=[]
    def split(left,right):
        wrap=column(sizing_mode='stretch_width')
        splits.append((wrap,left,right))
        return wrap
    map_source=ColumnDataSource(dict(x=ds.ra,y=ds.dec,group=ds.groups),name='pointing_centres')
    marker=ColumnDataSource(dict(x=[ds.ra[state['index']]],y=[ds.dec[state['index']]]),name='selected_point')
    visible=BooleanFilter(booleans=[True]*len(ds.ra));view=CDSView(source=map_source,filters=[visible]);active=np.ones(len(ds.ra),dtype=bool)
    protect=Spinner(title='Protected H I half-width ± km/s',value=ds.cfg['protect_kms'],low=0,high=500,step=25,width=230,name='protect_kms')
    local=Spinner(title='Local linear half-width ± km/s',value=ds.cfg['local_fit_kms'],low=550,high=3000,step=50,width=230,name='local_fit_kms')
    baseline=CheckboxGroup(labels=['Subtract baseline'],active=[0] if ds.cfg.get('baseline',True) else [],width=170,name='baseline_enabled')
    mode=Select(title='Spectral processing',value=ds.cfg['spectral_mode'],options=['rebin','gaussian'],width=180,name='spectral_mode')
    nb=Spinner(title='n_b · native channels/bin',value=ds.cfg['n_b'],low=1,high=128,step=1,width=190,name='n_b',disabled=ds.cfg['spectral_mode']!='rebin')
    ng=Spinner(title='n_g · final Gaussian FWHM',value=ds.cfg['n_g'],low=1,high=128,step=.5,width=210,disabled=ds.cfg['spectral_mode']!='gaussian',name='n_g')
    nf=Spinner(title='Assumed native Gaussian FWHM',value=ds.cfg['native_fwhm_chan'],low=.1,high=10,step=.1,width=230,disabled=ds.cfg['spectral_mode']!='gaussian',name='native_fwhm_chan')
    vmin=Spinner(title='Product v_min · km/s',value=ds.cfg['velocity_range_kms'][0],low=-3000,high=2999,step=25,width=180,name='velocity_min')
    vmax=Spinner(title='Product v_max · km/s',value=ds.cfg['velocity_range_kms'][1],low=-2999,high=3000,step=25,width=180,name='velocity_max')
    validation=CheckboxGroup(labels=['Optional direct-gridding validation'],active=[],width=300,name='optional_validation')
    make_cube=CheckboxGroup(labels=['Proceed to combined cube after resampling'],active=[0],width=350,name='make_cube')
    ref_channel=Spinner(title='Reference diagnostic frequency · MHz',value=1420.,low=1400,high=1440,step=.122,width=290,name='reference_frequency')
    nt_cal=Spinner(title='n_t_cal · ON−OFF timestamps (1 = off)',value=ds.cfg.get('n_t_cal',10),low=1,high=3000,step=1,width=360,name='n_t_cal',disabled=not (ds.cfg.get('extracted_cal_root') and ds.cfg.get('tcal_path')))
    nt=Spinner(title='n_t_apply · Tsys-reference timestamps (1 = off)',value=ds.cfg.get('n_t_apply',10),low=1,high=3000,step=1,width=360,name='n_t_apply')
    batch_status=Div(text='Batch idle. No reduction has been launched.',sizing_mode='stretch_width',name='batch_status')
    previous=Button(label='← Previous spectrum',width=170,name='previous_spectrum')
    next_=Button(label='Next spectrum →',width=170,name='next_spectrum')
    refit=Button(label='Refit selected preview',button_type='primary',width=190,name='refit_preview')
    undo=Button(label='Undo window edit',width=160,disabled=True,name='undo_windows')
    clear=Button(label='Clear extra exclusions',width=180,name='clear_exclusions')
    draw_mode=Select(title='Draw a velocity interval on either spectrum',value='Inspect',options=['Inspect','Add exclusion','Set H I protection','Set local fit range'],width=350,name='draw_baseline_window')
    exclusions=TextAreaInput(title='Extra baseline exclusions · lo,hi km/s per line',value=format_windows(ds.cfg.get('baseline_exclusions_kms',[])),value_input=format_windows(ds.cfg.get('baseline_exclusions_kms',[])),max_length=4096,rows=2,width=350,name='baseline_exclusions')
    compare=CheckboxGroup(labels=['Overlay previous corrected mean (same spectrum)'],active=[],width=390,name='compare_previous_fit')
    workspace_help=Div(text='Draw mode selects the operation. Window edits stay in this session; Run all applies them to the full dataset. Purple marks extra fit exclusions; channels are retained.',sizing_mode='stretch_width',name='workspace_help')
    state.update(comparison=None)
    buttons=[]

    def params():
        return dict(baseline_exclusions_kms=parse_windows(exclusions.value_input),n_t_apply=nt.value,n_t_cal=nt_cal.value,baseline=bool(baseline.active),protect_kms=float(protect.value),local_fit_kms=float(local.value),
            spectral_mode=mode.value,n_b=int(nb.value),n_g=float(ng.value),native_fwhm_chan=float(nf.value),velocity_range_kms=[float(vmin.value),float(vmax.value)])

    def pointing_map(label,rectangle=False):
        tools='pan,wheel_zoom,box_zoom,reset,tap'+(',box_select' if rectangle else '')
        p=figure(width=440,height=530,title=label,x_axis_label='RA (deg)',y_axis_label='Dec (deg)',tools=tools,active_scroll='wheel_zoom',output_backend='webgl',sizing_mode='stretch_width')
        p.circle('x','y',source=map_source,view=view,size=3,alpha=.25,color=linear_cmap('group',Category10[5],0,4))
        p.circle('x','y',source=marker,size=11,color='#e45756',line_color='white')
        p.on_event(Tap,lambda event:clicked(event.x,event.y))
        if rectangle:p.on_event(SelectionGeometry,rectangle_changed)
        maps.append(p)
        return p

    def plot(title,ylabel,height=290,xlabel='Native frequency (MHz)',limits=None):
        p=figure(width=940,height=height,title=title,x_axis_label=xlabel,y_axis_label=ylabel,
            tools='pan,wheel_zoom,box_zoom,reset',active_scroll='wheel_zoom',x_range=Range1d(*(limits or (1404.77,1436.03))),sizing_mode='stretch_width')
        return p

    ref_time=plot('Reference time filter · selected frequency','Temperature (K)',260,'Seconds from selected spectrum',(-35,35))
    ref_freq=plot('Reference frequency spectrum · selected timestamp','Temperature (K)',260)
    rt_source=ColumnDataSource(dict(x=[],rawAA=[],smoothAA=[],rawBB=[],smoothBB=[]))
    rf_source=ColumnDataSource(dict(x=[],rawAA=[],smoothAA=[],rawBB=[],smoothBB=[]))
    for p,s in [(ref_time,rt_source),(ref_freq,rf_source)]:
        for key,color,dash in [('rawAA','#8cb5d7','solid'),('smoothAA','#1f77b4','solid'),('rawBB','#edbc8f','solid'),('smoothBB','#d97920','solid')]:
            p.line('x',key,source=s,color=color,line_width=2 if key.startswith('smooth') else 1,legend_label=key)
        p.legend.click_policy='hide'
    cal_plot=plot('Calibrated AA, BB and their mean','Temperature (K)',530)
    cal_source=ColumnDataSource(dict(x=[],AA=[],BB=[],average=[],AA_low=[],AA_high=[],BB_low=[],BB_high=[],average_low=[],average_high=[]),name='calibration_spectrum')
    for key,color in [('AA','#1f77b4'),('BB','#e58a2b'),('average','#202b3b')]:
        cal_plot.varea(x='x',y1=key+'_low',y2=key+'_high',source=cal_source,fill_color=color,fill_alpha=.05)
        cal_plot.line('x',key,source=cal_source,color=color,line_width=1.5,legend_label=key)
    cal_plot.legend.click_policy='hide'
    baseline_pol=Select(title='Baseline components to inspect',value='Mean',options=['Mean','AA','BB'],width=260,name='baseline_polarization')
    upper=plot('Observed spectrum and baseline components','Temperature (K)',300,'LSRK radio velocity (km/s)',(-3300,3300))
    upper.y_range=Range1d(35,45);upper.extra_y_ranges={'correction':Range1d(-1,1)};upper.add_layout(LinearAxis(y_range_name='correction',axis_label='Local correction (K)'),'right')
    upper_source=ColumnDataSource(dict(x=[],observed=[],broad=[],local=[],total=[],observed_low=[],observed_high=[]),name='baseline_components')
    upper.varea(x='x',y1='observed_low',y2='observed_high',source=upper_source,fill_color='#777777',fill_alpha=.07)
    for key,color in [('observed','#777777'),('broad','#dc8c20'),('total','#723fa1')]:upper.line('x',key,source=upper_source,color=color,legend_label=key,line_width=1.6)
    upper.line('x','local',source=upper_source,color='#299d69',y_range_name='correction',legend_label='local · right axis')
    upper.legend.click_policy='hide'
    lower=plot('Baseline-subtracted spectrum · initial ±3 K view clips the H I peak','Residual (K)',260,'LSRK radio velocity (km/s)',(-3300,3300))
    lower.y_range=Range1d(-3,3)
    residual_keys=('residual','residual_AA','residual_BB')
    lower_source=ColumnDataSource(dict(x=[],**{key:[] for k in residual_keys for key in (k,k+'_low',k+'_high')}),name='baseline_residual')
    for key,color,label in [('residual_AA','#1f77b4','AA corrected'),('residual_BB','#e58a2b','BB corrected'),('residual','#202b3b','Corrected mean')]:
        lower.varea(x='x',y1=key+'_low',y2=key+'_high',source=lower_source,fill_color=color,fill_alpha=.04)
        lower.line('x',key,source=lower_source,color=color,line_width=1.8 if key=='residual' else 1,legend_label=label)
    lower.legend.click_policy='hide'
    previous_source=ColumnDataSource(dict(x=[],y=[]),name='previous_corrected_mean')
    previous_line=lower.line('x','y',source=previous_source,color='#aa43a2',line_dash='dashed',line_width=2,visible=False,legend_label='Previous mean')
    manual_boxes=[]
    for plot_ in (upper,lower):
        boxes=[]
        for _ in range(64):
            box=BoxAnnotation(fill_color='#9c50bf',fill_alpha=.12,visible=False);plot_.add_layout(box);boxes.append(box)
        manual_boxes.append(boxes)
    protected=[];local_boxes=[]
    for p in (upper,lower):
        box=BoxAnnotation(left=-250,right=250,fill_color='#969da6',fill_alpha=.15);p.add_layout(box);protected.append(box)
        for a,b in [(-1500,-250),(250,1500)]:
            box=BoxAnnotation(left=a,right=b,fill_color='#42a675',fill_alpha=.06);p.add_layout(box);local_boxes.append(box)
    region_plot=plot('Draw a rectangle to preview its mean · no spatial gridding','Mean temperature (K)',530,'LSRK radio velocity (km/s)',(-250,250))
    region_plot.extra_y_ranges={'fraction':Range1d(0,1)};region_plot.add_layout(LinearAxis(y_range_name='fraction',axis_label='Mean marked contribution'),'right')
    region_source=ColumnDataSource(dict(x=[],mean=[],marked=[],mean_low=[],mean_high=[]),name='region_mean');region_plot.varea(x='x',y1='mean_low',y2='mean_high',source=region_source,fill_color='#1f77b4',fill_alpha=.07);region_plot.line('x','mean',source=region_source,color='#1f77b4',legend_label='Mean spectrum')
    region_plot.line('x','marked',source=region_source,color='#dc544d',y_range_name='fraction',legend_label='Marked fraction',alpha=.6)
    region_label=Div(text='Rectangle preview: choose a region with the box-selection tool.',sizing_mode='stretch_width',name='region_status')
    refine=Button(label='Preview entire rectangle',button_type='default',width=230,name='refine_rectangle')
    cancel=Button(label='Cancel preview',width=150,name='cancel_preview')

    def redraw():
        point=state['point']
        if point is None:return
        cal=point['cal'];fit=point['fit'];j=point['row'];f=cal['frequency_ghz'];scale=cal['scales'][j]
        if tabs.active==1:cal_source.data=display_data(f*1000,dict(AA=cal['pol'][j,:,0],BB=cal['pol'][j,:,1],average=cal['average'][j]),cal_plot.x_range.start,cal_plot.x_range.end)
        v=(1-f*scale/ds.cfg['rest_freq_ghz'])*C
        pol_index={'AA':0,'BB':1}.get(baseline_pol.value)
        observed=cal['average'][j] if pol_index is None else cal['pol'][j,:,pol_index]
        components={key:fit[key][j] if pol_index is None else fit[key+'_pol'][j,:,pol_index] for key in ('broad','local')}
        upper.title.text=('Mean of independent AA/BB baseline fits' if pol_index is None else baseline_pol.value+' · independently fitted baseline')
        if tabs.active==2:upper_source.data=display_data(v,dict(observed=observed,broad=components['broad'],local=np.where(np.abs(v)<=fit['local_fit_kms'],components['local'],np.nan),total=components['broad']+components['local']),upper.x_range.start,upper.x_range.end)
        if tabs.active==2:lower_source.data=display_data(v,dict(residual=fit['residual'][j],residual_AA=fit['residual_pol'][j,:,0],residual_BB=fit['residual_pol'][j,:,1]),lower.x_range.start,lower.x_range.end)
        for boxes in manual_boxes:
            windows=point.get('baseline_exclusions_kms',[])
            for k,box in enumerate(boxes):
                box.visible=k<len(windows)
                if box.visible:box.left,box.right=windows[k]
        comparison=state['comparison']
        if comparison is not None and comparison['index']==int(index.value):
            previous_source.data=display_data(comparison['velocity'],dict(y=comparison['mean']),lower.x_range.start,lower.x_range.end)
            previous_line.visible=bool(compare.active)
        else:previous_source.data=dict(x=[],y=[]);previous_line.visible=False
        ref=point['reference'];jj=int(np.argmin(abs(ref['frequency_ghz']*1000-ref_channel.value)));ii=int(np.argmin(abs(ref['mjd']-cal['mjd'][j])))
        if tabs.active==0:rt_source.data=dict(x=(ref['mjd']-cal['mjd'][j])*86400,rawAA=ref['raw'][:,0,jj],smoothAA=ref['smoothed'][:,0,jj],rawBB=ref['raw'][:,1,jj],smoothBB=ref['smoothed'][:,1,jj])
        if tabs.active==0:rf_source.data=dict(x=ref['frequency_ghz']*1000,rawAA=ref['raw'][ii,0],smoothAA=ref['smoothed'][ii,0],rawBB=ref['raw'][ii,1],smoothBB=ref['smoothed'][ii,1])
        if tabs.active==0:ref_time.title.text='Tsys-reference filter · {:.4f} MHz · {} timestamps'.format(ref['frequency_ghz'][jj]*1000,ref['info']['temperature_smoothing_rows'])
        for p in (cal_plot,upper,lower):
            if (p is cal_plot and tabs.active!=1) or (p is not cal_plot and tabs.active!=2):continue
            slots=annotation_slots[id(p)]
            for k,box in enumerate(slots):
                box.visible=k<len(cal['intervals'])
                if box.visible:
                    a,b=cal['intervals'][k]
                    if p is not cal_plot:a,b=sorted([(1-a/1000*scale/ds.cfg['rest_freq_ghz'])*C,(1-b/1000*scale/ds.cfg['rest_freq_ghz'])*C])
                    box.left=a;box.right=b
        P=float(protect.value);L=float(local.value)
        for box in protected:box.left=-P;box.right=P
        for k,box in enumerate(local_boxes):box.left,box.right=((-L,-P) if k%2==0 else (P,L))
        lc=components['local'];bound=max(.05,float(np.nanmax(abs(lc)))*1.15);upper.extra_y_ranges['correction'].start=-bound;upper.extra_y_ranges['correction'].end=bound

    annotation_slots={}
    for p in (cal_plot,upper,lower):
        annotation_slots[id(p)]=[]
        for _ in range(max([len(v) for v in ds.priors.values()]+[1])):
            box=BoxAnnotation(left=0,right=0,fill_color='#db4c45',fill_alpha=.12,visible=False);p.add_layout(box);annotation_slots[id(p)].append(box)

    def completed(future,token,kind,started):
        if state['closed'] or token!=state['generation']:return
        try:result=future.result()
        except Exception as e:status.text='<b style="color:#b52323">Preview error:</b> '+html.escape(str(e));return
        if result is None:return
        if kind=='point':
            old_point=state['point']
            if old_point is not None and old_point['index']==result['index'] and old_point['fit_key']!=result['fit_key']:
                jold=old_point['row'];cold=old_point['cal']
                state['comparison']=dict(index=int(index.value),velocity=(1-cold['frequency_ghz']*cold['scales'][jold]/ds.cfg['rest_freq_ghz'])*C,mean=old_point['fit']['residual'][jold].copy())
            state['point']=result;ctx=ds.files[result['file_id']];j=result['row'];report=result['fit']['report']
            status.text='{} {} · raw row {} · MJD {:.8f} · chunk {} spectra · {:.2f} s · RAM cache {:.0f} MiB{}'.format(
                ctx['scan'],ctx['beam'],ctx['coords']['ids'][result['chunk_id']*ds.cfg['chunk_rows']+result['chunk_row']],result['cal']['mjd'][j],result['chunk_size'],time.time()-started,ds.cache_bytes/1024**2,
                ' · <b>Broad period at search boundary</b>' if report.get('at_search_boundary') else '')
            if report.get('baseline',False):
                status.text+=' · broad periods AA {:.3f}, BB {:.3f} MHz'.format(*(report['polarizations'][p]['period_MHz'] for p in ('AA','BB')))
            redraw()
        else:
            state['region_result']=result
            region_source.data=display_data(result['velocity'],dict(mean=result['mean'],marked=result['rfi_fraction']),region_plot.x_range.start,region_plot.x_range.end)
            region_label.text='Rectangle: {} / {} spectra used · {} · {:.2f} s · aligned per spectrum, no spatial gridding'.format(result['used'],result['total'],'complete selection' if result['used']==result['total'] else 'CHUNK-SAMPLED PREVIEW; refine for all',time.time()-started)
            status.text='Rectangle preview ready; no products written.'

    def submit(kind,full=False):
        try:options=params()
        except ValueError as e:status.text=html.escape(str(e));return
        if on_settings is not None:
            try:on_settings(options,int(state['index']))
            except ValueError as e:status.text=html.escape(str(e));return
        state['generation']+=1;token=state['generation'];started=time.time()
        status.text='Computing read-only '+kind+' preview…'
        previous=state.get('future')
        if previous is not None:previous.cancel()
        if kind=='point':future=executor.submit(lambda:ds.point(state['index'],options,selected_only=True) if token==state['generation'] else None)
        else:
            ids=state['region'] if state['region'] is not None else np.array([state['index']])
            future=executor.submit(ds.region_spectrum,ids,options,0 if full else 128,lambda:state['closed'] or token!=state['generation'],
                lambda partial:doc.add_next_tick_callback(lambda:region_progress(partial,token)) if not state['closed'] else None)
        state['future']=future
        future.add_done_callback(lambda f:doc.add_next_tick_callback(lambda:completed(f,token,kind,started)) if not state['closed'] else None)

    def region_progress(partial,token):
        if token==state['generation']:
            state['region_result']=partial
            region_source.data=display_data(partial['velocity'],dict(mean=partial['mean'],marked=partial['rfi_fraction']),region_plot.x_range.start,region_plot.x_range.end)
            region_label.text='Rectangle preview: '+str(partial['used'])+' / '+str(partial['target'])+' preview spectra processed; '+str(partial['total'])+' in the full rectangle. Partial mean; no gridding or output writes.'

    def clicked(x,y):
        allowed=np.flatnonzero(active)
        if not len(allowed):return
        distance=((ds.ra[allowed]-x)*np.cos(np.deg2rad(y)))**2+(ds.dec[allowed]-y)**2
        chosen=int(allowed[np.argmin(distance)]);index.value=chosen

    def index_changed(attr,old,new):
        state['index']=max(0,min(len(ds.ra)-1,int(new)))
        marker.data=dict(x=[ds.ra[state['index']]],y=[ds.dec[state['index']]])
        submit('point')
    index.on_change('value',index_changed)

    def filters_changed(attr,old,new):
        active[:]=True
        for fi,ctx in enumerate(ds.files):
            if (scan.value!='All' and ctx['scan']!=scan.value) or (beam.value!='All' and ctx['beam']!=beam.value):active[ds.offsets[fi]:ds.offsets[fi+1]]=False
        visible.booleans=active.tolist();selection_info.text='{} / {} centres shown · full dataset still used by Run all'.format(int(active.sum()),len(active))
        clicked(ds.cfg['header']['CRVAL1'],ds.cfg['header']['CRVAL2'])
    scan.on_change('value',filters_changed);beam.on_change('value',filters_changed)

    def rectangle_changed(event):
        if not event.final or event.geometry.get('type')!='rect':return
        g=event.geometry;ids=ds.rectangle(g['x0'],g['x1'],g['y0'],g['y1']);state['region']=ids[active[ids]]
        region_label.text='Rectangle contains {} spectra; preparing sampled preview…'.format(len(state['region']));submit('region')

    pending={'handle':None}
    def parameters_changed(attr,old,new):
        state['generation']+=1;status.text='Settings changed; updating read-only preview…'
        ng.disabled=mode.value!='gaussian';nf.disabled=ng.disabled;nb.disabled=mode.value!='rebin'
        if pending['handle'] is not None:
            try:doc.remove_timeout_callback(pending['handle'])
            except ValueError:pass
        pending['handle']=doc.add_timeout_callback(lambda:submit('region' if tabs.active==3 else 'point'),300)
    exclusions.on_change('value_input',parameters_changed)
    for widget in (nt_cal,nt,protect,local,mode,nb,ng,nf,vmin,vmax):widget.on_change('value',parameters_changed)
    baseline.on_change('active',parameters_changed)
    ref_channel.on_change('value',lambda attr,old,new:redraw())
    baseline_pol.on_change('value',lambda attr,old,new:redraw())
    compare.on_change('active',lambda attr,old,new:redraw())
    redraw_timer={'handle':None}
    def queue_redraw(attr,old,new):
        if redraw_timer['handle'] is not None:
            try:doc.remove_timeout_callback(redraw_timer['handle'])
            except ValueError:pass
        redraw_timer['handle']=doc.add_timeout_callback(redraw,100)
    for p in (cal_plot,upper,lower):
        p.x_range.on_change('start',queue_redraw)
        p.x_range.on_change('end',queue_redraw)
    def region_zoom(attr,old,new):
        result=state['region_result']
        if result is not None:region_source.data=display_data(result['velocity'],dict(mean=result['mean'],marked=result['rfi_fraction']),region_plot.x_range.start,region_plot.x_range.end)
    region_plot.x_range.on_change('start',region_zoom);region_plot.x_range.on_change('end',region_zoom)
    refine.on_click(lambda:submit('region',True))
    def cancel_preview():
        state['generation']+=1
        if state.get('future') is not None:state['future'].cancel()
        status.text='Preview cancelled; existing curves retained. No products written.'
    cancel.on_click(cancel_preview)

    def alive(pid):
        try:
            pid=int(pid);cmd=(Path('/proc')/str(pid)/'cmdline').read_bytes()
            if b'fastfish.native_stages' not in cmd or str(root/'jobs').encode() not in cmd:return False
            os.kill(pid,0);return True
        except (OSError,ValueError,TypeError):return False

    def start_batch(stage,options=None,chain=False,checked=False):
        if not allow_batches:
            batch_status.text='Use STAGE run at the fish prompt to apply the current shared settings.';return
        # Product writes are confined to this explicit button path and its child.
        root.mkdir(parents=True,exist_ok=True)
        with (root/'.launch.lock').open('a+') as lock:
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            active_job=root/'active_job.json'
            if active_job.exists() and alive(read(active_job).get('pid')):
                batch_status.text='A batch job is already active. Wait or stop it before launching another.';return
            opts=params() if options is None else options
            recipe=dict(config=str(ds.config_path),selection=str(ds.selection_path),root=str(root),stage=stage,params=opts,
                        catalog_sha256=ds.catalog_sha256,source_sha256=source_hash(),validate=checked)
            key=str(time.time_ns());jobs=root/'jobs';jobs.mkdir(exist_ok=True);p=jobs/(key+'.json');dump(p,recipe)
            log=jobs/(key+'.log')
            with log.open('ab') as stream:
                process=subprocess.Popen([sys.executable,'-m','fastfish.native_stages','--recipe',str(p)],stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            dump(active_job,dict(pid=process.pid,stage=stage,log=str(log),recipe=str(p)))
            state['job']=process;state['chain']=(opts,checked) if chain else None
            for button in buttons:button.disabled=True
            batch_status.text='Started '+stage+' for all '+str(len(ds.ra))+' selected spectra. Log: '+str(log)

    def run_button(stage,label,chain=False):
        button=Button(label=label,button_type='primary',width=310,name='run_'+stage)
        button.on_click(lambda:start_batch(stage,chain=chain and bool(make_cube.active),checked=bool(validation.active)))
        buttons.append(button);return button

    stop=Button(label='Stop batch',button_type='warning',width=130,name='stop_batch')
    def stop_batch():
        if not allow_batches:
            batch_status.text='Use stop at the fish prompt.';return
        path=root/'active_job.json'
        if path.exists():
            pid=read(path).get('pid')
            if alive(pid):os.killpg(int(pid),signal.SIGTERM)
        state['chain']=None
    stop.on_click(stop_batch)

    def poll():
        process=state['job']
        if process is not None and process.poll() is not None:
            result=process.returncode;state['job']=None
            chain=state['chain'];state['chain']=None
            if result==0 and chain:
                start_batch('cube',options=chain[0],checked=chain[1]);return
            update_availability()
            if result!=0:batch_status.text='<b style="color:#b52323">Batch stopped/failed (exit '+str(result)+').</b> Inspect the job log; completed beam checkpoints remain reusable.'
        path=root/'status.json'
        if path.exists() and (state['job'] is not None or process is not None):
            report=read(path)
            if state['job'] is None and process is not None and process.returncode!=0:return
            batch_status.text='Batch '+html.escape(json.dumps(report))
        if state['job'] is None:update_availability()
    doc.add_periodic_callback(poll,2000)
    ref_button=run_button('reference','Run all · reference smoothing')
    cal_button=run_button('calibration','Run all · calibration')
    baseline_button=run_button('baseline','Run all · baseline fitting')
    sample_button=run_button('resample','Run all and proceed',True)
    def update_availability():
        if not allow_batches:
            for button in buttons:button.disabled=True
            stop.disabled=True
            batch_status.text='FISH command session: plot edits sync to the prompt. Use STAGE run there to process data.';return
        if state['job'] is not None:return
        try:options=params()
        except ValueError:
            for button in buttons:button.disabled=True
            return
        dependencies={cal_button:'reference',baseline_button:'calibration',sample_button:'baseline'}
        ref_button.disabled=False
        for button,dependency in dependencies.items():button.disabled=not (stage_path(root,dependency,ds,options)/'summary.json').exists()
    tabs=Tabs(tabs=[
        Panel(title='1 · Reference',child=column(split(column(nt_cal,nt,ref_channel,Div(text='Raw = Tsys calculated with the selected ON−OFF window. Smooth = Tsys-reference application filter, before deriving SCI/Tsys gain. ON−OFF control requires extracted CAL; precomputed Tsys embeds that step. Display is read-only.',sizing_mode='stretch_width'),sizing_mode='stretch_width'),column(ref_time,ref_freq,sizing_mode='stretch_width')),ref_button,sizing_mode='stretch_width')),
        Panel(title='2 · Calibration',child=column(split(pointing_map('Click a pointing centre'),cal_plot),cal_button,sizing_mode='stretch_width')),
        Panel(title='3 · Baseline',child=column(row(previous,next_,refit,undo,clear),row(protect,local),row(baseline,baseline_pol),row(draw_mode,exclusions),compare,workspace_help,split(pointing_map('Pointing selection'),column(upper,lower,sizing_mode='stretch_width')),Div(text='AA and BB are fitted independently before averaging. Lower panel: corrected AA, BB and mean. Grey: protected H I. Red: RFI fit exclusions. Pale green: local-fit window. Local correction uses the right axis and applies only inside the local-fit window. Wheel/box zoom to inspect the initially clipped H I peak.',sizing_mode='stretch_width'),baseline_button,sizing_mode='stretch_width')),
        Panel(title='4 · Resample / cube',child=column(row(mode,nb),row(ng,nf),row(vmin,vmax),split(pointing_map('Draw rectangle using box-select',True),region_plot),region_label,row(refine,cancel),make_cube,validation,sample_button,sizing_mode='stretch_width'))],name='processing_tabs',sizing_mode='stretch_width')
    def arrange(attr,old,new):
        stacked=layout_choice.value=='Stacked'
        for wrap,left,right in splits:
            if stacked:
                left.sizing_mode='stretch_width';right.sizing_mode='stretch_width'
                wrap.children=[column(left,right,sizing_mode='stretch_width')]
            else:
                left.sizing_mode='fixed';left.width=440
                right.sizing_mode='stretch_width'
                wrap.children=[row(left,right,sizing_mode='stretch_width')]
        for p in maps:p.height=350 if stacked else 530
    attach_workspace(upper,lower,draw_mode,exclusions,protect,local,index,active,previous,next_,refit,undo,clear,workspace_help,submit,state)
    layout_choice.on_change('value',arrange);arrange(None,None,None);update_availability()
    tabs.active={'reference':0,'calibration':1,'baseline':2,'resample':3}.get(initial_view,0)
    tabs.on_change('active',lambda attr,old,new:submit('region' if new==3 else 'point'))
    doc.add_root(column(title,note,controls,status,tabs,stop,batch_status,sizing_mode='stretch_width',max_width=1420))
    doc.title='FISH · '+str(ds.cfg.get('name','observation'))
    def close(session):state.update(closed=True,generation=state['generation']+1);executor.shutdown(wait=False)
    doc.on_session_destroyed(close)
    doc.add_next_tick_callback(lambda:submit('point'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config');p.add_argument('--selection');p.add_argument('--run-root');p.add_argument('--session-state');p.add_argument('--port',type=int,default=8765)
    args=p.parse_args()
    from bokeh.application import Application
    from bokeh.application.handlers.function import FunctionHandler
    from bokeh.server.server import Server
    def build(doc):
        if args.session_state:
            from fastfish.processing.session import read_session,write_session
            snapshot=read_session(args.session_state);revision=[snapshot['revision']]
            from fastfish.commands.shell import validate_params
            original_cfg=json.loads(Path(snapshot['config']).read_text())
            def sync(options,index):
                options=validate_params(original_cfg,options)
                value=write_session(args.session_state,dict(params=options,index=index),expected_revision=revision[0]);revision[0]=value['revision']
            build_document(doc,snapshot['config'],snapshot['selection'],snapshot['root'],snapshot['params'],snapshot['index'],sync,False)
        else:build_document(doc,args.config,args.selection,args.run_root)
    if not args.session_state and not all((args.config,args.selection,args.run_root)):p.error('Supply --config, --selection and --run-root, or --session-state')
    app=Application(FunctionHandler(build))
    server=Server({'/':app},address='127.0.0.1',port=args.port,allow_websocket_origin=['localhost:'+str(args.port),'127.0.0.1:'+str(args.port)])
    server.start();print('Interactive trial at http://127.0.0.1:'+str(args.port),flush=True);server.io_loop.start()


if __name__=='__main__':main()
