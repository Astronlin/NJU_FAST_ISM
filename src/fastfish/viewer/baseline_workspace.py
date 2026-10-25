"""Session-only baseline window editing and filtered spectrum navigation."""
from copy import deepcopy
from bokeh.events import SelectionGeometry
from bokeh.models import BoxSelectTool
from fastfish.baseline.windows import parse_windows,format_windows

def attach_workspace(upper,lower,draw_mode,exclusions,protect,local,index,active,
                     previous,next_,refit,undo,clear,help_,submit,state):
    controls=(protect,local,exclusions)
    history=[];guard={'restoring':False}
    def snapshot():return dict(protect=protect.value,local=local.value,exclusions=exclusions.value_input)
    saved=snapshot()
    def changed(attr,old,new):
        nonlocal saved
        now=snapshot()
        if not guard['restoring'] and now!=saved:
            history.append(deepcopy(saved));del history[:-30]
        saved=now;undo.disabled=not history
    for widget in controls:widget.on_change('value_input' if widget is exclusions else 'value',changed)
    def undo_edit():
        nonlocal saved
        if not history:return
        target=history.pop();guard['restoring']=True
        try:
            protect.value=target['protect'];local.value=target['local'];exclusions.value_input=target['exclusions']
        finally:guard['restoring']=False;saved=snapshot();undo.disabled=not history
    undo.on_click(undo_edit)
    clear.on_click(lambda:setattr(exclusions,'value_input',''))
    refit.on_click(lambda:submit('point'))
    def navigate(direction):
        import numpy as np
        ids=np.flatnonzero(active);ids=ids[ids>int(index.value)] if direction>0 else ids[ids<int(index.value)]
        if len(ids):index.value=int(ids[0] if direction>0 else ids[-1])
        else:help_.text='No further spectrum in the current scan/beam filter.'
    previous.on_click(lambda:navigate(-1));next_.on_click(lambda:navigate(1))
    tools=[]
    def selected(event):
        if not event.final or event.geometry.get('type')!='rect' or draw_mode.value=='Inspect':return
        g=event.geometry;lo,hi=sorted((float(g['x0']),float(g['x1'])))
        if hi<=lo:return
        try:
            if draw_mode.value=='Add exclusion':exclusions.value_input=format_windows(parse_windows(exclusions.value_input)+[[lo,hi]])
            else:
                width=max(abs(lo),abs(hi));widget=protect if draw_mode.value=='Set H I protection' else local
                if not widget.low<=width<=widget.high:raise ValueError('Draw half-width between {} and {} km/s'.format(widget.low,widget.high))
                widget.value=width
            help_.text='Window edit is session-only. Run all explicitly applies the current windows to the entire dataset.'
        except ValueError as e:help_.text=str(e)
    for plot in (upper,lower):
        tool=BoxSelectTool(dimensions='width');plot.add_tools(tool);plot.on_event(SelectionGeometry,selected);tools.append((plot,tool))
    def mode_changed(attr,old,new):
        for plot,tool in tools:plot.toolbar.active_drag=tool if new!='Inspect' else 'auto'
    draw_mode.on_change('value',mode_changed)
    return history
