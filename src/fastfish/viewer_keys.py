"""Read-only CLASS/GREG GO VIEW navigation for the native plotting canvas."""
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore,QtWidgets

HELP='''CLASS/GREG GO VIEW navigation (supported subset)

Put the mouse over a plot, then use:
Z — drag a rectangle to zoom. On spectra, only the horizontal range changes.
Shift+Z — drag a rectangle to change BOTH horizontal and vertical ranges.
B — previous view of the panel under the cursor.
F — full footprint or spectrum. Outside panels: reset both.
R / P — create an editable rectangle / polygon on the footprint.
V — report cursor coordinates in the status bar; no spectrum is fetched.
H — this help.
E / Q — close this viewer; terminal and reduction jobs continue.
Escape — cancel rectangle zoom.

Keys apply only while the plotting canvas has keyboard focus.
Click the canvas to focus it. Typing in parameter/text controls is unaffected.
Free dragging and wheel navigation are disabled; use Z, B and F.
Escape cancels an armed zoom even when another control has focus.
Cube movie, slices, pixel editing and export keys are not implemented.
Reference: IRAM GREG tutorial, Dealing with LMV cubes (2019).
'''

class ClassViewBox(pg.ViewBox):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.class_zoom=False;self.keep_y=False;self.original_y=None;self.history=[];self.setMouseEnabled(x=False,y=False)
    def remember_view(self):
        self.history.append([list(axis) for axis in self.viewRange()])
        self.history=self.history[-50:]
    def previous_view(self):
        self.cancel_zoom()
        if not self.history:return False
        limits=self.history.pop()
        self.setRange(xRange=limits[0],yRange=limits[1],padding=0)
        return True
    def arm_zoom(self,keep_y=False):
        self.class_zoom=True;self.keep_y=keep_y;self.original_y=None;self.setMouseEnabled(x=True,y=True);self.setMouseMode(self.RectMode)
    def cancel_zoom(self):
        self.class_zoom=False;self.setMouseEnabled(x=False,y=False);self.rbScaleBox.hide()
        if hasattr(self,"pending_view"):del self.pending_view
    def finish_zoom(self):
        if hasattr(self,"pending_view"):
            self.history.append(self.pending_view);self.history=self.history[-50:];del self.pending_view
        if self.keep_y and self.original_y is not None:self.setYRange(*self.original_y,padding=0)
        self.cancel_zoom()
    def mouseDragEvent(self,event,axis=None):
        active=self.class_zoom and event.button()==QtCore.Qt.LeftButton
        if not active:event.accept();return
        if active and event.isStart():
            self.original_y=list(self.viewRange()[1]);self.pending_view=[list(a) for a in self.viewRange()]
        super().mouseDragEvent(event,axis=axis)
        if active and event.isFinish():self.finish_zoom()

    def wheelEvent(self,event,axis=None):
        event.accept()

class ClassKeys(QtCore.QObject):
    def __init__(self,viewer):
        super().__init__(viewer);self.viewer=viewer;self.canvas=viewer.graph;self.viewport=viewer.graph.viewport();self.position=None;self.help_dialog=None
        viewer.graph.setFocusPolicy(QtCore.Qt.StrongFocus)
        QtWidgets.QApplication.instance().installEventFilter(self)
        viewer.graph.scene().sigMouseMoved.connect(self.moved)
    def moved(self,position):self.position=position
    def target(self):
        if self.position is not None:
            for plot in (self.viewer.map,self.viewer.top,self.viewer.bottom):
                if plot.isVisible() and plot.vb.scene() is not None and plot.vb.sceneBoundingRect().contains(self.position):return plot
        return None
    def eventFilter(self,obj,event):
        if event.type() in (QtCore.QEvent.ShortcutOverride,QtCore.QEvent.KeyPress) and event.key()==QtCore.Qt.Key_Escape:
            armed=any(p.vb.class_zoom for p in (self.viewer.map,self.viewer.top,self.viewer.bottom))
            inside=isinstance(obj,QtWidgets.QWidget) and (obj is self.viewer or self.viewer.isAncestorOf(obj))
            if armed and inside:
                if event.type()==QtCore.QEvent.KeyPress:self.handle('ESC')
                event.accept();return True
        if obj is not self.canvas and obj is not self.viewport:return False
        if event.type()==QtCore.QEvent.Leave:self.position=None
        if event.type()!=QtCore.QEvent.KeyPress:return False
        if event.modifiers() & (QtCore.Qt.ControlModifier|QtCore.Qt.AltModifier|QtCore.Qt.MetaModifier):return False
        key=event.text().upper()
        if event.key()==QtCore.Qt.Key_Escape:key='ESC'
        if key=='Z' and event.modifiers() & QtCore.Qt.ShiftModifier:key='ZXY'
        if key not in ('Z','ZXY','B','F','R','P','V','H','E','Q','ESC'):return False
        self.handle(key);event.accept();return True
    def handle(self,key):
        v=self.viewer;plot=self.target()
        if key in ('R','P'):
            if plot is not v.map:v.statusBar().showMessage('Move the cursor over the footprint, then press R or P');return
            v.select_region('rectangle' if key=='R' else 'polygon');return
        if key=='H':
            if self.help_dialog is None:
                self.help_dialog=QtWidgets.QMessageBox(v);self.help_dialog.setWindowTitle('FISH · CLASS-style keys');self.help_dialog.setText(HELP)
            self.help_dialog.show();self.help_dialog.raise_();return
        if key in ('E','Q'):v.close();return
        if key=='ESC':
            for p in (v.map,v.top,v.bottom):p.vb.cancel_zoom()
            v.statusBar().showMessage('Zoom selection cancelled');return
        if key in ('Z','ZXY'):
            if plot is None:v.statusBar().showMessage('Move the cursor over a map or spectrum, then press Z');return
            for p in (v.map,v.top,v.bottom):p.vb.cancel_zoom()
            plot.vb.arm_zoom(keep_y=plot is not v.map and key!='ZXY');v.statusBar().showMessage('Drag a rectangle to zoom BOTH axes; Escape cancels.' if key=='ZXY' else 'Drag a zoom rectangle; Escape cancels. Spectrum zoom keeps the intensity range. Shift+Z changes both axes.');return
        if key=='B':
            self.handle('ESC')
            restored=False
            for p in (v.map,v.top,v.bottom):
                if plot is None or plot is p:restored=p.vb.previous_view() or restored
            v.statusBar().showMessage('Previous view restored' if restored else 'No previous view')
            return
        if key=='F':
            self.handle('ESC')
            for p in (v.map,v.top,v.bottom):
                if plot is None or plot is p:p.vb.remember_view()
            if plot is None or plot is v.map:v.full_footprint()
            for p in (v.top,v.bottom):
                if plot is None or plot is p:p.enableAutoRange();p.autoRange()
            v.statusBar().showMessage('Full extent restored');return
        if key=='V':
            if plot is None:v.statusBar().showMessage('Move the cursor over a plot, then press V');return
            p=plot.vb.mapSceneToView(self.position)
            if plot is v.map:text='Cursor RA={:.8f} deg  Dec={:.8f} deg'.format(p.x(),p.y())
            else:text='Cursor {}={:.8g} {}  y={:.8g} K'.format(plot.getAxis('bottom').labelText,p.x(),plot.getAxis('bottom').labelUnits,p.y())
            v.statusBar().showMessage(text)
