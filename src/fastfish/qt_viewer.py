"""Native local PyQtGraph viewer backed by a persistent local/SSH worker."""
import argparse,base64,colorsys,hashlib,io,json,os,re,shlex,subprocess,sys
from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore,QtWidgets,QtGui
from .viewer_keys import ClassViewBox,ClassKeys
from .baseline_components import normalize_components,default_components
Signal=QtCore.pyqtSignal if hasattr(QtCore,'pyqtSignal') else QtCore.Signal

def scan_color(name,dark=False):
    """Stable name-based hues, including when only one scan is loaded."""
    match=re.match(r'^(.*?)(\d+)$',str(name))
    prefix=match.group(1) if match else str(name)
    number=int(match.group(2)) if match else 0
    seed=int.from_bytes(hashlib.sha256(prefix.encode()).digest()[:4],'big')/2**32
    hue=(seed+number*.6180339887498949)%1
    return tuple(int(round(255*x)) for x in colorsys.hls_to_rgb(hue,.64 if dark else .42,.72))

class WrapLayout(QtWidgets.QLayout):
    """Controls wrap to available width instead of forcing a wide window."""
    def __init__(self):
        super().__init__();self.items=[];self.setContentsMargins(0,0,0,0);self.setSpacing(6)
    def addItem(self,item):self.items.append(item)
    def count(self):return len(self.items)
    def itemAt(self,index):return self.items[index] if 0<=index<len(self.items) else None
    def takeAt(self,index):return self.items.pop(index) if 0<=index<len(self.items) else None
    def expandingDirections(self):return QtCore.Qt.Orientations(0)
    def hasHeightForWidth(self):return True
    def heightForWidth(self,width):return self.arrange(QtCore.QRect(0,0,width,0),True)
    def minimumSize(self):
        result=QtCore.QSize()
        for item in self.items:result=result.expandedTo(item.minimumSize())
        return result
    def sizeHint(self):return self.minimumSize()
    def setGeometry(self,rect):
        super().setGeometry(rect);self.arrange(rect,False)
    def arrange(self,rect,measure):
        x=rect.x();y=rect.y();row=0;gap=self.spacing()
        for item in self.items:
            size=item.sizeHint();width=min(size.width(),max(1,rect.width()))
            if x>rect.x() and x+width>rect.x()+rect.width():
                x=rect.x();y+=row+gap;row=0
            if not measure:item.setGeometry(QtCore.QRect(QtCore.QPoint(x,y),QtCore.QSize(width,size.height())))
            x+=width+gap;row=max(row,size.height())
        return y+row-rect.y()

class ImageLoader(QtCore.QThread):
    loaded=Signal(object);failed=Signal(str)
    def __init__(self,path,parent):super().__init__(parent);self.path=path
    def run(self):
        try:
            from .sky_image import load_sky_image
            self.loaded.emit(load_sky_image(self.path))
        except Exception as exc:self.failed.emit(str(exc))

class Worker(QtCore.QThread):
    ready=Signal(object);failed=Signal(str)
    def __init__(self,command):
        super().__init__();self.command=command;self.request=None;self.proc=None
    def run(self):
        try:
            if self.proc is None:self.proc=subprocess.Popen(self.command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
            request=dict(self.request);request.setdefault('compress',os.path.basename(self.command[0])=='ssh')
            self.proc.stdin.write(json.dumps(request)+'\n');self.proc.stdin.flush();line=self.proc.stdout.readline()
            if not line:raise RuntimeError('Worker disconnected; check SSH authentication and server paths')
            r=json.loads(line)
            if 'error' in r:raise RuntimeError(r['error'])
            with np.load(io.BytesIO(base64.b64decode(r['arrays'])),allow_pickle=False) as z:r['data']={k:z[k] for k in z.files}
            self.ready.emit(r)
        except Exception as e:self.failed.emit(str(e))
    def shutdown(self):
        if self.proc is not None:self.proc.terminate()
        self.wait(3000)

class SettingsWindow(QtWidgets.QDialog):
    def __init__(self,viewer):
        super().__init__(viewer,QtCore.Qt.Tool);self.viewer=viewer;self.setWindowTitle('FISH · parameter console');self.setMinimumSize(650,180)
        layout=QtWidgets.QVBoxLayout(self);layout.setContentsMargins(6,6,6,6);bar=QtWidgets.QHBoxLayout();layout.addLayout(bar)
        apply_button=QtWidgets.QPushButton('Apply / refit');apply_button.clicked.connect(viewer.refit);bar.addWidget(apply_button)
        attach_button=QtWidgets.QPushButton('Attach to viewer');attach_button.clicked.connect(viewer.attach_settings);bar.addWidget(attach_button);bar.addStretch(1)
    def reject(self):
        self.viewer.attach_settings();super().reject()
    def closeEvent(self,event):
        self.viewer.attach_settings();event.accept()

class Viewer(QtWidgets.QMainWindow):
    def __init__(self,command,view='baseline'):
        super().__init__();self.setWindowTitle('FISH · native spectral workspace');self.resize(1100,750);self.state=None;self.data=None;self.masks=[];self.last_stage=None;self.footprint_layout=False;self.region_preview_data=None;self.render_cache=OrderedDict();self.render_cache_bytes=0;self.grid_orders={};self.mask_signature=None;self.pending_request=None;self.data_view=None;self.syncing_spectrum_x=False;self.settings_detached=False;self.settings_window=None;self.closing=False;self.synced_params=None
        self.scan_curves={};self.scan_names=[];self.scan_groups=None;self.scan_legend=None
        self.worker=Worker(command);self.worker.ready.connect(self.received);self.worker.failed.connect(self.error);self.worker.finished.connect(self.dispatch_pending)
        root=QtWidgets.QWidget();self.setCentralWidget(root);outer=QtWidgets.QVBoxLayout(root);outer.setContentsMargins(6,6,6,6);outer.setSpacing(4);self.outer_layout=outer
        toolbar=QtWidgets.QWidget();toolbar_layout=QtWidgets.QVBoxLayout(toolbar);toolbar_layout.setContentsMargins(0,0,0,0);toolbar_layout.setSpacing(3);outer.addWidget(toolbar)
        bar=WrapLayout();toolbar_layout.addLayout(bar);bar.addWidget(QtWidgets.QLabel('View'))
        self.stage=QtWidgets.QComboBox();self.stage.addItems(['Footprint','Baseline','Apply calibration','Tsys calibration','Table']);self.stage.setMinimumWidth(170);self.stage.currentTextChanged.connect(self.stage_changed);bar.addWidget(self.stage)
        for name,fn in [('Previous',lambda:self.navigate(-1)),('Next',lambda:self.navigate(1)),('Refresh',lambda:self.fetch('point' if hasattr(self,'ra') and self.stage.currentText()!='Footprint' else 'map')),('Full footprint',self.full_footprint),('Baseline zoom',self.baseline_zoom),('Apply / refit',self.refit)]:
            button=QtWidgets.QPushButton(name);button.clicked.connect(fn);bar.addWidget(button)
            if name=='Baseline zoom':self.baseline_zoom_button=button
            elif name=='Apply / refit':self.refit_button=button
        self.pol=QtWidgets.QComboBox();self.pol.addItems(['AA','BB','Mean']);self.pol.currentTextChanged.connect(lambda _:self.draw());bar.addWidget(self.pol)
        self.settings_button=QtWidgets.QToolButton();self.settings_button.setText('Settings');self.settings_button.setCheckable(True);self.settings_button.setChecked(True);self.settings_button.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon);self.settings_button.toggled.connect(self.update_settings_layout);bar.addWidget(self.settings_button)
        settings_menu=QtWidgets.QMenu(self.settings_button);self.detach_settings_action=settings_menu.addAction('Detach parameter console');self.detach_settings_action.triggered.connect(self.toggle_settings_window);self.settings_button.setMenu(settings_menu);self.settings_button.setPopupMode(QtWidgets.QToolButton.MenuButtonPopup);self.settings_button.setToolTip('Show/hide settings. Use the arrow menu to detach the parameter console.')
        bar=WrapLayout();toolbar_layout.addLayout(bar)
        bar.addWidget(QtWidgets.QLabel('Style'));self.spectrum_style=QtWidgets.QComboBox();self.spectrum_style.addItems(['Centered steps','Lines']);self.spectrum_style.currentTextChanged.connect(lambda _:self.draw());bar.addWidget(self.spectrum_style)
        bar.addWidget(QtWidgets.QLabel('Display bin (ch)'));self.display_bin=QtWidgets.QSpinBox();self.display_bin.setRange(1,512);self.display_bin.setValue(32);self.display_bin.valueChanged.connect(lambda _:self.draw());bar.addWidget(self.display_bin)
        image_button=QtWidgets.QPushButton('Load 2D FITS');image_button.clicked.connect(self.load_image);bar.addWidget(image_button)
        self.image_visible=QtWidgets.QCheckBox('Show image');self.image_visible.setChecked(True);self.image_visible.toggled.connect(lambda on:self.sky_image.setVisible(on));bar.addWidget(self.image_visible)
        self.theme=QtWidgets.QComboBox();self.theme.addItems(['White / black','Black / white']);bar.addWidget(self.theme)
        self.scan_colors_button=QtWidgets.QToolButton();self.scan_colors_button.setText('Scan colors');self.scan_colors_button.setCheckable(True);self.scan_colors_button.setChecked(True);self.scan_colors_button.setToolTip('Color footprint points by scan; turn off for a single color. Open the arrow menu for the legend.');self.scan_colors_button.toggled.connect(self.update_footprint);bar.addWidget(self.scan_colors_button)
        color_menu=QtWidgets.QMenu(self.scan_colors_button);self.scan_legend_action=color_menu.addAction('Show scan legend');self.scan_legend_action.triggered.connect(self.show_scan_legend);self.scan_legend_action.setEnabled(False);self.scan_colors_button.setMenu(color_menu);self.scan_colors_button.setPopupMode(QtWidgets.QToolButton.MenuButtonPopup)
        self.region_button=QtWidgets.QPushButton('Region mean (≤128)');bar.addWidget(self.region_button);self.region_button.clicked.connect(self.region_mean)
        panel=QtWidgets.QWidget();layout=QtWidgets.QVBoxLayout(panel);layout.setContentsMargins(4,2,4,2);layout.setSpacing(3)
        scroll=QtWidgets.QScrollArea();scroll.setWidgetResizable(True);scroll.setWidget(panel);scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);outer.addWidget(scroll)
        self.control_scroll=scroll;self.settings_panel=panel;self.controls={}
        self.calibration_settings=QtWidgets.QWidget();calibration_bar=WrapLayout();self.calibration_settings.setLayout(calibration_bar);layout.addWidget(self.calibration_settings)
        self.baseline_settings=QtWidgets.QWidget();baseline_layout=QtWidgets.QVBoxLayout(self.baseline_settings);baseline_layout.setContentsMargins(0,0,0,0);baseline_layout.setSpacing(3);layout.addWidget(self.baseline_settings)
        baseline_bar=QtWidgets.QHBoxLayout();baseline_layout.addLayout(baseline_bar)
        self.resample_settings=QtWidgets.QWidget();resample_bar=WrapLayout();self.resample_settings.setLayout(resample_bar);layout.addWidget(self.resample_settings)
        for key,label,maximum in [('n_t_cal','ON−OFF smoothing',3000),('n_t_apply','Tpsr smoothing',3000),('protect_kms','Protect ±km/s',10000),('n_b','Rebin n_b',1000)]:
            group=QtWidgets.QWidget();pair=QtWidgets.QHBoxLayout(group);pair.setContentsMargins(0,0,0,0);pair.addWidget(QtWidgets.QLabel(label));w=QtWidgets.QDoubleSpinBox() if key=='protect_kms' else QtWidgets.QSpinBox();w.setRange(0 if key=='protect_kms' else 1,maximum);self.controls[key]=w;pair.addWidget(w)
            (calibration_bar if key in ('n_t_cal','n_t_apply') else baseline_bar if key=='protect_kms' else resample_bar).addWidget(group)
        baseline_bar.addWidget(QtWidgets.QLabel('Exclude km/s'))
        self.exclusions=QtWidgets.QLineEdit();self.exclusions.setPlaceholderText('[[700,900]]');self.exclusions.setMinimumWidth(150);baseline_bar.addWidget(self.exclusions,1)
        self.mode=QtWidgets.QComboBox();self.mode.addItems(['rebin','gaussian']);resample_bar.addWidget(QtWidgets.QLabel('Resampling'));resample_bar.addWidget(self.mode);self.ng=QtWidgets.QDoubleSpinBox();self.ng.setRange(1,1000);self.ng.setValue(2);resample_bar.addWidget(QtWidgets.QLabel('Gaussian n_g'));resample_bar.addWidget(self.ng)
        self.component_table=QtWidgets.QTableWidget(0,7);self.component_table.setHorizontalHeaderLabels(['Function','n','Fit & apply range (km/s)','Period MHz','↑','↓','×'])
        self.component_table.horizontalHeader().setSectionResizeMode(2,QtWidgets.QHeaderView.Stretch)
        for column,width in [(0,125),(1,50),(3,140),(4,30),(5,30),(6,30)]:self.component_table.setColumnWidth(column,width)
        self.component_table.setMinimumHeight(100);self.component_table.setMaximumHeight(180);baseline_layout.addWidget(self.component_table)
        self.add_baseline_button=QtWidgets.QPushButton('+ Add component');self.add_baseline_button.clicked.connect(lambda:self.add_component());baseline_layout.addWidget(self.add_baseline_button,0,QtCore.Qt.AlignLeft)
        self.component_table.setToolTip('Fit in list order; H I/RFI excluded. Range: full or min:max km/s. Period: auto, value, or min:max MHz.')
        self.mode.currentTextChanged.connect(self.update_resample_controls)
        self.component_curves=[]
        self.graph=pg.GraphicsLayoutWidget();self.graph.setMinimumSize(300,200);outer.addWidget(self.graph,1)
        self.map=self.graph.addPlot(row=0,col=0,rowspan=2,viewBox=ClassViewBox(),title='Pointing footprint · click to select');self.map.setLabel('bottom','RA',units='deg');self.map.setLabel('left','Dec',units='deg');self.map.invertX(True)
        self.sky_image=pg.ImageItem(axisOrder='row-major');self.sky_image.setZValue(-100);self.map.addItem(self.sky_image,ignoreBounds=True);self.image_loader=None
        self.footprint_curve=self.map.plot(pen=None,symbol='o',symbolSize=2,symbolBrush=(100,170,220,180),symbolPen=None)
        self.lod_timer=QtCore.QTimer(self);self.lod_timer.setSingleShot(True);self.lod_timer.setInterval(80);self.lod_timer.timeout.connect(self.update_footprint)
        self.map.vb.sigRangeChanged.connect(lambda *args:self.lod_timer.start());self.map.vb.sigResized.connect(lambda *args:self.lod_timer.start())
        self.top=self.graph.addPlot(row=0,col=1,viewBox=ClassViewBox());self.bottom=self.graph.addPlot(row=1,col=1,viewBox=ClassViewBox());self.top.vb.sigXRangeChanged.connect(lambda view,limits:self.sync_spectrum_x(self.bottom,limits));self.bottom.vb.sigXRangeChanged.connect(lambda view,limits:self.sync_spectrum_x(self.top,limits));self.graph.ci.layout.setColumnStretchFactor(0,35);self.graph.ci.layout.setColumnStretchFactor(1,65);self.top.addLegend();self.bottom.addLegend()
        self.top.getAxis('bottom').enableAutoSIPrefix(False);self.bottom.getAxis('bottom').enableAutoSIPrefix(False);self.top.getAxis('left').enableAutoSIPrefix(False);self.bottom.getAxis('left').enableAutoSIPrefix(False)
        self.top.setLabel('left','Temperature',units='K');self.bottom.setLabel('left','Temperature',units='K');self.top.getAxis('left').setWidth(70);self.bottom.getAxis('left').setWidth(70)
        self.curves=[self.top.plot(pen=pg.mkPen(c,width=1),name=n,autoDownsample=True,downsampleMethod='peak',clipToView=True) for c,n in [('#a4acb6','Observed'),('#f471a0','Total model')]]
        self.corrected=self.bottom.plot(pen='#79baff',name='Corrected',autoDownsample=True,downsampleMethod='peak',clipToView=True)
        self.protected=pg.LinearRegionItem([-250,250],brush=(130,140,155,35));self.top.addItem(self.protected)
        self.protected.sigRegionChangeFinished.connect(lambda:self.dragged('protect_kms',self.protected))
        self.region_shape=None;self.polygon_roi=None;self.roi=pg.RectROI([0,0],[0.1,0.1],pen='#f4a340');self.map.addItem(self.roi);self.roi.hide()
        self.marker=self.map.plot(pen=None,symbol='o',symbolSize=12,symbolBrush='#ffba55');self.map.scene().sigMouseClicked.connect(self.clicked)
        self.keyboard=ClassKeys(self)
        self.stage.setCurrentText({'where':'Footprint','reference':'Tsys calibration','calibration':'Apply calibration','resample':'Table','baseline':'Baseline'}.get(view,'Baseline'));self.stage_changed();self.theme.currentIndexChanged.connect(self.apply_theme);self.apply_theme();self.statusBar().showMessage('Connecting to computation worker…');self.fetch('map')
    def add_component(self,spec=None):
        spec=spec or dict(function='polynomial',degree=1,range_kms=[-1500.,1500.])
        table=self.component_table;row=table.rowCount()
        if row>=16:self.error('At most 16 baseline components');return
        table.insertRow(row);table.setFixedHeight(max(100,min(180,32+30*table.rowCount())));function=QtWidgets.QComboBox();function.addItems(['sin/cos','polynomial']);function.setCurrentText('sin/cos' if spec['function']=='sincos' else 'polynomial');table.setCellWidget(row,0,function)
        degree=QtWidgets.QSpinBox();degree.setRange(0,20);degree.setValue(spec.get('degree',1));table.setCellWidget(row,1,degree)
        window=QtWidgets.QLineEdit('full' if spec.get('range_kms') is None else ':'.join(str(x) for x in spec['range_kms']));table.setCellWidget(row,2,window)
        text=str(spec['period_MHz']) if spec.get('period_MHz') is not None else ':'.join(map(str,spec['period_range_MHz'])) if spec.get('period_range_MHz') else 'auto'
        period=QtWidgets.QLineEdit(text);table.setCellWidget(row,3,period)
        def changed(*args):degree.setEnabled(function.currentText()=='polynomial');period.setEnabled(function.currentText()=='sin/cos')
        function.currentTextChanged.connect(changed);changed()
        for col,label,delta in [(4,'↑',-1),(5,'↓',1),(6,'×',0)]:
            button=QtWidgets.QPushButton(label);table.setCellWidget(row,col,button)
            button.clicked.connect(lambda checked=False,b=button,d=delta:self.edit_component(b,d))
        self.update_settings_layout()
    def read_components(self):
        result=[];table=self.component_table
        for row in range(table.rowCount()):
            kind='sincos' if table.cellWidget(row,0).currentText()=='sin/cos' else 'polynomial';text=table.cellWidget(row,2).text().strip()
            spec=dict(function=kind,range_kms=None if text.lower()=='full' else list(map(float,text.split(':'))))
            if kind=='polynomial':spec['degree']=table.cellWidget(row,1).value()
            else:
                text=table.cellWidget(row,3).text().strip().lower()
                if ':' in text:spec['period_MHz']=None;spec['period_range_MHz']=list(map(float,text.split(':')))
                else:spec['period_MHz']=None if text=='auto' else float(text)
            result.append(spec)
        return normalize_components(result)
    def populate_components(self,specs):
        self.component_table.setRowCount(0);self.component_table.setFixedHeight(100)
        for spec in specs:self.add_component(spec)
        self.update_settings_layout()
    def edit_component(self,button,delta):
        row=self.component_table.indexAt(button.pos()).row()
        if row<0:return
        try:specs=self.read_components()
        except (ValueError,TypeError) as exc:self.error(str(exc));return
        if delta==0:specs.pop(row)
        elif 0<=row+delta<len(specs):specs[row],specs[row+delta]=specs[row+delta],specs[row]
        self.populate_components(specs)
    def load_image(self):
        if self.image_loader is not None and self.image_loader.isRunning():return
        path,_=QtWidgets.QFileDialog.getOpenFileName(self,'Load a 2D celestial FITS image','','FITS images (*.fits *.fit *.fts);;All files (*)')
        if not path:return
        self.image_loader=ImageLoader(path,self);self.image_loader.loaded.connect(self.image_loaded);self.image_loader.failed.connect(self.error)
        self.statusBar().showMessage('Loading and aligning reference image…');self.image_loader.start()
    def image_loaded(self,result):
        self.sky_image.setImage(result['image'],autoLevels=False,levels=result['levels'])
        self.sky_image.setRect(QtCore.QRectF(*result['rect']));self.sky_image.setOpacity(.65)
        self.image_visible.setChecked(True);self.sky_image.show()
        self.statusBar().showMessage('Reference image aligned to ICRS; display only, 1–99 percentile grayscale')
    def fit_screen(self):
        screen=self.screen() if hasattr(self,'screen') else QtWidgets.QApplication.primaryScreen()
        if screen is None:return
        available=screen.availableGeometry()
        width=max(400,min(1380,available.width()-40))
        height=max(360,min(850,available.height()-60))
        self.resize(width,height)
        self.move(available.x()+20,available.y()+20)
    def apply_theme(self):
        dark=self.theme.currentIndex()==1;ink='#ffffff' if dark else '#000000';background='#000000' if dark else '#ffffff'
        self.graph.setBackground(background)
        for plot in (self.map,self.top,self.bottom):
            for name in ('left','bottom','right','top'):
                axis=plot.getAxis(name);axis.setPen(ink);axis.setTextPen(ink);axis.label.setDefaultTextColor(pg.mkColor(ink))
            plot.titleLabel.setText(plot.titleLabel.text,color=ink)
            if plot.legend is not None:plot.legend.setLabelTextColor(ink)
        self.footprint_curve.setSymbolBrush(ink)
        for group,curve in self.scan_curves.items():curve.setSymbolBrush(scan_color(self.scan_names[group],dark))
        if self.scan_legend is not None and self.scan_legend.isVisible():self.populate_scan_legend()
        self.curves[0].setPen(pg.mkPen(ink,width=1))
        self.corrected.setPen('#79baff' if dark else '#005bb5')
    def plot_title(self,plot,text):
        plot.setTitle(text,color='#ffffff' if self.theme.currentIndex()==1 else '#000000')
    def toggle_settings_window(self):
        if self.settings_detached:self.attach_settings();return
        if self.settings_window is None:self.settings_window=SettingsWindow(self)
        self.settings_detached=True;self.outer_layout.removeWidget(self.control_scroll);self.settings_window.layout().addWidget(self.control_scroll,1)
        self.control_scroll.setMinimumHeight(0);self.control_scroll.setMaximumHeight(16777215);self.control_scroll.setSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
        self.detach_settings_action.setText('Attach parameter console');self.settings_button.setChecked(True);self.update_settings_layout()
        screen=self.screen() if hasattr(self,'screen') else QtWidgets.QApplication.primaryScreen()
        width,height=850,max(220,min(440,self.settings_panel.sizeHint().height()+70))
        if screen is not None:
            rect=screen.availableGeometry();width=min(width,max(650,rect.width()-40));height=min(height,max(180,rect.height()-60))
        self.settings_window.resize(width,height);self.settings_window.show();self.settings_window.raise_()
    def attach_settings(self):
        if not self.settings_detached:return
        self.settings_window.layout().removeWidget(self.control_scroll);self.outer_layout.insertWidget(1,self.control_scroll);self.settings_detached=False
        self.detach_settings_action.setText('Detach parameter console');self.settings_window.hide();self.update_settings_layout()
    def sync_spectrum_x(self,target,limits):
        # Sync the numerical range exactly. ViewBox.setXLink aligns scene pixels
        # and can introduce a small range mismatch when axis widths differ.
        if self.syncing_spectrum_x:return
        self.syncing_spectrum_x=True
        try:target.setXRange(float(limits[0]),float(limits[1]),padding=0)
        finally:self.syncing_spectrum_x=False
    def update_settings_layout(self,*args):
        if not hasattr(self,'control_scroll'):return
        stage=self.stage.currentText();footprint=stage=='Footprint'
        self.calibration_settings.setVisible(stage in ('Tsys calibration','Apply calibration'))
        self.baseline_settings.setVisible(stage=='Baseline')
        self.resample_settings.setVisible(stage=='Table')
        self.settings_button.setEnabled(not footprint);self.detach_settings_action.setEnabled(not footprint)
        self.settings_button.setArrowType(QtCore.Qt.DownArrow if self.settings_button.isChecked() else QtCore.Qt.RightArrow)
        visible=not footprint and self.settings_button.isChecked();self.control_scroll.setVisible(visible)
        if self.settings_detached:
            self.settings_window.setWindowTitle('FISH · parameter console · '+stage);self.settings_window.setVisible(visible)
        self.region_button.setVisible(stage=='Table')
        self.refit_button.setVisible(not footprint);self.baseline_zoom_button.setVisible(stage in ('Baseline','Apply calibration','Table'))
        self.settings_panel.layout().activate()
        # Fixed to content height, with scrolling only for a long component list.
        height=self.component_table.height()+75 if stage=='Baseline' else 42
        if not self.settings_detached:self.control_scroll.setFixedHeight(min(230,height))
    def update_resample_controls(self,*args):
        active=self.stage.currentText()=='Table'
        self.controls['n_b'].setEnabled(active and self.mode.currentText()=='rebin')
        self.ng.setEnabled(active and self.mode.currentText()=='gaussian')
    def stage_changed(self):
        self.region_preview_data=None;self.update_settings_layout()
        footprint=self.stage.currentText()=='Footprint'
        if footprint!=self.footprint_layout:
            self.graph.ci.removeItem(self.map)
            if footprint:
                self.graph.ci.removeItem(self.top);self.graph.ci.removeItem(self.bottom);self.graph.ci.addItem(self.map,row=0,col=0,rowspan=2,colspan=2)
            else:
                self.graph.ci.addItem(self.map,row=0,col=0,rowspan=2);self.graph.ci.addItem(self.top,row=0,col=1);self.graph.ci.addItem(self.bottom,row=1,col=1)
            self.footprint_layout=footprint
        self.top.setVisible(not footprint);self.bottom.setVisible(not footprint);self.roi.setVisible(not footprint and self.region_shape=='rectangle')
        if self.polygon_roi is not None:self.polygon_roi.setVisible(not footprint and self.region_shape=='polygon')
        self.plot_title(self.map,'Pointing footprint');self.map.setToolTip('Click a pointing to inspect its spectrum. Zoom for more markers.')
        for control in list(self.controls.values())+[self.exclusions,self.mode,self.ng,self.pol,self.display_bin,self.spectrum_style,self.component_table,self.add_baseline_button]:control.setEnabled(not footprint)
        self.spectrum_style.setEnabled(not footprint and self.stage.currentText()!='Tsys calibration');self.display_bin.setEnabled(not footprint and self.stage.currentText()!='Tsys calibration')
        self.update_resample_controls()
        if not footprint and hasattr(self,'ra') and self.data_view!=self.stage.currentText():self.fetch('point')
        else:self.draw()
    def fetch(self,op,updates=None):
        if self.closing:return
        request=dict(op=op)
        if op=='point':request['view']=self.stage.currentText()
        if updates:request['updates']=updates
        if self.worker.isRunning():
            # Latest navigation/view wins; draft setting edits must survive it.
            if self.pending_request is not None and op=='point' and self.pending_request['op']=='point':
                previous=self.pending_request.get('updates',{});merged=dict(previous)
                if updates:
                    merged.update(updates)
                    if 'params' in previous and 'params' in updates:merged['params']=dict(previous['params'],**updates['params'])
                if merged:request['updates']=merged
            self.pending_request=request
            self.statusBar().showMessage('Computing; latest requested spectrum/view is queued.');return
        self.dispatch_request(request)
    def dispatch_request(self,request):
        if self.closing:return
        request=dict(request)
        if 'updates' in request:
            if self.state is None:return
            request['revision']=self.state['revision']
        self.worker.request=request
        self.statusBar().showMessage('Computing preview… Zoom and display controls remain local.');self.worker.start()
    def dispatch_pending(self):
        # A prior run's queued finished signal may arrive after a new start.
        if self.closing or self.pending_request is None or self.worker.isRunning():return
        request=self.pending_request;self.pending_request=None
        self.dispatch_request(request)
    def received(self,r):
        if self.closing:return
        self.state=r['state'];d=r['data']
        # Navigation and view responses must not erase uncommitted editor drafts.
        if self.synced_params!=self.state['params']:
            for k,w in self.controls.items():w.setValue(self.state['params'][k])
            self.populate_components(self.state['params'].get('baseline_components',default_components(self.state['params'])))
            self.exclusions.setText(json.dumps(self.state['params'].get('baseline_exclusions_kms',[])));self.mode.setCurrentText(self.state['params']['spectral_mode']);self.ng.setValue(self.state['params']['n_g'])
            self.synced_params=json.loads(json.dumps(self.state['params']))
        if 'region_mean' in d:
            self.region_preview_data=d
            o=np.argsort(d['region_velocity']);self.set_spectrum(self.corrected,d['region_velocity'][o],d['region_mean'][o]);self.plot_title(self.bottom,'Region preview · {} of {} spectra (chunk-sampled)'.format(self.state['region_used'],self.state['region_total']));self.statusBar().showMessage('Region mean only; no gridding or science products written');return
        if 'velocity' in d or 'reference_time' in d:self.region_preview_data=None
        if 'ra' in d:
            self.roi.setPos([float(np.median(d['ra'])),float(np.median(d['dec']))]);
            self.ra=d['ra'];self.dec=d['dec']
            self.scan_groups=d.get('scan_group',np.zeros(len(self.ra),dtype=np.uint8));self.scan_names=self.state.get('scan_names',['Selected input'])
            if len(self.scan_groups)!=len(self.ra) or (len(self.scan_groups) and int(self.scan_groups.max())>=len(self.scan_names)):raise ValueError('Footprint scan labels do not match the point coordinates')
            for curve in self.scan_curves.values():self.map.removeItem(curve)
            self.scan_curves={};self.scan_legend_action.setEnabled(True)
            if self.scan_legend is not None and self.scan_legend.isVisible():self.populate_scan_legend()
            self.full_footprint()
            if self.stage.currentText()!='Footprint' and self.pending_request is None:self.fetch('point')
            else:self.statusBar().showMessage('Footprint only · {} timestamps · no SCI arrays or calibration read'.format(len(self.ra)))
            return
        self.data=d;self.data_view=self.state.get('preview_view');self.render_cache.clear();self.render_cache_bytes=0;self.grid_orders.clear()
        if self.data_view is None:self.data_view=self.stage.currentText()
        if self.data_view==self.stage.currentText():self.draw()
        count=len(d.get('velocity',d.get('reference_time',[])));self.statusBar().showMessage('Spectrum {} · {} samples · revision {} · display changes use local arrays'.format(self.state['index'],count,self.state['revision']))
    def full_footprint(self):
        if not hasattr(self,'ra'):return
        good=np.isfinite(self.ra)&np.isfinite(self.dec)
        if not good.any():return
        self.map.vb.enableAutoRange(enable=False)
        self.map.setRange(xRange=(float(self.ra[good].min()),float(self.ra[good].max())),yRange=(float(self.dec[good].min()),float(self.dec[good].max())),padding=.04)
        self.update_footprint()
    def update_footprint(self,*args):
        if not hasattr(self,'ra'):return
        from .footprint_display import visible_indices
        ranges=self.map.vb.viewRange();rect=self.map.vb.sceneBoundingRect()
        self.displayed_indices=visible_indices(self.ra,self.dec,ranges[0],ranges[1],rect.width(),rect.height())
        ids=self.displayed_indices
        colored=self.scan_colors_button.isChecked() and self.scan_groups is not None
        self.footprint_curve.setVisible(not colored)
        active=set()
        if colored:
            groups=self.scan_groups[ids]
            for group in np.unique(groups):
                group=int(group);active.add(group);selected=ids[groups==group]
                if group not in self.scan_curves:
                    self.scan_curves[group]=self.map.plot(pen=None,symbol='o',symbolSize=2,symbolBrush=scan_color(self.scan_names[group],self.theme.currentIndex()==1),symbolPen=None)
                curve=self.scan_curves[group];curve.setData(self.ra[selected],self.dec[selected]);curve.show()
        else:self.footprint_curve.setData(self.ra[ids],self.dec[ids])
        for group,curve in self.scan_curves.items():
            if group not in active:curve.hide()
        self.plot_title(self.map,'Pointing footprint');self.map.setToolTip('{:,} of {:,} markers shown. Zoom for detail; all pointings remain selectable.'.format(len(ids),len(self.ra)))

    def populate_scan_legend(self):
        self.scan_legend_list.clear()
        counts=np.bincount(self.scan_groups.astype(np.int64),minlength=len(self.scan_names))
        for group,count in enumerate(counts):
            if not count:continue
            name=self.scan_names[group];pixmap=QtGui.QPixmap(16,16);pixmap.fill(pg.mkColor(scan_color(name,self.theme.currentIndex()==1)))
            self.scan_legend_list.addItem(QtWidgets.QListWidgetItem(QtGui.QIcon(pixmap),'{}  ·  {:,} spectra'.format(name,int(count))))

    def show_scan_legend(self):
        if self.scan_groups is None:return
        if self.scan_legend is None:
            self.scan_legend=QtWidgets.QDialog(self,QtCore.Qt.Tool);self.scan_legend.setWindowTitle('FISH · scan colors');layout=QtWidgets.QVBoxLayout(self.scan_legend)
            label=QtWidgets.QLabel('Colors identify scans; they do not change the selection.');label.setWordWrap(True);layout.addWidget(label)
            self.scan_legend_list=QtWidgets.QListWidget();layout.addWidget(self.scan_legend_list);self.scan_legend.resize(430,360)
        self.populate_scan_legend();self.scan_legend.show();self.scan_legend.raise_()
    def error(self,msg):self.statusBar().showMessage(msg)
    def navigate(self,delta):
        if self.stage.currentText()=='Footprint':self.statusBar().showMessage('Choose a processing view to inspect spectra');return
        if self.state:
            request=self.pending_request if self.pending_request is not None else self.worker.request if self.worker.isRunning() else None
            index=request.get('updates',{}).get('index',self.state['index']) if request is not None else self.state['index']
            self.fetch('point',dict(index=max(0,min(len(self.ra)-1,index+delta))))
    def clicked(self,event):
        if self.stage.currentText()=='Footprint' or self.map.vb.class_zoom:return
        if event.button()==QtCore.Qt.LeftButton and self.state and self.map.sceneBoundingRect().contains(event.scenePos()):
            p=self.map.vb.mapSceneToView(event.scenePos());dist=((self.ra-p.x())*np.cos(np.deg2rad(p.y())))**2+(self.dec-p.y())**2;self.fetch('point',dict(index=int(np.argmin(dist))))
    def dragged(self,key,region):
        self.controls[key].setValue(max(abs(x) for x in region.getRegion()));self.statusBar().showMessage('Draft window changed; Apply settings / refit commits it. Current spectrum still uses the previous fit.')
    def select_region(self,shape):
        if self.stage.currentText()=='Footprint':
            self.statusBar().showMessage('Choose a processing view to select a spectrum-mean region');return
        self.map.vb.cancel_zoom()
        (x0,x1),(y0,y1)=self.map.vb.viewRange()
        dx=(x1-x0)/4;dy=(y1-y0)/4;cx=(x0+x1)/2;cy=(y0+y1)/2
        self.region_shape=shape
        self.roi.hide()
        if self.polygon_roi is not None:
            self.map.removeItem(self.polygon_roi);self.polygon_roi=None
        if shape=='rectangle':
            self.roi.setPos([cx-dx/2,cy-dy/2]);self.roi.setSize([dx,dy]);self.roi.show()
        else:
            self.polygon_roi=pg.PolyLineROI([[cx-dx/2,cy-dy/2],[cx+dx/2,cy-dy/2],[cx+dx/2,cy+dy/2],[cx-dx/2,cy+dy/2]],closed=True,pen='#f4a340')
            self.map.addItem(self.polygon_roi)
        self.statusBar().showMessage(shape.title()+' selection: drag handles to reshape, then Preview region mean. Polygon edges can add vertices.')
    def region_mean(self):
        if self.stage.currentText()=='Footprint':self.statusBar().showMessage('Choose Table to preview calibrated region means');return
        if self.worker.isRunning():return
        if self.region_shape is None:self.statusBar().showMessage('Put the cursor over the map and press R or P to select a region');return
        if self.region_shape=='polygon':
            vertices=[self.polygon_roi.mapToParent(pos) for _,pos in self.polygon_roi.getLocalHandlePositions()]
            request=dict(op='region',vertices=[[p.x(),p.y()] for p in vertices])
        else:
            p=self.roi.pos();z=self.roi.size();request=dict(op='region',bounds=[p.x(),p.x()+z.x(),p.y(),p.y()+z.y()])
        self.worker.request=request;self.statusBar().showMessage('Computing region preview (at most 128 spectra)…');self.worker.start()
    def refit(self):
        if self.stage.currentText()=='Footprint':self.statusBar().showMessage('Choose a processing view before refitting');return
        if self.state:
            try:extra=json.loads(self.exclusions.text())
            except ValueError:self.error('Exclusions must be JSON, e.g. [[700,900]]');return
            try:components=self.read_components()
            except (ValueError,TypeError) as exc:self.error(str(exc));return
            params={k:w.value() for k,w in self.controls.items()};params['baseline_components']=components;params.update(baseline_exclusions_kms=extra,spectral_mode=self.mode.currentText(),n_g=self.ng.value());self.fetch('point',dict(params=params))
    def baseline_zoom(self):
        if self.data is None or self.stage.currentText() in ('Footprint','Tsys calibration') or self.data_view!=self.stage.currentText():return
        x=self.data['velocity'];y=self.choose(self.data['raw_pol'] if self.stage.currentText()=='Apply calibration' else self.data['pol'])
        windows=[c['range_kms'] for c in self.state['params'].get('baseline_components',[]) if c['range_kms'] is not None]
        low=min(w[0] for w in windows) if windows else float(np.min(x));high=max(w[1] for w in windows) if windows else float(np.max(x))
        xx,yy=self.display_values(x,y);good=(abs(xx)>self.state['params']['protect_kms'])&(xx>=low)&(xx<=high)&np.isfinite(yy)
        if not good.any():self.error('No unprotected finite values in the component ranges');return
        lo,hi=np.percentile(yy[good],[2,98]);pad=max(.1,(hi-lo)*.15);self.top.setXRange(low,high);self.top.setYRange(lo-pad,hi+pad);self.bottom.setXRange(low,high);
        if self.stage.currentText()=='Apply calibration':
            cx,cy=self.display_values(x,self.choose(self.data['pol']));ok=(abs(cx)>self.state['params']['protect_kms'])&(cx>=low)&(cx<=high)&np.isfinite(cy)
            if ok.any():
                ca,cb=np.percentile(cy[ok],[2,98]);cp=max(.1,(cb-ca)*.15);self.bottom.setYRange(ca-cp,cb+cp)
        else:self.bottom.setYRange(-3,3)
        self.statusBar().showMessage('Component detail; protected H I clipped by display only')
    def display_values(self,x,y):
        n=self.display_bin.value()
        if n==1:return x,y
        ids=np.arange(0,len(x),n);counts=np.minimum(n,len(x)-ids);good=np.isfinite(y);weight=np.add.reduceat(good.astype(float),ids);sums=np.add.reduceat(np.where(good,y,0.),ids);out=np.full(len(ids),np.nan);np.divide(sums,weight,out=out,where=weight>0);return np.add.reduceat(x,ids)/counts,out
    def spectrum_values(self,x,y):
        x=np.asarray(x);y=np.asarray(y);centers,values=self.display_values(x,y)
        if self.spectrum_style.currentText()=='Lines' or len(x)<2:return centers,values,None
        if np.any(np.diff(x)<=0):raise ValueError('Spectrum display coordinates must be strictly increasing')
        edges=np.empty(len(x)+1);edges[1:-1]=(x[:-1]+x[1:])/2;edges[0]=x[0]-(x[1]-x[0])/2;edges[-1]=x[-1]+(x[-1]-x[-2])/2
        starts=np.arange(0,len(x),self.display_bin.value())
        return np.r_[edges[starts],edges[-1]],values,'center'
    def set_spectrum(self,curve,x,y):
        xx,yy,mode=self.spectrum_values(x,y)
        # Older PyQtGraph versions do not preserve N+1 step edges through
        # PlotDataItem's clipping/peak-downsampling. Use the explicit display
        # bins for steps; native mode is bounded by one spectrum's channels.
        curve.setData(xx,yy,stepMode=mode,autoDownsample=mode is None,clipToView=mode is None,downsample=1,connect='finite')
    def choose(self,a):return a[:,self.pol.currentIndex()] if self.pol.currentIndex()<2 else np.mean(a,axis=1)
    def prepared_draw(self):
        """Cache display arrays for this point; fitting stays in the worker.

        UI-only changes often revisit the same polarization/style. Keep this
        small LRU bounded by both bytes and entries, and clear it on new data.
        """
        s=self.stage.currentText();key=(s,self.pol.currentIndex(),self.display_bin.value(),self.spectrum_style.currentIndex())
        if key in self.render_cache:
            self.render_cache.move_to_end(key);return self.render_cache[key][0]
        d=self.data;reference=s=='Tsys calibration'
        grid='reference_time' if reference else 'velocity';x=d[grid]
        if grid not in self.grid_orders:self.grid_orders[grid]=np.argsort(x)
        order=self.grid_orders[grid];xs=x[order]
        result=dict(reference=reference,components=[])
        if reference:
            raw=self.choose(d['reference_raw'])[order];smooth=self.choose(d['reference_smooth'])[order]
            result.update(observed=(xs,raw,None),model=(xs,smooth),corrected=(xs,smooth,None))
        else:
            observed=self.choose(d['raw_pol'] if s=='Apply calibration' else d['pol'])[order]
            result['observed']=self.spectrum_values(xs,observed)
            if s=='Apply calibration':result['model']=(np.empty(0),np.empty(0))
            else:
                total=self.choose(d['broad'])[order]+self.choose(d['local'])[order]
                result['model']=self.display_values(xs,total)
            if s in ('Baseline','Table') and 'component_models' in d:
                specs=self.state['params'].get('baseline_components',[])
                for i,model in enumerate(d['component_models']):
                    spec=specs[i];name=str(i+1)+': '+('polynomial n='+str(spec['degree']) if spec['function']=='polynomial' else 'sin/cos')
                    values=self.choose(model)[order];window=spec['range_kms']
                    if window is not None:
                        values=values.copy();values[(xs<window[0])|(xs>window[1])]=np.nan
                    result['components'].append((name,self.display_values(xs,values)))
            if s=='Table':
                if 'output_velocity' not in self.grid_orders:self.grid_orders['output_velocity']=np.argsort(d['output_velocity'])
                order2=self.grid_orders['output_velocity'];x2=d['output_velocity'][order2];y2=d['output'][order2]
            else:x2=xs;y2=self.choose(d['pol'] if s=='Apply calibration' else d['residual'])[order]
            result['corrected']=self.spectrum_values(x2,y2)
        arrays=[]
        for name in ('observed','model','corrected'):arrays.extend(a for a in result[name] if isinstance(a,np.ndarray))
        for _,pair in result['components']:arrays.extend(pair)
        seen=set();size=0
        for a in arrays:
            if id(a) not in seen:seen.add(id(a));size+=a.nbytes
        # One oversized display still renders, but is not retained in the LRU.
        if size<=32*1024**2:
            while self.render_cache and (len(self.render_cache)>=12 or self.render_cache_bytes+size>32*1024**2):
                _,(_,old_size)=self.render_cache.popitem(last=False);self.render_cache_bytes-=old_size
            self.render_cache[key]=(result,size);self.render_cache_bytes+=size
        return result
    def set_prepared_spectrum(self,curve,values):
        xx,yy,mode=values
        curve.setData(xx,yy,stepMode=mode,autoDownsample=mode is None,clipToView=mode is None,downsample=1,connect='finite')
    def update_masks(self,reference):
        windows=[] if reference else [(lo,hi,(235,90,90)) for lo,hi in self.data['rfi_velocity']]+[(lo,hi,(180,100,240)) for lo,hi in self.state['params'].get('baseline_exclusions_kms',[])]
        signature=tuple((float(lo),float(hi),tuple(color)) for lo,hi,color in windows)
        if signature!=self.mask_signature:
            for i,(lo,hi,color) in enumerate(windows):
                if i==len(self.masks):
                    region=pg.LinearRegionItem(movable=False);self.top.addItem(region);self.masks.append(region)
                region=self.masks[i];region.setRegion(sorted([lo,hi]));region.setBrush((*color,40))
                for line in region.lines:line.setPen((*color,80))
                region.show()
            for region in self.masks[len(windows):]:region.hide()
            self.mask_signature=signature
        self.protected.setVisible(not reference)
        if not reference:
            self.protected.blockSignals(True);v=self.state['params']['protect_kms'];self.protected.setRegion([-v,v]);self.protected.blockSignals(False)
    def draw(self):
        if self.data is None or self.stage.currentText()=='Footprint' or self.data_view!=self.stage.currentText():return
        if self.region_preview_data is not None:
            r=self.region_preview_data;o=np.argsort(r['region_velocity']);self.set_spectrum(self.corrected,r['region_velocity'][o],r['region_mean'][o]);return
        d=self.data;s=self.stage.currentText()
        if s=='Apply calibration' and 'raw_pol' not in d:self.error('Reopen the viewer to load raw SCI counts from the updated worker');return
        new_stage=s!=self.last_stage
        if new_stage:self.top.enableAutoRange();self.bottom.enableAutoRange(axis='y');self.last_stage=s
        prepared=self.prepared_draw();reference=prepared['reference'];self.update_masks(reference)
        label='MJD' if reference else 'LSRK velocity';units=None if reference else 'km/s'
        self.plot_title(self.top,'Tpsr · raw / smoothed' if reference else 'Raw SCI counts · '+self.pol.currentText() if s=='Apply calibration' else s+' · '+self.pol.currentText())
        self.top.setLabel('left','Raw SCI counts' if s=='Apply calibration' else 'Temperature',units=None if s=='Apply calibration' else 'K')
        self.bottom.setLabel('left','Temperature',units='K')
        for sample,label_item in self.top.legend.items:
            if sample.item is self.curves[0]:label_item.setText('Raw SCI counts' if s=='Apply calibration' else 'Raw Tpsr' if reference else 'Observed')
        for sample,label_item in self.bottom.legend.items:
            if sample.item is self.corrected:label_item.setText('Calibrated (K)' if s=='Apply calibration' else 'Smoothed Tpsr' if reference else 'Corrected (K)')
        self.set_prepared_spectrum(self.curves[0],prepared['observed'])
        self.curves[1].setData(*prepared['model'],stepMode=None);self.curves[1].setVisible(s!='Apply calibration')
        for sample,label_item in self.top.legend.items:
            if sample.item is self.curves[1]:
                sample.setVisible(s!='Apply calibration');label_item.setVisible(s!='Apply calibration');label_item.setText('Smoothed Tpsr' if reference else 'Total baseline')
        self.top.setToolTip('Band medians are diagnostic summaries; Tpsr smoothing is per frequency channel.' if reference else 'A component period at the search limit needs inspection.' if self.state.get('fit_report',{}).get('at_search_boundary') else '')
        # Retain graphics and legend entries when switching point or style.
        for i,(name,pair) in enumerate(prepared['components']):
            if i==len(self.component_curves):
                curve=self.top.plot(pen=pg.mkPen(pg.intColor(i,hues=max(6,len(prepared['components']))),width=1),name=name,autoDownsample=True,downsampleMethod='peak',clipToView=True)
                self.component_curves.append(curve)
            curve=self.component_curves[i];curve.setData(*pair);curve.show()
            curve.setPen(pg.mkPen(pg.intColor(i,hues=max(6,len(prepared['components']))),width=1))
            for sample,label_item in self.top.legend.items:
                if sample.item is curve:label_item.setText(name);sample.show();label_item.show()
        for curve in self.component_curves[len(prepared['components']):]:
            curve.hide()
            for sample,label_item in self.top.legend.items:
                if sample.item is curve:sample.hide();label_item.hide()
        self.set_prepared_spectrum(self.corrected,prepared['corrected'])
        if new_stage and s=='Table':
            output=d['output_velocity'];self.top.setXRange(float(np.min(output)),float(np.max(output)),padding=.02)
        self.top.setLabel('bottom',label,units=units);self.bottom.setLabel('bottom',label,units=units)
        self.plot_title(self.bottom,'Calibrated temperature · '+self.pol.currentText() if s=='Apply calibration' else 'Resampled product' if s=='Table' else 'Smoothed Tpsr' if reference else 'Baseline-subtracted spectrum')
        self.marker.setData([self.ra[self.state['index']]],[self.dec[self.state['index']]])
        if s in ('Baseline','Table') and self.state.get('fit_report',{}).get('at_search_boundary'):self.plot_title(self.top,self.top.titleLabel.text+' · period at limit')
    def closeEvent(self,event):
        if self.image_loader is not None and self.image_loader.isRunning():
            self.statusBar().showMessage('Wait for the FITS image to finish loading before closing');event.ignore();return
        self.closing=True;self.pending_request=None;self.attach_settings()
        if self.scan_legend is not None:self.scan_legend.hide()
        QtWidgets.QApplication.instance().removeEventFilter(self.keyboard)
        self.worker.shutdown();event.accept()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--remote');p.add_argument('--view',choices=['where','baseline','reference','calibration','resample'],default='where');p.add_argument('--session-state',required=True);p.add_argument('--remote-source',default='/bb8/lingrui/my_projects/FAST_pipeline_v0.0/src');a=p.parse_args()
    command=[sys.executable,'-m','fastfish.viewer_worker','--session-state',a.session_state]
    if a.remote and a.remote.startswith('-'):p.error('Invalid SSH host')
    if a.remote:
        remote='OPENBLAS_NUM_THREADS=1 PYTHONPATH={} /usr/bin/python3 -m fastfish.viewer_worker --session-state {}'.format(shlex.quote(a.remote_source),shlex.quote(a.session_state));command=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15',a.remote,remote]
    pg.setConfigOptions(background='#ffffff',foreground='#000000',antialias=False);app=pg.mkQApp('FISH');window=Viewer(command,a.view);window.fit_screen();window.show();QtCore.QTimer.singleShot(150,window.fit_screen);return app.exec_()
if __name__=='__main__':raise SystemExit(main())
