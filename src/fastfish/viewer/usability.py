"""Lightweight guidance and local controls for the Qt spectral workspace."""
import json
import shlex
import time
import numpy as np
from pyqtgraph.Qt import QtCore,QtWidgets

STEPS=['Footprint','Tsys calibration','Apply calibration','Baseline','Table']
GUIDANCE={
    'Footprint':'Check the sky coverage, then continue to Tsys calibration. Scan colors identify the loaded scans.',
    'Tsys calibration':'Click a pointing to inspect the temperature reference. Adjust the smoothing windows, then Update preview.',
    'Apply calibration':'Compare raw counts above with calibrated temperature below. Inspect AA and BB separately.',
    'Baseline':'Add functions in fit order. Inspect RFI/protected bands and residuals, then Update preview.',
    'Table':'Choose spectral smoothing and output velocity limits. R/P selects a region for a preview mean.',
}
COMMANDS=['tsyscal run','wait','applycal run','wait','baseline run','wait','table run','wait','xy_map run --validate','wait','status']

class UsabilityMixin:
    def setup_guidance(self,toolbar_layout,resample_bar):
        self.draft_dirty=False;self.syncing_controls=False;self.submitted_draft=None;self.last_preview_error='';self.info_dialogs=[]
        guide=QtWidgets.QWidget();row=QtWidgets.QHBoxLayout(guide);row.setContentsMargins(0,0,0,0)
        self.guide_label=QtWidgets.QLabel();self.guide_label.setWordWrap(True);row.addWidget(self.guide_label,1)
        self.step_next=QtWidgets.QPushButton('Next step →');self.step_next.clicked.connect(self.next_step);row.addWidget(self.step_next)
        self.commands_button=QtWidgets.QPushButton('Run commands…');self.commands_button.clicked.connect(self.show_run_commands);row.addWidget(self.commands_button)
        toolbar_layout.addWidget(guide)
        details=QtWidgets.QHBoxLayout();toolbar_layout.addLayout(details)
        self.selection_label=QtWidgets.QLabel('Connecting to data…');self.selection_label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse);details.addWidget(self.selection_label,1)
        self.draft_label=QtWidgets.QLabel();details.addWidget(self.draft_label)
        self.reset_draft_button=QtWidgets.QPushButton('Reset edits');self.reset_draft_button.clicked.connect(self.reset_draft);details.addWidget(self.reset_draft_button)
        self.notice=QtWidgets.QWidget();notice_row=QtWidgets.QHBoxLayout(self.notice);notice_row.setContentsMargins(6,2,6,2)
        self.notice_text=QtWidgets.QLabel();self.notice_text.setWordWrap(True);notice_row.addWidget(self.notice_text,1)
        self.fit_details_button=QtWidgets.QPushButton('Details…');self.fit_details_button.clicked.connect(self.show_fit_details);notice_row.addWidget(self.fit_details_button)
        toolbar_layout.addWidget(self.notice);self.notice.hide()
        for name in ('vmin','vmax'):
            spin=QtWidgets.QDoubleSpinBox();spin.setRange(-100000,100000);spin.setDecimals(3);spin.setSuffix(' km/s');setattr(self,name,spin)
        resample_bar.addWidget(QtWidgets.QLabel('Output velocity from'));resample_bar.addWidget(self.vmin);resample_bar.addWidget(QtWidgets.QLabel('to'));resample_bar.addWidget(self.vmax)
        tips={
            'n_t_cal':'Number of CAL timestamps used to smooth ON−OFF when deriving Tsys. 1 disables this averaging.',
            'n_t_apply':'Number of reference timestamps used to smooth Tpsr before matching SCI times. This is separate from ON−OFF smoothing.',
            'protect_kms':'H I half-width around 0 km/s in LSRK. These channels cannot anchor a baseline fit.',
            'n_b':'Number of native channels per output bin. This changes the science spectrum and output channel spacing.',
        }
        for key,w in self.controls.items():w.setToolTip(tips[key]);w.valueChanged.connect(self.update_draft_status)
        self.exclusions.setToolTip('Baseline-only exclusions in LSRK km/s, e.g. [[700,900],[-1200,-1000]]. Use [] for none.')
        self.exclusions.textChanged.connect(self.update_draft_status)
        self.mode.currentTextChanged.connect(self.update_draft_status)
        for w in (self.ng,self.vmin,self.vmax):w.valueChanged.connect(self.update_draft_status)
        self.ng.setToolTip('Target Gaussian FWHM in native-channel units. Output spacing equals this FWHM.')
        self.display_bin.setToolTip('Average this many channels only for drawing. 1 shows native channels; science products are unaffected.')
        self.refit_button.setToolTip('Save the edited settings to this session and recompute the selected preview. Full processing uses Run commands.')
        self.busy_label=QtWidgets.QLabel();self.statusBar().addPermanentWidget(self.busy_label)
        self.busy_bar=QtWidgets.QProgressBar();self.busy_bar.setRange(0,0);self.busy_bar.setFixedWidth(85);self.busy_bar.hide();self.statusBar().addPermanentWidget(self.busy_bar)
        self.busy_timer=QtCore.QTimer(self);self.busy_timer.setInterval(1000);self.busy_timer.timeout.connect(self.update_busy)
        self.worker.started.connect(self.preview_started);self.worker.finished.connect(self.preview_finished)
        menu=self.menuBar().addMenu('View')
        display=menu.addAction('Display tools (style, colors, images)');display.setCheckable(True);display.toggled.connect(self.display_tools.setVisible)
        menu.addAction('Plot ranges…',self.show_plot_ranges)
        menu.addAction('Keyboard / mouse help',lambda:self.keyboard.handle('H'))
        workflow=self.menuBar().addMenu('Workflow');workflow.addAction('Run commands…',self.show_run_commands);workflow.addAction('Fit details…',self.show_fit_details);workflow.addAction('RFI catalogue…',self.load_rfi_catalog)
        self.update_guidance();self.update_draft_status()

    def draft_parameters(self):
        try:extra=json.loads(self.exclusions.text())
        except ValueError:raise ValueError('Extra exclusions: use pairs such as [[700,900]], or [] for none.') from None
        from fastfish.baseline.windows import normalize_windows
        extra=normalize_windows(extra)
        if self.vmin.value()>=self.vmax.value():raise ValueError('Output velocity: the lower limit must be smaller than the upper limit.')
        p={k:w.value() for k,w in self.controls.items()}
        p.update(baseline_components=self.read_components(),baseline_exclusions_kms=extra,spectral_mode=self.mode.currentText(),n_g=self.ng.value(),velocity_range_kms=[self.vmin.value(),self.vmax.value()])
        return p

    def update_draft_status(self,*args):
        if not hasattr(self,'draft_label') or self.syncing_controls:return
        if self.state is None:self.draft_label.setText('');self.reset_draft_button.setEnabled(False);return
        try:
            p=self.draft_parameters();self.draft_dirty=any(p[k]!=self.state['params'].get(k) for k in p)
            text='Unapplied edits · Update preview' if self.draft_dirty else 'Settings applied'
        except (ValueError,TypeError,KeyError) as exc:self.draft_dirty=True;text='Check settings';self.draft_label.setToolTip(str(exc))
        else:self.draft_label.setToolTip('Edits are kept until you update the preview or reset them.')
        self.draft_label.setText(text);self.draft_label.setStyleSheet('color: #9a5800;' if self.draft_dirty else 'color: #25623a;')
        self.reset_draft_button.setEnabled(self.draft_dirty)
        self.commands_button.setEnabled(not self.draft_dirty and not self.worker.isRunning() and self.pending_request is None)

    def sync_parameter_controls(self,params):
        self.syncing_controls=True
        try:
            for k,w in self.controls.items():w.setValue(params[k])
            self.populate_components(params.get('baseline_components',[]))
            self.exclusions.setText(json.dumps(params.get('baseline_exclusions_kms',[])))
            self.mode.setCurrentText(params['spectral_mode']);self.ng.setValue(params['n_g'])
            self.vmin.setValue(params['velocity_range_kms'][0]);self.vmax.setValue(params['velocity_range_kms'][1])
        finally:self.syncing_controls=False
        self.update_draft_status()

    def reset_draft(self):
        if self.state is not None:
            self.sync_parameter_controls(self.state['params']);self.submitted_draft=None
            self.last_preview_error='';self.update_notice()
            self.update_masks(self.stage.currentText()=='Tsys calibration') if self.data is not None and 'rfi_velocity' in self.data else None
            self.statusBar().showMessage('Edits reset to the settings saved in this session.')

    def update_guidance(self):
        if not hasattr(self,'guide_label'):return
        stage=self.stage.currentText();i=STEPS.index(stage)
        self.guide_label.setText('{} / 5 · {}'.format(i+1,GUIDANCE[stage]));self.step_next.setEnabled(True)
        self.step_next.setText('Next step →' if i<len(STEPS)-1 else 'Region mean')
        self.step_next.setToolTip('Open '+STEPS[i+1] if i<len(STEPS)-1 else 'Preview the mean of the rectangle/polygon selected on the map (at most 128 spectra).')
        self.graph.ci.layout.setColumnMaximumWidth(0,16777215 if stage=='Footprint' else max(230,self.graph.width()*.37))
        self.update_notice()

    def next_step(self):
        i=STEPS.index(self.stage.currentText())
        if i+1<len(STEPS):self.stage.setCurrentText(STEPS[i+1])
        else:self.region_mean()

    def resizeEvent(self,event):
        super().resizeEvent(event)
        if hasattr(self,'graph'):
            self.graph.ci.layout.setColumnMaximumWidth(0,16777215 if self.footprint_layout else max(230,self.graph.width()*.37))

    def update_selection_info(self):
        if self.state is None:return
        info=self.state.get('selection_info',{});index=int(self.state.get('index',0))
        if info:text='{} / {} · raw row {} · spectrum {:,}'.format(info['scan'],info['beam'],info['raw_row'],index)
        elif hasattr(self,'ra'):text='{:,.0f} pointings · spectrum {:,} · RA {:.5f}° Dec {:.5f}°'.format(len(self.ra),index,self.ra[index],self.dec[index])
        else:text='Session loaded'
        self.selection_label.setText(text)
        host='remote server' if self.worker.command and self.worker.command[0]=='ssh' else 'local data host'
        self.selection_label.setToolTip('Data: {}\nFile: {}\nRun folder: {}\nSpectrum and raw-row indices start at 0.'.format(host,info.get('file','—'),self.state.get('root','—')))

    def preview_started(self):
        self.busy_since=time.monotonic();self.busy_bar.show();self.busy_timer.start();self.update_busy();self.update_draft_status()
    def update_busy(self):
        self.busy_label.setText('Preview · {} s'.format(int(time.monotonic()-self.busy_since)))
    def preview_finished(self):
        if self.worker.isRunning():return
        self.busy_timer.stop();self.busy_bar.hide();self.busy_label.setText('Ready');self.update_draft_status()

    def update_notice(self):
        if not hasattr(self,'notice'):return
        failed=self.last_preview_error
        warning=bool(self.state and self.state.get('fit_report',{}).get('at_search_boundary') and self.stage.currentText() in ('Baseline','Table') and self.data_view==self.stage.currentText())
        self.notice.setVisible(bool(failed or warning))
        self.notice.setStyleSheet('background: #fff1e6; color: #743b00;' if not failed else 'background: #fff0f0; color: #8c2020;')
        self.notice_text.setText('Preview could not be updated: '+failed if failed else 'Period at search limit. Inspect the residuals and open Details for the affected component / polarization.')

    def text_dialog(self,title,text):
        dialog=QtWidgets.QDialog(self);dialog.setWindowTitle(title);dialog.resize(680,450);layout=QtWidgets.QVBoxLayout(dialog)
        editor=QtWidgets.QPlainTextEdit();editor.setReadOnly(True);editor.setPlainText(text);layout.addWidget(editor)
        buttons=QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close);buttons.rejected.connect(dialog.close);layout.addWidget(buttons)
        copy=buttons.addButton('Copy text',QtWidgets.QDialogButtonBox.ActionRole);copy.clicked.connect(lambda:QtWidgets.QApplication.clipboard().setText(editor.toPlainText()))
        self.info_dialogs.append(dialog);dialog.finished.connect(lambda _:self.info_dialogs.remove(dialog) if dialog in self.info_dialogs else None);dialog.show()
        return dialog,editor,buttons

    def show_fit_details(self):
        from fastfish.viewer.qt_viewer import period_warning
        report=(self.state or {}).get('fit_report',{})
        text=self.last_preview_error or period_warning(report) or 'No period-boundary warning is reported for the displayed fit.'
        if report:text+='\n\nFit report (AA and BB):\n'+json.dumps(report,indent=2)
        else:text+='\nOpen Baseline and update the preview to inspect its fit report.'
        self.text_dialog('FISH · preview / fit details',text)

    def show_run_commands(self):
        if self.state is None:return
        if self.state.get('preview_only'):
            self.text_dialog('Compact example · previews only','This 32-spectrum teaching session supports interactive previews. Use a full observation project for batch reduction and cube gridding.');return
        if self.worker.isRunning() or self.pending_request is not None:
            self.error('Wait for the current preview to finish before preparing processing commands.');return
        if self.draft_dirty:self.error('Update preview to apply your edits before preparing processing commands.');return
        commands=list(COMMANDS)
        if self.state.get('omit_gridding'):
            commands=commands[:commands.index('xy_map run --validate')]+['status']
        if self.state.get('session_path'):commands.insert(0,'restore '+shlex.quote(self.state['session_path']))
        text='Run these commands in the FISH prompt on the data host.\nEach run processes the entire loaded selection with the saved settings.\nWait for a stage to finish before the next command; investigate any failure.\n\nRun folder: '+self.state.get('root','—')+'\n\n'+'\n'.join(commands)
        dialog,editor,buttons=self.text_dialog('FISH · process the loaded selection',text)
        copy=buttons.addButton('Copy commands only',QtWidgets.QDialogButtonBox.ActionRole);copy.clicked.connect(lambda:QtWidgets.QApplication.clipboard().setText('\n'.join(commands)))

    def show_plot_ranges(self):
        if self.stage.currentText()=='Footprint':self.error('Open a spectrum view to set its plot ranges.');return
        dialog=QtWidgets.QDialog(self);dialog.setWindowTitle('Spectrum display ranges');layout=QtWidgets.QFormLayout(dialog);fields=[]
        reference=self.stage.currentText()=='Tsys calibration'
        ranges=[('Time (MJD)' if reference else 'Velocity (km/s)',self.top.vb.viewRange()[0]),('Upper spectrum Y',self.top.vb.viewRange()[1]),('Lower spectrum Y',self.bottom.vb.viewRange()[1])]
        for label,limits in ranges:
            row=QtWidgets.QWidget();bar=QtWidgets.QHBoxLayout(row);bar.setContentsMargins(0,0,0,0);pair=[]
            for value in limits:
                w=QtWidgets.QLineEdit('{:.12g}'.format(value));bar.addWidget(w);pair.append(w)
            layout.addRow(label+' · min / max',row);fields.append(pair)
        info=QtWidgets.QLabel('Display only. The shared X range and independent Y ranges do not change the output cube.');info.setWordWrap(True);layout.addRow(info)
        error=QtWidgets.QLabel();layout.addRow(error)
        buttons=QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Apply|QtWidgets.QDialogButtonBox.Close);layout.addRow(buttons);buttons.rejected.connect(dialog.close)
        def apply():
            try:
                limits=[[float(w.text()) for w in pair] for pair in fields]
                if any(not np.isfinite(pair).all() or pair[0]>=pair[1] for pair in limits):raise ValueError()
            except ValueError:error.setText('Enter finite numbers with minimum < maximum for each range.');return
            for plot in (self.top,self.bottom):plot.vb.remember_view()
            self.top.setXRange(*limits[0],padding=0);self.top.setYRange(*limits[1],padding=0);self.bottom.setYRange(*limits[2],padding=0);error.setText('Display ranges applied. B restores the previous view.')
        buttons.button(QtWidgets.QDialogButtonBox.Apply).clicked.connect(apply)
        self.info_dialogs.append(dialog);dialog.finished.connect(lambda _:self.info_dialogs.remove(dialog) if dialog in self.info_dialogs else None);dialog.show()
