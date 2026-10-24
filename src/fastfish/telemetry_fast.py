"""Stream only FAST's five telemetry columns from XLSX worksheet XML."""
from pathlib import Path
import posixpath,re,zipfile
import xml.etree.ElementTree as ET
import numpy as np
from astropy.time import Time
import astropy.units as u
from .reference import sha

NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
KEYS=['SysTime','SDP_PhaPos_X','SDP_PhaPos_Y','SDP_PhaPos_Z','SDP_AngleM']


def telemetry_rows(path):
    with zipfile.ZipFile(path) as z:
        workbook=ET.fromstring(z.read('xl/workbook.xml'));views=workbook.find(NS+'bookViews');active=int(views[0].get('activeTab','0')) if views is not None and len(views) else 0
        sheet=list(workbook.find(NS+'sheets'))[active];rid=sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'));target=next(r.get('Target') for r in rels if r.get('Id')==rid)
        target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        props=workbook.find(NS+'workbookPr');date1904=props is not None and props.get('date1904') in ('1','true')
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            with z.open('xl/sharedStrings.xml') as stream:
                for event,node in ET.iterparse(stream,events=('end',)):
                    if node.tag==NS+'si':strings.append(''.join(t.text or '' for t in node.iter(NS+'t')));node.clear()
        def value(cell):
            kind=cell.get('t');v=cell.find(NS+'v')
            if kind=='inlineStr':return ''.join(t.text or '' for t in cell.iter(NS+'t'))
            if v is None or v.text is None:return None
            if kind=='s':return strings[int(v.text)]
            if kind in ('str','d','e'):return v.text
            return float(v.text)
        columns=None;sheetdata=None
        with z.open(target) as stream:
            for event,node in ET.iterparse(stream,events=('start','end')):
                if event=='start':
                    if node.tag==NS+'sheetData':sheetdata=node
                    continue
                if node.tag!=NS+'row':continue
                if columns is None:
                    head={}
                    for c in node:
                        if c.tag==NS+'c':head.setdefault(value(c),re.sub(r'\d+$','',c.get('r','')))
                    missing=[key for key in KEYS if key not in head]
                    if missing:raise ValueError('Missing telemetry columns: '+', '.join(missing))
                    columns={head[k]:k for k in KEYS};node.clear()
                    if sheetdata is not None:sheetdata.clear()
                    continue
                row={k:None for k in KEYS}
                for cell in node:
                    col=re.sub(r'\d+$','',cell.get('r',''))
                    if col in columns:row[columns[col]]=value(cell)
                node.clear()
                if sheetdata is not None:sheetdata.clear()
                if row['SysTime'] is None:continue
                if isinstance(row['SysTime'],float):
                    from openpyxl.utils.datetime import from_excel,CALENDAR_MAC_1904,CALENDAR_WINDOWS_1900
                    row['SysTime']=str(from_excel(row['SysTime'],epoch=CALENDAR_MAC_1904 if date1904 else CALENDAR_WINDOWS_1900))
                row['SysTime']=str(row['SysTime']).replace('T',' ')
                yield row


def telemetry_cache(scan,sheet,tmin,tmax,cache_dir):
    out=Path(cache_dir)/(scan+'.npz');out.parent.mkdir(parents=True,exist_ok=True)
    lower=(Time(tmin,format='mjd')+8*u.hr-5*u.s).iso;upper=(Time(tmax,format='mjd')+8*u.hr+5*u.s).iso
    data={k:[] for k in KEYS}
    rows=telemetry_rows(sheet)
    try:
        for row in rows:
            value=row['SysTime']
            if value<lower:continue
            if value>upper:break
            for k in KEYS:data[k].append(row[k])
    finally:rows.close()
    if len(data['SysTime'])<2:raise ValueError('No matching telemetry: '+scan)
    arrays={k:np.asarray(v,dtype=str if k=='SysTime' else float) for k,v in data.items()};np.savez_compressed(out,**arrays)
    print('TELEMETRY',scan,len(data['SysTime']),'(streamed five columns)',flush=True)
    return dict(path=str(out),rows=len(data['SysTime']),sha256=sha(out),source=str(sheet),source_bytes=Path(sheet).stat().st_size,system_time_first_last=[data['SysTime'][0],data['SysTime'][-1]])
