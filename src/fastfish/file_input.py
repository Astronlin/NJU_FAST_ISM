"""Resolve native FITS inputs without loading science arrays."""
from pathlib import Path
import copy,glob,json

def native_name(value):
    name=str(value).lower()
    return any(name.endswith(s) for s in ('.fits','.fit','.fits.gz','.fit.gz')) or any(c in name for c in '*?[')

def expand_inputs(words,base_dir=None):
    if words[0]=='/list':
        if len(words)!=2:raise ValueError('Usage: file in /list LISTNAME')
        listing=Path(words[1]).expanduser();listing=((Path(base_dir) if base_dir else Path.cwd())/listing).resolve() if not listing.is_absolute() else listing.resolve()
        patterns=[(line.strip(),listing.parent) for line in listing.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    else:patterns=[(word,Path(base_dir) if base_dir else Path.cwd()) for word in words]
    result=[]
    for pattern,base in patterns:
        path=Path(pattern).expanduser();path=path if path.is_absolute() else base/path
        matches=sorted(glob.glob(str(path)))
        if not matches:raise ValueError('No file matches '+str(path))
        for match in matches:
            p=Path(match).resolve()
            if not p.is_file() or not native_name(p) or any(c in p.name for c in '*?['):raise ValueError('Expected a FITS file: '+str(p))
            if p not in result:result.append(p)
    if not result:raise ValueError('Empty file list')
    return result

def absolute_config(config_path):
    path=Path(config_path).resolve();cfg=json.loads(path.read_text());base=path.parent
    def absolute(v):return str((base/Path(v).expanduser()).resolve())
    for key in ('input_root','extracted_cal_root','tcal_path'):
        if cfg.get(key):cfg[key]=absolute(cfg[key])
    if cfg.get('rfi_catalog','packaged')!='packaged':cfg['rfi_catalog']=absolute(cfg['rfi_catalog'])
    cfg['pointing']={k:absolute(v) for k,v in cfg.get('pointing',{}).items()}
    for item in cfg.get('inputs',[]):
        item['reference']=absolute(item['reference']);item['files']=[absolute(f) for f in item['files']]
    return cfg

def subset_config(config_path,files):
    cfg=absolute_config(config_path);wanted={str(p.resolve()) for p in files};known={f for i in cfg['inputs'] for f in i['files']}
    missing=wanted-known
    if missing:raise ValueError('Files have no calibration/pointing association in this observation: '+', '.join(sorted(missing))+'. Load the matching observation or use --config CONFIG.')
    cfg=copy.deepcopy(cfg);cfg['inputs']=[dict(i,files=[f for f in i['files'] if f in wanted]) for i in cfg['inputs'] if any(f in wanted for f in i['files'])]
    return cfg

def cached_subset(selection_path,files):
    s=json.loads(Path(selection_path).read_text());wanted={str(p.resolve()) for p in files};items=[];count=0
    import numpy as np
    for item in s['items']:
        contexts=[c for c in item['files'] if str(Path(c['entry']['path']).resolve()) in wanted]
        if contexts:
            items.append(dict(item,files=contexts))
            for c in contexts:
                with np.load(c['coordinates']) as z:count+=len(z['ids'])
    if not count:raise ValueError('These files have no spectra inside the current observation selection')
    return dict(s,items=items,spectra=count,input_files=[str(p) for p in files])
