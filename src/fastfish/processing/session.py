"""Atomic administrative session metadata shared by the FISH prompt and plots."""
from pathlib import Path
import fcntl,json,os,tempfile

def load_session_catalog(session_path,path,profile=None,expected_revision=None):
    """Select a validated, frozen reference without changing original configs/products."""
    import hashlib
    from fastfish.rfi.catalog import catalog_path,validate_catalog,profile_priors,load_published_statistics
    s=read_session(session_path)
    if expected_revision is not None and s['revision']!=expected_revision:
        raise ValueError('Session changed; refresh before loading a catalogue')
    root=Path(s['root']).resolve();job=root/'active_job.json'
    if job.is_file():
        record=json.loads(job.read_text());proc=Path('/proc')/str(record.get('pid',0))/'cmdline'
        if proc.is_file():
            command=proc.read_bytes().replace(b'\0',b' ').decode(errors='replace')
            if ('fastfish.native_stages' in command or 'fastfish.processing.stages' in command) and str(root) in command:
                raise ValueError('A stage is running. Wait for it before changing the RFI catalogue.')
    original=Path(s['config']).resolve();cfg=json.loads(original.read_text())
    source=catalog_path(path,original.parent);payload=source.read_bytes();catalog=validate_catalog(json.loads(payload))
    profile=profile.strip() if profile else None
    if profile:profile_priors(catalog,profile,cfg['scans'])
    published=None
    if 'published_statistics' in catalog:
        load_published_statistics(catalog,source.parent)
        published=(source.parent/catalog['published_statistics']['path']).read_bytes()
    # Preserve config-relative inputs when moving the edited config into this run.
    for key in ('input_root','extracted_cal_root','tcal_path','output_dir','fish_run_root'):
        if cfg.get(key):cfg[key]=str((original.parent/Path(cfg[key]).expanduser()).resolve())
    if cfg.get('pointing'):cfg['pointing']={k:str((original.parent/Path(v).expanduser()).resolve()) for k,v in cfg['pointing'].items()}
    for item in cfg.get('inputs',[]):
        if item.get('reference'):item['reference']=str((original.parent/Path(item['reference']).expanduser()).resolve())
        item['files']=[str((original.parent/Path(f).expanduser()).resolve()) for f in item['files']]
    digest=hashlib.sha256(payload).hexdigest();folder=root/'rfi_catalogues'/digest
    cfg['rfi_catalog']=str(folder/'rfi_catalog.json');cfg['rfi_profile']=profile
    cfg.update(rfi_catalog_path=cfg['rfi_catalog'],rfi_catalog_sha256=digest,rfi_catalog_version=catalog['catalog_version'])
    config_bytes=(json.dumps(cfg,sort_keys=True,indent=2)+'\n').encode()
    config=root/'rfi_catalogues'/('config_'+hashlib.sha256(config_bytes).hexdigest()+'.json')
    snapshots={folder/'rfi_catalog.json':payload,config:config_bytes}
    if published is not None:snapshots[folder/catalog['published_statistics']['path']]=published
    for target,data in snapshots.items():
        target.parent.mkdir(parents=True,exist_ok=True)
        try:
            with target.open('xb') as stream:stream.write(data)
        except FileExistsError:
            if target.read_bytes()!=data:raise ValueError('Conflicting catalogue snapshot: '+str(target))
    info=dict(source=str(source),snapshot=str(folder/'rfi_catalog.json'),version=catalog['catalog_version'],sha256=digest,profile=profile)
    return write_session(session_path,dict(config=str(config),rfi_catalog_info=info),s['revision'])

def read_session(path):
    value=json.loads(Path(path).read_text())
    if value.get('format')!='fastfish-session-v1':raise ValueError('Not a FISH session file')
    return value

def write_session(path,updates,expected_revision=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.with_suffix(path.suffix+'.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        current=read_session(path) if path.exists() else dict(format='fastfish-session-v1',revision=0)
        if expected_revision is not None and current['revision']!=expected_revision:
            raise ValueError('Session settings changed in another window or the terminal. Reopen this plot before editing.')
        if all(current.get(k)==v for k,v in updates.items()):return current
        current.update(updates);current['revision']+=1
        fd,name=tempfile.mkstemp(prefix='.fish-session-',dir=str(path.parent))
        try:
            with os.fdopen(fd,'w') as stream:json.dump(current,stream,indent=2);stream.write('\n')
            os.replace(name,str(path))
        finally:
            if os.path.exists(name):os.unlink(name)
        return current
