"""Validated, immutable coordinate products shared across project sessions."""
from pathlib import Path
import hashlib,json,os,tempfile


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def context(project,sheets):
    import astropy,numpy,erfa
    package=Path(__file__).parent
    sources={name:digest(package/name) for name in ('xyz2celestial.py','project_input.py','telemetry_fast.py','fits_stream.py','time_tables.py','iers_b.fits')}
    leap=Path.home()/'.astropy/Leap_Second.dat'
    return dict(schema=1,project=str(Path(project).resolve()),code=sources,
        versions=dict(astropy=astropy.__version__,numpy=numpy.__version__,erfa=erfa.__version__),
        leap_seconds_sha256=digest(leap) if leap.is_file() else None,
        telemetry={scan:dict(path=str(Path(p).resolve()),sha256=digest(p)) for scan,p in sheets.items()})


def candidate(project,scan,beam,entry,common):
    path=Path(entry['path']).resolve();stat=path.stat()
    signature=dict(schema=common['schema'],project=common['project'],code=common['code'],versions=common['versions'],leap_seconds_sha256=common['leap_seconds_sha256'],telemetry=common['telemetry'][scan],scan=scan,beam=beam,
        science=dict(path=str(path),size=stat.st_size,mtime_ns=stat.st_mtime_ns,ctime_ns=stat.st_ctime_ns,metadata=dict(entry,path=str(path))))
    key=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    folder=Path(project)/'fish_runs'/'cache'/'coordinates'/key[:2]
    return dict(path=str(folder/(key+'.npz')),manifest=str(folder/(key+'.json')),signature=signature)


def lookup(cache):
    locations=[(cache['path'],cache['manifest'])]
    for path,manifest in locations:
        try:
            p=Path(path);record=json.loads(Path(manifest).read_text())
            if record['signature']!=cache['signature'] or digest(p)!=record['sha256']:continue
            rows=cache['signature']['science']['metadata']['rows']
            if not isinstance(record['kept'],int) or not 0<=record['kept']<=rows or record['excluded']!=rows-record['kept']:continue
            if path!=cache['path']:
                # Keep the old file valid for existing sessions; promote only verified data.
                publish(cache,p,record['kept'],record['excluded']);p=Path(cache['path'])
            return dict(coordinates=str(p),kept=record['kept'],excluded=record['excluded'],cached=True,peak_rss_kib=0)
        except (OSError,ValueError,KeyError):continue
    return None


def publish(cache,source,kept,excluded):
    import shutil
    path=Path(cache['path']);path.parent.mkdir(parents=True,exist_ok=True)
    fd,temp=tempfile.mkstemp(prefix='.coordinates-',dir=str(path.parent));os.close(fd)
    try:
        shutil.copyfile(source,temp);checksum=digest(temp);os.replace(temp,str(path))
    finally:
        if os.path.exists(temp):os.unlink(temp)
    manifest=dict(signature=cache['signature'],sha256=checksum,kept=kept,excluded=excluded)
    fd,temp=tempfile.mkstemp(prefix='.manifest-',dir=str(path.parent))
    try:
        with os.fdopen(fd,'w') as f:json.dump(manifest,f,sort_keys=True)
        os.replace(temp,cache['manifest'])
    finally:
        if os.path.exists(temp):os.unlink(temp)
