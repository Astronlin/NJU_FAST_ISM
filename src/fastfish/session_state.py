"""Atomic administrative session metadata shared by the FISH prompt and plots."""
from pathlib import Path
import fcntl,json,os,tempfile

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
