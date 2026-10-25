"""Memory-admitted, one-task-per-process scheduling (no idle heavy worker pool).

Admission is an estimate, not an OS memory limit. Existing tasks finish when
external memory pressure rises; new work waits. No science parameters change.
"""
from pathlib import Path
import json,math,multiprocessing,os,resource,sys,time,traceback,re,subprocess
GIB=1024**3


def memory_state():
    """Linux MemAvailable, additionally bounded by visible cgroup headroom."""
    values={}
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            key,value=line.split(':',1);values[key]=int(value.split()[0])*1024
        total=values['MemTotal'];available=values.get('MemAvailable',values['MemFree'])
    except (OSError,KeyError,ValueError):
        try:
            total=os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')
            if sys.platform=='darwin':
                # macOS has no SC_AVPHYS_PAGES. Count only free/reclaimable
                # physical pages; exclude compressed/swapped memory.
                report=subprocess.check_output(['/usr/bin/vm_stat'],text=True,timeout=5)
                match=re.search(r'page size of (\d+) bytes',report)
                if match is None:raise ValueError('Cannot read macOS VM page size')
                page_size=int(match.group(1))
                counts={key:int(value) for key,value in re.findall(r'^(Pages [^:]+):\s+(\d+)\.',report,re.M)}
                available=page_size*sum(counts[key] for key in ('Pages free','Pages inactive','Pages speculative'))
            else:
                available=os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_AVPHYS_PAGES')
        except (ValueError,OSError,KeyError,subprocess.SubprocessError):raise RuntimeError('Cannot determine available memory; refusing automatic admission')
    bases={Path('/sys/fs/cgroup')}
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            _,controllers,name=line.split(':',2)
            if not controllers:
                leaf=Path('/sys/fs/cgroup')/name.lstrip('/')
                if '..' not in leaf.parts:
                    bases.update(p for p in (leaf,*leaf.parents) if str(p).startswith('/sys/fs/cgroup'))
            elif 'memory' in controllers.split(','):
                bases.add(Path('/sys/fs/cgroup/memory')/name.lstrip('/'))
    except (OSError,ValueError):pass
    for base in bases:
        for cap,usage in (('memory.max','memory.current'),('memory.limit_in_bytes','memory.usage_in_bytes')):
            try:
                limit=int((base/cap).read_text());used=int((base/usage).read_text())
                if limit<total*100:total=min(total,limit);available=min(available,max(0,limit-used))
            except (OSError,ValueError):pass
    return dict(total=total,available=min(total,available))


def process_rss(pid):
    try:
        for line in (Path('/proc')/str(pid)/'status').read_text().splitlines():
            if line.startswith('VmRSS:'):return int(line.split()[1])*1024
    except (OSError,ValueError):pass
    return 0


def admission_slots(available,reserve,estimate,reservations,active,maximum,budget):
    """Keep headroom for live tasks that have not reached their estimated peak."""
    future=sum(max(0,estimate-rss) for rss in reservations)
    headroom=max(0,available-reserve-future)
    own_slots=max(0,int(budget//estimate)-active)
    return min(maximum-active,own_slots,int(headroom//estimate))


def io_limit(cfg,stage,maximum):
    io_stages=('reference','calibration','baseline','resample','reduce','grid')
    value=cfg.get('stage_io_workers',{}).get(stage,cfg.get('io_workers',4)) if stage in io_stages else maximum
    if isinstance(value,bool) or int(value)!=value or not 1<=int(value)<=64:
        raise ValueError('I/O worker limit must be an integer in 1..64')
    return min(maximum,int(value))


def worker_estimate(cfg,stage,channels=65536,pixels=0,width=0):
    chunk=int(cfg.get('chunk_rows',32))*channels*2*8
    if stage=='reference':return int(cfg.get('reference_worker_gib',1.5)*GIB)
    if stage=='grid':return int(.75*GIB+pixels*width*3*40+int(cfg.get('grid_row_block',256))*width*3*16)
    factors={'calibration':32,'baseline':64,'resample':40,'reduce':72}
    return int(.75*GIB+chunk*factors.get(stage,64))


def safe_block(cfg,requested,pixels,factor=120):
    state=memory_state();reserve=max(2*GIB,int(state['total']*.1))
    budget=min(state['available']-reserve,state['total']*.6,float(cfg.get('memory_budget_gib',state['total']/GIB))*GIB)
    width=min(int(requested),max(0,int((budget-.75*GIB)//max(1,pixels*factor))))
    if width<1:raise MemoryError('Not enough available memory for even one gridding/combination channel; retry after memory is freed')
    return width


def _entry(function,item,initializer,initargs,connection):
    try:
        if initializer:initializer(*initargs)
        result=function(item)
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak=int(peak if sys.platform=='darwin' else peak*1024)
        connection.send(dict(ok=True,result=result,peak=peak))
    except BaseException:
        connection.send(dict(ok=False,error=traceback.format_exc()))
    finally:connection.close()


def memory_map(function,items,cfg,stage,estimate,initializer=None,initargs=(),status_path=None,on_result=None):
    """Dynamically admit tasks, returning results in the original input order."""
    items=list(items)
    if not items:return []
    maximum=int(cfg.get('max_workers',cfg.get('workers',64)));configured_maximum=maximum
    if not 1<=maximum<=64:raise ValueError('Worker maximum must be in 1..64')
    io_maximum=io_limit(cfg,stage,maximum)
    maximum=min(maximum,io_maximum,len(items));estimate=max(1,int(estimate))
    initial=memory_state();budget=min(initial['total']*.6,float(cfg.get('memory_budget_gib',initial['total']/GIB))*GIB)
    reserve=max(2*GIB,int(initial['total']*.1));timeout=float(cfg.get('memory_wait_seconds',300))
    if budget<=0 or timeout<=0:raise ValueError('Memory budget and wait timeout must be positive')
    context=multiprocessing.get_context('spawn');active={};results=[None]*len(items);next_index=0;completed=0;idle_since=None;last_report=0;last_active=-1
    names=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS');previous={k:os.environ.get(k) for k in names}
    for key in names:os.environ[key]='1'
    try:
        while completed<len(items):
            for index,(process,receiver) in list(active.items()):
                if receiver.poll():
                    try:message=receiver.recv()
                    except EOFError:raise RuntimeError('Worker exited without a result for task '+str(index))
                    process.join();receiver.close();del active[index]
                    if not message['ok']:raise RuntimeError(message['error'])
                    if process.exitcode!=0:raise RuntimeError('Worker exited with code '+str(process.exitcode))
                    results[index]=message['result'];completed+=1
                    estimate=max(estimate,int(math.ceil(message['peak']*1.3)))
                    if on_result:on_result(message['result'],completed,len(items))
                elif process.exitcode is not None:
                    raise RuntimeError('Worker exited without a result; exit code '+str(process.exitcode))
            state=memory_state()
            # One pilot task first; its measured peak raises the estimate before scaling up.
            limit=1 if completed==0 else maximum
            slots=admission_slots(state['available'],reserve,estimate,[process_rss(p.pid) for p,_ in active.values()],len(active),limit,budget)
            for _ in range(min(slots,len(items)-next_index)):
                receiver,sender=context.Pipe(duplex=False)
                process=context.Process(target=_entry,args=(function,items[next_index],initializer,initargs,sender))
                process.start();sender.close();active[next_index]=(process,receiver);next_index+=1
            now=time.monotonic()
            if status_path and (now-last_report>=2 or len(active)!=last_active or completed==len(items)):
                p=Path(status_path);p.parent.mkdir(parents=True,exist_ok=True);temp=p.with_suffix('.tmp')
                temp.write_text(json.dumps(dict(stage=stage,configured_maximum_workers=configured_maximum,maximum_workers=maximum,io_worker_limit=io_maximum,active_workers=len(active),completed=completed,total=len(items),queued=len(items)-next_index,estimated_worker_gib=estimate/GIB,available_gib=state['available']/GIB,reserve_gib=reserve/GIB,budget_gib=budget/GIB,state='complete' if completed==len(items) else 'running' if active else 'waiting_for_memory'),indent=2)+'\n');temp.replace(p)
                if len(active)!=last_active:
                    print('SCHEDULER {}: active {}/{}; I/O cap {}; complete {}/{}; estimate {:.2f} GiB/worker; available {:.1f} GiB'.format(stage,len(active),configured_maximum,io_maximum,completed,len(items),estimate/GIB,state['available']/GIB),flush=True)
                last_report=now;last_active=len(active)
            if completed==len(items):break
            if not active:
                if idle_since is None:idle_since=now;print('SCHEDULER waiting for memory:',stage,flush=True)
                if now-idle_since>=timeout:raise MemoryError('No task admitted for '+str(timeout)+' seconds; estimated worker '+str(round(estimate/GIB,2))+' GiB. Free memory or lower chunk/block sizes.')
            else:idle_since=None
            time.sleep(.1 if active else 1.)
        return results
    finally:
        for process,receiver in active.values():
            if process.is_alive():process.terminate()
        for process,receiver in active.values():process.join();receiver.close()
        for key,value in previous.items():
            if value is None:os.environ.pop(key,None)
            else:os.environ[key]=value
