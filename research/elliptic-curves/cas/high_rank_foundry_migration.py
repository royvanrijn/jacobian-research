"""Freeze and restore a stopped foundry without resetting its scientific history."""
import os
from pathlib import Path
import shutil


def freeze(parent, frozen):
    import run_high_rank_foundry as c
    parent=Path(parent).resolve()
    locks=[]
    try:
        for name in ('launch','guardian','controller'):
            locks.append(c.exclusive(parent/(name+'.lock')))
        oldroot=c.guard(parent);config=c.read(parent/'config.json')
        db=c.connect(parent)
        try:
            for role in ('guardian','controller'):
                p=c.meta(db,role,{}) or {}
                c.require(not c.same_process(p.get('pid'),p.get('token')),'parent is still running')
            cc,jj=c.curves(db),c.jobs(db)
            c.require(not any(j['state'] in ('RUNNING','PENDING') for j in jj),'parent jobs must drain before migration')
            saved_meta={k:c.meta(db,k,default) for k,default in
                        (('intake_cursor',{}),('fresh_completed',0),('last_revival',0))}
            old_status=c.meta(db,'status')
        finally:db.close()
        data=frozen/'foundry-inputs'
        for name in ('intake','catalogue'):
            shutil.copyfile(oldroot/config[name],data/(name+'.json'))
        # Canonical proof packets are copied byte-for-byte. Cached landscapes bind
        # old source hashes, so retire only those caches and select a fresh bank.
        transformations=[]
        for row in cc:
            before={k:row.get(k) for k in ('packet','head','bank_index','rank')}
            if row.get('packet'):
                src=oldroot/row['packet']
                c.require(c.sha(src)==row['packet_sha256'],'parent certificate changed')
                dest=frozen/'foundry-inherited'/f"{row['id']}.json"
                dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
                row['packet']=str(dest.relative_to(frozen))
            else:
                # An input placeholder is not a certified rank17 subgroup.
                row['rank']=None
            if row.get('head'):
                row['previous_runtime_head']=str(oldroot/row['head'])
                row['head']=None;row['bank_index']+=1
            transformations.append({'id':row['id'],'before':before,
                                    'after':{k:row.get(k) for k in before}})
        for j in jj:
            j['path']=str(oldroot/j['path'])
            if j.get('export'):
                c.require(c.sha(Path(j['export'])),'missing parent result export')
        state={'schema':'foundry-inherited-state.v1','parent':str(parent),'parent_status':old_status,
               'parent_manifest_sha256':c.sha(parent/'manifest.json'),
               'parent_config_sha256':c.sha(parent/'config.json'),
               'meta':saved_meta,'curves':cc,'jobs':jj,'transformations':transformations,
               'parent_exports':{j['export']:c.sha(Path(j['export'])) for j in jj if j.get('export')},
               'claim_boundary':'Source upgrade preserves equations, proofs, trajectories, cooling, intake cursor and charged work. Previous raw logs remain in the stopped parent. Source-bound landscape caches are replaced by a new parent bank.'}
        c.save(data/'inherited-state.json',state,immutable=True)
        c.save(data/'warm.json',{'rows':[],'source':str(parent),'inherited':True},immutable=True)
        return config
    finally:
        for fd in reversed(locks):os.close(fd)


def restore(frozen,db):
    import run_high_rank_foundry as c
    import high_rank_foundry_arithmetic as arithmetic
    state=c.read(frozen/'foundry-inputs/inherited-state.json')
    checked=0
    for row in state['curves']:
        if row.get('packet'):
            path=frozen/row['packet']
            c.require(c.sha(path)==row['packet_sha256'],'inherited packet changed')
            arithmetic.verify_packet(c.read(path),row);checked+=1
    for name,digest in state['parent_exports'].items():
        c.require(c.sha(Path(name))==digest,'parent export changed')
    with db:
        for row in state['curves']:c.putcurve(db,row)
        for j in state['jobs']:c.putjob(db,j)
        for key,value in state['meta'].items():c.putmeta(db,key,value)
        c.putmeta(db,'inherited_from',state['parent'])
        c.putmeta(db,'consecutive_failures',0)
    return checked
