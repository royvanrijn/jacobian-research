#!/usr/bin/env sage-python
"""Bounded complete generic17 centre-class geometry; no point search or oracle."""
import sys,json,hashlib,gzip,heapq
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
D=ROOT/'artifacts/local/elliptic-curves/curve302-complete-generic-parity-v1'
SEED=ROOT/'artifacts/local/elliptic-curves/curve302-focused-point-exposure-v2/curve302-generic17/seed.json'
READS=set()

def guard(event,args):
    if event!='open' or not isinstance(args[0],(str,bytes)):return
    path=Path(args[0]).resolve()
    if path.is_relative_to(ROOT/'artifacts'):
        if not (path.is_relative_to(D) or path==SEED):raise PermissionError('generic-only artifact guard: '+str(path))
        READS.add(str(path.relative_to(ROOT)))

sys.addaudithook(guard)
sys.path.insert(0,str(CAS))
from sage.all import matrix,ZZ,pari
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
geometry=SourceFileLoader('complete_generic_geometry',str(CAS/'prospective_half_lattice_v3.sage')).load_module()
mapper=SourceFileLoader('complete_generic_mapper',str(CAS/'factor_free_pari_mapping.sage')).load_module()

def main():
    p=cert.read(D/'protocol.json')
    assert all(cert.hashed(ROOT/n)==h for n,h in p['sources'].items())
    assert cert.hashed(SEED)==p['seed_sha256']
    assert p['classes']==131071 and p['charts']==49
    assert not (D/'maps.json').exists()
    seed=cert.read(SEED);assert len(seed['points'])==17
    model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    gram,asym=geometry.canonical_height_gram(model,points)
    g=matrix(ZZ,geometry.rounded_gram(gram,1000000));u=matrix(ZZ,pari(g).qflllgram()).transpose()
    assert abs(u.det())==1
    inverse=u.inverse();reduced=u*g*u.transpose();oracle=geometry.CosetOracle(reduced.rows())
    heap=[];chunk=[];blocks=[]
    checkpoint(D/'metric.json',dict(metric_gram=[[str(v) for v in row] for row in gram],
        maximum_gram_asymmetry=str(asym),rounded_gram=[list(map(int,r)) for r in g.rows()],
        change_of_basis=[list(map(int,r)) for r in u.rows()],reduced_gram=[list(map(int,r)) for r in reduced.rows()]))
    for mask in range(1,1<<17):
        residue=matrix(ZZ,1,17,[(mask>>j)&1 for j in range(17)])
        target=[int(v)%2 for v in (residue*inverse).row(0)]
        norm,rep,error=oracle.solve(target)
        word=list(map(int,(matrix(ZZ,1,17,rep)*u).row(0)))
        assert all((word[j]-(mask>>j))%2==0 for j in range(17))
        assert sum(word[j]*g[j,k]*word[k] for j in range(17) for k in range(17))==norm
        entry=dict(parity=mask,representative=word,metric_norm=int(norm),
                   cvp_error=error,reduced_representative=list(map(int,rep)))
        chunk.append(entry);key=(int(norm),-mask)
        if len(heap)<49:heapq.heappush(heap,(key,entry))
        elif key>heap[0][0]:heapq.heapreplace(heap,(key,entry))
        if len(chunk)==4096 or mask==131071:
            path=D/'parity-blocks'/('block-'+str(len(blocks)).zfill(3)+'.json.gz')
            path.parent.mkdir(exist_ok=True)
            assert not path.exists()
            path.write_bytes(gzip.compress(json.dumps(chunk,separators=(',',':')).encode(),mtime=0))
            blocks.append(dict(path=str(path.relative_to(ROOT)),sha256=cert.hashed(path),rows=len(chunk)))
            chunk=[]
            checkpoint(D/'progress.json',dict(status='RUNNING_GEOMETRY',completed_classes=mask,blocks=blocks))
            print('GENERIC PARITIES',mask,'/131071',flush=True)
    centres=[entry for key,entry in sorted(heap,reverse=True)]
    data=dict(status='RUNNING_MAPS',protocol_hash=digest(p),classes=131071,blocks=blocks,
              centres=centres,rows=[],**cert.read(D/'metric.json'))
    checkpoint(D/'maps.json',data)
    pari.allocatemem(256000000,silent=True)
    for centre in centres:
        data['rows'].append(mapper.mapping(model,points,centre));checkpoint(D/'maps.json',data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(D/'maps.json',data)
    checkpoint(D/'data-access.json',sorted(READS))
    print('PASS131071 generic parity classes,49 maps; no point search',flush=True)

if __name__=='__main__':main()
