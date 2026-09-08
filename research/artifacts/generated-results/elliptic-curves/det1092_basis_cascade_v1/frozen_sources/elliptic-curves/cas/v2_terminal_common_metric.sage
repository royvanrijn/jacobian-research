#!/usr/bin/env sage-python
"""Exact common-metric CVP for all 105 omission states and Agent2 M30 axes."""
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import ZZ, QQ, matrix, vector, pari
from research_runtime.store import checkpoint
from terminal_affine_cvp import AffineCVP
CAS=Path(__file__).resolve().parent
b=SourceFileLoader('common_terminal',str(CAS/'diagnose_v2_terminal.sage')).load_module()
D=b.D;plan=b.read(D/'state-plan.json');data=b.read(D/'coordinates.json')
h=matrix(ZZ,plan['common_ambient_metric']);full=matrix(ZZ,plan['B31_ambient_words']);identity=matrix.identity(ZZ,31)
out=D/'common-exact';out.mkdir(exist_ok=True)

def cost(rows,target):
    g=rows*h*rows.transpose();cross=rows*h*target
    u=matrix(ZZ,pari(g).qflllgram()).transpose();reduced=u*g*u.transpose()
    t=vector(QQ,g.solve_right(2*cross))*u.inverse()
    constant=4*(target*h*target)-t*reduced*t
    proof=AffineCVP(reduced.rows()).nearest(t,count=1,node_limit=30000000)
    value=QQ(proof['rows'][0]['distance'])+constant
    assert value.denominator()==1
    return {'minimum_numerator':int(value),'proof':proof,'closest_words':[list(map(int,vector(ZZ,r['word'])*u)) for r in proof['rows']]}

for state in plan['states']:
    path=out/(state['id']+'.json')
    if path.exists():continue
    rows=full[state['kept'],:]
    results={str(i):cost(rows,full.row(i)) for i in state['omitted']}
    for proxy in state['scores']:
        assert results[str(proxy['omitted'])]['minimum_numerator']<=proxy['fresh_common_metric_numerator']
    checkpoint(path,{'state':state,'targets':results})
    print('COMMON EXACT',state['id'],[r['minimum_numerator'] for r in results.values()],flush=True)

for omitted in range(17,31):
    path=out/f'agent2-omit-{omitted}.json'
    if path.exists():continue
    rows=identity[[i for i in range(31) if i!=omitted],:]
    checkpoint(path,{'state_mask':16383^(1<<(omitted-17)),'target_axis':omitted,'result':cost(rows,identity.row(omitted))})
    print('AGENT2 EXACT',omitted,flush=True)

normal=plan['V2_M30_primitive_normal']
core=sum(1<<(i-17) for i in range(17,31) if normal[i]==0)
equivalence={str(i):cost(full[:30,:],identity.row(i)) for i in range(17,31) if normal[i]}
assert len({r['minimum_numerator'] for r in equivalence.values()})==1
checkpoint(out/'equivalent-lifts.json',{'coordinate_core_mask':core,'primitive_normal':normal,'targets':equivalence})
