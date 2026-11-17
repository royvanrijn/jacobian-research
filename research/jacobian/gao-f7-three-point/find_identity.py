"""Weighted Macaulay certificate search with exact rational reconstruction.

Extends the retained Case-2 Macaulay method to weighted monomials. Modular
pivots select a small system only; a candidate is accepted solely after an
exact polynomial identity check. No finite-field absence claim is made.
"""
import hashlib
import json
from pathlib import Path
import sys
import time
import resource
import sympy as S
from flint import nmod_mat, fmpq_mat, fmpq_mpoly_ctx, fmpq

HERE=Path(__file__).resolve().parent
case=sys.argv[1] if len(sys.argv)>1 else 'cubic_cube'
data=json.loads((HERE/(case+'-input.json')).read_text())
variables=data['variables']
weights={'cubic_cube':[1,2,3,1,2],
         'linear_cube_cubic_square':[1,1,2,3,1,2],
         'cubic-delta-zero':[1,2,1,3],
         'cubic-delta-nonzero':[1,2,1,2],
         'cubic-raw-delta':[1,2,3,1,2],
         'cubic-raw-delta-batch':[1,2,3,1,2]}[case]
memory_mb=3072 if case=='cubic-raw-delta-batch' else 1536
resource.setrlimit(resource.RLIMIT_AS,(memory_mb*1024**2,memory_mb*1024**2))
target_exprs=[S.sympify(x) for x in data.get('targets',[data.get('target','d**3' if case=='cubic-delta-nonzero' else '1')])]
all_targets=[{m:int(c) for m,c in S.Poly(expr,*list(map(S.Symbol,variables)),domain=S.ZZ).terms()}
             for expr in target_exprs]
symbols=list(map(S.Symbol,variables))
generators=[{m:int(c) for m,c in S.Poly(S.sympify(row['polynomial']),*symbols,domain=S.ZZ).terms()}
            for row in data['equations']]
degrees=[max(sum(a*b for a,b in zip(m,weights)) for m in f) for f in generators]
prime=32003
started=time.monotonic()
rounds=[]


def monomials(bound,i=0):
    if i==len(weights):
        yield ()
        return
    if bound<0:return
    for a in range(bound//weights[i]+1):
        for rest in monomials(bound-a*weights[i],i+1):yield (a,)+rest


def pivot_columns(matrix,rank):
    out=[]
    for i in range(rank):
        out.append(next(j for j in range(matrix.ncols()) if matrix[i,j]))
    return out


for bound in ([int(sys.argv[2])] if len(sys.argv)>2 else [12,14,16,18,20,22]):
    columns=[(i,m) for i,deg in enumerate(degrees) for m in monomials(bound-deg)]
    sparse=[];rowset=set().union(*all_targets)
    for i,m in columns:
        col={tuple(a+b for a,b in zip(n,m)):c for n,c in generators[i].items()}
        sparse.append(col);rowset.update(col)
    rows=sorted(rowset,key=lambda m:(sum(a*b for a,b in zip(m,weights)),m))
    indices={m:i for i,m in enumerate(rows)}
    report={'weighted_bound':bound,'rows':len(rows),'columns':len(columns),'seconds_at_start':time.monotonic()-started}
    if len(rows)*(len(columns)+len(all_targets))>(45_000_000 if case=='cubic-raw-delta-batch' else 15_000_000):
        report['outcome']='matrix-entry cap';rounds.append(report);break
    matrix=nmod_mat(len(rows),len(columns)+len(all_targets),prime)
    for j,col in enumerate(sparse):
        for m,c in col.items():matrix[indices[m],j]=c%prime
    for k,target_terms in enumerate(all_targets):
        for m,c in target_terms.items():matrix[indices[m],len(columns)+k]=c%prime
    print(json.dumps(report),flush=True)
    reduced,rank=matrix.rref(inplace=True)
    pivots=pivot_columns(reduced,rank)
    base_pivots=[j for j in pivots if j<len(columns)]
    base_rank=len(base_pivots)
    consistent=[k for k in range(len(all_targets))
                if all(reduced[i,len(columns)+k]==0 for i in range(base_rank,rank))]
    report['consistent_targets']=[str(target_exprs[k]) for k in consistent]
    if not consistent:
        report['outcome']='no identity at this bound modulo prime';rounds.append(report)
        print(json.dumps(report),flush=True)
        del matrix,reduced
        continue
    chosen=min(consistent,key=lambda k:sum(bool(reduced[i,len(columns)+k]) for i in range(base_rank)))
    target_expr=target_exprs[chosen];target_terms=all_targets[chosen]
    support=[base_pivots[i] for i in range(base_rank) if reduced[i,len(columns)+chosen]]
    report['modular_rank']=base_rank;report['support_columns']=len(support)
    report['chosen_target']=str(target_expr)
    if len(support)>2500:
        report['outcome']='exact solve dimension cap';rounds.append(report);break
    # The support columns are independent. Choose independent original
    # rows by reducing their transpose, then solve that square system in Q.
    del matrix,reduced
    transposed=nmod_mat(len(support),len(rows),prime)
    for k,j in enumerate(support):
        for m,c in sparse[j].items():transposed[k,indices[m]]=c%prime
    rr,rk=transposed.rref(inplace=True);selected=pivot_columns(rr,rk)
    assert rk==len(support)
    del transposed,rr
    (HERE/(case+'-macaulay-checkpoint.json')).write_text(json.dumps({'weighted_bound':bound,
        'support_columns':support,'selected_rows':selected,'target':str(target_expr),
        'outcome':'exact rational solve not yet certified'})+'\n')
    square=fmpq_mat([[sparse[j].get(rows[i],0) for j in support] for i in selected])
    rhs=fmpq_mat([[target_terms.get(rows[i],0)] for i in selected])
    print('exact solve '+str(len(support)),flush=True)
    solution=square.solve(rhs)
    context=fmpq_mpoly_ctx.get(variables)
    multipliers=[{} for f in generators]
    for k,j in enumerate(support):
        if solution[k,0]:
            i,m=columns[j];multipliers[i][m]=solution[k,0]
    total=context.constant(0)
    for gen,mul in zip(generators,multipliers):
        total+=context.from_dict(mul)*context.from_dict({m:fmpq(c) for m,c in gen.items()})
    expected=context.from_dict({m:fmpq(c) for m,c in target_terms.items()})
    if total!=expected:
        report['outcome']='modular support did not reconstruct an exact identity'
        rounds.append(report);print(json.dumps(report),flush=True);break
    cert={'variables':variables,'source':case+'-input.json','target':str(target_expr),'weighted_bound':bound,
          'multipliers':[[{'exponents':list(m),'coefficient':str(c)} for m,c in mul.items()]
                         for mul in multipliers],
          'generator_indices':list(range(len(generators)))}
    (HERE/(case+'-macaulay-certificate.json')).write_text(json.dumps(cert,separators=(',',':'))+'\n')
    report['outcome']='EXACT IDENTITY PASS';report['multiplier_terms']=sum(len(m) for m in multipliers)
    rounds.append(report);print(json.dumps(report),flush=True);break

(HERE/(case+'-macaulay-search.json')).write_text(json.dumps({'prime':prime,'weights':weights,'rounds':rounds,
    'seconds':time.monotonic()-started,'source_sha256':hashlib.sha256((HERE/(case+'-input.json')).read_bytes()).hexdigest(),
    'boundary':'Only the exported exact rational polynomial identity is a characteristic-zero certificate; all failed finite bounds remain inconclusive.'},indent=2)+'\n')
