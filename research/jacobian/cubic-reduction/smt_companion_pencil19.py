#!/usr/bin/env python3
"""Exact LRA/cardinality queries on a changed companion pencil's finite pool."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,sys
import z3
HERE=Path(__file__).resolve().parent;OUT=HERE/'companion-pencil19'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ip=OUT/'input.json';bp=OUT/'basis.json';d=json.loads(ip.read_text());bd=json.loads(bp.read_text());assert bd['input_sha256']==sha(ip);A=[{int(k):Q(c) for k,c in row.items()} for row in bd['A']];B=[{int(k):Q(c) for k,c in row.items()} for row in bd['B']];basis=bd['basis'];m=len(d['forms']);nonbasis=[j for j in range(m) if j not in basis];positions={b:p for p,b in enumerate(basis)};sigma=z3.Real('sigma');P={(q,j):z3.Real('n_%d_%d'%(q,j)) for q in nonbasis for j in range(18)};active=[z3.Bool('active_'+str(j)) for j in range(m)];solver=z3.Solver();solver.set(timeout=20000);solver.add(sigma!=0,active[d['padding_index']]);expressions=[];started=time.monotonic()
 for q in range(m):
  if q in positions:
   p=positions[q];values=[z3.RealVal(str(B[p].get(j,Q(0))))+sigma*z3.RealVal(str(B[p].get(j+19,Q(0))))-z3.Sum([z3.RealVal(str(c))*P[v,j] for v,c in A[p].items() if v in nonbasis]) for j in range(18)]
  else:values=[P[q,j] for j in range(18)]
  expressions.append(values);solver.add(z3.Implies(z3.Not(active[q]),z3.And([value==0 for value in values])))
 tests=[]
 for bound in [88,97]:
  solver.push();solver.add(z3.PbLe([(a,1) for a in active],bound));t=time.monotonic();answer=solver.check();row=dict(bound=bound,status=str(answer),solver_seconds=round(time.monotonic()-t,3));print(row,flush=True)
  if answer==z3.sat:
   model=solver.model();sv=model.eval(sigma,model_completion=True);assert z3.is_rational_value(sv) and sv.as_fraction()!=0;rows=[]
   for i,values in enumerate(expressions):
    coeff=[]
    for value in values:
     v=model.eval(value,model_completion=True);assert z3.is_rational_value(v);coeff.append(str(v))
    coeff.append(str(int(i==d['padding_index'])));rows.append(coeff)
   selected=[i for i,row0 in enumerate(rows) if any(Q(v) for v in row0)];row.update(sigma=str(sv),size=len(selected));assert len(selected)<=bound;candidate=dict(status='solver-exact-candidate-needs-replay',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),sigma=str(sv),size=len(selected),forms=[d['forms'][i] for i in selected],B=[[rows[i][j] for i in selected] for j in range(19)],padding_index=selected.index(d['padding_index']),scope='Exact rational LRA assignment, not accepted until independent coefficient/Keller/collision witness verification.');(OUT/('smt-candidate-'+str(bound)+'.json')).write_text(json.dumps(candidate,indent=2,sort_keys=True)+'\n')
  elif answer==z3.unknown:row['reason']=solver.reason_unknown()
  solver.pop();tests.append(row);(OUT/'smt-checkpoint.json').write_text(json.dumps(dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),basis_sha256=sha(bp),tests=tests,solver_version=z3.get_version_string()),indent=2,sort_keys=True)+'\n')
  if answer==z3.sat and bound<89:break
 receipt=dict(status='complete-with-unknown' if any(t['status']=='unknown' for t in tests) else 'complete',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),basis_sha256=sha(bp),candidate_forms=m,cube_rank=d['rank'],free_kernel_columns=len(nonbasis),timeout_ms_per_query=20000,solver_version=z3.get_version_string(),seconds=round(time.monotonic()-started,3),tests=tests,scope='Finite277-form changed-tensor pencil with one common nonzero real parameter, mandatory pure-tau cube and normalized last coefficient row, and arbitrary real coefficients in the other18 outputs. SAT is rational but needs full replay; UNKNOWN is not UNSAT, and solver UNSAT alone is not an independent/global lower bound.');(OUT/'smt-checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('DONE',receipt['status'],flush=True)
if __name__=='__main__':main()
