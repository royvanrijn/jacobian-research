#!/usr/bin/env python3
"""Independent Fraction replay of the complete orthogonal-block support bound."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import ast,json,hashlib
from verify import add,mul,scale,variable,compose,decode,identity,mmul
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def parse(text,n,names):
 one={(0,)*n:Q(1)}
 def rec(a):
  if isinstance(a,ast.Constant) and isinstance(a.value,int):return scale(one,Q(a.value))
  if isinstance(a,ast.Name):return variable(names.index(a.id),n)
  if isinstance(a,ast.UnaryOp):return scale(rec(a.operand),-1 if isinstance(a.op,ast.USub) else 1)
  if isinstance(a,ast.BinOp):
   l=rec(a.left)
   if isinstance(a.op,ast.Pow):
    assert isinstance(a.right,ast.Constant) and isinstance(a.right.value,int) and a.right.value>=0
    out=one
    for _ in range(a.right.value):out=mul(out,l)
    return out
   r=rec(a.right)
   if isinstance(a.op,ast.Add):return add(l,r)
   if isinstance(a.op,ast.Sub):return add(l,r,Q(-1))
   if isinstance(a.op,ast.Mult):return mul(l,r)
   if isinstance(a.op,ast.Div):
    assert len(r)==1 and (0,)*n in r and r[(0,)*n]
    return scale(l,1/r[(0,)*n])
  raise AssertionError(ast.dump(a))
 return rec(ast.parse(text,mode='eval').body)
def exact_divide(p,f):
 lf=max(f);cf=f[lf];out={};rem=dict(p)
 while rem:
  lm=max(rem)
  if any(a<b for a,b in zip(lm,lf)):return None
  ex=tuple(a-b for a,b in zip(lm,lf));co=rem[lm]/cf;term={ex:co}
  out=add(out,term);rem=add(rem,mul(term,f),Q(-1))
 return out
def rep(a,b,c,d,sgn,n):
 L=identity(12,n)
 L[10]=[{} for _ in range(12)];L[10][10]=mul(a,a);L[10][1]=scale(mul(a,b),2);L[10][6]=scale(mul(b,b),Q(1,2))
 L[1]=[{} for _ in range(12)];L[1][10]=mul(a,c);L[1][1]=add(mul(a,d),mul(b,c));L[1][6]=scale(mul(b,d),Q(1,2))
 L[6]=[{} for _ in range(12)];L[6][10]=scale(mul(c,c),2);L[6][1]=scale(mul(c,d),4);L[6][6]=mul(d,d)
 for i in (1,6,10):L[i]=[scale(p,sgn) for p in L[i]]
 eye=identity(12,n);L[9]=[add(add(eye[9][j],eye[1][j]),L[1][j],Q(-1)) for j in range(12)]
 return L
def mvec(L,V):return [sum_poly(mul(p,v) for p,v in zip(row,V)) for row in L]
def sum_poly(items):
 out={}
 for p in items:out=add(out,p)
 return out
def audit(name,sgn,boundary):
 path=HERE/'orthogonal-surface12'/name/'model.json';model=json.loads(path.read_text());params=model['parameters'];np=len(params);n=12+np
 X=[variable(i,n) for i in range(n)];one={(0,)*n:Q(1)}
 if boundary:a,b,c,d=X[12],one,scale(one,-1),{}
 else:a,b,c,d=add(one,mul(X[12],X[13])),X[12],X[13],one
 L=rep(a,b,c,d,sgn,n);Li=rep(d,scale(b,-1),scale(c,-1),a,sgn,n)
 assert mmul(L,Li)==identity(12,n) and mmul(Li,L)==identity(12,n)
 forms=mvec(L,X[:12])
 assert add(scale(mul(forms[1],forms[1]),2),mul(forms[10],forms[6]),Q(-1))==add(scale(mul(X[1],X[1]),2),mul(X[10],X[6]),Q(-1))
 assert add(forms[1],forms[9])==add(X[1],X[9])
 raw=json.loads((HERE/'cubic12-integer-witness.json').read_text());F=[decode(row,12) for row in raw['coefficients']]
 FH=[compose(p,forms,n) for p in F];G=mvec(Li,FH)
 groups=Counter();factors=[parse(f,np,params) for f in model['factors']]
 for i,row in enumerate(G):
  grouped={}
  for ex,co in row.items():grouped.setdefault(ex[:12],{})[ex[12:]]=co
  primary={tuple(v['exponent']):parse(v['coefficient'],np,params) for v in model['coefficients'][i]}
  assert grouped==primary
  assert grouped.get(tuple(int(j==i) for j in range(12)))=={(0,)*np:Q(1)}
  for ex,coef in grouped.items():
   assert sum(ex)<=3
   if sum(ex)<2:continue
   labels=[];rem=coef
   for j,factor in enumerate(factors):
    used=False
    while True:
     quotient=exact_divide(rem,factor)
     if quotient is None:break
     rem=quotient;used=True
    if used:labels.append(j)
   assert len(rem)==1 and (0,)*np in rem and rem[(0,)*np]
   groups[tuple(labels)]+=1
 expected=Counter({tuple(row['factors']):row['multiplicity'] for row in model['factor_groups']});assert groups==expected
 assert model['source_sha256']==sha(HERE/'orthogonal_surface12.py') and model['input_sha256']==sha(HERE/'cubic12-integer-witness.json')
 result={'generic_support':sum(groups.values()),'constant_occurrences':groups[()], 'model_sha256':sha(path)}
 def mandatory(known):return sum(v for group,v in groups.items() if set(group)<=set(known))
 if boundary:
  nonzero=[j for j,factor in enumerate(factors) if len(factor)==1] # scalar*a^k.
  result['a_nonzero_lower_bound']=mandatory(nonzero)
  known0=[j for j,f in enumerate(factors) if f.get((0,),0)]
  result['a_zero_count']=mandatory(known0)
  if sgn==1:assert groups[()]>=38 and result['a_nonzero_lower_bound']>=47 and result['a_zero_count']==47
  else:assert result['a_nonzero_lower_bound']>=50 and result['a_zero_count']==36
 else:
  lookup=lambda text:next(j for j,f in enumerate(factors) if f==parse(text,np,params))
  it,iu,iv,iw=[lookup(t) for t in ['t','u','t*u+1','t*u+1/2']]
  bounds={'t_u_nonzero_v_w_nonzero':mandatory([it,iu,iv,iw]),'t_u_nonzero_v_zero':mandatory([it,iu,iw]),'t_u_nonzero_w_zero':mandatory([it,iu,iv]),'t_zero_u_nonzero':mandatory([iu,iv,iw])}
  known=[]
  for j,f in enumerate(factors):
   ff={e:c for e,c in f.items() if e[1]==0}
   if len(ff)==1 and next(iter(ff.values())):known.append(j)
  bounds['u_zero_t_nonzero']=mandatory(known)
  known0=[j for j,f in enumerate(factors) if f.get((0,0),0)];origin=mandatory(known0)
  assert all(v>=47 for v in bounds.values()),bounds
  assert origin==(36 if sgn==1 else 47),origin
  result.update(stratum_lower_bounds=bounds,origin_count=origin)
 return result

def main():
 results={name:audit(name,sgn,boundary) for name,sgn,boundary in [('open',1,False),('boundary',1,True),('open-minus',-1,False),('boundary-minus',-1,True)]}
 old=json.loads((HERE/'integer-validation12.json').read_text());assert old['witness_sha256']==sha(HERE/'cubic12-integer-witness.json') and old['nonlinear_total']==36
 receipt={'status':'pass','claim':'KCUB24','minimum_in_defined_orthogonal_block_group':36,'next_possible_nonlinear_count':47,'no_counts37_through46_in_this_group':True,'proof_sha256':sha(HERE/'ORTHOGONAL_BLOCK_MINIMUM.md'),'minimizing_chart_points':{'open':['0','0'],'boundary-minus':['0']},'all_four_charts_reconstructed':True,'independent_arithmetic':'stdlib Fraction and guarded AST; no CAS or primary factorization used','chart_results':results,'checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'producer_sha256':sha(HERE/'orthogonal_surface12.py'),'input_sha256':sha(HERE/'cubic12-integer-witness.json'),'base_validation_sha256':sha(HERE/'integer-validation12.json'),'proof_boundary':'Chart completeness and the conic description of the block orthogonal group are written. Exact full polynomial conjugacy, inverse matrices, coefficients, factor products and all support indicators are independently reconstructed. Not full affine-orbit, arbitrary nonlinear or global sparsity optimum.','positive_bounds_unchanged':True}
 (HERE/'orthogonal-surface12/validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
