#!/usr/bin/env python3
"""Independent Fraction audit of the collision-free 23-atom weighted class.

Written curve/pole and polynomial-degree arguments accompany the finite
identities. This is not proof-assistant formalization.
"""
from pathlib import Path
from fractions import Fraction as Q
import ast,hashlib,json,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode,scale
from weighted_verify_utils import substitution,expression,equal_generators,coefficient_buckets,power,normal
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def laurent(text,names):
 n=len(names);one=const(1,n)
 def visit(x):
  if isinstance(x,ast.Constant):assert type(x.value) is int;return const(x.value,n)
  if isinstance(x,ast.Name):return variable(names.index(x.id),n)
  if isinstance(x,ast.UnaryOp):return scale(visit(x.operand),Q(-1) if isinstance(x.op,ast.USub) else Q(1))
  assert isinstance(x,ast.BinOp)
  if isinstance(x.op,ast.Add):return add(visit(x.left),visit(x.right))
  if isinstance(x.op,ast.Sub):return add(visit(x.left),visit(x.right),Q(-1))
  if isinstance(x.op,ast.Mult):return mul(visit(x.left),visit(x.right))
  if isinstance(x.op,ast.Div):
   a,b=visit(x.left),visit(x.right);assert len(b)==1;e,c=next(iter(b.items()));return mul(a,{tuple(-k for k in e):1/c})
  assert isinstance(x.op,ast.Pow) and isinstance(x.right,ast.Constant) and type(x.right.value) is int and x.right.value>=0;return power(visit(x.left),x.right.value,n)
 return visit(ast.parse(text.replace('^','**'),mode='eval').body)
def clear_monomial_denominator(poly,n):
 if not poly:return {}
 shift=tuple(-min(0,min(e[i] for e in poly)) for i in range(n));return {tuple(a+b for a,b in zip(e,shift)):c for e,c in poly.items()}
def linear_audit(linear,oldnames,newnames,claimed,allowed):
 # Exact row operations over Laurent polynomials. Every pivot is a single
 # monomial in parameters explicitly assumed nonzero; no exceptional roots
 # are hidden by rational-function Gaussian elimination.
 nn=len(newnames);cols=[name for name in oldnames if name not in newnames or name in claimed];alltails=[name for name in oldnames if name[0] in ['g','h']];params=[name for name in oldnames if name not in alltails];rows=[]
 for f in linear:
  row=[{} for _ in range(len(alltails)+1)]
  for e,c in f.items():
   tail=[(j,e[oldnames.index(name)]) for j,name in enumerate(alltails) if e[oldnames.index(name)]];assert len(tail)<=1 and (not tail or tail[0][1]==1);idx=tail[0][0] if tail else len(alltails);pe=[0]*nn
   for name in params:pe[newnames.index(name)]=e[oldnames.index(name)]
   row[idx]=add(row[idx],{tuple(pe):c})
  rows.append(row)
 pivots=[];r=0;remaining=[(j,name) for j,name in enumerate(alltails) if name in claimed]
 def invertible(p):
  if len(p)!=1:return False
  e=next(iter(p));return all(k==0 or newnames[i] in allowed for i,k in enumerate(e))
 while remaining:
  choice=next(((k,j,name) for j,name in remaining for k in range(r,len(rows)) if invertible(rows[k][j])),None);assert choice is not None,('no allowed monomial pivot',[name for j,name in remaining])
  k,j,name=choice;remaining.remove((j,name));rows[r],rows[k]=rows[k],rows[r];e,c=next(iter(rows[r][j].items()));inv={tuple(-a for a in e):1/c};rows[r]=[mul(p,inv) for p in rows[r]]
  for k in range(len(rows)):
   if k!=r and rows[k][j]:
    factor=rows[k][j];rows[k]=[add(a,mul(factor,b),Q(-1)) for a,b in zip(rows[k],rows[r])]
  pivots.append((name,j));r+=1
 for k,(name,j) in enumerate(pivots):
  value=scale(rows[k][-1],Q(-1))
  for i,other in enumerate(alltails):
   if other not in claimed:value=add(value,mul(rows[k][i],variable(newnames.index(other),nn)),Q(-1))
  assert value==laurent(claimed[name],newnames),(name,'parameterization mismatch')
 assert all(not p for row in rows[r:] for p in row);return r

def main():
 started=time.monotonic();root=HERE/'weighted-unmarked7';source=root/'input.json';d=load(source);names=d['variables'];n=len(names);assert n==21 and d['source_sha256']==sha(HERE/'weighted_unmarked7.py') and d['input_sha256']==sha(HERE/'weighted-newfirst7/input.json');atoms=d['atoms'];expected=set()
 for row,target in enumerate([1,-1,-2]):
  for a in range(8):
   for b in range(8):
    for c in range(8):
     if 2<=a+b+c<=7 and a-b-2*c==target:expected.add((row,(a,b,c)))
 assert len(atoms)==23 and {(a['row'],tuple(a['exponent'])) for a in atoms}==expected
 N=3+n;F=[variable(i,N) for i in range(3)]
 for atom in atoms:
  e=list(atom['exponent'])+[0]*n
  if atom['coefficient'] in names:e[3+names.index(atom['coefficient'])]=1;c=Q(1)
  else:c=Q(atom['coefficient']);assert atom['row']==0 and c in [Q(-3),Q(-1)]
  F[atom['row']]=add(F[atom['row']],{tuple(e):c})
 det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),const(1,N),Q(-1));buckets=coefficient_buckets(det,3);prepared=[decode(p,n) for p in d['equations']];equal_generators(list(buckets.values()),prepared);assert all(a==b+2*c for a,b,c in buckets)
 qnames=['u','w']+names;nn=2+n;u,w=variable(0,nn),variable(1,nn);qgh=[expression('1-3*u-w',qnames),u,w]
 for atom in atoms:
  if 'quotient_exponent' not in atom:continue
  a,b=atom['quotient_exponent'];qgh[atom['row']]=add(qgh[atom['row']],mul(variable(qnames.index(atom['coefficient']),nn),mul(power(u,a,nn),power(w,b,nn))))
 f,g,h=qgh
 for p,text in zip(qgh,d['quotients']):assert p==expression(text,qnames)
 def jac(f,g,h,width):
  fu,fw=derivative(f,0),derivative(f,1);gu,gw=derivative(g,0),derivative(g,1);hu,hw=derivative(h,0),derivative(h,1)
  return add(add(mul(g,add(mul(fu,hw),mul(fw,hu),Q(-1))),scale(mul(h,add(mul(gu,fw),mul(gw,fu),Q(-1))),Q(2))),mul(f,add(mul(gu,hw),mul(gw,hu),Q(-1))))
 qj=add(jac(f,g,h,nn),const(1,nn),Q(-1));assert qj=={(b,c)+e:co for (a,b,c),poly in buckets.items() for e,co in poly.items()}
 assert f==expression('1-3*u-w+f1*u**3+f2*u**2+f3*u*w+f4*w**2',qnames)
 tails=[name for name in names if name[0]!='f'];branch_cache={};boundary_cache={};cases=[]
 def remap(newnames,mode):
  z={name:variable(i,len(newnames)) for i,name in enumerate(newnames)};zero={};R=z.get('f1',zero);P=z.get('f2',zero);Qp=z.get('f3',zero);T=z.get('f4',zero)
  if mode=='graph':vals=[R,P,Qp,zero]
  elif mode=='conic':vals=[zero,P,Qp,T]
  elif mode=='cubicgraph':vals=[R,P,zero,zero]
  elif mode=='quadraticgraph':vals=[zero,P,zero,zero]
  elif mode=='parabola':vals=[zero,expression('(kappa+3)**2*f4',newnames),expression('2*(kappa+3)*f4',newnames),T]
  elif mode=='parallel':vals=[zero,scale(T,Q(9)),scale(T,Q(6)),T]
  elif mode=='linear':vals=[zero]*4
  else:raise AssertionError(mode)
  mapping=dict(zip(['f1','f2','f3','f4'],vals));return [mapping[name] if name in mapping else z[name] for name in names]
 def homcomp(poly,values,den,width):
  degree=max(e[0]+e[1] for e in poly);lift={e+(degree-e[0]-e[1],):c for e,c in poly.items()};return substitution(lift,values+[den],width),degree
 def branch(mode):
  if mode in branch_cache:return branch_cache[mode]
  ip=root/mode/'input.json';inp=load(ip);assert inp['input_sha256']==sha(source);P=inp['parameters']+inp['variables'];width=len(P);ren=remap(P,mode);base=[substitution(p,ren,width) for p in prepared];base=[p for p in base if p];boundary=[]
  if mode not in ['graph','conic']:
   nt=(['t','l'] if mode=='parallel' else ['t'])+P;nnt=len(nt);z=const(1,nnt)
   if mode=='cubicgraph':U=expression('t',nt);B=expression('1-3*t+f2*t**2+f1*t**3',nt)
   elif mode=='quadraticgraph':U=expression('t',nt);B=expression('1-3*t+f2*t**2',nt)
   elif mode=='parabola':U=expression('-1+t-f4*t**2',nt);B=expression('kappa*t-(kappa+3)*(-1+t-f4*t**2)',nt);z=expression('kappa',nt)
   elif mode=='parallel':U=expression('t',nt);B=expression('l-3*t',nt)
   else:U=expression('t',nt);B=expression('1-3*t',nt)
   coords=[U,B]+remap(nt,mode);GH=[]
   for poly,limit in [(g,1),(h,2)]:
    val,degree=homcomp(poly,coords,z,nnt);GH.append(val)
    if mode=='parallel':
     for (k,l),c in coefficient_buckets(val,2).items():
      if k>limit:assert l<2;boundary.append(c)
    else:boundary += [p for (k,),p in coefficient_buckets(val,1).items() if k>limit]
   if mode!='parallel':
    vf=['t','v']+P;nv=len(vf);den=const(1,nv)
    if mode=='parabola':Uf=expression('v-1+t-f4*t**2',vf);Bf=expression('kappa*t-(kappa+3)*(v-1+t-f4*t**2)',vf);den=expression('kappa',vf)
    else:
     Uf=expression('t',vf);Bf=expression({'cubicgraph':'1-3*t+f2*t**2+f1*t**3-v','quadraticgraph':'1-3*t+f2*t**2-v','linear':'1-3*t-v'}[mode],vf)
    val,degree=homcomp(f,[Uf,Bf]+remap(vf,mode),den,nv);assert val==mul(power(den,degree,nv),variable(1,nv));assert add(mul(derivative(Uf,0),derivative(Bf,1)),mul(derivative(Uf,1),derivative(Bf,0)),Q(-1))==scale(den,Q(-1))
   equal_generators(boundary,[decode(p,width) for p in inp['boundary_equations']]);boundary_cache[mode]=boundary
  equal_generators(base+boundary,[decode(p,width) for p in inp['equations']]);branch_cache[mode]=(P,base+boundary);return P,base+boundary
 def fieldcert(mode,targettext):
  P,derived=branch(mode);width=len(P);folder=root/mode;cp=folder/'certificate.json';cert=load(cp);assert cert['source_sha256']==sha(HERE/'weighted_unmarked7_lift.py') and cert['input_sha256']==sha(folder/'input.json') and cert['singular_source_sha256']==sha(folder/'lift.sing') and cert['variables']==P;gen=[decode(p,width) for p in cert['generators']];equal_generators(derived,gen);target=decode(cert['target'],width);assert target==expression(targettext,P);multipliers=[decode(p,width) for p in cert['multipliers']];assert len(gen)==len(multipliers);result={}
  for a,b in zip(gen,multipliers):result=add(result,mul(a,b))
  assert result==target;cases.append(dict(mode=mode,certificate_sha256=sha(cp),input_sha256=sha(folder/'input.json'),target=targettext,multiplier_monomials=sum(map(len,multipliers)),max_multiplier_degree=max(sum(e) for p in multipliers for e in p)));print('REPLAYED FIELD',mode,flush=True)
 fieldcert('graph','f3**16');fieldcert('conic','f4**2*(f2*f4-f3**2/4)**7');fieldcert('cubicgraph','f2*f1**5')
 def reduced(mode):
  base='cubicgraph' if mode=='purecubic' else mode;P,orig=branch(base);ip=root/'reduced'/mode/'input.json';inp=load(ip);assert inp['input_sha256']==sha(root/base/'input.json') and inp['source_sha256']==sha(HERE/'weighted_unmarked7_reduce.py');N=inp['variables'];width=len(N);linear=boundary_cache[base];params=P[:len(P)-len(tails)]
  if mode=='purecubic':
   vals=[{} if name=='f2' else variable(i,len(P)) for i,name in enumerate(P)];orig=[substitution(p,vals,len(P)) for p in orig];linear=[substitution(p,vals,len(P)) for p in linear];P=[name for name in P if name!='f2'];remove=[{} if name=='f2' else variable(P.index(name),len(P)) for name in branch_cache[base][0]];orig=[substitution(p,remove,len(P)) for p in orig];linear=[substitution(p,remove,len(P)) for p in boundary_cache[base]]
  allowed={'purecubic':{'f1'},'quadraticgraph':{'f2'},'parabola':{'f4','kappa'},'parallel':{'f4'}}[mode];rank_=linear_audit(linear,P,N,inp['linear_solution'],allowed);assert rank_==inp['linear_rank'];vals=[laurent(inp['linear_solution'][name],N) if name in inp['linear_solution'] else variable(N.index(name),width) for name in P];expected=[clear_monomial_denominator(substitution(p,vals,width),width) for p in orig if p];expected=[p for p in expected if p];D=expression(inp['linear_denominator'],N);assert D==expression('f4' if mode=='parabola' else '1',N);sat={'purecubic':'f1','quadraticgraph':'f2','parabola':'kappa*f4**2','parallel':'f4*(4*f4-1)'}[mode];assert expression(inp['saturation'],N)==expression(sat,N);expected.append(expression('1-branch_inv*('+sat+')',N));gen=[decode(p,width) for p in inp['equations']];equal_generators(expected,gen)
  return ip,inp,gen,vals
 for mode in ['purecubic','quadraticgraph','parabola']:
  ip,inp,gen,vals=reduced(mode);width=len(inp['variables']);cp=ip.parent/'certificate.json';cert=load(cp);assert cert['source_sha256']==sha(HERE/'weighted_unmarked7_unit.py') and cert['input_sha256']==sha(ip) and cert['singular_source_sha256']==sha(ip.parent/'lift.sing');multipliers=[decode(p,width) for p in cert['multipliers']];assert len(multipliers)==len(gen);result={}
  for a,b in zip(gen,multipliers):result=add(result,mul(a,b))
  assert result==const(1,width);cases.append(dict(mode=mode,certificate_sha256=sha(cp),input_sha256=sha(ip),target='1',multiplier_monomials=sum(map(len,multipliers)),max_multiplier_degree=max(sum(e) for p in multipliers for e in p)));print('REPLAYED UNIT',mode,flush=True)
 # Parallel-line proof: after the complete seven-row reduction, G is linear
 # in t and H quadratic, with deg_L A<=2, deg_L C<=1, B(0)=E(0)=0,D(0)=-3.
 ip,inp,gen,vals=reduced('parallel');N=['t','l']+inp['variables'];width=len(N);W=branch_cache['parallel'][0];mapping=[expression(inp['linear_solution'][name],N) if name in inp['linear_solution'] else variable(N.index(name),width) for name in W];ren=remap(['t','l']+W,'parallel');embed=[variable(N.index(name),width) if name in N else mapping[W.index(name)] for name in ['t','l']+W];coords=[expression('t',N),expression('l-3*t',N)]+[substitution(p,embed,width) for p in ren];G,H=[substitution(p,coords,width) for p in [g,h]];phi=expression('1-l+f4*l**2',N);assert substitution(f,coords,width)==phi
 def coefft(poly,k):return {e:c for e,c in poly.items() if e[0]==k}
 def strip_t(poly,k):return {tuple([0]+list(e[1:])):c for e,c in poly.items() if e[0]==k}
 assert max(e[0] for e in G)<=1 and max(e[0] for e in H)<=2;A,B=strip_t(G,1),strip_t(G,0);C,D,E=strip_t(H,2),strip_t(H,1),strip_t(H,0);assert max(e[1] for e in A)<=2 and max(e[1] for e in C)<=1
 atzero=[variable(i,width) if i!=1 else {} for i in range(width)];assert not substitution(B,atzero,width) and not substitution(E,atzero,width) and substitution(D,atzero,width)==const(-3,width)
 phip=derivative(phi,1);Ap,Bp,Cp,Dp,Ep=[derivative(p,1) for p in [A,B,C,D,E]];Wr1=add(scale(mul(B,C),Q(2)),mul(A,D),Q(-1));Wr0=add(mul(B,D),scale(mul(A,E),Q(-2)));c2=mul(phi,add(mul(A,Cp),scale(mul(Ap,C),Q(-2))));c1=add(scale(mul(phip,Wr1),Q(-1)),mul(phi,add(add(mul(A,Dp),mul(Ap,D),Q(-1)),scale(mul(Bp,C),Q(-2)))));c0=add(scale(mul(phip,Wr0),Q(-1)),mul(phi,add(mul(A,Ep),mul(Bp,D),Q(-1))));J=add(add(mul(c2,power(variable(0,width),2,width)),mul(c1,variable(0,width))),c0);assert J==substitution(add(qj,const(1,nn)),coords,width)
 # The final written divisibility contradiction is checked as a universal
 # algebra identity after C=0,A=alpha*phi,D=d (constants alpha,d).
 sn=['l','T','alpha','d','B','Bp','E','Ep'];sw=len(sn);phi0=expression('1-l+T*l**2',sn);phip0=derivative(phi0,0);target=add(scale(mul(phip0,add(expression('B*d',sn),scale(mul(expression('alpha',sn),mul(phi0,expression('E',sn))),Q(-2)))),Q(-1)),mul(phi0,add(mul(mul(expression('alpha',sn),phi0),expression('Ep',sn)),expression('Bp*d',sn),Q(-1))));expanded=expression('alpha*(2*(1-l+T*l**2)*(-1+2*T*l)*E+(1-l+T*l**2)**2*Ep)-d*((-1+2*T*l)*B+(1-l+T*l**2)*Bp)',sn);assert target==expanded
 print('PARALLEL POLYNOMIAL IDENTITIES CHECKED',flush=True)
 # The fixed-point exceptions A=1 or B=0 are exact units. All other boundary
 # G=A*u+B have a nonzero fixed point, supplying KCUB11 without assuming one.
 P,orig=branch('linear');W=P;wn=['t']+W;G0=substitution(g,[expression('t',wn),expression('1-3*t',wn)]+remap(wn,'linear'),len(wn));Ac={e[1:]:c for e,c in G0.items() if e[0]==1};Bc={e[1:]:c for e,c in G0.items() if e[0]==0}
 for mode,extra in [('Aone',add(Ac,const(1,len(W)),Q(-1))),('Bzero',Bc)]:
  folder=root/'tails'/mode;ip=folder/'input.json';inp=load(ip);assert inp['input_sha256']==sha(root/'linear/input.json') and inp['source_sha256']==sha(HERE/'weighted_unmarked7_tail.py');N=inp['variables'];width=len(N);linear=boundary_cache['linear']+[extra];claim={name:text for name,text in inp['substitution'].items() if name not in N};linear_audit(linear,W,N,claim,set());vals=[expression(inp['substitution'][name],N) for name in W];expected=[substitution(p,vals,width) for p in orig];gen=[decode(p,width) for p in inp['equations']];equal_generators(expected,gen);cp=folder/'certificate.json';cert=load(cp);assert cert['source_sha256']==sha(HERE/'weighted_unmarked7_unit.py') and cert['input_sha256']==sha(ip) and cert['singular_source_sha256']==sha(folder/'lift.sing');multipliers=[decode(p,width) for p in cert['multipliers']];assert len(gen)==len(multipliers);result={}
  for a,b in zip(gen,multipliers):result=add(result,mul(a,b))
  assert result==const(1,width);cases.append(dict(mode=mode,certificate_sha256=sha(cp),input_sha256=sha(ip),target='1',multiplier_monomials=sum(map(len,multipliers)),max_multiplier_degree=max(sum(e) for p in multipliers for e in p)));print('REPLAYED TAIL',mode,flush=True)
 # Generic fixed-point algebra and diagonal normalization of nonzero first
 # coefficients; the actual prior tail theorem determines beta1 and map13.
 fn=['u','A','B','beta'];assert expression('A*(-beta)+B+beta',fn)==expression('B-(A-1)*beta',fn)
 assert laurent('mu*(-3/mu)',['mu','nu'])==const(-3,2) and laurent('nu*(-1/nu)',['mu','nu'])==const(-1,2)
 tail=load(HERE/'weighted-tail7/validation.json');assert tail['status']=='pass' and tail['source_sha256']==sha(HERE/'verify_weighted_tail7.py') and tail['unique_nonlinear_total']==13
 earlier=load(HERE/'weighted-newfirst7/validation.json');assert earlier['status']=='pass' and earlier['source_sha256']==sha(HERE/'verify_weighted_newfirst7.py') and earlier['completed_square_identity_checked']
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),input_sha256=sha(source),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),tail_classification_sha256=sha(HERE/'weighted-tail7/validation.json'),smooth_cubic_lemma_source='WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md',nonlinear_atoms=23,coefficient_parameters=21,assumed_collision_equations=0,cases=cases,cleared_parameter_identities=3,full_ring_unit_identities=5,multiplier_monomials=sum(c['multiplier_monomials'] for c in cases),max_multiplier_degree=max(c['max_multiplier_degree'] for c in cases),universal_three_variable_determinant_reconstructed=True,quotient_jacobian_checked=True,all_boundary_constraints_reconstructed=True,linear_parameterizations_certified_with_allowed_monomial_pivots=True,all_used_certificates_multiplied=True,parallel_line_polynomial_identities_checked=True,parallel_line_divisibility_argument='Written proof; not formalized.',all_exceptional_branches_covered=True,automatic_nonzero_collision_from_boundary_fixed_point=True,diagonal_normalization_nonzero_pair_checked=True,unique_nonlinear_total=13,seconds=round(time.monotonic()-started,3),scope='All23 degree<=7 negative-weight atoms, identity linear part and fixed first coefficients-3,-1, with det1 alone force the established13-term map. Diagonal conjugacy extends the13-term conclusion to any nonzero pair of those first coefficients. No collision assumed. Written smooth-cubic and parallel-line arguments complement eight exactly replayed certificates; global minima and zero first-coefficient branches remain open.')
 (root/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
