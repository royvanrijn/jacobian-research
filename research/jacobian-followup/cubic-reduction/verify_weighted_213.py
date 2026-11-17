#!/usr/bin/env python3
"""Independent arithmetic and parity audit of the complete mixed-weight proof.

The compact-curve pole and univariate parity arguments are written proof,
not formal proof-assistant artifacts. Plane degree5 is an attributed import.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,re,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode,scale
from weighted_verify_utils import substitution,expression,equal_generators,coefficient_buckets,power
from verify_weighted_unmarked7 import linear_audit
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def strip(poly,k):return {(0,)+e[1:]:c for e,c in poly.items() if e[0]==k}
def main():
 started=time.monotonic();source=HERE/'weighted-sparse7/2_-1_-3/input.json';d=load(source);names=d['variables'];n=len(names);assert n==16 and d['linear_rank']==0 and d['source_sha256']==sha(HERE/'weighted_sparse7_search.py');atoms=d['atoms'];expected=set()
 for row,target in enumerate([2,-1,-3]):
  for a in range(8):
   for b in range(8):
    for c in range(8):
     if 2<=a+b+c<=7 and 2*a-b-3*c==target:expected.add((row,(a,b,c)))
 assert {(a['row'],tuple(a['exponent'])) for a in atoms}==expected;N=3+n;F=[variable(i,N) for i in range(3)]
 for atom in atoms:
  e=list(atom['exponent'])+[0]*n;e[3+names.index(atom['coefficient'])]=1;F[atom['row']]=add(F[atom['row']],{tuple(e):Q(1)})
 det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),const(1,N),Q(-1));buckets=coefficient_buckets(det,3);prepared=[decode(p,n) for p in d['equations']];equal_generators(list(buckets.values()),prepared);assert all(2*a==b+3*c for a,b,c in buckets)
 qn=['u','v']+names;nn=len(qn);fgh=[const(1,nn),variable(0,nn),variable(1,nn)]
 for atom in atoms:
  b,c=atom['exponent'][1:];coef=variable(qn.index(atom['coefficient']),nn);fgh[atom['row']]=add(fgh[atom['row']],mul(coef,expression('u**%d*v**%d'%(b,c),qn)))
 f,g,h=fgh;assert f==expression('1+c1*u**2+c2*u*v+c3*u**4+c4*v**2',qn)
 for poly,parity in [(f,0),(g,1),(h,1)]:assert all((e[0]+e[1])%2==parity for e in poly)
 fu,fv,gu,gv,hu,hv=[derivative(p,j) for p in fgh for j in range(2)];J=add(add(mul(f,add(mul(gu,hv),mul(gv,hu),Q(-1))),scale(mul(g,add(mul(fu,hv),mul(fv,hu),Q(-1))),Q(1,2))),scale(mul(h,add(mul(gu,fv),mul(gv,fu),Q(-1))),Q(3,2)));assert add(J,const(1,nn),Q(-1))=={(b,c)+e:co for (a,b,c),poly in buckets.items() for e,co in poly.items()}
 # Even quartic chart and sign symmetry for the two-infinity-point pole lemma.
 yy=expression('2*c4*v+c2*u',qn);quartic=expression('-4*c3*c4*u**4+(c2**2-4*c1*c4)*u**2-4*c4',qn);assert add(power(yy,2,nn),quartic,Q(-1))==scale(mul(variable(qn.index('c4'),nn),f),Q(4));assert fv==yy
 Z=['z','Y','mu','nu','P','T'];zz=len(Z);chart=expression('Y**2+4*P*T-(nu**2-4*mu*T)*z**2+4*T*z**4',Z);assert {e:c for e,c in chart.items() if e[0]==0}==expression('Y**2+4*P*T',Z);assert expression('(8*m+2)-6*m',['m'])==expression('2*m+2',['m'])
 root=HERE/'weighted-213';certs=[]
 for mode,target in [('conic','(c1*c4-c2**2/4)**4'),('linearv','c2**6')]:
  folder=root/mode;inp=load(folder/'input.json');P=inp['parameters']+inp['variables'];width=len(P);zero='c3' if mode=='conic' else 'c4';vals=[{} if name==zero else variable(P.index(name),width) for name in names];derived=[substitution(p,vals,width) for p in prepared];gen=[decode(p,width) for p in inp['equations']];equal_generators(derived,gen);cp=folder/'certificate.json';c=load(cp);assert c['source_sha256']==sha(HERE/'weighted_213_lift.py') and c['input_sha256']==sha(folder/'input.json') and c['singular_source_sha256']==sha(folder/'lift.sing') and c['variables']==P;assert gen==[decode(p,width) for p in c['generators']];t=decode(c['target'],width);assert t==expression(target,P);mult=[decode(p,width) for p in c['multipliers']];assert len(gen)==len(mult);res={}
  for a,b in zip(gen,mult):res=add(res,mul(a,b))
  assert res==t;certs.append(dict(mode=mode,input_sha256=sha(folder/'input.json'),certificate_sha256=sha(cp),target=target,multiplier_monomials=sum(map(len,mult)),max_multiplier_degree=max(sum(e) for p in mult for e in p)))
 # Global transverse source-degree constraints in the two residual branches.
 tails=names[4:];lemma_cases=[]
 for mode in ['linearv','parallel']:
  ip=root/'lemmas'/mode/'input.json';inp=load(ip);assert inp['source_sha256']==sha(HERE/'weighted_213_lemmas.py') and inp['input_sha256']==sha(source);P=inp['parameters']+tails;nt=['t','l']+P;nnt=len(nt);U=expression('l' if mode=='linearv' else 't',nt);V=expression('t' if mode=='linearv' else 'l-sigma*t',nt);vals=[{}]*4+[variable(nt.index(name),nnt) for name in tails];GH=[substitution(p,[U,V]+vals,nnt) for p in [g,h]];linear=[]
  for poly,limit in zip(GH,[1,3]):
   for (t,l),co in coefficient_buckets(poly,2).items():
    if t>limit:assert l<2;linear.append(co)
  equal_generators(linear,[decode(p,len(P)) for p in inp['linear_equations']]);new=inp['parameters']+inp['free_tail_variables'];rename=lambda name:('g'+name[1:] if name.startswith('c') else name);replace=lambda text:re.sub(r'\bc[0-9]+\b',lambda m:rename(m.group()),text);claims={rename(k):replace(v) for k,v in inp['solution'].items()};assert linear_audit(linear,list(map(rename,P)),list(map(rename,new)),claims,set())==inp['linear_rank'];nr=['t','l']+new;nrn=len(nr);sv=[expression(inp['solution'][name],nr) if name in inp['solution'] else variable(nr.index(name),nrn) for name in P];G,H=[substitution(poly,[variable(0,nrn),variable(1,nrn)]+sv,nrn) for poly in GH];assert max(e[0] for e in G)<=1 and max(e[0] for e in H)<=3;A,B=strip(G,1),strip(G,0);C,D,E,F0=[strip(H,k) for k in [3,2,1,0]];assert max(e[1] for e in A)<=2 and max(e[1] for e in C)==0 and max(e[1] for e in D)<=1
  for poly,parity in [(A,0),(B,1),(C,0),(D,1),(E,0),(F0,1)]:assert all(e[1]%2==parity for e in poly)
  atzero=[{} if i==1 else variable(i,nrn) for i in range(nrn)];assert not substitution(B,atzero,nrn) and not substitution(D,atzero,nrn) and not substitution(F0,atzero,nrn)
  if mode=='linearv':assert substitution(E,atzero,nrn)==const(1,nrn)
  else:assert substitution(E,atzero,nrn)==expression('-sigma',nr) and substitution(A,atzero,nrn)==expression('1-sigma*c5',nr) and substitution(derivative(F0,1),atzero,nrn)==const(1,nrn)
  lemma_cases.append(dict(mode=mode,input_sha256=sha(ip),linear_rank=inp['linear_rank'],parity_and_source_bounds_checked=True))
 pairs=[(m,n) for m in range(6) for n in range(6) if m+n>0 and (m+n==1 or n==3*m)];assert pairs==[(0,1),(1,0),(1,3)]
 # Generic coefficient identities for the written parity/univariate proof.
 X=['t','A','B','C','D','E','F','Ap','Bp','Dp','Ep','Fp','phi','phip'];xn=len(X);G=expression('A*t+B',X);H=expression('C*t**3+D*t**2+E*t+F',X);Gu=expression('Ap*t+Bp',X);Hu=expression('Dp*t**2+Ep*t+Fp',X);phi,phip=expression('phi',X),expression('phip',X);jj=add(scale(mul(phip,add(mul(G,derivative(H,0)),scale(mul(derivative(G,0),H),Q(-3)))),Q(1,2)),mul(phi,add(mul(Gu,derivative(H,0)),mul(derivative(G,0),Hu),Q(-1))));assert strip(jj,3)==expression('3*phi*Ap*C',X);assert strip(jj,2)==expression('phip*(-A*D+3*B*C)/2+phi*(2*Ap*D+3*Bp*C-A*Dp)',X)
 j2=expression('-phip*A*D/2+phi*(2*Ap*D-A*Dp)',X);assert expression('A*(phip*D**2+2*phi*D*Dp)-4*phi*D**2*Ap',X)==scale(mul(expression('D',X),j2),Q(-2))
 # Constant first component: plane maps over K(a) have degree<=5 and a
 # unimodular linear matrix over K[a], so their formal inverse has no poles.
 fn=['a','b','c']+tails;fnn=len(fn);sv=[variable(i,fnn) for i in range(3)]+[{}]*4+[variable(fn.index(name),fnn) for name in tails];FF=[substitution(p,sv,fnn) for p in F];assert FF[0]==variable(0,fnn);assert max(e[1]+e[2] for p in FF[1:] for e in p)<=5
 linear_tail=[{e:c for e,c in p.items() if e[1]+e[2]==1} for p in FF[1:]];assert linear_tail==[expression('b+c5*a*c',fn),expression('c',fn)];assert all(e[1]+e[2]>0 for p in FF[1:] for e in p)
 impath=root/'plane-degree5-import.json';imp=load(impath);old=HERE/'weighted-zero7/sources/plane-quartic-import.json';assert imp['source_import_sha256']==sha(old) and imp['primary_pdf_sha256']==sha(old.with_name('nakai-baba-1977.pdf')) and imp['used_degrees']==[1,2,3,4,5] and not imp['reproved_here']
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),input_sha256=sha(source),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),linear_audit_sha256=sha(HERE/'verify_weighted_unmarked7.py'),coefficients=16,nonlinear_slots=16,universal_three_variable_determinant_reconstructed=True,double_cover_quotient_checked=True,even_first_odd_tails_checked=True,completed_even_quartic_checked=True,certificates=certs,certificate_count=2,multiplier_monomials=sum(c['multiplier_monomials'] for c in certs),max_multiplier_degree=max(c['max_multiplier_degree'] for c in certs),both_identities_multiplied=True,lemma_cases=lemma_cases,global_transverse_bounds_and_parity_checked=True,univariate_coefficient_identities_checked=True,plane_family_degree_bound=5,plane_family_unimodular_linear_part_checked=True,plane_import_sha256=sha(impath),plane_theorem_reproved=False,written_two_infinity_pole_argument_formalized=False,written_univariate_argument_formalized=False,complete_weights=[2,-1,-3],ordinary_degree_bound=7,all_maps_automorphisms=True,seconds=round(time.monotonic()-started,3),scope='All weight(2,-1,-3) Keller maps through degree7 are polynomial automorphisms. Written sign-involution elliptic-boundary pole counting gives P*T=0; two exact identities give the residual conic discriminant or cross coefficient zero. Necessary global transverse degree constraints and a written parity-univariate proof exclude every nonlinear-first branch. The linear-first family uses the attributed plane degree<=5 theorem and formal inverse over K[a]. Neither written curve theory nor the plane theorem is formalized/reproved by the arithmetic replay. The original broad query stays UNKNOWN; global minima and weight(1,-1,-3) remain open.')
 (root/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
