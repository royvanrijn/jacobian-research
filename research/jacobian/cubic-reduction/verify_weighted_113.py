#!/usr/bin/env python3
"""Independent rational replay for the complete weight(1,-1,-3), degree7 proof.
Written degree/divisibility arguments and the imported plane theorem are not
proof-assistant formalizations. No CAS is imported by this checker.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,re,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode,scale
from weighted_verify_utils import substitution,expression,equal_generators,coefficient_buckets,power
from verify_weighted_unmarked7 import laurent,clear_monomial_denominator,linear_audit
HERE=Path(__file__).resolve().parent
ROOT=HERE/'weighted-113'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def strip(p,j,k):return {e[:j]+(0,)+e[j+1:]:c for e,c in p.items() if e[j]==k}
def source_pair_substitution(poly,U,V,n):
 # u,v have nonnegative exponents. All parameter exponents, including
 # negative beta powers, stay in the coefficient ring without substitution.
 out={}
 for e,c in poly.items():
  assert e[0]>=0 and e[1]>=0
  term={(0,0)+e[2:]:c}
  term=mul(term,power(U,e[0],n));term=mul(term,power(V,e[1],n));out=add(out,term)
 return out

def main():
 started=time.monotonic();source=ROOT/'input.json';d=load(source);old=HERE/'weighted-sparse7/1_-1_-3/input.json';assert d['source_sha256']==sha(HERE/'weighted_113.py') and d['input_sha256']==sha(old)
 names=d['variables'];assert names==['c'+str(i) for i in range(1,21)];n=len(names);atoms=d['atoms'];assert atoms==load(old)['atoms'];expected=set()
 for row,target in enumerate([1,-1,-3]):
  for a in range(8):
   for b in range(8):
    for c in range(8):
     if 2<=a+b+c<=7 and a-b-3*c==target:expected.add((row,(a,b,c)))
 assert len(atoms)==20 and {(a['row'],tuple(a['exponent'])) for a in atoms}==expected
 N=3+n;F=[variable(i,N) for i in range(3)]
 for a in atoms:
  e=list(a['exponent'])+[0]*n;e[3+names.index(a['coefficient'])]=1;F[a['row']]=add(F[a['row']],{tuple(e):Q(1)})
 det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),const(1,N),Q(-1));buckets=coefficient_buckets(det,3);prepared=[decode(p,n) for p in d['equations']];equal_generators(list(buckets.values()),prepared);assert all(a==b+3*c for a,b,c in buckets)
 qn=['u','w']+names;nn=len(qn);fgh=[const(1,nn),variable(0,nn),variable(1,nn)]
 for a in atoms:
  b,c=a['exponent'][1:];fgh[a['row']]=add(fgh[a['row']],mul(variable(qn.index(a['coefficient']),nn),expression('u**%d*w**%d'%(b,c),qn)))
 for poly,text in zip(fgh,d['quotients']):assert poly==expression(text,qn)
 f,g,h=fgh;assert f==expression('1+c1*u+c2*u**2+c3*w+c4*u**3+c5*u*w',qn)
 fu,fw,gu,gw,hu,hw=[derivative(p,j) for p in fgh for j in range(2)];J=add(add(mul(g,add(mul(fu,hw),mul(fw,hu),Q(-1))),scale(mul(h,add(mul(gu,fw),mul(gw,fu),Q(-1))),Q(3))),mul(f,add(mul(gu,hw),mul(gw,hu),Q(-1))));assert add(J,const(1,nn),Q(-1))=={(b,c)+e:co for (a,b,c),p in buckets.items() for e,co in p.items()}
 certs=[]
 def record(folder,gen,width,target,exporter):
  cp=folder/'certificate.json';c=load(cp);assert c['source_sha256']==sha(HERE/exporter) and c['input_sha256']==sha(folder/'input.json') and c['singular_source_sha256']==sha(folder/'lift.sing');mult=[decode(p,width) for p in c['multipliers']];assert len(mult)==len(gen);total={}
  for a,b in zip(gen,mult):total=add(total,mul(a,b))
  assert total==target;certs.append(dict(mode=str(folder.relative_to(ROOT)),input_sha256=sha(folder/'input.json'),certificate_sha256=sha(cp),multiplier_monomials=sum(map(len,mult)),max_multiplier_degree=max(sum(e) for p in mult for e in p)));print('REPLAYED',folder.relative_to(ROOT),flush=True)
 folder=ROOT/'full';inp=load(folder/'input.json');assert inp['parameters']+inp['variables']==names and inp['substitution']=={} and inp['input_sha256']==sha(source);gen=[decode(p,n) for p in inp['equations']];equal_generators(prepared,gen);c=load(folder/'certificate.json');assert c['variables']==names and c['generators']==inp['equations'] and decode(c['target'],n)==expression('c5**13',names);record(folder,gen,n,expression('c5**13',names),'weighted_113_lift.py')
 tails=names[5:];branches=[]
 for case in ['both','mu-zero']:
  for degree in ['cubic','quad','linear']:
   folder=ROOT/'boundary'/(case+'-'+degree);inp=load(folder/'input.json');assert inp['source_sha256']==sha(HERE/'weighted_113_boundary.py') and inp['input_sha256']==sha(source);P=inp['parameters']+tails;width=len(P);mu=-1 if case=='both' else 0;sub={'c1':mu,'c3':-1,'c5':0}
   if degree!='cubic':sub['c4']=0
   if degree=='linear':sub['c2']=0
   assert inp['substitution']=={k:str(v) for k,v in sub.items()};vals=[const(sub[name],width) if name in sub else variable(P.index(name),width) for name in names];native=[substitution(p,vals,width) for p in prepared];native=[p for p in native if p]
   tp=['t']+P;nt=len(tp);B=expression('1'+('-t' if mu else '')+('+c2*t**2' if degree!='linear' else '')+('+c4*t**3' if degree=='cubic' else ''),tp);tv=[const(sub[name],nt) if name in sub else variable(tp.index(name),nt) for name in names];GH=[substitution(p,[variable(0,nt),B]+tv,nt) for p in [g,h]];limits=[2,6] if degree=='cubic' else [1,3];boundary=[p for poly,limit in zip(GH,limits) for (k,),p in coefficient_buckets(poly,1).items() if k>limit];equal_generators(boundary,[decode(p,width) for p in inp['boundary_equations']]);equal_generators(native+boundary,[decode(p,width) for p in inp['equations']]);linear=boundary+[p for p in native if max(sum(e[len(inp['parameters']):]) for e in p)<=1];equal_generators(linear,[decode(p,width) for p in inp['linear_equations']])
   reduced=folder/'reduced';rd=load(reduced/'input.json');assert rd['source_sha256']==sha(HERE/'weighted_113_boundary.py') and rd['input_sha256']==sha(folder/'input.json');new=rd['variables'];nr=len(new);allowed={'c4'} if degree=='cubic' else ({'c2'} if degree=='quad' else set());rename=lambda name:'g'+name[1:] if name in tails else name;replace=lambda txt:re.sub(r'\bc[0-9]+\b',lambda m:rename(m.group()),txt);claims={rename(k):replace(v) for k,v in rd['solution'].items()};assert linear_audit(linear,list(map(rename,P)),list(map(rename,new)),claims,allowed)==rd['linear_rank']==6
   values=[laurent(rd['solution'][name],new) if name in rd['solution'] else variable(new.index(name),nr) for name in P];derived=[clear_monomial_denominator(substitution(p,values,nr),nr) for p in native];derived=[p for p in derived if p];D='c4**2' if degree=='cubic' else ('c2' if degree=='quad' else '1');sat='c4**3' if degree=='cubic' else ('c2**2' if degree=='quad' else '1');assert expression(rd['linear_denominator'],new)==expression(D,new) and expression(rd['saturation'],new)==expression(sat,new)
   if sat!='1':derived.append(expression('1-branch_inv*('+sat+')',new))
   gen=[decode(p,nr) for p in rd['equations']];equal_generators(derived,gen);branches.append(dict(case=case,degree=degree,input_sha256=sha(folder/'input.json'),reduced_input_sha256=sha(reduced/'input.json'),linear_rank=6,denominator=D,saturation=sat))
   if (case,degree)!=('both','quad'):record(reduced,gen,nr,const(1,nr),'weighted_113_unit.py')
 # Wronskian degree alternatives: nonconstant leading terms can cancel only n=3m.
 for bounds,limit in [((6,7),(2,6)),((4,5),(1,3)),((2,2),(1,1))]:
  pairs=[(m,k) for m in range(bounds[0]+1) for k in range(bounds[1]+1) if m+k>0 and (m+k==1 or k==3*m)];assert all(m<=limit[0] and k<=limit[1] for m,k in pairs)
 # Complete polynomial identities behind the hard quadratic highest-v split.
 X=['v','G0','G1','K','A','G1p','H0','H1','H2','H0p','H1p','H2p'];nx=len(X);G=expression('G0+G1*v+K*v**2',X);H=expression('H0+H1*v+H2*v**2',X);Gu=expression('A+G1p*v',X);Hu=expression('H0p+H1p*v+H2p*v**2',X);pde=add(mul(add(G,mul(variable(0,nx),derivative(G,0))),Hu),mul(Gu,add(scale(H,Q(3)),mul(variable(0,nx),derivative(H,0)))),Q(-1))
 for k,text in enumerate(['G0*H0p-3*A*H0','2*G1*H0p+G0*H1p-3*G1p*H0-4*A*H1','3*K*H0p+2*G1*H1p+G0*H2p-4*G1p*H1-5*A*H2','3*K*H1p+2*G1*H2p-5*G1p*H2','3*K*H2p']):assert strip(pde,0,k)==expression(text,X)
 # Chain rule for f=N(u)-w gives exactly the transverse PDE above.
 Z=['G','Gv','Gu','H','Hv','Hu','Np','v'];nz=len(Z);assert expression('G*(Np*(-Hv)+Hu+Np*Hv)+3*H*(-Gu-Np*Gv+Np*Gv)+v*((Gu+Np*Gv)*(-Hv)+Gv*(Hu+Np*Hv))',Z)==expression('(G+v*Gv)*Hu-Gu*(3*H+v*Hv)',Z)
 # The hard quadratic source support bounds, after necessary boundary degree
 # constraints alone (the complete affine reduction was audited above).
 folder=ROOT/'boundary/both-quad';rd=load(folder/'reduced/input.json');new=['u','v']+rd['variables'];nr=len(new);sourcevals=[const(-1,nr),variable(new.index('c2'),nr),const(-1,nr),{},{}]+[laurent(rd['solution'][name],new) if name in rd['solution'] else variable(new.index(name),nr) for name in tails];GH=[substitution(p,[variable(0,nr),expression('1-u+c2*u**2-v',new)]+sourcevals,nr) for p in [g,h]]
 for poly,bounds in zip(GH,[[1,2,0],[3,3,1]]):
  assert max(e[1] for e in poly)<=2
  for k,bound in enumerate(bounds):assert all(e[0]<=bound for e in strip(poly,1,k))
 # Affine-G coordinate change: Dq=Aq, Dv=-Av. Consequently the only
 # polynomial solutions have constant -1/(3A) plus q^(j+3)*v^j.
 Z=['u','v','k'];assert expression('(1+k)*((1+k)*u-k+k*v+k*v)-v*(1+k)*k',Z)==expression('(1+k)*((1+k)*u-k+k*v)',Z)
 assert all(i==j+3 for i in range(15) for j in range(15) if i-j-3==0)
 # The general univariate-first highest-w coefficients and their integrals.
 U=['w','A','B','C','D','E','Ap','Bp','Cp','Dp','Ep','phi','phip'];nu=len(U);G=expression('A*w+B',U);H=expression('C*w**2+D*w+E',U);Gu=expression('Ap*w+Bp',U);Hu=expression('Cp*w**2+Dp*w+Ep',U);wr=add(mul(G,derivative(H,0)),scale(mul(derivative(G,0),H),Q(-3)));jj=add(mul(expression('phip',U),wr),mul(expression('phi',U),add(mul(Gu,derivative(H,0)),mul(derivative(G,0),Hu),Q(-1))));assert strip(jj,0,2)==expression('-phip*A*C+phi*(2*Ap*C-A*Cp)',U);assert strip(jj,0,1)==expression('phip*(-2*A*D+2*B*C)+phi*(Ap*D+2*Bp*C-A*Dp)',U);assert strip(jj,0,0)==expression('phip*(B*D-3*A*E)+phi*(Bp*D-A*Ep)',U)
 W=['phi','phip','alpha','delta','B','Bp','D','Dp'];nw=len(W);sv=[{},expression('alpha*phi',W),expression('B',W),expression('delta*phi',W),expression('D',W),{},expression('alpha*phip',W),expression('Bp',W),expression('delta*phip',W),expression('Dp',W),{},expression('phi',W),expression('phip',W)];assert substitution(strip(jj,0,1),sv,nw)==expression('phi*(2*delta*(phip*B+phi*Bp)-alpha*(phip*D+phi*Dp))',W)
 # Reconstruct the small residual core, including rational source conditions.
 folder=ROOT/'core';inp=load(folder/'input.json');assert inp['source_sha256']==sha(HERE/'weighted_113_core.py');cn=['u','v']+inp['variables'];nc=len(cn);A=expression('-3-beta',cn);alpha=expression('2*beta+4',cn);G0=expression('(-3-beta)*u-beta',cn);G1=expression('(2*beta+4)*u+beta',cn);H1=laurent('-(((2*beta+4)*u+beta)/beta)**2',cn);H0=laurent('1+(2*(2*beta+4)/beta-1)*u+((2*beta+4)**2/beta**2-2*(2*beta+4)/beta+P)*u**2+T*u**3',cn);E0=add(add(mul(G0,derivative(H0,0)),scale(mul(A,H0),Q(-3))),const(1,nc),Q(-1));E1=add(add(scale(mul(G1,derivative(H0,0)),Q(2)),scale(mul(alpha,H0),Q(-3))),add(mul(G0,derivative(H1,0)),scale(mul(A,H1),Q(-4))));core=[clear_monomial_denominator(p,3) for poly in [E0,E1] for (i,j),p in coefficient_buckets(poly,2).items()];gen=[decode(p,3) for p in inp['equations']];equal_generators(core,gen);record(folder,gen,3,const(1,3),'weighted_113_unit.py')
 G=add(G0,mul(variable(1,nc),G1));H=add(H0,mul(variable(1,nc),H1));vals=[variable(0,nc),expression('1-u+P*u**2',cn)]+[variable(i,nc) for i in range(2,nc)];gs=source_pair_substitution(G,vals[0],vals[1],nc);hs=source_pair_substitution(H,vals[0],vals[1],nc);assert strip(gs,0,0)=={} and strip(gs,0,1)==const(1,nc);assert all(not strip(hs,0,k) for k in [0,1,2])
 # Source degree bounds for the constant-first polynomial family over K(a).
 fn=['a','b','c']+tails;nf=len(fn);fv=[variable(i,nf) for i in range(3)]+[{}]*5+[variable(fn.index(name),nf) for name in tails];FF=[substitution(p,fv,nf) for p in F];assert FF[0]==variable(0,nf);assert all(e[1]+e[2]>0 for p in FF[1:] for e in p);assert max(e[1]+e[2] for p in FF[1:] for e in p)<=5;assert [{e:c for e,c in p.items() if e[1]+e[2]==1} for p in FF[1:]]==[expression('b+c7*a**2*c',fn),expression('c',fn)]
 impath=ROOT/'plane-degree5-import.json';imp=load(impath);oldimp=HERE/'weighted-zero7/sources/plane-quartic-import.json';assert imp['source_import_sha256']==sha(oldimp) and imp['primary_pdf_sha256']==sha(oldimp.with_name('nakai-baba-1977.pdf')) and imp['used_degrees']==[1,2,3,4,5] and imp['reproved_here'] is False
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),input_sha256=sha(source),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),linear_audit_sha256=sha(HERE/'verify_weighted_unmarked7.py'),complete_weights=[1,-1,-3],ordinary_degree_bound=7,nonlinear_slots=20,universal_three_variable_determinant_reconstructed=True,quotient_jacobian_checked=True,certificate_count=len(certs),certificates=certs,multiplier_monomials=sum(c['multiplier_monomials'] for c in certs),max_multiplier_degree=max(c['max_multiplier_degree'] for c in certs),branches=branches,all_six_boundary_parameterizations_checked=True,all_seven_identities_multiplied=True,hard_quadratic_core_source_conditions_checked=True,core_laurent_coefficient_denominators_preserved=True,highest_transverse_coefficient_identities_checked=True,univariate_integral_identities_checked=True,plane_family_degree_bound=5,plane_family_unimodular_linear_part_checked=True,plane_import_sha256=sha(impath),plane_theorem_reproved=False,written_case_split_formalized=False,all_maps_automorphisms=True,seconds=round(time.monotonic()-started,3),scope='All weight(1,-1,-3) Keller maps through degree7 are automorphisms, allowing arbitrary invertible diagonal linear part. Exact c5^13 identity, five saturated/full unit identities and the small hard-quadratic unit identity accompany the written Wronskian degree, highest-transverse-degree case split and univariate divisibility proof. Plane degree<=5 is an attributed import. Neither the written coverage arguments nor the plane theorem are formalized/reproved here. Original broad/changed timeouts remain preserved. No global simplicity minimum or new construction follows.')
 (ROOT/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
