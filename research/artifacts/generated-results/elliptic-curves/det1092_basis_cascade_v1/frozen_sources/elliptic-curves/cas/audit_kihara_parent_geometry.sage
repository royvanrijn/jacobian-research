#!/usr/bin/env sage-python
"""Exact geometry of Kihara's path and one fixed six-root K3 parent."""
import argparse,json,hashlib,sys
from pathlib import Path
from fractions import Fraction
from sage.all import QQ,PolynomialRing,EllipticCurve,prod
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';sys.path.insert(0,str(ROOT/'elliptic-curves'))
from ecsearch.kihara import kihara_rank14_replay,verify_kihara_rank14_manifest
SOURCE=ART/'curve302_inverse_kihara_and_rank14_16_intake_v1.json'
def geometry(A,B):
 R=A.parent();common=A.gcd(B);scale=R(1);removed=[]
 for f,n in common.factor():
  va=A.valuation(f);vb=B.valuation(f);k=min(va//4,vb//6)
  if k:scale*=f**k;removed.append({'factor':list(map(str,f.list())),'power':int(k)})
 A,ar=A.quo_rem(scale**4);B,br=B.quo_rem(scale**6);assert not ar and not br
 delta=-16*(4*A**3+27*B**2);assert delta.gcd(A)==1
 chi=max((A.degree()+3)//4,(B.degree()+5)//6);infinity=12*chi-delta.degree();assert A.degree()==4*chi and B.degree()==6*chi and infinity>0
 factors=[]
 for f,n in delta.squarefree_decomposition():
  if n==1:factors.append({'degree':int(f.degree()),'multiplicity':1,'coefficients':list(map(str,f.list()))})
  else:
   for g,k in f.factor():
    assert k==1;factors.append({'degree':int(g.degree()),'multiplicity':int(n),'coefficients':list(map(str,g.list()))})
 euler=infinity+sum(r['degree']*r['multiplicity'] for r in factors);assert euler==12*chi
 rootrank=infinity-1+sum(r['degree']*(r['multiplicity']-1) for r in factors)
 return {'minimal_A':list(map(str,A.list())),'minimal_B':list(map(str,B.list())),'removed_scale':list(map(str,scale.list())),'removed_factors':removed,'chi':int(chi),'minimal_degrees':[int(A.degree()),int(B.degree()),int(delta.degree())],'infinity_type':'I'+str(infinity),'finite_factors':factors,'geometric_root_rank':int(rootrank),'h11':int(10*chi),'geometric_MW_Hodge_upper_bound':int(10*chi-2-rootrank),'c4_delta_coprime':True}
def global_path():
 d=json.loads(SOURCE.read_text());R=PolynomialRing(QQ,'t');e,dd,c,b,a=[R(coeff) for coeff in d['quartic_scaled_coefficients_low_to_high']]
 I=12*a*e-3*b*dd+c*c;J=72*a*c*e+9*b*c*dd-27*a*dd*dd-27*b*b*e-2*c**3
 result=geometry(-27*I,-27*J);print('global path',result['chi'],result['minimal_degrees'],result['infinity_type'],[(r['degree'],r['multiplicity']) for r in result['finite_factors']],flush=True)
 return {'schema':'elliptic-curves.kihara-global-geometry.v1','status':'PASS','geometry':result,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'scope':'Minimal elliptic-surface geometry of the existing Kihara rank>=14 path. No exact MW rank or K3 parent interpretation assumed, and no point search.'}
def fixed_parent():
 R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();S=PolynomialRing(R,'x');x=S.gen();p=QQ(5);q=QQ(-18)
 aa=[QQ(0),(2*p*p+p*q+2*q*q)**2,2*(p+q)**2*(2*p*p+p*q+q*q),q*q*(4*p*p-p*q+4*q*q),p*(2*p-q)*(2*p*p+4*p*q+5*q*q),4*p**4+8*p**3*q+9*p*p*q*q-2*p*q**3+2*q**4]
 normalizer=aa[1];roots=[a/normalizer for a in aa];shifted=[r+T for r in roots]+[r-T for r in roots]
 F=prod(x-r for r in shifted);G=x**6
 for j in range(5,-1,-1):G+=R((F[6+j]-(G*G)[6+j])/2)*x**j
 rem=G*G-F;assert rem.degree()==4
 coeff=[]
 for i in range(5):v,r=rem[i].quo_rem(T*T);assert not r;coeff.append(v)
 quartic=S(coeff);points=[(K(r),K(G(r)/T)) for r in shifted]
 den=2*p*p+2*p*q+3*q*q
 constant=(8*p**6+28*p**5*q+58*p**4*q*q+69*p**3*q**3+76*p*p*q**4+40*p*q**5+22*q**6)/(den*normalizer)
 xx=K(constant+(2*p*p+4*p*q+5*q*q)*T/den);yy=K(quartic(xx)).sqrt();assert yy in K
 controlT=QQ(275424)/normalizer
 if yy(controlT)<0:yy=-yy
 points.append((xx,yy));assert all(Y*Y==quartic(X) for X,Y in points)
 e,d,c,b,a=coeff;I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
 A=-27*I;B=-27*J;result=geometry(A,B);assert result['chi']==2
 # Use P1 as quartic origin and compute a degree-one pointed map.
 X0,Y0=points[0];translated=quartic(x+X0);e0,d0,c0,b0,a0=[K(translated[i]) for i in range(5)];assert e0==Y0*Y0
 E=EllipticCurve(K,[A,B]);images=[]
 for X,Y in points[1:]:
  xx=X-X0;wx=(2*Y0*(Y+Y0)+d0*xx)/(xx*xx);wy=(2*(wx*wx-4*Y0*Y0*a0)*xx-2*d0*wx-4*Y0*Y0*b0)/(4*Y0)
  images.append(E([9*wx+3*c0,27*wy]))
 control=kihara_rank14_replay(Fraction(2));Ec=EllipticCurve(QQ,list(map(QQ,control.weierstrass_coefficients)));old=[Ec([QQ(c) for c in P]) for P in control.weierstrass_points]
 Et=EllipticCurve(QQ,[A(controlT),B(controlT)]);special=[Et([c(controlT) for c in P.xy()]) for P in images]
 expected=[P-old[0] for P in old[1:13]]
 matches=[iso for iso in Et.isomorphisms(Ec) if all(iso(P)==Q for P,Q in zip(special,expected))];assert len(matches)==1
 manifest=json.loads((ART/'kihara_rank14_t2_v1.json').read_text());verify_kihara_rank14_manifest(manifest)
 print('fixed parent',result['minimal_degrees'],result['infinity_type'],[(r['degree'],r['multiplicity']) for r in result['finite_factors']],'12 generic sections proved',flush=True)
 return {'schema':'elliptic-curves.kihara-fixed-parent.v1','status':'PASS','p':'5','q':'-18','normalized_roots':list(map(str,roots)),'normalizer':str(normalizer),'control_T':str(controlT),'quartic_coefficients':[list(map(str,f.list())) for f in coeff],
 'quartic_points':[[str(c) for c in P] for P in points],'raw_A':list(map(str,A.list())),'raw_B':list(map(str,B.list())),'generic_sections':[[str(c) for c in P.xy()] for P in images],
 'control_isomorphism':list(map(str,matches[0].tuple())),'control_differences':'P2-P1,...,P13-P1 in the certified fourteen-dimensional group with origin P15','geometry':result,
 'generic_Q_rank_lower_bound':12,'generic_Q_rank_upper_bound':'UNKNOWN','geometric_Picard_rank':'UNKNOWN','full_NS_lattice':'UNKNOWN','different_NS_from_948_and_468':'UNKNOWN',
 'scope':'One fixed Kihara six-root K3 parent and twelve exact rational sections. Their images at the published t2 control are twelve independent differences in its certified rank14 subgroup. No rank14 generic parent claim, exact generic rank, NS identification, point search or novelty claim.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['global','fixed'],required=True);a=p.parse_args();r=global_path() if a.mode=='global' else fixed_parent();r['checker_sha256']=hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest();out=ART/('kihara_'+a.mode+'_parent_geometry_v1.json');assert not out.exists();out.write_text(json.dumps(r,indent=2)+'\n')
