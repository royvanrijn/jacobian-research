#!/usr/bin/env sage-python
"""Six fixed outer parameters: exact parent equations, sections and rank intake.

No Nagao scoring or elliptic point enumeration. A parent parameter and a
fibre parameter are distinct. Good-prime surface counts are computed later
from the exported degree-(8,12) equations under explicit semistability gates.
"""
import sys
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve,prod
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from icarm_curve245_mestre import fermigier_roots
from probe_mestre_fermigier_two_section_local_continuation import reconstructed_second_line,normalized_data
from screen_mestre_fermigier_two_section_height_triage import specialized_points
from research_runtime.store import checkpoint,digest
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from memory_rank_certificate import checked_rank
D=ROOT/'artifacts/local/elliptic-curves/mestre-parent-portfolio-intake-v1'

def main(index):
    protocol=cert.read(D/'protocol.json')
    assert all(cert.hashed(ROOT/n)==h for n,h in protocol['bindings'].items())
    row=protocol['rows'][index];u=cert.F(row['outer_u']);out=D/row['id']/'intake.json'
    if out.exists():raise FileExistsError('preserve fixed parent intake')
    data={'status':'RUNNING','protocol_hash':digest(protocol),'outer_u':str(u),'fibre_T':'1','stage':'symbolic_parent'};checkpoint(out,data)
    v,c2,m2=reconstructed_second_line(u);c1,m1=normalized_data(u,v)[4:6]
    original=fermigier_roots(u,v);roots=tuple(QQ(str((r-original[0])/(original[1]-original[0]))) for r in original)
    if len(set(roots))!=6:raise ArithmeticError('six distinct normalized roots required')
    R=PolynomialRing(QQ,'T');T=R.gen();Sring=PolynomialRing(R,'x');x=Sring.gen()
    F=prod((x-r-T)*(x-r+T) for r in roots);S=x**6
    for j in range(5,-1,-1):S+=R((F[6+j]-(S*S)[6+j])/2)*x**j
    remainder=S*S-F
    if remainder.degree()>4:raise ArithmeticError('quartic identity failed')
    quartic=[]
    for j in range(5):
        q,res=remainder[j].quo_rem(T*T)
        if res:raise ArithmeticError('square parameter content failed')
        quartic.append(q)
    extras=[]
    for c,m in ((c1,m1),(c2,m2)):
        xx=R(QQ(str(c))+QQ(str(m))*T);yy2=sum(q*xx**j for j,q in enumerate(quartic))
        if not yy2.is_square():raise ArithmeticError('extra ordinate not a rational polynomial square')
        yy=R(yy2.sqrt());assert yy*yy==yy2
        extras.append({'x_coefficients':list(map(str,xx.list())),'y_coefficients':list(map(str,yy.list()))})
    e,d,c,b,a=quartic;I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c*c*c
    A=-27*I;B=-27*J;delta=-16*(4*A**3+27*B**2);c4=-48*A
    if A.degree()>8 or B.degree()>12 or delta.degree()!=20 or delta.gcd(c4).degree()!=0 or delta.gcd(delta.derivative()).degree()!=0:raise ArithmeticError('expected20 I1 and infinity I4 geometry differs')
    data.update(stage='specialized_subgroup',roots=list(map(str,roots)),A_coefficients=list(map(str,A.list())),B_coefficients=list(map(str,B.list())),quartic_coefficients=[list(map(str,q.list())) for q in quartic],extra_sections=extras,discriminant_degree=20);checkpoint(out,data)
    model,images=specialized_points(u,cert.F(1));E=EllipticCurve(QQ,[QQ(str(v)) for v in model]);surface=EllipticCurve(QQ,[0,0,0,A(1),B(1)])
    if not E.is_isomorphic(surface):raise ArithmeticError('intake model is not the exported parent fibre')
    P=[E([QQ(str(z)) for z in p]) for p in images];points=[]
    for Q in P[1:]:
        halves=(Q-P[0]).division_points(2)
        if len(halves)!=1:raise ArithmeticError('unique rational half of covariant difference required')
        if halves[0].is_zero():raise ArithmeticError('nonzero marked divisor section required')
        points.append(tuple(cert.F(str(z)) for z in halves[0].xy()))
    cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];torsion=None
    from dataclasses import asdict
    for prime in _primes_up_to(997):
        if prime==2:continue
        if torsion is None and cert.short_curve_has_no_rational_2_torsion_modular_certificate(model,prime):torsion=prime
        try:sig=signature(cache,model,points,prime)
        except ValueError:continue
        before=len(pivots)
        for bits in sig.rows:insert(pivots,bits)
        if len(pivots)>before:signatures.append(asdict(sig))
        if len(pivots)==13 and torsion is not None:break
    if torsion is None:raise ArithmeticError('rational2-torsion exclusion not certified')
    selected=sorted(pivots);basis=[points[i] for i in selected];proof=checked_rank(model,basis,[s['prime'] for s in signatures],torsion)
    data.update(status='PASS_RANK13_PARENT' if len(basis)==13 else 'PARTIAL_SUBGROUP',stage='complete',curve=list(map(str,model)),points=[list(map(str,p)) for p in points],independent_column_indices=selected,independent_points=[list(map(str,p)) for p in basis],rank_certificate=proof,rank_lower_bound=len(basis),covariant_images=[list(map(str,p)) for p in images],argument='The14 rational quartic sections give13 divisor differences from the first section on their Jacobian. The covariant map sends these to doubled differences; unique rational halves at the smooth T=1 fibre are exactly checked. Independence at this fibre proves independence of the marked generic divisor sections. No generic-rank upper bound or full saturation claim.');checkpoint(out,data)
    print(row['id'],data['status'],len(basis),max(abs(q.numerator).bit_length() for q in model),flush=True)

if __name__=='__main__':main(int(sys.argv[1]))
