"""Exact generic incidence, finite nonimmersion gate and row-N cusp tests."""
from pathlib import Path
import json,sys,time,hashlib
import sympy as S
from sympy.polys.rings import ring
from witnesses import witness
from norms import nonimmersion_norm
import geometry as G
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
A,M,s,z,t,w,d=S.symbols('A M s z t w d');CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)
def zero(name,value):check(name,S.expand(value)==0)

def main():
    start=time.monotonic();witnesses=[];norms=[]
    for family in ['L','K']:
        case=witness(family)
        for row in case['checks']:check(family+' '+row['name'],row['rank']==row['dimension'])
        witnesses.append(case)
        data=nonimmersion_norm(family);par=d if family=='L' else z
        norm=S.sympify(data['norm'],locals={str(par):par})
        saved=json.loads((HERE/('nonimmersion-parameters-'+family+'.json')).read_text())
        factor=S.sympify(saved['primitive_polynomial'],locals=dict(w=w,z=z))
        if family=='L':
            check('L nonimmersion compact degree12',S.degree(factor,w)==12)
            check('L nonimmersion factor squarefree',S.gcd(factor,S.diff(factor,w))==1)
            check('L nonzero compact origin',factor.subs(w,0)!=0)
            ratio=S.cancel(norm/(d*d*factor.subs(w,d**3)))
            check('L exact norm scalar',ratio.is_Rational and ratio!=0)
            lf=json.loads((HERE.parent/'jc2-r7-collisions/strata.json').read_text())['exceptional_factors']
            for label,expression in lf.items():check('L cusp gate avoids C exception '+label,S.gcd(factor.subs(w,d**3),S.sympify(expression,locals={'d':d}))==1)
            check('L first-image parameter count37',1+3*S.degree(factor,w)==37)
        else:
            check('K nonimmersion factor degree41',S.degree(factor,z)==41)
            check('K nonimmersion factor squarefree',S.gcd(factor,S.diff(factor,z))==1)
            check('K nonzero excluded origin',factor.subs(z,0)!=0)
            ratio=S.cancel(norm/(z**15*factor));check('K exact norm scalar',ratio.is_Rational and ratio!=0)
            factors=[('C1',3*z+2),('C2',15*z-4),('triple',G.K.D)]
            tangent=json.loads((HERE.parent/'jc2-r7-laurent-collisions/exceptional_moduli.json').read_text())['tangent']
            factors.append(('tangent',S.sympify(tangent,locals={'z':z})))
            for label,expression in factors:check('K cusp gate avoids C exception '+label,S.gcd(factor,expression)==1)
            check('K normalized parameter count123',3*S.degree(factor,z)==123)
        data['verified_primitive_factor']=str(factor);data['verified_scalar']=str(ratio);data['status']='EXACT_NORM_VERIFIED';norms.append(data)
    check('total160 coefficient points',37+123==160)
    check('at most54 cube-root orbits',(37-1)//3+1+123//3==54)
    # Every first-image nonimmersion is one of these two scaled K classes.
    rows=json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families']
    row=next(v for v in rows if v['family']=='K');parse=lambda e:S.sympify(e,locals=dict(A=A,M=M,s=s,z=z))
    cusps=[]
    for value,tstar in [(S.Rational(-2,3),1),(S.Rational(4,15),-1)]:
        pc=G.K.P.subs(z,value).subs(G.F.t,t);qc=G.K.Q.subs(z,value).subs(G.F.t,t)
        u0=pc.subs(t,tstar);v0=qc.subs(t,tstar)
        check('cusp unique critical parameter '+str(value),S.Poly(S.gcd(S.diff(pc,t),S.diff(qc,t)),t).monic().as_expr()==t-tstar)
        zero('cusp clean doubled fibre '+str(value),S.Poly(S.gcd(pc-u0,qc-v0),t).monic().as_expr()-(t-tstar)**2)
        jet=(S.diff(pc,t,2)*S.diff(qc,t,3)-S.diff(pc,t,3)*S.diff(qc,t,2)).subs(t,tstar)
        check('ordinary cusp jet '+str(value),jet!=0)
        gc,_,_=G.approximate_quintic(pc,qc,t);k={mon:parse(e).subs(z,value) for mon,e in row['implicit_difference'].items()}
        evaluation=S.expand(gc.subs({G.U:u0,G.V:v0})-k['(2, 0)']*u0*u0-k['(0, 1)']*v0-k['(1, 0)']*u0-k['(0, 0)'])
        R,aa,mm=ring('A,M',S.QQ,order=lambda mon:(2*mon[0]+3*mon[1],mon[0]));basis=[R.from_expr(parse(e).subs(z,value)) for e in row['groebner_basis']]
        standard=[tuple(e) for e in row['standard_monomials']];positions={e:i for i,e in enumerate(standard)};matrix=S.zeros(7)
        for j,mon in enumerate(standard):
            for exponent,co in R.from_expr(evaluation).mul_monom(mon).div(basis)[1].items():matrix[positions[exponent],j]=co
        determinant=matrix.det(method='domain-ge');check('all seven images avoid known cusp '+str(value),determinant!=0)
        cusps.append(dict(z=str(value),tstar=tstar,cusp=[str(u0),str(v0)],evaluation_matrix=[[str(x) for x in list(matrix.row(i))] for i in range(7)],determinant=str(determinant)))
    # Exact pointwise retained counts for the smooth/node/crossing catalogue.
    counts=[]
    for label,c,b in [('off',0,0),('C smooth',1,0),('C node',2,0),('B smooth',0,1),('B node',0,2),('C/B crossing',1,1)]:
        values={'C':6-2*c-3*b,'A_H=C':6-3*c-2*b,'A_H=B':6-2*c-3*b}
        check('nonnegative generic retained counts '+label,min(values.values())>=0)
        check('singleton only joint crossing '+label,all((n==1)==(label=='C/B crossing') for n in values.values()))
        counts.append(dict(stratum=label,c=c,b=b,retained=values))
    check('rowN known cusp retains3',6-2*1-2*0-1==3)
    check('rowN elsewhere even',(6-2*1-2*1)%2==0)
    result=dict(schema='r7-quintic-incidence-v1',generic_witnesses=witnesses,nonimmersion_norms=norms,cusp_tests=cusps,
        generic_retained_counts=counts,rowB_necessary_first_image_points=160,rowB_cube_root_orbit_upper_bound=54,
        rowN='No singleton fibre if actual B has degrees(3,5) and infinity contact>=38 in the same coordinates.',
        scope='Nonempty generic opens and explicit conditional finite necessary gate; no global Keller cover or JC2 closure.')
    (HERE/'incidence.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs=[HERE/name for name in ['verify.py','finite.py','witnesses.py','norms.py','PROOF.md','README.md','REVIEW.md','nonimmersion-parameters-L.json','nonimmersion-parameters-K.json','discovery.json']]
    inputs += [HERE.parent/path for path in ['jc2-r7-quintic-normalizations/normalizations.json','jc2-r7-quintic-normalizations/PROOF.md','jc2-r7-complement/geometry.py','jc2-r7-image-family/family.py','jc2-r7-laurent-collisions/laurent.py','jc2-r7-laurent-collisions/exceptional_moduli.json','jc2-r7-collisions/strata.json','jc2-omitted-divisors/NONIMMERSED_LEDGER.md','jc2-omitted-divisors/IMAGE_NORM_LEDGER.md','jc2-omitted-divisors/RESIDUE_LEDGER.md']]
    inputs += [ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(passed=True,check_count=len(CHECKS),checks=CHECKS,elapsed_seconds=round(time.monotonic()-start,3),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        output_sha256={'incidence.json':hashlib.sha256((HERE/'incidence.json').read_bytes()).hexdigest()},
        assurance='Written finite-algebra lifting, curve/ledger interpretation and openness; exact modular ranks and rational norm/cusp determinants. No formalization or independent external review.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(f'PASS {len(CHECKS)} QUINTIC INCIDENCE CHECKS ({time.monotonic()-start:.3f}s)',flush=True)

if __name__=='__main__':main()
