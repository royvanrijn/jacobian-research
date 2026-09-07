#!/usr/bin/env sage-python
"""Picard ranks and full rational MW saturation of the six existing parents."""
import argparse,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,lcm,gcd
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves'
COUNTS=ROOT/'artifacts/local/elliptic-curves/mestre-picard-count-v1'
HEIGHTS=ROOT/'artifacts/local/elliptic-curves/mestre-generic-height-v1'
PILOT=ROOT/'artifacts/local/elliptic-curves/mestre-parent-calibration-v1'
OUT=ART/'mestre_parent_picard_and_saturation_v1.json'

def finite_rank(model,points,primes,ell):
    E=EllipticCurve(QQ,model);P=[E([QQ(x),QQ(y)]) for x,y in points];projective=[]
    for point in P:
        den=lcm([c.denominator() for c in point]);v=[ZZ(c*den) for c in point];g=gcd(v);projective.append([c//g for c in v])
    rows=[];records=[]
    for prime in primes:
        F=GF(prime);e=EllipticCurve(F,[F(c) for c in model])
        if not e.discriminant():raise ArithmeticError('good finite specialization required')
        key=lambda point:tuple(int(c) for c in point)
        elements=e.points();multiples={key(ell*P):ell*P for P in elements};mask={key(P):0 for P in multiples.values()};reps=[e(0)]
        while len(mask)<len(elements):
            P=next(P for P in elements if key(P) not in mask);old=list(reps);size=len(old)
            for digit in range(1,ell):
                for i,R in enumerate(old):
                    rep=R+digit*P;reps.append(rep)
                    for T in multiples.values():
                        k=key(rep+T)
                        if k in mask:raise ArithmeticError('quotient cosets overlap')
                        mask[k]=i+size*digit
        dimension=ZZ(len(reps)).valuation(ell)
        if ell**dimension!=len(reps) or dimension>2:raise ArithmeticError('elliptic quotient dimension differs')
        reduced=[e([F(c) for c in P]) for P in projective]
        for j in range(dimension):rows.append([(mask[key(P)]//ell**j)%ell for P in reduced])
        records.append({'prime':prime,'group_order':len(elements),'ell_multiple_subgroup_order':len(multiples),'quotient_dimension':int(dimension)})
    rank=int(matrix(GF(ell),rows).rank())
    if rank!=11:raise ArithmeticError('selected generic seed is not injective in finite ell quotients')
    return {'modulus':ell,'rank':rank,'groups':records}

def frobenius(record):
    p=ZZ(record['prime']);n1,n2=[ZZ(c['surface_point_count']) for c in record['counts']]
    for c in record['counts']:
        if len(c['weierstrass_fibre_counts'])!=c['field_order']+1 or c['surface_point_count']!=sum(c['weierstrass_fibre_counts'])+sum(r['correction'] for r in c['resolution_corrections']):raise ArithmeticError('surface count sum differs')
    a=n1-1-p*p-18*p;b=QQ(a*a-(n2-1-p**4-18*p*p))/2
    if b.denominator()!=1:raise ArithmeticError('nonintegral residual coefficient')
    # Orthogonal reciprocity forces sign+ if b!=0. For b=0 the
    # sign+ Weil condition is |a|<=p, and sign- requires |a|<=2p.
    if b:sign=1
    elif p<abs(a)<=2*p:sign=-1
    else:raise ArithmeticError('this fixed audit requires a unique determinant sign')
    R=PolynomialRing(QQ,'z');z=R.gen();f=z**4-a*z**3+b*z*z-sign*p*p*a*z+sign*p**4;remaining=f;cycles=[];rational=0
    for s in (1,-1):
        multiplicity=0
        while remaining(s*p)==0:remaining=remaining//(z-s*p);multiplicity+=1
        if multiplicity:cycles.append({'sign':s,'multiplicity':multiplicity})
        if s==1:rational=multiplicity
    if remaining.degree()!=2 or remaining[2]!=1 or remaining[0]!=p*p:raise ArithmeticError('two algebraic eigenvalues and a reciprocal pair required')
    trace=-remaining[1]
    if trace in (-2*p,-p,0,p,2*p):raise ArithmeticError('remaining pair might be algebraic')
    if abs(trace)>=2*p:raise ArithmeticError('remaining pair violates strict Weil bound')
    return {'prime':int(p),'surface_counts':[int(n1),int(n2)],'residual_polynomial':list(map(str,f.list())),
      'algebraic_residual_eigenvalues':cycles,'nonalgebraic_pair_trace':str(trace),
      'arithmetic_NS_upper_bound':18+rational,'geometric_reduction_NS_rank':20,
      'artin_tate_over_Fp2_discriminant_squareclass_representative':str(trace*trace-4*p*p)}

def compute():
    if cert.read(COUNTS/'ledger.json')['status']!='PASS' or cert.read(HEIGHTS/'ledger.json')['status']!='PASS':raise ArithmeticError('all fixed computations must be terminal and successful')
    independent=ART/'mestre_generic_height_independent_v1.json';involution=ART/'mestre_base_involution_v1.json'
    odd=ART/'mestre_parent_calibration_seed_clouds_v1.json';geometry=ART/'mestre_parent_fibre_geometry_v1.json'
    paths=[Path(__file__).resolve(),COUNTS/'protocol.json',COUNTS/'ledger.json',HEIGHTS/'protocol.json',HEIGHTS/'ledger.json',independent,involution,odd,geometry]
    for p in (independent,involution,odd,geometry):
        if cert.read(p)['status']!='PASS':raise ArithmeticError('input theorem required')
    rows=[]
    for u in cert.read(HEIGHTS/'protocol.json')['outer_parameters']:
        hpath=HEIGHTS/('u'+u)/'height.json';seed_path=PILOT/('u'+u+'-unit')/'seed.json';h=cert.read(hpath);seed=cert.read(seed_path);paths += [hpath,seed_path]
        determinant=QQ(h['known_rank18_divisor_lattice_absolute_determinant'])
        if determinant!=468:raise ArithmeticError('known rank18 divisor determinant differs')
        witnesses=[]
        for r in cert.read(COUNTS/'protocol.json')['rows']:
            if r['outer_u']==u:
                path=COUNTS/r['id']/'counts.json';paths.append(path);witnesses.append(frobenius(cert.read(path)))
        if len(witnesses)!=2 or min(r['arithmetic_NS_upper_bound'] for r in witnesses)!=18:raise ArithmeticError('matching rational NS upper bound required')
        ratio=QQ(witnesses[0]['artin_tate_over_Fp2_discriminant_squareclass_representative'])/QQ(witnesses[1]['artin_tate_over_Fp2_discriminant_squareclass_representative'])
        if ratio.is_square():raise ArithmeticError('two-prime geometric rank20 exclusion failed')
        sigma=next(r for r in cert.read(involution)['rows'] if r['outer_u']==u)
        M=matrix(QQ,sigma['action_matrix']);identity=matrix.identity(QQ,11)
        if M*M!=identity or 11-(M-identity).rank()!=5 or 11-(M+identity).rank()!=6 or sigma['geometric_NS_lower_bound']!=19:raise ArithmeticError('involution geometric lower bound differs')
        finite=[];proof=seed['rank_certificate']
        finite.append(finite_rank(seed['curve'],seed['points'],[r['prime'] for r in proof['signatures']],2))
        source=next(r for r in cert.read(odd)['rows'] if r['id']=='u'+u+'-unit');third=next(a for a in source['audits'] if a['modulus']==3)
        finite.append(finite_rank(seed['curve'],seed['points'],[r['prime'] for r in third['signatures']],3))
        possible=[d for d in range(1,22) if ZZ(determinant)%(d*d)==0]
        if possible!=[1,2,3,6]:raise ArithmeticError('possible saturation indices differ')
        # The full arithmetic MW group has rank11 by arithmetic Shioda-Tate.
        # Nonzero geometric torsion is excluded by h>=4-(1/2+1/2+1)=2.
        # The finite2/3 injections exclude every nontrivial possible index.
        selected=h['seed_indices'];C=matrix(QQ,h['divisor_change_matrix']).matrix_from_rows_and_columns(selected,selected)
        action=C.transpose().inverse()*M*C.transpose()
        if any(c.denominator()!=1 for c in action.list()):raise ArithmeticError('full arithmetic basis action must be integral')
        G=matrix(QQ,h['seed_height_gram'])
        plus=matrix(ZZ,action-identity).right_kernel().basis_matrix();minus=matrix(ZZ,action+identity).right_kernel().basis_matrix()
        plus_gram=plus*G*plus.transpose();minus_gram=minus*G*minus.transpose()
        q=next(r for r in cert.read(geometry)['rows'] if r['outer_u']==u)
        R=PolynomialRing(QQ,'x');x=R.gen();zero_cubic=x**3+QQ(q['quotient_A_coefficients'][0])*x+QQ(q['quotient_B_coefficients'][0])
        if not zero_cubic.is_irreducible():raise ArithmeticError('fixed quotient zero cubic was expected irreducible')
        rows.append({'outer_u':u,'frobenius_witnesses':witnesses,'artin_tate_ratio_is_square':False,'artin_tate_ratio':str(ratio),
          'arithmetic_NS_rank':18,'geometric_NS_rank':19,'arithmetic_generic_MW_rank':11,'geometric_generic_MW_rank':12,
          'finite_saturation_injections':finite,'initial_possible_arithmetic_lattice_indices':possible,'certified_index':1,
          'full_arithmetic_NS_absolute_determinant':468,'arithmetic_MW_height_determinant':'117/4','geometric_torsion_order':1,
          'arithmetic_MW_base_involution':[[int(c) for c in row] for row in action.rows()],
          'arithmetic_invariant_rank':5,'arithmetic_anti_invariant_rank':6,
          'invariant_height_gram':[[str(c) for c in row] for row in plus_gram.rows()],
          'invariant_height_determinant':str(plus_gram.det()),'anti_invariant_height_determinant':str(minus_gram.det()),
          'all_Q_fibrations_MW_rank_upper_bound':16,'full_geometric_NS_determinant':'UNKNOWN',
          'quadratic_twist_quotient_generic_MW_rank':6,'quadratic_twist_quotient_arithmetic_NS_rank_upper_bound':17,
          'quadratic_twist_zero_cubic_coefficients':list(map(str,zero_cubic.list()))})
        print('u'+u,'PICARD Q/Qbar18/19; MW11/12; SATURATED Q NS DET468; + DET',plus_gram.det(),flush=True)
    return {'schema':'elliptic-curves.mestre-parent-picard-saturation.v1','status':'PASS','rows':rows,
      'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
      'scope':'Frobenius fixes exactly18 rational divisor directions at a certified good prime on each parent. Two incompatible Artin-Tate squareclasses bound geometric Picard rank by19, and the exact5+6 base-involution decomposition supplies the matching lower bound via the rank6 geometric rational quotient. Arithmetic/geometric generic MW ranks are11/12. Complete finite groups show the actual generic basis injects modulo2 and3 at T=1, excluding every possible index in the rank18 arithmetic NS lattice of determinant468. The full geometric NS lattice is not identified. Every Q-fibration on each parent has MW rank at most16; this does not bound specialized elliptic ranks. The twist quotient has an I0* with irreducible cubic at zero, forcing at least two nonrational divisor directions and arithmetic NS rank at most17. No new fibration, rootless-frame search or foundry admission.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve Picard/saturation certificate')
        checkpoint(OUT,result)
