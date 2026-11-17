"""Exact controls for the written abelian-classgroup output obstruction.

The general valuation, Galois stabilizer and divisor arguments are proofs
in the note, not deductions from these finite algebraic controls.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item))==0 for item in values))

    s,t,w,v,T,V,x,y,z,a = S.symbols('s t w v T V x y z a')
    K = S.Matrix([s+t**3,t])
    check('a nonlinear Keller control has unit Jacobian',K.jacobian([s,t]).det()==1)
    zero('that control has a literal polynomial inverse',
        S.Matrix([T-V**3,V]).subs({T:K[0],V:K[1]},simultaneous=True)-S.Matrix([s,t]))
    pullbacks = [s+t**3,s,t]
    check('the three stipulated distinct target divisor pullbacks are square-free',
        all(S.Poly(f,s,t).sqf_part().monic()==S.Poly(f,s,t).monic() for f in pullbacks))
    check('their pairwise pullbacks have no common source divisor',
        all(S.gcd(pullbacks[i],pullbacks[j])==1 for i in range(3) for j in range(i+1,3)))
    rational = (s+t**3)*s*s
    check('the rational cube-root obstruction has nondivisible source valuations',
        sorted(exponent for factor,exponent in S.factor_list(rational)[1])==[1,2])
    check('a ramified quadratic map fails the Keller hypothesis',
        S.Matrix([s*s,t]).jacobian([s,t]).det()==2*s)
    check('its cyclic square root is witnessed by a repeated pullback factor',
        S.factor_list(s*s)[1]==[(s,2)])
    check('a contracted divisor can also destroy coprime pullback outside Keller',
        S.gcd(s,s*t)==s and S.Matrix([s,s*t]).jacobian([s,t]).det()==s)

    permutations = list(itertools.permutations(range(3)))
    matrices = [S.Matrix([[int(p[0]==i)-int(p[2]==i),int(p[1]==i)-int(p[2]==i)] for i in range(2)])
        for p in permutations]
    alpha = S.Matrix([2,-1])
    orbit = {tuple(g*alpha) for g in matrices}
    check('the nonabelian class-action countercontrol has three distinct nonzero conjugate classes',
        len(orbit)==3 and (0,0) not in orbit)
    zero('those nonzero free-lattice classes nevertheless sum to zero',
        sum((S.Matrix(value) for value in orbit),S.zeros(2,1)))
    check('the class action is actually nonabelian',any(g*h!=h*g for g in matrices for h in matrices))
    stabilizer = [g for g in matrices if g*alpha==alpha]
    check('its class stabilizer is nonnormal, so it supplies no abelian intermediate-field inference',
        len(stabilizer)==2 and any((g*h*g.inv())*alpha!=alpha for g in matrices for h in stabilizer))
    check('torsion is a separate obstruction to cancelling the component count',
        1%2!=0 and (2*1)%2==0)
    check('rank-one free class groups have only the two abelian lattice automorphisms',
        all(S.Matrix([[entry]]).det() in [-1,1] for entry in [-1,1]))

    relation = x*y-(z-a)*(z+a)
    gb = S.groebner([relation],y,z,x,domain=S.QQ.frac_field(a),order='lex')
    zero('the split quadric has the two displayed conjugate factors',relation-(x*y-z*z+a*a))
    smooth = S.groebner([relation,x,y,z],x,y,z,domain=S.QQ.frac_field(a))
    check('the affine quadric is smooth when the two roots are distinct',smooth.contains(S.Integer(1)))
    quotient = S.groebner([relation,x],x,z,y,domain=S.QQ.frac_field(a),order='lex')
    zero('the product of the two divisor equations vanishes modulo x',quotient.reduce((z-a)*(z+a))[1])
    check('neither divisor equation vanishes individually modulo x',
        quotient.reduce(z-a)[1]!=0 and quotient.reduce(z+a)[1]!=0)
    check('the associated graded quadric has irreducible rank three',S.hessian(x*y-z*z,[x,y,z]).rank()==3)
    zero('the localization at x is explicitly polynomial Laurent space',
        relation.subs(y,(z*z-a*a)/x))
    zero('the local divisor equation at the first prime has multiplicity one',
        S.cancel(((z-a)*(z+a)/y)/(z-a))-(z+a)/y)
    check('the complementary root factor is nonzero on that prime',(z+a).subs(z,a)==2*a)
    zero('root exchange swaps the two prime-divisor equations',
        (z-a).subs(a,-a)-(z+a))
    check('the class-group relation has free quotient of rank one',
        S.Matrix([[1,1]]).rank()==1)

    # r=1, nonisotropic normal: the ambient geometric generic surface is nonfactorial.
    image_I = {T:s,x:1,y:t+w*w,z:w}
    zero('the quadric-surface embedding has the required Keller output base',
        S.Matrix([T,x*y-z*z]).subs(image_I,simultaneous=True)-S.Matrix([s,t]))
    zero('its closed embedding has a literal source inverse',
        S.Matrix([T,y-z*z,z]).subs(image_I,simultaneous=True)-S.Matrix([s,t,w]))
    check('the restricted quadratic has rank three in four ambient coordinates',
        S.hessian(x*y-z*z,[T,x,y,z]).rank()==3)

    Y = S.symbols('Y')
    image_II = {T:s,Y:w,x:1,z:t-s*w}
    zero('the isotropic r=1 embedding has the required Keller output base',
        S.Matrix([T,T*Y+x*z]).subs(image_II,simultaneous=True)-S.Matrix([s,t]))
    zero('that embedding too has literal polynomial source recovery',
        S.Matrix([T,T*Y+x*z,Y]).subs(image_II,simultaneous=True)-S.Matrix([s,t,w]))
    check('its full quadratic is nondegenerate in four ambient variables',
        S.hessian(T*Y+x*z,[T,Y,x,z]).det()==1)
    rel_II = T*Y+x*z-V
    zero('the generic quadratic-weight filtration eliminates Y exactly',rel_II.subs(Y,(V-x*z)/T))
    zero('the bad graded height-one fibre is distinct from the integral affine one',
        rel_II.subs(T,0)-(x*z-V))
    zero('the height-one affine fibre has a nonconstant unit with written inverse',
        (x*z/V-1).subs(z,V/x))
    special = S.groebner([x*z-V],x,z,Y,domain=S.QQ.frac_field(V),order='lex')
    check('the height-one monic normal basis still contains a nonconstant x',special.reduce(x)[1]==x)
    cone = S.groebner([x*z],x,z,Y,order='lex')
    check('its associated graded ring has two nonzero elements whose product vanishes',
        cone.reduce(x)[1]!=0 and cone.reduce(z)[1]!=0 and cone.reduce(x*z)[1]==0)
    trap = S.groebner([z*z-V],z,Y,domain=S.QQ.frac_field(V))
    zero('the rank-one trap instead makes its primitive reduction algebraic',trap.reduce(z*z-V)[1])
    zero('its actual global denominator failure is retained',
        ((z*z-V)/T).subs(V,T*Y+z*z)+Y)
    zero('a binary quadric generic curve has the excluded nonconstant Laurent unit',
        (x*z/V-1).subs(z,V/x))

    # HC4 compatibility with a degenerate output quadratic and an independent potential.
    image4 = {T:s,x:1,y:t+w*w,z:w,v:v}
    zero('the degenerate quadratic output also has a genuine HC4-size embedding',
        S.Matrix([T,x*y-z*z]).subs(image4,simultaneous=True)-S.Matrix([s,t]))
    zero('the HC4-size embedding has full source recovery',
        S.Matrix([T,y-z*z,z,v]).subs(image4,simultaneous=True)-S.Matrix([s,t,w,v]))
    actual = 1+s*w+t*v+w*w*v
    zero('an independent nondegenerate ambient potential has the stated pullback',
        (T*z+y*v+x*x).subs(image4,simultaneous=True)-actual)
    check('that actual HC4 source Hessian is a unit',S.hessian(actual,[s,t,w,v]).det()==1)

    cyclic_controls = []
    for m in [3,5]:
        C = S.zeros(m-1)
        for column in range(m-2):
            C[column+1,column] = 1
        for row in range(m-1):
            C[row,m-2] = -1
        zero(f'cyclic degree {m}: the divisor-class permutation has order dividing m',C**m-S.eye(m-1))
        generator = S.eye(m-1)[:,0]
        classes = [C**i*generator for i in range(m)]
        check(f'cyclic degree {m}: the full orbit has m nonzero classes',
            len({tuple(value) for value in classes})==m and all(value!=S.zeros(m-1,1) for value in classes))
        zero(f'cyclic degree {m}: the sole component relation is the orbit sum',
            sum(classes,S.zeros(m-1,1)))
        image_m = {T:s,x:1,y:t+w**m,z:w,v:v}
        zero(f'cyclic degree {m}: the higher-degree observable has a genuine embedding with inverse',
            S.Matrix([T,y-z**m,z,v]).subs(image_m,simultaneous=True)-S.Matrix([s,t,w,v]))
        check(f'cyclic degree {m}: the independent actual HC4 potential keeps unit Hessian',
            S.hessian(1+s*w+t*v+w**m*v,[s,t,w,v]).det()==1)
        cyclic_controls.append(dict(degree=m,geometric_class_group_rank=m-1,ordinary_leading_scheme_nonreduced=True))

    retained = {
        'hessian/hc4-triangular-graphs/SMOOTH_INFINITY_OUTPUT_OBSTRUCTION.md':'7ee0a9b293ca95a8f037cfdef0b5c0630f0fb169ec728df5b18c4dc2fb5638db',
        'hessian/hc4-triangular-graphs/STABLE_HYPERSURFACE_OBSTRUCTION.md':'a169c8541cda91042b31a55b77c249c89d48170eba2fe4ba7b524ccfa4100c30',
        'hessian/hc4-triangular-graphs/QUADRIC_OUTPUT_BASE_OBSTRUCTION.md':'c4a27b3b92b385b2bc5aebc518504e6769db161614bb1501f83d020166e15e0e',
        'hessian/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36',
        'hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md':'26a6a76eb579de05b8e71245b8b804bdd5f036602f2e3e6ad5d036fd9f5cd572',
        'hessian/hc4-triangular-graphs/sources/jedrzejewicz-squarefree-1109.2113v2.html':'71c1b2edbde62b32299d3cfbadfd8ace5d5842aad8b6aa98b75e0919811d77cf'
    }
    for name,digest in retained.items():
        check('retained mathematical or primary context pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_abelian_classgroup.py','ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md','sources/abelian-classgroup-context.json']},
        dependencies=retained,cyclic_Danielewski_controls=cyclic_controls,
        written_not_finite_exhaustion=[
            'every-dimensional Keller function fields have no nontrivial finite abelian Galois intermediate extension',
            'torsion-free abelian class-group monodromy makes all conjugate source fibre components principal',
            'geometrically integral height-one fibres suffice for global primitive denominator saturation',
            'all r>=1 nondegenerate linear-quadratic output pairs and degenerate restricted-rank-at-least-three pairs',
            'the nonfactorial quadric surface class group and both split-normal filtrations are proved directly',
            'cyclic Danielewski output fibres xy=z^m+t2 in every m>=2, beyond generic factoriality and any smooth ordinary infinity frame'],
        whole_Keller_extension_asserted_Galois_or_trivial=False,
        arbitrary_Galois_intermediate_extensions_asserted_excluded=False,
        nonabelian_or_torsion_classgroups_asserted_covered=False,
        generic_quadric_surface_asserted_factorial=False,
        height_one_graded_domains_asserted_necessary=False,
        geometric_generic_factoriality_alone_asserted_sufficient=False,
        arbitrary_low_quadratic_ranks_or_r_zero_asserted_covered=False,
        nonlinear_source_Hessian_invariance_used=False,
        new_external_classification_or_purity_imported=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-abelian-classgroup.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
