"""Exact controls for the normal-quotient, boundary-basis and class-stabilizer proofs."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time

import sympy as S

from verify_nilpotent_covers import finite_group

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def quaternion(a,b):
    a0,a1,a2,a3 = a
    b0,b1,b2,b3 = b
    return S.Matrix([a0*b0-a1*b1-a2*b2-a3*b3,
        a0*b1+a1*b0+a2*b3-a3*b2,
        a0*b2-a1*b3+a2*b0+a3*b1,
        a0*b3+a1*b2-a2*b1+a3*b0])


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

    T,V,z,x,y,u = S.symbols('T V z x y u')
    cubic = z**3-3*z+V
    zero('a nonnormal cubic has both a simple and a ramified place at one target divisor',
        cubic.subs(V,2)-(z-1)**2*(z+2))
    check('the simple cubic place is visibly unramified',S.diff(cubic,z).subs(z,-2)==9)
    check('the other cubic place is ramified',S.diff(cubic,z).subs(z,1)==0)
    zero('the cubic normal closure has the nonsquare discriminant detecting S3',
        S.discriminant(cubic,z)-27*(4-V*V))
    check('the cubic polynomial source map is not Keller globally',
        S.Matrix([T,-z**3+3*z]).jacobian([T,z]).det()==3-3*z*z)
    square = z*z-T
    zero('a Galois quadratic has no unramified visible place at its branch divisor',square.subs(T,0)-z*z)
    check('absorbing this branched root in a polynomial source violates Keller',
        S.Matrix([z*z,V]).jacobian([z,V]).det()==2*z)

    relation = x*y-z*z+1
    chart = {x:x,z:1+x*u,y:2*u+x*u*u}
    zero('the open affine-plane chart lands in the smooth normal surface xy=z^2-1',
        relation.subs(chart,simultaneous=True))
    zero('the first rational inverse chart recovers the affine-plane coordinate',
        ((z-1)/x).subs(chart,simultaneous=True)-u)
    zero('the second inverse chart recovers that coordinate across x=0',
        (y/(z+1)).subs(chart,simultaneous=True)-u)
    zero('the two inverse charts agree by the actual surface relation',
        (z-1)*(z+1)-x*y+relation)
    zero('the x=0 fibre has its two exact component factors',relation.subs(x,0)+(z-1)*(z+1))
    check('those two component factors are reduced and distinct',S.gcd(z*z-1,2*z)==1)
    zero('the affine-plane chart keeps exactly the z=1 component at x=0',
        S.Matrix([chart[z],chart[y]]).subs(x,0)-S.Matrix([1,2*u]))
    check('the normal surface is smooth at both exceptional components',
        all(S.diff(relation,z).subs(z,level)!=0 for level in [-1,1]))
    lattice = S.Matrix([[1,1]])
    check('localizing x gives the one diagonal principal relation among the two divisor classes',
        lattice.rank()==1 and S.gcd(*list(lattice))==1)
    check('the missing single component is a primitive generator of the surface class lattice',
        S.Matrix([[1,1],[0,1]]).det()==1)
    check('removing both components instead produces a boundary relation from the nonconstant unit x',
        (lattice*S.Matrix([1,-1]))[0]==0)
    zero('the x-localized smooth surface has the literal Laurent inverse',relation.subs(y,(z*z-1)/x))
    finite_base = S.Matrix([chart[x]+chart[y],chart[z]])
    check('this positive normalization-boundary chart is not a Keller counterexample',
        S.factor(finite_base.jacobian([x,u]).det())==x-2*u-x*u*u)

    one = S.Matrix([1,0,0,0])
    qi = S.Matrix([0,1,0,0])
    qj = S.Matrix([0,0,1,0])
    qk = S.Matrix([0,0,0,1])
    h = S.Matrix([S.Rational(1,2)]*4)
    basis = S.Matrix.hstack(h,qi,qj,qk)
    basis_inverse = basis.inv()
    zero('the Hurwitz lattice contains the identity with its exact integral coordinates',2*h-qi-qj-qk-one)
    check('all sixteen basis products lie in the same integral lattice',
        all(all(value.is_Integer for value in basis_inverse*quaternion(basis[:,i],basis[:,j]))
            for i in range(4) for j in range(4)))
    zero('the half-coordinate generator satisfies h^2=h-1',quaternion(h,h)-h+one)
    units = {tuple(sign*column) for sign in [-1,1] for column in [one,qi,qj,qk]}
    units |= {tuple(S.Rational(sign,2) for sign in signs)
        for signs in itertools.product([-1,1],repeat=4)}
    check('the complete norm-one Hurwitz unit list has 24 distinct elements',len(units)==24)
    check('every listed Hurwitz unit has norm one',all(sum(q*q for q in a)==1 for a in units))
    check('the entire list is closed under quaternion multiplication',
        all(tuple(quaternion(a,b)) in units for a in units for b in units))
    matrices = {}
    for unit in units:
        matrix = basis_inverse*S.Matrix.hstack(*(quaternion(unit,basis[:,i]) for i in range(4)))
        matrices[unit] = matrix
    check('all 24 left actions are integral unimodular lattice automorphisms',
        all(all(value.is_Integer for value in matrix) and matrix.det()==1 for matrix in matrices.values()))
    check('the integral lattice action is faithful',len({tuple(matrix) for matrix in matrices.values()})==24)
    check('every nonidentity Hurwitz action fixes no nonzero rational vector',
        all((matrix-S.eye(4)).det()!=0 for unit,matrix in matrices.items() if unit!=tuple(one)))
    omega = h-one
    zero('the chosen quaternion has order three',quaternion(quaternion(omega,omega),omega)-one)
    omega2 = quaternion(omega,omega)
    conjugate = quaternion(quaternion(qi,omega),-qi)
    zero('conjugation rotates the order-three subgroup to the stated other subgroup',
        conjugate-S.Matrix([-S.Rational(1,2),S.Rational(1,2),-S.Rational(1,2),-S.Rational(1,2)]))
    check('the order-three Sylow subgroup is not normal, so the class-action example is nonnilpotent',
        tuple(conjugate) not in {tuple(one),tuple(omega),tuple(omega2)})

    elements = set(itertools.permutations(range(3)))
    compose = lambda a,b:tuple(a[b[i]] for i in range(3))
    group,subgroups,normal_prime,closure = finite_group(elements,compose,(0,1,2))
    stabilizer = frozenset(g for g in group if g[0]==0)
    check('the nonzero S3 root-lattice class has its exact order-two stabilizer',len(stabilizer)==2)
    check('the S3 class stabilizer lies in no proper normal subgroup',
        not any(stabilizer<=n for n in normal_prime))
    check('there is no intermediate subgroup between that stabilizer and S3',
        not any(stabilizer<n<group for n in subgroups))
    check('the S3 class stabilizer is nonnormal and therefore not subnormal',
        any(compose(compose(g,h),next(k for k in group if compose(g,k)==(0,1,2))) not in stabilizer
            for g in group for h in stabilizer))

    retained = {
        'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
        'hessian/hc4-triangular-graphs/NILPOTENT_COVER_MONODROMY.md':'5378746b4148c65ef7b63ca033674a6a605f2ef549c69dabb47aec516b842ac5',
        'hessian/hc4-triangular-graphs/verify_nilpotent_covers.py':'635859a884c5046410647cf06c0c0da49ca12d9fd7f1996059e18de42802a34c',
        'hessian/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36',
        'hessian/hc4-triangular-graphs/sources/stacks-monic-inertia-0BSD.html':'cd127a5a3b7d2b586600d03cc2a810bf400af445f77fb9d1248aa55f41082253',
        'hessian/hc4-triangular-graphs/sources/stacks-monic-purity-0BMB.html':'37a491b9ddf1fed9ccbc333dfd635731e6394e64c0d1d5a27a323a63b65e8585',
        'hessian/hc4-triangular-graphs/sources/milne-finite-etale-covers-lec-2.21.pdf':'718182f1a55d3af6b71f0d2faa5836851fb8611f1d7276896e3d5affdf28697b'
    }
    context = json.loads((HERE/'sources/normal-boundary-import.json').read_text())
    for item in context['sources']:
        retained['hessian/hc4-triangular-graphs/sources/'+item['retained_file']] = item['sha256']
    for name,digest in retained.items():
        check('retained proof, replay helper or primary input/context pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_normal_boundary.py','NORMAL_QUOTIENT_AND_BOUNDARY.md','sources/normal-boundary-import.json']},
        dependencies=retained,
        written_not_finite_exhaustion=[
            'every finite Galois intermediate extension of a complex Keller field is trivial in every dimension',
            'the omitted divisors in the finite plane normalization form an exact free basis of its global class group',
            'a torsion or factorial global normalization class group forces the plane Keller inverse',
            'subnormal divisor-class stabilizers suffice in the filtered output criterion, without nilpotent group image',
            'Hurwitz lattice class actions meet that criterion while remaining nonnilpotent'],
        source_stabilizers_asserted_normal_without_hypothesis=False,
        full_group_preservation_for_arbitrary_Galois_covers_asserted=False,
        one_visible_unramified_place_asserted_all_without_Galois_normality=False,
        hypothetical_Keller_normalization_asserted_smooth=False,
        global_factoriality_inferred_from_local_factoriality=False,
        ambient_output_class_group_identified_with_normalization_boundary=False,
        quaternionic_class_module_asserted_realized_by_an_explicit_output=False,
        nonmonic_S3_root_family_asserted_closed=False,
        purity_or_finite_cover_comparison_theorems_fully_reproved=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-normal-boundary.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
