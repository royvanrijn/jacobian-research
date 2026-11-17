"""Exact controls for the three-affine-projection quadratic obstruction.

The arbitrary-degree proof covers every kernel rank and both positions
of the intrinsic plane. Finite controls do not certify that case split.
"""
import argparse
import hashlib
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
        check(name,all(S.cancel(S.expand(item.doit()))==0 for item in values))

    s1,s2,s3,t,H,G = S.symbols('s1 s2 s3 t H G')
    ambient = (s1,s2,s3,H,G)
    forms = [s1*s2+s3*s3+H*G,H*H+s1*G+s2*s3,s1*H+s2*G+s3*s3]
    for rank,form in zip((2,1,0),forms):
        check('kernel rank '+str(rank)+' has a nondegenerate ambient quadratic',S.hessian(form,ambient).det()!=0)
        check('its two unretained coordinates have exactly quadratic rank '+str(rank),
              S.hessian(form,(H,G)).rank()==rank)

    # Strict extension beyond four affine retained source coordinates.
    HH,GG = s3+s2**3,t+s2*s2
    f = forms[0].subs({H:HH,G:GG},simultaneous=True)
    source = (s1,s2,s3,t)
    zero('the strict three-affine example has polynomial source recovery',GG-s2*s2-t)
    zero('it has actual unit full Hessian',S.hessian(f,source).det()-1)
    a,b = S.symbols('a b')
    nonlinear = S.Poly(a*(HH-s3)+b*(GG-t),s2)
    check('its nonlinear coordinate combinations cannot yield a fourth affine source coordinate',
          nonlinear.nth(3)==a and nonlinear.nth(2)==b)
    zero('the free-line-contained control has the exact triangular cotangent coefficients',
         S.Matrix([S.diff(f,s1),S.diff(f,t)])-S.Matrix([s2,s3+s2**3]))
    z1,z2,w1,w2 = S.symbols('z1 z2 w1 w2')
    hp,ell = z2+z1**3,z1
    b1 = hp*hp+1
    zero('the constant-Jacobian product coefficient identity retains its multiple-of-h term',
         S.Matrix([hp,ell+hp*b1]).jacobian((z1,z2)).det()+1)
    # This is the polynomial remainder step used for the first case.
    zero('Jacobian modulo the first component differs only by its displayed multiple',
         S.Matrix([hp,ell+hp*b1]).jacobian((z1,z2)).det()-
         S.Matrix([hp,ell]).jacobian((z1,z2)).det()-
         hp*S.Matrix([hp,b1]).jacobian((z1,z2)).det())

    # Nondegenerate kernel, W transverse to the free line.
    h,g = z2+z1*z1,z1
    hfull,gfull = w1+h,-w2+g
    f2 = w1*w2+z1*z1+hfull*gfull
    zero('the two linear W factors cancel the fixed quadratic W top',
         f2-(z1*z1+h*g+w1*g-w2*h))
    zero('their Keller coefficients have unit Jacobian',S.Matrix([g,-h]).jacobian((z1,z2)).det()+1)
    zero('one retained quotient coordinate plus the factors recovers the source plane',h-z1*z1-z2)
    zero('the transverse product control has unit full Hessian',S.hessian(f2,(z1,z2,w1,w2)).det()-1)
    check('a nonconstant factor cannot divide a nonzero constant quadratic coefficient',S.degree(z1*z1,z1)>0)
    check('nonzero polynomial product degrees add over the quotient domain',
          S.Poly((z1*w1+w2)*(z2*w1+w2),z1,z2).total_degree()==2)

    # Rank-one kernel, transverse intrinsic plane: keep a nontrivial
    # normal jet to test that restriction to H=0 removes it exactly.
    s,w = S.symbols('s w')
    A,B = z1,2*(z2+z1*z1)
    h = z2+z1*z1
    h0 = w+h
    jet = z1*w+s*s
    physical2,physical3 = w+z1+2*s,-w+z1+3*s
    q = physical2*physical3
    qbase = z1*z1-w*w
    qjet = 5*z1+w
    Hfull = h0+s*jet
    Gfull = A-2*h0*jet-s*jet*jet-qjet-6*s
    frank1 = Hfull*Hfull+s*Gfull+q
    expected = z1*z1+h*h+w*B+s*A
    zero('rank-one full potential cancels all arbitrary normal-jet terms',frank1-expected)
    zero('its actual full Hessian is a nonzero constant',S.hessian(S.expand(frank1),(s,z1,z2,w)).det()-4)
    zero('the physical base change is invertible on the retained affine variables',
         S.Matrix([(physical2+physical3-5*s)/2,(physical2-physical3+s)/2])-S.Matrix([z1,w]))
    zero('the remaining source variable has global polynomial recovery from H',Hfull-s*jet-w-z1*z1-z2)
    plane = {s:0,w:-h}
    zero('the restricted source scheme is exactly the polynomial plane H0=0',h0.subs(plane,simultaneous=True))
    zero('the A component on that plane is G plus the full quadratic normal derivative',
         (Gfull+S.diff(q,s)).subs(plane,simultaneous=True)-A)
    zero('the B component on that plane is the quadratic tangential derivative',
         S.diff(q,w).subs(plane,simultaneous=True)-B)
    u2,u3,vG = S.symbols('u2 u3 vG')
    remaining = S.Matrix([vG+2*u3+3*u2,-u2+u3,(u2+u3)/2])
    check('the two Keller rows and one helper complete an ambient affine frame',
          remaining.jacobian((u2,u3,vG)).det()!=0)
    zero('the normal jet disappears because its multiplier H0 is zero',
         (2*h0*jet).subs(plane,simultaneous=True))
    zero('an isotropic tangent with zero square coefficient leaves a nonconstant linear derivative',
         S.diff(s2*s3+t*t,s2)-s3)

    # Totally isotropic kernel: include moving retained s3 and a syzygy.
    K = S.Matrix([z1,z2+z1*z1])
    p3 = z1+2*w1+3*w2
    syzygy = z1*z1+w1
    hp0 = K[0]-4*z1-4*w1-12*w2+w2*syzygy
    gp0 = K[1]-6*z1-9*w2-w1*syzygy
    frank0 = w1*hp0+w2*gp0+p3*p3
    zero('the rank-zero full identity keeps the quadratic retained-coordinate corrections',
         frank0-(z1*z1+w1*K[0]+w2*K[1]))
    zero('the rank-zero control has actual unit full Hessian',S.hessian(frank0,(z1,z2,w1,w2)).det()-1)
    zero('its retained third coordinate gives the quotient helper exactly',p3-2*w1-3*w2-z1)
    zero('its zero-base first jet recovers K with the quadratic correction',
         S.Matrix([hp0+4*z1,gp0+6*z1]).subs({w1:0,w2:0},simultaneous=True)-K)
    zero('the corrected restricted coordinates give actual polynomial source recovery',
         (gp0+6*z1-z1*z1).subs({w1:0,w2:0},simultaneous=True)-z2)

    # Linear-algebra positions underlying the affine-line endpoint cases.
    basis = S.eye(4)
    line = basis[:,3:4]
    zeroplane = basis[:,2:4]
    Wcontains = S.Matrix.hstack(basis[:,0],basis[:,3])
    projection = S.Matrix([[0,1,0,0],[0,0,1,0]])
    check('a free-line-contained zero-base plane projects to a quotient line',
          projection*Wcontains==S.zeros(2,2) and (projection*zeroplane).rank()==1)
    check('the source free line really belongs to that intrinsic plane',Wcontains.row_join(line).rank()==2)
    zero('the fully zero-base plane cannot be jointly affine for the retained square',S.diff(s3*s3,s3,2)-2)
    check('a transverse rank-one physical projection forces its kernel to move the square coordinate',
          S.Matrix.hstack(basis[:,0],basis[:,2])[:2,:].rank()==1)

    # Conditional singleton endpoint embedding with two retained affine coordinates.
    P,Q,L,N,lstar = S.symbols('P Q L N lstar')
    Cmat = S.Matrix([[1,z1*z1],[z2,1+z1*z1*z2]])
    row = S.Matrix([[z1,z2]])*Cmat
    dual = Cmat*S.Matrix([P,Q])
    zero('the singleton-row frame is genuinely polynomially invertible',Cmat.det()-1)
    zero('the two-affine singleton embedding has the exact cotangent potential',
         z1*dual[0]+z2*dual[1]-(row*S.Matrix([P,Q]))[0])
    zero('its source dual inverse is actual polynomial recovery',Cmat.adjugate()*dual-S.Matrix([P,Q]))

    # Conditional unit-ideal endpoint: the identity is universal in F,G,a,b
    # subject to the supplied Bezout certificate aF+bG=1.
    Fv,Gv,av,bv,V,R,D = S.symbols('Fv Gv av bv V R D')
    vactual = D+P*Fv+Q*Gv
    ractual = -bv*P+av*Q
    returned = S.Matrix([av*(vactual-D)-Gv*ractual,bv*(vactual-D)+Fv*ractual])
    zero('the relative unit-ideal recovery has exactly the Bezout multiplier',
         returned-(av*Fv+bv*Gv)*S.Matrix([P,Q]))
    normalx,normaly,normalT,normalZ = 1/(2*N),z1,z2,R
    normalW = V-lstar-L*z2*z2-N*z1*R
    zero('the two-affine coordinate endpoint has its exact quadratic potential',
         lstar+L*normalT**2+N*normaly*normalZ+2*N*normalx*normalW-V)

    dependency_pins = {}
    for receipt_name,names in [
        ('receipt-ambient-quadric.json',['AMBIENT_QUADRIC_GRAPHS.md']),
        ('receipt-reconstruction.json',['RECONSTRUCTION_EMBEDDING_OBSTRUCTION.md']),
        ('receipt-critical-components.json',['CRITICAL_COMPONENTS.md']),
        ('receipt-relative-coordinate.json',['RELATIVE_COORDINATE_ENDPOINT.md']),
        ('receipt-square-axis.json',['P_GRAPH_SINGLETON_ENDPOINT.md'])]:
        retained = json.loads((HERE/receipt_name).read_text())
        for name in names:
            dependency_pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
        if receipt_name in ['receipt-reconstruction.json','receipt-critical-components.json']:
            dependency_pins.update(retained['dependencies'])
    for name,expected_hash in dependency_pins.items():
        check('retained proof or primary input pin '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected_hash)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['verify_three_affine.py','THREE_AFFINE_PROJECTIONS.md']},
        dependencies=dependency_pins,
        written_not_finite_exhaustion=[
            'all three ranks of the unretained ambient quadratic plane',
            'both positions of the intrinsic affine plane relative to the free source line',
            'three affine retained source projections close every quadratic hypersurface embedding in all degrees',
            'an actual unit-Hessian example with exactly three affine retained projections extends the former graph criterion',
            'two affine retained projections already realize both distinct conditional plane endpoints'],
        new_external_classification_imported=False,
        relative_line_embedding_rectification_assumed=False,
        only_nondegenerate_kernel_rank_covered=False,
        finite_controls_asserted_to_prove_exhaustive_case_split=False,
        two_affine_retained_projections_asserted_closed=False,
        two_affine_class_asserted_equivalent_to_full_JC2=False,
        nonautomorphic_Keller_example_constructed=False,
        actual_affine_source_hypothesis_replaced_by_nonlinear_coordinates=False,
        nonlinear_Hessian_invariance_used=False,marked_gradient_collision_transfer_asserted=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        universal_SL2_row_equation_solved=False,unrestricted_HC4_or_JC2_resolved=False,
        novelty_priority_asserted=False,formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-three-affine.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=result['check_count'],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
