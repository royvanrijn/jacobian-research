"""Exact controls for written all-degree quadratic-symmetry saturation."""
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

    x,y,z,w,t,h = S.symbols('x y z w t h')
    ambient = (x,y,z,w,t)

    def quadratic(values):
        xx,yy,zz,ww,tt = values
        return xx*yy+zz*ww+tt*tt

    def shear(values,index,coefficient):
        values = list(values)
        mate = {0:1,1:0,2:3,3:2}[index]
        fixed,old_t = values[index],values[4]
        values[mate] -= 2*old_t*coefficient+fixed*coefficient**2
        values[4] += fixed*coefficient
        return tuple(values)

    for index in range(4):
        changed = shear(ambient,index,h)
        zero('generic elementary symmetry '+str(index)+' preserves the entire quadratic',
            quadratic(changed)-quadratic(ambient))
        zero('generic elementary symmetry '+str(index)+' has its literal polynomial inverse',
            S.Matrix(shear(changed,index,-h))-S.Matrix(ambient))
        zero('generic elementary symmetry '+str(index)+' has unit ambient Jacobian',
            S.Matrix(changed).jacobian(ambient).det()-1)
    coefficient = x*x+z*w+z**3
    changed = shear(ambient,0,coefficient)
    reverse_coefficient = changed[0]**2+changed[2]*changed[3]+changed[2]**3
    zero('a fixed-coordinate polynomial coefficient remains unchanged',reverse_coefficient-coefficient)
    zero('the variable-coefficient inverse is a genuine polynomial automorphism',
        S.Matrix(shear(changed,0,-reverse_coefficient))-S.Matrix(ambient))
    zero('variable-coefficient symmetry preserves the full potential',
        quadratic(changed)-quadratic(ambient))
    zero('variable-coefficient symmetry still has unit ambient Jacobian',
        S.Matrix(changed).jacobian(ambient).det()-1)
    for index,coefficient in [(0,z),(1,y),(2,z),(3,w)]:
        changed = shear(ambient,index,coefficient)
        zero('four-word step '+str(index)+' preserves the entire potential',
            quadratic(changed)-quadratic(ambient))
        zero('four-word step '+str(index)+' inverts with its fixed coefficient',
            S.Matrix(shear(changed,index,-coefficient))-S.Matrix(ambient))
    a,b,c,d,T = S.symbols('a b c d T')
    source = (a,b,c,d)
    f = a*b+c*d
    zero('the zero-affine-pullback family has actual unit source Hessian',S.hessian(f,source).det()-1)
    zero('its source gradient has an explicit polynomial inverse',
        S.Matrix([S.diff(f,b),S.diff(f,a),S.diff(f,d),S.diff(f,c)])-S.Matrix(source))

    # Small degree controls supplement the written arbitrary-integer proof.
    for m,n,r in [(1,1,1),(2,1,2),(1,3,1)]:
        A0 = 2*m+1
        B0 = (n+1)*A0
        degrees = [A0,B0+1,(2*n+1)*A0,(r+1)*(B0+1),(2*r+1)*(B0+1)]
        check('strict degree ordering in control '+str((m,n,r)),
            degrees==sorted(set(degrees)) and degrees[0]>1)
        check('unique leading summands in control '+str((m,n,r)),
            B0>m+1 and (2*n+1)*A0>m+1+n*A0 and B0+1>3 and
            (r+1)*(B0+1)>B0 and (2*r+1)*(B0+1)>B0+r*(B0+1))
    # An actual one-variable restriction sees every top term for m=n=r=1;
    # these lower bounds equal the proven multivariate upper bounds.
    zeta = S.symbols('zeta')
    values = (zeta,zeta,zeta,zeta,S.Integer(0))
    for index in [0,1,2,3]:
        coefficient = values[2] if index in [0,2] else values[index]
        values = tuple(S.expand(v) for v in shear(values,index,coefficient))
    check('the explicit nonisotropic control has degrees 9,3,21,7,14',
        [S.degree(v,zeta) for v in values]==[9,3,21,7,14])
    zero('its full word retains the quadratic identity',quadratic(values)-2*zeta*zeta)
    recovered = values
    for index in [3,2,1,0]:
        coefficient = recovered[2] if index in [0,2] else recovered[index]
        recovered = tuple(S.expand(v) for v in shear(recovered,index,-coefficient))
    zero('the actual inverse word recovers every source coordinate',
        S.Matrix(recovered)-S.Matrix([zeta,zeta,zeta,zeta,0]))
    for epsilon in [0,1]:
        values = (S.Integer(epsilon),zeta,zeta,zeta,zeta)
        for index in [0,1,2,3]:
            coefficient = values[2] if index in [0,2] else values[index]
            values = tuple(S.expand(v) for v in shear(values,index,coefficient))
        check('isotropic offset '+str(epsilon)+' has degrees 6,2,15,5,10',
            [S.degree(v,zeta) for v in values]==[6,2,15,5,10])
        zero('isotropic offset '+str(epsilon)+' retains its exact potential',
            quadratic(values)-(epsilon*zeta+2*zeta*zeta))
    iso_source = (b,c,d,T)
    zero('the critical-line model retains its degenerate Hessian',
        S.hessian(c*d+T*T,iso_source).det())
    check('the coordinate model is not asserted to have unit Hessian in these parameters',
        S.diff(b+c*d+T*T,b)==1 and S.hessian(b+c*d+T*T,iso_source).det()==0)
    scale = S.symbols('scale',nonzero=True)
    zero('a constant similitude scales the whole quadratic',
        quadratic((scale**2*x,y,scale**2*z,w,scale*t))-scale**2*quadratic(ambient))
    zero('actual source Hessian determinants scale by the fourth power',
        S.hessian(scale*f,source).det()-scale**4)
    check('a nonlinear source change need not preserve the actual Hessian',
        S.hessian(a+c*d,source)!=S.hessian(a+b*b+c*d,source))
    pins = {}
    for receipt_name,name in [
        ('receipt-three-affine.json','THREE_AFFINE_PROJECTIONS.md'),
        ('receipt-two-affine-geometry.json','TWO_AFFINE_ENDPOINT_GEOMETRY.md'),
        ('receipt-ambient-quadric.json','AMBIENT_QUADRIC_GRAPHS.md'),
        ('receipt-global-value.json','GLOBAL_VALUE_ORBITS.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    for name,digest in pins.items():
        check('retained proof pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ['verify_quadratic_symmetry.py','QUADRATIC_ISOMETRY_SATURATION.md']},
        dependencies=pins,written_not_finite_exhaustion=[
            'all polynomial quadratic-similitude saturations of the closed actual-source classes',
            'all hyperplane-image saturations retain exactly the previous endpoint classification',
            'any polynomial reparameterization of q induced by an ambient automorphism is a similitude',
            'arbitrary positive m,n,r in the four-word family have no affine pullbacks'],
        explicit_nonisotropic_control_degrees=[9,3,21,7,14],
        explicit_isotropic_control_degrees=[6,2,15,5,10],
        nonlinear_Hessian_invariance_used=False,
        elementary_words_asserted_to_generate_all_symmetries=False,
        all_embeddings_asserted_rectifiable=False,
        arbitrary_source_reparameterizations_asserted_to_preserve_three_affine_gate=False,
        raw_nonlinear_image_asserted_new_mechanism=False,
        unit_ideal_Keller_pair_constructed=False,singleton_zero_endpoint_resolved=False,
        unit_ideal_endpoint_resolved=False,entire_two_affine_class_asserted_closed=False,
        unrestricted_HC4_or_JC2_resolved=False,new_external_classification_imported=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-quadratic-symmetry.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
