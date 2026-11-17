"""Exact controls for isotropic-ideal graphs and positive weighted images."""
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

    R,x,y,z,w,c = S.symbols('R x y z w c')
    q = R*R+x*y+z*w
    variables = (R,x,y,z,w)
    # Generic multi-isotropic shear, including its actual inverse.
    a,b = S.symbols('a b')
    h = x*a+z*b
    outR = R-h
    outy = y+2*R*a-h*a
    outw = w+2*R*b-h*b
    zero('the generic nonisotropic graph shear preserves the full quadratic',
        outR*outR+x*outy+z*outw-q)
    zero('its inverse recovers the normal coordinate',outR+h-R)
    zero('its inverse recovers the first dual coordinate',outy-2*outR*a-h*a-y)
    zero('its inverse recovers the second dual coordinate',outw-2*outR*b-h*b-w)
    zero('the entire graph goes to the specified hyperplane',outR.subs(R,c+h)-c)
    # Coefficients can depend on unchanged residual, nonisotropic variables.
    P = x+z*z+w**3
    hres = x*P
    transformed = S.Matrix([R-hres,x,y+2*R*P-hres*P,z,w])
    zero('arbitrary residual-variable dependence still preserves q',
        transformed[0]**2+transformed[1]*transformed[2]+transformed[3]*transformed[4]-q)
    zero('its full ambient Jacobian determinant is one',transformed.jacobian(variables).det()-1)
    zero('the nonlinear residual coefficient is fixed by the shear',
        P.subs(dict(zip(variables,transformed)),simultaneous=True)-P)
    zero('the residual-variable graph has literal source recovery',
        S.Matrix([transformed[0]+hres,
            transformed[2]-2*transformed[0]*P-hres*P])-S.Matrix([R,y]))
    # Isotropic normal: h may use the fixed normal mate and residual R.
    Piso = z*z+y**3+R*R
    hiso = z*Piso
    iso = S.Matrix([R,x-hiso,y,z,w+y*Piso])
    zero('the isotropic graph shear preserves the whole quadratic',
        iso[0]**2+iso[1]*iso[2]+iso[3]*iso[4]-q)
    zero('the isotropic shear has ambient Jacobian one',iso.jacobian(variables).det()-1)
    zero('its coefficient remains unchanged even with normal-mate dependence',
        Piso.subs(dict(zip(variables,iso)),simultaneous=True)-Piso)
    zero('the isotropic shear has an exact polynomial inverse',
        S.Matrix([iso[1]+hiso,iso[4]-y*Piso])-S.Matrix([x,w]))
    zero('the whole isotropic graph goes to x=c',iso[1].subs(x,c+hiso)-c)
    # Hyperplane types, used after the literal ambient inverse.
    check('the nonisotropic normal has nonzero inverse-metric norm',
        S.hessian(q,variables).inv()[0,0]!=0)
    check('the isotropic normal has zero inverse-metric norm',
        S.hessian(q,variables).inv()[1,1]==0)
    zero('the nonisotropic hyperplane restriction is Q4 plus its critical value',
        q.subs(R,c)-(x*y+z*w+c*c))
    check('the through-center isotropic restriction has rank three',
        S.hessian(q.subs(x,0),(R,y,z,w)).rank()==3)
    zero('that isotropic restriction has its entire critical line',
        S.Matrix([S.diff(q.subs(x,0),v) for v in (R,y,z,w)]).subs({R:0,z:0,w:0}))
    value = S.symbols('value')
    zero('an off-center isotropic hyperplane potential is an actual coordinate',
        q.subs({x:c,y:(value-R*R-z*w)/c},simultaneous=True)-value)

    def euler(poly,weights):
        return sum(weight*v*S.diff(poly,v) for v,weight in weights.items())

    # A compatible positive image whose defining graph has two nonlinear terms.
    weights = {R:3,x:1,y:5,z:2,w:4}
    hp = x**3+x*z
    ap = S.diff(hp,x)/3
    bp = 2*S.diff(hp,z)/3
    zero('positive ambient weights give the whole quadratic degree six',euler(q,weights)-6*q)
    zero('the image equation has the required middle weight three',euler(R-hp,weights)-3*(R-hp))
    zero('weighted Euler supplies the literal isotropic-ideal certificate',x*ap+z*bp-hp)
    positive = S.Matrix([R-hp,x,y+2*R*ap-hp*ap,z,w+2*R*bp-hp*bp])
    for i,(coordinate,weight) in enumerate(weights.items()):
        zero('the rectifying polynomial coordinate respects weight '+str(i),
            euler(positive[i],weights)-weight*positive[i])
    zero('the positive rectification still preserves the full potential',
        positive[0]**2+positive[1]*positive[2]+positive[3]*positive[4]-q)
    zero('the positive rectification has global Jacobian one',positive.jacobian(variables).det()-1)
    H = S.hessian(q,variables)
    low = S.Matrix([[0,1,0,0,0],[0,0,0,1,0]])
    zero('the lower-weight covectors form a totally isotropic support',low*H.inv()*low.T)
    zero('that support is inverse-metric orthogonal to the middle normal',
        low*H.inv()*S.Matrix([1,0,0,0,0]))
    native = q.subs(R,hp)
    zero('the compatible native graph has actual unit source Hessian',
        S.hessian(native,(x,y,z,w)).det()-1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    hy = hp.subs({x:p2,z:p4},simultaneous=True)
    inverse = {x:p2,z:p4,y:p1-2*hy*(3*p2*p2+p4),w:p3-2*hy*p2}
    zero('its actual positive-weight gradient has a polynomial inverse',
        S.Matrix([S.diff(native,v) for v in (x,y,z,w)]).subs(inverse,simultaneous=True)-S.Matrix([p1,p2,p3,p4]))
    # The coefficient space can distinguish weight mismatch from isotropy.
    check('a non-middle normal is isotropic in the inverse metric',H.inv()[1,1]==0 and weights[x]*2!=6)
    # Negative controls for dropping one of the two image hypotheses.
    Ebad = R-x*y
    badweights = {R:2,x:1,y:1,z:1,w:1}
    zero('the incompatible graph nevertheless has a positive homogeneous equation',euler(Ebad,badweights)-2*Ebad)
    check('the ambient quadratic is not homogeneous for those weights',euler(q,badweights)!=2*q and euler(q,badweights)!=4*q)
    fbad = q.subs(R,x*y)
    rho = S.symbols('rho',nonzero=True)
    badcurve = {x:rho,y:-1/(2*rho),z:0,w:0}
    zero('the incompatible positive graph has a whole full critical curve',
        S.Matrix([S.diff(fbad,v) for v in (x,y,z,w)]).subs(badcurve,simultaneous=True))
    zero('that curve and the origin have different critical values',
        fbad.subs(badcurve,simultaneous=True)-fbad.subs({x:0,y:0,z:0,w:0})+S.Rational(1,4))
    mixed = {R:0,x:1,y:-1,z:1,w:-1}
    Emixed = R-x*y+1
    zero('without positivity both image equation and quadratic can have degree zero',euler(Emixed,mixed))
    zero('the same nonpositive action preserves q',euler(q,mixed))
    check('this invariant graph omits the ambient critical point',Emixed.subs({R:0,x:0,y:0,z:0,w:0})==1)
    fmixed = q.subs(R,x*y-1)
    mixedcurve = {x:rho,y:1/(2*rho),z:0,w:0}
    zero('it too has a full critical curve',
        S.Matrix([S.diff(fmixed,v) for v in (x,y,z,w)]).subs(mixedcurve,simultaneous=True))
    zero('its two critical values also differ',
        fmixed.subs(mixedcurve,simultaneous=True)-fmixed.subs({x:0,y:0,z:0,w:0})+S.Rational(1,4))
    # Nonlinear source conjugacy does not preserve actual full Hessians.
    a,b,c0,d0 = S.symbols('a b c0 d0')
    q4 = (a*a+b*b+c0*c0+d0*d0)/2
    conjugated = q4.subs(a,a+b*b)
    zero('the source-conjugated quadratic records its nonunit Hessian explicitly',
        S.hessian(conjugated,(a,b,c0,d0)).det()-(1+2*a+2*b*b))
    check('the unaltered source quadratic has a unit Hessian',S.hessian(q4,(a,b,c0,d0)).det()==1)
    pins = {}
    for receipt_name,name in [
        ('receipt-quadratic-symmetry.json','QUADRATIC_ISOMETRY_SATURATION.md'),
        ('receipt-two-affine-geometry.json','TWO_AFFINE_ENDPOINT_GEOMETRY.md'),
        ('receipt-global-value.json','GLOBAL_VALUE_ORBITS.md'),
        ('receipt-ambient-quadric.json','AMBIENT_QUADRIC_GRAPHS.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    context = json.loads((HERE/'sources/positive-weight-images-context.json').read_text())
    pins[str((HERE/'sources'/context['retained_file']).relative_to(ROOT))] = context['sha256']
    for name,digest in pins.items():
        check('retained mathematical or primary context pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_positive_weight_images.py','POSITIVE_WEIGHT_IMAGE_RECTIFICATION.md','sources/positive-weight-images-context.json']},
        dependencies=pins,written_not_finite_exhaustion=[
            'unrestricted-degree isotropic-ideal graph shears, both normal types and literal inverses',
            'all smooth positive-weight invariant images compatible with a nondegenerate quadratic',
            'the isotropic normal excludes unit Hessian and the nonisotropic case constructs a grading-compatible ambient isometry',
            'arbitrary source parameterizations retain only the singleton endpoint for compatible positive images',
            'a positive linear weighting on actual source variables makes the formal gradient inverse polynomial in every dimension'],
        image_equivariance_asserted_actual_source_linear_equivariance=False,
        nonlinear_source_Hessian_invariance_used=False,
        arbitrary_image_rectification_or_group_generation_assumed=False,
        nonpositive_weight_images_asserted_classified=False,
        positive_image_without_quadratic_compatibility_asserted_classified=False,
        all_source_parameterizations_asserted_unit_Hessian=False,
        singleton_zero_endpoint_resolved=False,unit_ideal_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,new_external_classification_imported=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-positive-weight-images.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
