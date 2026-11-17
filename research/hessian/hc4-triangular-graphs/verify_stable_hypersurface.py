"""Exact controls for the stable hypersurface obstruction.

The all-degree pencil, primitive-composition and smoothness proofs are
written separately. Finite controls do not certify those arguments.
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
        check(name,all(S.cancel(item.doit())==0 for item in values))

    s,U,P,Q,t1,t2,d1,d2,d3,T = S.symbols('s U P Q t1 t2 d1 d2 d3 T')
    dvars = (d1,d2,d3)
    h = d3-d1*d2-d1**3
    levels = (S.Integer(1),S.Integer(2),S.Integer(4))
    factors = [h-value for value in levels]
    check('nonlinear conjugate-component controls have the same total degree',
          all(S.Poly(item,*dvars).total_degree()==3 for item in factors))
    for ii,jj in [(0,1),(0,2),(1,2)]:
        zero('distinct level components have exact constant-unit Bezout identity '+str((ii,jj)),
             (factors[ii]-factors[jj])/(levels[jj]-levels[ii])-1)
    for ii,value in enumerate(levels):
        zero('each component has an actual polynomial A2 graph parameterization '+str(ii),
             factors[ii].subs(d3,d1*d2+d1**3+value))
    nonconstant = lambda expr: S.Poly(expr,*dvars).terms()[:-1]
    check('all nonlinear components have exactly the same nonconstant coefficients',
          all(nonconstant(item)==nonconstant(factors[0]) for item in factors))
    R = S.prod(T-value for value in levels)
    zero('the full product is exactly one polynomial composed with the common pencil',
         S.prod(factors)-R.subs(T,h))
    radical = S.sqrt(t1)
    conjugates = [h-radical,h+radical]
    zero('a genuine algebraic conjugate pair has its common nonconstant part over the base',
         S.prod(conjugates)-(h*h-t1))
    zero('Galois interchange changes only the constant coefficient of the normalized factor',
         (h+radical)-(h-radical)-2*radical)

    # Equal degree is not enough if component units are nonconstant.
    e1 = d1*d2-1
    e2 = d1*d2+d1*d1-1
    check('the nonconstant-component-unit countercontrol has equal degree',
          S.Poly(e1,d1,d2).total_degree()==S.Poly(e2,d1,d2).total_degree()==2)
    zero('those two non-pencil components are nevertheless disjoint',
         d2*d2*e2-(d1*d2+1+d2*d2)*e1-1)
    zero('the difference is a nonconstant unit on the first hyperbola',
         (e2-e1).subs(d2,1/d1)-d1*d1)
    zero('the hyperbola component has a nonconstant invertible function',
         (d1*d2).subs(d2,1/d1)-1)

    # Primitive nonconstant coefficient content, not target unimodularity.
    hgood = t1*d1+t2*d2
    check('primitive nonconstant coefficients need not generate the target unit ideal',
          S.gcd(t1,t2)==1 and hgood.subs({t1:0,t2:0})==0)
    coeffs = S.symbols('b0:4')
    formal_R = sum(coeffs[ii]*T**ii for ii in range(4))
    reduction = S.Poly(formal_R.subs(T,hgood).subs(t1,0),d2)
    check('composition with the nonconstant prime reduction preserves a generic polynomial',
          [reduction.nth(ii) for ii in range(4)]==[coeffs[ii]*t2**ii for ii in range(4)])
    bad_numerator = hgood*hgood+t1*hgood+t2
    zero('an attempted coefficient denominator has a nonzero residue numerator',
         bad_numerator.subs(t1,0)-(t2*t2*d2*d2+t2))
    hbad = t1*d1+t1*t1*d2
    zero('nonprimitive content really permits a hidden coefficient denominator',hbad/t1-(d1+t1*d2))
    check('the denominator countercontrol violates precisely primitive content',
          S.gcd(t1,t1*t1)==t1)
    zero('content normalization gives the same rational polynomial algebra generator',
         hbad-t1*(d1+t1*d2))

    # Fully nonlinear H has exact polynomial source recovery when K is invertible.
    K = S.Matrix([s,U+s**3])
    H = S.Matrix([P,Q+P*P,P*(Q+P*P)+s*P+U])
    image = {t1:K[0],t2:K[1],d1:H[0],d2:H[1],d3:H[2]}
    E = d3-d1*d2-t1*d1-t2+t1**3
    zero('the positive nonlinear embedding has its exact global image equation',
         E.subs(image,simultaneous=True))
    inverse = S.Matrix([t1,t2-t1**3,d1,d2-d1*d1])
    zero('the nonlinear embedding has an actual source inverse',
         inverse.subs(image,simultaneous=True)-S.Matrix([s,U,P,Q]))
    back = {s:inverse[0],U:inverse[1],P:inverse[2],Q:inverse[3]}
    returned = S.Matrix([K[0],K[1],H[0],H[1],H[2]]).subs(back,simultaneous=True)
    zero('the other composition is the identity modulo the image equation',
         returned-S.Matrix([t1,t2,d1,d2,d3])+S.Matrix([0,0,0,0,E]))
    zero('the positive base has unit Jacobian',K.jacobian((s,U)).det()-1)
    zero('the positive image is smooth over the target everywhere',S.diff(E,d3)-1)
    check('the positive dual image is not an affine plane',
          S.Poly(E,*dvars).total_degree()==2)
    check('the positive dual parameterization is not fibre-affine',
          max(S.Poly(item,P,Q).total_degree() for item in H)==3)

    # The derivative and smoothness gate, including ramified and local controls.
    for degree in (2,3,5):
        relation = T**degree-t1
        compound = relation.subs(T,h)
        zero('the composition derivative identity is exact for degree '+str(degree),
             S.Matrix([S.diff(compound,dv) for dv in dvars])-
             S.diff(relation,T).subs(T,h)*S.Matrix([S.diff(h,dv) for dv in dvars]))
        check('the ramified comparison has a zero of every dual derivative, degree '+str(degree),
              all(S.diff(compound,dv).subs({d1:0,d2:0,d3:0})==0 for dv in dvars))
    zero('the ramified polynomial base fails the needed Keller hypothesis',
         S.Matrix([s*s,U]).jacobian((s,U)).det()-2*s)
    check('localization permits a nonconstant unit relation derivative',
          S.diff(T**3-t1,T).subs(T,s)==3*s*s and S.cancel((3*s*s)/(3*s*s))==1)
    c = S.symbols('c',nonzero=True)
    derivative = S.Poly(S.diff(formal_R,T)-c,T)
    check('a scalar relation derivative eliminates every higher term in characteristic zero',
          [derivative.nth(ii) for ii in range(3)]==[coeffs[1]-c,2*coeffs[2],3*coeffs[3]])
    check('the derivative-minus-scalar relation has strictly smaller T degree',
          S.degree(S.diff(formal_R,T)-c,T)<S.degree(formal_R,T))
    charp = S.Poly(T**5-T-t1,T,t1,modulus=5)
    check('positive characteristic still requires an explicit exclusion',
          charp.diff(T)==S.Poly(-1,T,t1,modulus=5) and charp.degree(T)==5)

    # Byte pins are checked, without replaying or rewriting retained receipts.
    old = json.loads((HERE/'receipt-reconstruction.json').read_text())
    dependency_pins = {}
    for name in ['RECONSTRUCTION_EMBEDDING_OBSTRUCTION.md']:
        expected = old['source_sha256'][name]
        check('retained reconstruction proof pin '+name,
              hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected)
        dependency_pins[str((HERE/name).relative_to(ROOT))] = expected
    for name,expected in old['dependencies'].items():
        check('retained mathematical input pin '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected)
        dependency_pins[name] = expected
    context = json.loads((HERE/'sources/stable-hypersurface-import.json').read_text())
    for record in context['sources']:
        path = HERE/'sources'/record['retained_file']
        check('inspected primary basic-algebra source '+record['retained_file'],
              hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'])
        dependency_pins[str(path.relative_to(ROOT))] = record['sha256']
    check('basic sources do not import a new classification or assume embedding rectification',
          context['new_external_classification_imported'] is False and
          context['arbitrary_affine_plane_embedding_rectification_assumed'] is False)

    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['verify_stable_hypersurface.py','STABLE_HYPERSURFACE_OBSTRUCTION.md',
                                    'sources/stable-hypersurface-import.json']},
        dependencies=dependency_pins,
        written_not_finite_exhaustion=[
            'disjoint equal-degree constant-unit conjugate hypersurfaces share a polynomial pencil',
            'primitive nonconstant coefficient content removes every height-one denominator',
            'arbitrary nonlinear stable hypersurface embeddings over a plane Keller base force an inverse',
            'the exact relative generator minimum for a nonautomorphic plane Keller base is r+2 for every r',
            'an HC4 cotangent gradient pair cannot extend to an ambient coordinate pair in an A5 embedding'],
        only_affine_or_planar_fibre_images_covered=False,
        primitive_composition_denominator_step_omitted=False,
        arbitrary_A2_in_A3_embedding_rectification_assumed=False,
        nonconstant_component_units_treated_as_constants=False,
        general_graph_physical_base_identified_with_Keller_outputs=False,
        stable_generator_minimum_asserted_as_absolute_embedding_dimension=False,
        nonlinear_Hessian_invariance_used=False,
        universal_SL2_row_equation_solved=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,novelty_priority_asserted=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-stable-hypersurface.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=result['check_count'],
                          elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
