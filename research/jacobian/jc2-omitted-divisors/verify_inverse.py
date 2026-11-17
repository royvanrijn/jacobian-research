"""Source-finite inverse equations, Newton edges and exact fibre controls.

The general normalization and flat-rank arguments are written in
ADMISSIBLE_INVERSE.md. No unseen Keller polynomial or braid is computed.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from verify import quotient_length, checks as local_checks

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def zero(name, expr):
    check(name, S.cancel(expr) == 0)


def valuation(expr, prime, variables):
    if expr == 0:
        return None
    p, h = S.Poly(expr,*variables,domain=S.QQ), S.Poly(prime,*variables,domain=S.QQ)
    order = 0
    while True:
        quotient, remainder = p.div(h)
        if not remainder.is_zero:
            return order
        p = quotient
        order += 1


def newton_polygon(poly, prime, variables):
    vals = [valuation(poly.nth(i),prime,variables) for i in range(poly.degree()+1)]
    hull = []
    for i, order in enumerate(vals):
        if order is None:
            continue
        point = (i,order)
        while len(hull) >= 2:
            a,b = hull[-2:]
            cross = (b[0]-a[0])*(point[1]-b[1])-(b[1]-a[1])*(point[0]-b[0])
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    edges = []
    for a,b in zip(hull,hull[1:]):
        width = b[0]-a[0]
        slope = S.Rational(b[1]-a[1],width)
        edges.append(dict(start=list(a),end=list(b),width=width,slope=str(slope),escaping=bool(slope>0)))
    return dict(coefficient_valuations=vals,vertices=[list(p) for p in hull],edges=edges)


def root_fibre(P,t):
    degree = int(S.degree(P,t))
    simple = 0
    omitted = []
    for factor,k in S.factor_list(P,t)[1]:
        root_count = int(S.degree(factor,t))
        if k == 1:
            simple += root_count
        else:
            omitted.extend([2*k]*root_count)
    infinity = 4-degree
    if infinity == 1:
        simple += 1
    elif infinity >= 2:
        omitted.append(2*infinity)
    zero_root = int(P.subs(t,0) == 0)
    omitted.extend([1]*zero_root)
    return 2*simple-zero_root,sorted(omitted)


def main():
    started = time.monotonic()
    t,T,q,r,b = S.symbols('t T q r b')
    P = q*t**4+r*t**3+t-b
    D = S.diff(P,t)
    marking_elimination = (t**2-D)**2-D*(T*t-1)**2
    raw = S.Poly(S.resultant(P,marking_elimination,t),T)
    content = S.gcd_list(raw.all_coeffs())
    zero('exact resultant coefficient content',content-b*q**2)
    R = S.Poly(S.cancel(raw.as_expr()/content),T)
    delta = S.discriminant(P,t)
    check('primitive inverse degree8',R.degree() == 8)
    zero('leading coefficient records both images',R.LC()-b*delta)
    check('degree8 inverse monomial count',len(S.Poly(R.as_expr(),T,q,r,b).terms()) == 218)
    check('discriminant image irreducible',len(S.factor_list(delta,q,r,b)[1]) == 1)
    G = S.groebner(R.all_coeffs(),q,r,b,domain=S.QQ)
    check('coefficient ideal unit, no vertical target fibre',list(G) == [1])
    witness = S.Poly(R.as_expr().subs({q:1,r:0,b:1}),T)
    check('degree8 squarefree primitive-element witness',witness.degree()==8 and S.gcd(witness,witness.diff()).degree()==0)
    zero('degree-drop line',R.as_expr().subs({q:0,r:0})-(-T+b)*(b*T+b**2-2))
    zero('singleton target inverse exactly2T',R.as_expr().subs({q:0,r:0,b:0})-2*T)
    newton = {}
    for name,prime,width,slope in [('zero hyperplane',b,1,'1'),('binary discriminant',delta,4,'1/4')]:
        polygon = newton_polygon(R,prime,(q,r,b))
        escaping = [edge for edge in polygon['edges'] if edge['escaping']]
        check(name+' one escaping Newton edge',len(escaping)==1)
        check(name+' escaping horizontal length',escaping[0]['width']==width)
        check(name+' escaping valuation slope',escaping[0]['slope']==slope)
        check(name+' leading pole-pushforward multiplicity',polygon['coefficient_valuations'][-1]==1)
        newton[name] = polygon
    targets = [
        ('generic quartic',{q:1,r:0,b:1},8,[]),
        ('zero root simple',{q:1,r:0,b:0},7,[1]),
        ('simple infinity',{q:0,r:1,b:1},8,[]),
        ('double finite with simple infinity',{q:0,r:-S.Rational(1,3),b:S.Rational(2,3)},4,[4]),
        ('triple infinity',{q:0,r:0,b:2},2,[6]),
        ('singleton zero',{q:0,r:0,b:0},1,[1,6]),
        ('double root plus zero',{q:2,r:-3,b:0},3,[1,4]),
        ('triple finite root',{q:S.Rational(1,2),r:-1,b:S.Rational(1,2)},2,[6]),
        ('two double finite roots',{q:1,r:-2,b:-S.Rational(1,4)},0,[4,4]),
    ]
    target_rows = []
    for label,value,expected_retained,expected_missing in targets:
        fibre_poly = R.as_expr().subs(value)
        check(label+' no identically zero specialization',fibre_poly != 0)
        inverse_degree = int(S.degree(fibre_poly,T))
        retained,missing = root_fibre(P.subs(value),t)
        check(label+' exact inverse degree',inverse_degree==retained==expected_retained)
        check(label+' complete omitted local lengths',missing==expected_missing)
        check(label+' full normalization length',retained+sum(missing)==8)
        target_rows.append(dict(label=label,target={str(k):str(v) for k,v in value.items()},
                                inverse_polynomial=str(S.factor(fibre_poly)),retained_length=retained,
                                omitted_local_lengths=missing,total_length=8))
    all_root_partitions = [(1,1,1,1),(2,1,1),(2,2),(3,1),(4,)]
    possible = {2*part.count(1)-is_zero for part in all_root_partitions
                for is_zero in (0,1) if is_zero <= part.count(1)}
    check('complete degree8 cardinality upper list from projective root partitions',
          possible == {0,1,2,3,4,7,8})
    check('every listed degree8 cardinality has an exact target witness',
          {row['retained_length'] for row in target_rows} == possible)
    collision = S.Poly(R.as_expr().subs({q:0,r:0,b:1}),T)
    zero('two retained points can share one inverse root',collision.as_expr()+(T-1)**2)
    x,y,z3,z4 = S.symbols('x y z3 z4')
    p = 1+x*y
    Phi = ((1-y*t)**3*(p*t-x)+(p*t-x)**2
           *(4*y*(y*t+S.Rational(1,2))+z3*(p*t+x/2)+z4*t**2))
    outputs = S.Matrix([S.expand(Phi).coeff(t,3),S.expand(Phi).coeff(t,4),
                        S.expand(Phi).coeff(t,2),-S.expand(Phi).coeff(t,0)])
    J = outputs.jacobian([x,y,z3,z4])
    for label,point in [('ordinary',{x:1,y:0,z3:0,z4:0}),('escaping',{x:-1,y:2,z3:12,z4:-8})]:
        zero(label+' exact target t-1',Phi.subs(point)-(t-1))
        zero(label+' same T value1',(x+y).subs(point)-1)
        zero(label+' etale source point',J.subs(point).det()-1)
    # A source-finite nonprimitive coordinate scales aggregate degree by
    # k. This finite non-Keller map also tests a reduced-root collision.
    u,v = S.symbols('u v')
    bad = S.Poly(u*T-v,T)
    check('gcd-one coefficients alone do not exclude a vertical fibre',
          S.gcd_list(bad.all_coeffs())==1 and bad.as_expr().subs({u:0,v:0})==0)
    check('bad coefficient ideal is not the unit ideal',list(S.groebner(bad.all_coeffs(),u,v))!=[1])
    check('rejected birational plane control is not quasifinite',
          not S.groebner([x,x*y],x,y).is_zero_dimensional)
    finite_R = S.Poly(T**4-2*(u+v)*T**2+(u-v)**2,T)
    zero('finite plane inverse identity',finite_R.as_expr().subs({u:x**2,v:y**2,T:x+y}))
    beta=S.symbols('beta')
    zero('monic source elimination gives primitive inverse',
         S.resultant((T-beta)**2-u,beta**2-v,beta)-finite_R.as_expr())
    source_matrix=S.Matrix([[u-v-T**2,2*T*(u-T**2)],[2*T,3*T**2+u-v]])
    zero('source multiplication determinant',source_matrix.det()-finite_R.as_expr())
    smith=smith_normal_form(source_matrix.subs({u:1,v:1}),domain=S.QQ.poly_ring(T))
    smith_factors=[S.Poly(smith[i,i],T).monic().as_expr() for i in range(2)]
    check('source Smith factors at primitive collision',smith_factors==[T,T**3-4*T])
    check('source Smith degrees count four retained lengths',sum(S.degree(f,T) for f in smith_factors)==4)
    collision_at_11 = S.Poly(finite_R.as_expr().subs({u:1,v:1}),T)
    check('four-sheet inverse can have only three distinct roots',
          S.degree(S.div(collision_at_11.as_expr(),S.gcd(collision_at_11.as_expr(),S.diff(collision_at_11.as_expr(),T)),T)[0],T)==3)
    four = quotient_length([x**2-1,y**2-1],(x,y))
    nonreduced_zero = quotient_length([x**2,y**2],(x,y))
    check('nonprimitive source-finite coordinate needs factor k2',four==2*S.degree(T**2-1,T)==4)
    check('one geometric source point need not have length1',nonreduced_zero==4)
    for xx in (-1,1):
        for yy in (-1,1):
            check('finite collision points etale '+str((xx,yy)),4*xx*yy!=0)
    # Complete non-Keller plane controls, all with source A2 and a single
    # omitted e1 prime. The symbolic proof covers every integer a>=1.
    u,v,z,s = S.symbols('u v z s')
    plane_controls = []
    for a in range(1,5):
        inv = S.Poly(u**a*T**2-T-v,T)
        zero(f'a{a} retained ramification image',
             (1+4*u**a*v).subs({u:x,v:x**a*y**2-y})-(2*x**a*y-1)**2)
        zero(f'a{a} escaping path',(x**a*y**2-y).subs({x:s,y:s**(-a)}))
        zero(f'a{a} infinity-chart derivative unit',S.diff(u**a-z-v*z**2,z).subs(z,0)+1)
        check(f'a{a} origin ideal',(S.groebner([x,x**a*y**2-y],x,y).polys)==S.groebner([x,y],x,y).polys)
        retained_zero = quotient_length([x,x**a*y**2-y],(x,y))
        omitted_zero = quotient_length([u,u**a-z],(u,z))
        ramified_value = quotient_length([x-1,x**a*y**2-y+S.Rational(1,4)],(y,x))
        check(f'a{a} origin lengths1+1',retained_zero==omitted_zero==1)
        check(f'a{a} retained ramification length2',ramified_value==2)
        polygon = newton_polygon(inv,u,(u,v))
        edge = [e for e in polygon['edges'] if e['escaping']]
        check(f'a{a} one omitted sheet despite pole order a',len(edge)==1 and edge[0]['width']==1 and edge[0]['slope']==str(a))
        W=S.symbols('W')
        bad_primitive=W**2-W-u**a*v
        zero(f'a{a} other primitive hides omission',bad_primitive.subs({u:x,v:x**a*y**2-y,W:x**a*y}))
        zero(f'a{a} hidden primitive regular on omitted chart',
             bad_primitive.subs(W,1+v*z)-v*(z+v*z**2-u**a))
        zero(f'a{a} escaping path keeps hidden primitive finite',(x**a*y).subs({x:s,y:s**(-a)})-1)
        check(f'a{a} hidden primitive leading coefficient unit',S.Poly(bad_primitive,W).LC()==1)
        plane_controls.append(dict(a=a,inverse=str(inv.as_expr()),leading_coefficient=str(u**a),
                                   omitted_prime=dict(image='u=0',ramifies=False,omitted=True,e=1,f=1,pole_order=a),
                                   retained_ramification=dict(image=str(1+4*u**a*v),ramifies=True,omitted=False,e=2,f=1),
                                   origin_lengths=[1,1],retained_ramified_fibre_length=2,Keller=False,newton=polygon,
                                   source_coordinate_fibration=dict(finite_T_fibre='A1',infinity_fibre='a*D',multiplicity=a,
                                                                     distinct_from_target_radial_pencil=True),
                                   inadmissible_primitive=dict(function='x^a*y',equation=str(bad_primitive),
                                                              source_finite=False,value_on_omitted_prime=1,
                                                              leading_coefficient=1),
                                   complete_target_strata=[
                                       dict(condition='u!=0 and discriminant!=0',retained_lengths=[1,1],omitted_lengths=[]),
                                       dict(condition='discriminant=0',retained_lengths=[2],omitted_lengths=[]),
                                       dict(condition='u=0',retained_lengths=[1],omitted_lengths=[1])]))
    equations = dict(as_of='2026-11-14',degree8_slice='n4 marked-root lambda2=0; not a plane Keller map',
                     target_polynomial=str(P),marking_elimination=str(S.expand(marking_elimination)),
                     raw_content=str(content),inverse_polynomial=str(R.as_expr()),
                     coefficients_low_to_high=[str(R.nth(i)) for i in range(9)],
                     leading_coefficient=str(S.factor(R.LC())),binary_discriminant=str(delta),
                     coefficient_ideal_basis=['1'],newton_polygons=newton,exceptional_targets=target_rows,
                     degree8_occurring_cardinalities=sorted(possible),
                     degree8_all_fibre_degree_strata='Retained length j iff a_8=...=a_(j+1)=0 and a_j!=0; omitted length 8-j',
                     nonKeller_plane_controls=plane_controls,
                     actual_r7_inverse_polynomial=None,actual_singleton_Keller_inverse_polynomial=None)
    (HERE/'inverse_equations.json').write_text(json.dumps(equations,indent=2)+'\n')
    manifest=json.loads((HERE/'sources/manifest.json').read_text())
    for source in manifest['sources']:
        check('primary source pin '+source['file'],
              hashlib.sha256((HERE/'sources'/source['file']).read_bytes()).hexdigest()==source['sha256'])
    inputs = ['ADMISSIBLE_INVERSE.md','verify.py','sources/manifest.json']
    receipt = dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),local_staircase_checks=len(local_checks),
                   checks=checks,runtime_seconds=round(time.monotonic()-started,6),python=platform.python_version(),sympy=S.__version__,
                   checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   equations_sha256=hashlib.sha256((HERE/'inverse_equations.json').read_bytes()).hexdigest(),
                   input_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs},
                   general_inverse_theorem='Written source-finiteness, normalization and flat-rank proof, not proved by finite examples',
                   complete_control_ledgers=True,global_r7_images_determined=False,singleton_zero_properness_proved=False,
                   goal_complete=False,formal_verification=False,external_review=False)
    (HERE/'inverse_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','local_staircase_checks','runtime_seconds']}))


if __name__=='__main__':
    main()
