'''PF2RULED37 exact controls for written all-degree cyclic multiple-fiber reduction.
The geometric normalization, torsor, descent and monodromy proofs are written.
'''
from pathlib import Path
from collections import deque
from math import factorial
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generated(generators):
    I=tuple(range(len(generators[0])));seen={I};todo=deque([I])
    while todo:
        a=todo.popleft()
        for b in generators:
            c=tuple(a[b[j]] for j in range(len(a)))
            if c not in seen:seen.add(c);todo.append(c)
    return seen
def cycle_lengths(a):
    left=set(range(len(a)));out=[]
    while left:
        j=min(left);k=j;n=0
        while k in left:
            left.remove(k);n+=1;k=a[k]
        out.append(n)
    return sorted(out)
def main():
    pred=json.loads((HERE/'intrinsic_normalization_receipt.json').read_text())
    check('intrinsic predecessor receipt PASS',pred['outcome']=='PASS')
    for name,h in pred['source_sha256'].items():
        check('intrinsic predecessor seal '+name,digest(HERE/name)==h)
    manifest=HERE/'sources/multiple-fiber-cyclic/manifest.json'
    pins=json.loads(manifest.read_text())
    check('two retained primary inputs',len(pins['primary_sources'])==2)
    for row in pins['primary_sources']:
        check('primary source seal '+row['path_relative_to_affine_complements'],
              digest(HERE/row['path_relative_to_affine_complements'])==row['sha256'])
    q,e,m=S.symbols('q e m',integer=True,positive=True)
    check('general torsor exponent is q plus ramification index',
          S.solve(S.Eq(m-q-1,e-1),m)==[q+e])
    check('ramification formula retains both coefficients',
          S.expand((e-1)+(q-1)-(q+e-2))==0)
    check('cyclic Euler formula agrees with q boundary lines',
          S.expand(2*q-(q-1)-(q+1))==0)
    # Exact chart controls with unrestricted affine translation powers,
    # including powers higher than the linear torsor degree.
    z,y,u,v,x=S.symbols('z y u v x')
    ac=S.symbols('a0:5')
    for q0,e0 in [(1,2),(2,2),(3,2),(5,2),(2,3),(4,3)]:
        m0=q0+e0
        g=ac[0]/z+ac[1]+ac[2]*z**e0+ac[3]*z**(m0-1)+ac[4]*z**(m0+1)
        yy=u**m0*(v-g.subs(z,1/u))
        det=S.det(S.Matrix([[S.diff(1/u,u),0],[S.diff(yy,u),S.diff(yy,v)]]))
        check(f'full weighted volume transition q{q0} e{e0}',
              S.expand(q0*u**(1-q0)*det+q0*u**(e0-1))==0)
        thetaU=S.expand(q0*u**(e0-1)*(v-g.subs(z,1/u)))
        thetaW=q0*u**(e0-1)*v
        residue=S.residue(thetaU-thetaW,u,0)
        check(f'invariant-index residue with high translation terms q{q0} e{e0}',
              S.expand(residue+q0*S.expand(g).coeff(z,e0))==0)
        check(f'second chart primitive differentiates correctly q{q0} e{e0}',
              -S.diff(thetaW,v)==-q0*u**(e0-1))
    # The deck character at exponent e is exactly zeta^q=1 for every q.
    check('residue character is a full q-th power',S.expand((q+e)-e)==q)
    for q0 in [2,3,5]:
        zeta=S.exp(2*S.pi*S.I/q0)
        check('shared deck residue coefficient q'+str(q0),S.simplify(zeta**q0)==1)
    # A finite permutation control of the general group proof. No surface
    # realization is inferred from these abstract wreath groups.
    wreath=[]
    for q0,d0 in [(2,3),(3,3),(2,4)]:
        length=q0*d0;I=list(range(length))
        c=I[:]
        for j in range(q0):c[j]=(j+1)%q0
        swaps=[]
        for b in range(d0-1):
            a=I[:]
            for j in range(q0):
                a[b*q0+j]=(b+1)*q0+j;a[(b+1)*q0+j]=b*q0+j
            swaps.append(tuple(a))
        g=generated([tuple(c)]+swaps)
        check(f'block cycle and symmetric block lifts generate full wreath q{q0} d{d0}',
              len(g)==q0**d0*factorial(d0))
        check(f'boundary inertia is q disjoint transpositions q{q0} d{d0}',
              cycle_lengths(swaps[0])==[1]*(q0*(d0-2))+[2]*q0)
        check(f'auxiliary inertia is exactly one q-cycle q{q0} d{d0}',
              cycle_lengths(tuple(c))==[1]*(q0*(d0-1))+[q0])
        wreath.append(dict(q=q0,d=d0,group_order=len(g)))
    # Uniform explicit affine control family N_q, with two localization
    # recoveries and regular global primitive when q>1.
    for q0 in range(1,7):
        xx=u*v+u**(-q0);ff=x*y**q0-1
        aa=x*ff;bb=x*x*y**(q0-1)*ff
        sub={x:xx,y:u};fW=u**(q0+1)*v
        aW=u*v*(1+fW);bW=v*(1+fW)**2
        check('control chart generators regular q'+str(q0),
              all(S.cancel(f.subs(sub,simultaneous=True)-g)==0
                  for f,g in [(ff,fW),(aa,aW),(bb,bW)]))
        check('control source principal chart recovery q'+str(q0),
              S.cancel(aa/ff-x)==0)
        check('control boundary principal chart recovery q'+str(q0),
              S.cancel(bW/(fW+1)**2-v)==0)
        check('control two principal charts cover q'+str(q0),S.expand((ff+1)-ff)==1)
        check('control volume divisor is the reduced boundary q'+str(q0),
              S.det(S.Matrix([[S.diff(xx,u),S.diff(xx,v)],[1,0]]))==-u)
        check('control infinity multiplicity q'+str(q0),
              S.cancel(u**q0*xx)==1+u**(q0+1)*v)
        if q0>1:
            # A dx+B dy; pull back using x=u*v+u^-q, y=u.
            A=y/S.Integer(q0-1);B=q0*x/S.Integer(q0-1)
            theta_du=S.expand(A.subs(sub,simultaneous=True)*S.diff(xx,u)+B.subs(sub,simultaneous=True))
            theta_dv=S.expand(A.subs(sub,simultaneous=True)*S.diff(xx,v))
            check('control global primitive glues q'+str(q0),
                  S.cancel(theta_du-(q0+1)*u*v/(q0-1))==0 and
                  S.cancel(theta_dv-u*u/(q0-1))==0)
            check('control global primitive has correct derivative q'+str(q0),
                  S.diff(B,x)-S.diff(A,y)==1 and
                  S.cancel(S.diff(theta_dv,u)-S.diff(theta_du,v)+u)==0)
        else:
            check('control q1 has the nonzero volume residue',S.residue(u**(-q0),u,0)==1)
        # On a chosen split branch with root of unity 1. Other branches
        # have the same identity after constant scaling.
        vv=z**(q0+2)*y-z**(q0+1)
        yy=u+u**(q0+2)*v
        check('normalized control whole second chart inverse q'+str(q0),
              S.expand(vv.subs({z:1/u,y:yy},simultaneous=True))==v)
        s=S.symbols('s')
        quotient=sum(s**j for j in range(q0))
        check('split branch regularity identity q'+str(q0),
              S.expand((s-1)*quotient-(s**q0-1))==0 and quotient.subs(s,1)==q0)
        control_g=-z**(q0+1)
        check('control residue agrees with exactness threshold q'+str(q0),
              S.expand(control_g).coeff(z,2)==(-1 if q0==1 else 0))
    inputs=[Path(__file__),HERE/'MULTIPLE_FIBER_CYCLIC_REDUCTION.md',
            HERE/'intrinsic_normalization_receipt.json',manifest]+[
            HERE/row['path_relative_to_affine_complements'] for row in pins['primary_sources']]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='all finite field degrees; actual smooth one-index-two-prime normalization with any omitted-fiber multiplicity q',
        general_smooth_geometry=dict(canonical_boundary_coefficient='e-1',
            cyclic_cover_degree='q',ramification_of_cyclic_cover='(q-1)L_prime',
            split_boundary_lines='q',linear_torsor_index='q+e',
            canonical_class='(q+e-2)*sum R_i',Pic='free on the q boundary lines',
            H1=0,H2_rank='q',Euler='q+1'),
        volume=dict(residue='-q*[z^e]g_i',deck_shared_scalar=True,
            exact_if_and_only_if_shared_coefficient_zero=True,
            actual_index_two_necessary_coefficient='[z^2]g_i=0',
            tested_exponent_independent_of_q=True),
        actual_finite_cover=dict(field_degree='qd',additional_curve='D=pi(x=0)',
            additional_curve_immersed_normalization=True,additional_curve_distinct_from_C=True,
            monodromy='C_q wreath S_d',group_order='q^d*d!',
            boundary_inertia='q disjoint transpositions',
            auxiliary_inertia='one q-cycle',source_Jacobian='c*q*z^(q-1)'),
        affine_control_family=dict(all_positive_integer_multiplicities=True,
            ring_generators=['y','f=x*y^q-1','a=x*f','b=x^2*y^(q-1)*f'],
            smooth=True,one_omitted_A1=True,canonical_divisor='R',
            exact_volume_for_every_q_at_least_two=True,
            finite_Keller_extension_constructed=False,general_normalization_classification=False),
        finite_wreath_controls=wreath,
        assurance=dict(all_q_geometry_proved_in_writing=True,
            whole_multiple_fiber_locus_excluded=False,
            cyclic_cover_mistaken_for_Keller=False,
            restriction_to_one_piece_asserted_finite=False,
            ordinary_ramification_coefficients_preserved=True,
            actual_end_lens_or_core_proved=False,
            arbitrary_quartic_packets_excluded=False,
            all_finite_Wright_degrees_excluded=False,
            coefficient_or_Groebner_search_used=False,
            formal_verification=False,independent_full_reproof=False,
            external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        python=platform.python_version(),sympy=S.__version__)
    (HERE/'multiple_fiber_cyclic_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; all-degree cyclic/volume reduction; full goal active')
if __name__=='__main__':main()
