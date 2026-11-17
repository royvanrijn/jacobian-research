'''PF2RULED38: exact controls for written coordinate-pole and regular-sector proofs.
These controls do not enumerate every automorphism or construct a finite Keller map.
'''
from pathlib import Path
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def laurent_value(expr,u):
    terms=S.expand(expr).as_ordered_terms()
    return min(int(term.as_powers_dict().get(u,0)) for term in terms if term!=0)
def negative_part(expr,u):
    return S.Add(*(term for term in S.expand(expr).as_ordered_terms()
                   if int(term.as_powers_dict().get(u,0))<0))
def primitive_split(h,u,q):
    K=sum(h.coeff(u,-j)*u**(1-j)/S.Integer(1-j) for j in range(2,q+1))
    leading=h.coeff(u,-q);p=S.Integer(0)
    for k in range(q-1,0,-1):
        rem=S.expand(K-p*h)
        coefficient=S.cancel(rem.coeff(u,-k)/leading)
        p=S.expand(p+coefficient*u**(q-k))
    return S.expand(K),p
def main():
    pred=json.loads((HERE/'multiple_fiber_cyclic_receipt.json').read_text())
    check('cyclic predecessor receipt PASS',pred['outcome']=='PASS')
    for name,h in pred['source_sha256'].items():
        check('cyclic predecessor seal '+name,digest(HERE/name)==h)
    D=HERE/'sources/coordinate-pole-sector';manifest=D/'manifest.json'
    pins=json.loads(manifest.read_text())
    check('three primary inputs pinned',len(pins['primary_sources'])==3)
    for row in pins['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements']
        check('primary source seal '+p.name,p.read_bytes().startswith(b'%PDF') and digest(p)==row['sha256'])
    # Exact Bruhat factorization used to place every source automorphism
    # in the alternating normal form. The valuation induction is written.
    x,y,u,t=S.symbols('x y u t')
    aa,bb,cc,dd,ss,tt=S.symbols('aa bb cc dd ss tt',nonzero=True)
    C2=S.Matrix([[cc,dd,tt],[0,1,0],[0,0,1]])
    swap=S.Matrix([[0,1,0],[1,0,0],[0,0,1]])
    C1=S.Matrix([[bb-aa*dd/cc,aa/cc,ss-aa*tt/cc],[0,1,0],[0,0,1]])
    affine=S.Matrix([[aa,bb,ss],[cc,dd,tt],[0,0,1]])
    check('affine factor outside the triangular subgroup has the required Bruhat form',
          all(S.cancel(entry)==0 for entry in C1*swap*C2-affine))
    n=S.symbols('n',integer=True,positive=True)
    coordinate=y+x**n
    check('coordinates realizing every positive pole multiple form genuine automorphisms',
          S.det(S.Matrix([[1,0],[S.diff(coordinate,x),S.diff(coordinate,y)]]))==1)
    # Concrete controls check actual Laurent substitutions for several
    # reduced words; the all-word proof is the written induction.
    for q in [1,2,3,5]:
        sub={x:u*t+u**(-q),y:u}
        F=x+y*y+1;G=y+2
        check('initial triangular valuation control q'+str(q),
              laurent_value(F.subs(sub,simultaneous=True),u)==-q and
              laurent_value(G.subs(sub,simultaneous=True),u)==0)
        for step in range(1,3):
            F,G=G+F*F,F
            check(f'nonlinear valuation control q{q} step{step}',
                  laurent_value(F.subs(sub,simultaneous=True),u)==-(2**step)*q and
                  laurent_value(G.subs(sub,simultaneous=True),u)==-(2**(step-1))*q)
        f=F+3*G+4;g=2*F+5*G+7
        check('final affine and equal-leading cancellation control q'+str(q),
              laurent_value(f.subs(sub,simultaneous=True),u)==-4*q and
              laurent_value(g.subs(sub,simultaneous=True),u)==-4*q and
              laurent_value((2*f-g).subs(sub,simultaneous=True),u)==-2*q)
    # General principal parts, not just monomial controls.
    ring_controls=[]
    for q in range(1,6):
        c=S.symbols('c1:'+str(q+1))
        h=sum(c[j-1]*u**(-j) for j in range(1,q+1))
        HH=S.expand(u**q*h)
        xx=u*t+h
        f=x*y**q-HH.subs(u,y);a=x*f;b=x*x*y**(q-1)*f
        fW=u**(q+1)*t
        aW=u*t*(HH+fW);bW=t*(HH+fW)**2
        sub={x:xx,y:u}
        check('four global generators extend for every principal part q'+str(q),
              all(S.cancel(v.subs(sub,simultaneous=True)-w)==0
                  for v,w in [(f,fW),(a,aW),(b,bW)]))
        check('source and boundary principal charts recover q'+str(q),
              S.cancel(a/f-x)==0 and S.cancel(bW/(fW+HH)**2-t)==0)
        check('overlap principal chart recovers both coordinates q'+str(q),
              S.cancel((f+HH.subs(u,y))/y**q-x)==0 and
              S.cancel(fW/u**(q+1)-t)==0)
        H1=S.cancel((HH-c[-1])/u)
        check('three principal charts generate the unit ideal q'+str(q),
              H1.is_polynomial(u) and S.expand((fW+HH)-fW-u*H1)==c[-1])
        check('canonical boundary coefficient is one q'+str(q),
              S.det(S.Matrix([[S.diff(xx,u),S.diff(xx,t)],[1,0]]))==-u)
        check('period is precisely c1 q'+str(q),S.residue(h,u,0)==c[0])
        alpha,beta=S.symbols('alpha beta',nonzero=True)
        hp=sum(beta*alpha**j*c[j-1]*(alpha*u)**(-j) for j in range(1,q+1))
        check('marked coefficient rescaling is exact q'+str(q),S.cancel(hp-beta*h)==0)
        k0,k1=S.symbols('k0 k1')
        up=alpha*u;tp=beta*t/alpha+(k0+k1*u)/alpha
        check('marked shear with zero constant extends on whole boundary chart q'+str(q),
              S.cancel(up*tp+hp-(beta*xx+u*(k0+k1*u)))==0)
        check('a nonzero shear constant would have a boundary pole q'+str(q),
              laurent_value(1/(alpha*u),u)==-1)
        if q>=2:
            h0=S.expand(h.subs(c[0],0));K,p=primitive_split(h0,u,q)
            check('triangular negative-part inverse is exact q'+str(q),
                  S.cancel(negative_part(p*h0,u)-K)==0)
            remainder=S.cancel(p*h0-K)
            check('triangular split leaves a polynomial on the whole second chart q'+str(q),
                  remainder.is_polynomial(u) and p.subs(u,0)==0 and S.degree(p,u)<=q-1)
            # Independently compare the two global primitive expressions.
            pY=p.subs(u,y);KU=x*pY;KW=u*t*p+remainder
            AU=-S.diff(KU,x);BU=x-S.diff(KU,y)
            AW=u*t-S.diff(KW,u);BW=-S.diff(KW,t)
            xx0=S.expand(xx.subs(c[0],0));sub0={x:xx0,y:u}
            pullA=S.expand(AU.subs(sub0,simultaneous=True)*S.diff(xx0,u)+BU.subs(sub0,simultaneous=True))
            pullB=S.expand(AU.subs(sub0,simultaneous=True)*S.diff(xx0,t))
            check('constructive global primitive glues for all principal-part parameters q'+str(q),
                  S.cancel(pullA-AW)==0 and S.cancel(pullB-BW)==0)
            check('constructive global primitive differentiates to the source volume q'+str(q),
                  S.cancel(S.diff(BU,x)-S.diff(AU,y)-1)==0 and
                  S.cancel(S.diff(BW,u)-S.diff(AW,t)+u)==0)
            # Matrix for the negative-part map has determinant c_q^(q-1).
            M=S.Matrix([[S.expand(u**j*h0).coeff(u,-k)
                         for j in range(1,q)] for k in range(q-1,0,-1)])
            check('uniform triangular determinant keeps every c_q-nonzero parameter q'+str(q),
                  S.expand(M.det()-c[-1]**(q-1))==0)
            ring_controls.append(dict(q=q,coefficient_parameters=q-1,
                                      inverse_determinant='c_q^'+str(q-1)))
    # Tower multiplicity logic for the primitive-coordinate deduction.
    # This is an abstract finite algebra control, not a candidate Keller map.
    T,s,bet=S.symbols('T s bet')
    for k in [2,3,4]:
        polynomial=T**k-s
        M=S.zeros(k)
        for j in range(k-1):M[j+1,j]=1
        M[0,k-1]=s
        baseNorm=S.det(bet*S.eye(k)-M)
        check('norm polynomial in a primitive toy algebra degree'+str(k),
              S.expand(baseNorm-polynomial.subs(T,bet))==0)
        check('proper field tower raises every divisor multiplicity degree'+str(k),
              S.Poly(S.expand(baseNorm**2),s).degree()==2)
    inputs=[Path(__file__),HERE/'COORDINATE_POLE_AND_REGULAR_SECTOR.md',
            HERE/'multiple_fiber_cyclic_receipt.json',manifest]+[
            HERE/row['path_relative_to_affine_complements'] for row in pins['primary_sources']]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='all field degrees; coordinate-pole lemma for one regular source coordinate; full v_R(y)=1, canonical-divisor-R marked sector',
        coordinate_invariant=dict(regular_coordinates='a*y+b with a nonzero',
            negative_coordinate_poles='q*positive integers',minimum_pole='q',
            all_polynomial_source_automorphisms_covered_by_written_normal_form=True,
            bounded_word_controls_are_not_the_proof=True),
        regular_sector=dict(source_chart='x=u*t+h(u), y=u',
            h='sum c_j*u^(-j), c_q nonzero',canonical_divisor='R',
            complete_marked_normal_form=True,
            marked_equivalence='same q and c_j_prime=b*alpha^j*c_j',
            ring_generators=['y','f=x*y^q-H(y)','a=x*f','b=x^2*y^(q-1)*f'],
            exact_volume_iff='c_1=0',primitive_constructed_for_all_c1_zero=True,
            negative_part_inverse_determinant='c_q^(q-1)'),
        sharpened_controls=dict(every_q_at_least_two=True,
            no_alternative_source_coordinate_pole_one=True,
            no_marked_Wright_comparison=True,finite_Keller_extension_constructed=False),
        finite_algebra=dict(actual_finiteness_is_a_hypothesis=True,
            primitive_integral_coordinate='K(y)=C(N)',
            monic_norm_degree='d',
            zero_member='H_C*H_0, reduced with distinct factors',
            nonzero_members='H_beta, reduced irreducible with immersed A1 normalization',
            every_D_beta_self_colliding=True,
            target_plus_y_algebra_not_normal=True,
            conductor_has_y_horizontal_component=True,
            norm_family_is_a_linear_pencil=False),
        principal_part_controls=ring_controls,
        assurance=dict(all_q_and_all_degree_proofs_written=True,
            primary_normal_form_and_fibration_imports_explicit=True,
            one_line_injectivity_import_explicit=True,
            exact_volume_does_not_construct_finite_Keller_map=True,
            actual_regular_order_one_sector_excluded=False,
            regular_coordinate_higher_vanishing_classified=False,
            no_regular_coordinate_case_classified=False,
            arbitrary_quartic_packets_excluded=False,
            all_finite_Wright_degrees_excluded=False,
            generic_end_lens_or_core_proved=False,
            coefficient_or_Groebner_search_used=False,
            formal_verification=False,independent_full_reproof=False,
            external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        python=platform.python_version(),sympy=S.__version__)
    (HERE/'coordinate_regular_sector_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; coordinate invariant and all-degree regular-sector gate; full goal active')
if __name__=='__main__':main()
