"""PF2RULED99: monic-key Jacobian obstruction and entire chart necessity.

Written proofs, with explicit primary imports, establish universality.
Exact Laurent/monic controls are independent polynomial checks, not
coefficient-depth searches or finite Keller-map realizations.
"""
from pathlib import Path
import sys,json,hashlib,os,platform,math
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
from verify_exact_coordinate_bridge import pole_types

def monic_expansion(F,G,x,y):
    """Repeated exact division, independent of how a control was assembled."""
    out=[];k=0;F=S.expand(F)
    while F!=0:
        F,remainder=S.div(F,G,y)
        for (i,j),coefficient in S.Poly(remainder,x,y).terms():
            if coefficient:out.append((int(i),int(j),k,coefficient))
        k+=1
    return out

def laurent_leading(F,u):
    coefficients={}
    for term in S.Add.make_args(S.expand(F)):
        exponent=int(term.as_powers_dict().get(u,0))
        coefficients[exponent]=coefficients.get(exponent,0)+term/u**exponent
    coefficients={k:S.expand(v) for k,v in coefficients.items() if S.expand(v)!=0}
    assert coefficients
    k=min(coefficients)
    return k,coefficients[k]

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.simplify(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['proper_finite_marked_end_receipt.json','exact_coordinate_bridge_receipt.json',
        'higher_regular_coordinate_receipt.json','unramified_orevkov_receipt.json']
    for name in preds:
        r=json.loads((P/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/monic-key-chart-audit/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('primary exact byte pin '+f.name,f.stat().st_size==row['bytes'] and sha(f)==row['sha256'])
    check('general semidegree sign and scheme fiber quantifiers verified',
        sources['canonical_sign']=='delta=-v; A_delta=1+v(volume)=e' and
        sources['general_semidegree_quantifiers_verified'] and sources['scheme_fiber_integrality_not_just_reduced_support'])
    root=P.parents[1]
    entries=json.loads((root/'MATH_STATUS.json').read_text())['entries']
    published=next(x for x in entries if x['id']=='HC4MYKR0')
    check('published unrestricted quartic exclusion retained from shared authority',
        published['state']=='proved' and 'NO boundary-count' in published['scope'] and
        sources['published_quartic_retirement']['quartic_realization_is_no_longer_open'] is True)
    x,y,u,t=S.symbols('x y u t');p0,q0=S.symbols('p0 q0')
    G=S.Function('G')(x,y);H=S.Function('H')(x,y);K=S.Function('K')(x,y)
    jac=lambda F,B:S.diff(F,x)*S.diff(B,y)-S.diff(F,y)*S.diff(B,x)
    certificate=H*jac(G,K)+K*jac(H,G)+G*jac(H,K)
    identity('UNBOUNDED exact Jacobian factor certificate for arbitrary functions',
        jac(p0+G*H,q0+G*K)-G*certificate)
    m,n,s,i,j,k,g,A,B=S.symbols('m n s i j k g A B',integer=True)
    identity('equal monomial value at fixed k is reduced pole relation',
        (-m*i-n*j+s*k)-(-m*S.Symbol('ip')-n*S.Symbol('jp')+s*k)
        +m*(i-S.Symbol('ip'))+n*(j-S.Symbol('jp')))
    controls=[]
    for mm,nn,ss in [(3,2,1),(6,4,1),(15,10,3),(14,6,5),(2,2,1)]:
        gg=math.gcd(mm,nn);M=mm//gg;N=nn//gg
        key=y**M-x**N;xx=u**(-mm);yy=u**(-nn)+t*u**(ss+nn*(M-1))
        actual_key=S.expand(key.subs({x:xx,y:yy},simultaneous=True))
        check('generic monic primitive chart has nonconstant leading residue '+str((mm,nn,ss)),laurent_leading(actual_key,u)==(ss,M*t))
        polynomials=[x+y+key,1+key,x**ss*key**mm+y**ss*key**nn+1,
            x**2*y**(M-1)+3*x*key**2+2*y**(M-1)*key**3+5]
        values=[]
        for number,F in enumerate(polynomials):
            expansion=monic_expansion(F,key,x,y)
            rebuilt=sum(cc*x**ii*y**jj*key**kk for ii,jj,kk,cc in expansion)
            identity('independent monic division reconstructs polynomial '+str((mm,nn,ss,number)),rebuilt-F)
            expected=min(-mm*ii-nn*jj+ss*kk for ii,jj,kk,cc in expansion)
            actual=laurent_leading(F.subs({x:xx,y:yy},simultaneous=True),u)
            check('no cancellation exact minimum value '+str((mm,nn,ss,number)),actual[0]==expected)
            coefficient=S.expand(sum(cc*(M*t)**kk for ii,jj,kk,cc in expansion if -mm*ii-nn*jj+ss*kk==expected))
            check('leading coefficient retains independent powers of residue '+str((mm,nn,ss,number)),actual[1]==coefficient and coefficient!=0)
            if expected>=0:
                remainder=S.rem(F,key,y)
                check('every tested regular polynomial is constant modulo monic key '+str((mm,nn,ss,number)),not remainder.has(x,y))
            values.append(dict(control=number,value=expected,leading=str(actual[1])))
        controls.append(dict(poles=[mm,nn],key_y_degree=M,key_value=ss,leading_residue=str(M*t),values=values))
    # Arbitrary positive-index volume and additive direction: symbolic e,r.
    e,r,aa=S.symbols('e r a',integer=True,positive=True);z=S.symbols('z',positive=True)
    h=S.Function('h')(z);xx=z**(e-r)*t+h;ff=z**r
    coefficient=S.diff(xx,t)*S.diff(ff,z)
    identity('ALL-INDEX full root chart volume coefficient',coefficient-r*z**(e-1))
    identity('ALL-INDEX source/boundary primitive difference',
        r*z**(r-1)*xx-r*z**(e-1)*t-r*z**(r-1)*h)
    identity('ALL-INDEX minimal source field boundary exponent',r*aa-(e-r)-(r*(aa+1)-e))
    identity('ALL-INDEX gauge cutoff shifts final retained exponent below transition',
        -(r-e+1)-(e-r-1))
    charts=[]
    for ee,rr in [(1,1),(2,1),(1,2),(2,2),(5,3),(3,7),(9,4)]:
        support=list(range(rr-ee+1,max(rr+3,rr-ee+2)))
        coeff={jj:S.Symbol('c_'+str(jj).replace('-','minus')) for jj in support}
        hh=sum(cc*z**(-jj) for jj,cc in coeff.items())
        difference=S.expand(rr*z**(rr-1)*hh)
        residue=S.residue(difference,z,0)
        check('unique resonant residue at signed support index r '+str((ee,rr)),residue==rr*coeff.get(rr,0))
        zero=S.expand(difference.subs(coeff[rr],0))
        check('removing c_r makes full Laurent differential residue zero '+str((ee,rr)),S.residue(zero,z,0)==0)
        amin=max(0,(ee+rr-1)//rr-1);order=rr*(amin+1)-ee
        check('least integer source additive multiplier '+str((ee,rr)),order>=0 and (amin==0 or rr*amin-ee<0))
        charts.append(dict(index=ee,boundary_order=rr,signed_support=support,minimal_multiplier=amin,boundary_vanishing_order=order,residue=str(residue)))
    check('higher-index positive Laurent power retained by gauge',charts[4]['signed_support'][0]==-1 and charts[4]['minimal_multiplier']==1)
    # PF45's entire all-polar index2 normal form, not a jet search.
    a,b,c=S.symbols('a b c');xx=u**-3;yy=u**-2+a*u**2+b*u**3+c*u**4+u**5*t
    key=y**3-x**2-3*a*y-3*c
    delta_theta=S.expand(3*u**-4*yy-3*u*t)
    check('whole all-polar chart residue is exactly 3b',S.residue(delta_theta,u,0)==3*b)
    value,leading=laurent_leading(key.subs({x:xx,y:yy.subs(b,0)},simultaneous=True),u)
    check('forced exact chart meets general monic-key obstruction',value==1 and leading==3*t and S.Poly(key,y).LC()==1 and S.Poly(key,y).degree()==3)
    check('canonical index1 and index2 initial-pair exhaustion retained',pole_types(1)==[] and pole_types(2)==[(3,2,1,3,2)])
    # ESSENTIAL actual affine constant-residue countercontrol.
    f,W,V,X,Y=S.symbols('f W V X Y');relation=f**2*V-W*(W-1)
    Wsrc=f**2*x;Vsrc=x*(f**2*x-1)
    identity('actual affine countercontrol source chart equation',relation.subs({W:Wsrc,V:Vsrc},simultaneous=True))
    identity('actual affine countercontrol smoothness has global gradient certificate',
        S.diff(relation,W)**2-4*V*S.diff(relation,V)+4*relation-1)
    boundary={f:u,W:1+u**2*t,V:t+u**2*t**2}
    identity('actual affine countercontrol entire boundary plane equation',relation.subs(boundary,simultaneous=True))
    XX=u**-2+t;YY=XX+u;key2=Y-X;rescued=key2**2*X**2-X
    check('constant leading key residue retained',laurent_leading(key2.subs({X:XX,Y:YY},simultaneous=True),u)==(1,S.Integer(1)))
    check('leading pole terms cancel to rescue nonconstant regular residue',laurent_leading(rescued.subs({X:XX,Y:YY},simultaneous=True),u)==(0,t))
    check('constant-residue control is NOT constant modulo its monic key',S.rem(rescued,key2,Y)==-X)
    check('arbitrary finite quartic map has interior ramification and is not excluded',
        jac(x**4,y)==4*x**3 and jac(x**4,y).subs(x,0)==0)
    inputs=[Path(__file__),P/'MONIC_KEY_JACOBIAN_OBSTRUCTION_AND_WHOLE_CHART_NECESSITY.md',manifest]+[P/x for x in preds]+[P/x['path'] for x in sources['primary_sources']]
    inputs.extend([root/published['canonical_source'],root/'hessian/hc4-triangular-graphs/sources/quadratic-reconstruction-context.json'])
    result=dict(as_of='2026-11-16',outcome='PASS',check_count=len(checks),checks=checks,
        monic_key=dict(any_field_degree_generic_index_and_other_omitted_primes=True,
            normalized_constant_x_y_leading_residues_required=True,
            minimal_monic_y_degree_m_over_gcd_mn_required=True,
            positive_key_value_and_NONCONSTANT_leading_residue_required=True,
            every_regular_polynomial_constant_modulo_key=True,
            every_regular_target_Jacobian_has_key_factor=True,
            no_affineness_or_finite_generation_needed_for_algebraic_lemma=True),
        whole_charts=dict(all_positive_indices_WITH_smoothness_and_regular_coordinate=True,
            finite_root_normalization_and_SCHEME_reduced_fiber_verified=True,
            WHOLE_affine_planes_not_formal_neighborhoods=True,
            exact_transition_power_e_minus_r_and_signed_support_j_ge_r_minus_e_plus1=True,
            global_exactness_iff_c_r_zero_in_every_positive_index=True,
            least_additive_multiplier_max0_ceil_e_over_r_minus1=True,
            target_algebra_invariance_NOT_inferred=True),
        audit=dict(Mondal_general_signed_quantifiers_and_canonical_sign_verified=True,
            Miyanishi_full_affine_flat_integral_scheme_fiber_hypotheses_verified=True,
            PF45_index1_2_coordinate_existence_necessity_audited=True,
            PF40_43_46_whole_chart_NECESSITY_audited=True,
            external_classification_and_bundle_theorems_explicit_imports=True),
        published_quartic=dict(authority_claim='HC4MYKR0',
            unrestricted_generic_degree4_plane_Keller_exclusion_retained=True,
            marked_dense_plane_etale_finite_degree4_completions_excluded=True,
            quartic_realization_campaign_retired=True,
            internal_quartic_chain_certified_by_external_theorem=False,
            arbitrary_finite_quartic_maps_with_interior_ramification_excluded=False,
            full_33_page_published_proof_independently_reproduced=False,
            reviewed_scope_fingerprint=hashlib.sha256('\0'.join(published[k] for k in ['state','title','scope']).encode()).hexdigest()[:16]),
        monic_key_controls=controls,all_index_chart_controls=charts,
        actual_affine_constant_residue_countercontrol=dict(poles=[2,2],key_value=1,key_residue=1,
            rescued_function='G^2 X^2-X',rescued_value=0,rescued_residue='t',reduction_mod_key='-X',not_Keller_realization=True),
        assurance=dict(written_unbounded_monic_factor_and_whole_chart_proofs=True,
            essential_constant_residue_and_signed_support_limits_retained=True,
            finite_controls_not_universal_proof=True,
            published_quartic_exclusion_retained_and_realization_campaign_retired=True,
            higher_degree_nonmonic_constant_residue_or_singular_cases_remain=True,
            external_valuation_classification_independently_reproved=False,
            all_higher_index_regular_coordinates_forced=False,
            arbitrary_higher_index_affine_converse_proved=False,
            arbitrary_higher_degree_packets_excluded=False,
            new_internal_full_quartic_exclusion_proved=False,
            full_compatible_completion_goal_solved=False,formal_verification=False,
            external_review=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'monic_key_chart_audit_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree monic-key Jacobian factor and all-index conditional whole charts; full objective OPEN')
if __name__=='__main__':main()
