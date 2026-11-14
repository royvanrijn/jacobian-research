"""Complete conductor models and central inverse pole/index controls."""
from pathlib import Path
import hashlib,json,math,platform,time
import sympy as s
import flint
import central_conductor as C

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression(poly,variable):
    return s.Add(*(s.Rational(str(co))*variable**i for i,co in enumerate(poly.coeffs())))


def order(poly,variable):
    return min(mon[0] for mon,co in s.Poly(poly,variable).terms())


def main():
    started=time.monotonic()
    models=C.build_models(check)
    t,v,T,a=s.symbols('t v T a')
    for family,value in [('L',0),('L',1),('K',1)]:
        model=next(row for row in models if row['family']==family)
        p=expression(C.specialize(C.decode(model['normalization_p'],C.PT),value),t)
        q=expression(C.specialize(C.decode(model['normalization_q'],C.PT),value),t)
        conductor=expression(C.specialize(C.decode(model['conductor'],C.PT),value),t)
        dp=s.div(p.subs(t,v)-p,v-t,v)[0]
        dq=s.div(q.subs(t,v)-q,v-t,v)[0]
        check(f'{family}{value} independent divided-difference resultant',s.expand(s.resultant(dp,dq,v)-conductor)==0)
    model=next(row for row in models if row['family']=='K')
    for value,cusp in [(s.Rational(-2,3),1),(s.Rational(4,15),-1)]:
        conductor=expression(C.specialize(C.decode(model['conductor'],C.PT),str(value)),t)
        p=expression(C.specialize(C.decode(model['normalization_p'],C.PT),str(value)),t)
        q=expression(C.specialize(C.decode(model['normalization_q'],C.PT),str(value)),t)
        check('scaled K cusp'+str(value)+' derivative parameter',s.diff(p,t).subs(t,cusp)==0 and s.diff(q,t).subs(t,cusp)==0)
        check('scaled K cusp'+str(value)+' conductor order2',conductor.subs(t,cusp)==0 and s.diff(conductor,t).subs(t,cusp)==0 and s.diff(conductor,t,2).subs(t,cusp)!=0)
    prior=json.loads((HERE/'generic_line_ledger.json').read_text())
    weights=dict(nodes=2,triples=6,tacnodes=4,cusps=2,A5=6)
    for i,row in enumerate(prior['first_image_catalogue']):
        check(f'conductor catalogue{i} total24',sum(weights[k]*value for k,value in row['singularity_counts'].items())==24)
    pR,pS,pD,b,u,j,l,W=s.symbols('p_R p_S p_D b u j l W',positive=True)
    rows=[]
    for grouping,degree,conductor_exponent,B_exponent,H_exponent,degree_formula in [
        ('A_distinct',4,'p_R','p_S','p_D',24*pR+b*pS+u*pD),
        ('A_H=C',3,'p_R+p_D','p_S','0',24*(pR+pD)+b*pS),
        ('A_H=B',4,'p_R','p_S+p_D','0',24*pR+b*(pS+pD)),
        ('B',4,'p_R','p_S','0',24*pR+b*pS),
        ('N',4,'p_R','p_S','0',24*pR+b*pS),
        ('C',4,'p_R','p_S','0',24*pR+b*pS)]:
        rows.append(dict(grouping=grouping,central_T_degree=degree,
            central_leading_coefficient_mod_H_C=dict(conductor_adjoint_exponent=conductor_exponent,
                H_B_exponent=B_exponent,H_H_exponent=H_exponent,nonzero_constant_factor_unknown=True),
            exact_pullback_degree=str(degree_formula),actual_poles=None,actual_inverse_coefficient=None))
    check('general F2 central-degree budget',s.expand((24*pR+30*j*pS+30*l*pD)-(30*(2*pR+j*pS+l*pD)-36*pR))==0)
    check('shared C F2 central-degree budget',s.expand(24*(pR+pD)+30*j*pS-(30*(2*pR+2*pD+j*pS)-36*(pR+pD)))==0)
    check('F2 caps684 and648',30*24-36==684 and 30*24-36*2==648)
    # Local models have one fold component w^2=t and one unramified
    # omitted component, with T=(1+w)/w^p_R and T=1/(3t^p_D).
    controls=[]
    for pr in range(1,5):
        for pd in range(1,4):
            quadratic=s.resultant(v*v-t,v**pr*T-1-v,v)
            linear=1-3*t**pd*T
            inverse=s.Poly(s.expand(linear*quadratic),T)
            disc=s.discriminant(inverse.as_expr(),T)
            valuations=sorted([s.Rational(pr,2),s.Rational(pr,2),s.Integer(pd)])
            check(f'node poles{pr},{pd} constant unit and leading charge',inverse.nth(0).subs(t,0)!=0 and order(inverse.LC(),t)==pr+pd)
            check(f'node poles{pr},{pd} reciprocal coefficient bounds',all(order(inverse.nth(i),t)>=s.ceiling(sum(valuations[:i])) for i in range(1,4)))
            bound=s.Integer(pr)+4*min(s.Rational(pr,2),s.Integer(pd))
            disc_order=order(disc,t)
            index=(disc_order-1)//2
            expected=int(s.ceiling((bound-1)/2))
            check(f'node poles{pr},{pd} exact local index bound control',disc_order%2==1 and index==expected)
            check(f'node poles{pr},{pd} discriminant product identity',
                  s.expand(disc-s.discriminant(quadratic,T)*s.resultant(linear,quadratic,T)**2)==0)
            controls.append(dict(p_R=pr,p_D=pd,leading_pole_charge=pr+pd,raw_discriminant_order=disc_order,
                normalized_different=1,normalization_index_order=index,bound_attained_in_local_model=True,
                original_ramified_prime_count=1,unramified_omitted_prime_count=1,global_Keller_realization=False))
    check('unit pole node forces extra index',controls[0]['raw_discriminant_order']==3 and controls[0]['normalization_index_order']==1)
    primitive=s.resultant(v*v-t,v*T-1-v,v)
    raw=s.factor(s.discriminant((1-3*t*T)*primitive,T))
    check('explicit coordinate-collision factor',s.expand(raw-4*t**3*(1-15*t+9*t*t)**2)==0)
    check('collision factor does not add normalized ramification',(1-15*t+9*t*t).subs(t,0)==1)
    # Direct rational-root factor control, with all denominator roots kept.
    denom=t*t-1
    numerator=t**3+2
    quadratic=(t-2)*T*T+T+1
    inverse=s.expand((denom*T-numerator)*quadratic)
    check('split section has reduced denominator',s.gcd(denom,numerator)==1)
    check('split section exact inverse factor',s.cancel(inverse.subs(T,numerator/denom))==0)
    check('split resultant-square discriminant',s.expand(s.discriminant(inverse,T)-s.discriminant(quadratic,T)*s.resultant(denom*T-numerator,quadratic,T)**2)==0)
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT20',conductor_models=models,
        conductor_degree=24,target_conductor_representative_degree_at_most=7,
        conductor_is_full_affine_normalization_conductor=True,
        generic_admissible_T_separates_retained_C_embeddings=True,
        central_inverse_coefficient_ideal_is_unit=True,central_root_cover_normalization='E_C of PF2OMIT19',
        exact_leading_coefficient_rows=rows,
        conditional_F2_central_pullback_degree_caps=dict(general=684,A_H_C=648,requires_source_degrees=[75,125],geometric_sharpness_claimed=False),
        clean_central_normalized_discriminant='c_C*R_CB^(e_S-1)',
        clean_central_raw_discriminant='constant*c_C*R_CB^(e_S-1)*J_C^2',
        restricting_global_degree6_discriminant_recovers_central_discriminant=False,
        shared_H_C_node_raw_discriminant_bound='p_R+4min(p_R/2,p_D)',
        shared_H_C_node_index_bound='ceil((p_R+4min(p_R/2,p_D)-1)/2)',
        local_pole_index_controls=controls,
        clean_C_split_central_factor='T-A(t), A polynomial',
        clean_A_H_C_split_central_factor='c_C^p_D*T-A(t), gcd(A,c_C)=1',
        actual_R6=None,actual_central_inverse=None,actual_first_image_parameter=None,actual_p_R=None,actual_p_S=None,actual_p_D=None,
        actual_B_equation=None,actual_H_equation=None,actual_central_index=None,actual_section_numerator=None,
        actual_global_source_realization=None,unrestricted_singleton_zero_properness_proved=False,full_goal_complete=False)
    ledger_path=HERE/'central_inverse_ledger.json'
    ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=['jacobian-followup/jc2-omitted-divisors/'+name for name in ['CENTRAL_INVERSE_CONDUCTOR.md',
        'central_conductor.py','CENTRAL_IMAGE_FIBRE.md','central_image_ledger.json','ADMISSIBLE_INVERSE.md',
        'DISCRIMINANT_LEDGER.md','GENERIC_LINE_ACCOUNTING.md','generic_line_ledger.json']]
    inputs+=['jacobian-followup/jc2-r7-image-family/family.py','jacobian-followup/jc2-r7-elliptic-contact/PROOF.md',
        'jacobian-followup/jc2-r7-elliptic-contact/models.json','archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={path:digest(ROOT/path) for path in inputs},ledger_sha256=digest(ledger_path),
        runtime_seconds=round(time.monotonic()-started,6),python=platform.python_version(),sympy=s.__version__,flint=flint.__version__,
        proof_boundary='Complete chart conductor/gradient/adjoint identities and independent divided-resultant/cusp controls exact. Central inverse normalization, conductor interpretation, pole-divisor identities and trace-index application are written. Exact source-degree budget algebra and twelve local pole/index models support, but do not realize, global covers. Actual inverse, poles, other images and source remain absent.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'central_inverse_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({key:receipt[key] for key in ['outcome','check_count','runtime_seconds']}))


if __name__=='__main__':
    main()
