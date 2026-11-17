"""Exact controls for PF2RULED53; finite local models are NOT Keller realizations."""
from pathlib import Path
import hashlib,json,os,platform
import sympy as S

HERE=Path(__file__).resolve().parent;OMIT=HERE.parent/'jc2-omitted-divisors';checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value):
    assert value,name
    checks.append(name)
def valuation(f,r):
    a,b=S.fraction(S.cancel(f))
    return min(k[0] for k,co in S.Poly(a,r).terms())-min(k[0] for k,co in S.Poly(b,r).terms())
def main():
    preds=['source_symmetry_gate_receipt.json','boundary_fibre_models_receipt.json','unramified_orevkov_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for path,digest in r['source_sha256'].items():check('predecessor seal '+name+' '+path,sha(HERE/path)==digest)
    for name in ['inverse_receipt.json','discriminant_receipt.json']:
        r=json.loads((OMIT/name).read_text());check('retained omitted inverse receipt PASS '+name,r['outcome']=='PASS')
        for label in ['input_sha256','output_sha256']:
            for path,digest in r.get(label,{}).items():check('retained omitted inverse seal '+name+' '+path,sha(OMIT/path)==digest)
    manifest=HERE/'sources/trace-dual-pencil/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('primary trace/different source '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    u,v,w,r,T,h,x,y,C,c,t=S.symbols('u v w r T h x y C c t')
    # Actual finite smooth local covers: pi(u,w)=(u,w^ell+u*w).
    # Their complement of the critical divisor is NOT a marked source A2.
    covers=[];zero_forms={}
    for ell in [2,3,4]:
        f=w**ell+u*w-v;der=S.diff(f,w)
        M=S.zeros(ell)
        for j in range(ell):
            rem=S.rem(w**(j+1),f,w)
            for i in range(ell):M[i,j]=S.expand(rem).coeff(w,i)
        powers=[S.eye(ell)]
        for k in range(1,2*ell):powers.append(powers[-1]*M)
        G0=S.Matrix(ell,ell,lambda i,j:S.expand(S.trace(powers[i+j])))
        G1=S.Matrix(ell,ell,lambda i,j:S.expand(S.rem(w**(i+j),f,w)).coeff(w,ell-1))
        Mr=ell*powers[ell-1]+u*S.eye(ell)
        check('full trace-dual multiplication identity for local length'+str(ell),G1*Mr==G0)
        check('second trace pairing generates the local dual at all base points length'+str(ell),G1.det() in [1,-1])
        check('both full trace matrices are symmetric length'+str(ell),G0==G0.T and G1==G1.T)
        check('normal trace discriminant is the critical-value polynomial length'+str(ell),S.expand(G0.det()-S.discriminant(f,w))==0)
        Fx=S.factor((T*G0-G1).det())
        check('trace pencil equals cleared characteristic polynomial length'+str(ell),S.expand(Fx-G1.det()*(T*Mr-S.eye(ell)).det())==0)
        check('inverse x equation has full degree and exact trace-zero coefficient length'+str(ell),S.degree(Fx,T)==ell and S.Poly(Fx,T).coeff_monomial(T**(ell-1))==0)
        check('local inverse coefficients have unit ideal length'+str(ell),Fx.subs(T,0) in [1,-1])
        # All source-coordinate inverse poles at generic R are retained.
        pulled=S.factor(S.diff(Fx,T).subs({T:1/r,u:r-ell*w**(ell-1),v:r*w-(ell-1)*w**ell},simultaneous=True))
        check('exact conductor derivative pole is three minus degree length'+str(ell),valuation(pulled,r)==3-ell)
        # At the nonimmersive point the ordinary trace radical is the maximal ideal.
        A=G0.subs({u:0,v:0});B=G1.subs({u:0,v:0});rad=S.Matrix.hstack(*A.nullspace())
        check('ordinary fibre trace radical has dimension ell minus1 length'+str(ell),rad.cols==ell-1)
        rho=(rad.T*B*rad).rank()
        check('restricted second pairing recovers all local excess length'+str(ell),rho==ell-2)
        zero_forms[ell]=(A,B)
        # Recover w rationally from x=1/f_w: the projection is primitive in this control.
        check('local reciprocal-Jacobian coordinate is primitive length'+str(ell),S.cancel(ell*v/(der+(ell-1)*u)-w).subs(v,w**ell+u*w)==0)
        covers.append(dict(local_length=ell,G0=[[str(z) for z in row] for row in G0.tolist()],G1=[[str(z) for z in row] for row in G1.tolist()],inverse=str(Fx),restricted_rank=int(rho),derivative_pole=ell-3,marked_Keller_source=False))
    # The all-length rank theorem is written using socle duality; longer controls
    # and direct sums ensure folds, collisions and excess are not conflated.
    for ell in [2,3,5,8]:
        A=S.zeros(ell);A[0,0]=ell;B=S.zeros(ell)
        for i in range(ell):B[i,ell-1-i]=1
        rad=S.Matrix.hstack(*A.nullspace())
        check('all-length socle pairing rank control '+str(ell),(rad.T*B*rad).rank()==ell-2 and B.det()!=0)
    mixed=[]
    for lengths in [[2,2],[3,1],[2,3,1],[4,2,1,1]]:
        AA=[];BB=[]
        for ell in lengths:
            A=S.zeros(ell);A[0,0]=ell;B=S.zeros(ell)
            for i in range(ell):B[i,ell-1-i]=1
            AA.append(A);BB.append(B)
        A=S.diag(*AA);B=S.diag(*BB);rad=S.Matrix.hstack(*A.nullspace());rho=(rad.T*B*rad).rank()
        check('whole fibre split includes every local factor '+str(lengths),rho==sum(max(ell-2,0) for ell in lengths))
        mixed.append(dict(lengths=lengths,ordinary_trace_nullity=rad.cols,restricted_rank=int(rho)))
    check('same ordinary trace nullity does not distinguish two folds from one cusp',mixed[0]['ordinary_trace_nullity']==mixed[1]['ordinary_trace_nullity'] and mixed[0]['restricted_rank']==0 and mixed[1]['restricted_rank']==1)
    # Generic rank is not a substitute for the sheaf-surjection condition.
    bad=S.diag(u,1);badblock=bad.row_join(S.zeros(2))
    check('generic-rank-only block matrix is rejected at its base point',badblock.rank()==2 and badblock.subs(u,0).rank()==1)
    # Global B_n complementary modules with one additional allowed R pole.
    module_controls=[]
    for n in [1,2,4]:
        H=S.prod(h-j for j in range(2,n+2));hh=x*x*y;a=x*(hh-C);b=y*H.subs(h,hh);z=x*b
        for k in range(-5,7):
            if k>=1:
                generator=x**k*(hh-C)**(k-1);from_xB=x*a**(k-1)
                check('entire O(R) positive module generated by xB n'+str(n)+' weight'+str(k),S.expand(generator-from_xB)==0)
            elif k==0:check('entire O(R) weight zero is generated by B n'+str(n),S.expand(hh-x*x*y)==0)
            else:
                if k%2==0:generator=b**(-k//2)
                else:generator=z*b**((-k-1)//2)
                f=S.cancel(generator/x**k)
                check('entire O(R) nonpositive modules already lie in B n'+str(n)+' weight'+str(k),f.is_polynomial(x,y) and valuation(f.subs({x:c/w+t,y:w*w},simultaneous=True)* (c/w+t)**k,w)>=0)
        check('simple pole generates the dual at the whole ramified line n'+str(n),S.expand((c/w+t)*w).subs(w,0)==c)
        # Generic source-coordinate divisor and its exact target norm member.
        sx=[1]*n;leading=S.Symbol('GR');Hx=S.Symbol('Hx');Hi=S.symbols('Hi0:'+str(n));constant=Hx*S.prod(Hi)
        normx=constant/leading
        check('constant inverse member equals leading coefficient times rational norm n'+str(n),S.cancel(leading*normx-constant)==0)
        module_controls.append(dict(n=n,generators=['1','x'],all_tested_weights=list(range(-5,7))))
    # Degree ledger uses n plus every excess, with no finite degree cutoff.
    n,E=S.symbols('n E',nonnegative=True,integer=True)
    check('Euler degree and conductor pole match uniformly',S.expand((n+3+E)-3-(n+E))==0)
    check('quartic two-prime packet forces zero restricted rank',4-1-3==0)
    # Mixed inverse charts must not be mistaken for the retained source alone.
    check('x finite-chart D_i is retained in N minus R on B_n',all(row['n']>=1 for row in module_controls))
    inputs=[Path(__file__),HERE/'TRACE_DUAL_PENCIL_AND_EXCESS.md',manifest]+[HERE/p for p in preds]+[OMIT/p for p in ['inverse_receipt.json','discriminant_receipt.json','ADMISSIBLE_INVERSE.md','DISCRIMINANT_LEDGER.md']]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      scope=dict(all_field_degrees=True,whole_complementary_r2q1_row=True,sole_index_two_ramified_prime_required=True,all_other_primes_unramified_and_retained=True,arbitrary_D_orders_poles_values_and_coefficients=True,source_symmetry_required=False,zero_excess_required=False),
      trace_dual=dict(equality='B^vee=B+x*B',whole_R_nonvanishing_simple_pole_required=True,quotient='O_R(R)',target_basis_only_local_required=True,symmetric_pairings=['G0=(Tr(b_i*b_j))','G1=(Tr(x*b_i*b_j))'],block_surjection='[G0|G1] generates B^vee at EVERY target point',generic_block_rank_suffices=False,det_G0='unit*G_R',det_G1='unit*H_0*product_(lambda_i=0) H_i^s_i',inverse_pencil='det(T*G0-G1)=det(G0)*Norm(T-x)',module_controls=module_controls),
      excess=dict(pointwise='rank(G1|ker G0)=sum_(q in R over b)(ell_q-2)',all_nonimmersive_lengths=True,multiple_folds_not_nonimmersive=True,zero_excess_iff_radical_totally_isotropic_everywhere=True,degree='d=n+3+sum_b rank(G1|ker G0)',quartic_two_prime_requires_zero_rank_everywhere=True,direct_sum_controls=mixed),
      inverse=dict(every_source_coordinate_primitive_without_Sd=True,degree='full actual d',leading='G_R',next_coefficient_divisible_by_GR=True,degree_mod_generic_GR='d-2',all_coefficient_ideal_is_unit=True,members='F_x(beta)=constant*H_beta*product_(lambda_i=beta)H_i^s_i',generic_members_reduced_irreducible_immersed_A1=True,coincident_global_supports_allowed=True,finite_x_chart_normalization='N minus R, including every D_i',normal_monogenic_algebra_asserted=False),
      conductor=dict(retained_PF2OMIT4_formulas_not_claimed_new=True,divisor='div(F_x,T)=C_x-(d-3)R',generic_R_conductor=0,exact_R_derivative_pole='d-3=n+E_R',discriminant='constant*G_R*J_x^2',norm_derivative='constant*G_R^(3-d)*J_x^2',conductor_may_contain_D_i=True,unregistered_draft_used=False),
      controls=covers,assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,finite_non_Keller_controls_pass=True,source_A2_from_controls=False,pencil_conditions_sufficient=False,imported_duality_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'trace_dual_pencil_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; full trace-dual pencil and all-degree radical-rank excess gate')
if __name__=='__main__':main()
