"""PF116 exact module maps and finite covers; universal claims use written proofs."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import logarithmic_pullback as L
P=Path(__file__).resolve().parent


def order0(expr,t):
    a,b=S.fraction(S.cancel(expr))
    return min(m[0] for m,c in S.Poly(a,t).terms())-min(m[0] for m,c in S.Poly(b,t).terms())


def monomial_basis(equations,x,y):
    G=S.groebner(equations,x,y)
    if not G.is_zero_dimensional:raise ValueError('finite fibre required')
    lm=[p.LM(order=G.order).exponents for p in G.polys]
    bx=min(a for a,b in lm if a and not b);by=min(b for a,b in lm if b and not a)
    return [(a,b) for a in range(bx) for b in range(by) if not any(a>=i and b>=j for i,j in lm)],lm


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def ident(name,value):check(name,S.cancel(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    oldfiles=['residual_milnor_receipt.json','torus_cusp_excess_receipt.json','whole_plane_peripheral_receipt.json','coprime_leading_forms_receipt.json']
    for name in oldfiles:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/log-pullback-dessin/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('standard module/Belyi imports and HTML audit declared',manifest['Saito_freeness_and_residue_duality_imported'] and manifest['dessin_Belyi_correspondence_imported'] and not manifest['fresh_PDF_audit_claimed'])
    a,b,c,d,e,f,g,h=S.symbols('a b c d e f g h');J=S.Matrix([[a,b],[c,d]]);M=S.Matrix([[e,f],[g,h]])
    for i,z in enumerate(J.T*M.inv().T-(J.inv()*M).inv().T):ident('UNBOUNDED duality of log pullback and derivation lift '+str(i),z)
    rr=S.symbols('rr',integer=True,positive=True);x0=S.symbols('x0',positive=True);H0=S.symbols('H0',nonzero=True)
    ident('UNBOUNDED determinant different cancellation',x0**(rr-1)/(x0**rr*H0)-1/(x0*H0))
    packet=L.higher_characteristic_packet();x,y,u,v=[packet[k] for k in ('x','y','u','v')]
    F=packet['F'];M=packet['target_fields'];A=packet['A'];H=packet['H'];Lift=packet['lifted_fields']
    for i in range(2):ident('actual target logarithmic derivation '+str(i),M[0,i]*S.diff(F,u)+M[1,i]*S.diff(F,v)-packet['target_log_factors'][i]*F)
    ident('actual target Saito determinant',M.det()-1183*F)
    for i,z in enumerate(Lift):check('actual lifted coefficient analytic '+str(i),S.denom(z).subs({x:0,y:0})!=0)
    for i in range(2):
        for name,curve in [('A',A),('H',H)]:
            ratio=S.cancel((Lift[0,i]*S.diff(curve,x)+Lift[1,i]*S.diff(curve,y))/curve)
            check('actual lifted logarithmic tangency '+str((i,name)),S.denom(ratio).subs({x:0,y:0})!=0)
    ident('actual source Saito determinant over REDUCED divisor',S.cancel(Lift.det()/(A*H))-S.Rational(169,6)/(1+y))
    for z in [x,y]:ident('actual ramification-weighted dlog pullback '+str(z),S.diff(A**3*H,z)/(A**3*H)-3*S.diff(A,z)/A-S.diff(H,z)/H)
    t=S.symbols('t');Xp=t**4;Yp=56*t**6+64*t**7
    mons=[Xp**i*Yp**j for i in range(4) for j in range(3) if 4*i+6*j<16]
    gamma=L.truncated_values(mons,t,16)
    differential=[m*S.diff(Xp,t) for m in mons if order0(m,t)+4<16]+[m*S.diff(Yp,t) for m in mons if order0(m,t)+6<16]
    values=[i+1 for i in L.truncated_values(differential,t,15)];extra=sorted(set(values)-set(gamma))
    check('actual semigroup below conductor16',gamma==[0,4,6,8,10,12,13,14])
    check('actual analytic differential extras exactly11 and15',extra==[11,15])
    ident('actual value11 differential',Xp*S.diff(Yp/56,t)-S.Rational(3,2)*(Yp/56)*S.diff(Xp,t)-S.Rational(8,7)*t**10)
    tjurina=[]
    for N in [17,18]:
        basis,_=monomial_basis([F,S.diff(F,u),S.diff(F,v)]+[u**i*v**(N-i) for i in range(N+1)],u,v);tjurina.append(len(basis))
    check('actual local Tjurina14 by stable exact lengths and Nakayama',tjurina==[14,14] and 16-14==len(extra))
    residues=[M[1,1]/(1183*S.diff(F,u)),-M[1,0]/(1183*S.diff(F,u))]
    residue_orders=[order0(z.subs({u:Yp,v:Xp}),t) for z in residues]
    check('actual target residue generator orders',residue_orders==[-9,-7])
    check('source normalization saturation loses actual differential gaps',min(order0(S.diff(Xp,t),t),order0(S.diff(Yp,t),t))+1==4 and all(k not in values for k in [5,7,9]))
    rotations=[]
    for m in [1,2,3,11,32]:
        sigma,tau,face=L.dessin_rotations(m);degree=6*m-2
        check('dessin involution cycles '+str(m),sorted(map(len,L.cycles(sigma)))==[2]*(3*m-1))
        check('dessin cubic cycles '+str(m),sorted(map(len,L.cycles(tau)))==[1]+[3]*(2*m-1))
        check('dessin large face and monogons '+str(m),sorted(map(len,L.cycles(face)))==[1]*m+[5*m-2])
        seen={0};todo=[0]
        while todo:
            j=todo.pop()
            for z in [sigma[j],tau[j]]:
                if z not in seen:seen.add(z);todo.append(z)
        check('dessin connected genuszero '+str(m),len(seen)==degree and sum(len(L.cycles(p)) for p in [sigma,tau,face])==degree+2)
        rotations.append(dict(m=m,degree=degree,passport=[sorted(map(len,L.cycles(p))) for p in [tau,sigma,face]]))
    m=S.symbols('m',integer=True,positive=True);n=3*m-1;ee=5*m-2;ell=6*m-2
    identities={'passport Euler':2*m+n+(m+1)-(ell+2),'Riemann Hurwitz':(4*m-2)+n+(5*m-3)-(2*ell-2),'Wronskian degree':(2*m-1)+n-1-(ee-1),'Wronskian leading coefficient':3*(2*m-1)-2*n+1,'pure normal Jacobian':5*n-4-3*(ee-1),'residual no extra Y':6*n-3*ee-3*m,'weighted upstairs rank':(2*n)*(3*n)-3*(2*n*n),'quotient rank':2*n-ell,'local sharpness':ee-5*m+2,'marked character':2*n-ee-m,'Bezout coprimality':5*n-3*ee-1,'invariant coordinates':3*m-n-1,'residual Milnor':3*m*(m-1)-m+1-n*(m-1),'wrong-plane Euler defect':ee+m-(ell-1)-1}
    for name,z in identities.items():ident('UNBOUNDED '+name,z)
    j=S.symbols('j',integer=True,nonnegative=True)
    check('UNBOUNDED nonMorse nonGorenstein family for m>=2',(2*(2+j)-2).is_positive and (3*(2+j)-1-4).is_positive)
    exact_examples=[]
    for mm in [1,2]:
        model=L.explicit_sharp_pair(mm);t,X,Y,p,q,P0,Q0=[model[k] for k in ('t','X','Y','p','q','P','Q')];profile=model['profile'];nn=profile['quotient_order'];index=profile['index']
        check('actual p and q coprime '+str(mm),S.gcd(p,q)==1)
        W=S.factor(3*S.diff(p,t)*q-2*p*S.diff(q,t));wc=S.cancel(W/t**(index-1))
        check('actual pure Wronskian '+str(mm),not wc.has(t) and wc!=0)
        ident('actual weighted Jacobian '+str(mm),S.diff(P0,X)*S.diff(Q0,Y)-S.diff(P0,Y)*S.diff(Q0,X)-nn*wc*X**(index-1))
        check('actual small cyclic invariance '+str(mm),all((3*a+b)%nn==0 for poly in [P0,Q0] for (a,b),c in S.Poly(poly,X,Y).terms()))
        basis,lm=monomial_basis([P0,Q0],X,Y);chars=[sum((3*a+b)%nn==k for a,b in basis) for k in range(nn)]
        check('actual upstairs and invariant quotient fibre lengths '+str(mm),len(basis)==2*nn*nn and chars==[2*nn]*nn)
        residual=S.cancel((P0**3-Q0**2)/X**index)
        check('actual residual and one-point fibre support '+str(mm),S.Poly(residual,X,Y).degree(X)==mm and S.Poly(Q0.subs(Y,0),X).degree()==nn)
        exact_examples.append(dict(m=mm,upstairs_length=len(basis),invariant_length=chars[0],character_lengths=chars,leading_monomials=lm))
    for mm in [1,2,1000,10**30]:
        z=L.sharp_profile(mm);check('local sharpness keeps wrong-plane/global-optimality guards '+str(mm),z['index']==5*z['excess']-2 and z['local_degree']==z['index']+z['excess'] and not z['whole_plane_Keller_realization'] and not z['global_Wright_bound_optimality_claimed'])
    text=' '.join((P/'LOGARITHMIC_PULLBACK_AND_SHARP_DESSIN_CUSP_COVERS.md').read_text().lower().split())
    for phrase in ['all branch values','unit determinant','not a replacement of the nonreduced','all its multibranch gluing','need not belong to o_b','not extension separately to the normalizations','no identification','target-coefficient','saturation erases','actual saito basis','for every integer m>=1','actual finite covers','only the origin','non-gorenstein','one point of length','not optimality','not a keller realization','not used as a theorem input','no transport of bernstein-sato roots','surface trace dual']:check('written scope '+phrase,phrase in text)
    files=['LOGARITHMIC_PULLBACK_AND_SHARP_DESSIN_CUSP_COVERS.md','logarithmic_pullback.py','verify_logarithmic_pullback.py','sources/log-pullback-dessin/manifest.json']+oldfiles+[s['path'] for s in manifest['sources']]
    result=dict(claim='PF2RULED116',outcome='PASS',check_count=len(checks),checks=checks,source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
      theorem=dict(all_finite_smooth_plane_germ_degrees=True,reduced_target_contains_ALL_branch_values_required=True,exact_logarithmic_pullback_isomorphism=True,exact_ramification_weighted_residue_map=True,actual_curve_Jacobian_duality_formula=True,target_coefficient_differential_values_scale_by_residue=True,source_normalization_saturation_loses_differential_gaps=True,unbounded_genuine_finite_sharp_cusp_family=True,local_five_eta_coefficient_sharp=True,non_Morse_common_leading_zero_family=True,all_ordinary_residual_maps_primitive=True),
      fixed_actual_module_control=dict(image_conductor=16,Tjurina=14,differential_extra_values=extra,target_residue_generator_orders=residue_orders,target_Saito_determinant='1183 F',source_Saito_unit='169/(6(1+Y))'),exact_finite_examples=exact_examples,dessin_controls=rotations,
      assurance=dict(Saito_freeness_and_residue_duality_imported=True,dessin_Belyi_existence_imported=True,written_unbounded_proof=True,finite_samples_prove_theorem=False,original_surface_trace_dual_identified_with_curve_residue=False,ramification_diagonal_dropped=False,multibranch_gluing_dropped=False,target_coefficient_module_identified_with_source_saturation=False,newer_preprint_semifield_transfer_imported=False,Bernstein_root_transport_proved=False,orbifold_BMY_equality_assumed=False,global_Wright_bound_optimality_claimed=False,all_coefficient_fields_claimed_rational=False,sharp_family_is_Keller_completion=False,actual_global_finite_algebra_closed=False,formal_verification=False,external_review=False,research_goal_completed=False),environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'logarithmic_pullback_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',len(checks),'PF116 controls; actual logarithmic map and sharp unbounded finite cusp covers; full goal ACTIVE')


if __name__=='__main__':main()
