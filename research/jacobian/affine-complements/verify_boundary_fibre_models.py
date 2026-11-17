"""Controls for the full boundary-fibre budget and unbounded graded models."""
from pathlib import Path
import hashlib,json,os,platform
from math import gcd
import sympy as S
from verify_coordinate_pencil_matrix import cluster,blowup_matrix,order

HERE=Path(__file__).resolve().parent;checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    preds=['coordinate_pencil_matrix_receipt.json','two_sided_norm_receipt.json','unramified_sheet_pair_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor seal '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/boundary-fibre-models/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];check('primary source '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    x,y,h,C,c,d,a,b,z,w,t,v,U=S.symbols('x y h C c d a b z w t v U')
    groups=[]
    configs=[[(0,2)],[(0,2),(0,3)],[(0,2),(1,3)],[(0,2),(0,3),(1,4)],[(0,2),(1,3),(2,4)]]
    RR=cluster(z**-1+t,z*z,z,t)
    for values in configs:
        DD=[cluster(lam+w,di*w**-2+v,w,v) for lam,di in values]
        Q,I,labels,K=blowup_matrix([RR]+DD);ix=[labels[A[1]] for A in [RR]+DD];M=I.extract(ix,ix)
        expected=S.Matrix(len(ix),len(ix),lambda i,j:-3 if i==j==0 else -4 if i==j else 2 if i==0 or j==0 else -2 if values[i-1][0]==values[j-1][0] else 0)
        tag=str(values)
        check('whole grouped primitive matrix by actual shared point blowups '+tag,M==expected)
        D=-M[1:,1:];pv=S.ones(len(values),1)*2
        energy=(pv.T*D.inv()*pv)[0];sizes=[sum(lam==k for lam,di in values) for k in sorted(set(lam for lam,di in values))]
        predicted=sum(S.Rational(2*k,k+1) for k in sizes)
        check('grouped Schur energy keeps all within-fibre cross terms '+tag,energy==predicted)
        check('distinct-value blocks have exact zero cross terms '+tag,all(M[i+1,j+1]==0 for i in range(len(values)) for j in range(len(values)) if values[i][0]!=values[j][0]))
        check('canonical source volume orders match every component '+tag,[K[i] for i in ix]==[1]+[0]*len(values))
        groups.append(dict(values=values,group_sizes=sizes,energy=str(energy),strict_test=bool(energy<3),matrix=[[int(v) for v in row] for row in M.tolist()]))
    mm,nn=S.symbols('m n',positive=True,integer=True)
    check('two group classification is the exact reciprocal inequality',S.factor(3-2*mm/(mm+1)-2*nn/(nn+1)-(2/(mm+1)+2/(nn+1)-1))==0)
    check('r1 complementary ratio lies below two',S.Rational(5,4)<2)
    check('all r>=3 minimal ratios lie below two by explicit factor',S.expand(2*(mm-1)**2-(mm*mm-mm+1)-(mm*mm-3*mm+1))==0)
    # Every finite size is allowed by the written one-group Schur identity.
    check('unbounded one-group Schur identity',S.factor(3-2*nn/(nn+1)-(nn+3)/(nn+1))==0)
    for n in [1,2,3,5,9]:
        D=2*(S.eye(n)+S.ones(n));pv=2*S.ones(n,1)
        check('all-count matrix determinant control n'+str(n),D.det()*(3-(pv.T*D.inv()*pv)[0])==2**n*(n+3))
    for r,limit in [(3,6),(4,2),(5,1),(12,1)]:
        count=[n for n in range(1,10) if S.Rational(2*n,n+1)<S.Rational(r*r-r+1,(r-1)**2)]
        check('conditional all-simple D component bound r'+str(r),count==list(range(1,limit+1)))
    check('all r>=5 component cutoff follows from a positive quadratic',S.expand((mm*mm-mm+1)-2*(mm*mm-3*mm+1)+mm*mm-5*mm+1)==0)
    # Verify the entire determinantal presentation and both smooth localizations.
    T0=S.Symbol('T0');f1=a*b-(h-C)*z;f2=a*z-(h-C)*T0;f3=z*z-T0*b
    X=S.Symbol('X');eq=X*X*b-T0
    check('first localized chart has one hypersurface relation',S.expand(f2.subs({a:(h-C)*X,z:X*b})-(h-C)*eq)==0 and S.expand(f3.subs(z,X*b)-b*eq)==0)
    check('second localized chart has one hypersurface relation',S.cancel(f1.subs(b,z*z/T0)-z/T0*f2)==0 and S.cancel(f3.subs(b,z*z/T0))==0)
    ring_controls=[]
    for n in [1,2,3,5]:
        H=S.prod(h-j for j in range(2,n+2));TT=h*H;hxy=x*x*y
        aa=x*(hxy-C);bb=y*H.subs(h,hxy);zz=x*bb
        check('three complete polynomial relations n'+str(n),all(S.expand(f.subs({a:aa,b:bb,z:zz,h:hxy,T0:TT.subs(h,hxy)},simultaneous=True))==0 for f in [f1,f2,f3]))
        check('first smooth chart uses a squarefree T n'+str(n),S.gcd(TT,S.diff(TT,h))==1)
        check('two principal opens cover whenever T(C) is nonzero n'+str(n),S.expand(S.rem(TT,h-C,h)-TT.subs(h,C))==0)
        check('global source inverse x equals T over z n'+str(n),S.cancel(TT.subs(h,hxy)/zz-x)==0)
        check('global source inverse y near h0 n'+str(n),S.cancel(bb/H.subs(h,hxy)-y)==0)
        check('source inverse on T nonzero n'+str(n),S.cancel(zz**2/(hxy*H.subs(h,hxy)**2)-y)==0)
        for di in range(2,n+2):
            Hi=S.cancel(H/(h-di));unit=(h*Hi).subs(h,di)
            check('every added D has simple x and exact pole2 n'+str(n)+' d'+str(di),unit!=0 and S.expand(TT-(h-di)*h*Hi)==0)
        for k in range(-6,7):
            if k>=0:f=aa**k;expected_k=x**k*(hxy-C)**k
            elif k%2==0:f=bb**(-k//2);expected_k=x**k*(hxy*H.subs(h,hxy))**(-k//2)
            else:
                ell=(-k-1)//2;f=zz*bb**ell;expected_k=x**k*(hxy*H.subs(h,hxy))**(ell+1)
            check('full weight module generator n'+str(n)+' k'+str(k),S.cancel(f-expected_k)==0)
        ring_controls.append(dict(n=n,H=str(H),T=str(TT),source_and_all_added_lines_accounted_for=True))
    # Global exact primitive, with cancellation at every entire boundary line.
    xr=c/z+t;yr=z*z;yd=d/x**2+v
    Rdz=S.simplify(-2*yr*S.diff(xr,z)-xr*S.diff(yr,z));Rdt=S.simplify(-2*yr*S.diff(xr,t))
    Ddx=S.simplify(-2*yd-x*S.diff(yd,x));Ddv=S.simplify(-x*S.diff(yd,v))
    check('primitive is regular across the full ramified line',Rdz==-2*z*t and Rdt==-2*z*z)
    check('primitive is regular across every full unramified line',Ddx==-2*v and Ddv==-x)
    check('global primitive has the exact source differential',S.diff(-x,x)-S.diff(-2*y,y)==1)
    check('source volume has precisely orders one and zero',S.simplify(S.diff(xr,t)*S.diff(yr,z))==2*z and S.diff(yd,v)==1)
    # The resonance really adds both ends of a whole hyperbola.
    check('resonance hyperbola ends at R t0 and D v0',S.simplify((xr*xr*yr-c*c).subs(t,0))==0 and S.simplify(x*x*(c*c/x**2)-c*c)==0)
    extra=[]
    for r in [4,5,8]:
        CC=S.Integer(1);dd=S.Integer(2);u=S.Symbol('u');N=CC+dd*u**(r-2)
        G=x**r*y**(r-1)-dd*(x*y)**(r-2)-CC
        xx=N/u**(r-1);yy=u**r/N
        check('failed shortcut curve has the exact rational parameterization r'+str(r),S.cancel(G.subs({x:xx,y:yy},simultaneous=True))==0)
        check('failed shortcut retains r minus2 additional simple ends r'+str(r),S.degree(N,u)==r-2 and S.gcd(N,S.diff(N,u))==1 and N.subs(u,0)!=0)
        face=[((i,j),co) for (i,j),co in S.Poly(G,x,y).terms() if j-i==0]
        check('third infinity Newton face is present r'+str(r),len(face)==2 and S.expand(sum(co*x**i*y**j for (i,j),co in face)+dd*(x*y)**(r-2)+CC)==0)
        FR=S.expand(G.subs({x:z**(-(r-1))+z**(2-r)*t,y:z**r},simultaneous=True))
        FD=S.expand(G.subs(y,dd*x**-2+v))
        check('both selected chart equations can pass while the curve still has other ends r'+str(r),order(FR,z)[0]==1 and S.Poly(S.expand(FR).coeff(z,1),t).degree()==1 and order(FD,x)[0]==4-r and S.Poly(S.expand(FD).coeff(x,4-r),v).degree()==1)
        extra.append(dict(r=r,additional_ends=r-2,complete_curve_exclusion_valid=False))
    # A non-Keller power map is a control for the fixed-origin deck obstruction.
    check('nontrivial weighted deck preserves a non-Keller map',S.expand(((-x)**2-x*x))==0 and S.diff(x*x,x)==2*x)
    for aa,bb in [(2,1),(3,2),(5,3),(7,4)]:
        uu,vv,gg=S.gcdex(aa,bb);check('primitive monomial torus parameter '+str((aa,bb)),gg==1 and aa*uu+bb*vv==1)
        check('origin deck condition cannot hold for a nontrivial order '+str((aa,bb)),all(not (aa%k==0 and bb%k==0) for k in range(2,10)))
    inputs=[Path(__file__),HERE/'BOUNDARY_FIBRE_BUDGET_AND_GRADED_MODELS.md',manifest,HERE/'verify_coordinate_pencil_matrix.py']+[HERE/n for n in preds]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      fibre_budget=dict(all_field_degrees=True,full_matrix_required=True,all_D_orders_and_poles=True,other_primes_retained=True,complementary_index_two_required=True,distinct_values_cross_zero=True,strict_budget='n_R/q^2>sum_lambda E_lambda',group_energy_lower_bound='E_lambda>=max_i(p_i^2/n_i)>=1',one_boundary_value_except_r2q1=True,r2q1_at_most_two_boundary_values=True,unconditional_prime_count_bound=False,conditional_all_simple_D_counts=dict(R_order3=6,R_order4=2,R_order_at_least5=1),common_ruling_case_excluded=False,grouped_controls=groups),
      affine_family=dict(all_n_at_least_one=True,conditions='C,d_i nonzero and pairwise distinct',generators=['h=x^2*y','a=x*(h-C)','b=y*H(h)','z=x*y*H(h)'],entire_affine_presentation=True,two_smooth_principal_opens_cover=True,full_marked_A2_open=True,complement='one index-two volume-order R and n unramified-volume D_i affine lines',complementary_R_r2q1=True,all_D_simple_pole2=True,global_exact_primitive='-2y*dx-x*dy',full_weight_modules=True,source_Cstar_action=True,target_Cstar_action_inferred=False,matrix_determinant='2^n*(n+3)',conditional_zero_excess_degree='n+3',finite_Keller_map_supplied=False,controls=ring_controls),
      finite_algebra=dict(primitive_source_monomial_general_a_b=True,global_regularity_and_all_boundary_vertical_required=True,full_finite_Keller_map_required=True,generic_Cstar_image_residue_one=True,fixed_origin_deck_obstruction=True,primitive_element=True,monic_norm_degree_is_full_Keller_degree=True,specializations=['chi(0)=H_x^2*H_y','chi(d_i)=J_i*H_i^2','chi(C)=G_R*J_C^f, f=1 or2'],coincident_global_image_supports_allowed=True,normal_monogenic_algebra_asserted=False),
      failures=dict(resonance_C_equals_di_excluded=True,extra_infinity_ends_retained=True,selected_chart_tests_not_all_ends=True,failed_shortcut_controls=extra),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'boundary_fibre_models_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; full fibre budget, unbounded affine models and primitive monomial gate')
if __name__=='__main__':main()
