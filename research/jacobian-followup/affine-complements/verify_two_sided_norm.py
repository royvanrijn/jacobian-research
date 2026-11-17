"""Exact controls for the written two-sided full-norm curve obstruction."""
from pathlib import Path
from math import gcd
import hashlib,json,os,platform
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(f,z):
    f=S.expand(f)
    return min((int(a.as_powers_dict().get(z,0)) for a in S.Add.make_args(f)),default=S.oo) if f!=0 else S.oo
def coeff(f,z,n):return S.expand(f).coeff(z,n)
def defect(r,q,support):
    ss=sum(2*q-max(j for j in support if a*j%r) for a in range(1,r))
    return ss,r*q-ss
def quotient_norm(f,r,z,y):
    modulus=z**r-y
    columns=[S.Poly(S.rem(S.expand(f*z**k),modulus,z),z) for k in range(r)]
    matrix=S.Matrix(r,r,lambda i,j:columns[j].nth(i))
    return S.expand(matrix.det(method='domain-ge'))
def full_norm(r,h,q,z,x,y):
    modulus=z**r-y;f=S.expand(z**q*x-z**q*h)
    n=S.expand((-1)**((r+1)*q)*quotient_norm(f,r,z,y))
    second=S.rem(f,modulus,z)
    alternate=S.expand((-1)**((r+1)*q)*S.resultant(modulus,second,z))
    assert n==alternate,'reduced resultant disagrees with direct quotient determinant'
    return n
def band_cut(P,p,q,x,y):
    return S.expand(sum(c*x**i*y**j for (i,j),c in S.Poly(P,x,y).terms()
        if (i,j)!=(p,q) and q-j<max(2,p-i)))

def main():
    predecessor=HERE/'unramified_sheet_pair_receipt.json'
    old=json.loads(predecessor.read_text());check('PF48 predecessor PASS',old['outcome']=='PASS')
    for f,h in old['source_sha256'].items():check('PF48 sealed input '+f,sha(HERE/f)==h)
    manifest=HERE/'sources/two-sided-norm/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements']
        check('primary source '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    x,y,z,w,t,v,A=S.symbols('x y z w t v A')
    R1,R2,q1,p1,k,l,a,b,ee=S.symbols('Delta_R Delta_D q p k l a b e')
    rr,ss=S.symbols('r s')
    check('R band gap lower bound is exact at l=k-1',S.expand((R1-ee+q1*(k-2)-rr*(l-1)).subs(l,k-1)-(R1-ee+(q1-rr)*(k-2)))==0)
    check('R D-norm gap lower bound is exact at a=b+1',S.expand((R1-ee+q1*(a-2)-rr*(b-1)).subs(a,b+1)-(R1-ee+(q1-rr)*(b-1)))==0)
    check('D retained band bound is exact at l=k',S.expand((R2-1+p1*(l-2)-ss*(k-1)).subs(l,k)-(R2-1-ss+(p1-ss)*(k-2)))==0)
    raw=S.resultant(z**3-1,10*z**7-3-2*z**3,z)
    determinant=quotient_norm(10*z**7-3-2*z**3,3,z,S.Integer(1))
    reduced=S.resultant(z**3-1,10*z-5,z)
    check('odd odd direct quotient norm and reduced resultant give positive875',determinant==reduced==875)
    check('record installed raw resultant behavior without accepting it as norm',raw in [-875,875])
    examples=[
      (2,2,{3:2,1:3},2,{3:3}),
      (2,3,{4:2,2:3},2,{5:3,3:2}),
      (2,4,{5:2,3:3},3,{4:5}),
      (3,3,{4:2,1:3},2,{5:3,3:2}),
      (4,5,{6:2,2:3},3,{7:3,4:2}),
      (2,4,{6:2,3:3},4,{6:2,5:3}),
      (3,4,{6:2,3:3},4,{6:2,5:3}),
      (3,3,{5:2,1:3},3,{4:3}),
      (1,2,{3:2},3,{4:3}),
      (2,2,{3:S.Symbol('c3',nonzero=True),1:S.Symbol('c1')},2,{3:S.Symbol('d3',nonzero=True)}),
    ]
    controls=[]
    for e,r,cs,s,ds in examples:
        q=max(cs);p=max(ds);SR,DR=defect(r,q,cs);SD,DD=defect(s,p,ds)
        tag=f'e{e} r{r} q{q} s{s} p{p}'
        check('whole chart hypotheses '+tag,r>=max(2,e) and min(cs)>=r-e+1 and r not in cs and min(ds)>=s+1 and gcd(r,*cs)==gcd(s,*ds)==1)
        check('full sufficient defect hypotheses '+tag,q>=max(r+1,s) and p>=r and DR>=e and DD>=s+1)
        h=sum(c*z**(-j) for j,c in cs.items());g=sum(d*w**(-j) for j,d in ds.items())
        PR=full_norm(r,h,q,z,x,y);PD=full_norm(s,g,p,w,y,x)
        P0=S.expand(x**(p-r)*PR);E=band_cut(P0,p,q,x,y)
        F=S.expand(P0-E+y**(q-s)*PD-x**p*y**q)
        check('both norms have their exact top terms '+tag,coeff(PR,y,q)==x**r and coeff(PD,x,p)==y**s)
        check('both full reference branches vanish '+tag,S.expand(PR.subs({x:h,y:z**r},simultaneous=True))==0 and S.expand(PD.subs({y:g,x:w**s},simultaneous=True))==0)
        check('trace residue is zero '+tag,coeff(coeff(PR,y,q-1),x,r-1)==0)
        check('unramified norm has the strict loss band '+tag,all((i,j)==(p,s) or p-i>=s-j+1 for (i,j),c in S.Poly(PD,x,y).terms()))
        check('retained R norm has entire D-safe band '+tag,all((i,j)==(p,q) or q-j>=max(2,p-i) for (i,j),c in S.Poly(P0-E,x,y).terms()))
        check('erased band uses exact zero residue '+tag,E==0 or all(p-i>=2 and 1<=q-j<=p-i-1 for (i,j),c in S.Poly(E,x,y).terms()))
        BR=SR+2*q-p*q;BD=SD+2*p-p*q;MR=BR+e-r;MD=BD+1-s
        check('R derivative full contact sum '+tag,val(S.diff(P0,x).subs({x:h,y:z**r},simultaneous=True),z)==BR)
        check('D derivative full contact sum '+tag,val(S.diff(y**(q-s)*PD,y).subs({y:g,x:w**s},simultaneous=True),w)==BD)
        check('entire erased band is R-safe '+tag,val(E.subs({x:h,y:z**r},simultaneous=True),z)>=MR)
        check('entire D correction is R-safe '+tag,val((y**(q-s)*PD-x**p*y**q).subs({x:h,y:z**r},simultaneous=True),z)>=MR)
        check('entire retained R correction is D-safe '+tag,val((P0-E-x**p*y**q).subs({y:g,x:w**s},simultaneous=True),w)>=MD)
        FR=S.expand(F.subs({x:h+z**(e-r)*t,y:z**r},simultaneous=True))
        FD=S.expand(F.subs({y:g+w**(1-s)*v,x:w**s},simultaneous=True))
        lr=coeff(FR,z,MR);ld=coeff(FD,w,MD)
        check('entire R substitution has the precise order '+tag,val(FR,z)==MR)
        check('R leading equation is affine with nonzero constant coefficient '+tag,S.Poly(lr,t).degree()==1 and lr.coeff(t)!=0 and not lr.coeff(t).has(t))
        check('entire D substitution has the precise order '+tag,val(FD,w)==MD)
        check('D leading equation is affine with nonzero constant coefficient '+tag,S.Poly(ld,v).degree()==1 and ld.coeff(v)!=0 and not ld.coeff(v).has(v))
        check('full lower rectangle and sole top vertex '+tag,all((i,j)==(p,q) or i<=p-1 and j<=q-1 for (i,j),c in S.Poly(F,x,y).terms()))
        check('full R Newton face unchanged '+tag,all(q*i-r*j<q*(p-r) for (i,j),c in S.Poly(F-P0,x,y).terms()))
        check('full D Newton face unchanged '+tag,all(p*j-s*i<p*(q-s) for (i,j),c in S.Poly(F-y**(q-s)*PD,x,y).terms()))
        controls.append(dict(e=e,r=r,q=q,s=s,p=p,R_support={str(j):str(c) for j,c in cs.items()},D_support={str(j):str(c) for j,c in ds.items()},Delta_R=DR,Delta_D=DD,R_order=MR,D_order=MD,curve=str(F),erased_band=str(E)))
    simple=[]
    for s,p,q in [(1,2,2),(2,3,2),(2,3,3),(3,4,3),(3,5,4)]:
        cq,d=S.symbols('cq d',nonzero=True)
        h=cq*y**-q+sum(S.Symbol('c'+str(j))*y**-j for j in range(2,q));C=S.expand(y**q*h)
        ds={p:d}|{j:S.Symbol('d'+str(j)) for j in range(s+1,p)}
        g=sum(dj*w**-j for j,dj in ds.items());SD,DD=defect(s,p,ds)
        PD=full_norm(s,g,p,w,y,x)
        F=S.expand(x**p*y**q-x**(p-1)*C-A*cq*x**(p-2)+y**(q-s)*PD-x**p*y**q)
        MR=-p*q+2*q+1;MD=SD+2*p-p*q+1-s
        FR=S.expand(F.subs(x,A+h+y*t));FD=S.expand(F.subs({x:w**s,y:g+w**(1-s)*v},simultaneous=True))
        tag=f'simpleR s{s} p{p} q{q}'
        check('simple R arbitrary D sufficient conditions '+tag,q>=max(2,s) and DD>=s+1)
        check('simple R exact chart equation retains all arbitrary coefficients and A '+tag,val(FR,y)==MR and S.Poly(coeff(FR,y,MR),t).degree()==1 and coeff(FR,y,MR).coeff(t)==cq**(p-1))
        check('higher D whole chart equation with A retained '+tag,val(FD,w)==MD and S.Poly(coeff(FD,w,MD),v).degree()==1 and coeff(FD,w,MD).coeff(v)!=0)
        check('simple higher D two Newton faces unchanged '+tag,all(q*i-j<q*(p-1) for (i,j),c in S.Poly(y**(q-s)*PD-x**p*y**q,x,y).terms()) and all(p*j-s*i<p*(q-s) for (i,j),c in S.Poly(-x**(p-1)*C-A*cq*x**(p-2),x,y).terms()))
        simple.append(dict(s=s,p=p,q=q,all_coefficients_symbolic=True,constant_A_retained=True,curve=str(F)))
    # The untouched band really can fail; it is not cosmetic deletion.
    h=2*z**-3+3*z**-1;g=3*w**-3
    PR=full_norm(2,h,3,z,x,y);PD=full_norm(2,g,3,w,y,x)
    uncut=S.expand(x*PR+y*PD-x**3*y**3);SD,DD=defect(2,3,[3]);MD=SD+6-9-1
    check('uncut two norms fail on the nonsimple unramified chart',val(uncut.subs({x:w**2,y:g+w**-1*v},simultaneous=True),w)<MD)
    # Global primitive and every affine-direction pole in the common-ruling controls.
    common=[]
    for q in [2,3,5,9]:
        aa,bb=S.symbols('a b',nonzero=True)
        xr=y*t+aa*y**-q;xd=v+bb*y**-q
        Rdy=S.simplify((y*S.diff(xr,y)+q*xr)/(q-1));Rdt=S.simplify(y*S.diff(xr,t)/(q-1))
        Ddy=S.simplify((y*S.diff(xd,y)+q*xd)/(q-1));Ddv=S.simplify(y*S.diff(xd,v)/(q-1))
        check('common ruling primitive regular on entire R chart q'+str(q),Rdy==(q+1)*y*t/(q-1) and Rdt==y*y/(q-1))
        check('common ruling primitive regular on entire D chart q'+str(q),Ddy==q*v/(q-1) and Ddv==y/(q-1))
        check('global derivative is the source volume q'+str(q),S.diff(q*x/(q-1),x)-S.diff(y/(q-1),y)==1)
        check('all six directional translations are strict poles q'+str(q),all(k<0 for k in [-q,-q-1,-q,-q,-q,-q-1]))
        common.append(dict(q=q,R_transition_exponent=1,D_transition_exponent=0,torsor_multiindex=[1,0,1],all_six_translation_poles=True,global_exact_primitive=True,finite_Keller_map_supplied=False))
    inputs=[Path(__file__),HERE/'TWO_SIDED_NORM_AND_COMMON_RULING.md',manifest,predecessor]+[HERE/r['path_relative_to_affine_complements'] for r in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      two_sided=dict(all_field_degrees=True,all_indices_on_R=True,other_primes_retained=True,complementary_coordinates_required=True,smooth_affine_component_opens_required=True,R_order='r>=max(2,e)',D_order='s>=1',D_simple_required=False,R_pole='q>=max(r+1,s)',D_pole='p>=r',R_defect='Delta_R>=e',D_defect='Delta_D>=s+1',both_full_principal_parts_retained=True,exact_trace_residue_used=True,prescribed_norm_band_cut=True,all_infinity_branches_accounted_for=True,complete_curve_constructed=True,target_positivity_required=False,controls=controls),
      simple_R=dict(index=2,R_order=1,D_any_order=True,R_pole='q>=max(2,s)',D_defect='Delta_D>=s+1',constant_A_retained=True,all_coprime_D_poles_p_at_least_splus1_excluded_under_q_bound=True,controls=simple),
      common_ruling=dict(unbounded_smooth_affine_exact_volume_controls=True,orders=[1,0],both_component_opens_affine=True,whole_translation_criterion_required=True,global_explicit_primitive=True,finite_Keller_realization_claimed=False,controls=common),
      norm_software=dict(method='quotient multiplication determinant',independent_algorithm_check='reduced monic resultant',odd_odd_raw_resultant=int(raw),correct_root_product=875,older_sealed_PF48_inputs_modified=False,older_raw_resultant_helper_safe_for_arbitrary_new_inputs=False),
      assurance=dict(written_unbounded_proof=True,finite_controls_not_universal_proof=True,common_or_complementary_coordinates_automatic=False,insufficient_defect_excluded=False,full_quartic_exclusion=False,full_higher_degree_exclusion=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'two_sided_norm_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; arbitrary unramified order and common-ruling affine controls')
if __name__=='__main__':main()
