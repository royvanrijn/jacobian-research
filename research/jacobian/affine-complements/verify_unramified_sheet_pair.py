"""Controls for the universal marked-sheet and complete-curve proofs."""
from pathlib import Path
from math import gcd
import json,hashlib,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(expr,z):
    expr=S.expand(expr)
    if expr==0:return S.oo
    return min(int(t.as_powers_dict().get(z,0)) for t in S.Add.make_args(expr))
def coeff(expr,z,n):return S.expand(expr).coeff(z,n)
def norm_reference(r,h,q,z,x,y):
    return S.expand((-1)**((r+1)*q)*S.resultant(z**r-y,z**q*x-S.expand(z**q*h),z))
def defect(r,q,support):
    contacts=[2*q-max(j for j in support if a*j%r) for a in range(1,r)]
    ss=sum(contacts);return ss,r*q-ss
def main():
    preds=['lefschetz_centralizer_receipt.json','unramified_orevkov_receipt.json','universal_regular_end_receipt.json','higher_regular_coordinate_receipt.json','coordinate_regular_sector_receipt.json']
    for name in preds:
        d=json.loads((HERE/name).read_text());check('predecessor PASS '+name,d['outcome']=='PASS')
        for f,h in d['source_sha256'].items():check('predecessor seal '+name+' '+f,sha(HERE/f)==h)
    manifest=HERE/'sources/unramified-sheet-pair/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];check('primary bytes '+f.name,sha(f)==row['sha256'] and len(f.read_bytes())==row['bytes'])
    d,eta,delta=S.symbols('d eta delta')
    check('point stabilizer has the exact stronger degree inequality',S.expand(d*eta-(eta+2*delta)-((d-1)*eta-2*delta))==0)
    check('quartic unramified normal bundle bound is three times eta',S.expand(((d-1)*eta-2*delta).subs(d,4)-(3*eta-2*delta))==0)
    x,y,z,t,v,A=S.symbols('x y z t v A')
    simple=[]
    # Arbitrary coefficient controls; the written top-row/column proof is unbounded.
    for p,q in [(2,2),(2,3),(3,2),(3,3),(4,3)]:
        cc={j:S.Symbol('c'+str(j)) for j in range(2,q+1)};dd={j:S.Symbol('d'+str(j)) for j in range(2,p+1)}
        hR=sum(cc[j]*y**(-j) for j in cc);hD=sum(dd[j]*x**(-j) for j in dd)
        C=S.expand(y**q*hR);D=S.expand(x**p*hD)
        F=x**p*y**q-x**(p-1)*C-y**(q-1)*D-A*cc[q]*x**(p-2)
        FR=S.expand(F.subs(x,A+hR+y*t));MR=-q*(p-2)+1
        FD=S.expand(F.subs(y,hD+v));MD=-p*(q-2)
        check('all-coefficient simple R chart minimum p'+str(p)+' q'+str(q),val(FR,y)==MR)
        check('simple R leading chart equation is affine with exact unit coefficient p'+str(p)+' q'+str(q),S.Poly(coeff(FR,y,MR),t).degree()==1 and coeff(FR,y,MR).coeff(t)==cc[q]**(p-1))
        check('all-coefficient simple D chart minimum p'+str(p)+' q'+str(q),val(FD,x)==MD)
        check('simple D leading chart equation is affine with exact unit coefficient p'+str(p)+' q'+str(q),S.Poly(coeff(FD,x,MD),v).degree()==1 and coeff(FD,x,MD).coeff(v)==dd[p]**(q-1))
        check('simple polynomial has only the full lower rectangle and its top vertex p'+str(p)+' q'+str(q),all((i,j)==(p,q) or i<=p-1 and j<=q-1 for (i,j),c in S.Poly(F,x,y).terms()))
        simple.append(dict(p=p,q=q,R_order=MR,D_order=MD,all_coefficients_retained=True,simultaneous_constant_retained=True))
    Fbad=x*x*y*y-x-y
    check('discarding the simultaneous R constant fails the chart test',val(S.expand(Fbad.subs(x,1+y**-2+y*t)),y)==0)
    Ffixed=Fbad-1
    check('the retained constant correction restores the whole chart',val(S.expand(Ffixed.subs(x,1+y**-2+y*t)),y)==1)
    Fres=x*x*y*y-x-y*(1+x)
    check('nonzero unramified residue breaks the simple compact-curve construction',val(S.expand(Fres.subs(x,y**-2+y*t)),y)<1)
    controls=[]
    examples=[(2,2,{3:2,1:3},2),(2,3,{4:2,2:3},3),(2,4,{6:2,3:3},4),(2,4,{7:2,5:3,3:5},5),(2,6,{8:2,5:3},6),(1,2,{3:2},2),(3,3,{4:2,1:3},3),(3,4,{6:2,3:3},4),(4,5,{6:2,2:3},5)]
    for e,r,cs,p in examples:
        q=max(cs);h=sum(c*z**(-j) for j,c in cs.items());ss,Delta=defect(r,q,list(cs));tag='e'+str(e)+' r'+str(r)+' q'+str(q)+' p'+str(p)
        check('full-support data meet explicit sufficient hypotheses '+tag,gcd(r,*cs)==1 and r not in cs and min(cs)>=r-e+1 and r>=max(2,e) and p>=r and Delta>=e)
        P=norm_reference(r,h,q,z,x,y);P0=S.expand(x**(p-r)*P);AA=coeff(P0,y,q-1)
        DD=2+(x**(p-2) if p>2 else 0)
        F=S.expand(P0-y**(q-1)*(AA+DD))
        check('norm reference is a genuine polynomial '+tag,S.Poly(P,x,y) is not None and coeff(P,y,q)==x**r and P.subs(y,0)!=0)
        check('norm polynomial vanishes on the entire reference branch '+tag,S.expand(P.subs({x:h,y:z**r},simultaneous=True))==0)
        check('exact residue eliminates the top right-column coefficient '+tag,S.Poly(AA,x).degree()<=p-2)
        check('whole repaired top column preserves every D principal-part coefficient '+tag,coeff(F,y,q-1)==-DD)
        B=ss+2*q-p*q;M=B+e-r
        deriv=S.expand(S.diff(P0,x).subs({x:h,y:z**r},simultaneous=True))
        check('full conjugate derivative order equals contact sum '+tag,val(deriv,z)==B)
        err=S.expand((y**(q-1)*(AA+DD)).subs({x:h,y:z**r},simultaneous=True))
        check('exact characteristic defect controls the entire correction '+tag,val(err,z)>=M+Delta-e)
        FR=S.expand(F.subs({x:h+z**(e-r)*t,y:z**r},simultaneous=True));lead=coeff(FR,z,M)
        check('full R-root substitution is regular after the predicted factor '+tag,val(FR,z)==M)
        check('full R-root boundary equation has nonzero linear coefficient '+tag,S.Poly(lead,t).degree()==1 and lead.coeff(t)!=0)
        hD=DD/x**p;FD=S.expand(F.subs(y,hD+v));MD=-p*(q-2);ld=coeff(FD,x,MD)
        check('full D plane has its one finite implicit root '+tag,val(FD,x)==MD and S.Poly(ld,v).degree()==1 and ld.coeff(v)==2**(q-1))
        check('R infinity face survives the whole top-column repair '+tag,max(q*i-r*j for (i,j),c in S.Poly(F-P0,x,y).terms())<q*(p-r))
        check('D infinity face has exactly a degree-one jump '+tag,max(p*j-i for (i,j),c in S.Poly(F,x,y).terms())==p*(q-1) and 2*q-r>0)
        controls.append(dict(e=e,r=r,p=p,q=q,support=cs,S=ss,Delta=Delta,norm_reference=str(P),top_column=str(AA),curve=str(F),R_chart_order=M,D_chart_order=MD))
    # Trace is the algebraic residue identity, not a coordinate heuristic.
    for r in range(2,5):
        cr=S.Symbol('cr');q=r+1;h=2*z**(-q)+cr*z**(-r)+3*z**(-(r-1));P=norm_reference(r,h,q,z,x,y)
        check('general norm trace coefficient is exactly minus r*c_r r'+str(r),coeff(coeff(P,y,q-1),x,r-1)==-r*cr)
    r,q,p=4,12,4;h=2*z**-12+3*z**-3;ss,Delta=defect(r,q,[12,3]);P=norm_reference(r,h,q,z,x,y);DD=x*x+2;F=S.expand(P-y**(q-1)*(coeff(P,y,q-1)+DD));M=ss+2*q-p*q+2-r
    check('failed-defect exact control is outside the sufficient range',Delta<2)
    check('failed-defect control actually loses the R implicit chart',val(S.expand(F.subs({x:h+z**(2-r)*t,y:z**r},simultaneous=True)),z)<M)
    for r in range(2,10):
        rows=[q for q in range(max(2,r-1),r+6) if gcd(r,q)==1]
        check('all coprime numerical controls have defect exactly q r'+str(r),all(defect(r,q,[q])[1]==q for q in rows))
    inputs=[Path(__file__),HERE/'UNRAMIFIED_SHEET_AND_PAIRWISE_CURVE.md',manifest]+[HERE/r['path_relative_to_affine_complements'] for r in m['primary_sources']]+[HERE/n for n in preds]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        marked_sheet=dict(all_field_degrees=True,full_finite_normalization_required=True,unramified_omitted_prime=True,residue_one_immersive_image=True,shared_ramified_image_allowed=True,other_local_groups_may_be_nonabelian=True,tube_image='L subset G_a',index='[G:G_a]=d',positive_normal_degree_required=True,bound='(d-1)*eta_H<=2*delta_H',quartic_bound='3*eta_H<=2*delta_H',retained_unramified_coordinate_line_corollary=True,eta_nonpositive_excluded=False),
        simple_pair=dict(all_pole_orders=True,all_lower_coefficients_retained=True,simultaneous_constant_A_retained=True,source_coordinates_complementary_required=True,R_regular_coordinate_simple=True,D_regular_coordinate_simple=True,indices=[2,1],exact_residues_required=True,whole_sector_excluded=True,other_primes_allowed=True,finite_map_not_required_for_curve_obstruction=True,polynomial='x^p*y^q-x^(p-1)*C(y)-y^(q-1)*D(x)-A*c_q*x^(p-2)',complete_curve_constructed=True,controls=simple),
        higher_pair=dict(all_generic_indices=True,R_component_open_smooth_required=True,regular_coordinate_exists_required=True,source_coordinates_complementary_required=True,D_regular_coordinate_simple_required=True,R_order='r>=max(2,e)',D_pole='p>=r',support='r-e+1<=j<=q',residues='c_r=d_1=0',contact_sum='S=sum_a(2*q-j(a))',sufficient_defect='r*q-S>=e',full_norm_resultant=True,entire_top_column_repaired=True,polynomial='x^(p-r)*P-y^(q-1)*(A0+D)',complete_curve_constructed=True,index_two_all_coprime_q_r_with_q_at_least_two_and_p_at_least_r_excluded=True,index_two_borderline_q_rminus1_with_r_at_least_three_and_p_at_least_r_excluded=True,p_less_than_r_excluded=False,index_two_defect_at_most_one_excluded=False,r2_q1_uniformly_excluded=False,controls=controls,failed_defect_control=dict(r=4,q=12,p=4,S=ss,Delta=Delta,affine_realization_claimed=False)),
        assurance=dict(written_unbounded_marked_sheet_proof=True,written_unbounded_complete_curve_proof=True,all_infinity_branches_accounted_for=True,actual_source_chart_embeddings_retained=True,no_deep_coefficient_search=True,finite_controls_not_universal_proof=True,common_complementary_coordinates_automatic=False,deleting_other_primes_preserves_finiteness=False,all_D_vanishing_orders_excluded=False,full_quartic_exclusion=False,higher_degree_exclusion=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'unramified_sheet_pair_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(checks),'controls; all-degree marked sheet bound and complementary two-prime complete-curve exclusion')
if __name__=='__main__':main()
