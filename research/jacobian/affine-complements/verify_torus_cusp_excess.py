"""PF114 exact arithmetic/finite-map controls; universal proof is the written argument."""
from pathlib import Path
from fractions import Fraction
from math import gcd
import hashlib,json,platform
import sympy as S
import torus_cusp_excess as T

P=Path(__file__).resolve().parent


def cycles(p):
    seen=set();out=[]
    for i in range(len(p)):
        if i in seen:continue
        c=[];j=i
        while j not in seen:
            seen.add(j);c.append(j);j=p[j]
        out.append(c)
    return out


def mul(p,q):return tuple(p[q[i]] for i in range(len(p)))


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    oldfiles=['local_monodromy_branch_receipt.json','triple_boundary_residue_receipt.json','cyclic_quotient_morse_receipt.json','coprime_leading_forms_receipt.json']
    for name in oldfiles:
        old=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/local-cusp-excess/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('complete primary Brin pages freshly rendered and read',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[26,27,86,87,145,146,147])
    check('classification/topology imports declared and no finite-sample theorem',manifest['Brin_topological_inputs_imported'] and manifest['Palka_classification_imported'] and not manifest['universal_claims_inferred_from_finite_examples'])
    a,b,w,z,alpha=S.symbols('a b w z alpha',positive=True)
    eta=((a*b-a-b)*w+a*b*z/alpha)/(a+b)
    ident('UNBOUNDED orbifold Euler excess identity',(w+eta)*(1-1/a-1/b)-(eta-z/alpha))
    A=5*a*b-6*a-6*b;C=5*a*b-3*a-8*b
    ident('UNBOUNDED first positive decomposition',A-((5*a-6)*(b-3)+9*(a-2)))
    ident('UNBOUNDED second positive decomposition',C-((5*a-8)*(b-3)+12*(a-2)))
    ident('UNBOUNDED complete positive excess certificate',(a+b)*(5*eta-w-2*z)-(A*(w-z)+C*z+5*a*(b/alpha-1)*z))
    x,y,t,u=S.symbols('x y t u',nonnegative=True)
    # a=2+x,b=3+y,w=z+t,b/alpha=1+u covers every allowed parameter.
    cert=S.expand((A*t+C*z).subs({a:2+x,b:3+y})+5*(2+x)*u*z)
    check('UNBOUNDED certificate coefficients all nonnegative',all(c>=0 for mon,c in S.Poly(cert,x,y,t,u,z).terms()))
    d,e,r,h=S.symbols('d e r h',integer=True)
    ident('UNBOUNDED index ledger with ceil excess',(5*(d-1-r)-2)- (5*d-7-5*r))
    ident('UNBOUNDED r>=2 universal numerator',(5*d-7-5*(2+h))-(5*d-17-5*h))
    F=S.symbols('F',integer=True,positive=True)
    ident('UNBOUNDED shared mass numerator including residues',(5*(d-2-r+h)-2*F)-(5*d-10-5*r+5*h-2*F))
    ident('UNBOUNDED minimal outside ramified prime numerator',(5*d-10-5*(h+1)+5*h-2*F)-(5*d-15-2*F))
    q=S.symbols('q',positive=True)
    ident('UNBOUNDED q divides delta strengthens local residue charge',(5*eta-w-2*q)-((5*eta-w-2*z)+2*(z-q)))
    j=S.symbols('j',integer=True,nonnegative=True)
    check('UNBOUNDED low even degrees d>=10 satisfy new cap',(5*(10+2*j)-17-6*(5+j)).is_nonnegative)
    check('UNBOUNDED low odd degrees d>=9 satisfy new cap',(5*(9+2*j)-17-6*(4+j)).is_nonnegative)
    ident('UNBOUNDED old d-4 strict high profile has eta1 contradiction',(d-4)-3-(d-7))
    s,m,k=S.symbols('s m k',integer=True,positive=True)
    # Genuine all-s A1 quotient family from sealed PF105, including q>1 and delta>1.
    w_s=(2*s+1)*s;ell_s=4*s*s;eta_s=ell_s-w_s
    ident('UNBOUNDED A1 family exact local Euler identity',ell_s*(1-S.Rational(1,2)-1/(4*s-1))-(eta_s-s/(4*s-1)))
    ident('UNBOUNDED A1 family universal bound slack',5*eta_s-w_s-2*s-8*s*(s-1))
    ident('UNBOUNDED source weighted action target V degree',S.Rational(1,2)/s*(4*s-1+1)-2)
    ident('UNBOUNDED source weighted action target U degree',S.Rational(1,2)/s*2*s*(4*s-1)-(4*s-1))
    # Genuine corank-one polynomial family: Q has constant leading v coefficient.
    ell_m=m*(k-1)+1;eta_m=ell_m-k
    ident('UNBOUNDED corank-one family exact local Euler identity',ell_m*(1-1/m-1/ell_m)-(eta_m-1/m))
    controls=[]
    certified=dict(cyclic_axis_pair_certified=True,sole_critical_branch_certified=True,torus_cusp_certified=True)
    for args,label in [((2,3,2,1,3,1,2),'genuine smooth cubic'),((2,3,3,1,4,1,3),'genuine singular A1 eta1'),((2,3,9,1,12,1,1),'genuine PF107 finite A1')]:
        g=T.local_gate(*args,**certified);check('allowed local control '+label,g['necessary_arithmetic_passes'] and not g['hypotheses_contradicted']);controls.append(dict(label=label,data=args))
    for ss in [1,2,7,101]:
        g=T.local_gate(2,4*ss-1,2*ss+1,ss,4*ss*ss,ss,4*ss-1,**certified)
        check('actual A1 family including longitude joining '+str(ss),g['necessary_arithmetic_passes'] and g['fibre_order_multiple_of_local_residue'] and g['base_degree']==4*ss and g['unramified_boundary_count']==2*ss-1)
        check('weighted lift order modulo actual quotient is s '+str(ss),Fraction(2,2*ss).denominator==ss)
    U,V,Z=S.symbols('U V Z')
    for ee,mm in [(2,2),(3,2),(5,3)]:
        Q=S.integrate((U-Z**mm)**(ee-1),(Z,0,V));jac=S.diff(Q,V)
        ident('actual finite smooth family Jacobian '+str((ee,mm)),jac-(U-V**mm)**(ee-1))
        poly=S.Poly(Q,V);rank=mm*(ee-1)+1
        check('actual finite smooth family monic after constant rescaling '+str((ee,mm)),poly.degree()==rank and not poly.LC().has(U) and poly.LC()!=0)
        branch=S.expand(Q.subs(U,V**mm));check('nonzero cusp coefficient and exact orders '+str((ee,mm)),S.Poly(branch,V).degree()==rank and gcd(mm,rank)==1)
        check('actual finite smooth family allowed by local gate '+str((ee,mm)),T.local_gate(mm,rank,ee,1,rank,1,mm,**certified)['necessary_arithmetic_passes'])
    # Explicit covers of the cusp exterior: the orders of the two cone loops
    # and the actual meridian give the predicted unfilled solid-torus base.
    for Aperm,Bperm,ee,expected_alpha in [((1,0,2),(1,2,0),2,2),((1,0,3,2),(1,2,0,3),3,3)]:
        n=len(Aperm);identity=tuple(range(n));mu=mul(Aperm,Bperm)
        check('exact torus group relation '+str(n),mul(Aperm,Aperm)==identity and mul(mul(Bperm,Bperm),Bperm)==identity)
        check('exact cusp meridian cycle type '+str(n),sorted(map(len,cycles(mu)))==[1,ee])
        cones=[2//len(c) for c in cycles(Aperm) if len(c)<2]+[3//len(c) for c in cycles(Bperm) if len(c)<3]
        check('exact source base sole cone and genus zero '+str(n),cones==[expected_alpha] and len(cycles(Aperm))+len(cycles(Bperm))+len(cycles(mu))==n+2)
    for ee in [4,5,9,1001]:
        g=T.local_gate(2,3,ee,1,ee+1,1,3,**certified)
        check('entire eta1 mass>=4 class excluded '+str(ee),not g['universal_excess_bound'] and g['hypotheses_contradicted'])
    uncertified=T.local_gate(2,3,5,1,6,1,3)
    check('arbitrary source or unverified cusp not excluded by failed arithmetic',not uncertified['source_hypotheses_certified'] and not uncertified['hypotheses_contradicted'])
    cone=T.local_gate(2,3,4,1,5,1,6,**certified)
    check('Euler equality alone does not allow inherited order6 cone',cone['exact_euler_balance'] and not cone['cone_order_allowed'])
    check('regular fibre cannot ignore a joined local residue orbit',not T.local_gate(2,3,3,2,8,1,3,**certified)['fibre_order_multiple_of_local_residue'])
    for dd in [9,10,14,31,1000,10**30]:
        cap=T.universal_index_cap(dd)
        check('universal proportional cap '+str(dd),cap==(5*dd-17)//6 and cap>=dd//2)
        check('old entire d-4 profile excluded '+str(dd),not T.high_index_gate(dd,dd-4)['necessary_budget_passes'])
        if 2*cap>dd:check('current high-index endpoint remains arithmetically possible '+str(dd),T.high_index_gate(dd,cap)['necessary_budget_passes'])
        check('above current high-index endpoint fails '+str(dd),not T.high_index_gate(dd,cap+1)['necessary_budget_passes'])
    check('small-degree bounds and strict-half boundary retained',T.universal_index_cap(7)==3 and T.universal_index_cap(8)==4 and T.high_index_gate(8,4)['necessary_budget_passes'])
    check('same-image proportional mass gate excludes old profile',not T.image_mass_gate(14,[(10,1)])['necessary_budget_passes'])
    check('same-image residues are counted and not deleted',T.image_mass_gate(30,[(7,2),(3,1)])['mass']==17 and not T.image_mass_gate(30,[(7,2),(3,1)])['unknown_primes_deleted'])
    check('outside ramified image really costs another prime',not T.image_mass_gate(30,[(7,2),(3,1)],2)['necessary_budget_passes'])
    text=' '.join((P/'TORUS_CUSP_SOLID_TORUS_EXCESS_AND_FIVE_SIXTHS_BOUND.md').read_text().lower().split())
    for phrase in ['whole marked u=a2','all omitted primes','one joint','ordinary unramified preimage curves','not additional omitted primes','intersection with the actual filling meridian m is one','nonidentity regular fibre','base degree ell/delta','a local link statement','classification is imported','not a cover realization certificate','e <= floor((5*d-17)/6) for d>=9','no scalar-morse','not replaced by a proposed census','remain open']:
        check('written scope '+phrase,phrase in text)
    files=['TORUS_CUSP_SOLID_TORUS_EXCESS_AND_FIVE_SIXTHS_BOUND.md','torus_cusp_excess.py','verify_torus_cusp_excess.py','sources/local-cusp-excess/manifest.json']+oldfiles+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED114',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_field_degrees=True,all_normal_actual_sources=True,
          original_finite_map_whole_plane_all_primes_required=True,
          cyclic_axis_torus_cusp_mass_at_most_five_excess_minus_two_fibre_order=True,
          nontrivial_local_residue_and_longitude_included=True,
          exact_orbifold_Euler_identity=True,actual_unramified_meridian_fillings_retained=True,
          universal_index_cap_floor_five_d_minus17_over6_for_d_at_least9=True,
          entire_prior_e_equals_d_minus4_strip_excluded=True,
          shared_image_all_residue_cusp_cost=True),
        allowed_genuine_finite_controls=controls,
        assurance=dict(Brin_topological_inputs_imported=True,Palka_classification_imported=True,
          written_unbounded_proof=True,finite_samples_prove_theorem=False,
          arbitrary_plane_branches_covered=False,scalar_Morse_assumed=False,
          local_residue_one_assumed=False,regular_fibre_identity_assumed=False,
          local_quotient_identifies_global_end=False,arithmetic_gate_realizes_cover=False,
          actual_wrong_plane_finite_models_excluded=False,source_symmetry_descent_assumed=False,
          global_finite_algebra_closed=False,formal_verification=False,external_review=False,
          research_goal_completed=False),environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'torus_cusp_excess_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF114 controls; all-degree torus-cusp excess and proportional index cap; full goal ACTIVE')


if __name__=='__main__':main()
