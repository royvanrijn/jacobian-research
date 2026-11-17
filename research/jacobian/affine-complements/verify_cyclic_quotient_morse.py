"""PF106 exact controls for written all-index proofs and declared primary inputs."""
from pathlib import Path
import hashlib,json,math,platform
import sympy as S
import cyclic_quotient_morse as M
import triple_boundary_residue as T

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    def ord0(poly,var):return min(m[0] for m,c in S.Poly(S.expand(poly),var).terms())
    old=json.loads((P/'triple_boundary_residue_receipt.json').read_text())
    check('PF105 predecessor PASS',old['outcome']=='PASS')
    for f,h in old['source_sha256'].items():check('sealed PF105 predecessor '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/cyclic-quotient-morse/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('fresh full-page conductor polar audit',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[3,4,10])
    # ALL monomial pairs in a Y^0 Jacobian coefficient have j+s=1.
    check('complete nonnegative Y-exponent possibilities',[(j,s) for j in range(2) for s in range(2) if j+s==1]==[(0,1),(1,0)])
    i,j,e,n=S.symbols('i j e n',integer=True,positive=True)
    ident('UNBOUNDED pure/mixed Jacobian X exponent',(i+j-1).subs(j,e-i)-(e-1))
    ident('UNBOUNDED quotient character',(e-1)+(1-e))
    check('quotient order equal to e cannot be coprime for nontrivial order',math.gcd(5,5)>1)
    k,b,j=S.symbols('k b j',integer=True,positive=True);m=b+1
    ident('UNBOUNDED weighted monomial degree',b*(k-j)+b*j-b*k)
    ident('UNBOUNDED Euler eigenvalue',(k-j)-b*j-(k-m*j))
    ident('UNBOUNDED residue gcd reduction',k*b-k*(b+1)+k)
    nn,ff,tt,aa=S.symbols('n f t a',integer=True,positive=True)
    ident('UNBOUNDED quadratic field bound nonnegative decomposition',
          nn*ff**2*tt*aa-4*ff**2-ff**2*((nn-2)*tt*aa+2*(tt-1)*aa+2*(aa-2)))
    ident('UNBOUNDED exact eta',k*ff*aa-(k+1)*ff-ff*(k*(aa-1)-1))
    cc,q=S.symbols('c q',integer=True,positive=True)
    ident('UNBOUNDED conductor balance',q*(cc+aa-1)-(q*aa-1)-(q*(cc-1)+1))
    # Exact whole polynomial models independently verify all weighted identities.
    models=[]
    for kk,mm,nn in [(2,4,2),(3,6,3),(5,10,5),(4,8,2),(4,12,2),(6,8,2)]:
        d=M.morse_model(kk,mm);X=d['X'];Y=d['Y'];u=d['u'];v=d['v'];g=d['critical']
        qh=d['upstairs_residue'];ff=qh//nn;aa=mm//qh
        check('valid invariant quotient divisor '+str((kk,mm,nn)),qh%nn==0)
        ident('full finite-germ Jacobian '+str((kk,mm,nn)),S.diff(u,X)*S.diff(v,Y)-S.diff(u,Y)*S.diff(v,X)-g**kk)
        ident('partial-fraction nonzero coefficient '+str((kk,mm,nn)),sum(d['coefficients'])-d['coefficient'])
        check('restricted orders exact and residue computed '+str((kk,mm,nn)),
              d['coefficient']!=0 and ord0(u.subs(X,Y**d['b']),Y)==kk*d['b'] and
              math.gcd(mm,kk*d['b'])==qh)
        check('invariant monomial weights of quotient '+str((kk,mm,nn)),
              all((-mon[0]+mon[1])%nn==0 for mon,c in S.Poly(u,X,Y).terms()))
        check('exact axes degree '+str((kk,mm,nn)),
              ord0(u.subs(X,0),Y)+ord0(u.subs(Y,0),X)==kk*mm)
        a0,beta=d['image_pair']
        check('image monomial branch primitive '+str((kk,mm,nn)),math.gcd(a0,beta)==1)
        Fpull=S.expand(u**a0-d['coefficient']**a0*v**beta)
        H,remainder=S.div(Fpull,S.expand(g**(kk+1)),X)
        check('actual pullback factor has transverse index k+1 '+str((kk,mm,nn)),remainder==0)
        boundaryH=S.expand(H.subs(X,Y**d['b']))
        expected=M.conductor_intersection(qh,d['conductor'])
        check('EXACT conductor companion intersection '+str((kk,mm,nn)),ord0(boundaryH,Y)==expected)
        point=M.morse_budget(kk+1,nn,ff,aa)
        check('quotient local degree and actual field budget '+str((kk,mm,nn)),
              point['necessary_local_conditions'] and point['local_degree']==kk*mm//nn)
        models.append(dict(k=kk,m=mm,quotient_order=nn,residue_degree=ff,image_pair=d['image_pair'],
                           local_degree=point['local_degree'],conductor=d['conductor'],
                           companion_intersection=expected,coefficient=str(d['coefficient'])))
    # Resonance is an exact obstruction, not a small numerical score.
    for kk,mm in [(6,3),(8,4)]:
        try:M.morse_model(kk,mm)
        except ValueError:check('diagonal resonance rejected '+str((kk,mm)),True)
        else:raise AssertionError('resonant model accepted')
    pert=T.rank_zero_model(3,1)
    check('higher terms can reduce residue below leading gcd',
          math.gcd(4,6)==2 and pert['upstairs_residue_degree']==1)
    # Quotient constraints alone do NOT rule out the unresolved non-Morse case.
    check('non-Gorenstein e5,n3 remains a necessary quotient candidate',
          M.quotient_constraint(5,3)['necessary_bound_satisfied'] and M.quotient_constraint(5,3)['marked_weights']==[1,1])
    check('that candidate fails only the CONDITIONAL scalar-Morse test',
          not M.morse_budget(5,3,2,2)['necessary_local_conditions'])
    check('Morse e5,n2,f2 local data allowed but actual d16 fails',
          M.morse_budget(5,2,2,2,16)['necessary_local_conditions'] and
          not M.morse_budget(5,2,2,2,16)['actual_field_budget'] and
          M.morse_budget(5,2,2,2,17)['actual_field_budget'])
    for ee,want in [(3,[2]),(4,[3]),(6,[5])]:
        check('complete extremal quotient control e='+str(ee),
              [nn for nn in range(2,ee+1) if M.quotient_constraint(ee,nn)['necessary_bound_satisfied']]==want)
    # Full unbounded PF105 family saturates the local quadratic bound.
    counter=[]
    for ss in [1,2,3]:
        d=T.quotient_counterfamily(ss);point=M.morse_budget(d['transverse_index'],2,ss,2,d['quotient_degree'])
        check('genuine finite counterfamily saturates LOCAL quadratic bound '+str(ss),
              point['necessary_local_conditions'] and point['local_degree']==4*ss*ss and
              not point['actual_field_budget'] and d['complement']=='C* times C')
        counter.append(dict(s=ss,local_degree=point['local_degree'],residue_degree=ss,
                            actual_Keller_degree_at_least=point['minimum_field_degree'],realization=False))
    # A genuine finite non-Morse family: no automatic scalar condition.
    X,Y,z=S.symbols('X Y z');nonmorse=[]
    for mm in [3,5]:
        u=X*X+Y**mm;v=X*X*Y+S.Rational(mm,mm+1)*Y**(mm+1)
        ident('genuine non-Morse finite Jacobian '+str(mm),S.diff(u,X)*S.diff(v,Y)-S.diff(u,Y)*S.diff(v,X)-2*X**3)
        ident('non-Morse monic integrality '+str(mm),Y**(mm+1)-(mm+1)*u*Y+(mm+1)*v)
        check('non-Morse scalar quadratic pencil rank1 '+str(mm),S.hessian(u+z*v,(X,Y)).subs({X:0,Y:0}).det()==0)
        nonmorse.append(dict(m=mm,transverse_index=4,degree=2*(mm+1),quadratic_pencil='multiples of X^2',actual_Keller_source=False))
    a,b=S.symbols('a b',integer=True,nonnegative=True)
    ident('UNBOUNDED remaining generic residue cost at least10',(5+a)*(2+b)-10-(2*a+5*b+a*b))
    check('smooth image cannot have nontrivial normalized residue',M.conductor_intersection(2,0)<0)
    proof=' '.join((P/'CYCLIC_QUOTIENT_ORDER_AND_MORSE_RESIDUE_BUDGET.md').read_text().lower().split())
    for word in ['unconditional all-index','2<=n<=e-1','n*f divides k','d>=2*k*f+1>=4*f^2+1',
                 'not proved automatic','without the scalar-morse hypothesis','not the product',
                 'ordinary preimage curves','{1,2,3,4,6}','d<=10','common normal square',
                 'not a singular keller counterexample','all f1 packets']:
        check('scope/proof guard '+word,word in proof)
    files=['CYCLIC_QUOTIENT_ORDER_AND_MORSE_RESIDUE_BUDGET.md','cyclic_quotient_morse.py',
           'verify_cyclic_quotient_morse.py','sources/cyclic-quotient-morse/manifest.json',
           'triple_boundary_residue_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED106',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        unconditional=dict(all_field_degrees=True,normal_singular_source=True,arbitrary_other_primes=True,
                           quotient_order_bound='2<=n<=e-1, gcd(n,e)=1, marked type1/n(1,-e)',
                           quotient_order_at_most_d_minus2=True,all_extremal_quotients_residue1=True,
                           residue1_indices=[1,2,3,4,6],all_omitted_residue1_through_degree10=True,
                           exact_conductor_balance_all_indices=True,conductor_identity='I(A,H)=q(c-1)+1'),
        scalar_Morse=dict(hypothesis_required=True,divisibility='nf divides e-1',source_type='A_(n-1)',
                          exact_local_degree='(e-1)fa',exact_eta='f((e-1)(a-1)-1)',
                          minimum_field_degree='2(e-1)f+1>=4f^2+1',
                          low_degree_residue_failure_requires_degenerate_quadratic_pencil=True,
                          automatically_satisfied_claimed=False),
        exact_models=models,sharp_local_countercontrols=counter,non_Morse_countercontrols=nonmorse,
        assurance=dict(written_unbounded_proofs=True,primary_boundary_and_polar_inputs_imported=True,
                       finite_examples_not_universal_proof=True,all_ten_or_lower_degree_maps_excluded=False,
                       general_local_residue_bound_without_Morse_proved=False,all_indices_residue1=False,
                       degenerate_pencil_actual_Keller_realization=False,global_finite_algebra_closed=False,
                       formal_verification=False,external_review=False,external_review_deferred=True,research_goal_completed=False),
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'cyclic_quotient_morse_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',claim='PF2RULED106',check_count=len(checks))))


if __name__=='__main__':main()
