"""PF107 controls support a written universal proof, never a finite sweep proof."""
from pathlib import Path
import hashlib,json,math,platform
import sympy as S
import coprime_leading_forms as M
import triple_boundary_residue as T

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    old=json.loads((P/'cyclic_quotient_morse_receipt.json').read_text())
    check('PF106 predecessor PASS',old['outcome']=='PASS')
    for f,h in old['source_sha256'].items():check('PF106 seal '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/coprime-leading-forms/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('source/context exact pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('full selected passport/context pages audited',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[1,2,3])
    check('bijective-curve germ classification not imported',manifest['germ_classification_with_bijectivity_imported'] is False)
    a,b,g,n,f,t=S.symbols('a b g n f t',integer=True,positive=True)
    r,s=g*a,g*b;e=r+s-1;N=r*s/g
    ident('UNBOUNDED projective degree',N-g*a*b)
    ident('UNBOUNDED e<=N rewritten',N-e-(g*(a*b-a-b)+1))
    ident('UNBOUNDED a1 forces g<=1',(e-N).subs(a,1)-(g-1))
    ident('UNBOUNDED e+1 divisible by residue',(e+1).subs(g,n*f*t)-n*f*t*(a+b))
    ident('UNBOUNDED actual downstairs localdegree',(r*s/n).subs(g,n*f*t)-n*f*f*t*t*a*b)
    ident('UNBOUNDED transverse lower bound gap',(e-(5*n*f-1)).subs(g,n*f*t)-n*f*(t*(a+b-5)+5*(t-1)))
    ident('UNBOUNDED ambient lower bound gap',n*f*f*t*t*a*b-6*n*f*f-n*f*f*(t*t*(a*b-6)+6*(t*t-1)))
    check('least distinct coprime normalized orders',math.gcd(2,3)==1 and 2+3==5 and 2*3==6)
    ident('UNBOUNDED q-family constant Wronskian degree cancellation',3*(2*g)-2*(3*g))
    ident('UNBOUNDED q-family passport counts',4*g+3*g+(5*g-2)-(12*g-2))
    ident('UNBOUNDED q-family excess',6*g*g/n-(5*g-1)*g/n-g*(g+1)/n)
    models=[]
    for q in [1,2]:
        d=M.explicit_pair(q);X,Y,t0=d['X'],d['Y'],d['t'];u,v=d['P'],d['Q']
        jac=S.diff(u,X)*S.diff(v,Y)-S.diff(u,Y)*S.diff(v,X)
        c=3 if q==1 else -S.Rational(9,4)
        ident('exact homogeneous Jacobian q'+str(q),jac-c*X**(5*q-2))
        check('no common projective zeros q'+str(q),S.resultant(d['p'],d['z'],t0)!=0 and
              S.degree(d['p'],t0)==2*q and S.degree(d['z'],t0)==3*q)
        check('exact critical restriction residue q'+str(q),math.gcd(S.degree(u.subs(X,0),Y),S.degree(v.subs(X,0),Y))==q)
        check('finite homogeneous degree q'+str(q),S.total_degree(u)*S.total_degree(v)==6*q*q)
        diff=S.expand(d['p']**3-d['z']**2)
        check('DZ exact difference degree q'+str(q),S.degree(diff,t0)==q+1)
        W=S.expand(3*S.diff(d['p'],t0)*d['z']-2*d['p']*S.diff(d['z'],t0))
        check('nonzero CONSTANT projective Wronskian q'+str(q),W!=0 and not W.has(t0))
        H,rem=S.div(S.expand(u**3-v**2),X**(5*q-1),X)
        check('single image branch pullback factor q'+str(q),rem==0 and H.subs(X,0)!=0)
        check('conductor companion intersection q'+str(q),S.degree(H.subs(X,0),Y)==q+1)
        if q==2:
            check('both complete quadratic parts zero',min(S.total_degree(u),S.total_degree(v))>=4)
            check('actual scalar order2 invariance',all((i+j)%2==0 for poly in [u,v] for (i,j),c0 in S.Poly(poly,X,Y).terms()))
            packet=M.leading_budget(4,6,2,1,12)
            check('genuine degree12 quotient fails whole-plane d13 budget',packet['necessary_local_conditions'] and packet['local_degree']==12 and not packet['actual_field_budget'])
            check('degree13 budget is necessary not realization',M.leading_budget(4,6,2,1,13)['actual_field_budget'] and not packet['realization_claimed'])
        models.append(dict(q=q,Jacobian=str(jac),degree=6*q*q,exact_residue=q,actual_Keller_source=False))
    controls=[]
    # Numerical controls of a theorem-based family; no coefficient search.
    for qq,nn in [(2,2),(3,3),(4,2),(6,3),(6,6),(10,5)]:
        d=M.cusp_family(qq,nn);budget=M.leading_budget(2*qq,3*qq,nn,qq//nn,d['degree'])
        check('passport sum/gcd/localdegree '+str((qq,nn)),3*(2*qq)==2*(3*qq) and math.gcd(3,2)==1 and 5*qq-1<=6*qq)
        check('sharp conditional projective budget '+str((qq,nn)),budget['necessary_local_conditions'] and budget['local_degree']==d['degree'] and not budget['actual_field_budget'])
        check('wrong complement retained '+str((qq,nn)),d['complement']=='C* times C' and not d['actual_Keller_source'])
        check('exact downstairs excess '+str((qq,nn)),d['degree']-d['transverse_index']*d['residue_degree']==d['finite_excess'])
        if nn>2:check('non-Gorenstein actually permitted '+str((qq,nn)),d['non_Gorenstein'] and 2%nn!=0)
        controls.append(d)
    check('non-Gorenstein q6n3 contradicts automatic Morse divisor',M.cusp_family(6,3)['transverse_index']==29 and 28%6!=0)
    check('nontrivial residue first conditional d49',M.leading_budget(8,12,2,2,49)['actual_field_budget'] and not M.leading_budget(8,12,2,2,48)['actual_field_budget'])
    check('e5n3 coprime-leading orders impossible only in this sector',not M.leading_budget(3,3,3,1)['necessary_local_conditions'])
    try:M.cusp_family(5,3)
    except ValueError:check('invalid quotient rejected',True)
    else:raise AssertionError('invalid quotient accepted')
    d=T.rank_zero_model(3,0);X,Y=d['X'],d['Y']
    initial_u=X**2/2;initial_v=X*Y
    check('genuine Morse countercontrol has common tangent factor',S.gcd(initial_u,initial_v)==X)
    check('common tangent invalidates naive degree product',d['local_surface_degree']!=4)
    proof=' '.join((P/'COPRIME_LEADING_FORMS_AND_DEGENERATE_RESIDUE_BUDGET.md').read_text().lower().split())
    for word in ['no common zero','q divides g','n*f divides e+1','d>=6*n*f^2+1>=12*f^2+1',
                 'd>=49','common projective zero','not keller completions','non-gorenstein',
                 'apparent value-label inconsistency','general inequality','still unproved','all f1 packets']:
        check('scope guard '+word,word in proof)
    files=['COPRIME_LEADING_FORMS_AND_DEGENERATE_RESIDUE_BUDGET.md','coprime_leading_forms.py',
           'verify_coprime_leading_forms.py','sources/coprime-leading-forms/manifest.json',
           'cyclic_quotient_morse_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED107',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={x:sha(P/x) for x in dict.fromkeys(files)},
        theorem=dict(all_degrees=True,all_indices=True,normal_quotient_source=True,
                     coprime_leading_forms_hypothesis_required=True,automatic_hypothesis=False,
                     residue_divides_e_plus1=True,scalar_marked_quotient=True,
                     bound='d>=6*n*f^2+1>=12*f^2+1',nontrivial_residue_minimum_degree=49),
        explicit_polynomial_models=models,passport_family_controls=controls,
        assurance=dict(written_unbounded_proof=True,passport_existence_imported=True,
                       imported_germ_bijectivity_classification=False,finite_sweep_proves_theorem=False,
                       all_degenerate_pencils_classified=False,general_local_residue_bound_proved=False,
                       all_indices_residue1=False,actual_Keller_controls=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'coprime_leading_forms_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF107 controls; full goal remains active')


if __name__=='__main__':main()
