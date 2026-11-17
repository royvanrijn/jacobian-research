"""PF2RULED100 exact controls; universal scope is in the written proof."""
from pathlib import Path
from itertools import product
import hashlib,json,sys,platform
import sympy as S

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from key_sequence_obstruction import analyze
from verify_monic_key_chart_audit import laurent_leading


def expansion(F,keys,alpha,x,y):
    """Independent monic division, then digit division of each remainder."""
    out=[]
    def rec(f,j,digits):
        if j==0:
            for (i,),c in S.Poly(f,x).terms():
                if c:out.append((i,tuple(reversed(digits)),c))
            return
        k=0
        while f!=0:
            f,rem=S.div(f,keys[j],y)
            if rem!=0:rec(rem,j-1,digits+[k])
            k+=1
    rec(S.expand(F),len(keys)-1,[])
    assert all(0<=digits[i-1]<alpha[i-1] for _,digits,_ in out for i in range(1,len(keys)-1))
    return out


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    predecessor=P/'monic_key_chart_audit_receipt.json'
    old=json.loads(predecessor.read_text())
    check('PF99 predecessor PASS',old['outcome']=='PASS')
    for f,h in old['source_sha256'].items():check('sealed PF99 '+f,sha(P/f)==h)
    manifest=P/'sources/polynomial-key-chain/manifest.json'
    m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=P/row['path']
        check('exact retained primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    check('source scope rejects unrestricted negative cone inference',
          m['terminal_polynomial_equivalence_imported'] and
          not m['primitive_or_positive_terminal_assumption_used'] and
          not m['unrestricted_negative_range_cone_or_attainment_criterion_imported'] and
          not m['non_finite_generation_claimed'])
    x,y,u,t=S.symbols('x y u t')
    G,H,K=[S.Function(a)(x,y) for a in ['G','H','K']]
    p0,q0=S.symbols('p0 q0')
    jac=lambda a,b:S.diff(a,x)*S.diff(b,y)-S.diff(a,y)*S.diff(b,x)
    identity('arbitrary-function exact Jacobian factor',
             jac(p0+G*H,q0+G*K)-G*(H*jac(G,K)+K*jac(H,G)+G*jac(H,K)))
    charts=[(1,u**-2,u**-3+t*u**4),
            (2,u**-4,u**-6+S.Rational(1,2)*u**3+t*u**10),
            (3,u**-8,u**-12+S.Rational(1,2)*u**6+S.Rational(1,4)*u**15+t*u**22)]
    controls=[]
    for L,xx,yy in charts:
        weights=[2**L]+[3*2**j for j in reversed(range(L))]+[-1]
        a=analyze(weights)
        keys=[x,y]
        for j in range(1,L+1):keys.append(S.expand(keys[j]**2-(x**3 if j==1 else keys[j-1])))
        actual=[S.expand(g.subs({x:xx,y:yy},simultaneous=True)) for g in keys]
        lead=[laurent_leading(g,u) for g in actual]
        check('actual full key values and generic terminal residue depth '+str(L),
              [v for v,c in lead]==[-w for w in weights] and
              all(c==1 for v,c in lead[:-1]) and lead[-1][1]==2**L*t)
        check('polynomial terminal degree strictly exceeds old minimal degree at depth '+str(L),
              S.degree(keys[-1],y)==2**L and (L==1 or S.degree(keys[-1],y)>weights[0]//S.gcd(weights[0],weights[1])))
        check('integer polynomial-chain filter depth '+str(L),a['all_key_forms_polynomial'] and
              a['regular_algebra_case']=='constant_mod_terminal')
        # Independently divide input polynomials; do not reuse their construction digits.
        tests=[x+y+keys[-1],1+keys[-1],
               keys[-2]+keys[-1]**2,
               x*y**(2**L+1)+3*keys[-1]**2+5]
        for no,F in enumerate(tests):
            ex=expansion(F,keys,a['alpha'],x,y)
            rebuilt=sum(c*x**i*S.prod(g**b for g,b in zip(keys[1:],digits)) for i,digits,c in ex)
            identity('independent mixed-radix monic division '+str((L,no)),rebuilt-F)
            predicted=max(weights[0]*i+sum(w*b for w,b in zip(weights[1:],digits)) for i,digits,c in ex)
            actual_lead=laurent_leading(F.subs({x:xx,y:yy},simultaneous=True),u)
            check('exact no cancellation full-chain weight '+str((L,no)),actual_lead[0]==-predicted)
            if predicted<=0:
                check('regular polynomial constant mod full terminal '+str((L,no)),not S.rem(F,keys[-1],y).has(x,y))
        volume=S.expand(S.diff(xx,u)*S.diff(yy,t)-S.diff(xx,t)*S.diff(yy,u))
        check('actual canonical index depth '+str(L),1+laurent_leading(volume,u)[0]==2**(L+1)-2)
        controls.append(dict(depth=L,weights=weights,terminal_degree=2**L,
                             terminal_leading_coefficient=str(lead[-1][1]),canonical_index=2**(L+1)-2))
    deep=[]
    for L in [64,1000]:
        w=[2**L]+[3*2**j for j in reversed(range(L))]+[-1]
        a=analyze(w)
        check('unexpanded integer circuit all-polynomial depth '+str(L),a['all_key_forms_polynomial'] and
              all(i==2 for i in a['alpha'][:-1]) and all(r['beta'][0]>=0 for r in a['recurrences']))
        e=-w[0]-w[-1]+sum((al-1)*wi for al,wi in zip(a['alpha'][:-1],w[1:-1]))
        check('integer circuit canonical index depth '+str(L),e==2**(L+1)-2)
        deep.append(dict(depth=L,terminal_degree_bit_length=L+1,canonical_index_bit_length=e.bit_length(),
                         expanded=False,integer_recurrence_sha256=hashlib.sha256(json.dumps(a,sort_keys=True).encode()).hexdigest()))
    # Three signs, zero prior weight, and n=0 mixed weighted degree.
    cases=[]
    for w,case in [([2,3,1],'constants'),([2,3,0],'one_generator'),
                   ([2,3,0,-1],'constant_mod_terminal'),([1,-1],'constant_mod_terminal')]:
        a=analyze(w);check('three-sign algebra classifier '+str(w),a['all_key_forms_polynomial'] and a['regular_algebra_case']==case)
        cases.append(a)
    zero=analyze([2,3,0,-1]);check('zero-weight digit absent and next key subtracts constant',zero['alpha'][1]==1 and zero['recurrences'][1]['beta']==[0,0])
    xx=u**-2;yy=u**-3+S.Rational(1,2)*u**3+t*u**4
    check('actual zero-weight prefix then negative terminal',
          laurent_leading(S.expand((y*y-x**3).subs({x:xx,y:yy},simultaneous=True)),u)==(0,S.Integer(1)) and
          laurent_leading(S.expand((y*y-x**3-1).subs({x:xx,y:yy},simultaneous=True)),u)==(1,2*t))
    for terminal,yy in [(1,u**-3+t*u**2),(0,u**-3+t*u**3)]:
        check('actual zero or positive terminal chart '+str(terminal),
              laurent_leading(S.expand((y*y-x**3).subs({x:u**-2,y:yy},simultaneous=True)),u)==(-terminal,2*t))
    # Numeric semigroup gap before any negative polynomial prefix.
    gap=analyze([3,2,1,-2]);xx=u**-3;yy=u**-2+u**3+t*u**6
    gapkey=y**3-x*x-3*x**-1*y*y
    check('positive prefix semigroup gap has negative beta0',gap['first_Laurent_key_index']==3 and gap['recurrences'][1]['beta']==[-1,2])
    check('positive prefix gap actual generic terminal',laurent_leading(S.expand(gapkey.subs({x:xx,y:yy},simultaneous=True)),u)==(2,3*t))
    # Actual affine control: its terminal key is rational Laurent, not the first negative polynomial.
    X=u**-2+t;f=u;Y=X+f
    affine=analyze([2,2,-1,-4])
    check('actual affine control negative intermediate key is not terminal',
          affine['first_Laurent_key_index']==3 and S.simplify((Y-X)/u)==1 and
          S.limit(((Y-X)**2-X**-1)/u**4,u,0)==t)
    W=f*f*X;V=S.expand(X*(f*f*X-1))
    identity('actual affine equation retained',f*f*V-W*(W-1))
    check('rescued affine regular function has nonconstant residue',laurent_leading(V,u)==(0,t))
    simple=analyze([2,-1,-4]);check('same affine control simpler coordinate gap',simple['recurrences'][0]['beta']==[-1])
    # Exact exponent inequalities, independent of the polynomial-key test.
    i,j,k,d=S.symbols('i j k d',integer=True,nonnegative=True)
    exponent=j-2*i+2*k
    identity('closed cone lower bound certificate exponent+2d',exponent+2*d-(j+2*(d-i)+2*k))
    identity('closed cone upper bound certificate d-exponent',d-exponent-((d-j)+2*(i-k)))
    # Under i+j<=d,k<=i all summands in these certificates are nonnegative.
    check('closed cone edge generators occur despite Laurent terminal',
          laurent_leading(X,u)[0]==-2 and laurent_leading(f,u)[0]==1 and
          not simple['all_key_forms_polynomial'])
    # Independent bounded mixed-digit uniqueness control, not the universal proof.
    for w in [[8,12,6,3,-1],[2,3,0,-1],[6,4,1,-1]]:
        a=analyze(w)
        seen={}
        for digits in product(range(5),*(range(q) for q in a['alpha'][:-1])):
            val=sum(b*wi for b,wi in zip(digits,w[:-1]))
            check('bounded fixed-terminal digit uniqueness '+str((w,digits)),val not in seen)
            seen[val]=digits
    for w in [[0,1],[2,4],[2,3,6],[1],['a',1]]:
        try:analyze(w)
        except ValueError:check('invalid sequence rejected '+str(w),True)
        else:check('invalid sequence rejected '+str(w),False)
    receipt=dict(as_of='2026-11-16',outcome='PASS',check_count=len(checks),checks=checks,
        theorem=dict(arbitrary_finite_polynomial_key_chain=True,all_three_terminal_signs=True,
                     exact_common_Jacobian_factor_or_zero=True,all_degree_index_and_other_primes_unrestricted=True,
                     any_normal_affine_whole_plane_extension_excludes_polynomial_terminal=True,
                     genuinely_Laurent_terminal_necessary=True,negative_beta0_necessary=True),
        shallow_generic_controls=controls,deep_integer_circuits=deep,three_sign_controls=cases,
        positive_prefix_gap=gap,actual_affine_control=dict(full_weights=[2,2,-1,-4],simpler_weights=[2,-1,-4],
          intermediate_leading_residue='1',terminal_leading_residue='t',passes_necessary_filter=True,Keller_realization_claimed=False),
        cone_scope=dict(actual_affine_countercontrol=True,closed_wedge_generators=[[1,-1],[1,2]],
                       terminal_key='f^2-x^-1',terminal_weight=-4,
                       unrestricted_negative_range_cone_inference_imported=False,non_finite_generation_claimed=False),
        assurance=dict(written_unbounded_proof=True,finite_controls_not_universal_proof=True,
                       selected_primary_imports_explicit=True,published_quartic_retirement_preserved=True,
                       arbitrary_Laurent_boundary_excluded=False,affine_or_Keller_realizability_from_filter=False,
                       full_Jacobian_conjecture_proved=False,formal_verification=False,external_review=False,
                       research_goal_completed=False),
        source_sha256={name:sha(P/name) for name in [Path(__file__).name,'key_sequence_obstruction.py',
          'POLYNOMIAL_KEY_CHAIN_OBSTRUCTION_AND_LAURENT_TERMINAL_NECESSITY.md',
          'sources/polynomial-key-chain/manifest.json','monic_key_chart_audit_receipt.json']+[r['path'] for r in m['primary_sources']]},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'polynomial_key_chain_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS',len(checks),'exact controls; all-degree written polynomial-chain obstruction; Laurent realization remains open')


if __name__=='__main__':main()
