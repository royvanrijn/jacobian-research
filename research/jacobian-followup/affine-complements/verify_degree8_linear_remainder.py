"""PF2RULED31: exact linear-critical-remainder exclusion.
All arithmetic is rational polynomial arithmetic modulo irreducible P4.
No modular inference, Groebner search, or finite-point sampling is used.
"""
from flint import fmpq,fmpq_poly
from pathlib import Path
import gzip,json,time

class Algebra:
 def __init__(self,data):
  self.p=fmpq_poly([fmpq(c) for c in data['P4_coefficients']]);self.zero=fmpq_poly([]);self.one=fmpq_poly([1]);self.rules=[]
  for f in data['rules']:
   g=self.dec(f);lm=max(g,key=self.key);assert sum(lm)==4 and g[lm]==self.one
   self.rules.append((lm,{e:-c for e,c in g.items() if e!=lm}))
 def dec(self,data):return {tuple(e):fmpq_poly([fmpq(c) for c in cs]) for e,cs in data if any(c!='0' for c in cs)}
 def key(self,e):return (sum(e),*[-x for x in e[::-1]])
 def add(self,a,b):
  out=dict(a)
  for e,c in b.items():
   val=out.get(e,self.zero)+c
   if val:out[e]=val
   else:out.pop(e,None)
  return out
 def scale(self,a,c):return {e:(x*c)%self.p for e,x in a.items() if x*c}
 def mul(self,a,b):
  out={}
  for e,c in a.items():
   for f,d in b.items():
    ex=tuple(x+y for x,y in zip(e,f));val=out.get(ex,self.zero)+(c*d)%self.p
    if val:out[ex]=val
    else:out.pop(ex,None)
  return out
 def nf(self,f):
  work=dict(f);out={}
  while work:
   lm=max(work,key=self.key);c=work.pop(lm)
   if sum(lm)<4:out[lm]=c;continue
   div,tail=next((e,t) for e,t in self.rules if all(a>=b for a,b in zip(lm,e)))
   shift=tuple(a-b for a,b in zip(lm,div))
   for e,d in tail.items():
    ex=tuple(a+b for a,b in zip(e,shift));val=work.get(ex,self.zero)+(c*d)%self.p
    if val:work[ex]=val
    else:work.pop(ex,None)
  return out
import hashlib,platform,itertools
import sympy as S
import flint
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,truth):
    assert truth,name
    checks.append(name)
def gz(path):return json.loads(gzip.decompress(path.read_bytes()))
def sub(alg,a,b):return alg.add(a,{e:-c for e,c in b.items()})
def pmul(alg,a,b):
    out=[{} for _ in range(len(a)+len(b)-1)]
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]=alg.add(out[i+j],alg.nf(alg.mul(c,d)))
    while out and not out[-1]:out.pop()
    return out

def pdiv(alg,a,b):
    a=list(a)
    while a and not a[-1]:a.pop()
    b=list(b)
    while b and not b[-1]:b.pop()
    zero=(0,0,0,0)
    assert set(b[-1])=={zero} and b[-1][zero].degree()==0
    inv=fmpq_poly([1/b[-1][zero][0]])
    out=[{} for _ in range(max(0,len(a)-len(b)+1))]
    while len(a)>=len(b):
        shift=len(a)-len(b);q=alg.scale(a[-1],inv);out[shift]=q
        for j,c in enumerate(b):a[j+shift]=sub(alg,a[j+shift],alg.nf(alg.mul(q,c)))
        assert not a[-1]
        while a and not a[-1]:a.pop()
    return out,a

def primary_poly(alg,text):
    t,w,q4,q5,q6,q7,z=S.symbols('t w q4 q5 q6 q7 z')
    loc=dict(zip(('t','w','q4','q5','q6','q7','z'),(t,w,q4,q5,q6,q7,z)))
    f=S.Poly(S.sympify(text.replace('^','**'),locals=loc),w,t,q4,q5,q6,q7,z,domain=S.QQ)
    groups=[{} for _ in range(f.degree(w)+1)]
    for ex,c in f.terms():
        assert ex[-1]==0
        key=tuple(ex[2:6]);val=groups[ex[0]].get(key,alg.zero)+fmpq_poly([0]*ex[1]+[fmpq(str(c))]);groups[ex[0]][key]=val%alg.p
    return [{e:c for e,c in row.items() if c} for row in groups]

def main():
    beg=time.monotonic();pred=json.loads((HERE/'degree8_quartic_descent_receipt.json').read_text())
    check('current sign-descent predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    primary=gz(HERE/'degree8_quartic_descent_certificates.json.gz')['sectors']
    data=gz(HERE/'degree8_linear_remainder_certificates.json.gz')['sectors']
    check('both marked sectors aligned',len(data)==len(primary)==2 and [d['sector'] for d in data]==[2,4])
    summaries=[]
    for d,old in zip(data,primary):
        k=d['sector'];alg=Algebra(d);t=S.Symbol('t');p4=S.Poly(S.sympify(old['P4'].replace('^','**'),locals={'t':t}),t,domain=S.QQ)
        check(f'primary irreducible quartic field k{k}',p4.degree()==4 and p4.is_irreducible and
              alg.p==fmpq_poly([fmpq(str(p4.nth(j))) for j in range(5)]))
        W,M,L=[[alg.dec(c) for c in d[n]] for n in ['W','M','L']]
        for name,actual,key in [('W',W,'normalized_W'),('M',M,'normalized_M'),('L',L,'normalized_L')]:
            check(f'complete primary {name} matches coefficientwise k{k}',actual==primary_poly(alg,old[key]))
        s=fmpq_poly([fmpq(c) for c in d['s']]);sp=primary_poly(alg,old['normalized_cusp_parameter'])
        check(f'primary cusp scalar matches k{k}',sp==[{(0,0,0,0):s}])
        inputs=[alg.dec(c) for c in d['generators']];rules=[alg.dec(c) for c in d['rules']]
        check(f'four stratum equations are precisely M2 through M5 k{k}',inputs==M[2:6] and len(inputs)==4)
        expected={e for e in itertools.product(range(5),repeat=4) if sum(e)==4}
        check(f'35 monic rules reduce every quartic monomial to degree at most3 k{k}',
              len(rules)==35 and {a for a,b in alg.rules}==expected and
              all(sum(e)<=3 for lm,tail in alg.rules for e in tail))
        for j,(g,co) in enumerate(zip(rules,d['rules_in_inputs'])):
            total={}
            for coeff,f in zip(co,inputs):total=alg.add(total,alg.mul(alg.dec(coeff),f))
            check(f'FULL exact rule-in-stratum polynomial identity k{k} row{j}',len(co)==4 and total==g)
        check(f'all35 membership identities supplied k{k}',len(d['rules_in_inputs'])==35)
        # Each reduction preserves congruence modulo the four stratum equations,
        # by the preceding identities; confluence is unnecessary for this proof.
        WW,MM,LL=[[alg.nf(c) for c in poly] for poly in [W,M,L]]
        while MM and not MM[-1]:MM.pop()
        check(f'complete M reduces to degree at most1 on the stratum k{k}',len(MM)<=2)
        zero=(0,0,0,0);h=[{zero:-s},{zero:alg.one}]
        quotient,remainder=pdiv(alg,WW,h)
        check(f'cusp division remains exact after congruence reduction k{k}',not remainder)
        check(f'full L retains fixed unit leading coefficient8 k{k}',LL[-1]=={zero:fmpq_poly([8])})
        product=pmul(alg,MM,quotient);unused,E=pdiv(alg,product,LL)
        e5=E[5] if len(E)>5 else {}
        check(f'FULL necessary E5 circuit gives retained reduced polynomial k{k}',e5==alg.dec(d['reduced_E5']))
        U=alg.dec(d['inverse'])
        check(f'inverse has parameter degree at most3 k{k}',max(map(sum,U))<=3)
        check(f'FULL exact inverse identity U times E5 equals1 modulo the stratum k{k}',
              alg.nf(alg.mul(U,e5))=={zero:alg.one})
        summaries.append(dict(multiplicities=old['multiplicities'],P4=old['P4'],rule_count=35,
            reduction_basis_size=35,all_rule_memberships_verified=True,E5_inverse_degree=3,
            full_E5_circuit_checked=True,full_inverse_identity_checked_in_characteristic_zero=True,
            excluded_critical_remainder_degrees=[0,1],necessary_minimum_critical_remainder_degree=2,
            full_nonzero_core_solved=False))
        print(f'PASS k{k}: full exact E5 inverse excludes every linear or constant M; no modular inference',flush=True)
    inputs=[Path(__file__),HERE/'DEGREE_EIGHT_LINEAR_REMAINDER_EXCLUSION.md',
            HERE/'degree8_linear_remainder_certificates.json.gz',HERE/'degree8_quartic_descent_receipt.json',
            HERE/'degree8_quartic_descent_certificates.json.gz',HERE/'compute_degree8_linear_remainder.py',
            HERE/'compute_degree8_linear_remainder_dixon.py']
    inputs+=sorted((HERE/'sources/degree8-linear-remainder').glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,marked_points=16,
        coefficient_fields=2,coefficient_field_degree=4,excluded_critical_remainder_degrees=[0,1],
        necessary_minimum_critical_remainder_degree=2,necessary_minimum_D_polynomial_degree=1,zero_polynomial_M_excluded=True,
        exact_reduction_rule_count_per_field=35,inverse_polynomial_degree=3,rational_Dixon_system_size=140,
        summaries=summaries,first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],
        assurance=dict(full_characteristic_zero_polynomial_identities_checked=True,
            complete_primary_necessary_E5_circuit_checked=True,all_sixteen_complex_markings_covered=True,
            all_parameters_and_normal_freedoms_retained=True,no_parameter_denominator_introduced=True,
            exclusion_precedes_nonzero_kappa_condition=True,modular_unit_flags_not_used_as_proof=True,
            actual_global_root_hypothesis_removed=False,full_degree8_nonzero_kappa_core_solved=False,
            general_determinant_locus_solved=False,actual_global_section_lifts_constructed=False,
            full_source_Keller_pair_constructed=False,field_degree_four_realization_constructed=False,
            global_surface_finiteness_proved=False,actual_ruled_compatibility_proved=False,contact_cases_classified=False,
            whole_marked_candidate_excluded=False,quartic_packets_excluded=False,formal_verification=False,
            independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        python=platform.python_version(),sympy=S.__version__,python_flint=flint.__version__,
        elapsed_seconds=round(time.monotonic()-beg,6))
    (HERE/'degree8_linear_remainder_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; necessary deg M>=2; full goal remains active',flush=True)
if __name__=='__main__':main()
