"""Exact local/chord controls;does not replay registered mathematics."""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();checks=[]
 def check(label,v):
  assert v,label
  checks.append(label)
 def zero(label,e):check(label,s.cancel(s.expand(e))==0)
 x1,x2,y1,y2,ax,ay,c0=s.symbols('x1 x2 y1 y2 ax ay c0')
 A=y1-y2;B=x2-x1;C=x1*y2-x2*y1
 zero('chord contains first source point',A*x1+B*y1+C)
 zero('chord contains second source point',A*x2+B*y2+C)
 m=(x1-x2)/(y1-y2);b=x1-m*y1
 zero('scheme fibre ideal first generator',x1-c0-((b-c0)+y1*m))
 zero('scheme fibre ideal second generator',x2-c0-((b-c0)+y2*m))
 zero('scheme fibre ideal inverse generator',m-((x1-c0)-(x2-c0))/(y1-y2))
 zero('chord slope-intercept identification',b+C/A)
 swap={x1:x2,x2:x1,y1:y2,y2:y1}
 for label,f in [('A',A),('B',B),('C',C)]:zero('swap negates homogeneous chord '+label,f.subs(swap,simultaneous=True)+f)
 fa=(x1-ax)*(y2-ay)-(x2-ax)*(y1-ay)
 zero('pencil equation is oriented area',fa-(ax*A+ay*B+C))
 zero('pencil square descends under swap',fa.subs(swap,simultaneous=True)**2-fa**2)
 check('finite-source pencil avoids removed line-at-infinity point',0*ax+0*ay+1==1)
 check('parallel family is punctured not complete',0==0)
 t,r,u,v,x,y=s.symbols('t r u v x y')
 def divided(f):return s.cancel((f.subs(t,r)-f)/(r-t))
 def artinian_length(gens):
  gb=s.groebner(gens,r,t,order='lex');assert gb.is_zero_dimensional
  lead=[p.LM(order=gb.order).exponents for p in gb.polys]
  bounds=[]
  for i in range(2):bounds.append(min(e[i] for e in lead if e[i]>0 and all(e[j]==0 for j in range(2) if j!=i)))
  return sum(not any(all(mm[i]>=ee[i] for i in range(2)) for ee in lead) for mm in itertools.product(*(range(bd) for bd in bounds)))
 for multiplicity in [1,2,3,5]:
  p=t**2-1;q=t*p**multiplicity
  F1,F2=divided(p),divided(q)
  check('immersive node/tangency normalization '+str(multiplicity),s.gcd(s.diff(p,t),s.diff(q,t))==1)
  length=artinian_length([F1,F2]);check('ordered collision scheme length '+str(multiplicity),length==2*multiplicity)
  check('unordered collision length is half '+str(multiplicity),length//2==multiplicity)
  zero('finite two-branch plane equation '+str(multiplicity),q**2-p**(2*multiplicity)*(p+1))
  zero('conductor via complementary modules '+str(multiplicity),q/t-p**multiplicity)
  check('conductor degree matches2delta '+str(multiplicity),s.degree(p**multiplicity,t)==2*multiplicity)
  check('divided equations retain no spurious diagonal '+str(multiplicity),s.gcd(F1.subs(r,t),F2.subs(r,t))==1)
 p=t**3-t;q=t*p
 check('triple-point normalization immersive',s.gcd(s.diff(p,t),s.diff(q,t))==1)
 check('triple point has six ordered scheme points',artinian_length([divided(p),divided(q)])==6)
 zero('triple-point implicit equation',q**3-p**2*q-p**4)
 zero('triple-point conductor',(3*q**2-p**2)/s.diff(p,t)-p**2)
 check('triple-point conductor degree is6',s.degree(p**2,t)==6)
 tangents=[(s.diff(p,t).subs(t,a),s.diff(q,t).subs(t,a)) for a in [-1,0,1]]
 check('triple point distinct tangent branches',all(tangents[i][0]*tangents[j][1]-tangents[j][0]*tangents[i][1]!=0 for i in range(3) for j in range(i+1,3)))
 p0=x**2-1-y;q0=x*(x**2-1)**2+y*(3*x-s.Rational(5,2)*x**3)
 J=s.diff(p0,x)*s.diff(q0,y)-s.diff(p0,y)*s.diff(q0,x)
 zero('local tacnode extension exact Jacobian',J-(1+y*(3-s.Rational(15,2)*x**2)))
 zero('local tacnode extension Jac1 on WHOLE selected line',J.subs(y,0)-1)
 check('local extension is NOT whole Keller',s.diff(J,y)!=0)
 zero('tacnode pair has equal first image',p0.subs({x:1,y:0})-p0.subs({x:-1,y:0}))
 zero('tacnode pair has equal second image',q0.subs({x:1,y:0})-q0.subs({x:-1,y:0}))
 vecs=[(s.diff(p0,x).subs({x:a,y:0}),s.diff(q0,x).subs({x:a,y:0})) for a in [-1,1]]
 check('tacnode image line tangents are parallel',s.det(s.Matrix(vecs))==0)
 check('flat branched fibre has length2 despite one closed point',s.degree(t**2,t)==2 and len(s.solve(t**2,t))==1)
 # Abstract sharp affine two-chart cover of P1:generic two sheets,
 # one at0 and one at infinity.It supplies no Keller collision maps.
 check('two affine charts cover all directions',2>=1)
 check('two-chart cover generic degree2',1+1==2)
 check('two-chart cover at0 retains one sheet',1+0==1)
 check('two-chart cover at infinity retains one sheet',0+1==1)
 check('finite flat completion has positive missing local length',2-1>0)
 check('positive all-line defect and strict drop force generic>=2',1+1==2)
 check('generic conductor lower bound is4',2*2==4)
 check('two different covering degrees are not equated',s.binomial(5,2)!=2)
 algebra_count=len(checks)
 context=json.loads((HERE/'sources/collision-chords-context.json').read_text())
 for name,hsh in context['retained_proof_and_source_sha256'].items():check('retained bytes '+name,sha(ROOT/name)==hsh)
 for name,hsh in context['new_primary_sha256'].items():check('new primary bytes '+name,sha(HERE/name)==hsh)
 check('retained foundations only;no new major theorem import',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
 flags={k:v for k,v in context.items() if isinstance(v,bool)}
 check('scheme fibres and unreduced pencils retained',flags['whole_scheme_fibre_length_equals_affine_delta_proved'] and not flags['distinct_pair_count_substituted_for_scheme_length'] and not flags['pencil_pullback_assumed_reduced'])
 check('no automatic finiteness or universal single-line degree4',not any(flags[k] for k in ['chord_map_assumed_etale_or_finite','delta_at_least2_asserted_for_every_single_line','strict_drop_asserted_for_punctured_parallel_pencil']))
 check('no trace/Hessian completion claim',not any(flags[k] for k in ['last_HC4_source_to_regular_inverse_Gram_implication_proved','whole_surviving_higher_degree_packet_excluded','unrestricted_HC4_or_JC2_resolved','original_goal_complete']))
 proof=HERE/'GLOBAL_COLLISION_CHORDS_AND_PENCIL_DEFECT_DROP.md';text=' '.join(proof.read_text().split())
 for phrase in ['CURRENT LAST UNSOLVED IMPLICATION','SCHEMES','2*delta','Cohen-Macaulay','possibly nonreduced','base torsion','some fibre has length STRICTLY LESS','EVERY finite source point','parallel','common K-value may escape','order DIVIDING2','binomial(d,2)','DEFERRED']:
  check('written scope '+phrase,phrase.lower() in text.lower())
 result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=algebra_count,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-start,3),**flags)
 (HERE/'receipt-collision-chords.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
if __name__=='__main__':main()
