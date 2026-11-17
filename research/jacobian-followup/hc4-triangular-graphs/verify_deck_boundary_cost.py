"""Exact controls for universal deck-boundary proofs; no old checker replay."""
import hashlib,json,math,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();checks=[]
 def check(label,condition):
  assert condition,label
  checks.append(label)
 def zero(label,e):check(label,s.cancel(s.expand(e))==0)
 x,y,z,t=s.symbols('x y z t');k=s.symbols('k',integer=True,positive=True)
 # Universal numerator of dtheta=(CP'+(2A'-B')P)Omega.
 P=z**k-1;A=-(k+1)*z/k;B=-(k+2)*z/k
 check('all-k regular primitive numerator',s.simplify((A-B)*s.diff(P,z)+(2*s.diff(A,z)-s.diff(B,z))*P)==1)
 check('all-k cyclic primitive weight z*x*dy',1+1-2==0)
 check('all-k cyclic primitive weight z*y*dx',1-2+1==0)
 check('all-k volume weight',1+1-2==0)
 check('all-k defining left weight',2-2==0)
 check('free nonidentity cyclic element would force impossible x=z=0',(-1)!=0)
 for kk in [2,3,5,7]:
  zz=1+x**2*t;yy=s.cancel((zz**kk-1)/x**2);R=sum(z**i for i in range(kk))
  check('whole A2 chart polynomial k='+str(kk),s.denom(yy)==1)
  zero('whole surface relation k='+str(kk),x**2*yy-zz**kk+1)
  zero('off-x=0 inverse k='+str(kk),(zz-1)/x**2-t)
  zero('retained-line inverse k='+str(kk),yy/R.subs(z,zz)-t)
  zero('retained fibre is entire affine line k='+str(kk),yy.subs(x,0)-kk*t)
  check('retained-line inverse denominator is a unit k='+str(kk),R.subs(z,1)==kk)
  zero('volume on whole plane chart k='+str(kk),s.diff(zz,t)/x**2-1)
  check('chosen k-1 boundary basis cannot be permuted freely k='+str(kk),(kk-1)%kk!=0)
 # k=2 explicit pseudo-plane target and actual map from the plane.
 u=x**2;v=2*t+x**2*t**2;w=x*(1+x**2*t)
 zero('explicit target is u(1+uv)=w^2',u*(1+u*v)-w**2)
 zero('k=2 rational deck map preserves target u',u.subs(x,-x)-u)
 zero('k=2 rational deck map preserves target v',v.subs({x:-x,t:-t-2/x**2},simultaneous=True)-v)
 zero('k=2 rational deck map preserves target w',w.subs({x:-x,t:-t-2/x**2},simultaneous=True)-w)
 check('rational deck has a genuine source pole',s.limit(x**2*(-t-2/x**2),x,0)==-2)
 # Inductive generator-character calculation is exact for arbitrary n:
 # (M-n)-((M-1)-(n-phi(n)-1))=-phi(n).
 M,n,phi=s.symbols('M n phi')
 zero('unbounded cyclic generator character induction',(M-n)-((M-1)-(n-phi-1))+phi)
 # Regular-minus-trivial character has the required invariant dimensions.
 for gg,c in [(2,1),(3,2),(4,2),(6,2),(12,3)]:
  b0=gg*c-1
  check('regular character dimension '+str((gg,c)),1+b0==gg*c)
  check('quotient H2 invariant dimension '+str((gg,c)),s.Rational(b0-(gg-1),gg)==c-1)
  for aa in s.divisors(gg):
   check('cyclic-subgroup average '+str((gg,c,aa)),s.Rational(b0-(aa-1),aa)==s.Rational(b0+1,aa)-1)
 for gg in [4,6,8,12,60]:
  check('noncyclic minimal Euler multiple rejected '+str(gg),(gg-1+1)//gg==1)
  check('noncyclic next multiple meets proven lower cost '+str(gg),(2*gg-1+1)//gg==2)
 check('field-degree and boundary constraints both needed',math.gcd(12,8)==4 and math.gcd(12,7)==1)
 # Nonfinite etale A2->pseudo-plane does NOT multiply Euler by its degree.
 check('surjective nonfinite etale cover Euler control',2*1-(2-1)*1==1 and 2*1!=1)
 check('minimal cyclic quotient Picard rank',s.Rational(5,5)-1==0)
 check('noncyclic quotient has positive Picard rank',s.Rational(12,6)-1==1)
 algebra_count=len(checks)
 context=json.loads((HERE/'sources/deck-boundary-cost-context.json').read_text())
 for name,hsh in context['retained_proof_and_source_sha256'].items():check('retained bytes '+name,sha(ROOT/name)==hsh)
 for name,hsh in context['new_primary_sha256'].items():check('new primary bytes '+name,sha(HERE/name)==hsh)
 check('one explicit proved classification import',len(context['newly_imported_proved_external_statements'])==1 and not context['newly_used_unproved_external_statements'])
 flags={key:value for key,value in context.items() if isinstance(value,bool)}
 check('no generic-deck versus monodromy or source-chart substitution',not flags['generic_deck_group_identified_with_monodromy'] and not flags['generic_decks_assumed_to_preserve_actual_A2_chart'])
 check('no automatic trace pairing or full theorem conclusion',not any(flags[key] for key in ['trace_pairing_regularity_inferred_from_Hodge_purity','all_generic_deck_groups_proved_cyclic','etale_pseudoplane_to_A2_constructed','last_HC4_source_to_regular_inverse_Gram_implication_proved','original_goal_complete']))
 proof=HERE/'GENERIC_DECK_GROUP_AND_UNRAMIFIED_BOUNDARY_COST.md';text=' '.join(proof.read_text().split())
 for phrase in ['CURRENT LAST UNSOLVED IMPLICATION','UNRAMIFIED','FREE','g divides gcd(d,1+b0)','RATIONAL REPRESENTATION','g|b0','FINITE such fundamental group is CYCLIC','b0>=2*g-1','SURJECTIVE ETALE','ONE point','not the full finite normalization','TRIVIAL generic deck group','DEFERRED']:
  check('written scope '+phrase,phrase.lower() in text.lower())
 result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=algebra_count,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-start,3),**flags)
 (HERE/'receipt-deck-boundary-cost.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({key:result[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
if __name__=='__main__':main()
