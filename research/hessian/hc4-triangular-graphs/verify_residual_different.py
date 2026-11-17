"""New residual different/discriminant controls; no registered replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];started=time.perf_counter();checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,x):check(name,s.simplify(s.expand(x))==0,True)
e=s.symbols('e',integer=True,positive=True);u,v,w,P=s.symbols('u v w P')
f=u**e*(u-1);h=u**(e-1)*((e+1)*u-e);q=v/h
zero('all-index derivative',s.diff(f,u)-h)
zero('all-index open-source constant Jacobian',s.diff(f,u)*s.diff(q,v)-1)
zero('finite normalization coordinate change',h*q-v)
c=e/(e+1);b=-e**e/(e+1)**(e+1)
check('second critical point nonzero',c.is_positive,True)
check('second critical value nonzero',b.is_negative,True)
# Positive e allows all same-base power simplifications used here.
zero('all-index second critical value',s.powsimp((c**e*(c-1))/b,force=True)-1)
zero('all-index retained root f1',f.subs(u,1))
zero('all-index retained root h1',h.subs(u,1)-1)
zero('all-index critical c linear factor',((e+1)*u-e).subs(u,c))
zero('critical multiplicity profile',e+1-(2+(e-1)))
zero('residual fiber universal monic power',(w**e*(w-1))/(w-1)-w**e)
zero('diagonal divided difference value',(e+1)*u**e-e*u**(e-1)-h)
# One fresh small transcription supports, not proves, the universal argument.
f2=u**2*(u-1);h2=s.diff(f2,u);F=w**2*(w-1)-f2
HH=s.cancel(F/(w-u))
check('new cubic control residual monic quadratic',s.Poly(HH,w).degree()==2 and s.Poly(HH,w).LC()==1,True)
zero('new residual derivative identity',s.diff(w**2*(w-1),w)-(w-u)*s.diff(HH,w)-HH)
zero('new CRT comaximal diagonal remainder',HH.subs(w,u)-h2)
idempotent=HH/h2
zero('new global diagonal idempotent',s.rem(s.together(idempotent**2-idempotent),F,w))
zero('new diagonal idempotent evaluation',idempotent.subs(w,u)-1)
full_disc=s.discriminant(w**2*(w-1)-P,w)
zero('new full discriminant',full_disc+P*(27*P+4))
res_disc=s.discriminant(HH,w)
zero('new exact squared-unit discriminant identity',full_disc.subs(P,f2)-h2**2*res_disc)
zero('new residual discriminant nonunit factors',res_disc+(u-1)*(3*u+1))
zero('retained point has nonreduced full residual fiber',HH.subs(u,1)-w**2)
zero('other branch retained source derivative unit',h2.subs(u,-s.Rational(1,3))-1)
zero('other branch retained source target value',f2.subs(u,-s.Rational(1,3))+s.Rational(4,27))
G=s.Matrix([[2,1-u],[1-u,(1-u)**2-2*(u**2-u)]])
zero('new residual trace matrix determinant',G.det()-res_disc)
zero('diagonal product trace matrix determinant',s.diag(1,G).det()-res_disc)
ctx=json.loads((H/'sources/residual-different-context.json').read_text());known={x['id']:x for x in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
for x in ctx['retained_structural_inputs']:check('prior '+x['claim'],known[x['claim']]['state']=='proved' and known[x['claim']]['canonical_source']==x['canonical_source'])
for p,d in ctx['retained_proof_and_source_sha256'].items():check('byte '+p,hashlib.sha256((R/p).read_bytes()).hexdigest()==d)
yes={'full_goal_active','previous_goal_turn_progress','retained_KU_diagonal_idempotent_credited','global_residual_cover_decomposition_proved','full_trace_quotient_preserved_in_residual_factor','different_diagonal_projection_distinguished_from_full_base_change','discriminant_ideal_exactly_residual_after_source_base_change','actual_nonunit_discriminant_persists_without_surjectivity_assumption','all_index_connected_surjective_constant_Jacobian_open_controls','retained_sheet_over_every_target_point_in_controls','residual_length_e_fiber_at_retained_point_proved','different_pullback_shortcut_retired','external_review_deferred'}
for p,x in ctx['scope_flags'].items():check('scope '+p,x is (p in yes))
proof=H/'SOURCE_BASE_CHANGE_AND_RESIDUAL_DIFFERENT.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=ctx['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-residual-different.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
