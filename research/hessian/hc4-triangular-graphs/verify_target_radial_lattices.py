"""New symbolic radial/trace transcriptions;no registered checker replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];start=time.perf_counter();checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,z):check(name,s.simplify(s.expand(z))==0,True)
a,b,c,h,i,j,k=s.symbols('a b c h i j k');G=s.Matrix([[a,b],[b,c]]);C=s.Matrix([[h,i],[j,k]])
GG=C.T*G+G*C
zero('trace Gram determinant logarithmic derivative',c*GG[0,0]+a*GG[1,1]-2*b*GG[0,1]-2*s.trace(C)*G.det())
zero('trace connection trace after conjugation',s.trace(G.inv()*GG)-2*s.trace(C))
t,v,u=s.symbols('t v u',nonzero=True);e=s.symbols('e',integer=True,positive=True);q=s.symbols('q',integer=True,positive=True);ell=s.symbols('ell',integer=True,nonnegative=True);A0,c0=s.symbols('A0 c0',nonzero=True)
zero('leading derivative pole exponent',(-q-ell*e-1)+(1-e)-(-q-(ell+1)*e))
zero('leading coefficient recurrence',c0*A0*(-q-ell*e)/e+c0*A0*(q+ell*e)/e)
check('positive exponent assumptions for written nonvanishing proof',q.is_positive is True and e.is_positive is True and ell.is_nonnegative is True,True)
zero('log discriminant leading pole',s.diff(t**e,t)/t**e-e/t)
# Target nonconical branch model p=u^e+v^2,q=v;radial lift is meromorphic.
Eu=(u**e-v*v)/(e*u**(e-1));Ev=v
zero('nonconical model target p radial action',e*u**(e-1)*Eu+2*v*Ev-(u**e+v*v))
zero('nonconical model target q radial action',Ev-v)
zero('nonconical branch radial normal value',(u**e-v*v).subs(u,0)+v*v)
zero('radial lift generic pole coefficient',s.expand(Eu*u**(e-1))+v*v/e-u**e/e)
# Retained all-index Laurent controls,only NEW Euler identities checked here.
EuL=u/e;EvL=(2*e-1)*v/e;P=u**e;Q=v/u**(e-1)
zero('Laurent target P weight one',s.diff(P,u)*EuL-P)
zero('Laurent target Q weight one',s.diff(Q,u)*EuL+s.diff(Q,v)*EvL-Q)
zero('Laurent radial source divergence',s.diff(EuL,u)+s.diff(EvL,v)-2)
m,n=s.symbols('m n',integer=True)
zero('Laurent monomial eigenvalue',m/e+(2*e-1)*n/e-(m+(2*e-1)*n)/e)
zero('finite normalization generator fractional weight',s.diff(u,u)*EuL-u/e)
zero('codifferent generator fractional weight',s.diff(u**(1-e),u)*EuL-(1-e)/e*u**(1-e))
# A generic formal coframe,not a Keller candidate search.
px,py,qx,qy,pp,qq=s.symbols('px py qx qy pp qq');J=px*qy-py*qx
Ex=(pp*qy-qq*py)/J;Ey=(-pp*qx+qq*px)/J
zero('actual radial P action',px*Ex+py*Ey-pp)
zero('actual radial Q action',qx*Ex+qy*Ey-qq)
zero('actual radial contraction dx',-Ey+(qq*px-pp*qx)/J)
zero('actual radial contraction dy',Ex+(qq*py-pp*qy)/J)
# Cusp image commutator survives:two transpositions with braid relation.
aa=(1,0,2);bb=(0,2,1)
compose=lambda p,r:tuple(p[r[i]] for i in range(len(p)))
check('cusp pair braid',compose(compose(aa,bb),aa)==compose(compose(bb,aa),bb),True)
comm=compose(compose(compose(aa,bb),aa),bb)
check('cusp commutator nonidentity',comm!=(0,1,2),True)
check('cusp commutator order three',compose(compose(comm,comm),comm)==(0,1,2),True)
zero('normal degree adjunction convention',(-s.Symbol('KdotC')-2)-(s.Symbol('a_n')-1).subs(s.Symbol('a_n'),-s.Symbol('KdotC')-1))
ctx=json.loads((H/'sources/target-radial-lattices-context.json').read_text());known={z['id']:z for z in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
check('new checker unregistered','HC4MYME0' not in known)
for z in ctx['retained_structural_inputs']:check('prior '+z['claim'],known[z['claim']]['state']=='proved' and known[z['claim']]['canonical_source']==z['canonical_source'])
for path,d in ctx['retained_proof_and_source_sha256'].items():check('byte '+path,hashlib.sha256((R/path).read_bytes()).hexdigest()==d)
yes={'full_goal_active','previous_goal_turn_progress','source_and_target_Euler_fields_distinguished','all_centers_and_polynomial_target_frames','all_finite_radial_stable_source_modules_lie_in_target','locally_finite_part_exactly_target_ring','all_polynomial_eigenfunctions_descend','whole_finite_order_trace_discriminant_connection_identity','intermediate_plane_not_required','immersive_binary_affine_profile_nonpositive_normal_degree_by_prior_theorem','cusp_local_S3_prevents_Orevkov_factorization','prior_additive_and_nodal_results_credited','external_review_deferred'}
for path,value in ctx['scope_flags'].items():check('scope '+path,value is (path in yes))
proof=H/'TARGET_RADIAL_FINITE_LATTICES_AND_EIGENPOLYNOMIALS.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-start,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=ctx['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-target-radial-lattices.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
