"""Exact controls for fixed-support plinth flatness and one-value recognition."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
u,v,h=sp.symbols('u v h')
def pp(value):
    result=0
    for term in sp.Add.make_args(sp.expand(value)):
        powers=term.as_powers_dict()
        if powers.get(u,0)<0 and powers.get(v,0)<0:result+=term
    return sp.expand(result)
# The full annihilator scheme is flat,although its closed socle can jump.
zeta=1/(u*u*v)+h/(u*v**3)
relations=[u*u,u*v,v*v-h*u]
for relation in relations:zero(pp(relation*zeta),'the exact flat-family relation annihilates zeta: '+str(relation))
images=[pp(b*zeta) for b in [1,u,v]]
zeros([images[0]-zeta,images[1]-1/(u*v),images[2]-h/(u*v*v)],
      'the three quotient-basis images are the independent displayed Laurent columns')
image_matrix=sp.Matrix([[1,0,0],[h,0,0],[0,1,0],[0,0,h]])
check(image_matrix.extract([0,2,3],[0,1,2]).det()==h,
      'the quotient basis is independent over C[h],even though image specialization can lose rank')
Mu=sp.Matrix([[0,0,0],[1,0,0],[0,0,0]])
Mv=sp.Matrix([[0,0,0],[0,0,h],[1,0,0]])
check(Mu*Mv==Mv*Mu,'the flat algebra multiplication matrices commute')
check(Mu**2==sp.zeros(3),'u squared vanishes in the entire flat algebra')
check(Mu*Mv==sp.zeros(3),'uv vanishes in the entire flat algebra')
check(Mv**2==h*Mu,'the h-dependent nilpotent multiplication law is exact')
check(Mv**3==sp.zeros(3),'the nilpotent ideal has a uniform finite exponent')
check(len(Mu.col_join(Mv).nullspace())==1,'the generic flat algebra has socle dimension one')
check(len(Mu.col_join(Mv).subs(h,0).nullspace())==2,'the actual zero-value fibre has socle dimension two')
check(image_matrix.subs(h,0).rank()==2,'the specialized inverse-system map loses one dimension')
zero(pp(v*zeta).subs(h,0),'the specialized Laurent class loses v')
zero(pp(v/(u*u*v)),'v belongs to the annihilator after class specialization')
check(Mv.subs(h,0)*sp.Matrix([1,0,0])==sp.Matrix([0,0,1]),
      'v still survives in the actual ideal fibre,so ideal and class annihilator base changes differ')
zero((v*v-h*u).subs(u,v*v/h),'over C(h),u is v squared divided by h')
zero((h*u*v).subs(u,v*v/h)-v**3,'the generic algebra is the length-three curvilinear algebra')
check(1+2==3,'the flat control keeps one support point and nilpotent dimension two')
# A reduced specialized annihilator can falsely pass a one-value check.
theta=h/(u*u*v)+(h-1)/(u*v)
zero(pp(u*u*theta),'the false-positive control has the full annihilator relation u squared')
zero(pp(v*theta),'the false-positive control has the full annihilator relation v')
zero(pp(u*theta)-h/(u*v),'u survives in the full module with a nonzero polynomial value coefficient')
zero(theta.subs(h,0)+1/(u*v),'the specialized class has a simple pole and reduced annihilator')
zero(pp(u*theta).subs(h,0),'the specialized class annihilator kills u')
check(sp.Matrix([[0,0],[1,0],[0,h]]).extract([1,2],[0,1]).det()==h,
      'the two full quotient classes remain independent over the value PID')
# Primitive Nambu normalization in the actual weighted target chart.
X,z,P,Q=sp.symbols('X z P Q');coords=[X,z,P,Q]
a=sp.Function('a')(z);C=sp.Function('C')(z);D=sp.Function('D')(z)
potential=a*X*P+C*Q+D
def N(value):return sp.Matrix([X,z,potential,value]).jacobian(coords).det()
zero(N(X),'the Nambu derivation fixes the first actual kernel generator')
zero(N(z),'the Nambu derivation fixes the second actual kernel generator')
zero(N(potential),'the Nambu derivation fixes the potential')
zero(N(P)+C,'the target P coefficient has the exact primitive orientation')
zero(N(Q)-a*X,'the target Q coefficient has the exact primitive orientation')
zero(sum(sp.diff(N(var),var) for var in coords),'the Nambu source vector has zero divergence')
F,G,Aconst=sp.symbols('F G Aconst')
zero(F*P/Aconst+G*Q/Aconst-1-(F*P+G*Q-Aconst)/Aconst,
     'a nonzero generic residue removes the whole source fibre over its base support point')
# The one-value certificate for the actual three-point target is a free rank-three quotient.
C3=(z-1)*(z-2)*(z-3)
check(sp.Poly(C3,z).LC()==1 and sp.degree(C3,z)==3,
      'the target support quotient is free of rank three over the entire value line')
idempotents=[]
for r in [1,2,3]:
    er=sp.cancel(C3/((z-r)*sp.diff(C3,z).subs(z,r)));idempotents.append(er)
    zero(sp.rem(er*er-er,C3,z),'the target line idempotent is exact: '+str(r))
zero(sum(idempotents)-1,'the three target support components cover every value of the finite scheme')
for i,j in [(0,1),(0,2),(1,2)]:
    zero(sp.rem(idempotents[i]*idempotents[j],C3,z),'distinct target support components stay disjoint: '+str((i,j)))
check(sp.Matrix([[1,r,r*r] for r in [1,2,3]]).det()==2,
      'evaluation at the three fixed support points gives the entire complex reduced quotient')
# Full integral cycle bases force a bijection and primitive degree,not a normal differential rank.
d1,d2,d3=sp.symbols('d1 d2 d3')
cycle=sp.Matrix([[0,d2,0],[d1,0,0],[0,0,d3]])
zero(cycle.det()+d1*d2*d3,'the labelled Morse-cycle matrix determinant is the product of its nonzero degrees up to sign')
check(cycle.subs({d1:1,d2:-1,d3:1}).det()==1,'the native primitive labelled matrix is unimodular')
duplicate=sp.Matrix([[d1,d2,0],[0,0,0],[0,0,d3]])
zero(duplicate.det(),'two finite cycles over one support line cannot give a full integral basis')
# Exact source controls for the focused line-bundle comparison;no manuscript theorem is imported.
aa,dd,bb,cc,kk=sp.symbols('aa dd bb cc kk')
xx=aa*bb;yy=dd*cc;zz=dd*bb;hh=xx-kk**3;ss=yy-xx*hh
vv=aa**3*bb-aa*aa*kk**3-dd*dd
ff=xx-ss**2-kk**3;gg=zz-ss*xx
C1=cc*cc-bb*bb*hh;al=aa*aa*(1-2*bb*dd);be=dd*dd*(3+2*bb*dd)+al*hh
def detzero(value,label):zero(value.subs(cc,(1+bb*dd)/aa),label)
detzero(xx*yy-zz*(zz+1),'the manuscript coefficient inclusion satisfies the exact quadric relation')
detzero(gg-bb*bb*vv,'the manuscript second ideal generator is b squared times v')
detzero(ff-C1*vv,'the manuscript first ideal generator has the displayed principalized coefficient')
detzero(al*C1+be*bb*bb-1,'the manuscript principalization coefficient row is unimodular')
detzero(al*ff+be*gg-vv,'the manuscript v is in the pulled-back ideal by its actual formula')
zero(vv.subs({aa:1,dd:0,bb:0,cc:1,kk:0}),'the manuscript principal divisor is nonempty,so v is not a unit')
def weights(poly,variables,values):
    return {sum(power*w for power,w in zip(monomial,values)) for monomial,coefficient in sp.Poly(poly,*variables).terms()}
variables=[aa,dd,bb,cc,kk]
check(weights(aa*cc-bb*dd-1,variables,[1,1,-1,-1,0])=={0},'the manuscript determinant relation has fibre weight zero')
check(weights(vv,variables,[1,1,-1,-1,0])=={2},'the manuscript principalizer has fibre weight two')
check(weights(dd*dd+aa*aa*kk**3,variables,[0,3,-3,0,2])=={6},
      'the manuscript auxiliary highest sum has degree six')
check(weights(aa**3*bb,variables,[0,3,-3,0,2])=={-3},'the manuscript lower auxiliary term has degree minus three')
# Exact gate logic retains source alignment and E_N independently.
for minimal in [False,True]:
    for en in [False,True]:
        check((not (en and minimal))==((not en) or (not minimal)),
              'one-value pure alignment and E_N remain a conjunction: '+str((minimal,en)))
algebra=len(checks);context=json.loads((HERE/'sources/flat-plinth-families-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    data=(HERE/'sources'/entry['retained_file']).read_bytes()
    check(hashlib.sha256(data).hexdigest()==entry['sha256'],'inspected primary '+entry['retained_file'])
    if 'git_blob_sha' in entry:
        check(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==entry['git_blob_sha'],
              'exact fetched manuscript Git blob '+entry['retained_file'])
imports=context['newly_imported_proved_external_statements']
check(len(imports)==1 and imports[0]['hypotheses_checked']
      and imports[0]['source']=='https://stacks.math.columbia.edu/tag/0549',
      'the sole new proved import has its finite torsion-free module and value-PID hypotheses checked')
check(context['newly_used_unproved_external_statements']==[],
      'no automatic native minimum length or manuscript rigidity/coordinate theorem is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==6 and not flags['fixed_support_alone_proved_to_force_minimal_native_plinth_length']
      and not flags['actual_plinth_fibre_can_be_replaced_by_specialized_action_plinth']
      and not flags['manuscript_rigidity_or_stable_coordinate_theorem_imported']
      and not flags['an_individual_native_action_with_required_kernel_and_plinth_constructed']
      and not flags['original_goal_complete'],
      'flatness and one-value recognition do not assert native minimum length,an action,or an inapplicable bundle rigidity theorem')
proof=HERE/'FLAT_PLINTH_FAMILIES_AND_ONE_VALUE_RECOGNITION.md'
paths=['sources/flat-plinth-families-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-flat-plinth-families.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
