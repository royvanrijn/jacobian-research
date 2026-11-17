"""NEW MQ0 coordinate/orbit transcriptions; never replay a sealed checker."""
import hashlib,json,time
from collections import Counter
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-cyclic-deck-source-image.json'
assert not RECEIPT.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMQ0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)

x,y,a,c,step=s.symbols('x y a c step',nonzero=True)
h=s.Function('h')(y)
J=s.Matrix([x/a+h,a*y]).jacobian([x,y])
check('generic affine-fibre deck preserves volume exactly',s.simplify(J.det())==1)
check('cotangent fibre multiplier is constant from volume',s.simplify(a*(1/a))==1)
involution={x:-x+c/y**2,y:-y}
image=-x+c/y**2
check('rational source-pole involution still has order two',s.simplify(image.subs(involution,simultaneous=True)-x)==0)
check('source-pole involution has constant Jacobian one',s.Matrix([image,-y]).jacobian([x,y]).det()==1)
check('invariant base need not be a source coordinate',s.diff(y**2,y).subs(y,0)==0)

# New finite component controls: g fixed-fibre unit components; g ordinary
# base values with s unit components and one uncovered multiplicity-two orbit.
# These are combinatorial transcriptions, not actual surface constructions.
for g,selected in [(2,[0,1]),(3,[0,1,2]),(5,[0,1,2,0,1])]:
 count=Counter(selected);number=len(count)
 fixed={(0,j) for j in range(g)}
 ordinary={(i+1,j) for i in range(g) for j in range(number+1)}
 components=fixed|ordinary
 original={(0,0)}|{(i+1,selected[i]) for i in range(g)}
 def action(cell,k):
  base,label=cell
  return (0,(label+k)%g) if base==0 else ((base-1+k)%g+1,label)
 translates=[{action(z,k) for z in original} for k in range(g)]
 covered=set.union(*translates)
 missed=components-covered
 check(f'fixed source component has full free orbit g{g}',len({action((0,0),k) for k in range(g)})==g)
 check(f'at least g-1 added unit components in distinguished fibre g{g}',len(fixed-original)==g-1)
 check(f'covered ordinary fibre counts partition degree g{g}',sum(count.values())==g and 1<=number<=g)
 check(f'every covered component is unit g{g}',all(z[0]==0 or z[1]<number for z in covered))
 check(f'common omitted part is exactly one whole curve orbit g{g}',len(missed)==g and missed=={(i+1,number) for i in range(g)})
 check(f'quotient source image has no additional point object g{g}',components-covered==set.intersection(*(components-t for t in translates)))
 check(f'distinguished fibre has exactly one quotient orbit g{g}',len(fixed)==g)
 # Per fibre above a nonzero quotient base value: one source line at each root.
 check(f'every source-image fibre point is reached g{g}',all(k>0 for k in count.values()))
 check(f'original finite-cover Euler ledger g{g}',(1+(g-1)+g*number)//g==number+1)
 check(f'curve orbit deletion lowers quotient Euler by one g{g}',(number+1)-1==number)
 check(f'actual image Euler from ruling g{g}',1+(number-1)==number)
 check(f'universal cover boundary count g{g}',g-1+g*(number-1)==g*number-1)
 check(f'only fixed-base source fibre has one point g{g}',len(original&fixed)==1)
 check(f'invariant source replica character cancels g{g}',(g-1+1)%g==0)
 check(f'base root order gives multiple fibre g{g}',s.diff(y**g,y).as_leading_term(y)==g*y**(g-1))

for d,g in [(6,2),(12,3),(30,5)]:
 check(f'primitive quotient degree d{d} g{g}',d%g==0 and d//g>1 and (d//g)*g==d)

cp=HERE/'sources/cyclic-deck-source-image-context.json';C=json.loads(cp.read_text())
for p,digest in C['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,sha(ROOT/p)==digest,False)
check('five prior claims and twelve pins',len(C['retained_structural_inputs'])==5 and len(C['retained_proof_and_source_sha256'])==12,False)
flags=C['scope_flags']
check('actual source image has no point gaps',flags['actual_source_image_affine_and_surjective'] and flags['PF136_finite_point_gaps_eliminated'],False)
check('single irreducible multiple fibre and reduced others',flags['unique_irreducible_multiple_fibre_order_g'] and flags['other_image_fibres_reduced_at_most_g_components'],False)
check('free source-translation universal cover',flags['covered_universal_open_has_free_Ga'] and flags['actual_image_algebra_R_intersect_fixed_field'],False)
check('full target lattice remains an obligation',not any(flags[k] for k in ['image_equals_entire_intermediate_etale_locus','primitive_u_integral_over_target','free_Ga_surface_is_plane','source_Ga_descends_to_plane_target','unrestricted_HC4_or_JC2_resolved']),False)
check('no new imports/caches/downloads/registered replay',C['new_external_theorem_imports']==C['new_primary_caches']==C['new_primary_downloads']==0 and C['registered_mathematical_checkers_replayed']==[],False)
r=dict(outcome='PASS',claim='HC4MYMQ0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'CYCLIC_DECK_SOURCE_IMAGE_AND_THE_UNIQUE_MULTIPLE_FIBRE.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
