"""Exact KMFI94--97 subgroup marks, image controls and arithmetic normalizer."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
from collections import Counter
from itertools import product
import json,hashlib,platform,time
import sympy as S
from psl11_residues import exceptional_actions,inverse
from psl11_curve_spectrum import integer_group,closure
from local_image import complete_subgroups,conjugate,fixed_sheets,orbit_lengths,subgroup_kind
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(x):return S.cancel(x)==0
data,matorders,triangles,a5classes,actions=exceptional_actions()
elements,I,T=integer_group(data);idx={g:i for i,g in enumerate(elements)}
orders=[matorders[g] for g in elements]
invs=[idx[data['canon'](inverse(g))] for g in elements]
perms=[[a['permutations'][g] for g in elements] for a in actions]
check("PSL2(11) group order660",len(elements)==660)
check("complete element order counts",Counter(orders)=={1:1,2:55,3:110,5:264,6:110,11:120})
groups,free,transcript,extra=complete_subgroups(I,T,invs,orders)
check("complete pruned C2 overgroup sizes",Counter(map(len,free))=={2:1,4:3,6:6,10:6,12:3,60:6})
check("all10244 pruned extension pairs tested",len(transcript)==10244 and
 {(H,g) for H,g,E,bad in transcript}==
 {(H,g) for H in free for g in range(660) if orders[g] in [1,2,3,5] and g not in H})
check("all8150 bad extensions rejected",sum(E is None for H,g,E,bad in transcript)==8150)
check("every accepted extension retained and every rejection has derangement witness",
 all((E in free and H<E and g in E) if E is not None else orders[bad] in [6,11]
   for H,g,E,bad in transcript))
check("pruned maximal groups all order60",
 all(len(H)==60 for H in free if not any(H<J for J in free)))
for o,sizes,count in [(6,[6,12,660],1302),(11,[11,55,660],1254)]:
 found,edges=extra[o]
 check("complete overgroups order"+str(o),sorted(map(len,found))==sizes and len(edges)==count)
check("full subgroup census620",len(groups)==620)
expected_counts={'1':1,'C2':55,'C3':55,'V4':55,'C5':66,'S3a':55,'S3b':55,
 'C6':55,'D10':66,'C11':12,'A4':55,'D12':55,'F55':12,'A5a':11,'A5b':11,'G':1}
expected_marks={'1':[11,11],'C2':[3,3],'C3':[2,2],'V4':[1,1],'C5':[1,1],
 'S3a':[2,0],'S3b':[0,2],'C6':[0,0],'D10':[1,1],'C11':[0,0],
 'A4':[1,1],'D12':[0,0],'F55':[0,0],'A5a':[1,0],'A5b':[0,1],'G':[0,0]}
ordered=sorted(groups,key=lambda H:(len(H),tuple(sorted(H))));gids={H:j for j,H in enumerate(ordered)}
marks={H:[len(fixed_sheets(H,p)) for p in perms] for H in groups}
kinds={H:subgroup_kind(H,orders,marks[H]) for H in groups}
check("all16 subgroup-class type counts",Counter(kinds.values())==expected_counts)
check("all620 marks in both actions",all(marks[H]==expected_marks[kinds[H]] for H in groups))
classes=[]
for kind,count in expected_counts.items():
 H=next(H for H in groups if kinds[H]==kind)
 conjugates={conjugate(H,g,T,invs) for g in range(660)}
 check("complete conjugacy class "+kind,len(conjugates)==count and
   conjugates=={J for J in groups if kinds[J]==kind})
 row={'kind':kind,'order':len(H),'count':count,'marks':marks[H],
  'representative_id':gids[H],'orbit_profiles':[list(orbit_lengths(H,p)) for p in perms]}
 classes.append(row)
 if kind in ['S3a','S3b']:
  good=0 if kind=='S3a' else 1
  check("S3 good orbit profile "+kind,orbit_lengths(H,perms[good])==(1,1,3,6))
  check("S3 bad orbit profile "+kind,orbit_lengths(H,perms[1-good])==(2,3,3,3))
  check("Burnside orbit counts agree while common fixed sheets differ "+kind,
   len(orbit_lengths(H,perms[0]))==len(orbit_lengths(H,perms[1]))==4 and marks[H] in [[2,0],[0,2]])
for j in [0,1]:
 bad=[H for H in groups if marks[H][j]==0]
 minimal=[H for H in bad if not any(J<H for J in bad)]
 target={'C6','C11','S3b' if j==0 else 'S3a'}
 check("all inclusion-minimal bad subgroups action"+str(j),
  len(minimal)==122 and Counter(kinds[H] for H in minimal)=={k:expected_counts[k] for k in target})
 check("all201 bad subgroups exactly contain a minimal obstruction action"+str(j),
  len(bad)==201 and all((marks[H][j]==0)==any(J<=H for J in minimal) for H in groups))
abelian=[H for H in groups if all(T[a][b]==T[b][a] for a,b in product(H,repeat=2))]
check("all abelian subgroup types",set(kinds[H] for H in abelian)=={'1','C2','C3','V4','C5','C6','C11'})
snc=[H for H in abelian if closure(tuple(g for g in H if orders[g] in [2,3,5]),I,T)==H]
check("complete possible affine normal-crossing local types",set(kinds[H] for H in snc)=={'1','C2','C3','V4','C5','C6'})
for xorder,yorder in [(2,5),(3,5)]:
 check("no commuting mixed order"+str((xorder,yorder)),not any(
  T[a][b]==T[b][a] for a in range(660) if orders[a]==xorder for b in range(660) if orders[b]==yorder))
# Exact automorphism bound and outer representative, no S11 enumeration.
a=next(i for i in range(660) if orders[i]==2)
candidates=[g for g in range(660) if orders[g]==3 and orders[T[a][g]]==11]
b=candidates[0]
check("twenty-four triangle images for fixed involution",len(candidates)==24)
check("all triangle image pairs generate full group",all(len(closure((a,g),I,T))==660 for g in candidates))
check("trivial group center", {g for g in range(660) if T[g][a]==T[a][g] and T[g][b]==T[b][g]}=={I})
outer=[idx[data['canon']((g[0],2*g[1]%11,6*g[2]%11,g[3]))] for g in elements]
check("outer matrix map bijective",len(set(outer))==660)
check("outer map preserves full multiplication table",all(outer[T[a][b]]==T[outer[a]][outer[b]] for a,b in product(range(660),repeat=2)))
for ka,kb in [('A5a','A5b'),('S3a','S3b')]:
 check("outer automorphism exchanges "+ka+" and "+kb,
  {frozenset(outer[g] for g in H) for H in groups if kinds[H]==ka}=={H for H in groups if kinds[H]==kb})
check("automorphism upper equals explicit inner-outer lower",55*24==2*660==1320)
check("point stabilizers fix exactly one sheet, trivial permutation centralizer",
 all(marks[frozenset(idx[g] for g in action['H'])][j]==1 for j,action in enumerate(actions)))
# Cubic actual-map control and a global polynomial omission certificate.
x,y,z,A,B,C,R=S.symbols('x y z A B C R')
u=1+x*y;gamma=1-S.Rational(3,2)*x*y+x*x*z
F=[S.cancel((u+u*u-2*u**3*gamma)/(x*x)),
 S.cancel((1+2*u-3*u*u*gamma)/x),x*gamma]
check("actual cubic coordinates polynomial",all(S.denom(f)==1 for f in F))
check("actual cubic full Jacobian1",zero(S.Matrix(F).jacobian([x,y,z]).det()-1))
P0=S.expand(3*F[1]*F[2]-1);Q0=S.expand(27*F[0]*F[2]**2-1)
SS=3*u*gamma-1;LL=u*SS-1
check("omission identity P0",zero(P0-3*gamma+SS**2))
check("omission identity Q0",zero(Q0-SS**3-3*P0*(SS+1)))
check("omission identity unimodular step",zero(1-u*P0-SS*LL))
co1=3*u-3*u*u*P0+u**3*P0**2-3*(SS+1)*LL**3
check("exact pulled-back complete-intersection Bezout certificate",zero(P0*co1+Q0*LL**3-1))
sigma={A:1/(27*C*C),B:1/(3*C)}
poly=C*R**3-R*R+B*R-A
check("omitted locus exactly triple-root coefficients",zero(poly.subs(sigma)-C*(R-1/(3*C))**3))
targetP=3*B*C-1;targetQ=27*A*C*C-1
check("omitted curve transverse complete intersection",
 zero(S.Matrix([targetP,targetQ]).jacobian([B,A]).det()-81*C**3))
pp=targetP/(3*C*C);qq=(3*targetP-targetQ)/(27*C**3)
tt=S.symbols('tt')
check("local depressed cubic coordinates",zero(poly.subs(R,tt+1/(3*C))/C-(tt**3+pp*tt+qq)))
check("C0 root at infinity simple",
 S.diff(C-R+B*R*R-A*R**3,R).subs({C:0,R:0})==-1)
# The two branch transpositions in the cubic standard action generateS3.
p12=(1,0,2);p23=(0,2,1);seen={tuple(range(3))};todo=list(seen)
while todo:
 p=todo.pop()
 for q in [p12,p23]:
  val=tuple(p[q[i]] for i in range(3))
  if val not in seen:seen.add(val);todo.append(val)
check("cubic cusp localS3 has no common fixed sheet",len(seen)==6 and not any(all(g[i]==i for g in seen) for i in range(3)))
manifest=json.loads((ROOT/'sources/local-image-manifest.json').read_text())
for source in manifest['sources']:
 raw=(ROOT/source['path']).read_bytes()
 check("primary bytes "+source['path'],len(raw)==source['bytes'])
 check("primary hash "+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
free_ordered=sorted(free,key=lambda H:(len(H),tuple(sorted(H))));freeids={H:j for j,H in enumerate(free_ordered)}
records=[[freeids[H],g,freeids[E] if E is not None else None,bad] for H,g,E,bad in transcript]
model={'claims':['KMFI94','KMFI95','KMFI96','KMFI97'],
 'matrix_labels':'sorted canonical determinant-one matrices modulo +/-I',
 'pure_normalization_boundary':True,'image_fiber_bound':'actual geometric fiber injects into common fixed sheets of strict local group',
 'fixed_sheet_sufficiency_without_boundary_condition':False,
 'equality_condition':'no omitted unramified prime divisors',
 'complete_subgroup_count':len(groups),'complete_conjugacy_classes':classes,
 'subgroups':[{'id':gids[H],'kind':kinds[H],'elements':sorted(H),'marks':marks[H]} for H in ordered],
 'pruned_C2_overgroups':[{'id':freeids[H],'generators':list(free[H]),'elements':sorted(H)} for H in free_ordered],
 'pruned_extensions':records,'all_other_seed_extension_counts':{o:len(edges) for o,(found,edges) in extra.items()},
 'minimal_bad_classes':{'action_a':['C6','C11','S3b'],'action_b':['C6','C11','S3a']},
 'same_element_characters_do_not_determine_local_fiber_marks':True,
 'every_subgroup_claimed_as_realizable_local_group':False,
 'normal_crossing_local_types':['1','C2','V4','C3','C5','C6'],
 'mixed_2_3_normal_crossing_strata':'omitted','mixed_5_with2_or3_normal_crossing_strata':'impossible for the G-cover',
 'cubic_omitted_target_curve':['3*B*C-1','27*A*C^2-1'],'cubic_target_curve':'smooth G_m, exact image complement',
 'arithmetic_normalizer_in_S11_order':660,'automorphism_group_order':1320,
 'geometric_exceptional11_forces_same_arithmetic_group_and_no_new_constants':True,
 'local_obstruction_promoted_to_absolute_action_exclusion':False,
 'exceptional11_absolute_action_decided':False,'full_action_classification_complete':False}
(ROOT/'local_image_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['LOCAL_IMAGE_AND_MARKS.md','local_image.py','verify_local_image.py','local_image_model.json',
 'sources/local-image-manifest.json','BOUNDARY_CLASS_OBSTRUCTION.md','psl11_curve_spectrum.py','psl11_residues.py','central_lift.py']+[s['path'] for s in manifest['sources']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written pure-boundary and strict-local geometry; exhaustive620 subgroup mark/minimal-obstruction census, complete finite automorphism bound and exact actual cubic omission controls. No formal verification, independent whole-proof review or absolute exceptional-action exclusion.',
 'goal_complete':False}
(ROOT/'local_image_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'subgroups':620,'classes':16,'S3_classes_distinguish_actions':True,
 'arithmetic_group_forced_PSL11':True,'absolute_exceptional_decided':False,'goal_complete':False}))
