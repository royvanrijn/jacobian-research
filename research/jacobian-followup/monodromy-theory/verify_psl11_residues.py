"""Complete finite action/decomposition and Schur-residue controls KMFI81--84.

The affine purity and transfer arguments are written in the proof note.
No passed local or global necessary gate is promoted to occurrence.
"""
from pathlib import Path
from itertools import product,combinations
from collections import Counter
from math import gcd
import hashlib,json,platform,time
from central_lift import generated
from psl11_residues import (exceptional_actions,inverse,centralizer,cyclic_overgroups,
 local_profile,schur_character,permutation_sign_on,residue_visibility_gate)
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def cycle_type(p):
    left=set(range(len(p)));out=[]
    while left:
        a=min(left);b=a;k=0
        while b in left:left.remove(b);k+=1;b=p[b]
        out.append(k)
    return tuple(sorted(out))
data,orders,triangles,classes,actions=exceptional_actions()
G=data['quotient'];I=data['qI'];mul=data['qmul'];smul=data['mul']
check("full PSL2(11) order660",len(G)==660)
check("full group element order counts",dict(Counter(orders.values()))=={1:1,2:55,3:110,5:264,6:110,11:120})
check("A5 triangle subgroups through fixed involution complete",len(triangles)==6)
check("two A5 conjugacy classes of size11",list(map(len,classes))==[11,11] and not classes[0]&classes[1])
check("all enumerated A5 triangle subgroups accounted for",triangles<=classes[0]|classes[1])
expected_cycles={1:(1,)*11,2:(1,1,1,2,2,2,2),3:(1,1,3,3,3),5:(1,5,5),6:(2,3,6),11:(11,)}
for action_index,action in enumerate(actions):
    perms=action['permutations'];H=action['H']
    check(f"exceptional action faithful degree11 #{action_index}",len(set(perms.values()))==660 and all(len(p)==11 for p in perms.values()))
    check(f"exceptional stabilizer A5 order60 #{action_index}",len(H)==60 and {g for g in G if perms[g][0]==0}==set(H))
    check(f"all660 exceptional cycle types #{action_index}",all(cycle_type(perms[g])==expected_cycles[orders[g]] for g in G))
    HH={a for a in data['elements'] if data['canon'](a) in H}
    comm={smul(smul(smul(a,b),inverse(a)),inverse(b)) for a,b in product(HH,repeat=2)}
    derived=generated(comm,data['I'],smul)
    check(f"lifted A5 subgroup order120 #{action_index}",len(HH)==120)
    check(f"lifted A5 subgroup perfect; no cyclic quotient #{action_index}",len(derived)==120)
check("both exceptional actions have identical complete cycle characters",
 all(cycle_type(actions[0]['permutations'][g])==cycle_type(actions[1]['permutations'][g]) for g in G))
expect={
 (2,2):{'tails':(1,1,1),'ram':[(2,1)]*4,'residue':False,'label':'C2'},
 (2,4):{'tails':(1,2),'ram':[(2,1),(2,1),(2,2)],'residue':True,'label':'V4'},
 (2,6):{'tails':(3,),'ram':[(2,1),(2,3)],'residue':False,'label':'C6'},
 (2,12):{'tails':(3,),'ram':[(2,1),(2,3)],'residue':True,'label':'D12 order12'},
 (3,3):{'tails':(1,1),'ram':[(3,1)]*3,'residue':False,'label':'C3'},
 (3,6):{'tails':(2,),'ram':[(3,1),(3,2)],'residue':False,'label':'C6'},
 (5,5):{'tails':(1,),'ram':[(5,1)]*2,'residue':False,'label':'C5'}}
rows=[];local_records=[];seen_rows=set();representatives={}
for inertia_order,centralizer_order,subgroups_count in [(2,12,55),(3,6,55),(5,5,66)]:
    a=min(g for g in G if orders[g]==inertia_order)
    inertia=generated([a],I,mul);representatives[inertia_order]=a
    C=centralizer(G,a,mul)
    check(f"centralizer exact order inertia{inertia_order}",len(C)==centralizer_order)
    allcyc={generated([g],I,mul) for g in G if orders[g]==inertia_order}
    conjugatecyc={frozenset(mul(mul(g,h),data['canon'](inverse(g))) for h in inertia) for g in G}
    check(f"all inertia subgroups covered by one conjugacy orbit order{inertia_order}",
      conjugatecyc==allcyc and len(allcyc)==subgroups_count)
    over=cyclic_overgroups(inertia,C,I,mul)
    expected_sizes={2:[2,4,4,4,6,12],3:[3,6],5:[5]}[inertia_order]
    check(f"complete centralizer overgroups order{inertia_order}",sorted(map(len,over))==expected_sizes)
    # Normalizer elements that invert odd inertia are not valid decomposition
    # elements over algebraically closed characteristic-zero constants.
    if inertia_order in [3,5]:
        N={g for g in G if frozenset(mul(mul(g,h),data['canon'](inverse(g))) for h in inertia)==inertia}
        check(f"odd inertia normalizer strictly larger than tame centralizer order{inertia_order}",len(N)==2*len(C))
        check(f"excluded normalizer elements invert inertia order{inertia_order}",
          all(mul(mul(g,a),data['canon'](inverse(g)))==data['canon'](inverse(a)) for g in N-set(C)))
    char=schur_character(data,a,C) if inertia_order==2 else {g:0 for g in C}
    if inertia_order==2:
        check("involution lifts have order4",smul(a,a)==data['Z'])
        check("Schur character homomorphism on full involution centralizer",
          all(char[mul(g,h)]==(char[g]+char[h])%2 for g,h in product(C,repeat=2)))
        check("Schur character trivial on inertia",all(char[g]==0 for g in inertia))
        check("Schur character kernel has cyclic order6",
          len([g for g in C if not char[g]])==6 and max(orders[g] for g in C if not char[g])==6)
    for action_index,action in enumerate(actions):
        perms=action['permutations']
        if inertia_order==2:
            fixed={s for s in range(11) if perms[a][s]==s}
            check(f"Schur residue is full fixed-tail discriminant character #{action_index}",
              all(char[g]==permutation_sign_on(perms[g],fixed) for g in C))
        for D in sorted(over,key=lambda H:(len(H),tuple(sorted(H)))):
            key=(inertia_order,len(D));expected=expect[key]
            profile=local_profile(D,inertia,perms);orbits=profile['tail_orbits']
            tail_degrees=tuple(sorted(map(len,orbits)));ram=[ef for ef in profile['ef'] if ef[0]>1]
            tag=f"I{inertia_order},D{len(D)},action{action_index},rep{min(D)}"
            check("complete unramified residue orbit degrees "+tag,tail_degrees==expected['tails'])
            check("complete missing ramified prime indices "+tag,ram==expected['ram'])
            check("local degree sum11 "+tag,sum(e*f for e,f in profile['ef'])==11)
            target=any(char[g] for g in D)
            check("exact nonzero target Schur residue "+tag,target==expected['residue'])
            visible=[]
            for O in orbits:
                s=min(O);stab={g for g in D if perms[g][s]==s}
                value=any(char[g] for g in stab);visible.append(value)
                check("source residue stabilizer index "+tag+str(s),len(D)//len(stab)==len(O))
                if target:
                    check("nonzero residue survives exactly on odd source orbits "+tag+str(s),value==(len(O)%2==1))
            if key==(2,4):
                HH={g for g in data['elements'] if data['canon'](g) in D}
                check("V4 preimage is quaternion "+tag,
                  len(HH)==8 and sum(smul(g,g)==data['I'] for g in HH)==2)
            if key==(2,12):
                fixed=profile['fixed']
                action_on_fixed={tuple(perms[g][s] for s in sorted(fixed)) for g in D}
                check("D12 residual action is S3 "+tag,len(action_on_fixed)==6)
            if (key,action_index) not in seen_rows:
                seen_rows.add((key,action_index))
                patterns=[]
                for mask in range(1,1<<len(orbits)):
                    ids=[j for j in range(len(orbits)) if mask>>j&1]
                    rec={'inertia_order':inertia_order,'decomposition':expected['label'],
                      'target_residue_nonzero':target,
                      'retained_residue_degrees':sorted(len(orbits[j]) for j in ids),
                      'source_residues_nonzero':[visible[j] for j in ids]}
                    patterns.append(rec)
                    if action_index==0:local_records.append(rec)
                rows.append({'action':action_index,'inertia_order':inertia_order,'D':expected['label'],
                  'unramified_degrees':list(tail_degrees),'ramified_ef':[list(ef) for ef in ram],
                  'Schur_residue_nonzero':target,'locally_permitted_source_patterns':patterns})
check("all seven local types verified in both exceptional actions",len(rows)==14)
hidden=next(r for r in local_records if r['decomposition']=='V4' and r['retained_residue_degrees']==[2])
visible=next(r for r in local_records if r['decomposition']=='V4' and r['retained_residue_degrees']==[1])
dihedral=next(r for r in local_records if r['decomposition']=='D12 order12')
cyclic=next(r for r in local_records if r['inertia_order']==2 and r['decomposition']=='C6')
check("single hidden V4 residue rejected globally",not residue_visibility_gate([hidden]))
check("any number of hidden V4 residues rejected globally",all(not residue_visibility_gate([hidden]*k+[cyclic]) for k in range(1,8)))
check("visible V4 residue meets this necessary gate",residue_visibility_gate([visible]))
check("D12 odd-degree source orbit meets this necessary gate",residue_visibility_gate([hidden,dihedral]))
check("all cyclic decomposition rows meet this necessary gate",residue_visibility_gate([cyclic]))
check("degree11 restriction injective on2torsion arithmetic condition",gcd(11,2)==1 and 11%2==1)
check("even degree cannot use same injectivity step",12%2==0)
# Explicit normalized Q8 central cocycle and its edge character.
V=list(product(range(2),repeat=2))
add=lambda x,y:((x[0]+y[0])%2,(x[1]+y[1])%2)
c=lambda x,y:(x[0]*y[0]+x[1]*y[1]+x[1]*y[0])%2
check("quaternion normalized central cocycle",all(c((0,0),x)==c(x,(0,0))==0 for x in V))
check("quaternion cocycle identity exhaustive",all(
 (c(x,y)+c(add(x,y),z)-c(y,z)-c(x,add(y,z)))%2==0 for x,y,z in product(V,repeat=3)))
check("quaternion mixed commutator character",
 all((c((1,0),x)-c(x,(1,0)))%2==x[1] for x in V))
check("quadratic residual extension kills mixed character",
 all((c((1,0),x)-c(x,(1,0)))%2==0 for x in V if x[1]==0))
b=lambda a:(a//2)%2
check("cyclic C4 inertia cocycle inflates to coboundary on procyclic quotient",
 all((b(a)+b(d)-b((a+d)%4))%2==(a%2)*(d%2) for a,d in product(range(4),repeat=2)))
# A (2,3,11) abstract curve passport passes the new residue gate.
# It is deliberately not an ordinary-affine realization.
a2=representatives[2];passport=None
for a3 in G:
    if orders[a3]==3 and orders[mul(a2,a3)]==11:
        subgroup=generated([a2,a3],I,mul)
        if len(subgroup)==660:passport=(a2,a3);break
check("abstract2311 passport exists and generates full group",passport is not None)
p2,p3=[actions[0]['permutations'][a] for a in passport]
p11=actions[0]['permutations'][mul(*passport)]
indices=[11-len(cycle_type(a)) for a in [p2,p3,p11]]
check("abstract2311 root passport genus0",sum(indices)==2*11-2)
manifest=json.loads((ROOT/'sources/psl11-residues-manifest.json').read_text())
for source in manifest['retained_primary_inputs']:
    path=ROOT/source['path']
    check("retained primary bytes "+path.name,path.stat().st_size==source['bytes'])
    check("retained primary hash "+path.name,hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256'])
model={'claims':['KMFI81','KMFI82','KMFI83','KMFI84'],
 'A5_actions_checked':2,'complete_A5_conjugacy_classes':[11,11],
 'local_tame_tables':rows,
 'universal_statements':{
 'prime_to_degree_Brauer_restriction_visible_on_affine_source':'proved for every prime l not dividing generic degree; over nonclosed fields exclude constant classes',
 'nonclosed_coefficients_require_constant_Brauer_qualification':True,
 'tame_decomposition_centralizes_inertia':True,
 'Schur_residue_support':'exactly I=C2 with D=V4 or D12 order12',
 'Schur_residue_character':'discriminant of full rank3 unramified residue algebra',
 'proper_SL2_11_field_lift_criterion':'all affine involution decomposition groups cyclic',
 'nonzero_Schur_residue_requires_some_present_odd_residue_source_prime':True,
 'only_quadratic_V4_source_orbits_with_no_D12':'impossible if any V4 divisor occurs',
 'every_local_quadratic_V4_choice_individually_excluded':False,
 'decomposition_group_itself_required_to_fix_sheet':False,
 'unramified_tails_assumed_all_present':False,
 'proper_SL2_lift_claimed_as_Keller_occurrence_or_etale_source_cover':False},
 'controls':{'locally_allowed_but_globally_rejected':hidden,
   'visible_pattern_meets_only_necessary_gate':visible,
   'abstract2311_passport_generators':[list(a) for a in passport],
   'abstract2311_genus0':True,'abstract2311_ordinary_affine_completion_provided':False,
   'lifted_A5_preimages_perfect':True},
 'exceptional_PSL2_11_degree11_absolute_occurrence_decided':False,
 'natural_PSL2_11_degree12_absolute_occurrence':'excluded by KMFI67',
 'full_action_classification_complete':False}
(ROOT/'psl11_residue_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['PSL11_RESIDUE_OBSTRUCTION.md','psl11_residues.py','verify_psl11_residues.py',
 'psl11_residue_model.json','sources/psl11-residues-manifest.json']+[s['path'] for s in manifest['retained_primary_inputs']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),
 'runtime_seconds':round(time.monotonic()-started,3),
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written all-dimension geometric visibility and Schur-residue proofs using retained primary inputs; complete finite PSL2(11) group, both A5 actions, tame overgroups, all residue orbits and local cocycle controls. No formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'psl11_residue_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'complete_tame_rows':7,'all_hidden_quadratic_V4_pattern_excluded':True,
 'exceptional_PSL11_absolute_decided':False,'goal_complete':False}))
