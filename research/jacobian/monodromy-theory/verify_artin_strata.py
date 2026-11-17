"""Exact KMFI98--99 all-A_k local representation classification."""
from pathlib import Path
from collections import Counter
import json,hashlib,time,platform
from psl11_residues import exceptional_actions,inverse
from psl11_curve_spectrum import integer_group,closure
from local_image import fixed_sheets,subgroup_kind
from artin_strata import pair_orbits,alternating,power,satisfies,spectrum,full_group_orders
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
data,matorders,triangles,a5classes,actions=exceptional_actions()
elements,I,T=integer_group(data);idx={g:i for i,g in enumerate(elements)}
orders=[matorders[g] for g in elements]
invs=[idx[data['canon'](inverse(g))] for g in elements]
perms=[[act['permutations'][g] for g in elements] for act in actions]
check('full660 group and all element orders',len(elements)==660 and
 Counter(orders)=={1:1,2:55,3:110,5:264,6:110,11:120})
allowed={a for a,o in enumerate(orders) if o in [2,3,5]}
check('429 genuine allowed meridians',len(allowed)==429 and I not in allowed)
rows=pair_orbits(I,T,invs,orders,perms)
check('295 simultaneous-conjugacy pair orbits',len(rows)==295)
covered=set();all_centers=set()
for j,row in enumerate(rows):
    a,b=row['a'],row['b']
    orbit={(T[T[g][a]][invs[g]],T[T[g][b]][invs[g]]) for g in range(660)}
    check('orbit coverage and size '+str(j),len(orbit)==row['orbit_size'] and not (covered&orbit))
    covered|=orbit
    H=closure((a,b),I,T)
    marks=[len(fixed_sheets(H,p)) for p in perms]
    check('subgroup and marks '+str(j),sorted(H)==row['subgroup'] and marks==row['marks'] and
          subgroup_kind(H,orders,marks)==row['kind'])
    c=T[a][b];d=T[b][a];o=orders[c]
    check('closed-form Artin words and period '+str(j),
      row['period']==2*o and 660%(2*o)==0 and power(c,o,I,T)==power(d,o,I,T)==I and
      all(alternating(a,b,m,I,T)==
          (power(c,m//2,I,T) if m%2==0 else T[power(c,m//2,I,T)][a])
          and alternating(b,a,m,I,T)==
          (power(d,m//2,I,T) if m%2==0 else T[power(d,m//2,I,T)][b])
          for m in range(2*o)) and
      row['residues']==[m for m in range(2*o)
           if alternating(a,b,m,I,T)==alternating(b,a,m,I,T)])
    if satisfies(row,3) and row['kind'] not in ['C2','C3','C5']:
        all_centers.add(tuple(g for g in H if all(T[g][h]==T[h][g] for h in H)))
check('all184041 ordered meridian pairs covered exactly',covered=={(a,b) for a in allowed for b in allowed})
check('all periods divide660',all(660%row['period']==0 for row in rows))
check('node spectrum matches abelian local types',set(spectrum(rows,2))=={'C2','C3','C5','V4','C6'})
cusp={'C2':55,'C3':110,'C5':264,'S3a':330,'S3b':330,'A4':1320,'A5a':1320,'A5b':1320}
check('complete ordinary cusp table',spectrum(rows,3)==cusp)
check('5049 cusp pairs',sum(cusp.values())==5049)
check('ordinary cusp meridian orders by image',
 all(row['orders']==([2,2] if row['kind'] in ['C2','S3a','S3b'] else
                       [3,3] if row['kind'] in ['C3','A4'] else [5,5])
     for row in rows if satisfies(row,3)))
check('all noncyclic ordinary cusp images centerless',all_centers=={(I,)})
check('odd Artin relation conjugates meridians',
 all(row['orders'][0]==row['orders'][1] for row in rows if any(m%2 for m in row['residues'])))
period=[]
for m in range(2,662):
    full=bool(full_group_orders(rows,m));small=bool(full_group_orders(rows,m,(2,3)))
    criterion=(m%5==0 or m%11==0) if m%2 else any(m%d==0 for d in [4,6,10,22])
    no5=m%11==0 if m%2 else any(m%d==0 for d in [10,12,22])
    check('full-G criterion residue '+str(m%660),full==criterion)
    check('no-inertia5 full-G criterion residue '+str(m%660),small==no5)
    period.append({'m':m,'residue_mod660':m%660,'spectrum':spectrum(rows,m),
       'full_G_orders':[[*k,v] for k,v in sorted(full_group_orders(rows,m).items())],
       'full_G_without5_orders':[[*k,v] for k,v in sorted(full_group_orders(rows,m,(2,3)).items())]})
check('first full-G A_k is tacnode A3',not full_group_orders(rows,2) and
 not full_group_orders(rows,3) and full_group_orders(rows,4)=={(5,5,2):2640})
check('first irreducible full-G cusp is A4',full_group_orders(rows,5)=={(5,5,5):2640})
check('no-inertia5 first full-G A9',all(not full_group_orders(rows,m,(2,3)) for m in range(2,10))
 and full_group_orders(rows,10,(2,3))=={(3,3,5):2640})
check('no-inertia5 first irreducible full-G A10',
 all(not full_group_orders(rows,m,(2,3)) for m in range(3,11,2))
 and full_group_orders(rows,11,(2,3))=={(3,3,11):1320})
for label,mat_a,mat_b,m,expected in [
 ('tacnode',(0,1,10,3),(0,1,10,8),4,[5,5,2]),
 ('ramphoid cusp',(0,1,10,3),(1,7,4,7),5,[5,5,5])]:
    a,b=idx[mat_a],idx[mat_b];H=closure((a,b),I,T)
    check('explicit '+label+' matrices',len(H)==660 and
      [orders[a],orders[b],orders[T[a][b]]]==expected and
      alternating(a,b,m,I,T)==alternating(b,a,m,I,T) and
      all(not fixed_sheets(H,p) for p in perms))
check('all generated subgroups agree with retained620 census',
 {frozenset(row['subgroup']) for row in rows}<=
 {frozenset(s['elements']) for s in json.loads((ROOT/'local_image_model.json').read_text())['subgroups']})
manifest=json.loads((ROOT/'sources/artin-strata-manifest.json').read_text())
for source in manifest['sources']:
    raw=(ROOT/source['path']).read_bytes()
    check('primary bytes '+source['path'],len(raw)==source['bytes'])
    check('primary hash '+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
model={'claims':['KMFI98','KMFI99'],'matrix_labels':'sorted canonical determinant-one matrices modulo +/-I',
 'group_order':660,'allowed_meridian_orders':[2,3,5],'ordered_pairs':184041,'pair_conjugacy_orbits':rows,
 'unbounded_relation_period':660,'all_period_spectra':period,'ordinary_cusp_spectrum':cusp,
 'ordinary_cusp_pair_count':5049,'first_full_group_A_k':3,
 'first_irreducible_full_group_A_k':4,'first_without5_full_group_A_k':9,
 'first_without5_irreducible_full_group_A_k':10,
 'actual_local_finite_group_representations':True,'local_quotient_is_Keller_realization':False,
 'every_quotient_extends_to_global_Keller_closure':False,'absolute_exceptional_action_decided':False,'goal_complete':False}
(ROOT/'artin_strata_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['ARTIN_BRANCH_STRATA.md','artin_strata.py','verify_artin_strata.py','artin_strata_model.json',
 'sources/artin-strata-manifest.json','local_image.py','local_image_model.json','psl11_residues.py',
 'psl11_curve_spectrum.py','central_lift.py']+[s['path'] for s in manifest['sources']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'runtime_seconds':round(time.monotonic()-started,3),
 'python':platform.python_version(),'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths},
 'assurance':'Written local braid/inverse-section proof; exact exhaustive295 pair orbits cover184041 pairs and all A_k via660-period. No global Keller realization, formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'artin_strata_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'pair_orbits':len(rows),'ordered_pairs':len(covered),'unbounded_A_k_period':660,
 'absolute_exceptional_decided':False,'goal_complete':False}))
