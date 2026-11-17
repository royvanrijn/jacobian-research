#!/usr/bin/env python3
"""Register complete weighted tail rigidity and the sharp marked degree bound."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;BASE='jacobian-followup/cubic-reduction/'
sys.path.insert(0,str(ROOT/'scripts'))
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
ledger=load('MATH_STATUS.json');known={e['id']:e for e in ledger['entries']}
v=load(BASE+'weighted-tail7/validation.json');w=load(BASE+'weighted-newfirst6/validation.json')
assert v['status']=='pass' and v['source_sha256']==sha(HERE/'verify_weighted_tail7.py') and v['branch_unit_identities']==3
assert w['status']=='pass' and w['source_sha256']==sha(HERE/'verify_weighted_newfirst6.py') and w['tail_classification_sha256']==sha(HERE/'weighted-tail7/validation.json')
def register(identifier,title,scope,source,checker,deps,files):
 assert identifier not in known or known[identifier]['canonical_source']==BASE+source
 e=copy.deepcopy(known['KCUB10']);e.update(id=identifier,title=title,scope=scope,canonical_source=BASE+source,checker=BASE+checker,dependencies=deps,
  artifact_hash='sha256:'+sha(HERE/checker),software_lock=[str((HERE/f).resolve().relative_to(ROOT)) for f in files],supersedes_notes=[],supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
 if identifier in known:known[identifier].clear();known[identifier].update(e)
 else:ledger['entries'].append(e);known[identifier]=e
 return e
files=['verify.py','verify_power_single_compress.py','weighted_verify_utils.py','weighted_tail7.py','weighted-tail7/input.json','weighted-tail7/validation.json',
 'weighted_tail7_small_h.py','weighted-tail7/small-h/input.json','weighted_tail7_high_w.py','weighted-tail7/high-w/generic/input.json',
 'weighted_tail7_mid_w.py','weighted-tail7/mid-w/input.json','certify_weighted_units.py','weighted-unit-certificates/manifest.json',
 'verify_raw13_rigidity.py','raw13-rigidity/validation.json','raw13-specialization/input.json','raw13-normalized/input.json','requirements.txt','graph-requirements.txt']
for tag in ['small-h','high-w','mid-w']:files += ['weighted-unit-certificates/'+tag+s for s in ['-input.json','.json','.sing']]
e11=register('KCUB11','Complete19-atom negative-weight tail slice is rigid',
 'Over every characteristic-zero field, fix F1=a-3a^2b-a^3c, identity linear part, source/output weights(1,-1,-2) and a collision between(0,beta,c0) and(-1,beta,1+3beta) with beta nonzero. Allow every compatible nonlinear tail monomial of ordinary degree at most7: eight in F2 and nine in F3, for19 possible nonlinear atoms in total. The unique map is the established raw13-term representative, with beta1 and c0=-8; all six added tail coefficients vanish. Setting u=ab,w=a^2c and v=1-3u-w reduces the full Jacobian equation to (G+v G_v)H_u-G_u(2H+v H_v)=1. An elementary polynomial Wronskian degree argument forces boundary degrees at most1 and2. The top-v equations give an exhaustive case split, with a degree-free Euler-vector-field exclusion of linear g and three rational unit identities for the remaining branches. A separate standard-library Fraction audit reconstructs the full three-variable universal determinant, boundary equations, complete branch linear parameterizations, all15665 multiplier monomials (maximum degree7), top-v identities, linear-g coordinate change and the reduction to independently checked KCUB10 equations. Direct generic Groebner and real sparse queries remain stored as timeout/UNKNOWN; they are not the proof. This is complete only for the stated fixed-first-output weight/collision slice. It does not cover changing the first output, arbitrary weights, collision configurations or global raw12-term existence.',
 'WEIGHTED_TAIL_CLASSIFICATION.md','verify_weighted_tail7.py',['KCUB10'],files)
e12=register('KCUB12','Degree seven is sharp in a marked negative-weight class',
 'For three-variable characteristic-zero Keller maps with identity linear part, source/output weights(1,-1,-2), first-output coefficients-3 on a^2b and-1 on a^3c, and a collision between(0,1,c0) and(-1,1,4), the minimum ordinary degree is exactly7. The established raw13-term map supplies the upper bound. The degree<=6 model includes all17 compatible nonlinear atoms. Its first collision equation sets F1=a-3a^2b-a^3c+theta*(4a^3b^2+a^4bc); two exact linear eliminations leave13 coefficient parameters. A lifted generic Q(theta) unit clears to the explicit identity theta^7=sum M_i I_i with51 multiplier monomials of maximum degree6. Independent Fraction replay reconstructs the full three-variable determinant/collision coefficient equations, the complete linear parameterization, and this cleared identity, so the sole exceptional parameter is theta0 over any characteristic-zero field. That zero slice embeds into KCUB11, which requires the nonzero coefficient8 on a^3b^3c in F3, absent at degree6. A direct unit lift at zero timed out and is not used; no timeout is treated as UNSAT. Sharpness is only in the stated marked weighted class, not a global three-variable degree minimum, arbitrary degree-six exclusion or raw support optimum. The separate ordinary-degree minimum across varying dimensions remains3, and the larger23-atom degree7 raw12-term query is still UNKNOWN.',
 'WEIGHTED_MARKED_DEGREE.md','verify_weighted_newfirst6.py',['KCUB11'],['verify.py','verify_power_single_compress.py','weighted_verify_utils.py','weighted_newfirst.py','weighted-newfirst6/input.json','weighted_newfirst6_parameter.py',
 'weighted_newfirst6_lift.py','weighted-newfirst6/parameter-field/lift.sing','weighted-newfirst6/parameter-field/certificate.json','weighted-newfirst6/validation.json','verify_weighted_tail7.py','weighted-tail7/validation.json','weighted-tail7/input.json','raw13-specialization/input.json','requirements.txt','graph-requirements.txt'])
save('MATH_STATUS.json',ledger)
lessons=load('knowledge/lessons.json')
def put(r):
 r['programme_status']='archived';r['reviewed_claims']={k:claim_fingerprint(known[k]) for k in r['claims']}
 old=next((x for x in lessons['lessons'] if x['id']==r['id']),None)
 if old is None:lessons['lessons'].append(r)
 else:old.clear();old.update(r)
put(dict(id='METHOD-KELLER-WEIGHTED-TAIL-CASE-SPLIT-20261114',area='geometry',title='Use boundary degrees and top-v branches before a dense weighted coefficient solve',
 when='Enlarging the established raw13 support while fixing its first output and negative source/output weights.',
 use='Enumerate every allowed tail atom. Reduce the determinant through u=ab,w=a^2c, then use v=1-3u-w. The boundary Wronskian forces low polynomial degrees; extract the top-v equations to split zero/nonzero leading coefficients. Handle the linear-g case by an Euler-vector-field argument, normalize each linear branch completely, and retain explicit unit identities. Independently reconstruct the full three-source-variable determinant rather than trusting only the quotient formula.',
 avoid='Repeating timed-out14-variable generic solves after exact structural cases are available; confusing real NLSat UNKNOWN with complex UNSAT; dropping beta under diagonal normalization; claiming global sparsity from one fixed-first-output weight slice.',
 boundary='KCUB11 classifies the full19-atom degree7 tail slice with first output a-3a^2b-a^3c and beta!=0 collision. All six additions vanish. Three independently replayed unit identities contain15665 multiplier monomials and degree at most7. The proof also checks the top-v coefficient identities and degree-free linear-g obstruction. Generic full-ring20-second and real40-second attempts remain timeout/UNKNOWN artifacts.',
 reopen_when='Change the first output, weights or collision configuration. The23-atom degree7 marked model with new first-output terms and raw support<=12 remains UNKNOWN; use parameter-field elimination in its three first-output parameters instead of repeating the same bounded real query.',
 claims=['KCUB11','KCUB10'],sources=[BASE+'WEIGHTED_TAIL_CLASSIFICATION.md',BASE+'EXPLORATION.md'],implementation=[BASE+n for n in ['weighted_tail7.py','weighted_tail7_reduced.py','weighted_tail7_nra.py','weighted_tail7_small_h.py','weighted_tail7_high_w.py','weighted_tail7_mid_w.py','certify_weighted_units.py','verify_weighted_tail7.py']],keywords=['Keller','weights1,-1,-2','19 atoms','Wronskian','boundary degrees','top-v equations','three unit identities']))
put(dict(id='METHOD-KELLER-PARAMETER-CERTIFICATE-EXCEPTIONS-20261114',area='geometry',title='Clear a parameter-field unit certificate and resolve every exceptional root',
 when='A full coefficient-ring solve stalls while rational parameter slices quickly give unit ideals.',
 use='Treat the parameter as a coefficient-field variable, lift the generic unit and export coefficients monomial by monomial. Clear the least common denominator over the rational polynomial ring and verify the resulting D(theta)=sum M_i I_i coefficient by coefficient. Only after checking every root of D may the parameter family be excluded. Reuse a proved containing class for an exceptional slice when a direct certificate lift is expensive.',
 avoid='Inferring global exclusion from a generic Q(theta) unit or sampled rational values; omitting algebraic exceptional roots; rerunning an unchanged timeout; treating the failed theta-zero lift as a mathematical gap after an independent containing-class proof supplies that case.',
 boundary='KCUB12 allows all17 degree<=6 negative-weight atoms under the stated coefficient/collision marks. The explicit cleared target is theta^7 with51 multiplier monomials, forcing theta0 even over complex fields. KCUB11 closes that zero slice through its forced degree7 atom. The known degree7 map gives matching sharpness only in this marked class. The larger degree7 raw12-term model remains UNKNOWN. The original degree6 Groebner and real queries and direct zero lift retain their timeout statuses.',
 reopen_when='Apply parameter-field continuation to the larger degree7 model in its three new first-output parameters, retain denominator factors, and solve every exceptional branch; or change the marks and weight pattern. Do not re-run the closed degree6 class or the six numerical theta slices as new evidence.',
 claims=['KCUB12','KCUB11'],sources=[BASE+'WEIGHTED_MARKED_DEGREE.md',BASE+'EXPLORATION.md'],implementation=[BASE+n for n in ['weighted_newfirst.py','weighted_newfirst6_slice.py','weighted_newfirst6_parameter.py','weighted_newfirst6_lift.py','verify_weighted_newfirst6.py']],keywords=['Keller','degree6','sharp degree7','parameter field','theta^7','denominator exceptions','marked collision']))
for r in lessons['lessons']:
 if r['id']=='METHOD-KELLER-RAW13-TORUS-RIGIDITY-20261114':
  addition=' KCUB11 subsequently allows all six extra degree7 tail atoms with the first output fixed and still forces the known representative. KCUB12 excludes the entire17-atom degree6 marked weight model via theta^7 and the certified zero slice; sharpness is confined to that marked class.'
  if addition not in r['boundary']:r['boundary']+=addition
  r['reopen_when']='Change first-output terms in the23-atom degree7 marked model, weights or collision marks. Adding only the six allowed degree7 tail atoms is closed by KCUB11, and the degree6 marked weight sector is closed by KCUB12; global raw12-term and normalized10-term questions remain open.'
  r['claims']=list(dict.fromkeys(r['claims']+['KCUB11','KCUB12']));r['sources']=list(dict.fromkeys(r['sources']+[BASE+'WEIGHTED_TAIL_CLASSIFICATION.md',BASE+'WEIGHTED_MARKED_DEGREE.md']));r.setdefault('reviewed_claims',{}).update({k:claim_fingerprint(known[k]) for k in ['KCUB11','KCUB12']})
save('knowledge/lessons.json',lessons)
work=load('knowledge/work_items.json')
for r in work['actions']:
 if r['id']=='WORK-GLDC':
  addition=' KCUB11 closes all19 negative-weight tail atoms through degree7 with the fixed first output and nonzero-beta collision: all six additions vanish. KCUB12 gives sharp degree7 only in the marked weight-(1,-1,-2) sector; its17-atom degree6 model is excluded by an independently checked theta^7 certificate and the certified zero slice. Next continue the23-atom degree7 raw12-term model with new first-output parameters using parameter-field elimination and exact exceptional-branch checks, or change the weights/marks. Its30-second real query remains UNKNOWN. No global3D degree6 exclusion, raw13-term minimum, ordinary cubic minimum dimension or complete normalized support optimum is claimed.'
  if addition not in r['next_step']:r['next_step']+=addition
  r['claims']=list(dict.fromkeys(r['claims']+['KCUB11','KCUB12']));r['sources']=list(dict.fromkeys(r['sources']+[BASE+'WEIGHTED_TAIL_CLASSIFICATION.md',BASE+'WEIGHTED_MARKED_DEGREE.md']));r.setdefault('reviewed_claims',{}).update({k:claim_fingerprint(known[k]) for k in ['KCUB11','KCUB12']})
save('knowledge/work_items.json',work)
print('Registered KCUB11 and KCUB12: complete weighted slices only. Global simplicity goal remains open.')
