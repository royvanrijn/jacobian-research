#!/usr/bin/env python3
"""Register arbitrary two-new-form exclusion and raw13-pattern coefficient rigidity."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;BASE='jacobian/cubic-reduction/'
sys.path.insert(0,str(ROOT/'scripts'))
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
ledger=load('MATH_STATUS.json');known={e['id']:e for e in ledger['entries']}
v=load(BASE+'power-secant-two/validation.json');r=load(BASE+'raw13-rigidity/validation.json')
assert v['status']=='pass' and v['unit_identities']==607 and v['source_sha256']==sha(HERE/'verify_power_secant_two.py')
assert v['certification_sha256']==sha(HERE/'power-secant-two/certification.json')
assert r['status']=='pass' and r['identities']==15 and r['source_sha256']==sha(HERE/'verify_raw13_rigidity.py')
assert r['certification_sha256']==sha(HERE/'raw13-rigidity/certification-direct.json')
def register(identifier,title,scope,source,checker,dependencies,files,notes):
 assert identifier not in known or known[identifier]['canonical_source']==BASE+source
 e=copy.deepcopy(known['KCUB8']);e.update(id=identifier,title=title,scope=scope,canonical_source=BASE+source,checker=BASE+checker,dependencies=dependencies,
  artifact_hash='sha256:'+sha(HERE/checker),software_lock=[str((HERE/n).resolve().relative_to(ROOT)) for n in files],supersedes_notes=notes,supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
 if identifier in known:known[identifier].clear();known[identifier].update(e)
 else:ledger['entries'].append(e);known[identifier]=e
 return e
secantfiles=['verify.py','verify_power_single_compress.py','verify_power_span_optimality.py','power-span-split/validation.json','power-basis-anneal/candidate.json','power_secant_two.py','certify_power_secant_two.py','power_single_compress.py','power_line_compress.py','power-secant-two/checkpoint.json','power-secant-two/certification.json','power-secant-two/validation.json','requirements.txt','graph-requirements.txt']
primary=load(BASE+'power-secant-two/checkpoint.json')
for item in primary['cases']:
 if item['status']=='unit':
  tag='power-secant-two/pair-%02d-%02d'%tuple(item['pair']);secantfiles += [tag+'.json',tag+'-cert.json',tag+'-cert.sing']
secantfiles += ['power-single-compress/case-%02d.json'%i for i in range(15)]
e9=register('KCUB9','Two arbitrary new cubes cannot shorten the saved89 identity',
 'For the exact saved89-form padded companion identity of KCUB5, no shorter representation can retain a subset of its selected forms, add at most two arbitrary new linear cubes and freely re-solve all retained output coefficient columns, over any characteristic-zero extension. Independent new output columns force both new cubes into the saved scalar span, excluded by KCUB8. Proportional columns reduce to one scalar sum of two cubes in that span. For each of the102 missing source co-occurrence edges, its mixed second derivative is zero; independence of the two source linear forms forces each form to have clique support. The15 maximal cliques give120 unordered support pairs. Seven pairs exclude all triples by zero coordinates; the other113 produce607 unit systems covering214 proportional-column deletion triples. A separate standard-library Fraction checker reconstructs the combined membership spaces with one common projection, every coordinate map, all support-pair and coefficient-triple coverage, and all607 rational unit identities (29023 multiplier monomials, maximum degree12). No finite pool or coefficient bound limits the theorem. This is a replacement obstruction relative to one saved identity, not global rank89 or a minimum Keller dimension; at least three genuinely new forms remain possible. Primary three-form tests in the930-form pool and a tangent-cubic search at930 selected directions find no smaller witness only within their recorded scopes. The fixed-tensor unrestricted interval[34,89] and all other global simplicity questions remain open.',
 'TWO_CUBE_OBSTRUCTION.md','verify_power_secant_two.py',['KCUB8'],secantfiles,
 ['This excludes the unrestricted two-new-form escape route left by KCUB8, and strengthens the finite930-form pair search. It does not exclude three or more new forms.'])
rawfiles=['verify.py','verify_power_single_compress.py','raw13_specialization.py','raw13_normalized.py','certify_raw13_rigidity.py','parse_raw13_certificates.py',
 'raw13-specialization/input.json','raw13-specialization/checkpoint.json','raw13-normalized/input.json','raw13-normalized/checkpoint.json','raw13-rigidity/certification-direct.json','raw13-rigidity/validation.json','requirements.txt','graph-requirements.txt']
for tag in ['raw-unit-0','raw-unit-1','normalized-squares']:rawfiles += ['raw13-rigidity/'+tag+suffix for suffix in ['-parsed.json','.sing','.log']]
rawfiles.append('../de-bondt-2610.24889/verify.py')
e10=register('KCUB10','Raw13-atom collision pattern is rigid up to diagonal conjugation',
 'In the exact13-nonlinear-atom support pattern of the established three-variable degree-seven map F=(a-3a^2b-a^3c,3v^2ac+(9v^2-6v-2)b,v^3c+(6v^2+2v)b^2), v=1+2ab, consider all coefficient specializations with identity linear part, determinant one and a collision between(0,1,c0) and(-1,1,d0). Over any characteristic-zero field, all such maps are F_w=D_w F(D_w^-1 X), with w nonzero and D_w=diag(1,1,-1/w). Their points are(0,1,8/w) and(-1,1,-4/w), with common image(0,1,0). Every coefficient is nonzero, so all13 atoms are forced within this support/collision family. Two unit identities exclude zero first-output coefficients; diagonal normalization retains the transformed common second coordinate beta, and thirteen squared-coordinate ideal memberships force every normalized coefficient to the known representative with beta=1. An independent Fraction implementation reconstructs the raw and normalized universal determinant/collision equations, verifies all15 identities (94324 multiplier monomials, maximum degree12), and checks the complete Laurent family determinant and collision. Seven initial unnormalized boundary timeouts remain historical inconclusive attempts; the normalized certificate route resolves this stated family separately. The raw13-term representative is pre-existing, not a new construction, support record or normalized collision-axis bound. Other raw12-term patterns, different collision families and global sparsity remain open. The ordinary cubic12/36 and normalized three-variable[10,48] counts concern different constraints.',
 'RAW13_PATTERN_RIGIDITY.md','verify_raw13_rigidity.py',['BCRDB1'],rawfiles,
 ['Documents the pre-existing raw13-term baseline and closes coefficient deletion only in the specified pattern/collision family. No new counterexample or global term optimum is claimed.'])
save('MATH_STATUS.json',ledger)
lessons=load('knowledge/lessons.json')
def put(record):
 record['programme_status']='archived';record['reviewed_claims']={k:claim_fingerprint(known[k]) for k in record['claims']}
 old=next((x for x in lessons['lessons'] if x['id']==record['id']),None)
 if old is None:lessons['lessons'].append(record)
 else:old.clear();old.update(record)
put(dict(id='METHOD-KELLER-SECANT-CLIQUE-20261114',area='geometry',title='Use mixed second derivatives to cover arbitrary two-new-form replacements',
 when='A fixed-span lower bound forces any smaller shared-cube identity to introduce cubes outside the old span.',
 use='Separate the rank of the new output coefficient columns. Two independent columns force inside-span cubes. Proportional columns reduce to ell^3+m^3 in S; mixed second derivatives at missing co-occurrence edges force both source supports to be cliques. Use a common quotient projection across source supports, normalize the old proportional coefficient directions, and export unit-ideal identities for all support-pair/triple systems. Independently reconstruct membership row spaces and coordinate maps and multiply every certificate.',
 avoid='Adding separately normalized residual equations from different source cliques; retaining a finite-pool restriction after a continuous clique cover is available; inferring arbitrary three-form exclusion from two-form closure or finite circuit searches.',
 boundary='KCUB9 closes all at-most-two-new-form shortenings of the saved89 identity over characteristic-zero extensions through120 support pairs and607 unit identities. Three-form rank2 residue tests examine50138 pool triples, rank1 tests123379 tasks, and tangent tests930 p directions with all q solved linearly; these are primary scoped misses, not a global three-form theorem. All930 tangent directions finish, with455 nonzero tangent spaces and no four potential cancellation coordinates in one old output direction.',
 reopen_when='Introduce at least three new forms outside the saved span, enlarge the current finite pool, solve unrestricted three-form residue equations, or change the companion. The930 fixed tangent directions and recorded finite circuit classes are terminal misses and must not be rerun unchanged.',
 claims=['KCUB9','KCUB8'],sources=[BASE+'TWO_CUBE_OBSTRUCTION.md',BASE+'EXPLORATION.md'],implementation=[BASE+n for n in ['power_secant_same_support.py','power_secant_two.py','certify_power_secant_two.py','verify_power_secant_two.py','power_three_rank2.py','power_three_rank1.py','power_tangent_search.py']],keywords=['Keller','two new cubes','secant','mixed derivative','source clique','607 identities','three-form exchange']))
put(dict(id='METHOD-KELLER-RAW13-TORUS-RIGIDITY-20261114',area='geometry',title='Keep the raw13-term baseline and carry the collision coordinate through torus normalization',
 when='Searching coefficient specializations of the established degree-seven three-variable map, or comparing raw support with collision-axis-normalized support and ordinary cubic support.',
 use='Record the raw row counts2,5,6 before optimizing. Separate the support pattern and collision family. Exclude zero normalization coefficients exactly, then normalize by diagonal conjugation while retaining the transformed common b-coordinate beta. Reconstruct all normalized equations independently, including that changed collision point, and verify ideal memberships coefficient by coefficient. Parse large expanded Singular certificates directly as sparse rational monomials rather than through general expression expansion.',
 avoid='Reporting the pre-existing raw13 map as a new record; transferring normalized-axis lower bounds to raw coordinates; fixing beta=1 before proving it; calling seven unnormalized5-second timeouts UNSAT; treating the stopped slow postprocessor as an unfinished algebraic calculation.',
 boundary='KCUB10 classifies precisely the13-atom coefficient pattern with collisions(0,1,c),(-1,1,d). Two unit and thirteen squared-coordinate certificates give the one-parameter diagonal family, with all13 atoms nonzero. It is not a global raw12-term exclusion. Complete Singular logs were retained when the slow symbolic parser was stopped; direct sparse parsing and independent replay verify94324 multiplier monomials. The original degree-seven map and its13 raw occurrences are established input, separate from12-variable cubic36 and normalized three-variable[10,48].',
 reopen_when='Change the raw monomial support, collision configuration, or polynomial coordinate transformation. Do not repeat coefficient-deletion searches in the closed pattern/family. Global raw12-term existence and normalized10-term attainment remain open.',
 claims=['KCUB10','BCRDB1'],sources=[BASE+'RAW13_PATTERN_RIGIDITY.md',BASE+'EXPLORATION.md'],implementation=[BASE+n for n in ['raw13_specialization.py','raw13_normalized.py','certify_raw13_rigidity.py','parse_raw13_certificates.py','verify_raw13_rigidity.py']],keywords=['Keller','raw13 terms','coefficient pattern','torus normalization','collision beta','sparse certificates']))
for item in lessons['lessons']:
 if item['id'] in ['METHOD-KELLER-CUBE-SPAN-CLASSIFICATION-20261114','METHOD-KELLER-ONE-CUBE-COMPRESSION-20261114','METHOD-KELLER-SHARED-CUBE-BASIS-EXCHANGE-20261114']:
  addition=' KCUB9 now closes all arbitrary two-new-form shortenings, not only finite-pool pairs, using the mixed derivative clique lemma and607 independent unit identities. Construction must introduce at least three new forms. The recorded930-form three-exchange and tangent classes finish without improvement but do not exhaust arbitrary triples.'
  if addition not in item['boundary']:item['boundary']+=addition
  item['reopen_when']='Use at least three genuinely new forms outside the saved scalar span, enlarge the finite pool beyond the closed classes, solve an unrestricted three-form family, or change the companion. Arbitrary two-new-form compression is closed by KCUB9.'
  item['claims']=list(dict.fromkeys(item['claims']+['KCUB9']));item['sources']=list(dict.fromkeys(item['sources']+[BASE+'TWO_CUBE_OBSTRUCTION.md']));item.setdefault('reviewed_claims',{})['KCUB9']=claim_fingerprint(e9)
save('knowledge/lessons.json',lessons)
work=load('knowledge/work_items.json')
for item in work['actions']:
 if item['id']=='WORK-GLDC':
  addition=' KCUB9 closes arbitrary at-most-two-new-form reductions of the saved89 identity with607 independent unit certificates; next seek at least three new forms or a changed companion. The930-form three-exchange classes and930 tangent p directions finish without a smaller witness, without a global triple exclusion. KCUB10 separately classifies the established raw degree-seven13-atom pattern in its stated collision family and forces every atom nonzero; next raw sparsity work must change the pattern, collision family or polynomial coordinates. Raw13, ordinary cubic12/36 and normalized three-variable[10,48] are distinct metrics; global degree/term/special-form minima are not closed.'
  if addition not in item['next_step']:item['next_step']+=addition
  item['claims']=list(dict.fromkeys(item['claims']+['KCUB9','KCUB10']));item['sources']=list(dict.fromkeys(item['sources']+[BASE+'TWO_CUBE_OBSTRUCTION.md',BASE+'RAW13_PATTERN_RIGIDITY.md']))
  item.setdefault('reviewed_claims',{}).update({k:claim_fingerprint(known[k]) for k in ['KCUB9','KCUB10']})
save('knowledge/work_items.json',work)
print('Registered KCUB9 and KCUB10. No smaller construction or global minimum claimed; full goal remains open.')
