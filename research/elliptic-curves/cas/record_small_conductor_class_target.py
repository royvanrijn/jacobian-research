#!/usr/bin/env python3
"""Render audited norm waves, point relations and the exact lower/anchor certificates."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
NOTE=ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_TARGET_2026-09-06.md'
ID='EC-SMALL-CONDUCTOR-CLASS-TARGET-20260906'
LOW_ID='EC-SMALL-CONDUCTOR-CLASS-LOWER16-20260906'
CHAR_ID='EC-SMALL-CONDUCTOR-CLASS-CHARACTERS-20260906'
KINDS=[('box','', 'small-conductor-class-target-v1',''),
       ('strip','strip_', 'small-conductor-class-target-strips-v1','_strips'),
       ('protected','protected_', 'small-conductor-class-target-protected-v1','_protected'),
       ('capped','capped_', 'small-conductor-class-target-capped-v1','_capped'),
       ('residual','residual_', 'small-conductor-class-target-residual-v1','_residual')]


def read(p):return json.loads(p.read_text())
def rel(p):return str(p.relative_to(ROOT))
def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def evidence_link(p):return '../../'+rel(p)


def main():
    waves=[]
    for family,prefix,directory,suffix in KINDS:
        paths=sorted(ART.glob('small_conductor_class_target_'+prefix+'wave_*_v1.json'))
        for index,p in enumerate(paths,1):
            r=read(p)
            assert r['wave']==index and r['status']=='PASS'
            for name,h in r['sources'].items():assert hashed(ROOT/name)==h,name
            waves.append((family,ROOT/'artifacts/local/elliptic-curves'/directory,suffix,p,r))
    if not waves:raise ValueError('no audited waves')
    family,directory,suffix,last_path,last=waves[-1]
    runner=ROOT/('elliptic-curves/cas/pursue_small_conductor_class_target'+suffix+'.sage')
    combined_path=ART/('small_conductor_class_target_combined_%s_%03d_v1.json'%(family,last['wave']))
    combined=read(combined_path) if combined_path.exists() else None
    if combined:
        assert combined['status']=='PASS'
        for name,h in combined['sources'].items():assert hashed(ROOT/name)==h,name
    dim=combined['matrix']['quotient_dimension'] if combined else last['matrix']['quotient_dimension']
    gain=combined['independent_gain_from_point_relations'] if combined else 0
    seeded_point_gain=0
    if family in ['capped','residual']:
        first_cap=next(r for f,d,s,p,r in waves if f=='capped')
        seed_protected=next(r for f,d,s,p,r in waves if f=='protected' and r['wave']==3)
        seeded_point_gain=seed_protected['matrix']['quotient_dimension']-first_cap['starting_dimension']
    checker=ROOT/'elliptic-curves/cas/certify_small_conductor_class_target_combined.sage' if combined else runner
    command='sage -python '+rel(checker)+(' --family '+family+' --wave '+str(last['wave'])+' --check' if combined else ' check --wave '+str(last['wave']))
    formal_dim=dim
    completion_path=ART/'small_conductor_class_completion_v1.json'
    completion=read(completion_path) if completion_path.exists() else None
    completion_proof=ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md'
    if completion:
        assert completion['status']=='PASS' and family=='residual' and last['wave']==completion['residual_wave']
        assert completion['formal_relation_matrix']==last['matrix']
        for name,h in completion['sources'].items():assert hashed(ROOT/name)==h,name
        dim=completion['conditional_on_grh_class_two_rank_upper_bound']
        checker=ROOT/'elliptic-curves/cas/certify_small_conductor_class_completion.sage'
        command='sage -python '+rel(checker)+' check'
    lower_path=ART/'small_conductor_class_lower16_v1.json'
    char_path=ART/'small_conductor_class_characters_v1.json'
    lower=read(lower_path);chars=read(char_path) if char_path.exists() else None
    assert lower['status']=='PASS' and lower['unconditional_class_two_rank_lower_bound']==16
    assert dim>=16,'upper bound contradicts certified class lower bound16'
    portable=ART/'small_conductor_class_target_portable_replay_v1.json'
    replay=read(portable) if portable.exists() else None
    isolated=bool(replay and replay['status']=='PASS' and replay.get('family')==family and replay.get('last_wave')==last['wave'])
    if completion:isolated=bool(isolated and replay.get('completion_certificate_sha256')==hashed(completion_path))
    exact=dim==16
    lines=['# MW16 at 3/17: '+('class-2-rank 16 and curve rank 22 under GRH' if exact else 'reducing the class-group quotient toward 16'),'',
        'Mathematical authority: `MATH_STATUS.json`, entries `'+ID+'`, `'+LOW_ID+'` and `'+CHAR_ID+'`.','',
        '**Current class-2-rank bounds: 16 <= g <= '+str(dim)+'.** The lower bound is unconditional;',
        'the upper bound depends on GRH for the stated ideal-class characters.',
        'The unconditional curve-rank lower bound is **22**. '+('Matching bounds prove **exact rank 22 under GRH**.' if exact else 'Exact curve rank remains unknown.'),'',
        'The [preceding descent study](SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md) establishes',
        'the cubic field, local correction `rank <= dim Sel_2 <= g+7`, proved even',
        'Selmer parity and the interval-certified generating cutoff 37,638 under GRH.',
        'Thus `g<=16` suffices for rank 22 under the same assumption.','',
        *(['The [completion proof](SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md)',
           'closes the bound using a quadratic-character explicit formula. The formal',
           'relation quotient remains 18; exact memberships of small prime classes in',
           'the sixteen-anchor span suffice to prove that span generates under GRH.',
           'At cutoff 50,000 the worst-case corrected margin exceeds 17.16. A second',
           'interval calculation with a weaker archimedean bound is also positive.',
           'The [certificate]('+evidence_link(completion_path)+') lists every membership',
           'and unresolved prime. No twenty-third point or missing norm relation is claimed.',''] if completion else []),
        '## Unconditional lower bound and independent ideal anchors','',
        'The [lower-bound certificate]('+evidence_link(lower_path)+') uses the known points',
        'to form `beta_i=4*x(P_i)-theta`, with norm `(8*y(P_i)+4*x(P_i))^2`.',
        'Residue characters at a fixed set of 128 rational primes prove independence',
        'of all 22 field square classes. Their valuation-parity matrix at every bad',
        'prime has rank 4. Its kernel supplies 18 independent products with even',
        'valuations everywhere: away from `2*Delta(E)`, separability and the square',
        'norm prove the assertion, including when x has a pole. Adjoining -1 adds',
        'an independent class because its cubic norm is negative. The field Selmer',
        'group has dimension at least 19; the unit square-class dimension is 3.',
        'The [field Selmer exact sequence](https://arxiv.org/html/1606.07178#S4.SS3)',
        'therefore proves **g >= 19-3 = 16** unconditionally.','',
        'The [character certificate]('+evidence_link(char_path)+') imposes dyadic unit',
        'Hilbert symbols and positive signs at the three real embeddings on these',
        '19 classes. The constraint rank is 3, leaving 16 independent everywhere-',
        'unramified quadratic extensions. Each dyadic residue field is F2; the',
        'units `1+pi^k`, `1<=k<=2e`, generate units modulo squares because successive',
        'unit-filtration quotients are F2 and `U^(2e+1)` consists of squares by Hensel.',
        'Class field theory gives 16 ordinary ideal-class characters. Their values',
        'on the listed 16 prime ideals form an invertible matrix, proving those',
        'ideal classes independent. Protected waves keep these anchors free during',
        'elimination and never target them. A purported relation supported only',
        'on the anchors causes an explicit failure.','',
        '## Authorized relation continuation','',
        'The user set the goal of reducing the quotient to 16. Each wave freezes a',
        'finite target list, candidate region, smoothness bound, one-worker time',
        'and memory limits, source hashes and checkpoint policy before execution.',
        'Targets are free prime-ideal columns after supported row reduction.',
        'Square boxes use exact Hessian-reduced index-prime lattices. Near-root',
        'strips use the three real norm roots transformed into those lattices:',
        'exact algebraic arithmetic floors each slope times `2^96`, and each',
        'positive denominator v tests the integer center and its two neighbors.',
        'Primitive-coordinate filters and previous search regions are replayed.',
        'Protected waves also skip targets that have become elimination pivots;',
        'the checker reconstructs every such adaptive decision. Increasing the',
        'smoothness cutoff permits reexamining a previously searched region.','',
        'Every retained witness has its norm and principal ideal checked exactly,',
        'including the nonmonic norm identity\'s fixed square factor. All prime',
        'coordinates above 37,638 remain until exact elimination cancels them.',
        'Supported rank is independently checked as `rank(all)-rank(outside)`.',
        'Rejected candidates are unnecessary for a rank upper bound; counts and',
        'digests describe the search, without an exhaustive-miss assertion.','',
        'Residual-representative waves also allow a pivot prime ideal, or an outside-',
        'base ideal whose outside normal form vanishes, when its reduced',
        'normal form contains unresolved nonanchor coordinates. A greedy spanning',
        'set is selected in ascending prime order, then filled with inexpensive',
        'eligible ideals. A target is skipped only when its normal form lies',
        'entirely in the certified anchor span. These normal forms guide selection;',
        'they are not treated as newly proved class-group relations.', '',
        'The capped implementation bounds every product-tree node by `M+1`, where M',
        'is the primorial. Capping commutes with positive multiplication, and a',
        'remainder at most M is unchanged modulo any product above M. Thus every',
        'leaf remainder stays exact. The [fixed benchmark](../../artifacts/generated-results/elliptic-curves/small_conductor_capped_remainders_benchmark_v1.json)',
        'checks agreement against both the full tree and scalar division on 58,819',
        'actual values plus 200 edge-case lists. Timing is specific to that target.', '',
        '## Audited waves','',
        '| Wave | Targets checkpointed | Region | Retention cutoff | Candidate occurrences | Relation occurrences | Independent supported gain | Norm-matrix dimension | Worker seconds |',
        '| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for f,d,s,p,r in waves:
        protocol=read(d/('wave_%03d'%r['wave'])/'protocol.json')
        seconds=sum(read(ROOT/c['path'])['wall_seconds'] for c in r['chunks'])
        region='box '+str(protocol['box']) if f=='box' else 'v <= '+str(protocol['vmax'])
        values=[len(r['chunks']),region,protocol['smooth_bound'],r['candidate_occurrences'],r['relation_occurrences'],r['independent_supported_gain'],r['matrix']['quotient_dimension'],round(seconds,2)]
        lines.append('| ['+f+' '+str(r['wave'])+']('+evidence_link(p)+') | '+' | '.join(map(str,values))+' |')
    if combined:
        lines+=['','The [combined certificate]('+evidence_link(combined_path)+') also appends the',
            'known points\' principal parity relations. Their odd valuations occur',
            'only at the explicitly checked bad primes; all omitted valuations',
            'are proved even. These contribute **'+str(gain)+' additional independent rows**,',
            'reducing the final norm-matrix dimension from **'+str(last['matrix']['quotient_dimension'])+' to '+str(dim)+'**.']
    if seeded_point_gain:
        lines+=['','The capped phase seeds the known-point principal parity relations directly.',
                'They supply **'+str(seeded_point_gain)+' independent rows** beyond its norm-only seed; this gain',
                'is already included in the displayed capped-phase dimensions.']
    lines+=['','From the starting bound 2,879, the formal relation quotient has fallen to **'+str(formal_dim)+'**.',
        ('The explicit-formula completion proves that the sixteen independent anchor classes generate, reaching **16 under GRH**.' if completion else '**'+str(dim-16)+' further independent supported rows** would reach 16.'),
        'The curve-rank upper bound is **'+str(2*((dim+7)//2))+' under GRH**.',
        'These are ideal relations, not new rational points or an algebraic-rank parity claim.','',
        '## Replay','','```bash',command,'```','',
        'The checker replays the inherited curve and field proofs, principal-ideal',
        'witnesses, target selection, matrix transitions and applicable point relations.',
        *(['It also verifies all prime memberships in the anchor span with a second',
           'elimination order and checks both positive explicit-formula margins.'] if completion else []),
        'Raw protocols and checkpoints are under the matching `small-conductor-class-target*`',
        'directories in `artifacts/local/elliptic-curves/`. Each generated certificate',
        'pins its exact sources and witness chunks.',
        ('The [portable replay]('+evidence_link(portable)+') passes from a fresh extracted directory.' if isolated else 'A new isolated portable replay is not yet claimed.'),'']
    NOTE.write_text('\n'.join(lines))
    status=ROOT/'MATH_STATUS.json';original=status.read_text();data=read(status)
    prior=next(e for e in data['entries'] if e['id']=='EC-SMALL-CONDUCTOR-SMALL-BASE-TARGETS-20260906')
    def make(id,title,scope,source,locks,deps):
        e=json.loads(json.dumps(prior));e.update({'id':id,'title':title,'scope':scope,
            'canonical_source':rel(NOTE),'checker':rel(source),'dependencies':deps,
            'independent_replay':isolated,'artifact_hash':'sha256:'+hashed(source),
            'software_lock':list(dict.fromkeys(rel(p) for p in [source,NOTE,*locks]))})
        return e
    lower_source=ROOT/'elliptic-curves/cas/certify_small_conductor_class_lower16.sage'
    lower_entry=make(LOW_ID,'MW16 at 3/17: unconditional class-2-rank lower bound 16',
        'A fixed128-prime residue-character matrix proves22 independent field square classes from known points. The exact bad-prime valuation-parity matrix has rank4, leaving18 independent everywhere-even-valuation products; separability and square norms certify all omitted places. Adding -1 gives19 field Selmer classes, and the unit square-class dimension3 proves class-2-rank at least16 unconditionally. No upper bound or exact curve rank follows.',lower_source,[lower_path],[prior['id']])
    char_source=ROOT/'elliptic-curves/cas/certify_small_conductor_class_characters.sage'
    char_entry=make(CHAR_ID,'MW16 at 3/17: sixteen certified unramified characters and ideal anchors',
        'On the19 certified field Selmer classes, dyadic unit Hilbert symbols and real sign conditions have rank3. Unit filtration and Hensel certify completeness of the local unit tests. The kernel gives16 independent everywhere-unramified quadratic extensions. By ordinary class field theory their residue characters form an invertible matrix on16 listed degree-one prime ideals. These ideal classes are independent unconditionally and may be protected during row elimination. This does not prove an upper bound.',char_source,[char_path,lower_path],[LOW_ID])
    locks=[Path(__file__).resolve(),runner,*[p for f,d,s,p,r in waves]]
    if combined:locks += [combined_path]
    if family in ['capped','residual']:locks += [ROOT/'elliptic-curves/cas/capped_primorial_remainders.py',ART/'small_conductor_capped_remainders_benchmark_v1.json']
    if completion:locks += [completion_path,completion_proof,ROOT/'artifacts/local/elliptic-curves/small-conductor-class-completion-v1/protocol.json']
    if isolated:locks += [portable,ART/'small_conductor_class_target_evidence_v1.json',ART/'small_conductor_class_target_evidence_v1.zip']
    scope=str(len(waves))+' bounded, checkpointed norm waves plus '+str(gain+seeded_point_gain)+' additional independent point-derived parity rows reduce the formal class quotient from2879 to'+str(formal_dim)+'. All retained principal ideals, source bindings, deterministic selections, adaptive skips and supported-intersection rank identities replay exactly. Outside prime coordinates are never discarded without cancellation. '+('A quadratic-character explicit-formula test at50000, using exact memberships in the16-anchor span and worst-case signs at all unresolved primes, has corrected margin>17.16. A second conservative interval test also passes, proving that the16 anchors generate Cl(K)/2 under GRH. The formal quotient remains18; no missing principal relations are asserted. ' if completion else '')+'The unconditional class-2-rank lower bound is16; the GRH-conditional upper bound is'+str(dim)+', yielding curve-rank upper bound'+str(2*((dim+7)//2))+'. '+('Matching bounds prove class-2-rank16,2-Selmer dimension22 and curve rank22 under GRH for the nontrivial quadratic ordinary ideal-class characters (only those trivial on the anchor span are needed).' if exact else 'Exact curve rank remains unknown.')+' No new point or unconditional upper bound is asserted. '+('Isolated portable replay passes.' if isolated else 'Local replay passes; isolated portable replay is not yet claimed.')
    entry=make(ID,'MW16 at 3/17: audited class-rank interval 16 to '+str(dim),scope,checker,locks,[prior['id'],LOW_ID,CHAR_ID])
    if completion:
        entry['title']='MW16 at 3/17: class-2-rank 16 and exact curve rank 22 under GRH'
        entry['canonical_source']=rel(completion_proof)
    for e in [lower_entry,char_entry,entry]:
        indices=[i for i,z in enumerate(data['entries']) if z['id']==e['id']]
        if indices:data['entries'][indices[0]]=e
        else:data['entries'].append(e)
    assert status.read_text()==original,'concurrent status update'
    status.write_text(json.dumps(data,indent=2)+'\n')
    print('Recorded',len(waves),'waves; class interval16 to',dim,'; point gain',gain+seeded_point_gain)


if __name__=='__main__':main()
