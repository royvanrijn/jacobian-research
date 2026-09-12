#!/usr/bin/env python3
"""Deterministic final audit/pilot report, without promoting deficiencies to g."""
import argparse
from pathlib import Path

import audit_wide_arithmetic_censoring as audit
import run_two_class_relation_pilot as pilot
import wide_arithmetic_profile_core as storage


def build(out):
    plan=pilot.check(out,complete=True)
    for name,text in audit.build(out).items():
        pilot.require((out/name).read_bytes()==text.encode(),'censoring audit changed')
    census=pilot.read(out/'censoring_audit.json')
    replay=pilot.read(out/'REPLAY.json')
    pilot.require(replay['plan_sha256']==pilot.sha(out/'plan.json'),'replay plan mismatch')
    pilot.require(replay['checker_sha256']==pilot.sha(Path(__file__).with_name('verify_two_class_relation_pilot.sage')),'replay checker changed')
    rows=[]
    for row,checked in zip(plan['rows'],replay['rows']):
        key=row['curve_key']
        pilot.require(key==checked['curve_key'],'replay curve mismatch')
        r=dict(row)
        supervisor=pilot.read(out/key/'supervisor.json')
        r['worker_status']=supervisor['status']
        if supervisor['status']=='PASS_BOUNDED_WORKER':
            pilot.require(checked['status']=='PASS_EXACT_BOUNDED_RELATION_REPLAY','missing exact replay')
            for name,digest in checked['bindings'].items():
                pilot.require(pilot.sha(out/key/name)==digest,'replay artifact changed')
            result=pilot.read(out/key/'result.json')
            setup=pilot.read(out/key/'setup.json')
            r.update(candidate_count=result['candidate_prefixes'][-1]['attempted'],
                     noncanonical_relations=result['noncanonical_relation_count'],
                     matrix=result['final_matrix'],baseline=setup['baseline'],
                     measurement_gate=result['measurement_gate'],
                     median_norm_numerator_bits=storage.median([int(k) for k,v in result['norm_numerator_bit_histogram'].items() for _ in range(v)]),
                     wall_seconds=result['wall_seconds'],
                     reduced_cubic_ascending=setup['form']['reduced_cubic_ascending'],
                     bach_GRH_bound_interval=setup['bach_GRH_bound_interval'],
                     factor_base_generation=result['factor_base_generation'])
        r['global_class_2rank_estimate']=None
        r['global_class_2rank_upper_bound']=None
        rows.append(r)
    pilot.require(len(rows)==len(plan['rows'])==len(replay['rows'])==2,'pilot cardinality mismatch')
    result={'schema':'elliptic-curves.two-class-pilot-report.v1',
            'status':'PASS_AUDIT_AND_BOUNDED_PILOT_NO_CLASS_RANK_RESULT','rows':rows,
            'point_searches_launched':0,'full_class_group_calls':0,
            'old_evidence_unchanged':True,'historical_excluded_from_prospective_audit':True,
            'plan_sha256':pilot.sha(out/'plan.json'),'replay_sha256':pilot.sha(out/'REPLAY.json'),
            'audit_sha256':pilot.sha(out/'censoring_audit.json'),'reporter_sha256':pilot.sha(Path(__file__)),
            'boundary':'No independently measured g or validated predictor. Matrix deficiencies concern the chosen factor-base image, not the full class group. Known-rank-derived g lower bounds are calibration metadata only.'}
    lines=['# Censoring audit and two-class relation pilot', '',
           'Status: completed bounded experiment; **no class-group 2-rank estimate or upper bound obtained**. '
           'No point searches or full class-group computations. Earlier frozen arithmetic evidence remains unchanged.', '',
           '## 1. Censoring and ramification', '',
           'The [full audit](CENSORING.md) uses only the original 2,080 prospective rows. '
           'This is retrospective analysis of that frozen cohort, not new prospective validation. Historical rows are excluded.', '',
           'LOCAL completion is 235/416 (56.5%) in the smallest BASE-discriminant quintile and '
           '127/416 (30.5%) in the largest. Parent composition matters too: the class1 parent has '
           'only 3/80 completed LOCAL rows. The available arithmetic is not a uniform subsample.', '',
           'Among completed cases, ramification-count versus rank-lower-bound Spearman is '
           f"{census['completed_case_ramification_association']['spearman_ramification_vs_rank_LB']:.4f} overall (898 cases) and "
           f"{census['completed_11952_ramification_association']['spearman_ramification_vs_rank_LB']:.4f} within11952 (145 cases). "
           'The within-family-and-size centered rank association is -0.0544 (883 informative cases); it is descriptive, not a causal adjustment.', '',
           'The >=23 tail fractions by known ramification count are 5/168 (0–4 primes), '
           '7/312 (5–6), 3/258 (7–8), and2/160 (9+). Within11952, none of the17 completed '
           '9+-prime fibres is in that tail. Thus the historical ramification excess does not reproduce '
           'as positive enrichment in the observed broad population. Missing-data worst-case intervals '
           'remain wide; this is not a disproof of an underlying relation.', '',
           '## 2. Frozen two-field relation prototype', '',
           'The two fields are user-selected retrospective engineering controls. Each receives the same '
           '20,000 rational-prime bound, mandatory denominator-prime support,16,384 primitive binary-form '
           'candidates, and180-second hard wall limit. Rank labels and points are excluded from worker inputs. '
           'All rational principal-ideal relations (p) are inserted before measuring new information.', '',
           'The maximal cubic orders are certified; exact multiplication tables yield binary cubic forms '
           'of discriminant D_K and explicitly transported generators. A fixed finite height descent is used, '
           '**not** production Julia reduction. This implementation enumerates smooth-norm candidates; '
           'it is not KSW’s optimized NFS sieve or large-prime relation collector.', '',
           audit.text_table(['11952 fibre','Known rank LB','Columns','Canonical rows','New smooth relations','Final mod-2 deficiency','Rank gain beyond canonical'],
                [[r['t'],r['final_rank_lower_bound'],r.get('matrix',{}).get('columns'),r.get('baseline',{}).get('rows'),
                  r.get('noncanonical_relations'),r.get('matrix',{}).get('deficiency'),r.get('matrix',{}).get('rank_gain_beyond_canonical')] for r in rows]), '',
           'Both fields complete all16,384 candidates. The historical fibre’s sole new relation is the '
           'selected field generator itself, (a,b)=(0,1), at attempt2. No later prefix adds a relation. '
           'Deficiencies1,905 and1,890 are dominated by missing principal relations and factor-base splitting, '
           'not measurements near the desired class ranks. Fixed-N noncanonical-relation comparisons were '
           'not achieved. The small control x^3−x−1 admitted33 noncanonical relations in64 candidates and '
           'closed its displayed factor-base quotient, validating the basic relation path but not high-discriminant coverage.', '',
           '## Interpretation gate', '',
           'Let H_F be the subspace of Cl(K)/2Cl(K) generated by the selected factor-base ideals. '
           'Every verified principal relation gives a row, hence a surjection from the presented '
           'mod-2 quotient onto H_F. Consequently its deficiency bounds dim(H_F) **above**. '
           'It is not a lower bound or an approximation to g. Without factor-base generation, '
           'it is not an upper bound on the full g either.', '',
           '[Klagsbrun–Sherman–Weigandt, §§4–5](https://arxiv.org/abs/1606.07178) '
           'use mod-2 relations to avoid full integral class-group structure, but also supply '
           'a GRH-dependent factor-base generation argument for their global upper bounds. '
           'Their production relation collection is substantially larger than this pilot. '
           'No generation argument has been certified here; the Bach bounds are approximately '
           '721,495 and1,027,301, whereas the small-prime bound in this pilot is20,000.', '',
           'Known rational-point lower bounds and the existing Brumer–Kramer local terms still imply '
           '**g>=16** at921/653 and **g>=18** at110314/102227. These are consequences of known rank, '
           'not independent discoveries or estimates g≈16 andg≈20. Larger class rank alone does not '
           'force higher rational rank or prove that the Selmer upper bound is attained.', '',
           'Next gate, before any broader relation-deficiency panel: improve polynomial/skew selection '
           'and sieving enough to obtain substantial verified noncanonical relation coverage on these '
           'same controls. No extension is launched or scheduled by this report.', '',
           '## Reproduction and provenance', '',
           '[Plan](plan.json), [machine report](REPORT.json), and [exact replay](REPLAY.json). '
           'Every candidate norm/cofactor checksum, accepted ideal identity, and final GF(2) rank replays. '
           'Sage supplies a separate matrix-rank implementation; number-field arithmetic shares Sage/PARI.', '',
           'The live inventory was expanded elsewhere during this turn. The original445-row inventory '
           'is resolved from Git commit1959f550ca43e5ff4492e8b5eb383f3104562a1a, '
           'with exactly the SHA256 pinned in the historical plan. No live inventory or old frozen '
           'hash was altered or silently rebound.', '',
           '```sh',
           'python3 research/elliptic-curves/cas/audit_wide_arithmetic_censoring.py check',
           'python3 research/elliptic-curves/cas/run_two_class_relation_pilot.py check',
           'timeout 120 sage -python research/elliptic-curves/cas/verify_two_class_relation_pilot.sage',
           'python3 research/elliptic-curves/cas/report_two_class_relation_pilot.py check',
           'PYTHONPATH=research/elliptic-curves/cas python3 -m unittest discover -s research/elliptic-curves/tests -p test_two_class_relation_pilot.py -v',
           '```','']
    return {'REPORT.json':storage.stable_json(result),'SUMMARY.md':'\n'.join(lines)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command',choices=['report','check'])
    ap.add_argument('--output',type=Path,default=pilot.OUT)
    args=ap.parse_args()
    outputs=build(args.output)
    for name,text in outputs.items():
        p=args.output/name
        if args.command=='check' or p.exists():
            pilot.require(p.read_bytes()==text.encode(),'report replay mismatch: '+name)
        else:
            p.write_bytes(text.encode())
    print('MOD2_REPORT|PASS|no_class_rank_result|original_cohorts_unchanged')


if __name__=='__main__':
    main()
