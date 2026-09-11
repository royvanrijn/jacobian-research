#!/usr/bin/env python3
"""Render verified event/subset evidence, basis sensitivity, and claim limits."""
import argparse
import gzip
from hashlib import sha256
import json
from math import log2
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elliptic-curves/rank_accessibility_subsets_v1'


def load(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)


def digest(path):return sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    verification=load(OUT/'verified.json');protocol=load(OUT/'protocol.json');events=load(OUT/'events.json')['events']
    assert verification['status']=='PASS'
    inputs={r['id']:r for r in load(OUT/'inputs.json.gz')['curves']}
    summaries={};envelopes={}
    for record in verification['results']:
        path=OUT/record['file'];assert digest(path)==record['sha256']
        data=load(path);s=data['summary'];key=(s['curve'],s['policy']);summaries[key]=s
        envelopes[key]=record['minimal_support_envelopes']
    lines=[
        '# Accessibility gains: exact events, all subsets, and basis sensitivity',
        '',
        'The original-basis contrast survives removal of acquisition order, but the large contrast in relative gains is sensitive to modest generic-basis changes. The old zero on 11952 was partly order-dependent. These are sparse effects in a fixed coordinate atlas, not evidence of an intrinsic avalanche or of a rank upper bound.',
        '',
        'All labels below are one-based. M1,...,M17 are the original specialized generic basis, and E1,... are the original displayed exceptional directions. Full integer words in both the public D basis and the (M,E) basis, exact points, and N/F values are retained in [events.json](events.json). Curve302 E1,...,E14 correspond to public points 1,3,4,6,9,13,14,17,19,22,26,29,30,31; 11952 E1,...,E8 are D18,...,D25.',
        '',
        '## Four successful events',
        '',
        'Every contribution is before minus after, in bits; positive values help the drop. The residual column is one quarter of the rational x-height on the frozen short Weierstrass model. The identity used is k = h_x/4 + a + c, with c = log2(gcd(N,F))/4 and a = k - log2(max(|N|,|F|))/4. All N/F pairs are fixed primitive integral homogeneous quartics before target evaluation.',
        '',
        '| Hidden target | Added direction | Q before | Q after | Total drop | Residual x-height term | Archimedean term a | Finite term c |',
        '|---|---|---|---|---:|---:|---:|---:|',
    ]
    for e in events:
        c=e['drop_components_bits']
        lines.append(f"| E{e['target']+1} | E{e['added']+1} | {e['before']['label']} | {e['after']['label']} | {e['drop_bits']:.6f} | {c['quarter_hx']:+.6f} | {c['a']:+.6f} | {c['c']:+.6f} |")
    lines += [
        '',
        '- E10 after E1: the residual x-height and archimedean terms worsen. A 24.028-bit fall in c more than offsets them. The exact parameter-height ratio is about 229,097; this is not a runtime measurement.',
        '- E14 after E1: the main favorable term is the 12.385-bit decrease in residual x-height/4. The finite term contributes another 5.366 bits, while the archimedean term offsets 3.983 bits.',
        '- E7 after E2: the 4.953-bit residual x-height improvement is partly offset by both remaining terms.',
        '- E12 after E6: residual x-height/4 changes by only 0.342 bits. The archimedean term improves by 22.901 bits, while c increases and offsets 16.470 bits. The new chart is Q=M17+E6.',
        '',
        'The signs matter: a larger c is not automatically favorable. These accounts are exact coordinate-height identities, not fitted explanations or canonical-height causal claims. The verifier checks the multiplicative identity H_before^4/H_after^4 = residual_ratio * archimedean_ratio * finite_ratio and checks x(2P-Q)=N/F using a separate rational group law.',
        '',
        '## Complete subset landscape in the original basis',
        '',
        '| Metric | Curve302 | 11952 |',
        '|---|---:|---:|',
    ]
    a,b=[summaries[c,'original'] for c in ('302','11952')]
    rows=[('Known exceptional directions',a['known_exceptional_dimensions'],b['known_exceptional_dimensions']),
          ('Held-out target/subset cells',f"{a['cells']:,}",f"{b['cells']:,}"),
          ('Directed subset edges',f"{a['edges']['count']:,}",f"{b['edges']['count']:,}"),
          ('Targets helped somewhere',f"{a['targets_with_any_native_gain']}/14",f"{b['targets_with_any_native_gain']}/8"),
          ('Largest single-addition drop (bits)',f"{a['edges']['max_bits']:.6f}",f"{b['edges']['max_bits']:.6f}"),
          ('Positive-G cells',f"{a['G']['positive_count']:,} ({100*a['G']['positive_count']/a['cells']:.3f}%)",f"{b['G']['positive_count']:,} ({100*b['G']['positive_count']/b['cells']:.3f}%)"),
          ('Median G over cells (bits)',f"{a['G']['median_bits']:.6f}",f"{b['G']['median_bits']:.6f}"),
          ('Mean positive part of G (bits/cell)',f"{a['G']['mean_positive_part_bits']:.6f}",f"{b['G']['mean_positive_part_bits']:.6f}"),
          ('Positive edges',f"{a['edges']['positive_count']:,}/{a['edges']['count']:,}",f"{b['edges']['positive_count']:,}/{b['edges']['count']:,}"),
          ('Positive-edge AUC / all edges (bits)',f"{a['edges']['mean_positive_part_bits']:.6f}",f"{b['edges']['mean_positive_part_bits']:.6f}")]
    lines += [f'| {label} | {v} | {w} |' for label,v,w in rows]
    lines += ['',
        'Each target has every subset of the other exceptional directions represented. At subset size k the control has exactly 2(17+k)^2 generic-only charts. Full discrete survival curves, their exact rational height ratios, and all winner indices are in the compressed per-policy files. Cell and edge means have explicitly different denominators; neither is an independent-sample statistic.',
        '',
        'The original-basis functions reduce to the following minimal offers. For any target, subtract the largest listed improvement whose required subset is supplied from its initial kappa; use zero when no offer is available. Every unlisted target stays at its initial kappa for every eligible subset. The independent verifier checks this compact description against every stored cell.',
        '',
        '| Fibre | Target | Required exceptional direction(s) | Winning anchor | Improvement over M17 (bits) |',
        '|---|---|---|---|---:|']
    for cid in ('302','11952'):
        for envelope in envelopes[cid,'original']:
            for offer in envelope['offers']:
                supplied=','.join(f'E{i+1}' for i in offer['required_exceptional_indices'])
                lines.append(f"| {cid} | E{envelope['target']+1} | {supplied} | {offer['anchor_label']} | {offer['improvement_over_initial_bits']:.6f} |")
    lines += [
        '',
        '11952 E2 gains 1.654608 bits from Q=E4. The old ordering supplied E2 before E4 and could not expose this held-out improvement. On 302, E2 also offers E10 a 3.849481-bit improvement, which the old path masked by supplying the better anchor -E1 first. All six minimal offers across both fibres have singleton exceptional support. Thus no pair of exceptional directions supplies an additional gain beyond these offers in the original full signed-unit/pair atlas.',
        '',
        '## Frozen basis sensitivity',
        '',
        'The five perturbations use disjoint consecutive pairs (a,b) -> (a+s*b,b), with s=+1 or -1. They are unimodular, preserve M17, and never mix generic and exceptional blocks. Generic transformations are shared across targets. Exceptional transformations pair the other exceptional indices separately for each held-out target, keeping that original target fixed and excluding its coordinates from every supplied generator. These are leave-one-target-out basis sensitivity probes, not a single simultaneous rebased acquisition path. All policies were pinned before their measurements.',
        '',
        '| Basis policy | 302 targets helped | 11952 targets helped | 302 largest drop / max G | 11952 largest drop / max G | 302 mean positive G | 11952 mean positive G |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for policy,_,_ in protocol['policies']:
        a,b=[summaries[c,policy] for c in ('302','11952')]
        lines.append(f"| {policy} | {a['targets_with_any_native_gain']}/14 | {b['targets_with_any_native_gain']}/8 | {a['edges']['max_bits']:.3f} / {a['G']['max_bits']:.3f} | {b['edges']['max_bits']:.3f} / {b['G']['max_bits']:.3f} | {a['G']['mean_positive_part_bits']:.3f} | {b['G']['mean_positive_part_bits']:.3f} |")
    lines += [
        '',
        'Relative drops are sensitive to a worse initial bank. For 11952 E2 under generic_plus, initial kappa increases from 72.929942 to 87.655814 bits, while the best held-out chart remains Q=E4 at 71.275334 bits. Its new 16.380480-bit drop is the original 1.654608-bit effect plus a 14.725872-bit degradation of the starting bank. The endpoint has not improved.',
        '',
        'More generally, the four helped 11952 targets under generic_plus end at equal or worse absolute costs than under the original full held-out bank. The rebase changes which good generic charts are available; exceptional anchors can restore accessibility relative to that altered baseline. The original-basis advantage therefore does not establish intrinsic exceptional subgroup geometry. Exceptional-block shears alone preserve the strongest original gains but alter how many subsets make the required anchor available.',
        '',
        'The supported conclusion is a small collection of useful anchor relations with presentation-dependent accessibility gains. Acquisition order explains the old 11952 zero; generic-basis choice explains much of the dramatic cross-fibre contrast in relative drops. The fixed raw-coordinate metric and one prescribed generic control policy remain limits. These results give no information about missing rank on 11952.',
        '',
        '## Reproduction and verification',
        '',
        'Run from the repository root:',
        '',
        '```bash',
        'sage -python research/elliptic-curves/cas/rank_accessibility_subsets.py prepare',
        'sage -python research/elliptic-curves/cas/rank_accessibility_subsets.py original',
        'sage -python research/elliptic-curves/cas/rank_accessibility_subsets.py all',
        'python3 research/elliptic-curves/cas/verify_rank_accessibility_subsets.py --check',
        'python3 research/elliptic-curves/cas/report_rank_accessibility_subsets.py --check',
        '```',
        '',
        'The producer resumes immutable hash-bound checkpoints. The separate verifier uses Python Fraction arithmetic for group operations and checks every distinct rational slope, every subset minimum via direct sparse support predicates, generic chart-prefix minima, and complete edge/G survival counts. It checks the four original exact event identities separately. Logs are descriptive float displays; all cost comparisons use exact integers.',
    ]
    measured=sum(r['new_rational_slopes_checked'] for r in verification['results'])
    cells=sum(r['cells_checked'] for r in verification['results'])
    edges=sum(r['edges_checked'] for r in verification['results'])
    lines.append(f'Across six policies on two fibres, verification replays {measured:,} distinct target–anchor slopes, {cells:,} held-out cells, and {edges:,} directed edges. Repeated states are not treated as independent observations.')
    lines.append('')
    basis_differences=[]
    for cid in ('302','11952'):
        base=summaries[cid,'original']
        for policy,_,_ in protocol['policies'][1:]:
            for t,b in zip(summaries[cid,policy]['targets'],base['targets']):
                basis_differences.append({'curve':cid,'policy':policy,'target':t['target'],
                                          'initial_shift_bits':t['initial_kappa_bits']-b['initial_kappa_bits'],
                                          'full_hidden_endpoint_shift_bits':t['full_hidden_kappa_bits']-b['full_hidden_kappa_bits']})
    output={'schema':'rank-accessibility.subsets-report.v1','claim':'Order-dependent and presentation-dependent sparse coordinate gains.',
            'protocol_sha256':digest(OUT/'protocol.json'),'verifier_output_sha256':digest(OUT/'verified.json'),
            'reporter_sha256':digest(Path(__file__)), 'policy_summaries':list(summaries.values()),
            'basis_endpoint_shifts':basis_differences,
            'verified_cells':cells,'verified_distinct_slopes':measured,'verified_edges':edges}
    for path,raw in ((OUT/'summary.json',(json.dumps(output,sort_keys=True,indent=2)+'\n').encode()),
                     (OUT/'REPORT.md','\n'.join(lines).encode())):
        if args.check:assert path.read_bytes()==raw, f'Report differs: {path}'
        else:path.write_bytes(raw)
    print('RANK_ACCESSIBILITY_SUBSETS_REPORT PASS')


if __name__=='__main__':main()
