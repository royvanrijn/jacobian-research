"""Join existing outcome summaries only after immutable M18 feature sealing."""
import csv
import time
from m18_landscape_comparison import OUT, ROOT, SNAP, read, write, sha, require

def joined():
    seal=read(OUT/'feature-seal.json')
    require(seal['protocol_sha256']==sha(OUT/'protocol.json'),'protocol changed')
    for path,digest in seal['features'].items():require(sha(OUT/path)==digest,'features changed before label join')
    rows=[]
    for case in read(OUT/'protocol.json')['cases']:
        f=read(OUT/'features'/f"{case['case']}.json")
        folder=ROOT/case['folder']
        candidates=[folder/'seeded-verified.json',folder/'trial-verified.json']
        path=next((p for p in candidates if p.exists()),None)
        if path:
            label=read(path)
            require(label['status'] in ['PASS_INDEPENDENT_SEEDED_V3_REPLAY','PASS_INDEPENDENT_DET1092_REPLAY'],'unverified outcome')
            require(label['initial_rank']==18,'label starts at wrong rank')
            rank=label['rank_lower_bound']
            outcome={'rank':rank,'gain':rank-18,'charts':label['charts'],'stop_reason':label['stop_reason'],
                     'source':str(path.relative_to(ROOT)),'sha256':sha(path)}
        else:outcome={'rank':None,'gain':None,'charts':None,'stop_reason':'UNKNOWN_NO_SEALED_REPLAY'}
        rows.append({'case':case['case'],'parameter':f['parameter'],'j':f['j'], 'outcome':outcome,'features':f['features']})
    return rows

def main():
    rows=joined()
    write(OUT/'outcomes.json',{'feature_seal_sha256':sha(OUT/'feature-seal.json'),'joined_unix':time.time(),'rows':rows})
    flat=[{'case':r['case'],'parameter':r['parameter'],**r['outcome'],**r['features']} for r in rows]
    fields=list(dict.fromkeys(k for r in flat for k in r))
    with (OUT/'comparison.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(flat)
    positives=[r for r in rows if r['outcome']['gain'] is not None and r['outcome']['gain']>0]
    nulls=[r for r in rows if r['outcome']['gain']==0]
    margins=[]
    for name in rows[0]['features']:
        a=[r['features'][name] for r in positives];b=[r['features'][name] for r in nulls]
        if any(v is None for v in a+b) or not a or not b:continue
        span=max(a+b)-min(a+b)
        if not span:continue
        lowgap=min(b)-max(a);highgap=min(a)-max(b)
        if max(lowgap,highgap)>0:
            margins.append({'feature':name,'positive_min':min(a),'positive_max':max(a),
              'null_min':min(b),'null_max':max(b),'direction':'lower' if lowgap>0 else 'higher',
              'gap_over_observed_range':max(lowgap,highgap)/span})
    margins.sort(key=lambda r:-r['gap_over_observed_range'])
    write(OUT/'separations.json',{'description':'Descriptive scan of all frozen scalar features, any certified gain vs bounded no-gain. No fitted model, no holdout validation.',
      'positive_states':len(positives),'positive_distinct_j':len({r['j'] for r in positives}),
      'null_states':len(nulls),'null_distinct_j':len({r['j'] for r in nulls}),
      'strict_separations':margins})
    lines=['# Frozen M18 landscape comparison','',
      'All CVPs use the retained rounded integer Gram. Norms are divided by the median of the first 17 diagonal entries. Outcomes are certified lower bounds after bounded V3; 18 is not an upper bound. Quartic sizes cover all 64 refined parities, irrespective of which charts ran.','',
      '| M18 state / parameter | Outcome | CVP min / median / max | ≤2 masks | Within 1.25×min | Anchors | Quartic bits median | Tail max/min |',
      '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        f=r['features'];label=r['parameter'] or '302 '+r['case'].split('__')[1]
        end=r['outcome']['rank'];outcome=f'18→{end}' if end is not None else 'pending'
        lines.append(f"| {label} | {outcome} | {f['cvp_min']:.3f} / {f['cvp_median']:.3f} / {f['cvp_max']:.3f} | {f['masks_le_2']} | {f['competitive_masks']} | {f['competitive_anchors']} | {f['quartic_bits_median']:.0f} | {f['tail_ratio']:.3f} |")
    lines+=['','## Strict retrospective separations','',
      '| Feature | Gain states range | Bounded null range | Gap / total observed range |','|---|---:|---:|---:|']
    for m in margins:lines.append(f"| {m['feature']} | {m['positive_min']:.5g}–{m['positive_max']:.5g} | {m['null_min']:.5g}–{m['null_max']:.5g} | {m['gap_over_observed_range']:.3f} |")
    lines+=['','This scan is exploratory and uses many correlated features. Curve 302 repeats one fibre and its seeds were obtained retrospectively. No separating threshold has been prospectively validated. Missing panel outcomes are excluded from the separation scan.','']
    (OUT/'comparison.md').write_text('\n'.join(lines))
    print('\n'.join(lines))

if __name__=='__main__':main()
