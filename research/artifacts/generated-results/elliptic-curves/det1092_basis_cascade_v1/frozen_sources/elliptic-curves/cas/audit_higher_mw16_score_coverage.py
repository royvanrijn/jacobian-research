#!/usr/bin/env python3
"""Independent model reconstruction and strict retained-prefix coverage proof."""
import argparse,json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'

def expected(broad):
    prefix='broad' if broad else 'joint'
    scan=LOCAL/(prefix+'-mw16-higher-annuli-v1')
    scores=LOCAL/(prefix+'-mw16-higher-scores-v1')
    protocol=cert.read(scores/'protocol.json');population=cert.read(scan/'result.json')
    replay=cert.read(scan/'replay.json');control=cert.read(scan/'controller/ledger.json')
    if population['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or control['status']!='PASS' or replay['status']!='PASS' or replay['result_sha256']!=cert.hashed(scan/'result.json'):
        raise ArithmeticError('complete fixed initial population and score replay required')
    if len(protocol['rows'])!=10240 or any(cert.hashed(ROOT/n)!=h for n,h in protocol['sources'].items()):
        raise ArithmeticError('frozen scalar inputs changed')
    families={r['fibration_id']:r for r in cert.read(spec.ATLAS)['families']}
    groups=defaultdict(list);by_slice=defaultdict(list);seen={};model_count=0
    for row in protocol['rows']:
        n,d=row['numerator'],row['denominator'];t=F(n,d)
        lo,hi={2:(16384,65536),3:(65536,262144)}[row['band']]
        if str(t)!=row['parameter'] or (t.numerator,t.denominator)!=(n,d) or not lo<max(abs(n),d)<=hi:
            raise ArithmeticError('primitive parameter or higher band differs')
        f=families[row['family']];model=[F(0)]*3
        for key,weight in [('A_coefficients_low_to_high',8),('B_coefficients_low_to_high',12)]:
            coefficients=list(map(F,f[key]))
            if len(coefficients)>weight+1:raise ArithmeticError('unexpected homogeneous degree')
            model.append(sum(c*n**i*d**(weight-i) for i,c in enumerate(coefficients)))
        if any(a.denominator!=1 for a in model) or list(map(str,model))!=row['model']:
            raise ArithmeticError('independent exact homogeneous model differs')
        a,b=model[3:];delta=-16*(4*a**3+27*b**2)
        if not delta:raise ArithmeticError('singular model')
        j=(-48*a)**3/delta
        if any(cert.isomorphic(model,other) for other in seen.get(j,[])):
            raise ArithmeticError('prospective equation duplicate')
        seen.setdefault(j,[]).append(model);model_count+=1
        groups[(row['band'],row['family'])].append(row);by_slice[row['slice_id']].append(row)
    if len(groups)!=10 or any(len(rows)!=1024 for rows in groups.values()):
        raise ArithmeticError('equal ten-group scalar allocation differs')
    def key(row):return (-row['combined_selection_units'],-row['combined_good'],row['denominator'])
    thresholds={g:max(key(r) for r in rows) for g,rows in groups.items()}
    boundaries=defaultdict(list);matched=set();paths=[]
    for manifest in population['shards']:
        if broad:
            path=ROOT/manifest['retained_path'];paths.append(path)
            if path!=scan/manifest['id']/'retained.json' or cert.hashed(path)!=manifest['retained_sha256']:
                raise ArithmeticError('immutable broader shard differs')
            shard=cert.read(path)
        else:shard=manifest
        if len(shard['rows'])!=4096:raise ArithmeticError('fixed slice retention differs')
        g=(shard['band'],shard['family']);boundary=key(shard['rows'][-1])
        boundaries[g].append({'slice':shard['id'],'last_retained_key':list(boundary),
                              'selected_threshold_strictly_better':thresholds[g]<boundary})
        index={r['parameter']:r for r in shard['rows']}
        for row in by_slice[shard['id']]:
            if row['id'] in matched or row['parameter'] not in index:raise ArithmeticError('selected address missing or repeated')
            original=index[row['parameter']]
            keys=('numerator','denominator','score_units','good_primes','extension_selection_units',
                  'extension_good','combined_selection_units','combined_good')
            if any(row[k]!=original[k] for k in keys):raise ArithmeticError('scalar selection changed a retained address or score')
            matched.add(row['id'])
    if len(matched)!=10240:raise ArithmeticError('selected address coverage incomplete')
    records=[]
    for g,rows in sorted(groups.items()):
        b=boundaries[g]
        if len(b)!=(32 if broad else 2):raise ArithmeticError('complete signed slice count differs')
        records.append({'band':g[0],'family':g[1],'scalar_candidates':len(rows),
            'parameter_height_minimum':min(max(abs(r['numerator']),r['denominator']) for r in rows),
            'parameter_height_maximum':max(max(abs(r['numerator']),r['denominator']) for r in rows),
            'last_selected_primary_key':list(thresholds[g]),'slice_boundaries':b,
            'strict_retention_coverage':all(r['selected_threshold_strictly_better'] for r in b)})
    coverage=all(r['strict_retention_coverage'] for r in records)
    inputs=[Path(__file__).resolve(),Path(cert.__file__),spec.ATLAS,
            scores/'protocol.json',scan/'protocol.json',scan/'result.json',scan/'replay.json',
            scan/'controller/ledger.json',*paths]
    return {'schema':'elliptic-curves.higher-mw16-model-and-retention-audit.v1',
        'status':'PASS' if coverage else 'PASS_MODELS_WITH_UNRESOLVED_RETENTION_COVERAGE',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in inputs},
        'variant':'broad' if broad else 'narrow','models_checked':model_count,
        'distinct_j_invariants':len(seen),'distinct_Q_isomorphism_classes':model_count,
        'groups':records,'strict_retention_coverage':coverage,
        'argument':'The scanner retains the best4096 of each signed slice, ordered by decreasing combined units, decreasing good count, denominator, absolute numerator. Every omitted address is no better than its slice boundary in the first three fields. If the last admitted scalar candidate in a band/family is strictly better than every slice boundary in those fields, no omitted address can precede that candidate, regardless of the signed-numerator tie convention. Thus truncation cannot remove any address before this selected prefix, including any earlier within-roster duplicate that the frozen scalar selection skips.',
        'claim_boundary':'Independent exact homogeneous models, primitive higher-band membership, within-roster Q-isomorphism distinctness, saved score provenance and finite strict-prefix coverage through prime32749. This preserves the existing selector. It does not prove the globally best65521 scores, a rank predictor, coverage outside the declared slices, catalogue absence, rational-point independence or a new rank.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--broad',action='store_true');parser.add_argument('--check',action='store_true');a=parser.parse_args()
    output=ART/(('broad' if a.broad else 'narrow')+'_higher_mw16_score_coverage_v1.json');data=expected(a.broad)
    if a.check:
        if cert.read(output)!=json.loads(json.dumps(data)):raise ArithmeticError('model/coverage proof replay differs')
    else:
        if output.exists():raise FileExistsError('preserve model/coverage proof')
        checkpoint(output,data)
    print('HIGHER MW16 MODEL/RETENTION AUDIT',data['variant'],data['status'],data['models_checked'],data['distinct_j_invariants'],flush=True)
