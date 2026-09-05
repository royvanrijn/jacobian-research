#!/usr/bin/env sage-python
"""Portable exact witness replay and joined diagnostic report.

Reads committed artifacts only. Replays charts without enumerating search boxes,
exact group words, and every prime table with Sage finite-field point counts.
Does not recompute numerical CVP proposals, timings, or a parameter census.
"""
import gzip,json
from hashlib import sha256
from pathlib import Path
import sys,zipfile
from collections import Counter
from statistics import mean
from sage.all import EllipticCurve,GF,QQ,ZZ,matrix,lcm

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint,digest
from pointed_quartic_search import PointedQuarticSearch
from search_observability import point_visibility
from compare_bounded_prime_selectors import short_mod,contributions,ranks
from mod_l_reduction_independence import ModLReductionSignature,mod_l_reduction_signature,combined_mod_l_rank
from fractions import Fraction
from audit_mw18_retained_visibility import transport
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS,POINTS

ART=ROOT/'artifacts/generated-results/elliptic-curves'


def read(path):
    data=path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix=='.gz' else data)


def point(E,p):return E(0) if p is None else E(QQ(p['x']),QQ(p['y']))


def word_point(E,points,word):return sum((ZZ(c)*p for c,p in zip(word,points) if c),E(0))


def check_independence(row,certificate):
    model=tuple(map(Fraction,row['search_model']))
    pts=tuple((Fraction(p['x']),Fraction(p['y'])) for p in row['subgroup'])
    last=certificate['attempts'][-1];mod=last['modulus'];signatures=[]
    for s in last['signatures']:
        calculated=mod_l_reduction_signature(model,pts,s['prime'],mod)
        if calculated.rows!=tuple(tuple(r) for r in s['rows']) or calculated.group_order!=s['group_order']:
            raise ArithmeticError('original independence witness failed')
        signatures.append(calculated)
    assert combined_mod_l_rank(signatures,len(pts),mod)==len(pts)
    torsion=mod_l_reduction_signature(model,(),last['no_torsion_prime'],mod)
    assert torsion.group_order==last['no_torsion_group_order'] and torsion.group_order%mod


def run():
    paths=[ART/name for name in ('fibre_height_population_v1.json.gz','fibre_height_search_v1.json.gz',
        'ordinary_masked_controls_v1.json.gz','ordinary_masked_relations_v1.json',
        'mw18_retained_visibility_v1.json','mw18_translated_visibility_v1.json',
        'bounded_prime_selector_traces_v1.json.gz','bounded_prime_selector_comparison_v1.json',
        'point_supplied_mw19_diagnostic_v1.json')]
    population,old,masked,relations,raw,translated,traces,comparison,construction=map(read,paths)
    rows={r['id']:r for f in population['families'].values() for r in f['rows']}
    checks=Counter();masked_summary={}
    assert set(masked['results'])==set(old['results'])==set(relations['results'])
    for identifier,c in masked['results'].items():
        row=rows[identifier];check_independence(row,old['results'][identifier]['generic_independence'])
        checks['original_independent_subgroups']+=1
        E=EllipticCurve(QQ,row['search_model']);basis=[point(E,p) for p in row['subgroup']]
        rw=relations['results'][identifier];covered=set()
        for ch in c['charts']:
            assert ch['input']['curve']==row['search_model'] and ch['input']['subgroup']==row['subgroup'][1:]
            assert ch['height_bound']==100000 and ch['status']=='bounded_search_complete'
            # Search records are self-contained, so no local blind-input/cache
            # checkpoint is required by this replay.
            PointedQuarticSearch(**ch['input']).verify_record(ch)
            checks['masked_chart_replays']+=1
            covered.update(tuple(p.values()) for p in ch['finite_curve_points'])
        assert len(covered)==len(rw['exact_relations']) and not rw['unknown_point_indices']
        for r in rw['exact_relations'].values():
            pt=E(list(map(QQ,r['point'])));coeff=list(map(QQ,r['coefficients']));d=lcm(x.denominator() for x in coeff)
            assert tuple(r['point']) in covered and d*pt==word_point(E,basis,[d*x for x in coeff])
            checks['masked_exact_group_relations']+=1
        assert rw['withheld_direction_recovered']==any(QQ(r['withheld_coefficient']) for r in rw['exact_relations'].values())
    for family in population['families']:
        ids=[i for i in masked['results'] if rows[i]['id'].startswith(family+'-')]
        masked_summary[family]=dict(curves=len(ids),withheld_directions_recovered=sum(relations['results'][i]['withheld_direction_recovered'] for i in ids),
            finite_points=sum(masked['results'][i]['finite_point_count'] for i in ids),completed_charts=12*len(ids),
            worker_seconds=sum(ch['worker_seconds'] for i in ids for ch in masked['results'][i]['charts']),
            original_new_directions=sum(old['results'][i]['certified_new_directions'] for i in ids))
    print('REPLAY_MASKED',dict(checks),flush=True)
    z=zipfile.ZipFile(ART/'mw18_deep_centre_comparison_v1.zip');protocol=json.loads(z.read('protocol.json'))
    public={r['id']:r for r in read(ART/'icarm_curve_refresh_475_573_overview_v1.json')['snapshot']['records']}
    visibility_summary={}
    for identifier,c in translated['cases'].items():
        case=protocol['cases'][identifier];E=EllipticCurve(QQ,case['curve']);known=[point(E,p) for p in case['points']]
        if identifier.startswith('historical'):model,pub=GENERAL_WEIERSTRASS_COEFFICIENTS,POINTS
        else:
            row=public[int(identifier.split('-')[1])];model,pub=row['ainvs'],row['points']
        pub=[point(E,p) for p in transport(model,pub,case['curve'])]
        words=matrix(QQ,c['generic_words_in_public_basis'])
        assert words.rank()==18
        for i,p in enumerate(known):assert word_point(E,pub,words.column(i))==p
        visibility_summary[identifier]={}
        for policy,p in c['policies'].items():
            cell=json.loads(z.read(f'anchor-trial/cells/{identifier}--{policy}.json'))
            got={tuple(pt.values()):pt for ch in cell['charts'] for pt in ch['search']['finite_curve_points']}
            assert len(got)==len(p['recovered_public_basis_words'])
            for pt,word in zip(got.values(),p['recovered_public_basis_words']):
                assert word_point(E,pub,word)==point(E,pt);checks['mw18_recovered_group_relations']+=1
            span=words.augment(matrix(QQ,p['recovered_public_basis_words']).transpose()) if got else words
            assert span.rank()-18==p['exact_rational_span_gain']
            for label,best in p['minima'].items():
                index=int(label[1:])-1;unit=matrix(QQ,28,1,[int(j==index) for j in range(28)])
                assert best['already_in_recovered_rational_span']==(span.augment(unit).rank()==span.rank())
                if best['minimum_affine_height'] is None:continue
                pt=point(E,best['translated_point'])
                assert pt==best['sign']*pub[index]+word_point(E,known,best['translation_in_generic_basis'])
                checked=point_visibility(cell['charts'][best['chart']]['search'],best['translated_point'])
                for key in checked:assert checked[key]==best[key]
                assert checked['status']!='VISIBLE_NOT_RECORDED';checks['mw18_minimum_witnesses']+=1
            extras={}
            for bound in (100000,200000,1000000,10000000):
                extra=span
                for label,best in p['minima'].items():
                    if best['minimum_affine_height'] is not None and best['minimum_affine_height']<=bound:
                        index=int(label[1:])-1;extra=extra.augment(matrix(QQ,28,1,[int(j==index) for j in range(28)]))
                extras[str(bound)]=int(extra.rank()-span.rank())
            assert extras==p['witnessed_additional_rational_span_at_height'] and extras['100000']==0
            visibility_summary[identifier][policy]=dict(original_gain=p['exact_rational_span_gain'],witnessed_additional_gain=extras)
    print('REPLAY_MW18',dict(checks),flush=True)
    pp=traces['protocol'];ps=pp['training_primes'];vs=sum(pp['validation_blocks'],[])
    assert len(ps)==256 and len(vs)==128 and not set(ps)&set(vs) and set(traces['tables'])==set(map(str,ps+vs))
    for prime,table in traces['tables'].items():
        prime=int(prime)
        for key,t in table['traces'].items():
            a,b=map(int,key.split(','))
            if (4*a**3+27*b*b)%prime==0:assert t is None
            else:assert prime+1-int(EllipticCurve(GF(prime),[a,b]).cardinality())==t
            checks['independent_local_point_counts']+=1
        for identifier,key in table['candidate_lookup'].items():
            assert key==','.join(map(str,short_mod(rows[identifier]['minimal_model'],prime)))
    print('REPLAY_PRIMES',checks['independent_local_point_counts'],flush=True)
    joined=[]
    for r in comparison['comparisons']:
        ids=[i for i,v in comparison['values'].items() if v['family']==r['family'] and
            (r['subset']=='all' or r['subset']=='measured' and i in masked['results'] or v['split']==r['subset'])]
        kind=pp['scores'].index(r['score']);score={}
        for i in ids:
            score[i]=sum(contributions(traces['tables'][str(p)]['traces'][traces['tables'][str(p)]['candidate_lookup'][i]],p)[kind]
                for p in ps[:r['prefix']])
            assert score[i]==comparison['values'][i]['scores'][r['score']][str(r['prefix'])][0]
        assert r['selected']==sorted(ids,key=lambda i:(-score[i],i))[:len(r['selected'])]
        if r['subset']=='measured':joined.append({**r,'exact_withheld_direction_recoveries':sum(relations['results'][i]['withheld_direction_recovered'] for i in r['selected'])})
    legacy=[]
    for family,f in population['families'].items():
        ids=[r['id'] for r in f['rows']];chosen=f['arms']['nagao']
        for name in pp['scores']:
            vals={i:sum(s[0] for s in comparison['values'][i]['validation'][name]) for i in ids}
            percentiles=dict(zip(ids,ranks(list(vals.values()))))
            legacy.append(dict(family=family,validation_score=name,selected=chosen,
                policy='Retained original 25-prime capped Pareto selector; no refit.',
                selected_mean_validation_percentile=mean(percentiles[i]/(len(ids)-1) for i in chosen),
                original_new_directions=sum(old['results'][i]['certified_new_directions'] for i in chosen),
                exact_masked_recoveries=sum(relations['results'][i]['withheld_direction_recovered'] for i in chosen)))
    base=EllipticCurve(QQ,construction['base_curve']);p=point(base,construction['base_point']);ann=0
    from math import gcd
    for row in construction['base_good_reductions']:
        order=int(EllipticCurve(GF(row['prime']),list(base.ainvs())).cardinality());assert order==row['order'];ann=gcd(ann,order)
    assert ann==construction['base_torsion_annihilator'] and ann*p==point(base,construction['nonzero_annihilator_multiple']) and ann*p
    gram=matrix(QQ,construction['generic_gram']);assert gram.is_positive_definite() and gram.rank()==19 and gram.det()==QQ(construction['generic_gram_determinant'])
    for sample in construction['height_samples']:
        E=EllipticCurve(QQ,sample['curve'])
        for pt in sample['supplied_points']:assert point(E,pt)
    checks['point_supplied_construction']=1
    output=dict(schema='elliptic-curves.rank-jump-diagnostics-replay.v1',status='PASS',checks=dict(checks),
        inputs={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths+[Path(__file__)]},
        visibility=visibility_summary,masked_controls=masked_summary,selector_measured_comparisons=joined,
        legacy_selector_validation=legacy,construction=dict(covers=construction['covers'],generic_rank_gain=1,
            positive_rank_base=True,height_samples=[{k:r[k] for k in ('multiple','parameter_height_bits','j_log2_height','raw_model_size','point_sizes')} for r in construction['height_samples']]),
        claim_boundary='Retained enumerator completeness remains a worker trust boundary. '
        'Visibility witnesses do not certify globally optimal coordinates. '
        'The original unmasked recovery endpoint is constant zero and masked recovery is constant one; neither validates selector enrichment.')
    checkpoint(ART/'rank_jump_diagnostics_replay_v1.json',output)
    print('PASS',dict(checks),flush=True)


if __name__=='__main__':run()
