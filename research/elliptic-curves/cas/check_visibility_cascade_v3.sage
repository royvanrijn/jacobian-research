#!/usr/bin/env sage-python
"""Read-only replay of full-extension coverage, exact CVP and chart/rank witnesses."""
import argparse
import time
from visibility_selection_v3 import cheap_shortlist, final_shortlist, chart_profile
import csv
import itertools
import json
import sys
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import numpy as np
from sage.all import ZZ, matrix
from visibility_lattice_v2 import ExactParity
from research_runtime.store import checkpoint

CAS = Path(__file__).resolve().parent
v2 = SourceFileLoader('replayed_v2',str(CAS/'adaptive_visibility_cascade_v3.sage')).load_module()


def unit_tests():
    count = 0
    for g in ([[10,3,2],[3,7,1],[2,1,6]], [[2,0,0],[0,2,0],[0,0,2]], [[4,3],[3,4]]):
        e = ExactParity(g)
        for parity in itertools.product(range(2),repeat=len(g)):
            seed, _ = e.babai([parity])
            result = e.solve(parity,seed[0])
            vectors = [v for v in itertools.product(range(-6,7),repeat=len(g)) if tuple(x%2 for x in v)==parity]
            best = min(map(e.norm,vectors))
            minima = sorted(v for v in vectors if e.norm(v)==best)
            assert result['norm']==best and result['minima']==minima
            count += 1
    print('EXACT CVP independent brute-force tests',count,flush=True)
    return count


def replay(progress=None, start=30):
    v2.v1.guard()
    try:
        (v2.ART/'curve302_recovered_mw17_parent_v1.json').open()
    except PermissionError:
        rejected_oracle = True
    else:
        raise ArithmeticError('artifact guard accepted unredacted parent')
    d = v2.D
    p = v2.read(d/'protocol.json')
    assert p['sources']==v2.sources()
    assert all(v2.sha(v2.ROOT/k)==h for k,h in p['inputs'].items())
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl
    import pari_pointed_backend as backend
    from pointed_quartic_search import PointedQuarticSearch
    from research_runtime.search_state import raw_state
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    if start==17:
        model, basis = v2.v1.seed({'parameter':'0','presentation':'normalized'})
    else:
        fixed=v2.read(d/'fixed-M30.json')
        model=tuple(map(F,fixed['curve']))
        basis=tuple(tuple(map(F,q)) for q in fixed['points'])
    mapper=v2.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000,silent=True)
    reports = []
    tested = set()
    prior = v2.read(progress) if progress and progress.exists() else {}
    if prior:
        assert prior['checker_sha256']==v2.sha(Path(__file__))
        assert prior['protocol_sha256']==v2.sha(d/'protocol.json')
    for wd in sorted((d/('terminal-M30' if start==30 else 'replay-M17')).glob('epoch-*')):
        if not (wd/'selection.json').exists():
            break
        index = len(reports)
        if index<len(prior.get('stages',[])):
            old = prior['stages'][index]
            assert all(v2.sha(v2.ROOT/k)==h for k,h in old['checkpoint_hashes'].items())
            assert old['before']==len(basis)
            for path in sorted(wd.glob('chart-*.json')):
                tested.add(tuple(v2.read(path)['centre']['point']))
            stage = v2.read(wd/'stage.json')
            basis = tuple(tuple(map(F,q)) for q in v2.read(wd/stage['audit'])['independent_points'])
            assert len(basis)==old['after']
            reports.append(old)
            print('HASH-VERIFIED PREVIOUS INDEPENDENT REPLAY',index,flush=True)
            continue
        selection = v2.read(wd/'selection.json')
        assert selection['basis']==[list(map(str,q)) for q in basis]
        g = matrix(ZZ,selection['rounded_gram'])
        u = matrix(ZZ,selection['LLL'])
        assert abs(u.det())==1
        e = ExactParity((u*g*u.transpose()).rows())
        count = 1 << (len(basis)-17)
        assert count==selection['extensions_per_anchor']
        columns = v2.fingerprints(model,basis)
        assert columns==selection['fingerprint_columns']
        ext_columns = columns[17:]
        rebased_columns = ext_columns[:]
        if len(rebased_columns)>1:
            rebased_columns[0] ^= rebased_columns[1]
            rebased_columns.reverse()
        def all_fingerprints(cols):
            values = {0}
            for col in cols:
                values |= {value^col for value in values}
            return values
        assert all_fingerprints(ext_columns)==all_fingerprints(rebased_columns)
        assert len(all_fingerprints(ext_columns))==count
        buckets = {8:[],10:[]}
        with v2.ORBITS.open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                shell = int(row['minimum_norm'])
                if shell not in buckets:
                    continue
                word = list(map(int,row['parent_MW17_w'].split()))+[0]*(len(basis)-17)
                buckets[shell].append({'orbit':int(row['orbit_mask']),'shell':shell,
                    'representative':word,'fingerprint':v2.fp(word,columns)})
        anchors, chosen = [], []
        for shell, rows in buckets.items():
            words = np.asarray([r['representative'] for r in rows],dtype=np.int64)
            norms = np.einsum('ij,jk,ik->i',words,np.asarray(g.rows(),dtype=np.int64),words,optimize=True)
            for row,norm in zip(rows,norms): row['metric_norm']=int(norm)
            rows.sort(key=lambda r:(-r['metric_norm'],r['fingerprint']))
            anchors += rows[:p['anchors_per_shell']]
            count_canonical = 0
            for row in rows:
                key,word = v2.point_key(model,basis,row['representative'])
                if tuple(map(str,key)) in tested: continue
                chosen.append({**row,'representative':word,'point':list(map(str,key)),'lane':'canonical'})
                count_canonical += 1
                if count_canonical==p['canonical_per_shell']: break
        assert anchors==[r['anchor'] for r in selection['anchors']]
        cosets = 0
        exact_count = 0
        for ai, record in enumerate(selection['anchors']):
            path = wd/f'anchor-{ai:02d}-full.npz'
            assert v2.sha(path)==record['full_scores_sha256']
            data = np.load(path)
            residues = data['residues']
            assert residues.shape==(count,len(basis))
            assert np.all((residues==0)|(residues==1))
            assert np.all(residues[:,:17]==np.asarray(record['anchor']['representative'][:17])%2)
            assert len(set(map(tuple,residues[:,17:])))==count
            rp = (residues @ (np.asarray(u.inverse().rows(),dtype=np.int64)%2))%2
            assert np.array_equal(rp,data['reduced_residues'])
            words, norms = e.babai(rp)
            assert np.array_equal(words,data['babai_words']) and np.array_equal(norms,data['babai_norms'])
            keys = [v2.fp(row,columns) for row in residues]
            order = cheap_shortlist(norms,keys)
            assert set(order)=={r['extension'] for r in record['refined']}
            for row in record['refined']:
                i = row['extension']
                result = e.solve(rp[i],words[i],p['exact_cvp_node_limit'])
                assert json.dumps(result,sort_keys=True)==json.dumps(row['cvp'],sort_keys=True)
                assert row['metric_norm']==result['norm'] and row['fingerprint']==keys[i]
                word = matrix(ZZ,1,len(basis),row['representative'])
                assert int((word*g*word.transpose())[0,0])==result['norm']
                assert [int(x)%2 for x in word.row(0)]==list(residues[i])
                key, _ = v2.point_key(model,basis,row['representative'])
                assert list(map(str,key))==row['point']
                minimizing_keys = []
                for minimum in result['minima']:
                    transported = list(map(int,(matrix(ZZ,1,len(basis),minimum)*u).row(0)))
                    minimizing_keys.append(v2.point_key(model,basis,transported)[0])
                assert key==min(minimizing_keys), 'minimum-centre semantic tie changed'
                profile=chart_profile(mapper.mapping(model,basis,row))
                assert all(row[k]==v for k,v in profile.items())
                assert row['multiplicity']==len(result['minima'])//2
                exact_count += 1
            ranked = sorted(record['refined'],key=lambda r:(-r['metric_norm'],r['fingerprint']))
            chosen += final_shortlist([r for r in ranked if tuple(r['point']) not in tested])
            cosets += count
        assert cosets==selection['full_cosets_scored']
        chosen.sort(key=lambda r:(-r['metric_norm'],r['fingerprint'],r['lane'],tuple(map(F,r['point']))))
        unique, seen = [], set(tested)
        for row in chosen:
            key = tuple(row['point'])
            if key not in seen:
                unique.append(row);seen.add(key)
        assert unique==selection['centres'], 'final schedule differs'
        print('REPLAYED FULL V3 LANDSCAPE',selection['rank'],cosets,exact_count,flush=True)
        while not (wd/'stage.json').exists():
            if (d/('terminal-M30' if start==30 else 'replay-M17')/'terminal.json').exists():
                raise ArithmeticError('terminal without completed stage')
            time.sleep(5)
        state = raw_state(model,basis,cache=Cache(MemoryFactStore()),prime_bound=1000)
        charts = sorted(wd.glob('chart-*.json'))
        for j,path in enumerate(charts):
            chart = v2.read(path)
            assert chart['centre']==selection['centres'][j]
            tested.add(tuple(chart['centre']['point']))
            search = PointedQuarticSearch(state=state,centre={'coefficients':chart['centre']['representative']},coordinate_policy=chart['mapping']['coordinate_policy'])
            backend.replay(search,chart['mapping'],chart['search'])
            audit = wd/f'mod2-{j:03d}.json'
            mod2.check(audit)
            if j<len(charts)-1:
                assert v2.read(audit)['rank_lower_bound']==len(basis), 'stale chart executed after gain'
        stage = v2.read(wd/'stage.json')
        assert stage['audit_sha256']==v2.sha(wd/stage['audit'])
        cloud = v2.read(wd/stage['audit'])
        enlarged = tuple(tuple(map(F,q)) for q in cloud['independent_points'])
        assert enlarged[:len(basis)]==basis and len(enlarged)==stage['after']
        if len(enlarged)>len(basis):
            modl.check(wd/'modl.json')
        modl_ranks = {str(a['modulus']):a['finite_column_rank'] for a in v2.read(wd/'modl.json')['audits']} if (wd/'modl.json').exists() else {}
        reports.append({'epoch':stage['epoch'],'before':len(basis),'after':len(enlarged),
                        'full_cosets_replayed':cosets,'exact_CVP_replayed':exact_count,'charts_replayed':len(charts),
                        'independent_modl_ranks':modl_ranks,
                        'extension_basis_shear_and_reversal_coverage':'passed',
                        'checkpoint_hashes':{v2.rel(path):v2.sha(path) for path in sorted(wd.iterdir()) if path.suffix in ('.json','.npz')},
                        'stale_charts_cancelled':stage['stale_charts_cancelled']})
        basis = enlarged
        print('REPLAYED V3 EPOCH',{k:v for k,v in reports[-1].items() if k!='checkpoint_hashes'},flush=True)
        if progress:
            checkpoint(progress,{'checker_sha256':v2.sha(Path(__file__)),
                'protocol_sha256':v2.sha(d/'protocol.json'),'stages':reports})
    return {'schema':'visibility-cascade-v3-replay','protocol_sha256':v2.sha(d/'protocol.json'),
            'checker_sha256':v2.sha(Path(__file__)),'stages':reports,'rank_lower_bound':len(basis),
            'guard_rejected_unredacted_parent':rejected_oracle,
            'claim_boundary':'Calibration execution with full extension coverage and exact rounded-metric CVP on a heuristic shortlist. Not full-quotient chart search or arbitrary-basis invariance.'}


if __name__=='__main__':
    ap = argparse.ArgumentParser();ap.add_argument('--unit-only',action='store_true');ap.add_argument('--output',type=Path);ap.add_argument('--progress',type=Path)
    ap.add_argument('--start',type=int,choices=[17,30],default=30)
    args = ap.parse_args();unit_tests()
    if not args.unit_only:
        result = replay(args.progress,args.start)
        if args.output:
            checkpoint(args.output,result)
        else:
            print(json.dumps(result,indent=2))
