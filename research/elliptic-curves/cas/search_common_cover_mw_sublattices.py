#!/usr/bin/env python3
"""Fixed second pass: bounded representatives of every short-ball parity coset.

Only --evaluate reads the known generic embedding. Uses both curves identically.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
import search_low_height_mw_sublattices as base
from latent_lattice import exact_rational_ranks

PREFIX='common_cover_mw_sublattices_v1'
POLICY={'cosets':'every primitive ray in the v1 short ball',
        'representatives':'c+2z with z=0 or a signed support<=2 vector in the ambient LLL basis',
        'ranks':list(range(8,21)), 'retain_per_coset':128,
        'independence':'first independent representatives in increasing height order',
        'shortlist_per_rank_per_score':8,
        'scores':['generated_log_determinant','largest_generator_height','height_ratio'],
        'boundary':'Numerical ordering; exact rank, same-parity identities and generated subgroup indices. No exhaustive CVP or optimality claim.'}


def candidates(doc):
    h=np.asarray(doc['height_gram'],dtype=float); n=len(h)
    u=np.asarray(doc['lll_transform_rows'],dtype=np.int64)
    half=np.asarray(list(base.shell(n,2)),dtype=np.int64)@u.T
    shifts=np.concatenate((np.zeros((1,n),dtype=np.int64),2*half,-2*half))
    shift_norms=np.einsum('ij,jk,ik->i',shifts,h,shifts)
    cholesky=np.linalg.cholesky(h)
    byrank=defaultdict(list); total=0
    for ci,c in enumerate(doc['ball']):
        c=np.asarray(c,dtype=np.int64)
        norms=shift_norms+2*shifts@(h@c)+c@h@c
        indices=np.argsort(norms,kind='stable')[:POLICY['retain_per_coset']]
        vectors=[]; heights=[]; q=[]; seen=set(); logdet=0.
        for idx in indices:
            v=c+shifts[idx]
            if v[np.flatnonzero(v)[0]]<0:
                v=-v
            key=tuple(v)
            if key in seen:
                continue
            seen.add(key)
            y=v@cholesky
            residual=y.copy()
            for _ in range(2):
                for direction in q:
                    residual-=np.dot(residual,direction)*direction
            rn=residual@residual
            if rn<1e-8:
                continue
            q.append(residual/np.sqrt(rn)); vectors.append(v.tolist());heights.append(float(norms[idx]));logdet+=float(np.log(rn))
            r=len(vectors)
            if r>=8:
                byrank[r].append({'rank':r,'channel':'common_cover','seed_index':ci,
                    'coset_representative':c.tolist(),'basis_rows':[x[:] for x in vectors],
                    'generated_log_determinant':logdet,'largest_generator_height':max(heights),
                    'height_ratio':max(heights)/min(heights)})
                total+=1
            if r==20:
                break
        if (ci+1)%500==0:
            print(f"COSETS|{doc['curve']}|{ci+1}/{len(doc['ball'])}",flush=True)
    retained={}
    for r,rows in byrank.items():
        for score in POLICY['scores']:
            for row in sorted(rows,key=lambda p:(p[score],p['seed_index']))[:8]:
                key=(r,row['seed_index'])
                retained.setdefault(key,{**row,'preselection_scores':[]})['preselection_scores'].append(score)
    return {'examined_representatives':len(shifts)*len(doc['ball']),
            'proposal_count':total,'proposals':list(retained.values()),
            'all_proposals':[p for rows in byrank.values() for p in rows]}


def run(curve):
    doc=base.cloud(curve)
    path=base.OUT/f'{PREFIX}_{curve}_proposals.json.gz'
    if path.exists():
        proposed=base.read(path)
    else:
        proposed=candidates(doc);base.save(path,proposed)
    cp=base.OUT/f'{PREFIX}_{curve}_candidates.json.gz'
    if cp.exists():
        scored=base.read(cp)
    else:
        scored=base.score_batch(doc,proposed['proposals']);base.save(cp,scored)
    selected=base.finalists(doc,scored)
    for p in selected:
        p['generated_determinant']=format(float(p['determinant'])*p['generated_index_in_primitive_closure']**2,'.17g')
        p['generated_basis_is_one_parity_coset']=all(
            all((a-b)%2==0 for a,b in zip(row,p['coset_representative'])) for row in p['basis_rows'])
        if not p['generated_basis_is_one_parity_coset']:
            raise ArithmeticError('cover parity failed')
    base.save(base.OUT/f'{PREFIX}_{curve}_selection.json',{
        'curve':curve,'policy':POLICY,'examined_representatives':proposed['examined_representatives'],
        'proposal_count':proposed['proposal_count'],'shortlisted_count':len(scored),'finalists':selected,
        'status':'BOUNDED_COMMON_COVER_SUBLATTICES_NO_GENERIC_IDENTIFICATION',
        'inputs':{str(p.relative_to(base.ROOT)):base.digest(p) for p in (path,cp,Path(__file__))}})
    print(f"COVER_SELECTED|{curve}|{len(selected)}",flush=True)


def evaluate():
    selected_path=base.OUT/f'{PREFIX}_245_selection.json'
    selected=base.read(selected_path)
    truth_path=base.OUT/'latent_lattice_calibration_truth_v1.json'
    truth=next(t for t in base.read(truth_path)['fermigier_family_controls'] if t['label'].startswith('ICARM_245_'))
    t=truth['embedding_matrix_columns']
    path=base.OUT/f'{PREFIX}_245_proposals.json.gz'
    # Test recall against ALL rank-12 proposals, not merely preselected winners.
    rows=[p for p in base.read(path)['all_proposals'] if p['rank']==12]
    ranks=exact_rational_ranks([p['basis_rows']+t for p in rows],timeout=120)
    exact=[p for p,r in zip(rows,ranks) if r==12]
    selected_seeds={p['seed_index'] for p in selected['finalists'] if p['rank']==12}
    actual=[]
    for p in exact:
        # Equality of row Z-lattices, not just Q-spaces or equal determinants.
        lines=base.run_gp(f'A={base.gp_matrix(p["basis_rows"])};T={base.gp_matrix(t)};print(mathnf(A~)==mathnf(T~));',timeout=120)
        if lines==['1']:
            actual.append(p['seed_index'])
    outcome={'status':'PASS_ACTUAL_SUBGROUP_RECOVERY' if set(actual)&selected_seeds else 'FAIL_BLIND_CALIBRATION',
        'best_rank12_intersection':max(24-r for r in ranks),
        'exact_primitive_space_seed_indices':[p['seed_index'] for p in exact],
        'exact_generated_subgroup_seed_indices':actual,
        'selected_exact_primitive_space_seed_indices':[p['seed_index'] for p in exact if p['seed_index'] in selected_seeds],
        'selected_exact_generated_subgroup_seed_indices':sorted(set(actual)&selected_seeds),
        'inputs':{str(p.relative_to(base.ROOT)):base.digest(p) for p in (selected_path,truth_path,path)}}
    base.save(base.OUT/f'{PREFIX}_calibration.json',outcome);print(json.dumps(outcome,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curve',type=int,choices=(245,302));parser.add_argument('--evaluate',action='store_true')
    args=parser.parse_args()
    base.save(base.OUT/f'{PREFIX}_protocol.json',{'policy':POLICY,'sha256':base.digest(Path(__file__)),
        'base_script_sha256':base.digest(Path(base.__file__))})
    if args.evaluate:
        evaluate()
    elif args.curve:
        run(args.curve)
    else:
        parser.error('choose --curve or --evaluate')
