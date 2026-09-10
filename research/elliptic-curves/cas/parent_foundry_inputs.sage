#!/usr/bin/env sage-python
"""Freeze generic geometry, separating invariants and deterministic parent proposals."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ, ZZ, PolynomialRing, matrix, vector

CAS=Path(__file__).resolve().parent; ROOT=CAS.parents[1]
sys.path.insert(0,str(CAS))
geo=SourceFileLoader('parent_inputs_geometry',str(CAS/'parent_foundry_geometry.sage')).load_module()
from research_runtime.store import checkpoint


def build(output):
    art=ROOT/'artifacts/generated-results/elliptic-curves'
    atlas=json.loads(geo.ATLAS.read_text())
    keys=('family','A_coefficients_low_to_high','B_coefficients_low_to_high','sections',
          'generic_height_gram','generic_rank_lower_bound')
    sources=[{k:f[k] for k in keys} for f in atlas['families']]
    for f in sources:
        f.update(source_surface='X948',generic_rank_upper_bound=17,
                 source_certificate=str(geo.ATLAS.relative_to(ROOT)),source_certificate_sha256=geo.digest(geo.ATLAS))
    parent_path=ROOT/'parent-inputs/det1092-gram.json'
    reduced_path=art/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
    reduced=json.loads(reduced_path.read_text())
    if parent_path.exists():parent=json.loads(parent_path.read_text())
    else:
        # Development preparation only. Production uses the narrow projection.
        parent_path=art/'curve302_recovered_mw17_parent_v1.json'
        parent=json.loads(parent_path.read_text())
    R=PolynomialRing(QQ,'t'); K=R.fraction_field()
    def decode(r):return K(R(r['numerator']))/R(r['denominator'])
    a1,a2,a3,a4,a6=map(decode,reduced['a_invariants'])
    if a1 or a2 or a3:raise ArithmeticError('reduced source must be short')
    f={'family':'det1092','source_surface':'X1092','generic_rank_lower_bound':17,
        'generic_rank_upper_bound':17,'A_coefficients_low_to_high':geo.record(R(a4)),
        'B_coefficients_low_to_high':geo.record(R(a6)),
        'sections':[{'basis_index':i,'X':geo.encode(decode(p[0])),
                     'Y':geo.encode(decode(p[1]))} for i,p in enumerate(reduced['basis_weierstrass_coordinates'])],
        'generic_height_gram':parent['generic_height_gram'],
        'source_certificate':str(reduced_path.relative_to(ROOT)),
        'source_certificate_sha256':geo.digest(reduced_path),
        'generic_gram_source_sha256':geo.digest(parent_path)}
    sources.append(f)
    for f in sources:
        g=matrix(QQ,f['generic_height_gram'])
        if len(f['sections'])!=17 or not g.is_positive_definite():raise ArithmeticError('invalid source')
        A,B=(R(f[k]) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
        D=4*A**3+27*B**2
        if D.degree()!=24 or not D.is_squarefree():raise ArithmeticError('source must have 24I1')
        for s in f['sections']:
            x,y=geo.decode(s['X'],R),geo.decode(s['Y'],R)
            if y*y!=x*x*x+A*x+B:raise ArithmeticError('source section failed')
    template_path=ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json'
    template=json.loads(template_path.read_text()); known=[]
    for p in template['presentations']:
        a,b=(R(p['pencil'][k]) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
        known.append({'presentation':p['presentation_id'],'family':p['fibration_id'],
                      'invariant':geo.affine_invariant(a,b)})
    # Mandatory equivalence regression: nine presentations must give five keys.
    byfamily={}
    for row in known:
        byfamily.setdefault(row['family'],set()).add(row['invariant']['key'])
    if any(len(v)!=1 for v in byfamily.values()) or len({r['invariant']['key'] for r in known})!=5:
        raise ArithmeticError('known A1 collision regression failed')
    baselines=list(sources)
    for p in json.loads((art/'compact_five_mw16_atlas_v1.json').read_text())['families']:
        baselines.append({**{k:p[k] for k in keys if k!='family'},'family':p['fibration_id'],
            'generic_rank_upper_bound':16,'source_surface':'X948','baseline':True})
    pool=[]
    for row in csv.DictReader(geo.TABLE.open(),delimiter='\t'):
        if row['minimal_unoriented_count']=='1':
            pool.append({'source_family':'11952','priority':int(row['priority_rank']),
                         'trace_word':list(map(int,row['section_basis_w'].split()))})
    detrows=[]
    dettable=art/'curve302_parent_degree2_multisection_orbits_v1.tsv'
    for row in csv.DictReader(dettable.open(),delimiter='\t'):
        if row['minimum_norm']=='8':
            word=list(map(int,row['parent_MW17_w'].split()))
            detrows.append({'source_family':'det1092','priority':int(row['orbit_mask']),
                            'trace_word':word})
    detrows.sort(key=lambda r:hashlib.sha256(f'parent-construction-v1/{r["priority"]}'.encode()).hexdigest())
    # The cheapest coordinate words overrepresent pencils with many I2 fibres.
    # Interleave the full population before using cost within small blocks.
    detrows=[r for i in range(0,len(detrows),64) for r in sorted(detrows[i:i+64],
        key=lambda r:(sum(abs(v) for v in r['trace_word']),r['priority']))]
    # Batch cheap exact lattice gates so a run does not spend one CAS startup
    # per rejected multi-I2 proposal before reaching its first X1092 parent.
    import numpy as np
    from visibility_lattice_fast import IntegerExactParity
    gram=matrix(QQ,sources[-1]['generic_height_gram']);change=gram.LLL_gram().transpose()
    inverse=change.inverse();solver=IntegerExactParity((2*change*gram*change.transpose()).rows())
    preflight=[];eligible=[]
    for row in detrows[:4096]:
        w=matrix(ZZ,1,17,row['trace_word']);residue=tuple(int(x)%2 for x in (w*inverse).row(0))
        starts,_=solver.babai(np.asarray([residue],dtype=np.int64))
        proof=solver.solve(residue,tuple(map(int,starts[0])),2000000)
        preflight.append({'priority':row['priority'],'norm':proof['norm'],'signed_minima':len(proof['minima'])})
        if proof['norm']==16 and len(proof['minima'])==2:eligible.append(row)
    original_det_count=len(detrows)
    detrows=eligible+detrows[4096:]
    # Alternating construction lanes. The 1092 proposals still need their exact
    # unique-I2 equation gate; a norm-eight label alone does not imply MW16.
    queue=[]
    for i in range(max(len(pool),len(detrows))):
        if i<len(pool):queue.append(pool[i])
        if i<len(detrows):queue.append(detrows[i])
    checkpoint(output,{'schema':'parent-foundry.generic-inputs.v1','sources':sources,
        'baselines':baselines,'known_a1':known,'proposals':queue,
        'counts':{'proved_x948_a1_proposals':len(pool),'det1092_norm8_proposals':original_det_count,
                  'det1092_prefiltered':len(preflight),'det1092_a1_in_prefilter':len(eligible)},
        'det1092_lattice_preflight':preflight,
        'bindings':{str(p.relative_to(ROOT)):geo.digest(p) for p in
                    (geo.ATLAS,geo.TABLE,parent_path,reduced_path,template_path,dettable)},
        'selection':'Generic norm-eight classes only; low group-law cost order, alternating X948 and X1092. No specialization score or record fibre enters selection.',
        'rank_completion_gate':'Exact A1 MW16 has no 17th section. Rank-completion is inadmissible here; fibration changes are explicit constructions.',
        'claim_boundary':'Proposals are not accepted parents. Exact unique-I2, section, independence and separating-invariant gates precede evaluation.'})
    print('GENERIC_INPUTS',len(pool),len(detrows),len(baselines),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args().output)
