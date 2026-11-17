#!/usr/bin/env python3
"""Seal the two precise classical mutation inputs and eight-name source cues."""
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile
from audit_knotinfo_snapshots import RELEASES

HERE=Path(__file__).resolve().parent


def main():
    release=next(r for r in RELEASES if r['version']=='2026.11.2')
    archive=HERE/'tmp/sources'/release['url'].rsplit('/',1)[1]
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==release['sha256']
    selection=json.loads((HERE/'mutation_cover_plan.json').read_text())['selection']
    names={name for block in selection.values() for name in block}
    columns=['name','jones_polynomial','signature','conway_polynomial','symmetry_type']
    csv.field_size_limit(10**9)
    with zipfile.ZipFile(archive) as z:
        member=next(m for m in z.namelist() if m.endswith('/knotinfo_data_complete.csv'))
        raw=z.read(member);reader=csv.reader(io.StringIO(raw.decode('utf-8-sig')),delimiter='|')
        header=next(reader);indices={k:header.index(k) for k in columns}
        records=[{k:row[i] for k,i in indices.items()} for row in reader if row[indices['name']] in names]
    assert {r['name'] for r in records}==names
    for block in selection.values():
        rows=[r for r in records if r['name'] in block]
        assert len(rows)==4 and all(len({r[k] for r in rows})==1 for k in columns[1:])
    subset={'release':release,'csv_member':member,'csv_sha256':hashlib.sha256(raw).hexdigest(),
        'columns_retained':columns,'records':records,
        'scope':'Pinned invariant screening cues only; graph certificates, not invariant equality, prove mutation.'}
    path=HERE/'sources/knotinfo/2026.11.2_mutation_subset.json'
    if path.exists():assert json.loads(path.read_text())==subset
    else:path.write_text(json.dumps(subset,indent=2)+'\n')
    additions=[{
        'id':'GRE13','title':'Lattices, graphs, and Conway mutation','authors':'Joshua Evan Greene','year':2013,
        'doi':'10.1007/s00222-012-0421-4','url':'https://doi.org/10.1007/s00222-012-0421-4',
        'snapshot':'sources/greene-mutation-v1.pdf',
        'location':'Author preprint arXiv1103.0487v1, Section4.1-4.2 printedpp19-23, Corollary4.3, Lemma4.4, Corollary4.5; published Invent.Math.192,717-750 (2013), online11August2012',
        'audit':'Primary graph/drawing/mutation arguments read, including parallel-edge qualification and chosen checkerboard convention; official Springer publication metadata checked. Classical Whitney2isomorphism and plane-drawing results remain explicit imports. Full lattice/Floer proof and published-version comparison not audited.',
        'scope':'Connected reduced alternating diagrams with2isomorphic Tait graphs in the same crossing/checkerboard convention are related by Conway mutations. Only this elementary Section4 implication is imported; no rank, bridge, meridional-rank or priority conclusion.'},
        {'id':'MOR14','title':'Mutant knots','authors':'H. R. Morton','year':2014,
        'url':'https://www.liverpool.ac.uk/~su14/papers/mutantchapter.pdf',
        'snapshot':'sources/morton-mutant-knots-2014.pdf',
        'location':'Author manuscript16March2014, definitions printedpp3-4 and Theorem3 with full short proof printedp6',
        'audit':'Read the torus-gluing proof: a chosen lift of each Conway halfturn is isotopic to identity on the branched-cover boundary torus. Credit Viro1976 retained. This survey supplies the audited argument, not a claim to have read the full original Russian proof or exhaustive mutation history.',
        'original_credit':'O.Ya.Viro, Nonprojecting isotopies and knots with homeomorphic coverings, Zap.Nauchn.Sem.LOMI66(1976),133-147; English J.SovietMath12(1)(1979),86-96, doi10.1007/BF01098418; primary MathNet metadata checked, original full proof not audited.',
        'scope':'Conway mutation preserves the homeomorphism type of the double branched cover, hence its fundamental group and group rank. No mutation invariance of bridge number or meridional rank is assumed.'}]
    catalogue=json.loads((HERE/'sources.json').read_text());by_id={r['id']:r for r in catalogue['sources']}
    for record in additions:
        assert (HERE/record['snapshot']).is_file()
        if record['id'] in by_id:assert by_id[record['id']]==record
        else:catalogue['sources'].append(record)
    assert len(catalogue['sources'])==len({r['id'] for r in catalogue['sources']})>=42
    (HERE/'sources.json').write_text(json.dumps(catalogue,indent=2,ensure_ascii=False)+'\n')
    files=['audit_mutation_sources.py','sources.json','sources/knotinfo/2026.11.2_mutation_subset.json',
        'sources/greene-mutation-v1.pdf','sources/morton-mutant-knots-2014.pdf']
    receipt={'status':'PASS scoped mutation source audit','catalogue_distinct_ids':len(catalogue['sources']),
        'new_records':additions,'selection':selection,
        'assurance':'Only Greene Cor4.5 graph-to-mutation and Morton Thm3/Viro cover-homeomorphism imports are used. Exact finite graph inputs are checked separately. Invariant agreement is a screening cue; no new rank or knot equality follows.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'mutation_source_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'catalogue_distinct_ids':len(catalogue['sources']),'new_ids':['GRE13','MOR14']}))


if __name__=='__main__':main()
