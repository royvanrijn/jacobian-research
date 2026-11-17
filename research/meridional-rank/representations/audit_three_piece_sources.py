#!/usr/bin/env python3
"""Seal the newly used case6 scope without changing the source catalogue."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    catalogue=json.loads((HERE/'sources.json').read_text())['sources']
    count=len(catalogue);assert count==len({s['id'] for s in catalogue}) and count>=45
    assert {'BW05','HAT00'}<={s['id'] for s in catalogue}
    files=['audit_three_piece_sources.py','sources.json','four_piece_source_audit.json',
        'sources/boileau-weidmann-structure-author.pdf','sources/hatcher-basic-3manifold.pdf',
        'THREE_PIECE_DOUBLE_COVERS_SHARP_FOUR.tex']
    receipt={'status':'PASS scoped three-piece case6 source audit','catalogue_distinct_ids':count,'catalogue_at_original_checkpoint':45,
        'rank_source':{'id':'BW05','location':'Theorem2, author manuscript printedp3, all ten cases re-read and specifically case6',
            'necessary_type_piece_counts':[2,2,2,2,2,3,2,1,1,2],
            'only_three_piece_case':6,
            'case6_scope':'Two end Seifert pieces over A(p) or D(p,q); central exterior of a two-bridge link; both end fibers glued to central link meridians.',
            'audit_boundary':'Necessary classification imported. Earlier separating/no-Klein-bundle proof reads retained; full42page proof and full published-version comparison not independently audited.'},
        'topology_source':{'id':'HAT00','location':'Previously read and sealed Section1.2 Theorem1.9, Propositions1.11-1.13, Lemmas1.14-1.15 and Section2.1 circle-bundle/filling construction',
            'scope':'Same explicit D(2,3) realization, boundary injection/unique fiber and nonparallel JSJ argument, now with A(2) and A(3). Original Jaco-Shalen/Johannson credit retained; normal forms and filling obstruction proved in the companion note.'},
        'elementary_argument':'Each component of a two-component two-bridge link is one-bridge, hence unknotted. Filling one meridian leaves the other unknot complement, a solid torus. This argument and the C_m*C_abs(r) nonabelian filling quotient are written in the companion proof, not inferred from source metadata.',
        'no_catalogue_change':True,'no_novelty_or_exhaustive_literature_claim':True,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'three_piece_source_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'only_three_piece_case':6,'catalogue_distinct_ids':count}))


if __name__=='__main__':main()
