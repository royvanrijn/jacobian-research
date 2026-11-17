#!/usr/bin/env python3
"""Seal reuse of the precise BW05 case6/Seifert scopes for five direct models."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    catalogue=json.loads((HERE/'sources.json').read_text())['sources']
    previous=json.loads((HERE/'three_piece_source_audit.json').read_text())
    count=len(catalogue);assert count==len({s['id'] for s in catalogue}) and count>=45
    assert previous['rank_source']['only_three_piece_case']==6
    files=['audit_final_five_sources.py','sources.json','three_piece_source_audit.json',
        'sources/boileau-weidmann-structure-author.pdf','sources/hatcher-basic-3manifold.pdf',
        'FINAL_FIVE_DOUBLE_COVERS_SHARP_FOUR.tex']
    report={'status':'PASS scoped final-five classical source reuse','catalogue_distinct_ids':count,'catalogue_at_original_checkpoint':45,'new_sources':0,
        'rank_input':'BW05 Theorem2, author printedp3, only three-piece type case6 requires central two-bridge link and adjoining end-fiber meridians. Necessary classification imported; full proof/published-version comparison not independently audited.',
        'topology_input':'Hatcher Sections1.2-2.1, retained Theorem1.9/Propositions1.11-1.13/Lemmas1.14-1.15 and circle-bundle construction. Negative-base boundary uniqueness and nonparallel canonical JSJ reused. Original JSJ credit retained.',
        'new_written_checks':'Direct source-specific Q_R(2/3),Q_L,Q_S inverse eliminations, conjugations, ks versus sk boundary base paths, one negative inner-fiber marking, exact source assembly and central filling quotients. Two-component two-bridge unknot/solid-torus filling argument written explicitly.',
        'scope':'Five direct new original-target models; no graph transfer, source-cover homeomorphism, rank=genus or global classification result claimed. Working proofs unformalized and externally unreviewed; fullMRC/exhaustive literature/worldwide novelty remain open.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'final_five_source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'catalogue_distinct_ids':count,'new_sources':0}))


if __name__=='__main__':main()
