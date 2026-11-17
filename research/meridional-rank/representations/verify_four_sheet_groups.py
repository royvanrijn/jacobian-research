#!/usr/bin/env python3
"""Independent planar graph and exact triviality proof replay for all21 covers."""
from collections import Counter
from itertools import combinations
import hashlib
import json
from pathlib import Path
from parabolic_lifts import oriented_diagram  # index transport only
from verify_three_sheet_splices import planar_arcs
from verify_three_sheet_structure import canonical
from verify_three_sheet_nonabelian import free,cyclic,inverse

HERE=Path(__file__).resolve().parent


def main():
    report=json.loads((HERE/'four_sheet_group_results.json').read_text());assert report['status']=='COMPLETE'
    for p,h in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    old=json.loads((HERE/'small_symmetric_full_panel_results.json').read_text())
    expected=[r['name'] for r in old['results'] if r['degrees'][2]['certificates']]
    assert expected==[r['name'] for r in report['results']] and len(expected)==21
    edges=list(combinations(range(4),2));summaries=[]
    for row in report['results']:
        name=row['name'];pd=panel[name]['pd'];new_cr,_=planar_arcs(pd);transport=[c[2] for c in oriented_diagram(pd)]
        labels=[row['arc_labels'][i] for i in transport];tree=row['tree_generators'];assert len(set(tree))==3
        tree_vertices={0}
        for _ in range(3):
            for i in tree:
                a,b=edges[row['arc_labels'][i-1]]
                if a in tree_vertices or b in tree_vertices:tree_vertices.update([a,b])
        assert tree_vertices==set(range(4));words=[]
        for o,a,b in new_cr:
            for start in range(4):
                sheet=start;word=[]
                for arc in [o,a,o,b]:
                    p,q=edges[labels[arc]];name_i=transport[arc]+1
                    if sheet not in [p,q]:continue
                    if name_i not in tree:word.append(name_i if sheet==p else -name_i)
                    sheet=q if sheet==p else p
                assert sheet==start;words.append(word)
        assert Counter(tuple(canonical(w)) for w in words)==Counter(tuple(canonical(w)) for w in row['relators'])
        assert set(row['generators'])==set(range(1,14))-set(tree)
        current=[cyclic(w) for w in row['relators'] if w];active=set(row['generators'])
        for step in row['eliminations']:
            i=step['relator_index'];r=current[i];g=step['generator'];replacement=step['replacement']
            assert r==step['relator'] and sum(abs(a)==g for a in r)==1 and g in active
            assert all(abs(a) in active-{g} for a in replacement)
            def substitute(w):return free(a for letter in w for a in (replacement if letter==g else inverse(replacement) if letter==-g else [letter]))
            assert substitute(r)==[];current.pop(i)
            current=[v for w in current if (v:=cyclic(substitute(w)))];active.remove(g)
        assert not active and not current
        assert row['remaining_generators']==row['remaining_relators']==[] and row['status']=='PROVED_TRIVIAL_BY_ELIMINATION'
        summaries.append({'name':name,'cover_group':'trivial','Tietze_eliminations':len(row['eliminations'])})
    paths=['verify_four_sheet_groups.py','four_sheet_group_plan.json','four_sheet_group_results.json','joint_a5_panel.json',
           'small_symmetric_full_panel_results.json','verify_three_sheet_splices.py','verify_three_sheet_structure.py','verify_three_sheet_nonabelian.py','parabolic_lifts.py']
    receipt={'status':'PASS exact four-sheet group triviality','covers':21,'all_groups_trivial':True,'results':summaries,
             'new_sharp_cases':0,'scope':'Complete fixed21 S4-transposition covers on the293-name panel. Independent under-to-under planar arc and folded branch-word multiset reconstruction; every generator elimination is replayed. No knot3bridge presentation follows. Source parser is used only for certificate index transport.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'four_sheet_group_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'covers':21,'all_groups_trivial':True}))


if __name__=='__main__':main()
