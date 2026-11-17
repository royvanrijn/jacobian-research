#!/usr/bin/env python3
"""Frozen exact three-piece inputs and five small graph-transfer searches."""
import hashlib
import json
from pathlib import Path
import time
from compact_cover_groups import free,inv,substitute
from mutation_graphs import checkerboards,search

HERE=Path(__file__).resolve().parent


def power(w,n):return (w if n>=0 else inv(w))*abs(n)


def main():
    plan=json.loads((HERE/'three_piece_cover_plan.json').read_text());start=time.monotonic()
    structure={r['name']:r for r in json.loads((HERE/'remaining_cover_structure_results.json').read_text())['results']}
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    files=['three_piece_cover_plan.json','audit_three_piece_covers.py','remaining_cover_structure_results.json',
        'compact_cover_groups.py','mutation_graphs.py','parabolic_panel.json','four_piece_cover_verification.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},'roots':[],'diagrams':{},'pairs':[]}
    specs=[(2,[1,-2],{1:[3],2:[1]},[1],1,[[2,3],[-1,-2]],403),
        (3,[2,-1],{1:[3],2:[1]},[], -1,[[-1,-2],[0,-1]],361),
        (2,[2,-1,-1],{1:[1],2:[3]},[],1,[[-1,-3],[1,2]],337)]
    h=[1]*3;t=[2];Pmap={1:[1],2:[-1,-1,2]};e=free(Pmap[2]+inv(h))
    for name,spec in zip(plan['roots'],specs):
        m,d,Qmap,conjugator,eta,matrix,D=spec;cut=structure[name]['cuts'][0]
        pieces=[];portwords=[]
        for frag in cut['fragments']:
            reduction=frag['formal_core_reduction'];rename={g:[i+1] for i,g in enumerate(reduction['final']['generators'])}
            words={int(g):substitute(w,rename) for g,w in reduction['original_generators_in_final'].items()};words[1]=[]
            pieces.append({'relator':substitute(reduction['final']['relators'][0],rename),'arc_words':words})
        for edge in cut['boundary']:
            words=[]
            for frag,piece in zip(cut['fragments'],pieces):
                arc=next(i+1 for i,edges in enumerate(frag['arc_edge_classes']) if edge in edges)
                words.append(piece['arc_words'][arc])
            portwords.append({'edge':edge,'P':words[0],'Q':words[1]})
        f=[3]*m;a=free(power([3],eta)+h)
        if name=='13a_233':
            expectedP=[h,free(e+h),free(power(h,2)+inv(e)),power(e,2)]
            expectedQ=[[],free(inv(f)+power(a,2)),free(inv(f)+a),a]
        elif name=='13a_239':
            expectedP=[h,e,free(h+inv(e)),[]]
            expectedQ=[[],free(f+a),a,free(f+power(a,2))]
        else:
            expectedP=[free(h+power(e,-2)),e,free(h+inv(e)),[]]
            expectedQ=[[],f,free(inv(a)+f),a]
        report['roots'].append({'name':name,'boundary':cut['boundary'],'pieces':pieces,'port_words':portwords,
            'P_images':Pmap,'Q_images':Qmap,'Q_conjugator':conjugator,'root_degree':m,'root_word_d':d,'central_outer_eta':eta,
            'P_ports_expected':expectedP,'Q_ports_expected':expectedQ,'central_outer_a_f':[a,f],
            'P_to_C_matrix_e_h_in_a_f':matrix,'C_to_Q_matrix_k_f_in_e_h':[[0,1],[1,-1 if name=='13a_253' else 0]],
            'expected_H1_order':D,'fiber_intersections':[abs(matrix[0][1]),1],
            'proposed_JSJ_bases':['D(2,3)',f'A({m})','D(2,3)']})
    path=HERE/'three_piece_cover_results.json'
    for names in plan['blocks'].values():
        for name in names:report['diagrams'][name]={'pd':panel[name]['pd'],**checkerboards(panel[name]['pd'])}
        for name in names[1:]:
            assert time.monotonic()-start<plan['limits']['seconds_total']
            A=report['diagrams'][names[0]]['selected'];B=report['diagrams'][name]['selected']
            certificate=search(A,B,plan['limits'])
            report['pairs'].append({'source':names[0],'target':name,'certificate':certificate})
            path.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'source':names[0],'target':name,'status':certificate['status'],'switches':len(certificate.get('switches',[]))}),flush=True)
    report['status']='COMPLETE_PERIPHERAL_AND_TRANSFER_INPUTS' if all(p['certificate']['status'].startswith('CERTIFIED') for p in report['pairs']) else 'PARTIAL_TRANSFER_INPUTS'
    report['seconds']=round(time.monotonic()-start,6)
    report['assurance']='Finite proposed algebra/peripheral inputs and bounded graph certificates only. Independent checks and written realization/JSJ/classification/filling argument required before a knot conclusion.'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'roots':len(report['roots']),'pairs':len(report['pairs']),'seconds':report['seconds']}))


if __name__=='__main__':main()
