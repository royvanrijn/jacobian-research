#!/usr/bin/env python3
"""Retain five frozen proposed peripheral models without numerical recognition."""
import hashlib
import json
from pathlib import Path
from compact_cover_groups import free,inv,substitute

HERE=Path(__file__).resolve().parent


def power(w,n):return (w if n>=0 else inv(w))*abs(n)


def main():
    plan=json.loads((HERE/'final_five_peripheral_plan.json').read_text())
    sources={r['name']:r for r in json.loads((HERE/'final_five_structure_results.json').read_text())['results']}
    files=['final_five_peripheral_plan.json','audit_final_five_peripherals.py','final_five_structure_results.json',
        'compact_cover_groups.py','three_piece_cover_verification.json','three_piece_source_audit.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},'targets':[]}
    h=[1]*3;Pmap={1:[1],2:[-1,-1,2]};e=free(Pmap[2]+inv(h))
    Pcoords={233:[[0,1],[1,1],[-1,2],[2,0]],239:[[0,1],[1,0],[-1,1],[0,0]],253:[[-2,1],[1,0],[-1,1],[0,0]]}
    for spec in plan['models']:
        name=spec['name'];cut=sources[name]['cuts'][0];pieces=[];ports=[]
        for frag in cut['fragments']:
            reduction=frag['formal_core_reduction'];rename={g:[i+1] for i,g in enumerate(reduction['final']['generators'])}
            words={int(g):substitute(w,rename) for g,w in reduction['original_generators_in_final'].items()};words[1]=[]
            pieces.append({'relator':substitute(reduction['final']['relators'][0],rename),'arc_words':words})
        for edge in cut['boundary']:
            words=[]
            for fragment,piece in zip(cut['fragments'],pieces):
                arc=next(i+1 for i,g in enumerate(fragment['arc_edge_classes']) if edge in g);words.append(piece['arc_words'][arc])
            ports.append({'edge':edge,'P':words[0],'Q':words[1]})
        k=power(h,spec['k_h_sign']);f=[3]*spec['m'];c=free(k+[3] if spec['outer_order']=='ks' else [3]+k)
        pc=Pcoords[spec['P_type']];qc=plan['proposed_Q_ports_c_f'][name]
        expectedP=[free(power(e,A)+power(h,B)) for A,B in pc]
        expectedQ=[free(power(c,A)+power(f,B)) for A,B in qc]
        M=spec['P_to_C'];N=[[0,1],[spec['k_h_sign'],-1 if spec['Q_type']=='root_squared' else 0]]
        base=0 if spec['P_type']==233 else 3;z=qc[base] if base==3 else [-M[0][1],-M[1][1]]
        row={**spec,'boundary':cut['boundary'],'pieces':pieces,'ports':ports,
            'P_images':Pmap,'outer_c_f':[c,f],'P_port_coordinates_e_h':pc,'Q_port_coordinates_c_f':qc,
            'P_ports_expected':expectedP,'Q_ports_expected':expectedQ,'C_to_Q':N,'z_in_c_f':z,
            'fiber_intersections':[abs(M[0][1]),1],'proposed_JSJ_bases':['D(2,3)',f'A({spec["m"]})','D(2,3)']}
        report['targets'].append(row)
        (HERE/'final_five_peripheral_results.json').write_text(json.dumps(report,indent=2)+'\n')
    report['status']='COMPLETE_FROZEN_PERIPHERAL_INPUTS'
    report['assurance']='Finite proposed exact models only; independent replay plus written source-group/manifold/JSJ/classification/filling argument required before a knot closure. No mutation or common-cover transfer is used.'
    (HERE/'final_five_peripheral_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'models':len(report['targets'])}))


if __name__=='__main__':main()
