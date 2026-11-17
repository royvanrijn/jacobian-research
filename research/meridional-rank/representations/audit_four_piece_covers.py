#!/usr/bin/env python3
"""Retain exact proposed peripheral models; topology/rank require written imports."""
import hashlib
import json
from pathlib import Path
from compact_cover_groups import inv,free,substitute,cyclic

HERE=Path(__file__).resolve().parent


def power(w,n):return (w if n>=0 else inv(w))*abs(n)


def main():
    plan=json.loads((HERE/'four_piece_cover_plan.json').read_text())
    assemblies=json.loads((HERE/'port_assembly_results.json').read_text())['targets']
    files=['four_piece_cover_plan.json','audit_four_piece_covers.py','port_assembly_results.json',
        'port_assembly_verification.json','compact_cover_verification.json','mutation_cover_verification.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},'models':[],'roots':[]}
    # The common model alphabet is a=1,t=2,s=3, with t^2=a^3 and s^2=d.
    P=assemblies[0]['pieces'][0]['relators'][0]
    Q=[[x-2 if x>0 else x+2 for x in assemblies[0]['pieces'][1]['relators'][0]]][0]
    P_forward={1:[1],2:[-2,1,1]};P_back={1:[1],2:[1,1,-2]}
    assert all(substitute(P_forward[i],P_back)==[i] and substitute(P_back[i],P_forward)==[i] for i in (1,2))
    tP=[-1,-1,2,2];targetP=free(tP+tP+[-1]*3)
    tQ=[2,2,1,1];targetQ=free(tQ+tQ+[-1]*3)
    assert cyclic(substitute(P,P_forward))==cyclic(targetP)
    assert cyclic(Q)==cyclic(targetQ)
    report['models']=[
        {'name':'P','source_relator':P,'free_basis_forward':P_forward,'free_basis_inverse':P_back,
            'power_relator_two_generator':targetP,'old_generators_in_amalgam':{1:[1],2:[-3,1,1]},
            'root_word_d':[1,1,2],'H_peripheral_e':[-1,2],'H_fiber_h':[1,1,1],
            'root_abelian_exponent':7,'root_boundary_matrix_in_e_h':[[0,1],[1,1]],
            'orbifold_bases':['D(2,3)','A(2)']},
        {'name':'Q','source_relator':Q,'power_relator_two_generator':targetQ,
            'old_generators_in_amalgam':{1:[1],2:[3]},'root_word_d':[2,-1,-1],
            'H_peripheral_e':[2,1],'H_fiber_h':[1,1,1],'root_abelian_exponent':-1,
            'orbifold_bases':['D(2,3)','A(2)']}]
    for assembly,epsilon in zip(assemblies,[1,-1]):
        pports=[];qports=[]
        for port in assembly['ports']:
            p=substitute(port['inner_word'],{1:[1],2:[-3,1,1]})
            q=substitute(port['outer_word'],{3:[1],4:[3]})
            pports.append({'edge':port['edge'],'word':p});qports.append({'edge':port['edge'],'word':q})
        h=[1,1,1];aP=[3]+h;f=[3,3];aQ=[3]+power(h,epsilon)
        expectedP=[inv(f),aP,free(aP+inv(f)),[]]
        expectedQ=[[],f,free(inv(aQ)+f),aQ]
        # Columns are the source ordered peripheral basis expressed in the target basis.
        matrices=[[[0,1],[1,1]],[[0,1],[1,1]],[[0,1],[epsilon,-1]]]
        row={'name':assembly['name'],'epsilon':epsilon,'P_ports_in_model':pports,'Q_ports_in_model':qports,
            'P_port_expected':expectedP,'Q_port_expected':expectedQ,
            'P_outer_basis_a_f':[aP,f],'Q_outer_basis_a_f':[aQ,f],
            'gluing_matrices':matrices,'fiber_intersections':[1,1,1],
            'JSJ_bases_proposed':['D(2,3)','A(2)','A(2)','D(2,3)'],
            'orbifold_Euler_characteristics':['-1/6','-1/2','-1/2','-1/6'],
            'abelian_gluing_matrix':[[14,-(1-12*epsilon)],[19,-(1+12*epsilon)]],
            'expected_H1_order':391 if epsilon==1 else 401,
            'assurance':'Finite algebraic/peripheral inputs only; constructed manifold realization, irreducibility/JSJ, rank2 exclusion and knot conclusion require independently reviewed written arguments and classical imports.'}
        report['roots'].append(row)
    report['status']='COMPLETE_PERIPHERAL_MODEL_INPUTS'
    (HERE/'four_piece_cover_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'roots':2,'formal_pieces':4,'claimed_rank_from_finite_inputs_alone':False}))


if __name__=='__main__':main()
