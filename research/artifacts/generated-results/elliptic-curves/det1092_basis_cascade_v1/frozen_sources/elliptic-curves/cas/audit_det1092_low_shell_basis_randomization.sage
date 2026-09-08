#!/usr/bin/env sage-python
"""Exact unimodular-basis sensitivity audit of the frozen low-shell schedule.

This intentionally audits the entire policy, including its SHA parity sample,
rather than merely checking the canonical-height norm under a basis change.
"""
import argparse, hashlib, json, sys
from pathlib import Path

from sage.all import ZZ, identity_matrix, matrix

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
sys.path.insert(0,str(CAS))
from research_runtime.store import digest
from importlib.machinery import SourceFileLoader

INPUT=ROOT/'artifacts/local/elliptic-curves/det1092-low-shell-cascade-v3/scale-0131232/wave-01/maps.json'
PROTOCOL=ROOT/'artifacts/local/elliptic-curves/det1092-low-shell-cascade-v3/scale-0131232/wave-01/protocol.json'
OUTPUT=ART/'det1092_low_shell_basis_randomization_v1.json'
geometry=SourceFileLoader('low_shell_basis_geometry',str(CAS/'prospective_half_lattice_v3.sage')).load_module()

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def masks(rank,floor,domain,count):
    out=[];i=0
    while len(out)<count:
        value=int(digest([domain,i]),16)%(1<<rank);i+=1
        if value>>floor and value not in out:out.append(value)
    return out
def parity(word):return sum((int(value)&1)<<index for index,value in enumerate(word))
def main(check):
    maps_data,protocol=read(INPUT),read(PROTOCOL);rank=len(maps_data['rounded_gram']);floor=protocol['rows'][0]['mask_floor'];domain=protocol['sample_domain'];count=protocol['sample_size']
    gram=matrix(ZZ,maps_data['rounded_gram']);basis_change=identity_matrix(ZZ,rank);basis_change[0,1]=1
    if abs(basis_change.det())!=1:raise ArithmeticError('test change is not unimodular')
    changed_gram=basis_change*gram*basis_change.transpose();lll=matrix(ZZ,changed_gram).LLL_gram().transpose()
    if abs(lll.det())!=1:raise ArithmeticError('randomized-basis LLL transport is not unimodular')
    reduced=lll*changed_gram*lll.transpose();oracle=geometry.CosetOracle(reduced.rows());inverse=lll.inverse();rows=[]
    for mask in masks(rank,floor,domain,count):
        residue=matrix(ZZ,1,rank,[(mask>>bit)&1 for bit in range(rank)]);reduced_residue=[int(value)%2 for value in (residue*inverse).row(0)
        ];norm,representative,error=oracle.solve(reduced_residue);new_word=(matrix(ZZ,1,rank,representative)*lll).row(0);old_word=(matrix(ZZ,1,rank,new_word)*basis_change).row(0)
        if int(matrix(ZZ,1,rank,new_word)*changed_gram*matrix(ZZ,rank,1,new_word)[0]) if False else False:pass
        new_vector=matrix(ZZ,1,rank,new_word)
        if any((int(new_word[bit])-((mask>>bit)&1))%2 for bit in range(rank)) or int((new_vector*changed_gram*new_vector.transpose())[0,0])!=norm:raise ArithmeticError('changed-basis representative identity differs')
        rows.append({'new_basis_parity':mask,'old_basis_parity':parity(old_word),'metric_norm':int(norm),'representative_in_old_basis':[int(value) for value in old_word],'cvp_error':error})
    selected=sorted(rows,key=lambda row:(row['metric_norm'],row['new_basis_parity']))[:49];baseline=[row['parity'] for row in maps_data['centres']];randomized=[row['old_basis_parity'] for row in selected];intersection=sorted(set(baseline)&set(randomized))
    payload={'schema':'elliptic-curves.det1092-low-shell-basis-randomization.v1','status':'PASS_COORDINATE_SENSITIVITY_AUDIT','inputs':{str(INPUT.relative_to(ROOT)):sha(INPUT),str(PROTOCOL.relative_to(ROOT)):sha(PROTOCOL),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},'unimodular_row_coordinate_change':[list(map(int,row)) for row in basis_change.rows()],'baseline_selected_old_basis_parities':baseline,'randomized_selected_old_basis_parities':randomized,'selected_intersection_size':len(intersection),'selected_intersection_old_basis_parities':intersection,'policy_is_basis_invariant':baseline==randomized,'interpretation':'The canonical-height coset norm is preserved by the unimodular transport. The SHA mask schedule is coordinate-indexed; a false value therefore rejects an intrinsic-policy claim and requires canonical parity addressing before such a claim can be made.','reproducing_command':'sage -python elliptic-curves/cas/audit_det1092_low_shell_basis_randomization.sage --check'}
    rendered=json.dumps(payload,indent=2,sort_keys=True)+'\n'
    if check:
        if OUTPUT.read_text()!=rendered:raise ArithmeticError('basis-randomization audit did not replay')
    else:
        if OUTPUT.exists():raise FileExistsError('preserve basis-randomization audit')
        OUTPUT.write_text(rendered)
    print('LOW-SHELL BASIS AUDIT|invariant={}|intersection={}/49'.format(payload['policy_is_basis_invariant'],len(intersection)),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();main(args.check)
