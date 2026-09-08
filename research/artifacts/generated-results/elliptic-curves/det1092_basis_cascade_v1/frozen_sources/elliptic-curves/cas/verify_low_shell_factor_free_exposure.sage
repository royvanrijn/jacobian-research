#!/usr/bin/env sage-python
"""Exact verifier for a frozen factor-free low-shell parity exposure."""
import argparse
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from types import SimpleNamespace

from sage.all import ZZ, matrix

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from replay_retention24_geometry import cloud_check, geometry, tuples
from research_runtime.store import digest


def check(directory, protocol_path, cloud_path):
    protocol=cert.read(protocol_path);rank=len(cert.read(directory/'seed.json')['points']);control=SimpleNamespace(SEED=directory/'seed.json')
    def masks():
        values=[];counter=0;floor=next(row['mask_floor'] for row in protocol['rows'] if row['id']==directory.name)
        while len(values)<protocol['sample_size']:
            mask=int(digest([protocol['sample_domain'],counter]),16)%(1<<rank);counter+=1
            if mask>>floor and mask not in values:values.append(mask)
        return values
    maps,data,seed=cert.read(directory/'maps.json'),cert.read(directory/'result.json'),cert.read(control.SEED);initial=tuples(seed['points'])
    if maps['protocol_hash']!=digest(protocol) or data['protocol_hash']!=digest(protocol) or len(initial)!=rank or tuples(data['initial_state']['state']['reductions']['points'])!=initial:raise ArithmeticError('low-shell subgroup/protocol differs')
    g=matrix(ZZ,maps['rounded_gram']);u=matrix(ZZ,maps['change_of_basis']);h=matrix(ZZ,maps['reduced_gram'])
    with localcontext() as context:
        context.prec=110;expected=[[int((Decimal(value)*1000000).to_integral_value()) for value in row] for row in maps['metric_gram']]
    if g!=matrix(ZZ,expected) or not g.is_symmetric() or not g.is_positive_definite() or abs(u.det())!=1 or h!=u*g*u.transpose():raise ArithmeticError('low-shell metric transport differs')
    if [row['parity'] for row in maps['sample']]!=masks():raise ArithmeticError('low-shell SHA parity sample differs')
    for row in maps['sample']:
        word=matrix(ZZ,1,rank,row['representative']);reduced=matrix(ZZ,1,rank,row['reduced_representative'])
        if word!=reduced*u or int((word*g*word.transpose())[0,0])!=row['metric_norm'] or any((int(word[0,index])-((row['parity']>>index)&1))%2 for index in range(rank)):raise ArithmeticError('low-shell parity/norm identity differs')
    selected=sorted(maps['sample'],key=lambda row:(row['metric_norm'],row['parity']))[:49]
    if maps['centres']!=selected or [row['centre'] for row in maps['rows']]!=selected:raise ArithmeticError('low-shell selected roster differs')
    geometry(data,maps,initial,{**protocol,'maps_path':directory/'maps.json'});cloud_check(data,data['charts'],cloud_path,directory/'result.json')
    if data['status']!='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' or len(data['charts'])!=49:raise ArithmeticError('low-shell worker did not complete')
    print('EXACT LOW-SHELL 2048 PARITIES;49 SELECTED;RATIONAL MAP/POINT PROVENANCE',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True);parser.add_argument('--protocol',type=Path,required=True);parser.add_argument('--cloud',type=Path,required=True);args=parser.parse_args();check(args.run,args.protocol,args.cloud)
