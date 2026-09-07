#!/usr/bin/env sage-python
"""Freeze301 exact PARI maps using only the newly discovered25-point subgroup."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import followup_higher_rank25 as follow
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
mapper=SourceFileLoader('paired_adaptive_mapper',str(CAS/'prepare_fresh_r17_pari_batch.sage')).load_module()

def main():
    p=follow.protocol();seed=cert.read(follow.SEED);out=follow.D/'maps.json'
    if out.exists():raise FileExistsError('preserve adaptive maps')
    model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['final_state']['state']['reductions']['points']);gram,asymmetry=mapper.geometry.canonical_height_gram(model,points);oracle=mapper.geometry.CosetOracle(mapper.geometry.rounded_gram(gram,1000000));words=sorted(range(1,256),key=lambda w:(w.bit_count(),w));centres=[]
    for i,g in enumerate(p['generic_pool']):
        word=words[i%len(words)];mask=g['mask']|(word<<17);residue=[(mask>>j)&1 for j in range(25)];norm,rep,error=oracle.solve(residue)
        if len(rep)!=25 or any((rep[j]-residue[j])%2 for j in range(25)):raise ArithmeticError('adaptive parity mismatch')
        centres.append({'generic_mask':g['mask'],'quotient_word':word,'parity':mask,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
    centres.sort(key=lambda c:(-c['metric_norm'],c['generic_mask'],c['quotient_word']));data={'status':'RUNNING','protocol_hash':digest(p),'centres':centres,'metric_gram':[[str(x) for x in row] for row in gram],'maximum_gram_asymmetry':str(asymmetry),'rows':[]};checkpoint(out,data);mapper.pari.allocatemem(256000000,silent=True)
    for c in centres:data['rows'].append(mapper.mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data);print('FROZEN NEW25 ADAPTIVE MAPS',len(data['rows']),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);a=p.parse_args();follow.configure(a.id);main()
