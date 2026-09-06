#!/usr/bin/env sage-python
"""Freeze49 rational charts from each previously audited65536-mask sample."""
import argparse,sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import larger_specialized_parity_trial as control
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
mapper=SourceFileLoader('larger_trial_mapper',str(CAS/'prepare_fresh_r17_pari_batch.sage')).load_module()
def main(index):
    control.configure(index);p=control.protocol();row=control.ROW;seed=cert.read(control.SEED);old=cert.read(ROOT/row['old_maps_path']);sample=cert.read(ROOT/row['sample_path']);out=control.D/'maps.json'
    if out.exists():raise FileExistsError('preserve larger sample maps')
    if sample['status']!='COMPLETE_FIXED_SAMPLE' or len(sample['sample'])!=65536 or len(sample['centres'])!=49:raise ArithmeticError('complete audited sample required')
    model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    data={'status':'RUNNING_MAPS','protocol_hash':digest(p),'metric_gram':old['metric_gram'],
        'sample':sample['sample'],'centres':sample['centres'],'rows':[]};checkpoint(out,data)
    mapper.pari.allocatemem(256000000,silent=True)
    for c in data['centres']:data['rows'].append(mapper.mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data);print('LARGER PARITY FROZEN49 MAPS',row['id'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);main(p.parse_args().index)
