#!/usr/bin/env python3
"""One bounded coordinate repair: t=1/u moves curve398's bad infinity to zero.

Protocol: same equation and retained witness, no new parameter. Reverse the
degree-(8,12) coefficient arrays, verify the new infinity is good, then reuse
the original 30-second geometry worker once. Preserve the unsupported v1 row.
No retries after this completion worker. No elliptic points are inspected.
"""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import generic_selmer_capacity as source
import fixed_field_transfer_geometry as geometry

INPUT=r.OUT/'rank_jump_generic_selmer_geometry_completion_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_generic_selmer_geometry_completion_v1.json'
FAMILY='curve398-p16875'


def export():
    f=next(x for x in r.read(source.INPUT)['families'] if x['family']==FAMILY)
    assert len(f['A'])==9 and len(f['B'])==13 and r.F(f['irreducibility_witness_parameter'])
    row={**f,'A':f['A'][::-1],'B':f['B'][::-1],
        'irreducibility_witness_parameter':str(1/r.F(f['irreducibility_witness_parameter']))}
    r.write_new(INPUT,{'families':[row],'original_base_coordinate':'t=1/u',
        'bindings':source.bind([source.INPUT,source.OUTPUT,Path(__file__),source.PROTOCOL])})


def worker():
    geometry.INPUT=INPUT;geometry.PROTOCOL=source.PROTOCOL
    row=geometry.worker(FAMILY)
    row['bindings']=source.bind([INPUT,Path(__file__),Path(geometry.__file__),source.PROTOCOL])
    r.write_new(OUTPUT,{'rows':[row],'bindings':row['bindings'],
        'boundary':'Exact base inversion of the same curve398 family; original geometry row remains UNKNOWN because that worker required good infinity.'})


def capture():
    with (source.WORK/'curve398-base-inversion.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=30)
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded worker timeout'
    if error:r.write_new(OUTPUT,{'rows':[{'family':FAMILY,'status':'UNKNOWN','reason':error}]})
    print(r.read(OUTPUT)['rows'][0]['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);args=p.parse_args();globals()[args.mode]()
