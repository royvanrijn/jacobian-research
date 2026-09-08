#!/usr/bin/env sage-python
"""Bounded old-census norm6 basis for normalized Brauer evaluation.

Read at most4096 saved norm6 rows. No new lattice enumeration,
point, exceptional coordinate, Selmer or Brauer class calculation.25s cap.
"""
import csv,hashlib,json,signal
from pathlib import Path
from sage.all import ZZ,GF,matrix,vector
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_brauer_section_gate_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==s
    else:p.write_text(s)
def main():
    parentpath=ART/'curve302_recovered_mw17_parent_v1.json'
    census=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
    local=ART/'det1092_seed_local_code_v5/protocol.json'
    paths=[parentpath,census,local,Path(__file__)]
    protocol={'classification':'generic-only bounded extraction of section intersection witnesses',
      'rule':'Filter the unchanged old census by norm6, then inspect at most4096 matching rows. Keep words that increase binary rank and stop at17. Read at most the existing131072-row file; retain skipped-norm counts and every eligible-row decision. The version1 raw prefix contained no norm6 rows because the saved file is category-grouped; retain that negative input-exposure result. Reuse only the unchanged nine-address21-place local panel; no exceptional coordinate is an input.',
      'limits':{'wall_seconds':25,'eligible_saved_rows':4096,'maximum_file_rows':131072,'retained_words':17,'new_lattice_nodes':0,
        'new_addresses':0,'exceptional_coordinate_inputs':0,'point_searches':0,
        'Brauer_group_computation':0,'Selmer_groups':0,'class_groups':0,'V3_inputs':0,'pilot_changes':0},
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(OUT/'protocol.json',protocol)
    G=matrix(ZZ,read(parentpath)['generic_height_gram']);chosen=[];ledger=[];skipped={};eligible=0;scanned=0
    with census.open() as handle:
        for n,row in enumerate(csv.DictReader(handle,delimiter='\t'),1):
            assert n<=131072;scanned=n
            norm=int(row['minimum_norm'])
            if norm!=6:
                skipped[str(norm)]=skipped.get(str(norm),0)+1;continue
            eligible+=1
            if eligible>4096:break
            trial={'saved_row':n,'orbit_mask':int(row['orbit_mask']),'norm':norm}
            if norm!=6:trial['status']='OTHER_NORM'
            else:
                w=list(map(int,row['parent_MW17_w'].split()));v=vector(ZZ,w);assert v*G*v==6
                rank=matrix(GF(2),[r['word'] for r in chosen]+[w]).rank()
                trial['status']='DEPENDENT_MOD2'
                if rank>len(chosen):
                    trial['status']='RETAINED';chosen.append({'saved_row':n,'orbit_mask':int(row['orbit_mask']),
                      'word':w,'height':6,'intersection_with_O':1})
            ledger.append(trial)
            if len(chosen)==17:break
    rank=len(chosen);det=int(matrix(ZZ,[r['word'] for r in chosen]).det()) if rank==17 else None
    result={'status':'NORM6_ODD_INDEX_SECTION_FRAME' if rank==17 else 'UNKNOWN_WITHIN_ROW_LIMIT',
      'classification':'exact lattice and section-intersection certificate','rank_mod2':rank,
      'integer_determinant':det,'chosen':chosen,'ledger':ledger,'scanned_file_rows':scanned,'skipped_norm_counts':skipped,
      'height_formula':'For a nonzero section on the proved24I1 K3, height=4+2(S.O); norm6 therefore gives a degree-one rational intersection with O.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [parentpath,census,OUT/'protocol.json']}}
    retain(OUT/'frame.json',result)
    print(result['status'],'eligible rows',len(ledger),'file rows',scanned,'rank',rank,'determinant',det,flush=True)
if __name__=='__main__':
    signal.alarm(25);OUT.mkdir(parents=True,exist_ok=True);main()
