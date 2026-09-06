#!/usr/bin/env python3
"""Bounded exact13-adic lifting audit of possible integral model scalings."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_r17_integral13_charts as charts
from research_runtime.store import checkpoint
ROOT=charts.ROOT;ART=charts.ART;D=ROOT/'artifacts/local/elliptic-curves/r17-13-scaling-classification-v1';OUT=ART/'r17_13_scaling_classification_v1.json'

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),charts.INPUT,charts.OUT,Path(charts.__file__).resolve()]}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve bounded13-adic classification')
    if cert.read(charts.OUT)['integral_charts']!=6:raise ArithmeticError('six exact integral residue charts required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-13-scaling-classification.v1','sources':sources(),'prime':13,'maximum_depth':6,'maximum_live_residues':4096,'wall_seconds':120,'rss_bytes':1073741824,'gate':'The exact84-cell coefficient audit finds one universal integral determinant13 chart in each compact R17 family. Test whether other original residue cells can admit13^4|A and13^6|B, and whether eligible primitive pairs in each integral chart admit another13 scaling. Only exact finite residue lifting is performed.','scope':'At depth k retain all13 lifts satisfying A=0 modulo13^min(k,4) and B=0 modulo13^min(k,6). Empty levels prove local impossibility; a surviving level6 or the fixed4096-residue cap leaves the relevant classification UNKNOWN. The kernel of a determinant13 chart modulo13 is excluded only for original primitive parameter pairs. No point search, parameter score, rank predictor or automatic deeper lifting.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen13-scaling sources differ')
    return p

def evaluate(coefficients,x,modulus):
    value=0
    for c in reversed(coefficients):value=(value*x+int(c))%modulus
    return value

def lifts(a,b,initial,p):
    residues=[initial];levels=[]
    for depth in range(1,p['maximum_depth']+1):
        mod=13**depth;ma=13**min(depth,4);mb=13**min(depth,6)
        candidates=residues if depth==1 else (r+13**(depth-1)*digit for r in residues for digit in range(13))
        current=[]
        for x in candidates:
            if evaluate(a,x,ma)==0 and evaluate(b,x,mb)==0:
                current.append(x)
                if len(current)>p['maximum_live_residues']:return {'status':'UNKNOWN_RESIDUE_CAP','levels':levels,'cap_depth':depth}
        current.sort();levels.append({'depth':depth,'modulus':mod,'A_modulus':ma,'B_modulus':mb,'residues':current});residues=current
        if not residues:return {'status':'EXCLUDED_BY_EMPTY_LIFT_LEVEL','levels':levels}
    return {'status':'UNKNOWN_SURVIVING_DEPTH6','levels':levels}

def expected():
    p=protocol();proof=cert.read(charts.OUT);rows=[]
    for f in cert.read(charts.INPUT)['families']:
        a=list(map(int,f['A_coefficients_low_to_high']));b=list(map(int,f['B_coefficients_low_to_high']));universal=[r for r in proof['rows'] if r['family']==f['family'] and r['integral_after_curve_scale13']]
        if len(universal)!=1:raise ArithmeticError('one universal integral residue cell required')
        chart=universal[0];other=[]
        for residue in [*range(13),'infinity']:
            if residue==chart['residue_mod13']:continue
            result=lifts(a[::-1] if residue=='infinity' else a,b[::-1] if residue=='infinity' else b,0 if residue=='infinity' else residue,p);other.append({'residue':residue,**result})
        aa=list(map(int,chart['A_coefficients_low_to_high']));bb=list(map(int,chart['B_coefficients_low_to_high']));m=chart['matrix'];further=[];bad=[]
        for residue in [*range(13),'infinity']:
            u,v=(1,0) if residue=='infinity' else (residue,1);n=m[0]*u+m[1]*v;den=m[2]*u+m[3]*v
            if n%13==0 and den%13==0:further.append({'residue':residue,'status':'EXCLUDED_NONPRIMITIVE_OLD_PARAMETER_KERNEL'});continue
            result=lifts(aa[::-1] if residue=='infinity' else aa,bb[::-1] if residue=='infinity' else bb,0 if residue=='infinity' else residue,p);further.append({'residue':residue,**result})
            A=aa[-1]%13 if residue=='infinity' else evaluate(aa,residue,13);B=bb[-1]%13 if residue=='infinity' else evaluate(bb,residue,13);bad.append({'residue':residue,'A_mod13':A,'B_mod13':B,'discriminant_mod13':(-16*(4*A**3+27*B**2))%13})
        complete=all(r['status']=='EXCLUDED_BY_EMPTY_LIFT_LEVEL' for r in other) and all(r['status'] in ('EXCLUDED_BY_EMPTY_LIFT_LEVEL','EXCLUDED_NONPRIMITIVE_OLD_PARAMETER_KERNEL') for r in further)
        rows.append({'family':f['family'],'universal_residue_mod13':chart['residue_mod13'],'matrix':m,'other_original_cells':other,'additional_scaling_in_integral_chart':further,'classification':'EXACT_ONE_SCALE_MAXIMUM' if complete else 'UNKNOWN','scaled_discriminant_residues':bad,'all_eligible_scaled_models_remain_bad_at13':complete and all(r['discriminant_mod13']==0 for r in bad)})
    return {'schema':'elliptic-curves.r17-13-scaling-classification-result.v1','status':'PASS_BOUNDED_EXACT_AUDIT','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'complete_classifications':sum(r['classification']=='EXACT_ONE_SCALE_MAXIMUM' for r in rows),'claim_boundary':p['scope']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('exact13-adic classification differs')
        else:
            if OUT.exists():raise FileExistsError('preserve13-adic classification')
            checkpoint(OUT,d)
        print('EXACT13 SCALING AUDIT',[(r['family'],r['classification'],r['all_eligible_scaled_models_remain_bad_at13']) for r in d['rows']],flush=True)
