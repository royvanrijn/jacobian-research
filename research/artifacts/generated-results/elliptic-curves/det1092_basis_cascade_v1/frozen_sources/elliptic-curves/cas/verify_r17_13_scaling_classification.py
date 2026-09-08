#!/usr/bin/env python3
"""Independent exhaustive small-residue verification of the13-scaling exclusions."""
from pathlib import Path
import certify_compact_r17_candidates as cert
import classify_r17_13_scaling as audit

def value(coefficients,x,modulus):return sum(int(c)*pow(x,i,modulus) for i,c in enumerate(coefficients))%modulus

def verify_cell(a,b,record):
    residue=record['residue'];initial=0 if residue=='infinity' else residue
    if residue=='infinity':a,b=a[::-1],b[::-1]
    levels=record['levels']
    if record['status']!='EXCLUDED_BY_EMPTY_LIFT_LEVEL' or not levels or levels[-1]['residues'] or len(levels)>3:raise ArithmeticError('fixed exact empty witness by depth3 required')
    for depth,r in enumerate(levels,1):
        modulus=13**depth;ma=13**min(depth,4);mb=13**min(depth,6)
        expected=[x for x in range(initial,modulus,13) if value(a,x,ma)==0 and value(b,x,mb)==0]
        if (r['depth'],r['modulus'],r['A_modulus'],r['B_modulus'],r['residues'])!=(depth,modulus,ma,mb,expected):raise ArithmeticError('independent exhaustive residue set differs')

def main():
    p=audit.protocol();d=cert.read(audit.OUT);c=cert.read(audit.charts.OUT);families=cert.read(audit.charts.INPUT)['families'];cells=[*range(13),'infinity']
    if d['status']!='PASS_BOUNDED_EXACT_AUDIT' or d['complete_classifications']!=6 or len(d['rows'])!=6:raise ArithmeticError('complete six-family classification required')
    for f,r in zip(families,d['rows']):
        if f['family']!=r['family']:raise ArithmeticError('complete family order differs')
        originalA=list(map(int,f['A_coefficients_low_to_high']));originalB=list(map(int,f['B_coefficients_low_to_high']));q=next(x for x in c['rows'] if x['family']==r['family'] and x['integral_after_curve_scale13']);a=list(map(int,q['A_coefficients_low_to_high']));b=list(map(int,q['B_coefficients_low_to_high']));m=q['matrix']
        if r['matrix']!=m or r['universal_residue_mod13']!=q['residue_mod13'] or [x['residue'] for x in r['other_original_cells']]!=[v for v in cells if v!=q['residue_mod13']] or [x['residue'] for x in r['additional_scaling_in_integral_chart']]!=cells:raise ArithmeticError('complete projective cell roster differs')
        for cell in r['other_original_cells']:verify_cell(originalA,originalB,cell)
        expected_bad=[]
        for cell in r['additional_scaling_in_integral_chart']:
            s=cell['residue'];u,v=(1,0) if s=='infinity' else (s,1);kernel=(m[0]*u+m[1]*v)%13==0 and (m[2]*u+m[3]*v)%13==0
            if kernel:
                if cell['status']!='EXCLUDED_NONPRIMITIVE_OLD_PARAMETER_KERNEL':raise ArithmeticError('exact nonprimitive kernel differs')
            else:
                verify_cell(a,b,cell);A=a[-1]%13 if s=='infinity' else value(a,s,13);B=b[-1]%13 if s=='infinity' else value(b,s,13);delta=(-16*(4*A**3+27*B**2))%13
                expected_bad.append({'residue':s,'A_mod13':A,'B_mod13':B,'discriminant_mod13':delta})
        if expected_bad!=r['scaled_discriminant_residues'] or not all(v['discriminant_mod13']==0 for v in expected_bad) or not r['all_eligible_scaled_models_remain_bad_at13'] or r['classification']!='EXACT_ONE_SCALE_MAXIMUM':raise ArithmeticError('minimal bad13 residue claim differs')
    print('INDEPENDENT EXHAUSTIVE13-ADIC CLASSIFICATION: ALL SIX PASS; WITNESSES THROUGH13^3 ONLY',flush=True)
if __name__=='__main__':main()
