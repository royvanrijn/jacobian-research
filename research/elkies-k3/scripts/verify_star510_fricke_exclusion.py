#!/usr/bin/env python3
"""Check finite table premises for a sourced lower-level Fricke exclusion."""
import argparse,hashlib,json,subprocess
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'artifacts/generated-results/elkies-k3-star510-fricke-exclusion-v1'
def run():
 t=S.symbols('t')
 rows=[{'quadratic_radicand':-2,'j':'8000','reported_CM_discriminant':-8},
       {'quadratic_radicand':-2,'j':'8000','reported_CM_discriminant':-8},
       {'quadratic_radicand':17,'j':'-671956992*t-2770550784','reported_CM_discriminant':-51}]
 allowed=list(range(1,20))+[21,25,27,37,43,67,163]
 assert 51 not in allowed and 510%51==0 and 510//51==10
 j=-671956992*t-2770550784;jpol=S.expand((S.Symbol('J')-j)*(S.Symbol('J')-j.subs(t,-t)));jpol=S.rem(jpol,t*t-17,t)
 assert jpol==S.Symbol('J')**2+5541101568*S.Symbol('J')+6262062317568
 code='if(polclass(-8)!=x-8000,error("D8"));if(polclass(-51)!=x^2+5541101568*x+6262062317568,error("D51"));print("CM_PASS");quit;\n'
 p=subprocess.run(['gp','-q','-f'],input=code,text=True,capture_output=True,timeout=10,check=True)
 assert p.stdout.strip()=='CM_PASS',p.stderr
 return {'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'external_quadratic_classification':'Ozman-Siksek arXiv:1806.08192v3, Main Theorem and Table8.6','table_representatives':rows,'rational_isogeny_degrees':allowed,'lower_level':51,'target_level':510,'cyclic_subgroup_order':51,'quadratic_nonCM_points_on_X0_51':0,'rational_nonCM_points_on_X0_51':0,'nonCM_rational_Fricke_lifts_at_level510':0,'star510_candidate_status':'EXCLUDED for the rational Fricke lift route','classification_recomputed':False,'positive_correlated_target_complete':False}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
 if args.write:
  PACKET.mkdir(exist_ok=True);(PACKET/'result.json').write_text(json.dumps(out,indent=2)+'\n')
 else:assert out==json.loads((PACKET/'result.json').read_text())
 print('PASS: CM polynomials and level51 divisor gate; star510 rational Fricke route excluded using sourced classification')
