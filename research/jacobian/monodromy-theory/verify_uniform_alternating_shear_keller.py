#!/usr/bin/env python3
"""Exact checks for the primary A12 model and uniform odd-shear Keller family.

The companion note proves all odd m and full geometric monodromy. Bounded
controls supplement those proofs. --full-jacobian additionally expands the
actual degree(33,159,3) three-variable A12 map and differentiates it directly.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from math import comb
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
from alternating_shear_keller_model import (
 binomial_terms,coefficient_pair,expanded_pair,expanded_map,laurent_primitive)

HERE=Path(__file__).resolve().parent
INPUT=HERE/'uniform_alternating_shear_keller_inputs.json'


def require(value,message):
 if not value:raise AssertionError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def literal_A12(data):
 a=data['A12'];c,d=coefficient_pair(3)
 expected_c=s.Rational(*a['c_real'])+s.I*s.Rational(*a['c_imag'])
 expected_d=s.Rational(*a['d_real'])+s.I*s.Rational(*a['d_imag'])
 require(s.simplify(c-expected_c)==0 and s.simplify(d-expected_d)==0,'literal primary coefficient field values')
 primary={tuple(row[1:]):s.Rational(*row[0])for row in a['B_terms']}
 formula={tuple(row[1:]):row[0]for row in binomial_terms(3)}
 require(primary==formula,'generic binomial formula matches all twelve literal primary B terms')
 pair=expanded_pair(3);L,w=pair['Lambda'],pair['w'];A,B=pair['A'],pair['B']
 require(A.total_degree()==a['A_degree']==17 and B.total_degree()==a['B_degree']==80,'actual A12 coefficient degrees')
 require(A.degree(w)==a['A_w_degree']==12 and B.degree(w)==a['B_w_degree']==54,'actual w degrees')
 one=s.Poly(1,L,w,domain=s.QQ_I)
 master=s.Poly(L,L,w,domain=s.QQ_I)*(B.diff(L)*A.diff(w)-B.diff(w)*A.diff(L))+B*A.diff(w)-A*B.diff(w)
 require(master==one,'expanded exact QQ(i) master Jacobian identity')
 require(A.eval({L:1,w:-1})==0 and B.eval({L:1,w:-1})==0,'both lift-point conditions')
 require(A.eval(L,0)==s.Poly(w,w,domain=s.QQ_I)and B.eval(L,0)==s.Poly(1,w,domain=s.QQ_I),
         'zero-Lambda boundary coefficients')
 # The critical value polynomial provides an independent branch-value collision control.
 p=1;m=3;Z=s.symbols('Z');Xi=s.Integer(1)
 for j in range(m+1):Xi-=s.Rational(2*comb(m,j)*(-1)**j,2*j-1)*Z**j*(1+Z)**(m-j)
 for j in range(2*m+1):Xi+=s.Rational(comb(2*m,j)*(-1)**j,2*j-3)*Z**j*(1+Z)**(2*m-j)
 theta=sum(s.Integer(v)*Z**j for j,v in enumerate(a['critical_theta_coefficients']))
 require(s.expand(Xi-s.Rational(8,315)*theta)==0,'literal A12 critical-curve coefficient identity')
 require(theta.subs(Z,0)==105 and theta.subs(Z,-1)==28,'critical origin and other Y-zero values nonzero')
 return {'literal_primary_B_terms':12,'A_terms':len(A.terms()),'B_terms':len(B.terms()),
  'A_B_degrees':[17,80],'A_B_w_degrees':[12,54],
  'expanded_master_identity_over_Qi':True,'polynomial_lift_conditions':True,
  'critical_branch_curve_birationality_controls':True,
  'geometric_monodromy_by_written_proof':'A12','coefficient_field':'Q(i)',
  'arithmetic_group_over_Qi_by_discriminant_proof':'S12','constant_field':'Q(i,sqrt(3))'}


def uniform_controls():
 checked=[]
 for m in range(3,16,2):
  p=(m-1)//2;N=m*(m+1);model=laurent_primitive(m)
  X,Y,T,c,d,P,Phi=[model[k]for k in('X','Y','T','c','d','P','Phi')]
  require(s.expand(X**(m+1)*s.diff(P,X)-Phi**2)==0,'uniform squared derivative identity')
  f=s.Poly(s.expand(X**m*(P-T)),X)
  require(f.degree()==N,'degree from the two pole orders')
  require(s.expand(X*f.diff().as_expr()-m*f.as_expr()-Phi**2)==0,'inverse polynomial derivative identity')
  aN=c*c*d**(2*m)/m**2;a0=-c*c*Y**(2*m)/m
  require(s.simplify(f.LC()-aN)==0 and s.simplify(f.TC()-a0)==0,'leading and constant coefficients')
  D=(-1)**((m-1)//2)*m
  require(s.simplify(((-1)**(N//2)*aN/a0)/D-(d**m/(m*Y**m))**2)==0,
          'constant discriminant square class, with no Y twist')
  for a,b,j,Lpow,Cpow,Wpow in binomial_terms(m):
   require(min(Lpow,Cpow,Wpow)>=0,'all substituted B monomials polynomial')
  require(all((p+1)*j-p!=0 for j in range(m+1))and
          all((p+1)*j-m!=0 for j in range(2*m+1)),'all Laurent primitive denominators nonzero')
  Z=s.symbols('Z');Xi=s.Integer(1)
  for j in range(m+1):Xi-=s.Rational(2*comb(m,j)*(-1)**j,(p+1)*j-p)*Z**j*(1+Z)**(m-j)
  for j in range(2*m+1):Xi+=s.Rational(comb(2*m,j)*(-1)**j,(p+1)*j-m)*Z**j*(1+Z)**(2*m-j)
  require(s.simplify(Xi.subs(Z,0)-s.Rational(2*(p+1)**2,p*m))==0,'critical parametrization order at origin')
  require(s.simplify(Xi.subs(Z,-1)-s.Rational((m*m-1)**2,m*m*(m*m+1)))==0,'unique critical-curve zero fibre')
  require(m-2*(p+1)==-1,'coprime critical valuations')
  checked.append({'m':m,'degree':N,'finite_3cycle_branches':N//2,'discriminant_squareclass':D})
 # These rational identities retain the all-m lift calculation exactly.
 m=s.symbols('m',positive=True);cc=m*(m-s.I)/(m*m-1);D1=(m*m+1)/2
 condition=(1-cc)**2+2*cc*(1-cc)/D1+cc**2/m**2
 require(s.factor(condition)==0,'symbolic all-m lift-root identity')
 return {'bounded_laurent_and_critical_controls':checked,
  'symbolic_all_m_lift_root_identity':True,
  'all_odd_m_statement_requires_written_polynomial_and_monodromy_proof':True}


def full_Jacobian(data):
 started=time.monotonic();variables,F=expanded_map(3)
 require([f.total_degree()for f in F]==data['A12']['F_component_degrees'],'actual polynomial component degrees')
 require([len(f.terms())for f in F]==data['A12']['F_component_term_counts'],'actual polynomial coefficient counts')
 D=[[f.diff(x)for x in variables]for f in F]
 J=D[0][0]*(D[1][1]*D[2][2]-D[1][2]*D[2][1])-D[0][1]*(D[1][0]*D[2][2]-D[1][2]*D[2][0])+D[0][2]*(D[1][0]*D[2][1]-D[1][1]*D[2][0])
 require(J==s.Poly(-1,*variables,domain=s.QQ_I),'direct actual three-variable polynomial determinant')
 return {'full_expanded_polynomial_components':True,'component_degrees':[33,159,3],
  'component_term_counts':[77,1484,2],'direct_determinant':-1,
  'elapsed_seconds':round(time.monotonic()-started,3)}


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--full-jacobian',action='store_true');parser.add_argument('--write-receipt',type=Path)
 args=parser.parse_args();data=json.loads(INPUT.read_text());source=HERE/data['primary_source']['path']
 require(sha(source)==data['primary_source']['sha256'],'exact retained v1 source pin')
 report={'schema_version':1,'python':platform.python_version(),'sympy':s.__version__,
  'checker_sha256':sha(__file__),'model_sha256':sha(HERE/'alternating_shear_keller_model.py'),
  'input_sha256':sha(INPUT),'primary_source_sha256':sha(source),
  'literal_A12_certificate':literal_A12(data),'uniform_controls':uniform_controls()}
 bad=deepcopy(data);bad['A12']['B_terms'][2][0]=[-1,2]
 try:literal_A12(bad)
 except AssertionError:pass
 else:raise AssertionError('corrupt primary B coefficient accepted')
 try:binomial_terms(4)
 except ValueError:pass
 else:raise AssertionError('even shear degree accepted')
 report['negative_controls_rejected']=2
 if args.full_jacobian:report['direct_three_variable_Jacobian']=full_Jacobian(data)
 report['scope']={'formal_verification':False,'independent_replay':False,'external_review':False,
  'no_Q_coefficient_polynomial_lift_claim':True,'no_plane_realization_claim':True}
 report['result']='PASS'
 if args.write_receipt:args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: exact primary A12 coefficients/master/lifts; bounded uniform Laurent controls; symbolic all-m lift identity.')
 if args.full_jacobian:print('PASS: direct expanded three-variable polynomial determinant=-1.')
 print('checker_sha256='+report['checker_sha256']);print('model_sha256='+report['model_sha256']);print('input_sha256='+report['input_sha256'])

if __name__=='__main__':main()