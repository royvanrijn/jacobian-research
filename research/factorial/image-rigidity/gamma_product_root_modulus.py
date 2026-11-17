#!/usr/bin/env python3
"""Exact controls for the written all-gap Gamma-product root modulus theorem."""
import hashlib,itertools,json,platform
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
a,N,j,q,z,g=S.symbols("alpha N j q z gamma")
checks=[]
def check(name,ok):
 if not bool(ok):raise AssertionError((name,ok))
 checks.append(name)
def same(name,left,right):check(name,S.simplify(S.cancel(S.expand(left-right)))==0)
co=(1-q*q)/(1+q*q);si=2*q/(1+q*q)
same("generic-coherent-full-row",
 (N+1-2*j)*co+si*((j-1)/q+(N-j)*q),N-1)
same("generic-coherent-last-deficit",a*si*q,a*(1-co))
same("generic-alternating-diagonal-sine-flip",(-1)*si,-si)
# Complete characteristic polynomials, including every coefficient.
char_coeffs=0
for nn in range(1,6):
 J=S.diag(*[a+nn+1-2*j0 for j0 in range(1,nn+1)])
 for j0 in range(1,nn):
  J[j0-1,j0]=J[j0,j0-1]=S.I*S.sqrt(j0*(a+nn-j0))
 actual=S.Poly(J.charpoly(z).as_expr().subs(z,-z)*(-1)**nn,z)
 expected=S.Poly(sum(S.binomial(nn,k)*S.rf(a,k)*z**(nn-k) for k in range(nn+1)),z)
 for k in range(nn+1):
  same(f"Gamma-matrix-n{nn}-coefficient{k}",actual.nth(k),expected.nth(k));char_coeffs+=1
# Direct actual positive-vector matrices and the displayed sum-of-squares decomposition.
frame_entries=0
for nn,aa,qq in [(1,S.Rational(1,2),S.Rational(1,2)),
 (2,S.Rational(1,2),2),(3,S.Rational(7,3),S.Rational(1,2)),(4,3,1)]:
 aa,qq=S.sympify(aa),S.sympify(qq)
 NN=aa+nn;cc=(1-qq*qq)/(1+qq*qq);ss=2*qq/(1+qq*qq)
 D=S.diag(*[NN+1-2*k for k in range(1,nn+1)])
 K=S.zeros(nn)
 for k in range(1,nn):K[k-1,k]=K[k,k-1]=S.sqrt(k*(NN-k))
 H=cc*D+ss*K
 vv=S.Matrix([S.sqrt(S.binomial(NN-1,k-1))*qq**(k-1) for k in range(1,nn+1)])
 expected=S.Matrix([(NN-1)*vv[k] for k in range(nn)])
 expected[nn-1]-=aa*(1-cc)*vv[nn-1]
 for k in range(nn):
  same(f"coherent-vector-n{nn}-row{k}",(H*vv)[k],expected[k]);frame_entries+=1
 B=S.diag(*[S.simplify(NN-1-(H*vv)[k]/vv[k]) for k in range(nn)])
 for k in range(nn-1):
  edge=S.zeros(nn,1);edge[k]=1/vv[k];edge[k+1]=-1/vv[k+1]
  B+=H[k,k+1]*vv[k]*vv[k+1]*(edge*edge.T)
 for k in range(nn):
  for ell in range(nn):
   same(f"positive-form-n{nn}-entry{k}-{ell}",B[k,ell],((NN-1)*S.eye(nn)-H)[k,ell]);frame_entries+=1
# Full polarization composition identity with algebraically independent roots.
polar_coeffs=0
for nn in range(1,5):
 ds=S.symbols("d0:"+str(nn));xs=S.symbols("x0:"+str(nn))
 G=sum(S.binomial(nn,k)*S.rf(a,k)*z**(nn-k) for k in range(nn+1))
 F=0
 for degree in range(nn+1):
  coeff=S.expand(G).coeff(z,degree)/S.binomial(nn,degree)
  F+=coeff*sum(S.prod(xs[k] for k in subset) for subset in itertools.combinations(range(nn),degree))
 left=S.cancel(S.prod(ds)*F.subs({xs[k]:z/ds[k] for k in range(nn)}))
 previous=S.Poly(S.prod(z+d for d in ds),z)
 right=sum(previous.nth(nn-k)*S.rf(a,k)*z**(nn-k) for k in range(nn+1))
 for k in range(nn+1):
  same(f"polarization-n{nn}-coefficient{k}",S.expand(left).coeff(z,k),S.expand(right).coeff(z,k));polar_coeffs+=1
 check(f"polarization-exact-degree-n{nn}",S.degree(G,z)==nn)
# Residue-class rising-factorial decomposition, reciprocal normalization and endpoints.
decomposition=0
for dd in range(1,6):
 for kk in range(6):
  same(f"residue-class-d{dd}-k{kk}",S.rf(g,dd*kk),
       dd**(dd*kk)*S.prod(S.rf((g+j0)/dd,kk) for j0 in range(dd)));decomposition+=1
 for nn in (1,2,3):
  C=sum(S.binomial(nn,k)*S.prod(S.rf((g+j0)/dd,k) for j0 in range(dd))*z**(nn-k) for k in range(nn+1))
  t=S.symbols("t")
  T=sum(S.binomial(nn,k)*S.rf(g,dd*k)*t**k for k in range(nn+1))
  same(f"reciprocal-normalization-d{dd}-n{nn}",(dd**dd*t)**nn*C.subs(z,1/(dd**dd*t)),T)
  same(f"upper-product-d{dd}-n{nn}",dd**dd*S.prod((g+j0)/dd+nn-1 for j0 in range(dd)),S.rf(g+dd*(nn-1),dd))
  same(f"lower-product-d{dd}-n{nn}",dd**dd*S.prod((g+j0)/dd-nn+1 for j0 in range(dd)),S.rf(g-dd*(nn-1),dd))
 same(f"one-index-equality-d{dd}",-1/(-1/S.rf(g,dd)),S.rf(g,dd))
# Strict quadratic upper comparison and deliberate orientation errors.
M=S.symbols("M")
same("quadratic-final-margin",M*M-(M-2)*(M-1),3*M-2)
wrong=(a+1)*(a+2)
right=a*(a+1)
check("wrong-Gamma-index-rejected",S.expand(wrong-right)!=0)
check("reciprocal-orientation-rejected",
 S.expand((z*z+2*a*z+a*(a+1))-(1+2*a*z+a*(a+1)*z*z))!=0)
# Generic shape-shift ordering and the linear finite-prefix specialization.
rr,ss,dd,Delta,aa,beta,NN=S.symbols("r s d Delta a beta N")
same("generic-shift-radius-argument-margin",
 (g+Delta-dd*(ss-1))-(g+dd*(rr-1)),Delta-dd*(rr+ss-2))
same("generic-shift-lower-applicability-at-threshold",
 g+dd*(rr+ss-2)-dd*(ss-1),g+dd*(rr-1))
same("physical-shape-shift-margin",
 (aa*ss+beta-dd*(ss-1))-(aa*rr+beta+dd*(rr-1)),
 aa*(ss-rr)-dd*(rr+ss-2))
same("finite-prefix-threshold-identity",
 (2*NN-3)*(ss-rr)-(rr+ss-2),
 2*((NN-2)*(ss-rr)-(rr-1)))
# The maximum ratio over a complete small prefix occurs at adjacent top orders.
for bound in range(2,8):
 for r0 in range(1,bound):
  for s0 in range(r0+1,bound+1):
   check(f"prefix-ratio-N{bound}-r{r0}-s{s0}",
         S.Rational(r0+s0-2,s0-r0)<=2*bound-3)
# Full exact rational-shape pairs, at equality and above, not an index theorem.
literal_pairs=0
for d0,r0,s0,gg,excess in [
 (1,1,2,S.Rational(1,2),0),(2,2,3,S.Rational(2,3),0),
 (3,1,3,S.Rational(3,2),1),(4,2,4,S.Rational(7,3),0)]:
 delta=d0*(r0+s0-2)+excess
 t=S.symbols("t")
 lower=sum(S.binomial(r0,k)*S.rf(gg,d0*k)*t**k for k in range(r0+1))
 upper=sum(S.binomial(s0,k)*S.rf(gg+delta,d0*k)*t**k for k in range(s0+1))
 check(f"literal-shape-separated-d{d0}-r{r0}-s{s0}",S.degree(S.gcd(lower,upper),t)==0);literal_pairs+=1
# Full phase normalizations at nonordinary rational shapes.
for d0,a0,bb,j0 in [(1,3,S.Rational(1,2),2),(2,6,S.Rational(3,2),3),
 (3,9,S.Rational(2,3),2)]:
 for k in range(j0+1):
  same(f"phase-normalization-a{a0}-d{d0}-j{j0}-k{k}",
       S.rf(bb,a0*j0+d0*k),S.rf(bb,a0*j0)*S.rf(bb+a0*j0,d0*k))
check("wrong-shift-threshold-rejected",
      S.expand((g+Delta-dd*(ss-1))-(g+dd*(rr-1))-(Delta-dd*(rr+ss-3)))!=0)
names=["GAMMA_PRODUCT_ROOT_MODULUS.tex","gamma_product_root_modulus.py",
 "CUBIC_ADJOINT_STABILITY_SFC.tex","sources/wagner-multivariate-stable-polynomials-2011.pdf"]
result={"kind":"internal exact replay; all-gap quantifiers are proved in the written comparison and circular-polarization argument",
 "passed":True,"controls":len(checks),"control_names":checks,
 "scope":{"all_positive_real_Gamma_shapes":True,"every_gap_d_ge1":True,"all_indices_n_ge1":True,
  "upper_reciprocal_root_bound":True,"lower_bound_when_gamma_gt_d_times_n_minus_one":True,
  "strict_bounds_when_n_ge2":True,"shape_separated_pair_coprimality":True,"linear_pairwise_prefix_criterion":True,
  "all_index_SFC_for_fixed_lower_exponent":False,"Strong_Factorial":False,
  "unrestricted_FC_or_SFC":False},
 "control_categories":{"Gamma_characteristic_coefficients":char_coeffs,
 "actual_coherent_vectors_and_positive_forms":frame_entries,"polarization_coefficients":polar_coeffs,
 "residue_class_factorial_identities":decomposition,"literal_shape_separated_pairs":literal_pairs},
 "source_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 "exact_imports":{"Gamma_matrix":"CUBIC_ADJOINT_STABILITY_SFC.tex lem:gamma",
 "circular_GWS":"Wagner Theorem4.1 PDFpage11; degree=n allows exterior circular domains"},
 "python":platform.python_version(),"sympy":S.__version__,
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/gamma_product_root_modulus.py"}
(HERE/"gamma_product_root_modulus_verification.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"passed":True,"controls":len(checks),"control_categories":result["control_categories"]},indent=2))
