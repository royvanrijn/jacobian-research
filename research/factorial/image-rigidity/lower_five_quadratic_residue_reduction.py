"""Exact internal replay of residue exclusions and the global gcd bound.

Written proofs cover14/16 index classes and all-index gcd degree0or2.
The degree bound also uses the retained squarefreeness theorem.
"""
from pathlib import Path
from hashlib import sha256
from math import comb
import json
import sympy as s
u,M,n=s.symbols("u M n")
checks=[]
def check(name,value):
 if not bool(value):raise AssertionError(name)
 checks.append(name)
def vp(a):
 a=int(a)
 if not a:return 10**9
 return (abs(a)&-abs(a)).bit_length()-1
def euler(v):
 r,e,k=v;g=M-2*n
 return [-2*u*s.diff(r,u)-g*u*k,
         r-e-2*u*s.diff(e,u)+u*k,
         e+M*k-2*u*s.diff(k,u)]
jet=[[s.Integer(1),s.Integer(0),s.Integer(0)]]
for j in range(7):
 ev=euler(jet[-1])
 jet.append([s.expand(ev[i]+j*jet[-1][i]) for i in range(3)])
ss=2*(n+1)
U=s.rf(M+2,5)+3*ss*(M+4)**2*u+ss*u*u
W=ss*u*((M+4)*(M+5)*(M+6)+(M+2*n+8)*u)
upper=[s.expand(u*jet[5][i]-jet[7][i]) for i in range(3)]
check("universal_shift_five_E_component",s.expand(upper[1]+W)==0)
check("universal_shift_five_K_component",s.expand(upper[2]+U)==0)
check("corrupted_shift_coefficient_rejected",s.expand(upper[2]+U+u)!=0)
check("ordinary_W_quadratic_factor",
      s.expand((M+2*n+8).subs(M,7*n+1)-9*(n+1))==0)
z,X=s.symbols("z X")
def bc(N):
 g=5*N+1;b=s.Integer(1)
 for j in range(N+1):
  c=s.expand((z+2*j+g+1)*b-z*s.diff(b,z))
  if j==N:return b,c
  b=c
def section(p,parity):
 return s.Poly(sum(c*u**((k-parity)//2)
                   for (k,),c in s.Poly(p,z).terms() if k%2==parity),u).as_expr()
for N in [1,2,3,4,5,6,7,8,9,10,15,17,19,21,23,31,33,35,37]:
 mm=7*N+1
 b,c=bc(N)
 rb=[s.Poly(section(b,i),u,modulus=2).eval(0) for i in range(2)]
 rc=[s.Poly(section(c,i),u,modulus=2).eval(0) for i in range(2)]
 check(f"complete_unit_cyclic_components_n{N}",
       rb==rc==([0,1] if N%2==0 else [1,1]))
 coeff=[];r=1;g=(5*N+1)//2
 for k in range(N+1):
  if k:r*=g+k-1
  if N%2:
   direct=comb(N,k)
   rr=s.rf(5*N+1,2*k)
   check(f"exact_rising_block_valuation_n{N}_k{k}",
         vp(direct*rr)==k+vp(direct)+vp(r))
   vv=vp(direct*rr)
   if N%4==1 and k>=2:
    check(f"strict_after_linear_bound_n{N}_k{k}",2*vv>3*k-1)
   if N%8==7 and k>=1:
    check(f"gamma_div_four_root_bound_n{N}_k{k}",vv>k)
   if N%16==3 and k>=1:
    check(f"gamma_div_eight_root_bound_n{N}_k{k}",2*vv>3*k)
 if N%2==0:
  aa=int(s.rf(mm+2,5));m=vp(aa);q=2**(m-1)
  uu=U.subs({M:mm,n:N,u:q*X})
  ww=W.subs({M:mm,n:N,u:q*X})
  F=s.Poly(s.expand((uu*section(c,1).subs(u,q*X)
                     +ww*section(b,1).subs(u,q*X))/2**m),X)
  check(f"full_even_scalar_integral_n{N}",F.domain==s.ZZ)
  check(f"full_even_scalar_affine_n{N}",
        s.Poly(F.as_expr()-1-X,X,modulus=2).is_zero)

# Complete modulo16 partition; no search of moment gcds.
def cls(N):
 if N%2==0:return "even"
 ell=vp(N+1)
 if ell>=3:return "ell>=3"
 if ell==1:
  b=vp(7*N+5)
  return "ell1_b2" if b==2 else "ell1_b3" if b==3 else "open13"
 a=vp(7*N+3)
 return "ell2_a3" if a==3 else "open11"
for residue in range(16):
 NN=residue if residue else 16
 check(f"all_mod16_classes_{residue}",
       (cls(NN).startswith("open"))==(residue in [11,13]))
# Actual full finite polynomials control the two quadratic component flags.
for N in [11,43]:
 b,c=bc(N);mm=7*N+1
 check(f"quadratic_profile_parameter_valuations_n{N}",
       vp(N+1)==2 and vp(mm+2)==4 and vp(s.rf(mm+2,5))==7)
 for i in [0,1]:
  H=s.expand((U.subs({M:mm,n:N,u:4*X})*section(c,i).subs(u,4*X)
              +W.subs({M:mm,n:N,u:4*X})*section(b,i).subs(u,4*X))/128)
  HH=s.Poly(H,X)
  check(f"full_quadratic_component_integrality_n{N}_i{i}",HH.domain==s.ZZ)
  check(f"full_quadratic_component_profile_n{N}_i{i}",
        s.Poly(H-1-X-X*X,X,modulus=2).is_zero)
check("quadratic_residual_irreducible",
      s.Poly(1+X+X*X,X,modulus=2).is_irreducible)
here=Path(__file__).resolve().parent
repo=here.parents[2]
prior=here/"LOWER_TWO_PRIME_GAP_SFC.tex"
simple=here/"CUBIC_ADJOINT_STABILITY_SFC.tex"
assert sha256(prior.read_bytes()).hexdigest()=="8453b5b0ccd220573ef106fa26790e1f691885970b8c2116f036ed893ab2c200"
assert sha256(simple.read_bytes()).hexdigest()=="ea2510e511545b9dd7a49a0796d75c45575687e60770c3359d38c5bc99cb4125"
check("quadratic_residual_derivative_one",
      s.Poly(s.diff(1+X+X*X,X)-1,X,modulus=2).is_zero)

# Written all-odd-index residual order2 theorem:72 supporting controls.
prior_control_count=len(checks)
assert prior_control_count==580
A0=s.Matrix([[1,0],[1,0]])
check("global_exact_reduced_idempotent_operator",A0*A0==A0)
check("global_positive_seed_component_unit",A0*s.Matrix([1,0])==s.Matrix([1,1]))
for NN in [1,3,5,7,9,11,13,15,19,27,29,43,45,59]:
 ell=vp(NN+1);mm=7*NN+1
 aa=vp(mm+2);bb=vp(mm+4)
 UU=s.Poly(U.subs({M:mm,n:NN}),u)
 WW=s.Poly(W.subs({M:mm,n:NN}),u)
 check(f"global_full_U_valuation_formula_n{NN}",
       [vp(UU.nth(i)) for i in range(3)]
       ==[ell+aa+bb,ell+1+2*bb,ell+1])
 check(f"global_full_W_valuation_formula_n{NN}",
       [vp(WW.nth(i)) for i in [1,2]]
       ==[2*ell+1+bb,2*ell+1])
 b,c=bc(NN)
 JJ=s.Poly(s.expand((UU.as_expr()*section(c,1)
                     +WW.as_expr()*section(b,1))/2**(ell+1)),u)
 check(f"global_whole_polynomial_integrality_n{NN}",JJ.domain==s.ZZ)
 check(f"global_whole_residual_identity_n{NN}",
       s.Poly(JJ.as_expr()-u*u*section(c,1),u,modulus=2).is_zero)
 check(f"global_exact_residual_order_two_n{NN}",
       JJ.nth(0)%2==JJ.nth(1)%2==0 and JJ.nth(2)%2==1)
global_control_count=len(checks)-prior_control_count
assert global_control_count==72

receipt={
 "passed":True,"checks":len(checks),"scope":
 "Written14/16 residue theorem, uniform11mod32 obstruction and global all-index gcd degree0or2; finite exact controls support written proofs",
 "remaining_lower_index_classes":"n11,13mod16",
 "remaining_terminal_classes":"N12,14mod16",
 "all_index_support_five_seven_sfc":False,
 "n11mod32_gcd_degree_zero_or_two":True,
 "global_all_index_gcd_degree_zero_or_two":True,
 "global_nonconstant_gcd":"irreducible quadratic of negative discriminant",
 "global_odd_component_residual_order":2,
 "global_root_count_proof":"written synthetic division preserving integrality; no Gauss theorem import",
 "control_sections":{"prior":prior_control_count,"global":global_control_count},
 "quadratic_bound_uses_written_squarefreeness_input":True,
 "command":"python research/factorial/image-rigidity/lower_five_quadratic_residue_reduction.py",
 "sympy_version":s.__version__,
 "no_moment_gcd_prefix_extrapolation":True,
 "source_sha256":{p.name:sha256(p.read_bytes()).hexdigest()
                  for p in [here/"LOWER_FIVE_QUADRATIC_RESIDUE_REDUCTION.tex",Path(__file__),prior,simple]},
 "independent_replay":False,"formal_verification":False,"external_review":False}
(here/"lower_five_quadratic_residue_reduction_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print("PASS",len(checks),"exact controls")
