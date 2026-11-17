"""Replay algebra supporting the universal quadratic-pair classification.

The written interlacing induction is the all-weight argument. Finite
Sturm controls here are regression evidence, not its quantifier.
"""
from pathlib import Path
import hashlib,json,platform,time
import sympy as S
from quadratic_scalar_obstruction import curved_sector
from quadratic_continuant import leading_matrix,bottom_continuants
from rational_wave import rational_chart
from polynomial_wave import polynomial_chart
from pole_direction import wave_tensor
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
k=S.symbols("kappa")
i,j=S.symbols("i j",integer=True,nonnegative=True)
D=(i+4)*(i+3)+2*k*(i+3)*(j+1)+k*k*(j+1)*(j+2)
edgeD=D.subs({i:i-2,j:j+1},simultaneous=True)
check("generic top-factor negative discriminant",
 S.expand(S.discriminant(D,k)+4*(i+3)*(j+1)*(i+j+5))==0)
check("generic edge-factor negative discriminant",
 S.expand(S.discriminant(edgeD,k)+4*(i+1)*(j+2)*(i+j+4))==0)
a=2*(i+1)*(j+2);anext=2*(i-1)*(j+3)
upper=(j+2)*(j+3);lower=i*(i-1)
q=S.cancel(upper*lower/(a*anext))
check("generic normalized leading edge identity",
 S.cancel(q-i/(4*(i+1)))==0)
check("generic strict quarter gap",
 S.cancel(S.Rational(1,4)-q-1/(4*(i+1)))==0)
check("recurrence needs leading positivity",
 S.expand(k*k-(2*k*k+1))==-k*k-1 and
 S.Poly(-k*k-1,k).count_roots(-S.oo,S.oo)==0)
controls=[]
for m in [0,1,2,3,4,5,6,7,8,12,20,21,24,31,40]:
    T=curved_sector(m,k);B=T[1:,:];L=leading_matrix(m)
    H=bottom_continuants(m,k);n=B.rows
    check(f"complete bottom tridiagonal structure m={m}",
      all(B[r,c]==0 for r in range(n) for c in range(n) if abs(r-c)>1))
    piv=S.Integer(1);leading_product=S.Integer(1)
    principal_coefficients=[]
    for ell in range(1,n+1):
        row=ell-1;ii=m-2*row
        if row:
            edge_i=ii+2
            edgeq=S.Rational(edge_i,4*(edge_i+1))
            piv=1-edgeq/piv
        check(f"exact leading LDL m={m},ell={ell}",piv>S.Rational(1,2))
        leading_product*=L[row,row]*piv
        check(f"exact leading coefficient m={m},ell={ell}",
          H[ell].degree()==ell and H[ell].LC()==leading_product and leading_product>0)
        check(f"exact real-root count m={m},ell={ell}",
          H[ell].count_roots(-S.oo,S.oo)==ell and S.gcd(H[ell],H[ell].diff()).degree()==0)
        if ell>=2:
            check(f"consecutive continuants coprime m={m},ell={ell}",
              S.gcd(H[ell],H[ell-1]).degree()==0)
        principal_coefficients.append(str(leading_product))
    if m<=8 or m==12:
        for ell in range(1,n+1):
            check(f"independent actual minor m={m},ell={ell}",
              S.Poly(B[:ell,:ell].det(method="domain-ge"),k,domain=S.QQ)==H[ell])
    check(f"top and bottom coprime m={m}",
      S.gcd(S.Poly(S.prod(T[r,r] for r in range(n)),k),H[-1]).degree()==0)
    controls.append(dict(weight=m,dimension=n,bottom_degree=H[-1].degree(),
        bottom_real_root_count=n,leading_principal_coefficients=principal_coefficients))
# Three distinct unbounded witness families, reconstructed from older exact solvers.
w,z=S.symbols("w z")
witnesses=[]
for rv,sv,RR,power in [(4,0,1/w,4),(6,1,1/w,4),(8,3,1/w,4),
                      (2,0,w*w,2),(4,1,w*w,2),(3,4,w*w,2)]:
    C=rational_chart(rv,sv,RR,w,z,k)
    scalar=S.cancel(C["lambda_polynomial"]/w**power)
    S.Poly(scalar,w,z)
    b,c=(w*w,k*w*z+1) if power==4 else (w,w*w+k*z)
    check(f"reconstructed quadratic family wave r={rv},s={sv},den={power}",
      C["red"](wave_tensor(scalar,b,c,w,z))==0)
    expected=(sv+2)*(rv-sv-1) if power==4 else (sv+2)*(rv+2*sv+2)
    check(f"reconstructed quadratic family degree r={rv},s={sv},den={power}",C["N"]==expected)
    check(f"reconstructed quadratic family group r={rv},s={sv},den={power}",
      C["geometric_group"]==("S" if rv%2==0 and sv%2 else "A")+str(expected))
    witnesses.append(dict(r=rv,s=sv,denominator_degree=power//2,
      scalar_degree=int(S.Poly(scalar,w,z).total_degree()),degree=C["N"],
      geometric_group=C["geometric_group"]))
for sv,ep in [(0,0),(0,1),(1,0),(2,1),(5,0)]:
    pairs=[(0,0)]*(sv+1);pairs[sv]=(1,0) if ep else (0,1)
    C=polynomial_chart(w*w,pairs,w=w,z=z)
    check(f"quadratic heat wave s={sv},epsilon={ep}",
      wave_tensor(C["lambda_polynomial"],1,w*w,w,z)==0)
    expected=(sv+1)*(3*sv+5) if not ep else (sv+2)*(3*sv+4)
    check(f"quadratic heat degree s={sv},epsilon={ep}",C["N"]==expected)
    group="SL2(F3) on8 letters" if sv==0 and ep==1 else "A"+str(expected)
    check(f"quadratic heat group s={sv},epsilon={ep}",C["geometric_group"]==group)
    witnesses.append(dict(s=sv,epsilon=ep,degree=C["N"],geometric_group=group))
model=dict(universal_claims=["KMFI54","KMFI55"],
 proof_mechanism="positive leading principal coefficients plus elementary strict interlacing of bottom continuants; top minor has no real roots",
 all_weights_injective=True,all_bottom_continuants_simple_real_rooted=True,
 every_quadratic_primitive_wave_pair_reduces_to_completed_chart=True,
 scalar_degree_cutoff=None,
 field_groups="natural A_N, natural S_N, or binary SL2(F3) on8",
 unbounded_action_families=[
   dict(pair="(1,h),deg h=2",parameters="s>=0,epsilon in{0,1}",degree="(s+1)(3s+5) or (s+2)(3s+4)",group="A_N except even-centered affine-lambda binary8"),
   dict(pair="(w,w^2+kappa*z)",parameters="r>=2,s>=0,C(r,s)=0",degree="(s+2)(r+2s+2)",group="S_N iff r even,s odd; otherwise A_N"),
   dict(pair="(w^2,kappa*w*z+1)",parameters="r>=s+4,s>=0,C(r,s)=0",degree="(s+2)(r-s-1)",group="S_N iff r even,s odd; otherwise A_N")],
 only_prime_degree_in_quadratic_class=5,
 finite_regression_controls=controls,field_witnesses=witnesses,
 arbitrary_Keller_PSL2_11_actions_decided=False,
 absolute_natural_A4_decided=False,absolute_binary_degree8_decided=False,
 higher_primitive_pair_degrees_classified=False,
 full_action_classification_complete=False)
(ROOT/"quadratic_continuant_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["QUADRATIC_CLASSIFICATION.md","quadratic_continuant.py",
 "verify_quadratic_continuant.py","quadratic_continuant_model.json",
 "quadratic_scalar_obstruction.py","rational_wave.py","polynomial_wave.py","pole_direction.py"]
rec=dict(result="pass",check_count=len(checks),checks=checks,
 python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),
 sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Written all-weight recurrence/positivity/interlacing proof with generic symbolic algebra and exact finite regression controls; no extrapolation from Sturm counts, formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"quadratic_continuant_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"],
 controls=len(controls),all_weight_theorem=True,goal_complete=False)))
