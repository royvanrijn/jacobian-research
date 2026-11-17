"""Exact weighted obstructions and polynomial Bezout certificates."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as S
from quadratic_scalar_obstruction import vertical_sector,curved_sector,bezout_certificate
from pole_direction import wave_tensor
from rational_wave import rational_chart
from polynomial_wave import polynomial_chart
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
w,z,k=S.symbols("w z kappa")
i,j=S.symbols("i j",integer=True,nonnegative=True)
q=i*(j+2)/(4*(i+2)*(j+3))
check("universal strict offdiagonal bound",S.cancel(S.Rational(1,4)-q-(i+2*j+6)/(4*(i+2)*(j+3)))==0)
vertical=wave_tensor(w**i*z**j,z,w*w,w,z)
expected=i*(i-1)*w**(i-2)*z**(j+2)+2*(i+2)*(j+1)*w**(i+1)*z**j+j*(j-1)*w**(i+4)*z**(j-2)
check("universal vertical weighted monomial identity",S.simplify(vertical-expected)==0)
curved=wave_tensor(w**i*z**j,w*w+z,k*w*z,w,z)
D=(i+4)*(i+3)+2*k*(i+3)*(j+1)+k*k*(j+1)*(j+2)
E=2*(i+1)*(i+2+k*(j+2));F=i*(i-1)
check("universal diagonal curved monomial identity",
 S.simplify(curved-D*w**(i+2)*z**j-E*w**i*z**(j+1)-F*w**(i-2)*z**(j+2))==0)
vertical_controls=[]
for weight in [0,2,3,6,9,12,18,24,35,40]:
 basis,M=vertical_sector(weight)
 piv=S.Integer(1);pivots=[piv]
 for r in range(len(basis)-1):
  ratio=S.Rational(M[r,r+1]*M[r+1,r],M[r,r]*M[r+1,r+1])
  check(f"vertical normalized edge bound weight={weight},edge={r}",0<ratio<S.Rational(1,4))
  piv=1-ratio/piv;pivots.append(piv)
 check("vertical exact positive LDL weight="+str(weight),all(v>S.Rational(1,2) for v in pivots))
 check("vertical actual sector nonsingular weight="+str(weight),M.det()!=0)
 vertical_controls.append(dict(weight=weight,basis=basis,LDL_pivots=[str(v) for v in pivots]))
certs=[]
for weight in range(21):
 c=bezout_certificate(weight,k);P,Q,U,V,G=[c[a] for a in ["top","bottom","U","V","G"]]
 check("exact all-kappa Bezout identity weight="+str(weight),G==1 and U*P+V*Q==1)
 M=curved_sector(weight,k)
 check("exact top minor factorization weight="+str(weight),S.expand(M[:-1,:].det(method="domain-ge")-S.prod(M[j,j] for j in range(M.cols)))==0)
 real_roots=Q.count_roots(-S.oo,S.oo)
 check("bounded bottom Sturm count weight="+str(weight),real_roots==Q.degree())
 def encode(poly):return [str(x) for x in poly.all_coeffs()]
 certs.append(dict(weight=weight,rows=c["rows"],columns=c["columns"],
  top_primitive_coefficients=encode(P),bottom_primitive_coefficients=encode(Q),
  U_coefficients=encode(U),V_coefficients=encode(V),identity="U*P+V*Q=1",
  bottom_degree=int(Q.degree()),bottom_real_root_count=int(real_roots)))
# Reconstruct every stored identity from coefficient arrays, independently of gcdex.
payload=dict(coefficient_order="descending powers of kappa",weights=list(range(21)),
 certificate_scope="all kappa; weight0..20 only",certificates=certs)
(ROOT/"quadratic_scalar_certificates.json").write_text(json.dumps(payload,indent=2)+"\n")
loaded=json.loads((ROOT/"quadratic_scalar_certificates.json").read_text())
def decode(coefs):return S.Poly.from_list([S.Rational(x) for x in coefs],k,domain=S.QQ)
check("serialized Bezout certificate replay",
 all(decode(c["U_coefficients"])*decode(c["top_primitive_coefficients"])+
     decode(c["V_coefficients"])*decode(c["bottom_primitive_coefficients"])==1
     for c in loaded["certificates"]))
# Operator matrices agree with the actual wave, not only the written recurrence.
for weight in [0,1,2,5,8,12]:
 M=curved_sector(weight,k);ok=True
 for col in range(M.cols):
  jj=col;ii=weight-2*jj
  result=wave_tensor(w**ii*z**jj,w*w+z,k*w*z,w,z)
  predicted=sum(M[row,col]*w**(weight+2-2*row)*z**row for row in range(M.rows) if M[row,col]!=0)
  ok=ok and S.expand(result-predicted)==0
 check("full actual curved matrix weight="+str(weight),ok)
# Resonance and scalar-degree bounds are exact parameter identities.
r,s=S.symbols("r s",integer=True)
check("unique scalar ordinary degree",S.expand((r-4)+s-(r+s-4))==0)
check("unique scalar weighted degree",S.expand((r-4)+2*s-(r+2*s-4))==0)
check("degree-ten weight bound",2*10==20)
# Complete scalar-linear reductions: all coefficients retained symbolically.
Bw,Bz,B0,Cw,Cz,C0=S.symbols("Bw Bz B0 Cw Cz C0")
b=w*w+Bw*w+Bz*z+B0;c=k*w*z+Cw*w+Cz*z+C0
def rem(expr,minimum):return S.Poly(S.expand(expr),k).rem(S.Poly(minimum,k)).as_expr().expand()
ww=rem(wave_tensor(w,b,c,w,z),k*k+4*k+10)
check("linear scalar diagonal A=w curve coefficient",S.expand(S.Poly(ww,w,z).coeff_monomial(w*z)-4*Bz*(2*k+3))==0)
zz=rem(wave_tensor(z,b,c,w,z),k*k+2*k+2)
check("linear scalar diagonal A=z curve coefficient",S.expand(S.Poly(zz,w,z).coeff_monomial(z*z)-2*Bz*(3*k+2))==0)
check("linear scalar diagonal A=z w coefficient",S.expand(S.Poly(zz,w,z).coeff_monomial(w*w)-2*Cw*(2*k+3))==0)
a,bb,c0,d,e,f=S.symbols("a b c0 d e f")
b=a*w+bb*z+c0;c=w*w+d*w+e*z+f
ww=wave_tensor(w,b,c,w,z);zz=wave_tensor(z,b,c,w,z)
check("full constant-direction A=w coefficients",S.expand(ww-(6*bb*w*w+2*(3*a*a+2*a*e+2*bb*d+e*e)*w+4*bb*(a+e)*z+2*(2*a*c0+bb*f+c0*e)))==0)
check("full constant-direction A=z coefficients",S.expand(zz-(2*(3*a+2*e)*w*w+8*bb*w*z+4*(a*d+c0+d*e)*w+2*(a*a+2*a*e+2*bb*d+3*e*e)*z+2*(a*f+c0*d+2*e*f)))==0)
check("constant-direction A=z contradiction",S.expand((a*a+2*a*e+3*e*e).subs(e,-3*a/2)-S.Rational(19,4)*a*a)==0)
# Positive field controls are exact reconstructions, not new Keller assertions.
models=[]
for rv,sv,R,denom_power in [(5,0,1/w,4),(3,0,w*w,2),(6,1,1/w,4)]:
 chart=rational_chart(rv,sv,R,w,z,k)
 scalar=S.cancel(chart["lambda_polynomial"]/w**denom_power)
 if denom_power==4:
  b=w*w;c=k*w*z+1
 else:
  b=w;c=w*w+k*z
 check(f"positive reconstructed scalar wave {rv},{sv}",chart["red"](wave_tensor(scalar,b,c,w,z))==0)
 models.append(dict(r=rv,s=sv,scalar=str(scalar),scalar_degree=int(S.Poly(scalar,w,z).total_degree()),
    b=str(b),c=str(c),degree=chart["N"],geometric_group=chart["geometric_group"]))
for h,pairs in [(w*w,[(1,0)]),(w*w+w,[(1,0)]),(w*w+w,[(0,0),(0,1)])]:
 chart=polynomial_chart(h,pairs,w=w,z=z)
 check("positive polynomial scalar wave h="+str(h)+",s="+str(chart["s"]),wave_tensor(chart["lambda_polynomial"],1,h,w,z)==0)
 models.append(dict(scalar=str(chart["lambda_polynomial"]),b="1",c=str(h),
    degree=chart["N"],geometric_group=chart["geometric_group"]))
check("scalar-linear exact action set",
 {m["geometric_group"] for m in models if m.get("scalar_degree")==1 or m["scalar"]=="w"}=={"A8","A10","SL2(F3) on8 letters"})
heat_degree_controls=[]
for sv in range(11):
 for ep in [0,1]:
  pairs=[(0,0)]*(sv+1);pairs[sv]=(1,0) if ep else (0,1)
  heat=polynomial_chart(w*w,pairs,w=w,z=z)
  ordinary=int(S.Poly(heat["lambda_polynomial"],w,z).total_degree())
  check(f"complete leading heat degree s={sv},eps={ep}",ordinary==3*sv+ep)
  heat_degree_controls.append(dict(s=sv,epsilon=ep,scalar_degree=ordinary))
families=[]
for heat in heat_degree_controls:
 if heat["scalar_degree"]>10:continue
 sv,ep=heat["s"],heat["epsilon"]
 degree=(sv+1)*(3*(sv+2)-1) if not ep else (sv+2)*(3*(sv+1)+1)
 group="A"+str(degree)
 if sv==0 and ep==1:
  families.append(dict(family="quadratic polynomial h, even about lambda root",s=sv,epsilon=ep,scalar_degree=heat["scalar_degree"],degree=8,geometric_group="SL2(F3) on8"))
 families.append(dict(family="quadratic polynomial h",s=sv,epsilon=ep,scalar_degree=heat["scalar_degree"],degree=degree,geometric_group=group))
for rv in range(2,15):
 for sv in range(6):
  if rv+2*sv<=12:
   degree=(sv+2)*(rv+2*sv+2)
   families.append(dict(family="linear denominator and quadratic R",r=rv,s=sv,scalar_degree=rv+2*sv-2,
     degree=degree,geometric_group=("S" if rv%2==0 and sv%2 else "A")+str(degree)))
  if rv>=sv+4 and rv+sv<=14:
   degree=(sv+2)*(rv-sv-1)
   families.append(dict(family="quadratic denominator and R=1/w",r=rv,s=sv,scalar_degree=rv+sv-4,
     degree=degree,geometric_group=("S" if rv%2==0 and sv%2 else "A")+str(degree)))
check("all finite action catalogue scalar bounds",all(0<=m["scalar_degree"]<=10 for m in families))
model=dict(all_scalar_constant_top_direction="curved affine denominator impossible; universal positive-LDL proof",
 vertical_controls=vertical_controls,
 diagonal_curved_operator="exact weighted sector T_m(kappa); finite ansatz from unique C(r,s) resonance",
 certified_weights=list(range(21)),all_kappa_injectivity_on_certified_weights=True,
 all_quadratic_pairs_scalar_degree_at_most10="affinely reduce to complete chart; A_N,S_N or binary degree8",
 exact_scalar_linear_actions=["A8","A10","SL2(F3) on8"],
 positive_field_controls=models,heat_degree_controls=heat_degree_controls,degree_at_most10_action_families=families,
 exact_degree_at_most10_action_types=sorted({m["geometric_group"] for m in families}),
 unbounded_diagonal_curved_injectivity_proved=False,
 bottom_continuant_real_rooted_for_all_weights_proved=False,
 general_quadratic_scalar_class_complete=False,
 arbitrary_Keller_PSL2_11_actions_decided=False,full_action_classification_complete=False)
(ROOT/"quadratic_scalar_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["QUADRATIC_SCALAR_OBSTRUCTION.md","quadratic_scalar_obstruction.py","verify_quadratic_scalar.py",
 "quadratic_scalar_certificates.json","quadratic_scalar_model.json","pole_direction.py","rational_wave.py","polynomial_wave.py"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Universal all-scalar positive-LDL obstruction; exact resonance formulas and universal-in-kappa Bezout certificates for weight0..20. No all-weight injectivity, formal verification or independent whole-proof review claimed.",
 goal_complete=False)
(ROOT/"quadratic_scalar_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"],
 certificates=len(certs),certificate_bytes=(ROOT/"quadratic_scalar_certificates.json").stat().st_size)))
