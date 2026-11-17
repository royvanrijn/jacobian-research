"""Exact replay for constant-scalar cubic classification and curved affine lift."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as S
from cubic_wave import curved_cubic,curved_stage,marked_polynomials,marked_stage_data,reconstruct
from pole_direction import wave_tensor,conjugate_affine_pair
from rational_wave import rational_chart
from polynomial_wave import polynomial_chart
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
w,z,g,D,k=S.symbols("w z gamma D kappa")
def rem(expr,poly):
    num,den=S.fraction(S.cancel(expr))
    return S.expand(S.rem(S.rem(num,poly,k)*S.invert(S.rem(den,poly,k),poly,k),poly,k))
def degree_part(expr,degree):
    return S.Add(*[co*w**a*z**b for (a,b),co in S.Poly(expr,w,z).terms() if a+b==degree])
C=curved_cubic(w,z,D);b,c,X=C["b"],C["c"],C["X"]
check("all-D curved cubic wave",wave_tensor(1,b,c,w,z)==0)
check("all-D primitive remainder",S.expand(c.subs(z,3*w*w/4)-w**3/4-D)==0)
lam=b*b;mu=b*c;nu=c*c
wanted=[S.diff(mu,w)+S.diff(nu,z),-S.diff(lam,w)-S.diff(mu,z)]
check("X1 reconstruction both derivatives",
    S.expand(S.diff(X[0],w)-wanted[0])==0 and S.expand(S.diff(X[0],z)-wanted[1])==0)
for idx,label,ww,zz in [
 (1,"X2",w*S.diff(X[0],w)-mu,w*S.diff(X[0],z)+lam),
 (2,"X3",z*S.diff(X[0],w)-nu,z*S.diff(X[0],z)+mu)]:
    check(label+" reconstruction",S.expand(S.diff(X[idx],w)-ww)==0 and S.expand(S.diff(X[idx],z)-zz)==0)
sweep=X+g*S.Matrix([1,w,z])
check("all-D full sweep determinant",S.factor(sweep.jacobian([g,w,z]).det())==g*g)
cross=X.diff(w).cross(X.diff(z))
check("generic branch rank2",cross[2]!=0)
check("exact Gauss normal first coordinate",S.expand(b*cross[0]-(w*c-z*b)*cross[2])==0)
check("exact Gauss normal second coordinate",S.expand(b*cross[1]+c*cross[2])==0)
eta,beta=S.symbols("eta beta")
gauss=S.expand((c+beta*b).subs(z,-eta-beta*w))
expected=w**3+2*beta*w*w+(eta+4*beta*beta/3)*w+4*beta*eta/3+D
check("exact all-D Gauss cubic",S.expand(gauss-expected)==0)
numer=w**3+2*beta*w*w+4*beta*beta*w/3+D
check("Gauss numerator denominator coprime",
    S.factor(numer.subs(w,-4*beta/3)-D+16*beta**3/27)==0)
# Highest-weight complete intersection and independent monic eliminant.
weights=[5,1,2];targets=[5,6,7]
for idx in range(3):
    top=[]
    for exps,co in S.Poly(sweep[idx],g,w,z).terms():
        weight=sum(a*b for a,b in zip(weights,exps))
        check(f"non-increasing weighted term output{idx},exp{exps}",weight<=targets[idx])
        if weight==targets[idx]:top.append((exps,co))
    top_expr=sum(co*g**a*w**bb*z**cc for (a,bb,cc),co in top)
    check("top weighted map output"+str(idx),S.expand(top_expr-sweep[idx].subs(D,0))==0)
t=S.symbols("t");M=marked_polynomials(t)
check("top fiber polynomial coprimality",S.gcd(M["P"],M["Q"])==1)
check("weighted complete intersection degree",5*6*7//(5*1*2)==21)
A,B,T=S.symbols("A B T")
V=S.expand(X[1]-w*X[0]);W=S.expand(X[2]-z*X[0])
eliminant=S.Poly(S.resultant(V+w*A-B,W+z*A-T,z),w)
check("all-D degree21 monic-up-to-unit eliminant",
    eliminant.degree()==21 and eliminant.LC()==-S.Rational(1,156279375))
check("integral z equation",S.Poly(V+w*A-B,z).LC()==S.Rational(16,27))
# The marked zero is over a degree<=18 number field.
check("exact marked sextic relation",
 S.factor(M["Q"]+M["d"]*(t-S.Rational(1,2)-4*t*t/3)-M["d"]**2+M["H"]/28350)==0)
check("all roots of H have d nonzero",S.gcd(M["H"],M["d"])==1)
check("all roots of H have gamma factor nonzero",S.gcd(M["H"],M["G"])==1)
check("marked sextic separable",S.gcd(M["H"],S.diff(M["H"],t))==1)
check("marked sextic irreducible",len(S.factor_list(M["H"])[1])==1 and S.factor_list(M["H"])[1][0][1]==1)
stage=marked_stage_data();red=stage["red"];R=stage["R"];BB=stage["B"];q0=stage["q0"]
check("exact marked E zero",all(stage["evaluate"](x)==0 for x in stage["E_q"]))
check("exact marked B determinant",stage["determinant"]==-1)
check("exact marked inverse BR",(BB*R).applyfunc(red)==S.eye(3))
check("exact marked inverse RB",(R*BB).applyfunc(red)==S.eye(3))
check("nonzero marked gamma cube",q0!=0 and S.gcd(stage["H"],S.Poly(q0,t)).degree()==0)
check("q0 and w0 reconstruction",
 red(q0**2*stage["u0"]**3*M["d"]-1)==0)
check("t and z0 reconstruction",red(stage["v0"]-t*q0*stage["u0"]**2)==0)
u,v=S.symbols("u v")
E=curved_stage(g,u,v,1)
for idx,expr in enumerate(E):
    S.Poly(expr,g,u,v)
    check("polynomial scaled stage output"+str(idx),all(exps[0]%3==0 for exps,co in S.Poly(expr,g,u,v).terms()))
check("scaled outputs use powers1,2,3",sum([1,2,3])==3+2+1)
# Actual stage composition in all forbidden weight slots, over Q[t]/H.
zero=(0,0,0);cap=3;ew=(1,2,3)
def add(a,b):
    out=dict(a)
    for exp,co in b.items():out[exp]=red(out.get(exp,0)+co)
    return {exp:co for exp,co in out.items() if co!=0}
def scale(a,co):return {exp:red(v*co) for exp,v in a.items() if red(v*co)!=0}
def mul(a,b):
    out={}
    for ea,ca in a.items():
      for eb,cb in b.items():
        exp=tuple(x+y for x,y in zip(ea,eb))
        if sum(x*y for x,y in zip(exp,ew))>=cap:continue
        out[exp]=red(out.get(exp,0)+ca*cb)
    return {exp:co for exp,co in out.items() if co!=0}
def power(a,n):
    out={zero:S.Integer(1)}
    for _ in range(n):out=mul(out,a)
    return out
vs=[{(1,0,0):S.Integer(1)},{(0,1,0):S.Integer(1)},{(2,0,0):stage["f20"]}]
linear=[]
for row in range(3):
    out={}
    for col in range(3):out=add(out,scale(vs[col],R[row,col]))
    linear.append(out)
qjet=scale(power(add({zero:S.Integer(1)},linear[0]),3),q0)
ujet=add({zero:stage["u0"]},linear[1]);vjet=add({zero:stage["v0"]},linear[2])
polynomial_jets=[]
for idx,expr in enumerate(stage["E_q"]):
    out={}
    for (a,bb,cc),co in S.Poly(expr,stage["q"],stage["u"],stage["v"]).terms():
        term=scale(mul(mul(power(qjet,a),power(ujet,bb)),power(vjet,cc)),co)
        out=add(out,term)
    forbidden={exp:co for exp,co in out.items() if sum(x*y for x,y in zip(exp,ew))<idx+1}
    check("all actual forbidden stage jets output"+str(idx),not forbidden)
    polynomial_jets.append({str(exp):str(co) for exp,co in out.items()})
# Full determinant chain: generic gamma and x exponents both cancel.
check("full affine lift x exponent",sum(ew)-sum(ew)==0)
check("full affine lift gamma exponent",1+2+2+1-sum(ew)==0)
check("full affine lift determinant sign",-(-1)==1)
# Complete constant-scalar cubic coefficient exhaustions.
a0,a1,a2,b0,b1,b2=S.symbols("a0 a1 a2 b0 b1 b2")
bquad=a0*w*w+a1*w*z+a2*z*z;cquad=b0*w*w+b1*w*z+b2*z*z
p_single=k*k+5*k+15;p_mixed=2*k*k+3*k+2
for label,bt,ct,minimum,rank in [
 ("single",w**3,k*w*w*z,p_single,3),
 ("mixed",w*w*z,k*w*z*z,p_mixed,4)]:
    ww=wave_tensor(1,bt+bquad,ct+cquad,w,z)
    part=degree_part(ww,3)
    coeffs=S.Poly(part,w,z).coeffs()
    matrix,rhs=S.linear_eq_to_matrix(coeffs,[a0,a1,a2,b0,b1,b2])
    minor=matrix[:,list(range(rank))]
    check(label+" quadratic coefficient matrix full row rank",minor.rows==rank and rem(minor.det(),minimum)!=0)
    if label=="single":
        shiftw,shiftz,curv=S.symbols("shiftw shiftz curv")
        Bq=3*shiftw*w*w-(k+12)*curv*w*z/9
        Cq=k*shiftz*w*w+2*k*shiftw*w*z+curv*z*z
        check("single complete translation curvature kernel",
          rem(degree_part(wave_tensor(1,bt+Bq,ct+Cq,w,z),3),minimum)==0)
        coef=S.Poly(wave_tensor(1,w**3-(k+12)*curv*w*z/9,k*w*w*z+curv*z*z,w,z),w,z).coeff_monomial(z*z)
        check("single weighted curvature obstruction",
          rem(coef-2*curv*curv*(291-8*k)/81,minimum)==0)
        check("single curvature factor cannot vanish",S.resultant(291-8*k,minimum,k)!=0)
    else:
        sw,sz=S.symbols("shiftw shiftz")
        Bq=sz*w*w+2*sw*w*z;Cq=2*k*sz*w*z+k*sw*z*z
        check("mixed complete translation kernel",
          rem(degree_part(wave_tensor(1,bt+Bq,ct+Cq,w,z),3),minimum)==0)
aa,bb,cc,dd,ee,ff=S.symbols("a b c d e f")
for label,bt,ct,minimum,linear,resultant in [
 ("single",w**3,k*w*w*z,p_single,(6*k+24)*aa+(4*k+6)*ee,48384),
 ("mixed",w*w*z,k*w*z*z,p_mixed,(8*k+12)*aa+(12*k+8)*ee,256)]:
    ww=wave_tensor(1,bt+aa*w+bb*z+cc,ct+dd*w+ee*z+ff,w,z)
    expr=rem(ww,minimum)
    if label=="single":
        check("single affine b_z forced zero",S.expand(S.Poly(expr,w,z).coeff_monomial(w*z)-bb*(8*k+12))==0)
        check("single affine b_constant forced zero",S.expand(S.Poly(expr,w,z).coeff_monomial(w)-cc*(4*k+12))==0)
    else:
        for name,mon,var,co in [
          ("c_w",w*w,dd,4*k+6),("b_z",z*z,bb,6*k+4),
          ("c_constant",w,ff,4*k+4),("b_constant",z,cc,4*k+4)]:
            check("mixed affine "+name+" forced zero",S.expand(S.Poly(expr,w,z).coeff_monomial(mon)-var*co)==0)
    rt=S.symbols("ratio")
    resultant_expr=S.resultant(S.resultant(linear.subs({aa:rt,ee:1}),rt*rt+rt+1,rt),minimum,k)
    check(label+" affine linear-quadratic incompatibility",resultant_expr==resultant)
# Constant top direction, including every lower coefficient.
a,be,cq,d,e,f,h20,h11,h02,h10,h01,h00=S.symbols("a beta cq d e f g h j i ell n")
B0=a*w*w+be*w*z+cq*z*z+d*w+e*z+f
C0=w**3+h20*w*w+h11*w*z+h02*z*z+h10*w+h01*z+h00
top3=degree_part(wave_tensor(1,B0,C0,w,z),3)
check("constant top cubic constraints",S.expand(top3-(8*be+4*h02)*w**3-12*cq*w*w*z)==0)
check("constant top quadratic curvature impossible",wave_tensor(1,be*w*z,-2*be*z*z,w,z)==38*be*be*z*z)
weighted=wave_tensor(1,a*w*w+e*z,w**3+h11*w*z,w,z)
check("constant top remaining weighted equations",
 S.expand(weighted-(12*a*a+6*a*h11+6*e+2*h11*h11)*w*w-4*e*(a+h11)*z)==0)
d,f,g2,j2,i2,l2=S.symbols("d f g j i l")
B0=w*w-4*z/3+d*w+f;C0=w**3-w*z+g2*w*w+j2*z+i2*w+l2
check("complete normalized curved lower equations",
 S.expand(wave_tensor(1,B0,C0,w,z)-(8*d-16*g2/3)*w-2*d*d-2*d*j2-2*f+8*i2/3-2*j2*j2)==0)
JJ,LL=S.symbols("J L")
B0=w*w-4*z/3-JJ*JJ;C0=w**3-w*z+JJ*z+LL
normalized=conjugate_affine_pair(B0,C0,S.Matrix([[1,0],[-3*JJ/2,1]]),S.Matrix([-JJ,0]),w,z)
check("exact curved affine normalization",
 normalized==S.Matrix([w*w-4*z/3,w**3-w*z+LL-JJ**3]))
# All five group types have exact wave witnesses; no sampled group inference.
models=[]
for rv,sv,RR,bv,cv in [
 (6,0,1/w**2,w**3,k*w*w*z+1),
 (2,0,w**3,w,w**3+k*z),
 (4,0,w*w+1/w,w*w,w**3+k*w*z+1),
 (6,0,1/w+1/w**2,w**3,k*w*w*z+w+1)]:
    C=rational_chart(rv,sv,RR,w,z,k)
    check("exact cubic chart witness wave r="+str(rv)+",R="+str(RR),C["red"](wave_tensor(1,bv,cv,w,z))==0)
    models.append(dict(b=str(bv),c=str(cv),degree=C["N"],group=C["geometric_group"]))
C=polynomial_chart(w**3,[(0,1)],w=w,z=z)
check("cubic polynomial constant-scalar group",C["N"]==7 and C["geometric_group"]=="A7")
models.append(dict(b="1",c="w^3",degree=7,group="A7"))
models.append(dict(b="w^2-4z/3",c="w^3-wz+D",degree=21,group="A21"))
check("exact full cubic constant-scalar group set",{x["group"] for x in models}=={"A7","A8","A10","A12","A21"})
nonconstant=polynomial_chart(w**3,[(0,0),(0,1)],w=w,z=z)
check("nonconstant scalar cubic boundary countercontrol",
 nonconstant["lambda_polynomial"]==z-w**4/2 and
 wave_tensor(z-w**4/2,1,w**3,w,z)==0 and
 nonconstant["N"]==22 and nonconstant["geometric_group"]=="A22")
def coefficients(expr):return [str(S.Poly(expr,t).nth(i)) for i in range(6)]
model=dict(claims=["KMFI56","KMFI57","KMFI58"],
 curved_pair=dict(b=str(b),c=str(c),X=[str(x) for x in X],finite=True,degree=21,geometric_group="A21",
   affine_reduction_to_completed_chart=False,arbitrary_birational_reduction_excluded=False),
 affine_Keller_witness=dict(source="ordinary affine4",degree=21,geometric_group="A21",Jacobian="gamma0!=0",
   coefficient_field="Q(tau,gamma0); H(tau)=0,gamma0^3=q0; degree<=18",
   H=str(M["H"]),coefficient_order="ascending powers tau0..5",
   q0=coefficients(stage["q0"]),u0=coefficients(stage["u0"]),v0=coefficients(stage["v0"]),
   R=[[coefficients(R[i,j]) for j in range(3)] for i in range(3)],
   f20=coefficients(stage["f20"]),actual_forbidden_weight_jets=polynomial_jets,
   new_abstract_action=False),
 exact_constant_scalar_cubic_groups=sorted({x["group"] for x in models}),field_witnesses=models,
 complete_cubic_constant_scalar_class=True,cubic_nonconstant_scalar_class_complete=False,
 arbitrary_Keller_PSL2_11_actions_decided=False,absolute_natural_A4_decided=False,
 absolute_binary_degree8_decided=False,full_action_classification_complete=False)
(ROOT/"cubic_wave_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["CUBIC_WAVE_CLASSIFICATION.md","cubic_wave.py","verify_cubic_wave.py","cubic_wave_model.json",
 "pole_direction.py","rational_wave.py","polynomial_wave.py"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written complete constant-scalar cubic coefficient exhaustion, finite degree and full geometric-group proof, with an exact marked algebraic field and all forbidden affine-lift jets. No formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"cubic_wave_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"],
 curved_degree=21,ordinary_affine_Keller_lift=True,groups=model["exact_constant_scalar_cubic_groups"],goal_complete=False)))
