from sage.all import *
import argparse


ap=argparse.ArgumentParser(description="Analyze and lift the unique GF(23) MW2 target.")
ap.add_argument("--show",action="store_true")
ap.add_argument("--hensel-digits",type=int,default=0)
ap.add_argument("--series-order",type=int,default=0)
ap.add_argument("--relation-degree",type=int,default=8)
ap.add_argument("--pair-relations",action="store_true")
ap.add_argument("--fixed-lam",default=None,help="rational lambda used for transverse Hensel lifting")
ap.add_argument("--scan-lam-height",type=int,default=0,
                help="scan reduced rational lambda values of bounded numerator/denominator height in the target residue disk")
ap.add_argument("--algdep-degree",type=int,default=0,
                help="run PARI p-adic algdep up to this degree on each lifted coordinate")
args=ap.parse_args()

p=23
names=["lam","mu","s","w","x1","r","n1","n2","n3","m2","m3","m4","m5"]
S=PolynomialRing(QQ,names,order="degrevlex")
SF=FractionField(S)
d=S.gens_dict()
T=PolynomialRing(SF,"t")
t=T.gen()


def V(name): return SF(d[name])
lam,mu,s,w,x1,r,n1,n2,n3,m2,m3,m4,m5=[V(name) for name in names]

# Rebuild the characteristic-zero normalized P1 chart using exactly the
# sparse eliminations used for the complete GF(23) scan.  A single temporary
# polynomial ring holds the active and eliminated coordinates.
aux_names=["C1","C2","A2","A3","A4","A5"]
E=PolynomialRing(QQ,names+aux_names,order="degrevlex")
EF=FractionField(E)
ed=E.gens_dict()
TE=PolynomialRing(EF,"tt")
tt=TE.gen()
elam,emu,es,ew,ex1=[EF(ed[name]) for name in ("lam","mu","s","w","x1")]
C1,C2,A2,A3,A4,A5=[EF(ed[name]) for name in aux_names]
ex2=ew**2-2*es-1-ex1
XU=1+ex1*tt+ex2*tt**2
YU=(tt-elam)*(tt-emu)*(1/(elam*emu)+C1*tt+C2*tt**2)
AU=A2*tt**2+A3*tt**3+A4*tt**4+A5*tt**5
BU=YU**2-XU**3-AU*XU


def linear_rhs(expression,variable):
    expression=EF(expression)
    coefficient=expression.derivative(variable)
    numerator=E(expression.numerator())
    if numerator.degree(E(variable))!=1:
        raise RuntimeError(f"expression is not linear in {variable}")
    coefficient=EF(numerator.derivative(E(variable)))
    return -EF(numerator.subs({E(variable):0}))/coefficient


c1U=linear_rhs(BU[1],C1)
BU1=TE([EF(coefficient).subs({E(C1):c1U}) for coefficient in BU.list()])
c2U=linear_rhs(BU1[2],C2)
subs12={E(C1):c1U,E(C2):c2U}
Xlam=XU(elam); Xmu=XU(emu)
rhs_lam=-3*Xlam**2-A2*elam**2-A5*elam**5
rhs_mu=-3*Xmu**2-A2*emu**2-A5*emu**5
node_matrix=matrix(EF,[[elam**3,elam**4],[emu**3,emu**4]])
node_solution=node_matrix.solve_right(vector(EF,[rhs_lam,rhs_mu]))
a3U,a4U=node_solution
subs1234={**subs12,E(A3):a3U,E(A4):a4U}


def settle_e(value,substitutions):
    value=EF(value)
    for _ in range(10):
        old=value
        value=EF(value.subs(substitutions))
        if value==old:
            break
    return value


Ystage=TE([settle_e(coefficient,subs1234) for coefficient in YU.list()])
a2U=linear_rhs(Ystage(1)-ew*(ew**2-3*es),A2)
subs12345={**subs1234,E(A2):a2U}
Astage=TE([settle_e(coefficient,subs12345) for coefficient in AU.list()])
a5U=linear_rhs(Astage(1)+3*es**2,A5)
full_subs={**subs12345,E(A5):a5U}

phi=E.hom([S(d[name]) for name in names]+[S(0)]*len(aux_names),S)


def down(value):
    value=settle_e(value,full_subs)
    return SF(phi(E(value.numerator())))/SF(phi(E(value.denominator())))


X1=T([down(coefficient) for coefficient in XU.list()])
Y1=T([down(coefficient) for coefficient in YU.list()])
A=T([down(coefficient) for coefficient in AU.list()])
B=Y1**2-X1**3-A*X1
Delta=-16*(4*A**3+27*B**2)

node_lam=X1(lam)
p1_residuals=[
    6*node_lam**2*A.derivative(t,2)(lam)
    +6*node_lam*B.derivative(t,2)(lam)-A.derivative(t)(lam)**2,
    s*A.derivative(t)(1)+B.derivative(t)(1),
    6*s**2*A.derivative(t,2)(1)+6*s*B.derivative(t,2)(1)-A.derivative(t)(1)**2,
]

# P2 has one pole, the nonzero D4 class, the inverse E6 class, and component
# 1 at the normalized I3.  Incidence at t=1 and the E6 leading sign are
# built into n4,m6,m7.
q=t-r
n4=s*(1-r)**2-n1-n2-n3
m7=-Y1[4]
m6=-m2-m3-m4-m5-m7
N=T([0,n1,n2,n3,n4])
M=T([0,0,m2,m3,m4,m5,m6,m7])
H=M**2-(N**3+A*N*q**4+B*q**6)
raw_residuals=p1_residuals+[H[i] for i in range(15)]
equations=[]
source_indices=[]
for index,residual in enumerate(raw_residuals):
    numerator=S(SF(residual).numerator())
    if numerator==0:
        continue
    equations.append(numerator)
    source_indices.append(index)

target={
    "lam":16,"mu":13,"s":16,"w":0,"x1":10,"r":8,
    "n1":15,"n2":6,"n3":19,"m2":22,"m3":14,"m4":13,"m5":2,
}
K=GF(p)
Sp=PolynomialRing(K,names,order="degrevlex")
point={Sp.gen(i):K(target[name]) for i,name in enumerate(names)}
mod_equations=[Sp(eq)*~Sp(eq.denominator()) if hasattr(eq,"denominator") else Sp(eq) for eq in equations]
# equations are already polynomial numerators over QQ; clear scalar
# denominators before reduction.
mod_equations=[]
integral_equations=[]
for equation in equations:
    denominator=lcm(QQ(c).denominator() for c in equation.coefficients())
    integral=S(denominator*equation)
    content=gcd(abs(ZZ(c)) for c in integral.coefficients())
    if content:
        integral=S(integral/content)
    integral_equations.append(integral)
    mod_equations.append(Sp(integral))

failures=[i for i,equation in enumerate(mod_equations) if equation.subs(point)]
if failures:
    raise RuntimeError(f"target fails equations {failures}")
jacobian=matrix(K,[
    [equation.derivative(variable).subs(point) for variable in Sp.gens()]
    for equation in mod_equations
])
rank=jacobian.rank()
print(
    f"MW2LIFT|p={p}|variables={len(names)}|equations={len(equations)}"
    f"|jacobian_rank={rank}|tangent_dimension={len(names)-rank}",flush=True,
)
print("MW2LIFT|sources="+",".join(map(str,source_indices)),flush=True)
print("MW2LIFT|target="+",".join(f"{name}:{target[name]}" for name in names),flush=True)
for index,vector0 in enumerate(jacobian.right_kernel().basis()):
    print(
        f"MW2LIFT|tangent={index}|"+
        ",".join(f"{name}:{int(value)}" for name,value in zip(names,vector0) if value),
        flush=True,
    )
if args.show:
    for row,(source,equation) in enumerate(zip(source_indices,integral_equations)):
        print(
            f"MW2LIFT_EQ|row={row}|source={source}|degree={equation.total_degree()}"
            f"|terms={len(equation.monomials())}",flush=True,
        )

fixed_index=names.index("lam")
unknown_indices=[i for i in range(len(names)) if i!=fixed_index]
unknown_jacobian=jacobian.matrix_from_columns(unknown_indices)
row_indices=list(unknown_jacobian.transpose().pivots())
if len(row_indices)!=len(unknown_indices):
    raise RuntimeError("fixing lam did not give a transverse square system")
square_jacobian=unknown_jacobian.matrix_from_rows(row_indices)


def structured_residuals(values,coefficient_ring):
    vals=dict(zip(names,values))
    l=vals["lam"]; m=vals["mu"]; ss=vals["s"]; ww=vals["w"]; xx1=vals["x1"]
    CT=PolynomialRing(coefficient_ring,"T")
    T0=CT.gen()
    xx2=ww**2-2*ss-1-xx1
    XX=1+xx1*T0+xx2*T0**2
    yy1=3*xx1/2
    cc1=(yy1+(l+m)/(l*m))/(l*m)
    target_y=ww*(ww**2-3*ss)
    base=(T0-l)*(T0-m)*(1/(l*m)+cc1*T0)
    c2_coefficient=(1-l)*(1-m)
    cc2=(target_y-base(1))/c2_coefficient
    YY=(T0-l)*(T0-m)*(1/(l*m)+cc1*T0+cc2*T0**2)
    provisional_B=YY**2-XX**3
    aa2=provisional_B[2]

    determinant=l**3*m**4-l**4*m**3
    def make_A(aa5):
        rhs_l=-3*XX(l)**2-aa2*l**2-aa5*l**5
        rhs_m=-3*XX(m)**2-aa2*m**2-aa5*m**5
        aa3=(rhs_l*m**4-l**4*rhs_m)/determinant
        aa4=(l**3*rhs_m-rhs_l*m**3)/determinant
        return aa2*T0**2+aa3*T0**3+aa4*T0**4+aa5*T0**5

    A0=make_A(coefficient_ring(0)); Aone=make_A(coefficient_ring(1))
    aa5=-(A0(1)+3*ss**2)/(Aone(1)-A0(1))
    AA=make_A(aa5)
    BB=YY**2-XX**3-AA*XX
    node=XX(l)
    P1=[
        6*node**2*AA.derivative(T0,2)(l)
        +6*node*BB.derivative(T0,2)(l)-AA.derivative(T0)(l)**2,
        ss*AA.derivative(T0)(1)+BB.derivative(T0)(1),
        6*ss**2*AA.derivative(T0,2)(1)
        +6*ss*BB.derivative(T0,2)(1)-AA.derivative(T0)(1)**2,
    ]
    rr=vals["r"]
    nn4=ss*(1-rr)**2-vals["n1"]-vals["n2"]-vals["n3"]
    NN=vals["n1"]*T0+vals["n2"]*T0**2+vals["n3"]*T0**3+nn4*T0**4
    mm7=-YY[4]
    mm6=-vals["m2"]-vals["m3"]-vals["m4"]-vals["m5"]-mm7
    MM=(vals["m2"]*T0**2+vals["m3"]*T0**3+vals["m4"]*T0**4
        +vals["m5"]*T0**5+mm6*T0**6+mm7*T0**7)
    qq=T0-rr
    HH=MM**2-(NN**3+AA*NN*qq**4+BB*qq**6)
    return P1+[HH[i] for i in range(15)]


if args.hensel_digits:
    if args.scan_lam_height:
        height=args.scan_lam_height
        fixed_lams=[]
        for denominator in range(1,height+1):
            if denominator%p==0:
                continue
            for numerator in range(-height,height+1):
                if gcd(numerator,denominator)!=1:
                    continue
                candidate=QQ(numerator)/denominator
                if K(candidate)==K(target["lam"]):
                    fixed_lams.append(candidate)
        fixed_lams=sorted(set(fixed_lams),key=lambda value:(max(abs(value.numerator()),value.denominator()),value))
        print(f"MW2HENSELSCAN|height={height}|candidates={len(fixed_lams)}",flush=True)
    else:
        fixed_lams=[QQ(args.fixed_lam) if args.fixed_lam is not None else QQ(target["lam"])]

    exact_hits=[]
    for fixed_lam in fixed_lams:
        if fixed_lam.denominator()%p==0 or K(fixed_lam)!=K(target["lam"]):
            raise RuntimeError("fixed lambda does not reduce to the target")
        coordinates=[ZZ(target[name]) for name in names]
        modulus=ZZ(p)
        for digit in range(1,args.hensel_digits):
            next_modulus=modulus*p
            coordinates[fixed_index]=ZZ(
                (fixed_lam.numerator()*inverse_mod(fixed_lam.denominator(),next_modulus))
                % next_modulus
            )
            rhs=[]
            for row in row_indices:
                value=ZZ(integral_equations[row](*coordinates))
                if value%modulus:
                    raise RuntimeError(f"row {row} fails modulo {modulus}")
                rhs.append(K(-(value//modulus)))
            correction=square_jacobian.solve_right(vector(K,rhs))
            for column,index in enumerate(unknown_indices):
                coordinates[index]+=modulus*ZZ(correction[column])
            modulus=next_modulus
            failures=[
                row for row,equation in enumerate(integral_equations)
                if ZZ(equation(*coordinates))%modulus
            ]
            if failures:
                raise RuntimeError(f"full Hensel system fails at digit {digit+1}: {failures}")
            if not args.scan_lam_height:
                print(f"MW2HENSEL|digits={digit+1}|modulus={modulus}",flush=True)
        reconstructions=[]
        for coordinate in coordinates:
            try:
                reconstructions.append(ZZ(coordinate%modulus).rational_reconstruction(modulus))
            except ArithmeticError:
                reconstructions.append(None)
        exact=all(value is not None for value in reconstructions)
        if exact:
            exact=all(equation(*reconstructions)==0 for equation in integral_equations)
        if not args.scan_lam_height or exact:
            print(
                "MW2HENSEL|reconstruction="+
                ",".join(f"{name}:{value}" for name,value in zip(names,reconstructions)),flush=True,
            )
            print(f"MW2HENSEL|fixed_lam={fixed_lam}|exact_rational_point={int(exact)}",flush=True)
        if args.algdep_degree:
            if args.scan_lam_height:
                raise RuntimeError("--algdep-degree is intended for one --fixed-lam lift")
            Kp=Qp(p,args.hensel_digits)
            ZP=PolynomialRing(ZZ,"z")
            for name,coordinate in zip(names,coordinates):
                value=Kp(ZZ(coordinate%modulus)).add_bigoh(args.hensel_digits)
                polynomial=ZP(pari(value).algdep(args.algdep_degree))
                if polynomial.leading_coefficient()<0:
                    polynomial=-polynomial
                evaluation=polynomial(value)
                valuation=(args.hensel_digits if evaluation==0 else min(args.hensel_digits,evaluation.valuation()))
                coefficient_height=max(abs(coefficient) for coefficient in polynomial)
                print(
                    f"MW2ALGDEP|fixed_lam={fixed_lam}|coordinate={name}"
                    f"|degree={polynomial.degree()}|height={coefficient_height}"
                    f"|valuation={valuation}|polynomial={polynomial}",flush=True,
                )
        if not exact:
            continue
        exact_hits.append((fixed_lam,tuple(reconstructions)))
        exact_point={S.gen(i):QQ(value) for i,value in enumerate(reconstructions)}
        def exact_polynomial(poly):
            QR=PolynomialRing(QQ,"T")
            return QR([QQ(SF(coefficient).subs(exact_point)) for coefficient in poly.list()])
        exact_A=exact_polynomial(A)
        exact_B=exact_polynomial(B)
        exact_X1=exact_polynomial(X1)
        exact_Y1=exact_polynomial(Y1)
        exact_N=exact_polynomial(N)
        exact_M=exact_polynomial(M)
        exact_Delta=-16*(4*exact_A**3+27*exact_B**2)
        exact_q=exact_A.parent().gen()-fixed_lam*0-QQ(reconstructions[names.index("r")])
        pair_gcd=gcd(exact_X1*exact_q**2-exact_N,exact_Y1*exact_q**3-exact_M)
        for label,poly in (
            ("A",exact_A),("B",exact_B),("P1X",exact_X1),("P1Y",exact_Y1),
            ("N",exact_N),("M",exact_M),
        ):
            print(f"MW2RATIONAL|{label}="+",".join(map(str,poly.list())),flush=True)
        print(f"MW2RATIONAL|Delta_factor={factor(exact_Delta)}",flush=True)
        print(f"MW2RATIONAL|pair_gcd={pair_gcd}",flush=True)
    if args.scan_lam_height:
        print(
            f"MW2HENSELSCAN|status=done|height={args.scan_lam_height}"
            f"|exact_hits={len(exact_hits)}|lambdas="
            +",".join(str(item[0]) for item in exact_hits),flush=True,
        )


if args.series_order:
    if args.series_order<5:
        raise SystemExit("--series-order must be at least 5")
    PS=PowerSeriesRing(K,"z",default_prec=args.series_order)
    z=PS.gen()
    series=[PS(target[name]) for name in names]
    series[fixed_index]+=z
    raw_jacobian_rows=[]
    for residual in raw_residuals:
        fraction=SF(residual)
        numerator=S(fraction.numerator()); denominator=S(fraction.denominator())
        den_value=Sp(denominator).subs(point)
        raw_jacobian_rows.append([
            Sp(numerator.derivative(variable)).subs(point)/den_value
            for variable in S.gens()
        ])
    raw_jacobian=matrix(K,raw_jacobian_rows)
    raw_unknown=raw_jacobian.matrix_from_columns(unknown_indices)
    raw_rows=list(raw_unknown.transpose().pivots())
    if len(raw_rows)!=len(unknown_indices):
        raise RuntimeError("structured residuals lost transverse rank")
    raw_square=raw_unknown.matrix_from_rows(raw_rows)
    for degree in range(1,args.series_order):
        residual_values=structured_residuals(series,PS)
        correction=raw_square.solve_right(vector(K,[-PS(residual_values[row])[degree] for row in raw_rows]))
        for column,index in enumerate(unknown_indices):
            series[index]+=PS(correction[column])*z**degree
    failures=[]
    final_residuals=structured_residuals(series,PS)
    for row,value in enumerate(final_residuals):
        value=PS(value)
        if value and value.valuation()<args.series_order:
            failures.append((row,value.valuation()))
    if failures:
        raise RuntimeError(f"formal series fails: {failures}")
    print(f"MW2SERIES|order={args.series_order}|verified={len(final_residuals)}",flush=True)
    parameter=series[fixed_index]
    for index,name in enumerate(names):
        if index==fixed_index:
            continue
        value=series[index]
        found=None
        for y_degree in range(1,args.relation_degree+1):
            for x_degree in range(0,args.relation_degree+1):
                count=(x_degree+1)*(y_degree+1)
                if count>args.series_order-4:
                    continue
                monomials=[parameter**i*value**j for j in range(y_degree+1) for i in range(x_degree+1)]
                coefficient_matrix=matrix(K,[
                    [term[k] for term in monomials] for k in range(args.series_order)
                ])
                kernel=coefficient_matrix.right_kernel()
                if kernel.dimension()!=1:
                    continue
                relation=kernel.basis()[0]
                if all(relation[y_degree*(x_degree+1)+i]==0 for i in range(x_degree+1)):
                    continue
                found=(x_degree,y_degree,relation)
                break
            if found:
                break
        if found:
            xd,yd,relation=found
            print(
                f"MW2SERIES|relation={name}|xdegree={xd}|ydegree={yd}|coeffs="
                +",".join(map(str,map(int,relation))),flush=True,
            )
    if args.pair_relations:
        for left in range(len(names)):
            for right in range(left+1,len(names)):
                xs=series[left]; ys=series[right]
                found=None
                for x_degree in range(1,args.relation_degree+1):
                    for y_degree in range(1,args.relation_degree+1):
                        count=(x_degree+1)*(y_degree+1)
                        if count>args.series_order-4:
                            continue
                        monomials=[xs**i*ys**j for j in range(y_degree+1) for i in range(x_degree+1)]
                        coefficient_matrix=matrix(K,[
                            [term[k] for term in monomials] for k in range(args.series_order)
                        ])
                        kernel=coefficient_matrix.right_kernel()
                        if kernel.dimension()!=1:
                            continue
                        relation=kernel.basis()[0]
                        if all(relation[j*(x_degree+1)+x_degree]==0 for j in range(y_degree+1)):
                            continue
                        if all(relation[y_degree*(x_degree+1)+i]==0 for i in range(x_degree+1)):
                            continue
                        found=(x_degree,y_degree,relation)
                        break
                    if found:
                        break
                if found:
                    xd,yd,relation=found
                    print(
                        f"MW2PAIR|left={names[left]}|right={names[right]}"
                        f"|xdegree={xd}|ydegree={yd}|coeffs="
                        +",".join(map(str,map(int,relation))),flush=True,
                    )
