from sage.all import *
from pathlib import Path
import argparse,time

ap=argparse.ArgumentParser(description="Track-B MW3 P1 hard-coded triangular stage 2.")
ap.add_argument("--p",type=int,default=31)
ap.add_argument("--checkpoint",default="/tmp/mw3-p1-stage2.txt")
ap.add_argument("--save-state",default=None)
ap.add_argument("--probe",action="store_true")
ap.add_argument("--probe-tag",action="append",default=[])
ap.add_argument("--show-probe",action="store_true")
ap.add_argument("--factor-probe",action="store_true")
ap.add_argument("--discriminant-var",default=None)
ap.add_argument("--show-discriminant",action="store_true")
ap.add_argument("--parametrize-p1-2",action="store_true")
ap.add_argument("--parametrize-p1-10",action="store_true")
ap.add_argument("--parametrize-i2-1",action="store_true")
ap.add_argument("--rank-i2-pairs",action="store_true")
ap.add_argument("--solve-i2-pair",action="store_true")
ap.add_argument("--preserve-a",action="store_true")
ap.add_argument("--keep-a6",action="store_true")
ap.add_argument("--rank-i2-five",action="store_true")
ap.add_argument("--eliminate-i2-b2",action="store_true")
ap.add_argument("--eliminate-i2-a3",action="store_true")
ap.add_argument("--i2-pivot",choices=["a3","a4","b2","b9","x2","y4","y5"],default=None)
ap.add_argument("--eliminate-i2-lam-b1",action="store_true")
args=ap.parse_args()

K=GF(args.p)

names=[]
names += [f"a{i}" for i in range(9)]
names += [f"b{i}" for i in range(13)]
names += ["lam","sinf","s0","s1","sl"]
names += ["x1","x2","x3"]
names += ["y1","y2","y3","y4","y5"]
if args.parametrize_p1_2:
    names += ["r0"]
if args.parametrize_p1_10:
    names += ["r1"]
if args.parametrize_i2_1:
    names += ["r2"]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
Rt=PolynomialRing(RF,"t"); t=Rt.gen()
V=lambda n: RF(d[n])

aa=[V(f"a{i}") for i in range(9)]
bb=[V(f"b{i}") for i in range(13)]
lam,sinf,s0,s1,sl=[V(n) for n in ["lam","sinf","s0","s1","sl"]]
x1,x2,x3=[V(n) for n in ["x1","x2","x3"]]
y1,y2,y3,y4,y5=[V(n) for n in ["y1","y2","y3","y4","y5"]]
r0=V("r0") if args.parametrize_p1_2 else None
r1=V("r1") if args.parametrize_p1_10 else None
r2=V("r2") if args.parametrize_i2_1 else None

A=sum(aa[i]*t**i for i in range(9))
B=sum(bb[i]*t**i for i in range(13))
X=s0+x1*t+x2*t**2+x3*t**3+sinf*t**4
Y=y1*t+y2*t**2+y3*t**3+y4*t**4+y5*t**5
S=Y**2-X**3-A*X-B
Delta=-16*(4*A**3+27*B**2)

eq={}
def add(tag,e):
    e=RF(e)
    if e!=0: eq[tag]=e

add("I3_0_A",A(0)+3*s0**2)
add("I3_0_B",B(0)-2*s0**3)
add("I3_0_d1",Delta[1])
add("I3_0_d2",Delta[2])

add("I2_1_A",A(1)+3*s1**2)
add("I2_1_B",B(1)-2*s1**3)
add("I2_1_d1",Delta.derivative(t)(1))

add("I2_lam_A",A(lam)+3*sl**2)
add("I2_lam_B",B(lam)-2*sl**3)
add("I2_lam_d1",Delta.derivative(t)(lam))
add("P1_lam_X",X(lam)-sl)
add("P1_lam_Y",Y(lam))

for k in range(24,13,-1):
    add(f"I11_inf_D{k}",Delta[k])

add("I11_inf_A",aa[8]+3*sinf**2)
add("I11_inf_B",bb[12]-2*sinf**3)

for k in range(S.degree()+1):
    if S[k]!=0: add(f"P1_{k}",S[k])

subs={}
history=[]
active=[V(n) for n in names]

def settle(fr,passes=20):
    fr=RF(fr)
    for _ in range(passes):
        old=fr
        fr=RF(fr.subs(subs))
        if fr==old:
            break
    return fr

def eliminate(vname,tag):
    t0=time.time()
    v=V(vname)
    vv=R(v)
    e=settle(eq[tag])
    num=R(e.numerator())

    deg=num.degree(vv)
    if deg!=1:
        raise RuntimeError(f"{tag} is degree {deg} in {vname}, expected 1")

    c=RF(num.derivative(vv))
    if c==0:
        raise RuntimeError(f"{tag} derivative wrt {vname} vanished")

    # Verify affine-linear.
    if R(c.numerator()).degree(vv)!=0 or R(c.denominator()).degree(vv)!=0:
        raise RuntimeError(f"{tag} coefficient still depends on {vname}")

    rhs=RF(-num.subs({vv:0})/c)
    subs[v]=rhs
    active.remove(v)
    history.append((vname,tag,rhs))

    dt=time.time()-t0
    print(
        f"MW3S2|elim={len(history)}|var={vname}|from={tag}"
        f"|degree={num.total_degree()}|terms={len(num.monomials())}"
        f"|rhs_num_degree={R(rhs.numerator()).total_degree()}"
        f"|rhs_num_terms={len(R(rhs.numerator()).monomials())}"
        f"|rhs_den_degree={R(rhs.denominator()).total_degree()}"
        f"|rhs_den_terms={len(R(rhs.denominator()).monomials())}"
        f"|seconds={dt:.3f}",
        flush=True
    )

# Exact stage-1 chain already verified experimentally.
stage1=[
 ("a0","I3_0_A"),
 ("a1","P1_1"),
 ("a8","I11_inf_A"),
 ("a7","P1_11"),
 ("b0","I3_0_B"),
 ("b12","I11_inf_B"),
 ("a2","I3_0_d2"),
 (("b10" if args.preserve_a else "a6"),"I11_inf_D22"),
 (("b9" if args.preserve_a else "a5"),"I11_inf_D21"),
 ("y1","P1_lam_Y"),
 ("sinf","P1_lam_X"),
]

# New deterministic section-equation chain.
stage2=[
 ("b3","P1_3"),
 ("b4","P1_4"),
 ("b5","P1_5"),
 ("b6","P1_6"),
 ("b7","P1_7"),
 ("b8","P1_8"),
]
if not args.preserve_a:
    stage2.append(("b10","P1_9"))
else:
    if not args.keep_a6:
        stage2.append(("a6","P1_9"))

for v,tag in stage1+stage2:
    eliminate(v,tag)

if args.parametrize_p1_2:
    p12_linear = lam*y2 + lam**2*y3 + lam**3*y4 + lam**4*y5
    subs[s0] = r0**2 / 3
    subs[V("b1")] = 6*s0**2*x1 + 2*r0*s0*p12_linear
    active.remove(s0)
    active.remove(V("b1"))
    history.append(("s0", "P1_2_param", subs[s0]))
    history.append(("b1", "P1_2_param", subs[V("b1")]))
    if settle(eq["P1_2"]) != 0:
        raise ArithmeticError("P1_2 norm parametrization failed")
    print(
        "MW3S2|param=P1_2_norm|removed=s0,b1|added=r0|verified=1",
        flush=True,
    )

if args.rank_i2_five:
    if not (args.preserve_a and args.parametrize_i2_1):
        raise SystemExit("--rank-i2-five requires --preserve-a and --parametrize-i2-1")
    x_at_one_equation=settle(X(1)-(r2**2-2*s1))
    y_at_one_equation=settle(Y(1)-r2*(r2**2-3*s1))
    a_prime=settle(A.derivative(t)(1))
    x_prime=settle(X.derivative(t)(1))
    y_prime=settle(Y.derivative(t)(1))
    tate_equation=RF(2*r2*y_prime-3*(r2**2-s1)*x_prime-a_prime)
    five_equations=[
        R(settle(eq["I2_1_A"]).numerator()),
        R(settle(eq["I2_lam_A"]).numerator()),
        R(x_at_one_equation.numerator()),
        R(y_at_one_equation.numerator()),
        R(tate_equation.numerator()),
    ]
    five_variables=[R(d[name]) for name in ("a3","a4","a5","b2","x2","y2")]
    for omitted in five_variables:
        pivot_variables=[variable for variable in five_variables if variable!=omitted]
        coefficient_matrix=matrix(
            R,
            [
                [equation.derivative(variable) for variable in pivot_variables]
                for equation in five_equations
            ],
        )
        constants=vector(
            R,
            [
                equation.subs({variable:0 for variable in pivot_variables})
                for equation in five_equations
            ],
        )
        determinant=coefficient_matrix.det()
        numerator_terms=[]
        numerator_degrees=[]
        for column in range(5):
            replaced=matrix(R,coefficient_matrix)
            replaced.set_column(column,-constants)
            numerator=replaced.det()
            numerator_terms.append(len(numerator.monomials()))
            numerator_degrees.append(numerator.total_degree())
        print(
            f"MW3S2FIVE|omit={omitted}|"
            f"det_degree={determinant.total_degree()}|"
            f"det_terms={len(determinant.monomials())}|"
            f"num_degrees={','.join(map(str,numerator_degrees))}|"
            f"num_terms={','.join(map(str,numerator_terms))}",
            flush=True,
        )

if args.parametrize_p1_10:
    if not args.parametrize_p1_2:
        raise SystemExit("--parametrize-p1-10 requires --parametrize-p1-2")
    b11_variable = R(d["b11"])
    subs[sl] = lam**3*x3 + lam**2*x2 + lam*x1 + (r0**2-r1**2)/3
    p110 = R(settle(eq["P1_10"]).numerator())
    if p110.degree(b11_variable) != 2:
        raise ArithmeticError("P1_10 is not quadratic in b11")
    quadratic_a = RF(p110.derivative(b11_variable, 2) / 2)
    quadratic_b = RF(p110.derivative(b11_variable).subs({b11_variable: 0}))
    discriminant = RF(p110.discriminant(b11_variable))
    square_unit = RF(
        discriminant / (y5**2 * lam**20 * r1**6)
    )
    if square_unit.denominator() not in K or square_unit.numerator() not in K:
        raise ArithmeticError("P1_10 discriminant quotient is not constant")
    square_unit = K(square_unit)
    if not square_unit.is_square():
        raise ArithmeticError("P1_10 discriminant unit is not a square")
    square_root_unit = square_unit.sqrt()
    quadratic_root = RF(
        (-quadratic_b + square_root_unit*y5*lam**10*r1**3)
        / (2*quadratic_a)
    )
    subs[V("b11")] = settle(quadratic_root)
    active.remove(sl)
    active.remove(V("b11"))
    history.append(("sl", "P1_10_param", subs[sl]))
    history.append(("b11", "P1_10_param", subs[V("b11")]))
    if settle(eq["P1_10"]) != 0:
        raise ArithmeticError("P1_10 norm parametrization failed")
    print(
        f"MW3S2|param=P1_10_norm|removed=sl,b11|added=r1|"
        f"square_unit={square_unit}|verified=1",
        flush=True,
    )

if args.rank_i2_pairs:
    if not args.parametrize_p1_10:
        raise SystemExit("--rank-i2-pairs requires both norm parametrizations")
    linear_variables=[R(d[name]) for name in ("a3","a4","b2","b9")]
    linear_equations=[
        R(settle(eq[tag]).numerator())
        for tag in ("I2_1_A","I2_lam_A")
    ]
    for first_index in range(len(linear_variables)):
        for second_index in range(first_index+1,len(linear_variables)):
            first_variable=linear_variables[first_index]
            second_variable=linear_variables[second_index]
            coefficients=matrix(
                R,
                [
                    [
                        equation.derivative(first_variable),
                        equation.derivative(second_variable),
                    ]
                    for equation in linear_equations
                ],
            )
            constants=vector(
                R,
                [
                    equation.subs({first_variable:0,second_variable:0})
                    for equation in linear_equations
                ],
            )
            determinant=coefficients.det()
            numerator_first=(
                -constants[0]*coefficients[1,1]
                + coefficients[0,1]*constants[1]
            )
            numerator_second=(
                -coefficients[0,0]*constants[1]
                + constants[0]*coefficients[1,0]
            )
            print(
                f"MW3S2PAIR|vars={first_variable},{second_variable}|"
                f"det_degree={determinant.total_degree()}|"
                f"det_terms={len(determinant.monomials())}|"
                f"num1_degree={numerator_first.total_degree()}|"
                f"num1_terms={len(numerator_first.monomials())}|"
                f"num2_degree={numerator_second.total_degree()}|"
                f"num2_terms={len(numerator_second.monomials())}",
                flush=True,
            )

if args.solve_i2_pair:
    if not args.parametrize_p1_10:
        raise SystemExit("--solve-i2-pair requires both norm parametrizations")
    pair_variables=[R(d["a3"]),R(d["b2"])]
    pair_equations=[
        R(settle(eq[tag]).numerator())
        for tag in ("I2_1_A","I2_lam_A")
    ]
    pair_coefficients=matrix(
        R,
        [
            [equation.derivative(variable) for variable in pair_variables]
            for equation in pair_equations
        ],
    )
    pair_constants=vector(
        R,
        [
            equation.subs({pair_variables[0]:0,pair_variables[1]:0})
            for equation in pair_equations
        ],
    )
    pair_determinant=pair_coefficients.det()
    pair_numerators=[
        -pair_constants[0]*pair_coefficients[1,1]
        + pair_coefficients[0,1]*pair_constants[1],
        -pair_coefficients[0,0]*pair_constants[1]
        + pair_constants[0]*pair_coefficients[1,0],
    ]
    pair_rhs=[RF(numerator/pair_determinant) for numerator in pair_numerators]
    simultaneous={variable:value for variable,value in zip(pair_variables,pair_rhs)}
    for equation in pair_equations:
        if RF(equation).subs(simultaneous) != 0:
            raise ArithmeticError("simultaneous I2 A solve failed")
    for variable,value in zip(pair_variables,pair_rhs):
        subs[RF(variable)]=value
        active.remove(RF(variable))
        history.append((str(variable),"I2_A_pair",value))
    print(
        f"MW3S2|solve=I2_A_pair|vars=a3,b2|"
        f"det_degree={pair_determinant.total_degree()}|"
        f"det_terms={len(pair_determinant.monomials())}|"
        f"num_terms={len(pair_numerators[0].monomials())},"
        f"{len(pair_numerators[1].monomials())}|verified=1",
        flush=True,
    )

if args.parametrize_i2_1:
    if not args.parametrize_p1_10:
        raise SystemExit("--parametrize-i2-1 requires both norm parametrizations")
    x2_variable=V("x2")
    y2_variable=V("y2")
    x_at_one_without_x2=settle(X(1)-x2_variable)
    y_at_one_without_y2=settle(Y(1)-y2_variable)
    subs[x2_variable]=r2**2-2*s1-x_at_one_without_x2
    subs[y2_variable]=r2*(r2**2-3*s1)-y_at_one_without_y2
    active.remove(x2_variable)
    active.remove(y2_variable)
    history.append(("x2","I2_1_point_param",subs[x2_variable]))
    history.append(("y2","I2_1_point_param",subs[y2_variable]))
    x_at_one=settle(X(1))
    y_at_one=settle(Y(1))
    if x_at_one != r2**2-2*s1 or y_at_one != r2*(r2**2-3*s1):
        raise ArithmeticError("I2@1 point parametrization failed")
    if y_at_one**2 != (x_at_one-s1)**2*(x_at_one+2*s1):
        raise ArithmeticError("I2@1 nodal cubic parametrization failed")
    print(
        "MW3S2|param=I2_1_identity_point|removed=x2,y2|added=r2|"
        "verified=1",
        flush=True,
    )

if sum(bool(value) for value in (args.eliminate_i2_b2, args.eliminate_i2_a3, args.i2_pivot, args.solve_i2_pair)) > 1:
    raise SystemExit("choose at most one I2_1_A pivot")
if args.eliminate_i2_b2:
    eliminate("b2", "I2_1_A")
if args.eliminate_i2_a3:
    eliminate("a3", "I2_1_A")
if args.i2_pivot:
    eliminate(args.i2_pivot, "I2_1_A")
if args.eliminate_i2_lam_b1:
    if not args.eliminate_i2_b2:
        raise SystemExit("--eliminate-i2-lam-b1 requires --eliminate-i2-b2")
    eliminate("b1", "I2_lam_A")

remaining=[str(v) for v in active]
print(f"MW3S2|CHECKPOINT|eliminated={len(history)}|remaining={len(remaining)}",flush=True)
print("MW3S2|remaining="+",".join(remaining),flush=True)

cp=Path(args.checkpoint)
with cp.open("w") as h:
    h.write(f"p={args.p}\n")
    h.write("remaining="+",".join(remaining)+"\n")
    h.write("\nELIMINATIONS\n")
    for i,(v,tag,rhs) in enumerate(history,1):
        h.write(f"{i}|{v}|{tag}|{rhs}\n")

print(f"MW3S2|saved={cp}",flush=True)
if args.save_state:
    save(
        {
            "ring":R,
            "equations":eq,
            "substitutions":subs,
            "remaining":remaining,
            "prime":args.p,
        },
        args.save_state,
    )
    print(f"MW3S2|state_saved={args.save_state}",flush=True)

# Probe only a small hand-picked set after all 18 eliminations.
# Never globally expand the entire residual system.
if args.probe:
    probe_tags=[
        "P1_2","P1_10",
        "I2_1_A","I2_1_B","I2_1_d1",
        "I2_lam_A","I2_lam_B","I2_lam_d1",
        "I11_inf_D20","I11_inf_D19","I11_inf_D18",
        "I11_inf_D17","I11_inf_D16","I11_inf_D15","I11_inf_D14",
    ]
    if args.probe_tag:
        probe_tags=args.probe_tag
    for tag in probe_tags:
        if tag not in eq: continue
        t0=time.time()
        z=settle(eq[tag])
        dt=time.time()-t0
        if z==0:
            print(f"MW3S2PROBE|tag={tag}|zero=1|seconds={dt:.3f}",flush=True)
            continue
        num=R(z.numerator())
        lin=[]
        vars_here=[]
        for n in remaining:
            vv=R(d[n])
            dg=num.degree(vv)
            if dg:
                vars_here.append(f"{n}:{dg}")
                if dg==1 and num.derivative(vv).degree(vv)==0:
                    lin.append(n)
        print(
            f"MW3S2PROBE|tag={tag}|seconds={dt:.3f}"
            f"|degree={num.total_degree()}|terms={len(num.monomials())}"
            f"|linear={','.join(lin)}|vars={','.join(vars_here)}",
            flush=True
        )
        if args.show_probe:
            print(f"MW3S2PROBE_EXPR|tag={tag}|expr={num}",flush=True)
        if args.factor_probe:
            factors=num.factor()
            print(
                f"MW3S2PROBE_FACTOR|tag={tag}|unit={factors.unit()}|"
                + "factors="
                + ";".join(
                    f"degree:{factor.total_degree()},terms:{len(factor.monomials())},exp:{exponent}"
                    for factor,exponent in factors
                ),
                flush=True,
            )
        if args.discriminant_var:
            discriminant_variable=R(d[args.discriminant_var])
            if num.degree(discriminant_variable)!=2:
                raise ArithmeticError("discriminant variable is not quadratic")
            discriminant=num.discriminant(discriminant_variable)
            discriminant_factors=discriminant.factor()
            print(
                f"MW3S2PROBE_DISCRIMINANT|tag={tag}|var={args.discriminant_var}|"
                f"degree={discriminant.total_degree()}|terms={len(discriminant.monomials())}|"
                f"unit={discriminant_factors.unit()}|factors="
                + ";".join(
                    f"degree:{factor.total_degree()},terms:{len(factor.monomials())},exp:{exponent}"
                    for factor,exponent in discriminant_factors
                ),
                flush=True,
            )
            if args.show_discriminant:
                print(
                    f"MW3S2PROBE_DISCRIMINANT_EXPR|tag={tag}|"
                    f"factorization={discriminant_factors}",
                    flush=True,
                )

print("MW3S2|done",flush=True)
