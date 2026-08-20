from sage.all import *
from pathlib import Path

OUT=Path("artifacts/local/elkies-k3")
OUT.mkdir(parents=True,exist_ok=True)

R.<t> = PolynomialRing(QQ)

def ord_at_zero(f):
    f=R(f)
    if f==0: return Infinity
    k=0
    while f[k]==0: k+=1
    return k

def ord_at(f,a):
    u=polygen(QQ,'u')
    return ord_at_zero(R(f(t+a)))

def infinity_order_homogeneous_degree(poly, total_degree):
    return total_degree - R(poly).degree()

print("X3VERIFY|stage=start")

# ------------------------------------------------------------
# Utsumi No.1: E8^2 A2, MW=0
#
# homogeneous:
#   y^2 = x^3 + t^5 s^5 (t-s)^2
#
# affine s=1:
#   y^2 = x^3 + t^5 (t-1)^2
# ------------------------------------------------------------
f1=R(0)
g1=t^5*(t-1)^2
Delta1=-16*(4*f1^3+27*g1^2)
print(f"X3VERIFY|no=1|f={f1}|g={g1}|Delta={factor(Delta1)}")
print(f"X3VERIFY|no=1|ordDelta_t0={ord_at_zero(Delta1)}|ordDelta_t1={ord_at(Delta1,1)}|ordDelta_inf={infinity_order_homogeneous_degree(Delta1,24)}")
assert ord_at_zero(Delta1)==10
assert ord_at(Delta1,1)==4
assert infinity_order_homogeneous_degree(Delta1,24)==10

# Fiber root lattices II*=E8, II*=E8, IV=A2.
# disc(E8)=1, disc(A2)=3, MW=0, torsion=0.
disc_trivial_1=ZZ(1*1*3)
mw_tors_1=ZZ(1)
disc_NS_1=-disc_trivial_1//(mw_tors_1^2)
assert abs(disc_NS_1)==3
print(f"X3VERIFY|no=1|fibers=II*,II*,IV|ADE=E8+E8+A2|MW=0|discNS={disc_NS_1}|discT={abs(disc_NS_1)}")

# ------------------------------------------------------------
# Utsumi No.2: D16 A2, MW torsion Z/2
#
# homogeneous equation from the literature:
# y^2=x^3
#   -3 t^2(t^6-16t^3s^3+16s^6)x
#   +2 t^3(t^3-2s^3)(t^6+32t^3s^3-32s^6)
#
# affine s=1.
# ------------------------------------------------------------
f2=-3*t^2*(t^6-16*t^3+16)
g2=2*t^3*(t^3-2)*(t^6+32*t^3-32)
Delta2=-16*(4*f2^3+27*g2^2)
print(f"X3VERIFY|no=2|f={f2}|g={g2}|Delta={factor(Delta2)}")
print(f"X3VERIFY|no=2|ordDelta_t0={ord_at_zero(Delta2)}|ordDelta_inf={infinity_order_homogeneous_degree(Delta2,24)}")
assert ord_at_zero(Delta2)==18
assert infinity_order_homogeneous_degree(Delta2,24)==3
q=Delta2//(t^18)
assert R(q).degree()==3
print(f"X3VERIFY|no=2|residualDelta={factor(q)}")

# D16 discriminant 4; A2 discriminant 3; torsion order 2.
disc_trivial_2=ZZ(4*3)
tors2=ZZ(2)
disc_NS_2=-disc_trivial_2//(tors2^2)
assert abs(disc_NS_2)==3
print(f"X3VERIFY|no=2|fibers=I12*,I3,3I1|ADE=D16+A2|MWtors=Z/2|discNS={disc_NS_2}|discT={abs(disc_NS_2)}")

# Explicit T target found independently in our Clifford transport.
Tcm=matrix(ZZ,[[2,1],[1,2]])
print(f"X3VERIFY|Tcm={Tcm}|det={Tcm.det()}")
assert Tcm.det()==3

summary=[
    "X3 singular K3, discriminant 3",
    "No1: y^2=x^3+t^5*s^5*(t-s)^2 ; fibers II*,II*,IV ; E8^2+A2 ; MW=0",
    "No2: y^2=x^3-3t^2(t^6-16t^3s^3+16s^6)x+2t^3(t^3-2s^3)(t^6+32t^3s^3-32s^6)",
    "     fibers I12*,I3,3I1 ; D16+A2 ; MW torsion Z/2",
    "Both give |disc(NS)|=3 and T=[[2,1],[1,2]].",
    "Recommendation: deform No1 first; use No2 as neighbor/cross-check."
]
(OUT/"disc3-x3-fibrations.txt").write_text("\n".join(summary)+"\n")
print("X3VERIFY|stage=done|status=OK|preferred=no1")
