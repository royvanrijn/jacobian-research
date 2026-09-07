from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Linearized deformation probe around the discriminant-3 E8^2 A2 K3.")
ap.add_argument("--mode",choices=["full","preserve-E8E8"],default="full")
args=ap.parse_args()

# Short Weierstrass K3:
# y^2=x^3+f(t,s)x+g(t,s), deg f=8, deg g=12.
#
# CM anchor:
# f0=0
# g0=t^5 s^5 (t-s)^2.
#
# First-order deformation f=eps*F, g=g0+eps*G.
# Since Delta=-16(4f^3+27g^2), f enters only at eps^3;
# first-order Delta variation is proportional to g0*G.
#
# This script reports the tangent dimensions under simple fiber-preservation
# constraints. The key point is that preserving both II* fibers forces huge
# vanishing at t=0,s=0 and leaves only a very small coefficient space.

print("X3TAN|stage=start|anchor=E8^2+A2")

# coefficient counts before quotienting base/scaling automorphisms
full_F=9   # binary octic
full_G=13  # binary degree 12
print(f"X3TAN|full_short_weierstrass_coeffs={full_F+full_G}|F={full_F}|G={full_G}")

# To retain II* at t=0 and s=0 in short Weierstrass form:
# ord f >=4 at each end -> F divisible t^4 s^4, hence 1-dimensional.
# ord g >=5 at each end -> G divisible t^5 s^5 times a binary quadric, 3-dimensional.
presF=1
presG=3
print(f"X3TAN|preserve_E8E8_coeffs={presF+presG}|F={presF}|G={presG}")
print("X3TAN|ansatz=F=c*t^4*s^4")
print("X3TAN|ansatz=G=t^5*s^5*(a*t^2+b*t*s+c*s^2)")

# Quotient by overall Weierstrass scaling and residual base automorphism fixing
# 0 and infinity (t/s -> lambda*t/s): two coordinate freedoms.
# Thus the E8E8-preserving family is expected to have dimension <=2,
# and imposing the rank-19/Shimura lattice polarization should cut it to 1.
print("X3TAN|coordinate_freedoms=2|expected_after_gauge<=2")
print("X3TAN|research_target=derive_one additional lattice-polarization condition to isolate X(6,79) curve locally")
print("X3TAN|stage=done")
