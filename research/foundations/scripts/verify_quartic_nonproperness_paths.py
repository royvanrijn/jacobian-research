#!/usr/bin/env python3
"""Exact escape paths for both components of the quartic nonproperness set."""

import sympy as sp


e = sp.symbols("e", nonzero=True)
A, B, C, W, r = sp.symbols("A B C W r")

Delta = (
    27*(B*C)**4 - 36*(B*C)**2*(A*C**2) - (B*C)**2
    + 16*(A*C**2)**3 + 8*(A*C**2)**2 + A*C**2
)
Q4 = (
    A - B**2
    + C**2*(27*B**4 - 36*A*B**2 + 8*A**2)
    + 16*A**3*C**4
)
assert sp.factor(Delta-C**2*Q4) == 0


def reconstruction(root, target_a, target_b, target_c):
    gamma = target_b*target_c - root + 2*root**3
    x = target_c/gamma
    u = root/gamma
    v = (u-1)/3
    y = v/x
    s = 1-4*v-gamma
    z = s/x**2
    inverse = root**2-root**4-2*target_b*target_c*root+target_a*target_c**2
    assert sp.factor(inverse) == 0
    return tuple(map(sp.factor, (x, y, z))), sp.factor(gamma)


# COMPONENT 1: the entire plane C=0. Follow the simple primitive root W=-1.
# The chosen second-order jet makes the moving target tend to arbitrary
# (A,B,0), while the reconstructed z coordinate grows like 2/e^2.
L = (3*B**2-A)/2
W_plane = -1-B*e+L*e**2
A_plane = sp.cancel((W_plane**4-W_plane**2+2*B*e*W_plane)/e**2)
source_plane, gamma_plane = reconstruction(W_plane, A_plane, B, e)
assert sp.limit(A_plane, e, 0) == A
assert sp.limit(gamma_plane, e, 0) == -1
assert sp.limit(source_plane[0], e, 0) == 0
assert sp.limit(e**2*source_plane[2], e, 0) == 2


# COMPONENT 2, C!=0: parameterize every repeated-root discriminant point by r.
# Perturb the root, change only A to keep it a root, and obtain x=C/gamma -> oo.
s_boundary = r-2*r**3
t_boundary = r**2-3*r**4
B_boundary = s_boundary/C
A_boundary = t_boundary/C**2
W_boundary = r+e
A_near = sp.cancel((W_boundary**4-W_boundary**2
                    + 2*B_boundary*C*W_boundary)/C**2)
source_boundary, gamma_boundary = reconstruction(
    W_boundary, A_near, B_boundary, C
)
assert sp.factor(Q4.subs({A: A_boundary, B: B_boundary})) == 0
assert sp.limit(A_near, e, 0) == A_boundary
assert sp.limit(gamma_boundary, e, 0) == 0
assert sp.factor(gamma_boundary) != 0
assert sp.limit(1/source_boundary[0], e, 0) == 0


# At the intersection C=0, Q4=0 means A=B^2. This path approaches the same
# point through the discriminant component (not merely through the plane).
C_intersection = e
r_intersection = B*e
B_intersection = sp.cancel(
    (r_intersection-2*r_intersection**3)/C_intersection
)
W_intersection = r_intersection+e**2
A_intersection = sp.cancel(
    (W_intersection**4-W_intersection**2
     + 2*B_intersection*C_intersection*W_intersection)
    / C_intersection**2
)
source_intersection, gamma_intersection = reconstruction(
    W_intersection, A_intersection, B_intersection, C_intersection
)
assert sp.limit(B_intersection, e, 0) == B
assert sp.limit(A_intersection, e, 0) == B**2
assert sp.limit(gamma_intersection, e, 0) == 0
assert sp.limit(1/source_intersection[0], e, 0) == 0


print("PASS: Delta(BC,AC^2)=C^2 Q4 identifies the saturated component")
print("PASS: W=-1 gives an exact escape path to every target on C=0")
print("PASS: repeated roots give exact escape paths over Q4=0, C!=0")
print("PASS: Q4 paths extend to the intersection A=B^2, C=0")
