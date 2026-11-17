#!/usr/bin/env python3
"""Symbolic audit of the quadratic-stabilized intersective transfer."""

import sympy as sp


A, B, C, u, v, d, e, k = sp.symbols("A B C u v d e k", nonzero=True)

# Endpoint data for g:
# g(0)=u, g(1)=v, g'(0)=d, g'(1)=e, g''(1)=k.
# For R=A W^2+B W+C, the tangent equation for P=gR is L=0.
L = v * A + (v - u) * B + (v - u - d) * C

# Solve the tangent equation for A and audit the binary discriminant.
A_on_L = sp.cancel(-((v - u) * B + (v - u - d) * C) / v)
Delta = sp.factor((B**2 - 4 * A * C).subs(A, A_on_L))
assert Delta == sp.factor(
    B**2 + 4 * (v - u) * B * C / v + 4 * (v - u - d) * C**2 / v
)

# The discriminant binary form is a square exactly on vu-vd-u^2=0.
mixed = sp.expand(4 * (v - u) / v)
constant = sp.expand(4 * (v - u - d) / v)
assert sp.factor(constant - mixed**2 / 4) == sp.factor(
    4 * (v * u - v * d - u**2) / v**2
)

# First weighted nondegeneracy:
# M=P'(1)-P'(0).
M = (
    (e + 2 * v) * A
    + (e + v - u) * B
    + (e - d) * C
)

# Second weighted nondegeneracy:
# N=P''(1)-2(P'(1)-P'(0)).
N = (
    (k + 2 * e - 2 * v) * A
    + (k - 2 * v + 2 * u) * B
    + (k - 2 * e + 2 * d) * C
)

# On the explicit direction C=0 in the tangent plane, the two restrictions
# reduce to the endpoint expressions used in the proof.
direction = {C: 0, A: -(v - u) * B / v}
assert sp.factor(M.subs(direction)) == sp.factor(
    B * (u * (v + e) - v**2) / v
)
assert sp.factor(N.subs(direction)) == sp.factor(
    B * (u * (2 * e + k) - 2 * v * e) / v
)

print("PASS: tangent normalization is one linear equation in A,B,C")
print("PASS: the auxiliary discriminant is nonsquare off vu-vd-u^2=0")
print("PASS: both weighted exceptional loci restrict to proper linear forms")
