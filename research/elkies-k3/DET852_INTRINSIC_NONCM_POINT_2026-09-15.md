# Determinant 852: an intrinsic rational non-CM period

The subsequent [global root obstruction](DET852_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md)
now excludes every MW17 fibration on this NS. Its full rational marking
and non-CM period remain valid; any rootless UNKNOWN below belongs to
this earlier gate alone.

The subsequent [full rational NS theorem](DET852_RATIONAL_MARKING_SOURCE_2026-09-15.md)
now completes primitive lattice and actual divisor descent. Remaining
admission gaps stated below belong to this earlier gate alone.

Let C=X_0^6(71)/<w426>, the full marked curve identified with426b1.
Choose any of its four rational CM points O of discriminant-67.
There is a unique rational CM point R of discriminant-163 such that
R-O belongs to2Jac(C)(Q). Then

    S = O + 2(R-O)

is a rational non-CM point of C. This is an intrinsic specification of a
particular point after choosing O, not an assertion that an unspecified
member of a finite set works. Every one of the four choices of O works.
Coordinates on the displayed426b1 equation remain uncomputed.

## Atkin–Lehner action and the two CM orbits

The [complete CM calculation](DET852_RATIONAL_CM_LOCUS_2026-09-15.md)
gives exactly four rational points of each discriminant-67 and-163.
The residual Atkin–Lehner group on C is a Klein four-group, with nonidentity
classes w2,w3,w6. The checker computes the fixed-point counts upstairs
for labels2,3,6,71,142,213. A residual w_m-fixed point lifts to a point
fixed by w_m or w_(426/m). The two fixed loci are disjoint: otherwise a
Klein four-group would fix a point on a smooth characteristic-zero curve,
whose finite stabilizer is cyclic. Dividing their total size by two gives
four fixed points for each of residual w2 and w3, and none for w6.
None of their fixed CM discriminants is-67 or-163. Consequently neither
w2 nor w3 has a rational fixed point, and the action on each four-point
rational CM locus is free, hence transitive. These nontrivial residual
actions also verify faithfulness of the quotient group.

Use O as elliptic origin. The fixed-point-free involution w6 is translation
by a nonzero rational point T of order two. The other involutions have
forms x->Q-x and x->Q+T-x. Since neither has a rational fixed point,
Q and Q+T are not in2C(Q).

Exact PARI descent on426b1 gives rank bounds[1,1]. Good-reduction orders
at5,7,11 bound rational torsion by two, and a nonzero rational2-torsion
point is displayed. Thus C(Q) is isomorphic to Z plus Z/2, and
C(Q)/2C(Q) has four elements. The class of T is nonzero; the preceding
nondivisibility conditions show that the classes of T and Q are independent.

Modulo2C(Q), the four orbit points

    x, x+T, Q-x, Q+T-x

therefore occupy all four cosets. Both rational CM orbits have this
property. In the zero coset the CM(-67) orbit has just O, while the
CM(-163) orbit has a unique point R. This proves the stated existence
and uniqueness without needing numerical CM coordinates.

## Doubling leaves the entire rational CM locus

The point R is nonzero, since its CM order differs from O's. It belongs to
2C(Q), whose torsion is trivial because C(Q)=Z plus Z/2. Thus R is
non-torsion. Its double is neither O nor R. But O and R are the only CM
points in2C(Q), and2R remains in that subgroup. Therefore2R, equivalently
O+2(R-O) in origin-free notation, is non-CM.

This closes the particular-period question by an intrinsic algebraic
definition. It does not provide its coordinates, a K3 equation or a
rational divisor basis. The remaining arithmetic-admission step is the
primitive NS complement and actual rational divisor descent. No frame
search is claimed or started here.

```sh
sage -python research/elkies-k3/scripts/certify_det852_intrinsic_noncm_point.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det852-cm-gate-v1/intrinsic-noncm.json)
checks the residual fixed loci, exact rank bounds, torsion, and the finite
coset action. The written genus-one automorphism and doubling arguments
are proof inputs, together with the previous full-curve and complete-CM
theorems and Ogg's formula. Independent implementation, formal verification,
external review and literature novelty are unclaimed.
