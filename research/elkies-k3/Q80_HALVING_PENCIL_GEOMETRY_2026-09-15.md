# A rational bisection does not supply a second rational pencil of branch halves

**Written geometric deduction with checked finite premises.** For the primitive
trace parity65538 of the retained Q80 rational bisection
[alternate-orbit-11ae6](Q80_ONE_COVER_HALVING_CONTROL_2026-09-14.md), the genus-nine
halving curve D_T has a unique degree-four map to a rational curve, up to a
change of target coordinate. It is the original parameter map D_T -> P1_t.
There are no maps of smaller degree to a rational curve.

Therefore a closed point P of degree2,3 or4 on D_T with

```
Q(P)=Q(t(P))
```

cannot come from the rational fibres of any pencil of that same degree on
D_T. This is a condition on equality of residue fields, not on whether the
number field has proper subfields. It is not a finiteness or emptiness theorem
for these points. No new two-gain cover is constructed.

## 1. The ruled quotient of a smooth rational bisection

The following argument applies to any24-I1 elliptic K3 with a smooth rational
bisection B whose trace T is geometrically primitive modulo2MW. Put
eta(P)=T-P on the generic elliptic fibre. Translation and inversion extend
to automorphisms of the minimal K3, so eta is an involution preserving B.
Its quotient Y has the induced ruling over the original t-line.

The quotient is smooth: eta acts by -1 on the holomorphic two-form, so has
no isolated fixed points. Every smooth elliptic fibre has rational quotient.
At an I1 fibre the involution swaps the two branches at the node; the local
model xy=t with x and y exchanged has smooth quotient with coordinates
x+y,t. Thus all quotient fibres are smooth rational curves. The image of B
is a section C of the ruling, and B is its full inverse image. Projection
and adjunction give

```
B^2=-2,   2*C^2=B^2,   C^2=-1.
```

Hence Y is the Hirzebruch surface F1 and C is its negative section. Over Q,
C maps isomorphically to the original rational base, so this is the split
ruled surface. Write F for a ruling fibre. The fixed curve D_T={2P=T} maps
isomorphically to the branch curve R on Y. The double-cover canonical formula
and K3 trivial canonical bundle give

```
C^2=-1, C.F=1, F^2=0,
K_Y=-2C-3F,   R=-2K_Y=4C+6F.
```

The [halving monodromy proof](JOINT_HALVING_EXCEPTIONAL_BRANCH_LOCUS_2026-09-14.md)
gives connected D_T of genus9, full S4 monodromy and a primitive degree-four
map f to P1_t. The same genus follows from adjunction on F1. Notice R.C=2.

For the actual control, the prior certificate verifies the point equation,
nontrivial trace parity, smooth branch fibres and rational covering conic.
The new checker verifies gcd(x1,y1)=1 in the retained expression
x=x0+x1*w, y=y0+y1*w. Thus conjugate branches never collide away from the
branch divisor, and at a branch at least one of their first derivatives in
w is nonzero. At infinity deg(x1)=3 and w/t is nonzero, so the two points
have distinct minimal-chart abscissas. The image is therefore an embedded
smooth rational bisection; smooth normalization alone would not suffice.

## 2. Uniqueness of the degree-four pencil

Work over an algebraic closure. Suppose g:D_T -> P1 is another map of degree
at most4. If f and g do not generate the function field of D_T, primitivity
of f forces g to factor through f. Its degree then forces it to be f up to
a target automorphism.

Otherwise the [Castelnuovo–Severi inequality, Theorem14](https://link.springer.com/article/10.1007/s40993-024-00543-4)
gives genus at most3*(deg(g)-1). Degree at most3 contradicts genus9. In degree4,
the image in P1 times P1 has bidegree(4,4) and arithmetic genus9. Equality of
arithmetic and normalization genus makes that image smooth. Adjunction yields

```
K_D = 2*f + 2*g
```

as divisor classes of fibres restricted to D=D_T.

But K_D-2f is the restriction of 2C+F on F1. All sections of 2C+F vanish on C:
its intersection with C is -1, and the residual system |C+F| is basepoint-free.
The restriction has exactly those sections. Indeed the restriction sequence
has kernel O_Y(K_Y-2F), whose H0 and H1 vanish; the latter follows by Serre
duality from H1(Y,2F)=0. Thus |K_D-2f| has the nonempty fixed divisor C.R,
of degree2, and dimension of its space of sections3.

On the proposed (4,4) model this would be |2g|, a basepoint-free system,
since it contains the pullbacks of all quadratic polynomials on P1. This
contradiction proves uniqueness. The argument includes tangency of C and R:
the fixed length-two divisor need not consist of distinct points.

## 3. What this changes for branch construction

A usable branch value b requires an actual half over Q(b), so a corresponding
point P of D_T must satisfy Q(P)=Q(b). A rational fibre of f instead supplies
points over extensions of Q while b=t(P) stays rational. Even when such a
point has degree4, its residue field is too large for this purpose.

Uniqueness proves that changing the rational pencil on this halving curve
cannot repair the mismatch: there is no second degree-four pencil. For points
of degree2 or3 there is no rational pencil of that degree at all. In the usual
terminology these usable points are P1-isolated. This terminology does not
assert isolation in the Jacobian or finiteness of their set.

Possible points not generated by rational pencils, and maps to positive-genus
curves, remain outside the exclusion. In particular no simplicity, rank or
elliptic-factor result for Jac(D_T) is proved here. Even a usable branch point
would still need a second independent halving condition or a compatible
nonzero branch class, and global trace-torsor solutions on the same cover.
The [requested MW17 two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open. The first retained rational cover was already proved to have
exact rank18; this result concerns possible new branch values for its trace
parity, not a second rank bound for that fixed cover.

## Evidence

[Checker](scripts/verify_q80_halving_pencil_geometry.py) and
[certificate](../artifacts/generated-results/elkies-k3-q80-halving-pencil-geometry-v1/result.json):

```
.venv/bin/python research/elkies-k3/scripts/verify_q80_halving_pencil_geometry.py
```

The short replay checks the bisection's embedding premises and F1 intersection
arithmetic. The quotient geometry, cohomology calculation and uniqueness proof
are written mathematics; no independent whole-proof replay, formal verification
or external review is claimed. No point search or research census was run.
