# A 48 term collision frame and two exact optimality results

There is a rational polynomial map K in three variables with

\[
K(0)=K(e_1)=0,\qquad JK(0)=I,\qquad\det JK=1,
\]

coordinate degrees (4,6,7) and nonlinear occurrence counts (6,19,23).
Thus the [ten-term exclusion](../sparse-nine/PROOF.md) and this map bound
the normalized raw-degree-at-most-seven support minimum by [10,48].

The map is the known foundational counterexample in a different affine
frame. It has invariant degree profile (7,6,4), so it does not give a
counterexample in a profile strictly below that flag.

## Explicit formula

Put

\[
Y=y-2,\qquad Z=z+14y+24x-16,\qquad U=1+xY,
\]

and set

\[
\begin{aligned}
K_1&=x-\frac32x^2Y-\frac12x^3Z,\\
K_2&=y+3xU^2Z+3xY^2(4+3xY),\\
K_3&=U^3Z+Y^2U(4+3xY)+2K_2-64K_1.
\end{aligned}
\]

The [coefficient witness](normalized-48-witness.json) expands this formula
exactly. An independent standard-library rational implementation checks all
48 nonlinear occurrences, the three degrees, the full determinant, the
linear part and both colliding points.

If H is the [49-term map](../sparse-nine/NORMALIZED_WITNESS.md), let
S=I+14 E32. Then K(X)=S^(-1) H(SX). This conjugation fixes zero and e1,
preserves the identity derivative at zero and does not raise degree.

In terms of the foundational map F, the original collision points remain

\[
a=(0,-2,-16),\qquad b=(1,-2,8).
\]

The new source completion is

\[
L=\begin{pmatrix}1&0&0\\0&1&0\\24&14&1\end{pmatrix},
\qquad
K(X)=(JF(a)L)^{-1}\bigl(F(a+LX)-F(a)\bigr).
\]

The extra transverse coefficient 14 lies outside the earlier frame class
with fixed transverse columns ey,ez, so this does not contradict its sharp
minimum 49.

## Family A: varying the boundary pair and one transverse shear

Use the complete boundary collision family

\[
a(w)=(0,4w+6,-63w^2-189w-142),\qquad
b(w)=(1,w,2-3w)
\]

and its normalized maps H_w from the
[previous frame theorem](../sparse-nine/NORMALIZED_WITNESS.md).
For every w,s define

\[
K_{w,s}(X)=S_s^{-1}H_w(S_sX),\qquad S_s=I+sE_{32}.
\]

This permits the transverse columns ey+s*ez and ez. Nonzero independent
column scalings preserve support, so the family also covers every invertible
lower triangular transverse completion after the stated torus normalization.

The complete two-parameter support minimum is **48**, attained at
w=-2,s=14. The generic coefficient support has 70 occurrences. Away from
the finite critical projection, the smallest root stratum has 65.

Both parameters range over the algebraic closure; this is not a rational
sample or a search in a coefficient box.

## Family B: every transverse basis at the fixed collision pair

Keep the collision pair a,b above. Permit arbitrary independent columns
in the original y,z-plane, represented by a matrix B in GL2. Complete the
collision direction by those two columns and apply the same derivative
normalization.

The support minimum over this entire class is also **48**.

To cover all bases, nonzero diagonal column scalings and a column swap can
be discarded because they induce support-preserving simultaneous source and
target changes. After a swap, write

\[
B=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix},
\quad\alpha\ne0,\quad\Delta=\alpha\delta-\beta\gamma\ne0.
\]

Then

\[
B\,\operatorname{diag}(\alpha^{-1},\alpha/\Delta)
=
\begin{pmatrix}1&r\\s&1+sr\end{pmatrix}
=
\begin{pmatrix}1&0\\s&1\end{pmatrix}
\begin{pmatrix}1&r\\0&1\end{pmatrix},
\]

where s=gamma/alpha and r=alpha*beta/Delta. This two-parameter product is
always invertible. It therefore covers all GL2 transverse bases up to the
discarded support-preserving changes.

The resulting family has generic support 111. Away from its critical
projection, the smallest root stratum has 52. The point s=14,r=0 gives
the displayed 48-term map.

This class does not permit transverse columns with nonzero x-component.

## Complete algebraic support stratification

For each family, expand every nonlinear source monomial coefficient as a
polynomial c_j(w,t) in its two parameters. The parameter names in machine
files are w,t in both cases; in Family B they denote s,r.

A complete pairwise coprime factor basis for these coefficients describes
all possible cancellations. There are 41 factor curves for Family A and
62 for Family B.

Separate vertical factors from moving factors of positive t-degree.
The finite critical w-projection includes:

1. every root of a vertical factor;
2. every root of a moving factor's leading t-coefficient;
3. every root of the resultant in t of two distinct moving factors.

Outside this projection, all moving factors retain positive degree and are
pairwise coprime after specialization. A root of one factor kills exactly
the coefficients containing that factor. This proves the stated generic
root-stratum minima.

The critical projection has 469 irreducible factors for Family A and 824
for Family B. For each factor q(w), work exactly in Q[w]/(q) and form a
coprime root basis in the remaining parameter t. Every coefficient is either
zero or a nonzero constant times a product of root-basis factors. Pairwise
coprimality makes the vanishing pattern constant on the roots of each basis
factor. Counting the remaining coefficients gives every specialized minimum.

All 1,293 specializations have minimum at least 48. The explicit rational
point attains 48, proving sharpness in both families. No numerical root
approximation or irreducibility assumption about the root basis in t is
used.

## Independent evidence and corrected exploratory code

The primary calculations use exact SymPy coefficient arithmetic and formal
number-field residues. The independent Singular replay verifies:

- irreducibility of every critical projection factor;
- complete coefficient factorization and factor coprimality;
- every resultant and every projection factor product;
- all specialized coefficient polynomials, complete root factorization and
  coprime root bases;
- the minimum in every one of the 1,293 exact specializations.

The [standard-library checker](verify.py) separately reconstructs both
parameter families by rational polynomial multiplication, checks every
coefficient input and projection factor product, and verifies the actual
48-term determinant and collision.

An earlier exploratory specialization reused one formal algebraic root
symbol across different minimal polynomials. It was rejected after a
known rational parameter contradicted its output. Its preserved ledgers are
marked invalid; the theorem uses only the v2 ledgers with distinct
generators, exact modulus guards and the independent Singular replay.
The [failure record](INVALIDATED_SPECIALIZATION.md) gives the precise
boundary. The predecessor sparse proofs and 49-term frame certificate
did not use that construction.

The bounded beam search is a discovery aid, not evidence that the whole
affine orbit has been searched.

## Reproduction and remaining scope

From the repository root:

~~~sh
.venv/bin/python research/jacobian/collision-frame-search/two_parameter.py
.venv/bin/python research/jacobian/collision-frame-search/specialize.py
.venv/bin/python research/jacobian/collision-frame-search/full_transverse.py
.venv/bin/python research/jacobian/collision-frame-search/verify_singular.py
.venv/bin/python research/jacobian/collision-frame-search/verify.py
.venv/bin/python research/jacobian/collision-frame-search/check_field_generators.py
~~~

The intermediate stages checkpoint exact completed cases and reject source
drift. Existing completed cases are reused with their source hashes; choose
a separate copy of the output directory for an entirely fresh replay.

Unrestricted transverse x-components, jointly varying the boundary pair
and both GL2 parameters, interior collision pairs, other Keller families,
and the true normalized support minimum remain open. No ordinary degree
improvement, global affine-orbit optimum, formal proof-assistant certificate
or external review is claimed.

