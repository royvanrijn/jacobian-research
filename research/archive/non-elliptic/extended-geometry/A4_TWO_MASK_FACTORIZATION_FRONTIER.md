# The two-mask \(A_4\) factorization frontier

## 1. How the construction continues

The [pure-target lift](A4_PURE_TARGET_LEDGER_LIFT.md) gives

\[
\widehat\Phi(U,V,W,z_1)
=\left(P,Q,R,S\right)
=\left(WN_1,WN_2,WH,\frac{WL}{4}z_1\right)
\]

with

\[
\det D\widehat\Phi=\mathcal B(P,Q,R).
\]

To obtain an ordinary Keller map without losing the \(A_4\) extension, the
cleanest sufficient mechanism is a polynomial factorization

\[
\widehat\Phi\times\operatorname{id}
=\beta\circ F                                      \tag{1.1}
\]

after adjoining a second mask variable, where

\[
\det D\beta=\mathcal B,\qquad \det DF=1.
\]

If \(\beta\) is generically birational and \(F\) is polynomial, then (1.1)
automatically gives:

- constant Jacobian for \(F\);
- generic degree four;
- the same natural \(A_4\) normal closure;
- one-to-one reconstruction of the auxiliary variables.

The problem is therefore polynomial divisibility through a birational
target blowdown, not a free-form search through all fivefold maps.

## 2. A two-mask target blowdown

Put

\[
X=P-Q
\]

and

\[
\mathcal U
=27R^2-X^2-3(X-3R)Q.                               \tag{2.1}
\]

The homogeneous target cubic has the exact Bezout presentation

\[
\boxed{
\mathcal B(P,Q,R)=27R^3-X\mathcal U.
}                                                    \tag{2.2}
\]

It also has the opposite-puncture presentation

\[
\mathcal B(P,Q,R)
=(X-3R)\mathcal V-27R^3,                            \tag{2.3}
\]

where

\[
\mathcal V=3XQ+X^2+3XR-18R^2.
\]

These are the homogeneous versions of the two incidence identities for
the punctures \(r=0\) and \(r=1\).

Equation (2.2) is a \(2\times2\) determinantal representation:

\[
\mathcal B
=\det
\begin{pmatrix}
27R^3&X\\
\mathcal U&1
\end{pmatrix}.                                      \tag{2.4}
\]

Consequently

\[
\boxed{
\begin{aligned}
\beta_2:\mathbb A^5_{P,Q,R,A,B}
&\longrightarrow\mathbb A^5,\\
(P,Q,R,A,B)
&\longmapsto
\left(
P,Q,R,\,
27R^3A+XB,\,
\mathcal U A+B
\right)
\end{aligned}
}                                                    \tag{2.5}
\]

has

\[
\boxed{\det D\beta_2=\mathcal B(P,Q,R).}             \tag{2.6}
\]

It is generically birational.  Away from \(\mathcal B=0\), its inverse
mask coordinates are

\[
\begin{aligned}
A&=\frac{S-XT}{\mathcal B},\\
B&=\frac{-\mathcal U S+27R^3T}{\mathcal B}.
\end{aligned}                                       \tag{2.7}
\]

This is the first target modification adapted simultaneously to both
finite puncture characters of \(L\).  Unlike the diagonal map
\((S,T)\mapsto(\mathcal B S,T)\), its two inverse numerators can cancel
different boundary residues.

## 3. Every one-primitive triangular rechart fails

Before using both masks, consider any target automorphism which fixes the
primitive coordinate up to a nonzero constant and changes the three base
coordinates only by multiples of it:

\[
\begin{aligned}
(P,Q,R,S)\longmapsto
\bigl(
P+Sh_1,\ Q+Sh_2,\ R+Sh_3,\ cS
\bigr),
\qquad c\ne0,                                       \tag{3.1}
\end{aligned}
\]

where the \(h_i\) are arbitrary polynomials for which (3.1) is an
automorphism.

After pullback by \(\widehat\Phi\), the prospective denominator is

\[
E=\mathcal B(P+Sh_1,Q+Sh_2,R+Sh_3).
\]

Since \(S=(WL/4)z_1\),

\[
E\bmod z_1
=\mathcal B(P,Q,R)
=W^3K^3L^2.                                        \tag{3.2}
\]

If the inverse primitive \(cS/E\) were polynomial, then \(E\mid cS\).
The polynomial \(E\) is coprime to \(z_1\), so \(E\mid WL\).  Reducing
modulo \(z_1\) would then give

\[
W^3K^3L^2\mid WL,
\]

which is impossible.  Therefore no triangular shear retaining one
primitive coordinate can factor the lift.  This is an all-degree
obstruction, not a bounded coefficient search.

## 4. Exhaustive coordinate-permutation screen for \(\beta_2\)

Adjoin a second source coordinate \(z_2\) and write the five outputs

\[
\mathcal Y=
\left(
WN_1,\ WN_2,\ WH,\ \frac{WL}{4}z_1,\ z_2
\right).                                            \tag{4.1}
\]

There are \(5!=120\) ways to assign these outputs, in order, to

\[
(P,Q,R,S,T)
\]

in (2.7).  For each assignment, form

\[
\begin{aligned}
D&=\mathcal B(P,Q,R),\\
n_1&=S-(P-Q)T,\\
n_2&=-\mathcal U(P,Q,R)S+27R^3T.                    \tag{4.2}
\end{aligned}
\]

Polynomial factorization through \(\beta_2\) requires simultaneously

\[
D\mid n_1,\qquad D\mid n_2.                         \tag{4.3}
\]

Exact multivariate division over \(\mathbb Q\) gives:

\[
\boxed{
\text{None of the \(120\) coordinate assignments satisfies (4.3).}
}                                                    \tag{4.4}
\]

This is exhaustive for coordinate permutations of (4.1).  It is not a
claim about arbitrary linear changes or nonlinear target shears.

The result has a useful interpretation.  Merely deciding that a different
existing output should serve as a mask does not align the two numerator
residues in (2.7).  The masks must alter the base coordinates before the
two-mask blowdown is inverted.

## 5. Every zero-section-preserving conjugation fails

The apparent next move is to conjugate \(\beta_2\) by an automorphism which
preserves the mask ideal

\[
I=(S,T).
\]

In fact the whole class is impossible.  Modulo \(I^2\), such an
automorphism acts on the normal module \(I/I^2\) by a matrix

\[
A(P,Q,R)\in\operatorname{GL}_2
\bigl(\mathbb Q[P,Q,R]\bigr),
\qquad \det A\in\mathbb Q^*.                         \tag{5.1}
\]

The inverse of \(\beta_2\) acts on the same normal module by

\[
C^{-1}
=\frac{\operatorname{adj}C}{\mathcal B},
\qquad
C=
\begin{pmatrix}
27R^3&P-Q\\
\mathcal U&1
\end{pmatrix}.                                      \tag{5.2}
\]

If the conjugated inverse were polynomial, then \(C^{-1}A\) would be a
polynomial matrix.  But

\[
\det(C^{-1}A)
=\frac{\det A}{\mathcal B},                          \tag{5.3}
\]

which is not polynomial.  Therefore:

\[
\boxed{
\text{No polynomial target automorphism preserving \(I=(S,T)\)
can make the two-mask inverse polynomial.}
}                                                    \tag{5.4}
\]

This closes every affine-linear or nonlinear base shear whose corrections
vanish on the mask zero section.  The next rechart must move that zero
section into the exceptional divisor \(\mathcal B=0\).

## 6. First singular-line incidence charts

The target cubic has the singular line

\[
\Sigma=\{P=Q,\ R=0\}\subset V(\mathcal B).
\]

This gives a linear way to move the mask zero section into the exceptional
divisor.  For example,

\[
(p,q,r,s,t)=(P,\ P+S,\ T,\ Q,\ R)                   \tag{6.1}
\]

is a linear automorphism of five-space and sends \(S=T=0\) to \(\Sigma\)
in the new base coordinates.

There are six coordinate versions of (6.1): choose one of \(P,Q,R\) as
the repeated singular-line parameter and order the other two as the mask
outputs.  Exact division of the two inverse numerators (2.7) gives

\[
\boxed{\text{All six singular-line incidence charts fail.}}   \tag{6.2}
\]

The failure of (6.1) itself is already visible in the first numerator:
after pullback it is not divisible by the new
\(\mathcal B(p,q,r)\).  Thus moving the zero section into the singular
line supplies the required vanishing order, but not the required residue
alignment.

This six-chart result is finite.  It does not exclude general linear
embeddings of the zero section into \(\Sigma\), or embeddings into the
smooth part of \(V(\mathcal B)\).

## 7. Normalized-boundary assembly audit

The normalization of the target cubic supplies the next, nonlinear
incidence model.  With parameters \((\lambda,r)\), put

\[
\begin{aligned}
R_0&=\lambda r(r-1),\\
Q_0&=\lambda(-r^3+3r-1),\\
P_0&=\lambda(2r^3-3r^2+3r-1).
\end{aligned}                                       \tag{7.1}
\]

Then

\[
\mathcal B(P_0,Q_0,R_0)=0.                           \tag{7.2}
\]

The [assembly audit](A4_NORMALIZED_BOUNDARY_ASSEMBLY_AUDIT.md) constructs an
explicit unimodular ambient completion of (7.1) and performs the full
factorization test.  Both inverse numerators fail divisibility, and the
rational Jacobian is a nonconstant quotient of the old and new
\(\mathcal B\)-pullbacks.

More generally, an automorphic rechart can yield constant Jacobian only if
it preserves \(\mathcal B\) up to a scalar.  It therefore cannot move the
generic old mask-zero section into \(\mathcal B=0\).  The surviving problem
is a nonautomorphic log-crepant incidence map \(\alpha\) satisfying

\[
\mathcal B(\operatorname{pr}_{1,2,3}\alpha)
=u\,\mathcal B(P,Q,R)\det D\alpha
\]

together with the two adjugate divisibilities.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_two_mask_factorization.py
```

The checker verifies both Bezout presentations, the determinant and
birational inverse of \(\beta_2\), the one-primitive divisibility gate, all
\(120\) exact coordinate-permutation divisions, the six singular-line
incidence charts, and the normalized cubic parametrization (7.1).  The
separate assembly checker verifies the unimodular completion and its
failure.
