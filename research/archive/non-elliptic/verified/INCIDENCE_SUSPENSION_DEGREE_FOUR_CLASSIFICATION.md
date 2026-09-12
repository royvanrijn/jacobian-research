# Incidence suspensions through horizontal degree four

This note classifies the natural bounded search suggested by the weighted and
quadratic marked-line constructions.  It has two parts:

1. an exact incidence-suspension criterion for a fixed reciprocal chart;
2. a complete classification, through `deg X<=4`, of root-preserving,
   `P`-fibration-preserving affine recharts of the quadratic reciprocal
   chart.

The result is a scoped classification, not a classification of every
birational chart of affine three-space.  Inside the stated class, the cubic
and quartic horizontal coordinates exist as reciprocal rational charts but
fail polynomial pullback for an unavoidable reason.  Only the quadratic
horizontal coordinate survives.

Work over a characteristic-zero field `k`.

## 1. Exact incidence-suspension criterion

Let

\[
 \rho:\mathbb A^3_{x,y,z}\dashrightarrow
 \mathbb A^3_{P,S,Q}
\]

be a birational chart with

\[
 \det D\rho=D^{-1},\qquad
 D=D_0(P,S)-\kappa(P,S)Q.                              \tag{1}
\]

Fix `lambda in k^times`.  For functions `X,Y,beta` depending only on
`(P,S)`, put

\[
 B=Q+\beta(P,S),\qquad C=Y(P,S)-B X(P,S).              \tag{2}
\]

Then

\[
 \det\frac{\partial(P,B,C)}{\partial(P,S,Q)}
 =-\bigl(Y_S-BX_S\bigr).                               \tag{3}
\]

Consequently

\[
 \det\frac{\partial(P,B,C)}{\partial(x,y,z)}
 =-\lambda
\]

if and only if

\[
 \boxed{
 X_S=\lambda\kappa,\qquad
 \beta=\frac{Y_S-\lambda D_0}{X_S}.
 }                                                       \tag{4}
\]

The second equality includes the regularity requirement on `beta`.  Thus,
for a fixed chart, the horizontal curve coordinate `X` is forced up to a
function of `P`; adding such a function merely shears `C` by a multiple of
`B`.  The remaining construction problem is linear in `Y`, followed by the
nonlinear requirement that (2) pull back to polynomial source coordinates.

This proves both necessity and sufficiency within the
coordinate-preserving, affine-linear-in-`Q` incidence ansatz.

## 2. The quadratic reciprocal chart

Put

\[
 t=1+xy,\qquad q_z=t^2,
\]

and use the chart

\[
 P=tq,\qquad S=\frac{x}{t},\qquad Q=y+xq.              \tag{5}
\]

The choice of the `z`-free part of `q` does not affect the following chart
identities.  One has

\[
 D=\frac1t=1-SQ+PS^2,\qquad
 \det\frac{\partial(P,S,Q)}{\partial(x,y,z)}=t=D^{-1}. \tag{6}
\]

Here `kappa=S`, so (4) forces

\[
 X=\frac{\lambda}{2}S^2+\text{function of }P.          \tag{7}
\]

After normalization `lambda=2`, this is exactly the quadratic gauge
`X=S^2`.  In particular, changing `Y` alone can never produce a cubic or
quartic tilt in this chart.

## 3. All root-preserving `P`-fibre-affine recharts

We now allow the smallest birational change capable of altering `kappa`.
Keep the marked root `S` and set

\[
 \widetilde P=A(S)P+E(S),\qquad
 \widetilde Q=A(S)^{-1}Q+H(\widetilde P,S),             \tag{8}
\]

where `A in k[S]` is nonzero and `E,H` are polynomial.  The two variable
change `(P,Q)->(tilde P,tilde Q)` has determinant one at fixed `S`.
Conversely, every root-preserving birational change which preserves the
`P`-fibration, is affine in `P` and `Q`, and has determinant one on each
`S`-fibre is of the form (8), up to constant source and target scalings.

Substitution in (6) gives

\[
 D=
 1-SA(S)\widetilde Q
 +SA(S)H(\widetilde P,S)
 +\frac{(\widetilde P-E(S))S^2}{A(S)}.                 \tag{9}
\]

Thus

\[
 \widetilde\kappa=SA(S),\qquad
 \widetilde X_S=\lambda SA(S).                         \tag{10}
\]

The coefficient of `tilde P` in `D_0` is

\[
 \frac{S^2}{A(S)}+SA(S)H_1(S),                         \tag{11}
\]

where `H_1` is the coefficient of `tilde P` in `H`.  The second summand is
polynomial.  Therefore polynomiality of `D_0` forces

\[
 \boxed{A(S)\mid S^2.}                                 \tag{12}
\]

Conversely, if (12) holds, taking `H=E=0` makes `D_0` polynomial.  Over a
field, the nonzero polynomial divisors of `S^2`, up to a scalar, are

\[
 A=1,\qquad A=S,\qquad A=S^2.                          \tag{13}
\]

They give, after normalizing `lambda`,

\[
\begin{array}{c|c|c|c}
A(S)&\widetilde\kappa&X(S)&D_0\text{ for }E=H=0\\ \hline
1&S&S^2&1+\widetilde P S^2\\
S&S^2&S^3&1+\widetilde P S\\
S^2&S^3&S^4&1+\widetilde P.
\end{array}                                             \tag{14}
\]

This is the complete reciprocal-chart list in the bounded root-preserving,
`P`-fibration-preserving affine search through horizontal degree four.

## 4. Why the cubic and quartic charts do not polynomialize

Formula (8) retains the marked-variable coefficient

\[
 \widetilde Q=A(S)^{-1}Q+\text{function of }
 (\widetilde P,S).                                     \tag{15}
\]

The incidence criterion always gives

\[
 \widetilde B=\widetilde Q+\widetilde\beta
 (\widetilde P,S).                                     \tag{16}
\]

Neither `H` nor `beta` depends on `Q`, so the coefficient `A^{-1}` in front
of the algebraically independent chart variable `Q` cannot cancel.
If `A` is nonconstant, its zero divisor meets the original affine source
chart, and (16) is not regular there.  Equivalently, for the two new cases,

\[
\begin{array}{c|c|c}
X&\widetilde P&\text{unavoidable term in }\widetilde B\\ \hline
S^3&SP&Q/S\\
S^4&S^2P&Q/S^2.
\end{array}                                             \tag{17}
\]

On the source, `S=x/t` and `Q=y+xq`.  At the generic point of `x=0`,
`Q` is a unit while `S` vanishes simply, so these are genuine poles.

No polynomial target automorphism can repair this.  If a polynomial
automorphism made the rational pair `(tilde B,tilde C)` polynomial, its
polynomial inverse would make `tilde B` and `tilde C` polynomial as well,
contradicting (16).

The divisibility problem for `Y` is not the obstruction.  For `A=S` and
`A=S^2`, equation (4) imposes only a finite Taylor congruence on `Y_S`, which
is solvable in characteristic zero.  Failure occurs one step later, at
polynomial pullback.

We have proved:

### Degree-four fibre-affine classification theorem

Among all root-preserving, `P`-fibration-preserving affine reciprocal
recharts (8) of the quadratic source chart with polynomial controlled
divisor and `deg X<=4`, the only chart producing a polynomial marked-line
Keller suspension of `A^3` is the quadratic chart `X=S^2`.

The cubic and quartic cases are valid rational constant-Jacobian
suspensions on localized charts, but each carries an additional vertical
pole divisor and does not extend to a polynomial map of affine three-space.

## 5. Relation to the weighted linear incidence

The weighted tangent construction with `X=S` is not a rechart of (5) inside
the class (8).  It is a polynomial vertical suspension of a different plane
incidence, rather than a reciprocal birational suspension.  Thus the two
surviving polynomial mechanisms remain

\[
\begin{array}{c|c|c}
\text{mechanism}&X&\text{suspension type}\\ \hline
\text{weighted tangent}&S&\text{polynomial vertical}\\
\text{quadratic gauge}&S^2&\text{reciprocal birational}.
\end{array}                                             \tag{18}
\]

The bounded theorem explains why the naive `S^3` and `S^4` attacks do not
produce further stable components: their reciprocal charts necessarily
introduce a second pole divisor before seed engineering begins.

## 6. What remains open

The theorem does not exclude:

1. a reciprocal chart not birationally `P`-fibre-affine over the known
   `(P,S,Q)` chart;
2. a chart which changes the marked-root coordinate nontrivially;
3. a multi-boundary ledger in which an additional chart factor cancels the
   pole of `A^{-1}`;
4. an incidence with more than one marked-line coefficient.

A global classification of incidence suspensions must either control these
possibilities intrinsically or exhibit a genuinely new chart.  The bounded
search shows that merely reweighting the known reciprocal chart cannot do
so.

## Exact regression

Run

```bash
.venv/bin/python scripts/verify_incidence_suspension_degree_four.py
```

The checker verifies the universal incidence determinant, the reciprocal
chart, all three divisor cases in (14), the cubic and quartic rational
charts, and their unavoidable pullback poles.
