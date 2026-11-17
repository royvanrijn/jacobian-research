# Unrestricted binary GVC by Hall-envelope separation

## 1. Statement and status

Let $k$ be a characteristic-zero field, let

\[
 \Lambda=\lambda(\partial_x,\partial_y)
\]

be a constant-coefficient differential operator in two variables, and let
$P,Q\in k[x,y]$.  This note proves the following statement.

> **Theorem 1.1 (binary Generalized Vanishing Conjecture).**  If
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
>  \tag{1.1}
> \]
> then
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
>  \tag{1.2}
> \]

The new step is a finite Newton-envelope argument.  It combines three
previously proved ingredients:

1. binary Hall localization at the lowest homogeneous symbol;
2. shifted-ray endpoint rigidity for unequal positive weights; and
3. terminality at an unequal common threshold.

No Hall carry promotion, factorial packet classification, bounded radical
calculation, or degree induction is used.  Those calculations remain useful
as diagnostics, but they are not needed for Theorem 1.1.

The proof is internal and has not been externally reviewed.

The [finite-certificate corollary](BINARY_GVC_FINITE_CERTIFICATE.md) extracts
a strict separator over the original field, an exact rational decision
procedure, and the degree-only cutoff
\(m>(\deg P+\deg\lambda)\deg Q\) from this proof.

## 2. The three inputs

If \(\lambda=0\), the conclusion is immediate.  If the constant term of
\(\lambda\) is nonzero, the highest ordinary-degree part of (1.1) is that
constant to the power \(m\), times \(P_d^m\); hence \(P=0\).  We may
therefore assume that \(\Lambda\) has positive lowest order.

Write $r$ for the lowest positive order of $\Lambda$, and put
$d=\deg P$.  The case $d\leq r$ is the general degree-cutoff theorem in
[`SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md`](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md),
so assume $d>r$.

Let $\lambda_r$ and $P_d$ be the lowest operator symbol and the top
homogeneous part of $P$.  Taking the highest ordinary-degree part of
(1.1) gives

\[
 \lambda_r(\partial)^m(P_d^m)=0\qquad(m\geq1).
 \tag{2.1}
\]

The uniform binary Hall-localization theorem gives a root of
$\lambda_r$, of full multiplicity $e$, and coordinates in which

\[
 \lambda_r=X^e C_{r-e}(X,Y),\qquad C_{r-e}(0,1)\ne0,
 \tag{2.2}
\]

and

\[
 P_d=y^{d-e+1}D_{e-1}(x,y).
 \tag{2.3}
\]

Thus the least $X$-exponent in $\lambda_r$ is exactly $e$, while the
largest $x$-exponent in $P_d$ is some $t\leq e-1$.

The second input is the shifted-ray endpoint theorem.  In the notation used
below it says the following.  Let $A$ and $B$ be nonzero polynomials
homogeneous for one positive weight $w$, of weights $W_A$ and
$W_B=W_A+w(\delta)$, where $\delta\in\mathbb Q_{\geq0}^2$.  If

\[
 [x^{m\delta_x}y^{m\delta_y}]
 A(\partial)^m(B^m)=0
 \tag{2.4}
\]

at every integral scale after clearing the denominator of $\delta$, then

\[
 \operatorname {Newt}(A)\cap
 \bigl(\operatorname {Newt}(B)-\delta\bigr)=\varnothing.
 \tag{2.5}
\]

This is Theorem 3.2 of
[`BINARY_GVC_UNIFORM_FACE_TERMINATION.md`](BINARY_GVC_UNIFORM_FACE_TERMINATION.md).
Its proof is the prime-dilated unique-minimum valuation argument, including
rational output rays.

The third input is the common-threshold corollary from the same note.  If a
positive unequal weight and a number $W$ satisfy

\[
 w(\operatorname {supp}\lambda)\geq W,
 \qquad
 w(\operatorname {supp}P)\leq W,
 \tag{2.6}
\]

then (1.1) implies (1.2).

## 3. Horizontal separation before a common threshold

The following observation is the missing link.

> **Lemma 3.1 (unequal-face no-overlap lemma).**  Let $w=(u,v)$ with
> $u,v>0$ and $u\ne v$.  Let $A,B\ne0$ be $w$-homogeneous of
> respective weights $W_A<W_B$, and suppose
> \[
>  A(\partial)^m(B^m)=0\qquad(m\geq1).
>  \tag{3.1}
> \]
> Then the projections of $\operatorname {Newt}(A)$ and
> $\operatorname {Newt}(B)$ to the $x$-axis are disjoint.

### Proof

Suppose that the two projected intervals overlap.  Choose a rational number
$c$ in their intersection.  Since both Newton polygons lie on affine
weight lines, there are unique points

\[
 \alpha=\left(c,\frac{W_A-uc}{v}\right)
 \in\operatorname {Newt}(A),
 \qquad
 \beta=\left(c,\frac{W_B-uc}{v}\right)
 \in\operatorname {Newt}(B).
 \tag{3.2}
\]

Put

\[
 \delta=\beta-\alpha
 =\left(0,\frac{W_B-W_A}{v}\right)
 \in\mathbb Q_{\geq0}^2.
 \tag{3.3}
\]

Equation (3.1) makes every coefficient on the rational ray $m\delta$
zero.  After clearing the denominators, the shifted-ray endpoint theorem
applies.  But
$\alpha\in\operatorname {Newt}(A)\cap
(\operatorname {Newt}(B)-\delta)$, contradicting (2.5).  Hence the
projected intervals are disjoint.  $\square$

The lemma is stronger than any bounded face-radical calculation needed in
the earlier degree-by-degree proofs.  It works with the actual nonzero
support on every coefficient stratum and at every unequal rational slope.

## 4. The two global envelopes

For $s\geq1$, use the weight

\[
 w_s(x)=w_s(X)=s,
 \qquad
 w_s(y)=w_s(Y)=1,
 \tag{4.1}
\]

and define the finite piecewise-linear envelopes

\[
 L(s)=\min_{\alpha\in\operatorname {supp}\lambda}w_s(\alpha),
 \qquad
 U(s)=\max_{\beta\in\operatorname {supp}P}w_s(\beta),
 \qquad
 \Delta(s)=U(s)-L(s).
 \tag{4.2}
\]

For sufficiently small rational $\epsilon>0$, (2.2)--(2.3) imply that
the unique lower operator endpoint at $s=1+\epsilon$ has $x$-exponent
$e$, while the unique upper polynomial endpoint has $x$-exponent $t$.
Consequently

\[
 e>t,
 \qquad
 \Delta(1+\epsilon)>0.
 \tag{4.3}
\]

At a rational slope $s>1$, let $A_s$ be the complete $L(s)$-face of
$\lambda$, and let $B_s$ be the complete $U(s)$-face of $P$.  The
part of largest $w_s$-weight in (1.1) is

\[
 A_s(\partial)^m(B_s^m)=0\qquad(m\geq1).
 \tag{4.4}
\]

Indeed, maximal output weight requires every operator selection to lie on
the minimum face and every polynomial selection to lie on the maximum
face.

Assume $\Delta(s)>0$.  Lemma 3.1 says that the two $x$-intervals of
$A_s$ and $B_s$ are disjoint.  Just above $s=1$, the operator interval
is strictly to the right of the polynomial interval by (4.3).  This order
cannot reverse while $\Delta>0$: at the first reversal slope, the old and
new endpoints belong to the complete equality faces, so their projected
intervals overlap, contrary to Lemma 3.1.  Therefore

\[
 \min\{\alpha_x:\alpha\in\operatorname {Newt}(A_s)\}
 >
 \max\{\beta_x:\beta\in\operatorname {Newt}(B_s)\}
 \tag{4.5}
\]

at every rational $s>1$ for which $\Delta(s)>0$.

Because both supports are finite, $\Delta$ has only finitely many linear
pieces.  It cannot remain positive for every $s$.  Otherwise, after the
last breakpoint, the minimizing operator monomial would have the globally
least $X$-exponent $a_{\min}$, the maximizing polynomial monomial would
have the globally largest $x$-exponent $i_{\max}$, and (4.5) would give
$a_{\min}>i_{\max}$.  On that final interval,

\[
 \Delta(s)=\text{constant}+s(i_{\max}-a_{\min}),
 \tag{4.6}
\]

whose slope is strictly negative.  It must cross zero, a contradiction.

Hence there is a first $s_*>1$ such that

\[
 \Delta(s_*)=0.
 \tag{4.7}
\]

All breakpoints and zeroes of the finite integer-exponent envelopes are
rational.  Multiplying $w_{s_*}$ by a common denominator gives distinct
positive integral weights.  With $W=L(s_*)=U(s_*)$, definitions (4.2)
give exactly

\[
 w_{s_*}(\operatorname {supp}\lambda)\geq W,
 \qquad
 w_{s_*}(\operatorname {supp}P)\leq W.
 \tag{4.8}
\]

The common-threshold corollary now proves (1.2).  This completes the proof
of Theorem 1.1 for $d>r$, and the degree-cutoff theorem supplies the
remaining case.  $\square$

## 5. Why the computed survivors have Ferrers radicals

The envelope proof also explains the degree-eight calculation which led to
Lemma 3.1.  At a slope-two wall, let the marked operator and polynomial
$x$-exponents differ by $g>0$.  Index tied operator corrections by their
leftward deficits $i$, and tied polynomial corrections by their rightward
advances $j$.  The two face intervals overlap exactly when

\[
 i+j\geq g.
 \tag{5.1}
\]

Lemma 3.1 says that no actual pure-zero support can contain such a pair.
Set-theoretically, the allowed coefficient supports are therefore the
independent sets of the Ferrers graph with edge set $i+j\geq g$.  Its
squarefree Stanley--Reisner ideal is

\[
 (a_i:i\geq g)+(b_j:j\geq g)
 +(a_i b_j:i,j<g, i+j\geq g).
 \tag{5.2}
\]

For the first nontrivial octic wall, $g=3$, exact moments through order six
give

\[
 (A,S,T,\ BP,\ BQ,\ CQ),
 \tag{5.3}
\]

whose three noncoordinate components are precisely the vertex covers of
the staircase $B\!-!P,B\!-!Q,C\!-!Q$.  At degree nine, the exact
characteristic-zero saturation checks also recover the gap-four staircase.
These computations are consequences and regressions of the horizontal
separation mechanism, not ingredients in its proof.

On an independent set, if $i_{\max}$ and $j_{\max}$ are the extreme
nonzero migrations, then

\[
 g'=g-i_{\max}-j_{\max}>0.
 \tag{5.4}
\]

This is the discrete gap descent seen in the septic and octic searches.
The envelope proof packages all such descents at once and avoids computing
their primary decompositions.

## 6. Consequences for the former frontier

The all-scale Hall/carry promotion problem was a sufficient route to
unrestricted binary GVC, but it is not necessary.  Its semigroup, character,
carry, and factorial results remain valid statements about that proposed
route.  Theorem 1.1 bypasses the route before a prime-dependent shell is
formed: global lower and upper Newton envelopes are already forced to a
common unequal-weight threshold.

This also explains why the bounded searches kept finding only coordinate
and Ferrers-type survivors.  They were resolving coefficient schemes for a
support-overlap event which shifted-ray endpoint rigidity rules out directly
on every actual pure-zero stratum.

The result is specific to two variables in two places: the Hall-deficient
set localizes at one repeated root direction, and every unequal-weight face
projects to an interval on one axis.  Neither step extends formally to the
higher-dimensional GVC.
