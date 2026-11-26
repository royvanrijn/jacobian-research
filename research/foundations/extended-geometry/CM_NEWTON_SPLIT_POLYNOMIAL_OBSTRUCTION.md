# A polynomial trace obstruction at split real-multiplication primes

This note concerns the one-bad-pair split pattern in the
[three-balanced-pair problem](CM_NEWTON_THREE_BALANCED_GLOBAL.md).
It proves an obstruction for fixed polynomial tests, including
on realizable finite-field data. The varieties in the latter
construction vary with the residue prime. The argument makes
no assertion that a fixed characteristic-zero variety realizes
the pattern on a set of positive density.

## The formal calculation

Give the three projective trace variables \(X_1,X_2,X_3\) the
three possible valuation vectors obtained by permuting
\((1,-1,-1)\). For a fixed vector use \(X_i=t^{w_i}u_i\)
in \(k(u_1,u_2,u_3)((t))\), with algebraically independent
\(u_i\) and coefficients in the ground field \(k\).
Then
\[
 v_t\bigl(P(t^{w_1}u_1,t^{w_2}u_2,t^{w_3}u_3)\bigr)
  =\min_{\alpha\in\operatorname{supp}P} w\cdot\alpha ,
\]
since the terms of that least weight give distinct monomials
in \(u_1,u_2,u_3\). If a polynomial is integral in all three
models, each of its supported exponents \((a,b,c)\) satisfies
\[
 a\ge b+c,\quad b\ge a+c,\quad c\ge a+b.
\]
Their sum forces \(a=b=c=0\). Thus the common integral
subring of the polynomial algebra is just the constants.
Algebraic independence here is a device for detecting
cancellation in the formal model. It is not an assumption
on the reductions of the traces of a fixed variety.

If the variables are inverted there are nonconstant
valuation-integral expressions. For example
\(1/(X_1X_2X_3)\) has valuation one in all three patterns.
Its archimedean value need not be bounded when a trace
approaches zero, and its denominator can introduce
nonintegrality at primes other than \(p\). The formal
polynomial obstruction alone says nothing about such
rational expressions with further global information.

## Avoiding fixed cancellations on actual Frobenius data

Fix a totally real cubic field \(F\). Restrict \(p\) to
sufficiently large primes split in the normal closure of
\(F(i)\), congruent to one modulo a fixed auxiliary level
\(M=\ell^{2N}\) with \(\ell\) odd. These conditions give an
infinite set by Chebotarev. Choose a place over \(p\) in that
normal closure to label the local embeddings by the three
real embeddings. At any selected one of the three places
\(\mathfrak q\) of \(F\) above \(p\), the following
construction realizes exactly one bad pair.

For any fixed target \(x\in(-2,2)^3\), there are algebraic
integers \(a\in\mathfrak q\) with \(a\equiv2\bmod M\)
whose normalized archimedean trace vectors
\((\sigma_i(a)/\sqrt p)_i\) tend to \(x\). Moreover one
can require simultaneously that any fixed finite family
of nonzero polynomials over \(\mathbb Q\) in the three
leading residues
\[
 u_i=\begin{cases}
       (a/p)\bmod\mathfrak q,& i\text{ labels }\mathfrak q,\\
       a\bmod\mathfrak q_i,& i\text{ labels another place}
     \end{cases}
\]
have nonzero values after reduction. These residues are
viewed in \(\mathbb F_p\).

Indeed, using finitely many ideal class representatives
and the log-unit lattice, choose a \(\mathbb Z\)-basis of
\(\mathfrak q\) with embedding matrix of entries
\(O_F(p^{1/3})\). Explicitly write
\(\mathfrak q=\beta J\) for a representative \(J\) and multiply
\(\beta\) by a unit so that each of its real absolute values
is \(O_F(|N_{F/\mathbb Q}(\beta)|^{1/3})\); use the scaled
fixed basis of \(J\). The congruence \(a\equiv2\bmod M\)
is a single coset modulo \(M\mathfrak q\) in \(\mathfrak q\).
A bounded fundamental parallelepiped for this lattice lets
us approximate \(\sqrt p\,x\) with error \(O_{F,M}(p^{1/3})\).
The displayed local residue map induces an isomorphism
\[
       \mathfrak q/p\mathfrak q \simeq \mathbb F_p^3.
\]
This is local multiplication by \(p^{-1}\) at the selected
place and the usual residue map at the other two places.
The basis of \(M\mathfrak q\) maps to an invertible
\(\mathbb F_p\)-basis since \(p\nmid M\).

Multiply the finitely many prescribed polynomials and the
three coordinate functions. For all large \(p\) the resulting
polynomial is nonzero. Substitute the affine change of coordinates
from the chosen coset point and basis. A nonzero polynomial of
total degree \(D\) cannot vanish everywhere on the grid
\(\{0,\ldots,D\}^3\) for \(p>D\), by iterated univariate
interpolation. Thus a bounded grid displacement gives all
required nonvanishing while the normalized archimedean error
remains \(O_{F,M,D}(p^{-1/6})\). In particular \(a\) has
valuation precisely one at \(\mathfrak q\) and is a unit
at the other two places. This argument does not prescribe
arbitrary residues modulo the varying prime \(p\).

For completeness these trace data are realized by actual abelian
varieties. The trace \(a\ne0\) generates \(F\), since a rational
integer in \(\mathfrak q\) and of magnitude \(<2\sqrt p\) must
be zero for \(p>4\), and \(F\) has prime degree. Each conjugate
of \(X^2-aX+p\) has negative discriminant. For a root \(\pi\),
the identity \(a=\pi+p/\pi\) gives \(F\subset\mathbb Q(\pi)\).
Thus \(\pi\) is a degree-six Weil \(p\)-integer. At \(\mathfrak q\) the quadratic
is Eisenstein with slopes \(1/2,1/2\). At the other two places
it splits by Hensel's lemma with slopes \(0,1\).
The local endomorphism invariants of its Honda–Tate class vanish:
at the selected place the ramified quadratic degree multiplies
\(1/2\) by two, and the other slopes are integral. Consequently
the simple threefold \(B_p/\mathbb F_p\) has commutative rational
endomorphism field \(F(\pi)\). Waterhouse lets us choose maximal
endomorphism order. In it \((\pi-1)/\ell^N\) is integral, by its
monic equation
\[
 Z^2+\frac{2-a}{\ell^N}Z
          +\frac{1-a+p}{\ell^{2N}}=0.
\]
Thus \(B_p[\ell^N]\) is rational. On \(B_p^2\) the diagonal \(F\)
action and the standard rotation of the two factors give the
Rosati-compatible rank-two action of \(F(i)\).
Its one-bad-pair reduction is the asserted actual model.
[Tate, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0/)
· [Waterhouse, Theorem 6.1](https://www.numdam.org/item/ASENS_1969_4_2_4_521_0/).

The projective trace on a conjugate pair of this model is \(a^2/p\).
In the selected local orientation its leading residues, after
separately removing \(p^{\pm1}\), are \(u_i^2\).
The formal initial form of any fixed projective polynomial stays
nonzero as a polynomial in the \(u_i\), since substituting squares
is injective. Apply the grid argument to these initial forms.
Thus nonconstant fixed projective polynomials actually have
negative valuation on suitable one-bad models, including models
at arbitrary fixed auxiliary torsion level. Avoiding cancellation
for a finite list by this construction is weaker than residue
independence on the Frobenius sequence of a fixed threefold.

## Polynomial tests with explicit powers of the residue prime

This also rules out a larger bounded-integer test. Fix
\(Q(U_1,U_2,U_3,S)\in\mathbb Q[U_1,U_2,U_3,S]\) and an integer
\(n\ge0\), and evaluate
\[
                   p^{-n}Q(\sigma_1(a),\sigma_2(a),\sigma_3(a),p).
\]
Unless the normalized expression is a rational constant,
it cannot be both uniformly archimedean bounded and integral
at every place above \(p\) on all the one-bad finite-field models
just constructed.

Give \(U_i\) weight one and \(S\) weight two. If the largest
weighted degree \(D\) exceeds \(2n\), evaluation of the
highest weighted part at \((x,1)\) is a nonzero polynomial
in \(x\). The grid and archimedean approximation with a target
where that polynomial is nonzero give growth of order
\(p^{D/2-n}\). If \(D\le2n\), consider any supported monomial
\(U^\alpha S^m\) other than \(S^n\).
Writing \(k=n-m\), the weight condition gives
\(k>0\) and \(|\alpha|\le 2k\).
In the model where label \(i\) is bad, its normalized
valuation is \(\alpha_i-k\). There must be a negative
such valuation, since requiring \(\alpha_i\ge k\) for
all three labels would contradict \(|\alpha|\le2k\).
The initial residue polynomial for that orientation is
nonzero: distinct \(\alpha\)'s give distinct monomials,
and different \(m\)'s with the same \(\alpha\) have
different valuations. The grid argument again gives
actual negative valuation on suitable models. This covers
polynomial constructions using a fixed finite list of
Frobenius power traces at fixed weights, since the traces
of powers satisfy the rank-two characteristic recurrence.

## What remains global

In the independent case, over an algebraic closure the
connected derived monodromy surjects onto three independent
\(\mathrm{PGL}_2\) factors. The three projective trace functions
therefore have Zariski-dense joint image in affine three-space.
An exact polynomial identity among connected Frobenius traces
would contradict this. It does not exclude congruences modulo
the varying residue prime. Conversely the finite-field
models above keep the one-bad pattern while approaching
arbitrary normalized archimedean traces and matching a
prescribed auxiliary trace congruence. The companion
[three-pair proof](CM_NEWTON_THREE_BALANCED_GLOBAL.md)
strengthens the latter to auxiliary matrix classes. This
is why bounded-integer polynomial tests, fixed-level
Chebotarev and weak archimedean equidistribution do not by
themselves produce the missing split ordinary density.
They do not prevent a theorem on correlations with the
varying residue prime of a fixed compatible system.

Fité's cubic real-multiplication proof uses a place of
inert coefficient degree, where a single local failure carries
three units of divisibility. The inert case of Wang–Zhang
uses the same kind of gain, and its equal-degree two-place
alternative does not apply to a cubic field. If the base
contains the coefficient field, an inert rational prime
cannot be of absolute degree one there. Passing to a larger
base preserves that obstruction to the inert method. The
linked projective-field theorem removes another global stratum;
its hypothesis fails in the independent split case.
[Fité, Theorem 3 and §6](https://arxiv.org/html/2005.10185v3)
· [Wang–Zhang, Theorems 1.3 and 1.7–1.8](https://doi.org/10.1112/mtk.70081)
· [Projective trace field theorem](CM_NEWTON_PROJECTIVE_TRACE_FIELD.md).

Suh's Hilbert modular comparison displays the same residual
case. His §5.3 treats cyclic cubic base and coefficient fields
that coincide. When the smallest Frobenius field still has
degree three, his bound there allows exactly one nonordinary
component. The Sato–Tate argument of Theorem 4.1.1(5) addresses
the equality threshold in even degree. Its bisection hypothesis
in part (6) requires two equally sized coefficient-field
orbits, so it does not cover three independent split factors.
[Suh, Theorem 4.1.1 and §5.3](https://arxiv.org/html/2410.01182v2).

## Rational functions of the trace and determinant

Fix \(F\) and an auxiliary torsion level as in the construction above.
Write \(S\) for the determinant variable and \(U_i\) for the three
labelled trace variables. One can instead use the coordinates of the
relative trace in an integral basis of \(F\). The change of variables
is invertible over the normal closure of \(F\). For a rational function
\(R\) in either set of coordinates, take its coefficients in a fixed
number field \(L\) containing the normal closure of \(F(i)\) and
restrict to \(p\) split completely in \(L\).

Suppose the defined values of \(R(p,a)\) on all sufficiently large
minimal split one-failure models are uniformly bounded at a fixed
archimedean embedding and integral at all places outside a fixed
finite set of rational primes. Then \(R\) is constant. Here minimal
means that the valuation of \(a\) at the failed coefficient prime
is one and its valuations at the other two are zero. Both
\(p\equiv1\bmod M\) and \(a\equiv2\bmod M\) may be retained for
the fixed level \(M=\ell_0^{2N}\).

First write \(R=P/Q\) in lowest terms. If \(Q\) has an irreducible
factor other than \(S\), choose a smooth algebraic point on a
geometric component of that factor with \(S\ne0\) and \(P\ne0\).
Such a point exists because \(P,Q\) are coprime, a property preserved
by extension of the coefficient field. Choose an auxiliary prime
\(\ell\) split in the normal closure of the field of the point and
the coefficients, outside the fixed exceptional primes and coprime
to \(M\). Exclude the finitely many primes where the point coordinates
or polynomial coefficients fail to be integral, or where \(SP\)
reduces to zero at the point. The point then prescribes a nonzero determinant
\(\bar s\in\mathbb F_\ell^\times\) and a trace residue class with
\(P\) a unit and \(Q\) zero modulo the selected place above \(\ell\).

Require \(p\equiv\bar s\bmod\ell\) and complete splitting in the
normal closure of \(L(\zeta_M)\). These conditions are compatible
by Chebotarev since \(\mathbb Q(\zeta_\ell)\) is totally ramified
at \(\ell\), while that normal closure is unramified there. The
trace residues and \(a\equiv2\bmod M\) combine into a single
class modulo \(\ell M\mathcal O_F\). In labelled coordinates this
uses the split isomorphism
\(\mathcal O_F/\ell\mathcal O_F\simeq\mathbb F_\ell^3\).
The trace lattice construction with the fixed modulus \(\ell M\)
and the Honda–Tate realization above preserve these congruences,
any selected failed label, and approximation to any fixed
normalized archimedean target in the open trace cube.
[Auxiliary pole construction](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#imposing-a-pole-on-finite-field-models).

To keep \(Q(p,a)\ne0\), take an archimedean target at which the
highest weighted part of \(Q\), evaluated at \(S=1\), is nonzero.
Use weights \(1\) for \(S\) and \(1/2\) for each \(U_i\). At any
fixed weight, the \(U\)-exponents determine the \(S\)-exponent, so
this highest part is a nonzero polynomial after setting \(S=1\).
It has a nonzero value somewhere in the open trace cube. Along
the corresponding approximations the weighted leading term
keeps \(Q(p,a)\) nonzero. Consequently the models give a defined
value with negative valuation above \(\ell\), contradicting the
integrality assumption. It follows that \(R\in L[S,S^{-1},U_1,U_2,U_3]\).

For a supported monomial \(S^kU^\alpha\) in this Laurent expression,
its valuation in the model with label \(i\) failed is \(k+\alpha_i\)
at sufficiently large \(p\). Among the monomials of minimum
valuation for label \(i\), the \(U\)-exponents are distinct,
since \(\alpha_i\) determines \(k\) at fixed valuation. Their
initial form is therefore a nonzero polynomial in the leading
residues \(u_i\) of the construction above. The grid argument
chooses a model where it survives, so integrality forces
\(k+\alpha_i\ge0\) for every supported monomial and all three labels.

Archimedean boundedness likewise forces \(k+|\alpha|/2\le0\)
for all supported monomials. If the maximum weight were positive,
its nonzero initial polynomial on the trace cube would give
unbounded growth along the same archimedean approximations.
If \(k<0\), the three valuation inequalities give
\(|\alpha|\ge-3k\), contrary to \(k+|\alpha|/2\le0\).
For \(k\ge0\), the weight inequality leaves just \(k=0\)
and \(\alpha=0\). Thus only constants satisfying the assumed
integrality condition remain. In particular, for \(R\) over
\(\mathbb Q\) these are rational constants, and with integrality
at all places they are integers.
