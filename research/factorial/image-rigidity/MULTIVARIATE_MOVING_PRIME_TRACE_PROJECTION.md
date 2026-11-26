# Moving prime depths in a rational trace projection

Fix a hypothetical ordinary factorial common factor and use an unordered
Galois-invariant multiset of its reciprocal traces
\(t_{ij}=\alpha_i/\alpha_j+\alpha_j/\alpha_i\). Let \(T\) be their monic
rational polynomial of degree \(m\) and \(D\) its least positive
denominator. The [denominator identity](MULTIVARIATE_TRACE_HERMITE_JET_HEIGHT.md)
gives \(v_p(D)\) as the sum of the absolute valuation spreads of
the selected root pairs. In particular \(T\) has integral coefficients at a prime at which
each selected pair has equal root valuations.

Assume the traces approach a fixed algebraic selection
\(T_*=\prod_a(Y-\theta_a)^{\nu_a}\), where the trace rates at \(\theta_a\)
are \(O(\eta^{\gamma_{a,i}})\) with \(\gamma_{a,i}>0\) in increasing
order. Here \(\eta\to0\) can be the maximum added relative gap at a
positive diffuse contact. By the
[Hermite jet bounds](MULTIVARIATE_TRACE_HERMITE_JET_HEIGHT.md), the forms
\[
 J_{a,k}(P)=P^{(k)}(\theta_a),\qquad
 w_{a,k}=\sum_{i=1}^{\nu_a-k}\gamma_{a,i}
 \qquad (0\le k<\nu_a)
\]
satisfy \(|J_{a,k}(DT)|=O(D\eta^{w_{a,k}})\). Their common kernel on
polynomials of degree at most \(m\) is the line through \(T_*\).

Take a fixed rational linear subspace \(U\) of that coefficient
space with \(T_*\notin U\otimes\overline{\mathbb Q}\). Let \(M\) be an
integer matrix with kernel \(U\) and full row rank, and put
\[
      y_* = M T_*,\qquad x=M(DT),\qquad g=\gcd(x)>0 .
\]
Use the positive gcd of the integer coordinates. The vector \(x\)
is nonzero near the reference.

For each weight \(w\) let \(E_w\) be the span, over the algebraic
numbers, of the forms \(J_{a,k}\) with \(w_{a,k}\ge w\). The forms
in \(E_w\) annihilating \(U\) descend to linear forms on the
quotient. For a rational subspace \(W\) of the quotient containing
\(y_*\) after extension of scalars, with \(r=\dim W\ge2\), let
\(S_U(W)\) be the largest sum of weights of \(r-1\) independent
forms on \(W\) obtainable in this way. A form may use weight \(w\)
if it descends from \(E_w\). Set
\[
 A_U=\min_{\substack{W\ni y_*\\ \dim W\ge2}}
                                \frac{S_U(W)}{\dim W}.
\]
The minimum belongs to a finite list of weight sums divided by
integers. The forms of smallest weight span the annihilator of
\(y_*\) in the quotient, so \(A_U>0\) if the quotient has
dimension at least two.

**Theorem.** For every \(\epsilon>0\), sufficiently close to \(T_*\),
either \(MT\) is projectively equal to \(y_*\), or
\[
                \log D \ \ge\ \log g+
                        (A_U-\epsilon)\log(1/\eta)-O(1).       \tag{1}
\]
There is no restriction on the primes supporting \(g\). If at a prime
\(p\) there is a coefficient vector \(T_p\in U\otimes
\overline{\mathbb Q}_p\) satisfying the coefficientwise congruence
\[
                          v_p(DT-T_p)\ge\kappa_p\ge0 ,
                                                                    \tag{2}
\]
then \(v_p(g)\ge\lceil\kappa_p\rceil\). In particular (1) can charge
the sum of these depths times \(\log p\) at arbitrarily moving primes.
The case of a one-dimensional quotient is exact projective attainment
throughout.

## Proof

Since \(x/D\to y_*\ne0\), the primitive vector \(x/g\) has
projective logarithmic height \(\log(D/g)+O(1)\). Every form
descending from \(E_w\) has magnitude \(O(D\eta^w)\) on \(x\),
hence its magnitude divided by \(\|x\|_\infty\) is \(O(\eta^w)\).
Apply the projective
[Subspace Theorem, equation (1.1)](https://arxiv.org/pdf/math/0408381)
as in the Hermite jet bound, on a rational subspace \(W\)
containing infinitely many violating points. It necessarily contains
\(y_*\). Take the \(r-1\) independent descended forms of maximum
weight and one additional form nonzero on \(y_*\).
The product of their projectively normalized absolute values is
\(O(\eta^{S_U(W)})\). If (1) fails along a sequence in \(W\), this
is at most \(H(x/g)^{-r-\delta}\) for some \(\delta>0\) and
all sufficiently close points on that sequence. A bounded height
sequence can instead take only finitely many projective values,
so convergence here would force exact attainment. The theorem
puts the violating points in finitely many proper subspaces of
\(W\). The rational spans of their rational points are still
proper. Descend in dimension, retaining an infinite subsequence
at each step. Dimension one gives exact attainment. This proves
(1). Formula (2) implies that every coordinate of
\(x=M(DT-T_p)\) has valuation at least \(\kappa_p\), because \(M\)
is integral. Taking the gcd proves the last assertion.

For an explicit weaker exponent, let \(r_0\) be the rational rank
of the coordinates of \(y_*\) and let \(w_0\) be the smallest
Hermite weight. Every \(W\) above has dimension at least
\(\max(2,r_0)\), so \(A_U\ge(1-1/\max(2,r_0))w_0\).
For rational projective \(y_*\), rational linear forms
annihilating \(y_*\) also give individual elementary bounds.
When such a form, in integral coordinates, descends from
\(E_w\), it either vanishes exactly or its nonzero integer
value at \(x/g\) gives \(\log D\ge\log g+w\log(1/\eta)-O(1)\).

## Local rational target jets

Suppose \(t_0\in\mathbb Q\) is a fixed target and the selected traces
are integral at each certified prime. At each such prime assume
\(t_0\) is integral there and at least \(\nu_p\) traces satisfy
\(v_p(t_{ij}-t_0)\ge\kappa_p>0\). The remaining differences are
integral at that prime. For a fixed integer \(2\le r\le m\) take
\(U=\{F:(Y-t_0)^r\mid F\}\) inside the coefficient space.
Use the first \(r\) Taylor coefficient evaluations at \(t_0\),
clearing their fixed rational denominators, for \(M\).
The product expansion of the \(k\)-th Taylor coefficient gives
\[
       v_p\bigl(T^{(k)}(t_0)/k!\bigr)
                      \ge (\nu_p-k)_+\kappa_p .
\]
Here \(p\nmid D\). The integer factors clearing rational
denominators cannot lower any of these valuations.
Thus \(v_p(g)\ge\lceil(\nu_p-r+1)_+\kappa_p\rceil\), at each of
the moving primes. An exact local source for \(\kappa_p\) is
the [truncated face test for the trace](MULTIVARIATE_MOVING_PRIME_TRACE_RESULTANT.md#a-local-certificate-from-arbitrary-truncated-faces)
with \(f(Y)=Y-t_0\) after clearing its denominator and accounting for
that factor at the prime. In other settings the same construction
works with a fixed span of monic local torsion polynomials. If each
unit trace is within \(p^{-\kappa_p}\) of a selected algebraic
integer, the elementary symmetric expansion makes \(T\) coefficientwise
close to the corresponding monic polynomial to at least that depth.
Its coefficient vector can lie in \(U\) even when the target itself
is outside \(\mathbb Q_p[Y]\), by using the rational span
containing all its algebraic conjugates.

## Comparison with separate bounds

For example, if a fixed monic rational polynomial \(T_0\) with
\(m\) distinct torsion trace roots is the complete local selection,
one can take \(U=\mathbb Q T_0\). The local depth in (1) is
at least \(L=\sum_p\kappa_p\log p\). The earlier resultant with
\(f=T_0\), scaled to be integral, uses \(m\) local trace factors
and divides the resulting depth by \(m\). Hence its depth term
is also \(L\), at primes away from the fixed denominator. If
\(T_0\) and \(T_*\) have disjoint roots, that resultant charges
no real approximation rate. The quotient in (1) still charges
a real rate when its projected rational dimension is at least
two.

Consider a degree-two trace polynomial with simple real limiting
traces \(\theta,-\theta/2\), where \(\theta=\sqrt[3]2\), both at
rate \(\gamma\). At the illustrative local target \(t_0=0\),
suppose both traces have depth at least \(\kappa_p\) at the included
primes and put \(L=\sum_p\kappa_p\log p\).
Projection to \((T(0),T'(0))\) has limiting ratio \(1/\theta\).
The theorem therefore requires
\[
             \log D\ge L+(\gamma/2-\epsilon)\log(1/\eta)-O(1).
\]
The full real coefficient rank is three, giving the separate
Hermite rate \(2\gamma/3\). The degree-one target resultant gives
the separate depth \(2L\). If \(L=(\gamma/4)\log(1/\eta)+O(1)\),
the projected rate \(3\gamma/4\) before \(\epsilon\) exceeds both.
This compares the inequalities under specified trace data and
makes no assertion that an ordinary factorial profile attains them.

Near a given reference selection, a uniform bound on
\(D/\prod_p p^{\lceil\kappa_p\rceil}\), using depths in (2),
already excludes the branch if \(MT_*\) has irrational projective
coordinates. Indeed the bound bounds the height of the primitive
projection, and only finitely many rational points have a bounded
height. This conclusion does not require the Subspace Theorem.
Inequality (1) allows that quotient to grow at a rate strictly below
\(A_U\) in \(\log(1/\eta)\). Both tests permit arbitrary moving
ordinary rows. Their applicability to every remaining tail would
require a uniform source of moving depths or an independent bound
on the projected denominator.
