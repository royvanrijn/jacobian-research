# Scalar gauges, the subprincipal corner, and polynomial admission

Scalar formal-unit gauges preserve local finite-order string equations
and shift the subprincipal corner. In contrast, a nonconstant rational
scalar gauge cannot carry both members of an actual polynomial Weyl
pair to polynomial operators. The latter is a genuine global admission
constraint, proved below by finite-dimensional traces. It does not force
the subprincipal corner of an unknown pair to vanish.

Work over a characteristic-zero field \(k\), with
\(A_1(k)=k[x][D]\), \([D,x]=1\). All operators have finite
differential order. Constants such as a nonzero Weyl parameter may be
included in \(k\).

## Formal units and the mixed corner

Let \(u=x-c\), choose a formal completion at a point where \(u\)
is nonzero, and normalize the unit's constant term to one. For every
\(\gamma\in k\), the binomial series defines a scalar unit
\(f_\gamma\) with \(f_\gamma'/f_\gamma=\gamma/u\).
For example, at \(x=0\) take \(u=1+x\) and
\(f_\gamma=(1+x)^\gamma\). Its conjugation is
\[
 \Gamma_\gamma(T)=f_\gamma^{-1}Tf_\gamma,\qquad
 \Gamma_\gamma(x)=x,\quad
                \Gamma_\gamma(D)=D+\frac\gamma u.       \tag{1}
\]
This is an automorphism of the finite-order formal differential algebra.
It preserves any scalar relation \([Q,P]=\kappa\), and for each
nonzero operator preserves its differential order and highest
coefficient. In particular it does so for every fixed noncommutative
chopping polynomial in \(P,Q\). If the highest coefficient is a
formal unit, the usual equality of differential order with the
\(\max(j-\operatorname{val}_x a_j)\) order is preserved as well.
This does not assert global polynomial coefficients after conjugation.

Suppose the normal PBW operator \(P\) is supported in the rectangle
\(0\leq i\leq a,0\leq j\leq b\), its corner coefficient at
\(x^aD^b\) is one, and its mixed coefficient at
\(x^{a-1}D^{b-1}\) is \(A\). Define
\[
                              A_W=A-\frac{ab}{2}.        \tag{2}
\]
This is the mixed coefficient of the Weyl symbol. The coefficient of
\(D^{b-1}\) in (1) is changed by \(b\gamma p_b(x)/u\), where
\(p_b\) is the highest coefficient of \(P\). Its expansion at
infinity therefore gives exactly
\[
                 A\longmapsto A+b\gamma,\qquad
                 A_W\longmapsto A_W+b\gamma.            \tag{3}
\]
The coefficient here is a Laurent coefficient at infinity if the
gauged operator is rational rather than polynomial. The leading
corner and the immediately adjacent edge coefficients are unchanged.

For comparison, let \(P,Q\) have monic rectangular corners
\(x^aD^b,x^cD^d\), where all four integers are positive,
\(ad=bc\), \(a\ne b\), and
\(\rho=c/a=d/b\). Assume the adjacent coefficients at offsets
\((1,0)\) and \((0,1)\) from each corner vanish. If the normal mixed
coefficients are \(A,B\), the coefficient of
\(x^{a+c-2}D^{b+d-2}\) in \([Q,P]\) is
\[
                (a-b)(\rho A_W-B_W),\qquad
                B_W=B-\frac{cd}{2}.                    \tag{4}
\]
Indeed the one-contraction terms contribute
\((c-d)A+(b-a)B\), and the two-contraction corner terms contribute
\(
 [d(d-1)a(a-1)-b(b-1)c(c-1)]/2
\).
Their sum is (4); larger contractions cannot reach that coefficient.
The adjacent-edge hypotheses remove the other one-contraction terms.

The target monomial in (4) is nonconstant, so a Weyl pair satisfies
\(B_W=\rho A_W\). Under (1), the two mixed coefficients shift by
\(b\gamma,d\gamma\), preserving this relation. It does **not**
give \(A_W=0\). This is the distinction relevant to the
[PBW-corrected final coefficient comparison](PBW_TRANSPORT_AND_FINAL_COEFFICIENT_COMPARISON.md).

## Rational scalar gauges cannot land both members globally

**Theorem.** Let \(P,Q\in A_1(k)\), \([Q,P]=\kappa\in k^\times\).
If \(f\in k(x)^\times\) satisfies
\[
                f^{-1}Pf,\quad f^{-1}Qf\ \in A_1(k),    \tag{5}
\]
then \(f\in k^\times\).

For a nonconstant polynomial \(f\), (5) makes the ideal \(fk[x]\)
invariant under both \(P,Q\), since \(P(fv)=f(f^{-1}Pf)v\).
Their induced operators on \(k[x]/fk[x]\) satisfy the same scalar
commutator. Its trace is simultaneously zero and
\(\kappa\deg f\ne0\), a contradiction.

For the full rational statement write \(f=g/h\) with nonzero coprime
polynomials. In \(k(x)\), both lattices
\(O=k[x]\) and \(L=fk[x]\) are invariant. Their intersection is
\(gk[x]\): if \(hv=gw\) with \(v\in k[x]\), coprimality forces
\(g\mid v\). Consequently the two invariant finite-dimensional
quotients have dimensions
\[
 \dim_k O/(O\cap L)=\deg g,\qquad
 \dim_k L/(O\cap L)=\deg h.                             \tag{6}
\]
For the second equality use \(L/gk[x]\simeq k[x]/hk[x]\).
If either dimension is positive, the same trace contradiction applies.
Thus \(g,h\) are constant. This proof does not assume the image pair
generates the entire Weyl algebra.

The reverse conjugation \(fPf^{-1},fQf^{-1}\) is also excluded for
nonconstant rational \(f\), by interchanging the two pairs. In
particular, every nonzero integer shift \(f=u^m\), positive or
negative, fails to keep both members of an actual polynomial pair in
the polynomial chart. This theorem does not cover nonrational algebraic
or formal \(u^\gamma\), and gives no classification for such gauges.

The same argument gives a one-position matrix statement. If
\(P,Q\in\operatorname{Mat}_r(A_1(k))\),
\([Q,P]=\kappa I_r\), and \(G\in\operatorname{GL}_r(k(x))\)
has both \(G^{-1}PG,G^{-1}QG\) polynomial, then
\[
                             G\in\operatorname{GL}_r(k[x]).
                                                               \tag{7}
\]
Here the conclusion means a polynomial matrix with polynomial inverse;
it permits nonconstant unimodular matrices. Indeed \(O=k[x]^r\)
and \(L=GO\) are invariant. A common polynomial denominator for
\(G,G^{-1}\) shows that \(O\cap L\) has finite codimension in
both lattices. The two quotient traces force both quotients to be
zero, so \(O=L\), proving (7). This uses a single position variable:
in several variables the analogous polynomial quotients need not be
finite-dimensional. It makes no claim about gauges in the full Weyl
division ring.

## Exact controls and scope

[The checker](check_scalar_gauge_admission.py) uses independent finite
Ore multiplication over Laurent polynomial coefficients. It verifies
the mixed-corner shift for several integral and fractional \(\gamma\),
the convention \(A_W=A-ab/2\), and (4). Its genuine Weyl-pair control is
\([D,u]=1\): conjugation gives \(D+\gamma/u,u\), still satisfying
the relation, while the first operator loses polynomial admission for
every tested nonzero \(\gamma\). Matrix controls retain the allowed
unimodular gauge \(\left(\begin{smallmatrix}1&u\\0&1\end{smallmatrix}\right)\)
and show the loss of polynomial admission for
\(\operatorname{diag}(u,u^{-1})\) on the scalar matrix Weyl pair.

The corner controls are not asserted to be Weyl pairs; the full
rational-gauge theorem follows from the written trace proof, not finite
sampling. [The receipt](scalar_gauge_admission_controls.json) pins this
note and checker. No general matrix-gauge no-go, absolute \(A_W=0\)
theorem, noninteger-gauge classification, new finite lift, or DC1
resolution is claimed.
