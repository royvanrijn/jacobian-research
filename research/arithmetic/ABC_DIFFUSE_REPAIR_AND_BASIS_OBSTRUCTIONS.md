# Repair costs and valuation bases for broad abc profiles

Two scalar tests describe the limits of transferring a broad valuation
profile to a powerful-number count. The first optimizes the actual
mandatory common multiplier in the bounds (GU2) and (WH1). The second
quantifies the cost of encoding many exact valuation levels in a
fixed number of short roots.

## The mandatory multiplier and unequal fullness orders

Use the height-first normalization \(n=X^L\). Consider the formal
profile in which all three numerical heights are one, and every
exact scaled valuation measure tends to
\[
                         \nu=\tfrac13\operatorname{Unif}[2,4].
\]
Its mass is \(d=1/3\) and its first moment is one. The following
calculation tests the bounds on those prescribed measures. An actual
family would also have to satisfy the original primitive additive
relation with all three radical caps.

Write \(m\ge x\ge y>0\) for the three scaled orders \(Ls_i\),
so that \(M=m/L\) is the maximum prescribed order. Define
\[
 R(v)=\nu([0,v)),\qquad
 Q(m,x)=\frac1m\int_{t<m}\max(m-t,x)\,d\nu(t).
\]
The normalized contribution of the mandatory multiplier is
\[
                \frac{\log_X Z}{LM}=Q(m,x)+R(x)+R(y).
                                                               \tag{DB1}
\]
To check this, at a deficient prime on the first role the
deficiency on its own role is \(M-e\), and the other two
primitive valuations are zero. The required exponent of the
common multiplier is therefore \(\max(M-e,x/L)\).
At a deficient prime on either other role it is \(M\),
since the first role is absent at that prime. The formula
also handles ties: if \(x=m\), then \(Q(m,x)=R(m)\).
The cutoff has no atom in the profile, so floors and ceilings
have no effect on the limiting integrals. This calculation
uses the actual common multiplier that makes all three
multiplied values full.
[Mandatory multiplier and primitive quotient bins (GU1–GU2)](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#an-arbitrary-order-upper-bound-with-common-divisors)

**DB1.** In the exponent supplied by (GU2) put \(a=1/3-\delta\).
The infimum over all three scaled orders is
\[
 \inf_{m\ge x\ge y>0}
    \left\{a\left(\frac1m+\frac1x+\frac1y\right)
                       +Q(m,x)+R(x)+R(y)\right\}
                         =\frac{3a}{2}.
\]
It is attained at \(m=x=y=2\). Thus this scalar transfer on the
formal profile gives an exponent below \(L/3\) exactly when the
gap in (GU2) satisfies \(\delta>1/9\).

It suffices to prove the lower bound at \(a=d\). For
\(0\le a\le d\), the score is \(a/d\) times the score at \(d\)
plus \(1-a/d\) times the nonnegative repair cost. If \(m\le2\),
the result follows from the three reciprocals. Suppose \(m\ge2\).
For \(x\ge2\), both nonmaximal terms \(d/v+R(v)\) are minimized
at \(v=2\), which is compatible with their ordering, and
\(Q(m,x)\) is nondecreasing in \(x\). Indeed
\(R(v)=(v-2)/6\) for \(2\le v\le4\) and \(R(v)=d\) for
\(v\ge4\). For \(x\le2\) the best \(y\le x\) is \(y=x\).
Here \(2d/x+Q(m,x)\) decreases as \(x\) increases up to two.
Where derivatives exist, that of \(Q\) with respect to \(x\)
is at most \(R(m)/m\le1/12\), whereas that of \(2d/x\)
is at most \(-d/2=-1/6\). Thus take \(x=y=2\).

The remaining term \(d/m+Q(m,2)\) equals
\[
 \begin{cases}
    d-d/m,&2\le m\le4,\\
    m/12-2/3+(d+2)/m,&4\le m\le6,\\
    d+(d-1)/m,&m\ge6 .
 \end{cases}
\]
Every branch is at least \(d/2\). For the middle branch the
arithmetic-geometric mean inequality gives
\(\sqrt{(d+2)/3}-2/3>d/2\). Equality in the full bound holds
at all three orders equal to two.

## The weighted-height sector with finite slack

**DB2.** For (WH1), write \(0\le\tau\le1/8\) for the
weighted-height tolerance, \(0<\kappa<1\) for its gap, and
\(R_0=1+\tau/2\). On the same formal profile the optimal
coefficient in base \(n\) is
\[
       \min_{2\le m\le 2R_0}
               \left(\frac{1-\kappa}{m}
                         +\frac{m-2}{6R_0}\right).
\]
The minimum is attained by clipping
\(\sqrt{6R_0(1-\kappa)}\) to \([2,2R_0]\). A coefficient
strictly below \(1/3\) is possible exactly when
\[
                 \kappa>\frac13-\frac{(R_0-1)^2}{6R_0}.
\]
In particular the required limiting gap at zero slack is
\(\kappa>1/3\), and small positive slack moves the threshold
only quadratically.
[Weighted-height quotient estimate (WH1)](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-quotient-bound-in-the-weighted-height-sector)

The constraint is \(m/x+m/y\le2+\tau\). Both nonmaximal
orders are therefore at least \(m/(1+\tau)\). If \(m\ge4\),
their costs \(R(x)+R(y)\) alone are at least \(14/27>1/2\),
since \(m/(1+\tau)\ge32/9\). The choice of all three orders
equal to two has score \((1-\kappa)/2\le1/2\).
For \(m\le4\), the function
\(t\mapsto R(m/t)=(m/t-2)_+/6\) on \(t\ge1\) is convex
and nonincreasing. Jensen gives
\(R(x)+R(y)\ge2R(m/R_0)\). Also \(x\ge m/R_0\) because
\(x\ge y\). Thus the same choice \(x=y=m/R_0\) minimizes
the nonmaximal costs and \(Q(m,x)\). For \(m\le4\) and the
stated slack one has \(m/R_0\ge m-2\), so
\(Q(m,m/R_0)=R(m)/R_0\). The score decreases up to \(m=2\)
on \(m\le2\). On \(2R_0\le m\le4\) it is increasing,
since its derivative is
\(1/(2R_0)-(1-\kappa)/m^2>0\).
This leaves exactly the displayed interval. There the score
is below \(1/3\) for some \(m\) precisely when
\[
       1-\kappa<
       \max_{2\le m\le2R_0}\frac{m(2R_0+2-m)}{6R_0}
                       =\frac{(R_0+1)^2}{6R_0}.
\]

Separate role-local multipliers would instead charge
\(m-t\) on every deficient prime at equal order \(m\).
Their combined numerical penalty is \((m-2)^2/4\) on
\(2\le m\le4\). That penalty belongs to a weighted additive
equation with three separately changed terms. Applying a
powerful-number bound to that equation requires control of
the coefficient records in addition to the numerical penalty.

## Entropy and the cost of an adaptive valuation basis

The next test applies to actual integer factorizations. Let
\(T>1\) have normalized radical mass \(r\), defined by
\(\log\operatorname{rad}(T)=rL\log X\). Choose a prime
of \(T\) with probability \(\log p/(rL\log X)\), and write
\(E=v_p(T)\) for its integer degree. Assume
\[
                         H(E)\ge\log(1/L)-C
\]
in natural entropy units. For example this holds if each exact
degree has probability \(O(L)\). A diffuse weak scaled profile
by itself is a weaker condition on the discrete degree law.

**DB3.** Consider any exact representation
\[
              T=V\prod_{j=1}^m z_j^{q_j}
\]
by positive rationals, with fixed integer orders \(q_j\) for
the represented value. The root supports may overlap, and
the orders and root choices may depend on the value. Put
\(B_j=\mathbb E|v_p(z_j)|\) under the prime law above. If the
ordinary rational height of the correction satisfies
\(H(V)\le X^\omega\), then
\[
      \left(L+\frac{2\omega}{r}\right)
                       \prod_{j=1}^m(1+B_j)
                  \ \ge\ e^{-C}(2e)^{-(m+1)}.                 \tag{DB3}
\]
Here \(e\) in the numerical constant denotes the base of natural
logarithms.

Indeed set \(K_j=v_p(z_j)\) and \(D=v_p(V)\). The identity
\(E=\sum_j q_jK_j+D\) and entropy subadditivity give
\(H(E)\le\sum_j H(K_j)+H(D)\). A signed integer variable
\(J\) of mean absolute size \(b\) has
\(H(J)\le\log(2e(1+b))\): take its magnitude, compare with
the geometric law of mean \(b\) by nonnegativity of relative
entropy, and use at most \(\log2\) for its sign. The sum of
the numerator and denominator logarithms of \(V\) is at most
\(2\omega\log X\), so \(\mathbb E|D|\le2\omega/(rL)\).
Substitution gives the displayed bound.

The order in \(L+\omega\) is sharp for fixed rank and bounded
literal degrees \(E\le C_1/L\), with \(r\) bounded away from
zero. For \(0\le\omega\le1\), take
\(R_1=1+\lfloor\omega/(rL)\rfloor\), and put \(E\bmod R_1\)
into the integral correction. Write \(\lfloor E/R_1\rfloor\)
in \(m\) digits to base
\[
              b=\max\!\left(2,
                     \left\lceil(1+\lfloor C_1/(LR_1)\rfloor)^{1/m}
                                                              \right\rceil\right).
\]
Use the digit in position \(j\) as the prime multiplicity in
\(z_j\), with order \(q_j=R_1b^{j-1}\).
The correction has numerical height exponent at most
\(rL(R_1-1)\le\omega\), every \(B_j\le b-1\), and
\(b^m=O_{C_1,r,m}(1/(L+\omega))\).

For a full modulus of numerical height \(w\) and scaled
fullness \(\ell L\to\kappa_0>0\), a valid root-height
coefficient in (UW) satisfies
\(C_j\ge (1+o(1))\kappa_0rB_j/(2w)\), as the radical-weighted
sum of absolute valuations of \(z_j\) is bounded by its
two rational height logarithms. Thus the known product
dependence in (UW) is relevant to this tradeoff. When the
correction width is at least \(L\), (DB3) bounds
\(\omega\prod_j C_j\) away from zero for fixed rank and fixed
positive modulus parameters, using \(C_j\ge1\).
At widths \(o(L)\), the same tradeoff forces a height product
of order at least \(1/L\). A comparison with the actual admissible
constants in a determinant theorem is a separate question.
[Adaptive finite-basis approximation (AB1–AB2)](ABC_WIDE_COFACTOR_CONGRUENCES.md#adaptive-finite-basis-approximation-of-prime-valuations)
· [Unequal root-height coefficients (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)
