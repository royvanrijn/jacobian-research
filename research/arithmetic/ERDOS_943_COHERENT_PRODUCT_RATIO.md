# Singular six-cycles with small coherent product ratio

Let \(N\ge3\), and form a simple bipartite graph from actual representations
\[
 a^2b^3+c^2d^3=N,\qquad
 a,b,c,d\in\mathbb Z_{>0},\qquad b,d\ \hbox{squarefree},\qquad
 \gcd(a^2b^3,c^2d^3)=1.                                      \tag{1}
\]
Keep one selected representation for each represented kernel pair \((b,d)\).
Assume that all kernel vertices are at most \(T\), where \(T\ge1\). No bound
on the selected roots \(a,c\), other than that forced by (1), is needed.

Count cyclically ordered simple six-cycles with left vertices
\(b_1,b_2,b_3\), right vertices \(d_1,d_2,d_3\), and edges
\((b_i,d_i),(b_{i+1},d_i)\), with indices modulo three. Such ordered counts
differ from the usual unoriented count by an absolute factor.

A cycle is called singular when its selected edge ratios admit positive
rational vertex labels \(p_i,q_j\) such that
\[
                  a_{ij}/c_{ij}=p_i/q_j
                  \quad\hbox{on every edge}.                  \tag{2}
\]
Equivalently, the alternating product of the six selected edge ratios is
one. Normalize these labels to positive integers of global gcd one, and set
\[
 P=p_1p_2p_3,\quad Q=q_1q_2q_3,\qquad
 B=b_1b_2b_3,\quad D=d_1d_2d_3.                                \tag{3}
\]
The ratio \(P/Q\) is independent of a common rescaling of the labels.
In particular, the condition \(P=Q\) is intrinsic.

## 1. The counted sectors

**Theorem 1.** For every \(\varepsilon>0\), uniformly in \(N\ge3\) and
\(T\ge1\), the number of actual singular cycles satisfying \(P=Q\) is
\[
 \ll_\varepsilon N^\varepsilon
       \left(\frac{T^3}{N^{1/3}}+\frac{T^6}{N}\right).           \tag{4}
\]
If all variables lie in \([T,2T)\) and \(N\asymp T^5\), this is
\(O_\varepsilon(T^{4/3+\varepsilon})\). The balance constants are fixed.
For a maximum box with all variables at most \(T\) and \(16T^4<N\),
(4) is \(O_\varepsilon(T^2N^\varepsilon)\).

**Theorem 2.** If \(N=p\) is prime and \(p>2T^3\), the right-hand side
of (4) can be replaced by
\[
 \ll_\varepsilon p^\varepsilon
       \left(\frac{T^3}{\sqrt p}+\frac{T^6}{p}\right).           \tag{5}
\]
If in addition \(p\equiv2\pmod3\), there are no such cycles.
Consequently the common dyadic sector at prime \(p\asymp T^5\)
has \(O_\varepsilon(T^{1+\varepsilon})\) cycles.

**Theorem 3.** Let \(R\ge1\). The number of actual singular cycles with
\[
       \frac P Q=\left(\frac r s\right)^3,\qquad
       r,s\in\mathbb Z_{>0},\quad \gcd(r,s)=1,\quad
       \max(r,s)\le R                                        \tag{6}
\]
is
\[
 \ll_\varepsilon N^\varepsilon
       \left(\frac{R^4T^3}{N^{1/3}}+
                         \frac{R^6T^6}{N}\right).             \tag{7}
\]
For prime \(p>2R^2T^3\), this improves to
\[
 \ll_\varepsilon p^\varepsilon
       \left(\frac{R^4T^3}{\sqrt p}+
                         \frac{R^6T^6}{p}\right),             \tag{8}
\]
and the sector is empty when \(p\equiv2\pmod3\).
In balanced boxes \(N\asymp T^5\), (7) gives
\(O_\varepsilon(T^{2+\varepsilon})\) whenever \(R\le T^{1/6}\).

The assertions are frequency bounds for actual cycles. They do not bound
the complementary singular sector, nonsingular cycles, or the full
six-cycle count.

## 2. Coherent units and a matching congruence

Connectedness of a cycle makes a rational labeling (2) unique up to one
common positive rational scalar. Clear denominators and divide their
global gcd. Since \(\gcd(a_{ij},c_{ij})=1\), each edge has
\[
 h_{ij}=\gcd(p_i,q_j),\qquad
 p_i=h_{ij}a_{ij},\quad q_j=h_{ij}c_{ij}.
\]
For a prime \(\ell\mid N\), every selected root and kernel is an
\(\ell\)-unit by (1). Thus the valuations of the two coherent labels at
every edge agree. Connectedness makes all six valuations equal, and
primitive normalization makes their common value zero. All \(p_i,q_j\)
are therefore units modulo \(N\). So are \(B,D,P,Q\).

Multiplying (1) by \(h_{ij}^2\) gives the exact edge equation
\[
                   b_i^3p_i^2+d_j^3q_j^2=N h_{ij}^2.           \tag{9}
\]
Multiply the three modular equations on the disjoint matching
\((b_i,d_i)\). The minus signs occur three times, giving
\[
                       B^3P^2+D^3Q^2\equiv0\pmod N.          \tag{10}
\]
For \(P=Q\), its square is a unit, so
\[
                               N\mid B^3+D^3.                \tag{11}
\]
This argument uses actual edge equations; an arbitrary choice of six
kernel vertices or six congruence witnesses is not substituted for them.

## 3. Primitive fractions in one residue class

**Lemma 4.** Fix a unit residue \(\lambda\) modulo \(N\), and \(H>0\).
The number of positive coprime pairs \((u,v)\), with \(u,v\le H\),
\(\gcd(v,N)=1\), and \(u/v\equiv\lambda\pmod N\), is
\[
                              O(1+H^2/N).                    \tag{12}
\]
The implied constant is absolute.

Distinct such primitive pairs are nonproportional positive vectors.
Their determinant is a nonzero multiple of \(N\). As each vector has
Euclidean norm at most \(\sqrt2 H\), the sine of the angular gap between
any two is at least \(N/(2H^2)\). If there are two vectors this ratio is
at most one; their angular gap is at least that same ratio, because
\(\sin x\le x\) for \(x\ge0\). Order the vectors in the interval of angles
\((0,\pi/2)\). Summing consecutive gaps gives at most
\(1+\pi H^2/N\) vectors. If there is at most one, (12) is immediate.

**Lemma 5.** For every unit residue \(\eta\) modulo \(N\), the number
of unit solutions to \(z^3\equiv\eta\pmod N\) is at most
\(3^{\omega(N)}\).

A nonempty fiber of cubing in the unit group is a translate of its
kernel. It is therefore enough to count roots of \(z^3=1\).
For a prime power \(\ell^e\) with \(\ell\ne3\), every unit root is simple
modulo \(\ell\): the derivative \(3z^2\) is a unit. Each root modulo
\(\ell\) has one lift at each successive power, and a polynomial of
degree three has at most three roots over the prime field. At
\(\ell=2\) there is only one unit root modulo two. At \(\ell=3\), a root
is \(1\) modulo three, and
\[
                v_3(z^3-1)=v_3(z-1)+1
                \quad\hbox{for }z\equiv1\pmod3,\ z\ne1.
\]
This follows by writing \(z=1+3u\) and factoring. Hence there are at most
three roots modulo \(3^e\), also for the two smallest exponents.
CRT proves the assertion. In particular \(3^{\omega(N)}\ll_\varepsilon
N^\varepsilon\).

We will also use the elementary divisor estimate
\(d_3(n)\ll_\varepsilon n^\varepsilon\), where \(d_3\) counts ordered
factorizations into three positive integers. Both subpower estimates
follow prime by prime: for sufficiently large primes the local factors
are bounded by the chosen positive prime power, and the finitely many
small primes are absorbed into a constant.

## 4. Counting the positive product pairs

Let \(\mathcal A_N(H)\) be the set of ordered positive unit pairs
\((L,M)\), with \(L,M\le H\), satisfying
\[
                               N\mid L^3+M^3.                \tag{13}
\]
Write \(g=\gcd(L,M)\), \(u=L/g\), \(v=M/g\). Since \(g\) is a unit
modulo \(N\), (13) implies \(u^3+v^3\equiv0\pmod N\).
The ratio \(u/v\) is one of the classes in Lemma 5. Positivity and
nonvanishing of the integer \(u^3+v^3\) imply
\[
            N\le u^3+v^3\le2(H/g)^3,\qquad
            g\le 2^{1/3}H/N^{1/3}.                            \tag{14}
\]
For each fixed \(g\), Lemmas 4 and 5 bound the number of pairs by
\[
                O\left(3^{\omega(N)}
                              (1+H^2/(Ng^2))\right).
\]
A primitive pair and \(g\) determine \(L,M\) uniquely. Summing over
positive integers in (14), and using the convergence of
\(\sum_{g\ge1}g^{-2}\), gives
\[
 \#\mathcal A_N(H)\ll 3^{\omega(N)}
                       \left(H/N^{1/3}+H^2/N\right).          \tag{15}
\]
If (14) has no positive integer solution, the set is empty, so no
additional constant term is required.

For (11), take \(H=T^3\). A fixed ordered pair \((B,D)\) has at most
\(d_3(B)d_3(D)\) ordered factorizations into the six kernels. Actual
equation (1) gives \(b_i^3<N\), \(d_j^3<N\), and therefore \(B,D<N\).
The divisor multiplicity is consequently \(O_\varepsilon(N^\varepsilon)\),
regardless of the size of the stated cutoff \(T\). For a fixed ordered
kernel sextuple there is at most one cycle with these positions in the
simple selected-edge graph. It either has all six required edges or
does not; their selected representations have already been fixed.
Thus (15) counts actual cycle frequency, with only a subpower
factorization multiplicity. Reallocate \(\varepsilon\) to obtain (4).

This last argument needs no codegree estimate, bound for fixed-label
cubic completions, or restriction \(16T^4<N\). Under that maximum-box
restriction, however,
\[
 T^3/N^{1/3}<16^{-1/3}T^{5/3},\qquad T^6/N<T^2/16,
\]
so the asserted \(T^2N^\varepsilon\) comparison follows for \(T\ge1\).
For a common dyadic box replace its cutoff by \(2T\) in (4).

## 5. The prime refinement

Suppose \(N=p\) is prime and \(p>2H\). For any primitive pair in the
preceding argument,
\[
                       0<u+v\le2H/g<p.
\]
The factorization
\[
                     u^3+v^3=(u+v)(u^2-uv+v^2)
\]
therefore forces \(p\mid u^2-uv+v^2\). This quadratic expression is
positive and at most \(\max(u,v)^2\), so
\[
                                  g\le H/\sqrt p.             \tag{16}
\]
The same fraction-packing sum, now over (16), gives
\[
 \#\mathcal A_p(H)\ll_\varepsilon p^\varepsilon
                                  (H/\sqrt p+H^2/p).          \tag{17}
\]
If \(p\equiv2\pmod3\), cubing is a bijection on
\(\mathbb F_p^\times\). The congruence \(u^3+v^3=0\) then forces
\(u+v=0\) modulo \(p\), contradicting the strict inequality above.
Apply (17) with \(H=T^3\), followed by the same kernel factorizations,
to prove Theorem 2. No global sign choice for composite moduli is used.

## 6. Bounded cubic product ratios

Under (6), the positive primitive pair \((r,s)\) is unique. Since
\(P,Q\) are integers, their ratio in lowest terms shows
\[
                  P=Jr^3,\qquad Q=Js^3
                  \quad\hbox{for some }J\in\mathbb Z_{>0}.
\]
All \(J,r,s\) are units modulo \(N\). Substitute in (10) and cancel
\(J^2\), obtaining
\[
                    N\mid (r^2B)^3+(s^2D)^3.                 \tag{18}
\]
For fixed \(r,s\), the products \(L=r^2B\), \(M=s^2D\) have height at
most \(R^2T^3\). The map \((B,D)\mapsto(L,M)\) is injective for this
fixed pair. Equation (15) therefore gives at most
\[
 \ll_\varepsilon N^\varepsilon
            \left(R^2T^3/N^{1/3}+R^4T^6/N\right)
\]
product pairs. Each pair again has subpower many kernel factorizations,
because the original \(B,D\) are less than \(N\). There are at most
\(R^2\) choices of the positive reduced center pair. This proves (7).
If \(p>2R^2T^3\), the prime argument applies uniformly to every center
and gives (8) and its emptiness assertion.

For a balanced common dyadic box, (7), with the kernel cutoff changed
to \(2T\), is
\[
       O_\varepsilon\left(
            T^\varepsilon(R^4T^{4/3}+R^6T)\right).
\]
Both terms are at most \(T^{2+\varepsilon}\) for \(R\le T^{1/6}\).
This is a counted restriction on the intrinsic label-product ratio,
not a bound on all coherent labels or all singular cycles.

## 7. Scope with positive coefficients

For an equation with fixed positive coefficients
\(\kappa_1a^2b^3+\kappa_2c^2d^3=N\) and coprime summands, the exact
matching congruence is
\[
              N\mid \kappa_1^3B^3P^2+\kappa_2^3D^3Q^2.        \tag{19}
\]
When (6) holds, the transformed products are
\(\kappa_1r^2B,\kappa_2s^2D\), not \(r^2B,s^2D\).
Thus the same calculation with \(C=\max(\kappa_1,\kappa_2)\) proves
the coefficient-dependent upper bound
\[
 \ll_\varepsilon N^\varepsilon
       \left(CR^4T^3/N^{1/3}+C^2R^6T^6/N\right).               \tag{20}
\]
All coefficients and transformed products are units modulo \(N\) on a
nonempty coprime graph. The original kernel products are still below
\(N\), so the factorization bound remains uniform. Formula (20) does
not imply the same \(T^2\) sector bound uniformly for arbitrary
coefficients. In particular, this argument does not discharge a full
coefficient-uniform maximum-box cycle hypothesis.

The following section also counts noncubic ratios with small reduced
height. Outside its range and the bounded cubic-center range, no
estimate beyond the general cycle bounds is asserted. These restrictions
are separate from the intrinsic minimum edge gcd and the two kernel-side
lcms. No implication from a general high-gcd cycle to either ratio
restriction has been established.


## 8. All small coherent product ratios

**Theorem 6.** Let \(X\ge1\). The number of actual singular cycles with
\[
             P/Q=r/s,\qquad r,s\in\mathbb Z_{>0},\quad
             \gcd(r,s)=1,\quad \max(r,s)\le X                  \tag{21}
\]
is
\[
 \ll_\varepsilon N^\varepsilon
       \left(\frac{T^3X^{8/3}}{N^{1/3}}+
                         \frac{T^6X^2}{N}\right).             \tag{22}
\]
No perfect-power assumption on \(P/Q\) is imposed. In a balanced common
dyadic box \(N\asymp T^5\), this is \(O_\varepsilon(T^{2+\varepsilon})\)
whenever \(X\le T^{1/4}\). The theorem concerns equation (1) with both
coefficients equal to one.

For a fixed reduced pair \((r,s)\) in (21), write \(P=Jr,Q=Js\).
As in Section 6, \(J,r,s\) are positive integers and units modulo \(N\).
The matching equation (10), after canceling \(J^2\), is
\[
                             N\mid r^2B^3+s^2D^3.             \tag{23}
\]
Keep the gcd of the original kernel products:
\[
                      g=\gcd(B,D),\qquad u=B/g,\quad v=D/g.
\]
It is a unit modulo \(N\), and \(u,v\) are positive coprime integers,
at most \(T^3/g\). Thus
\[
     (u/v)^3\equiv -s^2/r^2\pmod N,\qquad
     N\le r^2u^3+s^2v^3\le2X^2T^9/g^3.                       \tag{24}
\]
The integer in the middle is positive, so it is nonzero. Consequently
\[
                         g\le 2^{1/3}T^3X^{2/3}/N^{1/3}.      \tag{25}
\]
For fixed \(r,s,g\), Lemmas 4 and 5 give
\[
               O\left(3^{\omega(N)}
                               (1+T^6/(Ng^2))\right)
\]
primitive pairs. Summation over (25) gives, for this fixed center,
\[
 \ll_\varepsilon N^\varepsilon
           \left(T^3X^{2/3}/N^{1/3}+T^6/N\right).
\]
The original products \(B,D\) and their ordered kernel factorizations
are recovered just as in Section 4. The factorization cost remains
\(O_\varepsilon(N^\varepsilon)\), because \(B,D<N\). Summing over at
most \(X^2\) reduced pairs, and reallocating \(\varepsilon\), proves
(22). If the cutoff in (25) is less than one, there are no such products.

For a common dyadic kernel cutoff \(2T\), the same argument gives a
fixed constant times
\[
               N^\varepsilon\left(
                          T^{4/3}X^{8/3}+TX^2\right)
               \quad\hbox{when }N\asymp T^5.
\]
For \(X\le T^{1/4}\), the two terms are at most \(T^2\) and \(T^{3/2}\),
respectively. No codegree estimate, bound on coherent-label height,
or individual lower bound on selected roots is used.

The corresponding coefficient equation has
\(\kappa_1^3r^2B^3+\kappa_2^3s^2D^3\) in (23).
With \(C=\max(\kappa_1,\kappa_2)\), (25) acquires a factor \(C\).
Thus the calculation gives
\[
 \ll_\varepsilon N^\varepsilon
       \left(\frac{CT^3X^{8/3}}{N^{1/3}}+
                         \frac{T^6X^2}{N}\right),             \tag{26}
\]
with no assertion that its first term is coefficient-uniform in a
maximum-box application. This separates the unweighted actual
frequency theorem from a full coefficient-uniform cycle hypothesis.


## 9. A refined residual-cycle equivalence

In this section only, work in the common dyadic box
\([T,2T)\) for all four variables, with \(T\ge2\) and
\(c_0T^5\le N\le C_0T^5\), where \(c_0,C_0>0\) are fixed.
The cycle graph still uses one selected actual representation per edge.
For a singular cycle use the primitive labels of (2), and put
\[
                      g_*(C)=\min_{(i,j)\ {\rm edge}}
                                          \gcd(p_i,q_j).
\]
For a positive reduced rational number \(u/v\), let
\(H(u/v)=\max(u,v)\). If \(P/Q\) is a rational cube, its positive reduced
cube root \(r/s\) is unique.

Let \(\mathcal R_N(T)\) be the actual cycles for which both
\[
       \operatorname{lcm}(b_1,b_2,b_3)>8T^2,\qquad
       \operatorname{lcm}(d_1,d_2,d_3)>8T^2
\]
hold, and which are either nonsingular, or singular with
\(g_*(C)>T^{1/2}\). Define \(\mathcal R_N^{\rm pr}(T)\) by keeping a nonsingular cycle
in \(\mathcal R_N(T)\) only if \(U/V\) is a non-power rational with
reduced height greater than \(T^{5/2+\eta/4}\), using the absolute
\(\eta\) of [FixedReducedCycleRatio](ERDOS_943_SIX_CYCLES.md#fixed-reduced-ratios-and-proper-powers-for-every-target).
Here \(U/V\) is the alternating product of the six selected root
ratios in cycle order. Keep a singular cycle in the residual only if
\[
 \begin{split}
                    H(P/Q)&>T^{1/4},\\
 H(r/s)&>T^{1/6}
      \quad\hbox{whenever }P/Q=(r/s)^3 .
 \end{split}                                                   \tag{27}
\]
The second condition is vacuous when the product ratio is not a
rational cube.

**Proposition 7 (RefinedProductRatioResidualEquivalence).** The uniform
bound
\[
                       \#C_6(G)\ll_\varepsilon T^{2+\varepsilon}
                                                               \tag{28}
\]
for all such actual balanced graphs is equivalent to
\[
               \#\mathcal R_N^{\rm pr}(T)
                                  \ll_\varepsilon T^{2+\varepsilon}.
                                                               \tag{29}
\]
Both unrestricted statements remain open.

The forward implication is the subset inclusion. For the converse use
the established actual low-kernel-lcm and small-intrinsic-gcd sector
counts: cycles with either kernel-side lcm at most \(8T^2\), or singular
cycles with \(g_*\le T^{1/2}\), together contribute
\(O_\varepsilon(T^{2+\varepsilon})\).
These are the counts in LowKernelLcmSector and
SingularWeightedAndSmallGcdSector in
[Six-cycle counts and their residual domain](ERDOS_943_SIX_CYCLES.md).

Every remaining cycle is in \(\mathcal R_N(T)\).
FixedReducedCycleRatio bounds the excluded nonsingular cycles by
\(O_\varepsilon(T^{2+\varepsilon})\) for every \(N\).
An excluded singular cycle has either
\(H(P/Q)\le T^{1/4}\), or a reduced cubic center of height at most
\(T^{1/6}\). Theorem 6 and Theorem 3 bound these two sectors by
\(O_\varepsilon(T^{2+\varepsilon})\), respectively. In applying those
theorems the maximum kernel cutoff is \(2T\); this changes their
constants, not the displayed exponents or the thresholds in (27).
The union bound and (29) prove (28).

The old sector theorems are premises only for this equivalence; none
is used in Theorems 1--3 or 6. In particular the new frequency counts
do not acquire a codegree or fixed-label cubic-completion dependency.
Neither implication supplies the still missing bound (29).


## 10. A coefficient-uniform weighted product-ratio sector

Fix any positive integers \(\kappa_1,\kappa_2\), and now use the
simple selected-edge graph of coprime positive representations
\[
                    \kappa_1a^2b^3+\kappa_2c^2d^3=N.          \tag{30}
\]
The kernels are squarefree and at most \(T\). As before, no separate
root cutoff is needed. The constants below are independent of both
coefficients; the assertion is for each fixed coefficient pair, not
a count aggregated over different pairs.

For a singular cycle form \(P,Q,B,D\) as in (3). Define the positive
intrinsic weighted ratio
\[
                       \rho=
                    \frac{\kappa_1^3P^2}{\kappa_2^3Q^2}.       \tag{31}
\]
It is unchanged by a common rescaling of the coherent labels. Let
\(H(\rho)\) be its reduced positive rational height.

**Theorem 8 (WeightedProductRatioFrequency).** For every
\(\varepsilon>0\), uniformly in \(N\ge3\), \(T\ge1\), \(X\ge1\), and
all positive coefficient pairs, the number of actual singular cycles
in (30) with \(H(\rho)\le X\) is
\[
 \ll_\varepsilon N^\varepsilon
       \left(\frac{T^3X^{4/3}}{N^{1/3}}+
                         \frac{T^6X}{N}\right).               \tag{32}
\]
In balanced boxes \(N\asymp T^5\), the sector \(H(\rho)\le T^{1/2}\)
has \(O_\varepsilon(T^{2+\varepsilon})\) cycles. For maximum boxes,
the same target bound holds under either of the exact choices
\[
 \begin{array}{c|c}
      T=N^{11/50} & H(\rho)\le T^{17/44}\\
      T=N^{17/80} & H(\rho)\le T^{29/68}.
 \end{array}                                                   \tag{33}
\]
These are counted filtered sectors. The full coefficient-uniform
six-cycle bound remains open.

We first count the possible reduced centers without a coefficient
loss. Let \(K\) be the squarefree part of \(\kappa_1\kappa_2\), namely
the product of its prime factors with odd valuations. If
\(\rho=r/s\) is in lowest positive terms, its rational squareclass
is that of \(\kappa_1\kappa_2\). Equivalently,
\[
                              \operatorname{sf}(rs)=K.
\]
Indeed, (31) has squareclass \(\kappa_1^3/\kappa_2^3\), which is
\(\kappa_1\kappa_2\); multiplication by \(s^2\) identifies the
squareclass of \(r/s\) with that of \(rs\).

There are therefore unique squarefree coprime integers \(A,C\)
and positive integers \(u,v\) with
\[
                 r=A u^2,\qquad s=C v^2,\qquad AC=K.           \tag{34}
\]
For each of the \(2^{\omega(K)}\) ordered decompositions \(AC=K\),
the number with \(r,s\le X\) is at most
\[
             \left\lfloor\sqrt{X/A}\right\rfloor
             \left\lfloor\sqrt{X/C}\right\rfloor
                              \le X/\sqrt K .
\]
Ignoring any further coprimality restrictions only increases this
upper bound. Furthermore
\[
 \frac{2^{\omega(K)}}{\sqrt K}
    =\prod_{\ell\mid K}\frac2{\sqrt\ell}
    \le\frac4{\sqrt6}<2.                                     \tag{35}
\]
Only the primes two and three give factors greater than one; every
prime at least five gives a factor less than one. The empty product
is one. Thus there are fewer than \(2X\) possible reduced centers
\(r/s\) of height at most \(X\), uniformly in the coefficient sizes.

Coprimality of the two summands in (30) makes every coefficient,
root and kernel a unit modulo \(N\) on a nonempty graph. The
connected-label argument in Section 2 again makes every primitive
coherent label a unit. On each selected edge,
\[
       \kappa_1b_i^3p_i^2+\kappa_2d_j^3q_j^2
                                      =N\gcd(p_i,q_j)^2.
\]
Multiplication along the matching gives the weighted congruence (19).
For the reduced ratio \(r/s\) in (31), write
\[
        \kappa_1^3P^2=Jr,\qquad \kappa_2^3Q^2=Js .
\]
Here \(J,r,s\) are positive integers and units modulo \(N\).
Canceling \(J\) yields
\[
                              N\mid rB^3+sD^3.                \tag{36}
\]
Keep \(g=\gcd(B,D)\), and write \(a=B/g\), \(b=D/g\).
All are positive units modulo \(N\), and \(\gcd(a,b)=1\).
For a fixed center with \(r,s\le X\),
\[
 \begin{split}
               (a/b)^3&\equiv-s/r\pmod N,\\
               N\le ra^3+sb^3&\le2XT^9/g^3 .
 \end{split}
\]
Positivity is essential in the second line. It gives
\[
                       g\le 2^{1/3}T^3X^{1/3}/N^{1/3}.        \tag{37}
\]
The cubic fiber has at most \(3^{\omega(N)}\) unit classes by
Lemma 5. Lemma 4, applied at height \(T^3/g\), counts at most
\[
               O\left(3^{\omega(N)}
                                 (1+T^6/(Ng^2))\right)
\]
primitive pairs for each \(g\). Summing over (37), with the same
empty-cutoff convention as before, bounds the product pairs for
this fixed center by
\[
 \ll_\varepsilon N^\varepsilon
          \left(T^3X^{1/3}/N^{1/3}+T^6/N\right).
\]
There are fewer than \(2X\) centers by (34)--(35). Each recovered
\(B,D\) has subpower many ordered kernel factorizations: (30)
implies \(\kappa_1B<N\), \(\kappa_2D<N\), so in particular
\(B,D<N\). A fixed kernel sextuple gives at most one cycle in the
fixed selected-edge graph. Summing the center bounds and reallocating
\(\varepsilon\) proves (32). No root-height estimate, codegree
bound, or coefficient-dependent divisor constant is used.

For the balanced comparison, the two terms in (32), before the
subpower factor, become at most
\(T^{4/3}X^{4/3}\) and \(TX\), up to fixed balance constants.
They are at most \(T^2\) and \(T^{3/2}\) when \(X\le T^{1/2}\).
A common dyadic kernel cutoff \(2T\) only changes constants.

For the first choice in (33), the two exponents of \(T\) are
\[
 3-\frac{50}{33}+\frac43\frac{17}{44}=2,\qquad
 6-\frac{50}{11}+\frac{17}{44}=\frac{81}{44}<2.
\]
For the second choice they are
\[
 3-\frac{80}{51}+\frac43\frac{29}{68}=2,\qquad
 6-\frac{80}{17}+\frac{29}{68}=\frac{117}{68}<2.
\]
Since these are fixed powers of \(N\), the \(N^\varepsilon\)
factor can be expressed as \(T^\varepsilon\) after reallocation.
More generally, at \(T=N^\alpha\), \(0<\alpha<1/4\), (32) gives
the target \(T^{2+\varepsilon}\) whenever
\[
 H(\rho)\le T^\eta,\qquad
 0\le\eta\le
       \min\left\{\frac{1/\alpha-3}{4},\,
                                \frac1\alpha-4\right\}.        \tag{38}
\]

When \(\kappa_1=\kappa_2=1\), one has
\(H(\rho)=H(P/Q)^2\). Substitution in (32) recovers Theorem 6:
there is no additional unweighted saving hidden in this change
of parameter. The uniformity in arbitrary positive coefficients
comes from the fixed-squareclass center count, not from removing
the coefficient factor in (26) while keeping its unweighted
height parameter. No estimate for the complementary weighted-ratio
sector, or a full cycle-frequency estimate, follows here.
