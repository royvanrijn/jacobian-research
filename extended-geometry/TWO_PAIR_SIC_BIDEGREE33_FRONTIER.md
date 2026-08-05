# The two-pair bidegree-\((3,3)\) frontier

## 1. Status

This note treats the first genuinely new two-pair bidegree after the
complete bidegree-\((2,2)\) theorem.  Unequal bidegrees are already
one-sided for the total dual-minus-coordinate grading, so the next balanced
case is \((3,3)\).

The moment--nullcone question in this degree is now settled negatively.
The explicit six-entry form in the
[Rodrigues-survivor theorem](TWO_PAIR_SIC_BIDEGREE33_RODRIGUES_SURVIVOR.md)
has every pure contraction moment zero and has invertible coefficient
matrix.  It is nevertheless SIC-safe: for every fixed multiplier \(Q\),
\(\mathcal E_2(QF^m)=0\) once
\(m>3\deg_{Z,Y}Q\).  Thus \(\mathrm{MN}_3\) is false, but this does not
produce a bidegree-\((3,3)\) SIC counterexample.

The ten earlier exact characteristic-zero results are:

1. the full pair-linear one-sided nullcone in the sixteen-dimensional
   coefficient space has an exact seven-parameter elimination, dimension
   seven, projective degree \(20\), and a Gröbner basis of size \(148\);
2. the first thirteen full moments are algebraically independent, proved
   by a nonzero exact \(13\times13\) Jacobian minor;
3. nevertheless, invariant-ring Hilbert coefficients prove that degrees
   \(1,\ldots,13\) cannot be the degrees of a homogeneous system of
   parameters; replacing \(13\) by \(14\) gives an algebraically
   independent, Hilbert-compatible minimal candidate;
4. on each pure irreducible summand, the moments cut out exactly the
   one-sided nullcone: orders \(2,4,6,10\) for binary sextics, orders
   \(2,3\) for binary quartics, and order \(2\) for binary quadratics.
5. on the mixed
   \(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) branch, a non-null
   quadratic component is excluded by moments through order six.
6. after normalizing a non-null quadratic, moments \(2,\ldots,12\) have
   exact full Jacobian rank eleven on every residual-torus chart orbit in
   the complete
   \(\operatorname{Sym}^6\oplus\operatorname{Sym}^4\oplus
   \operatorname{Sym}^2\) branch.
7. on those chart orbits, \(\mu_2\) eliminates one opposite-weight
   variable with a constant nonzero pivot; on the first chart, the
   reduced \(\mu_3\) gives two explicit principal-open eliminations and
   one common boundary.
8. a natural two-parameter plane in that common boundary is excluded
   exactly by moments \(3\) and \(4\), with a displayed unit certificate.
9. the moment map \((\mu_2,\ldots,\mu_{12})\) has maximal differential
   rank on all three \(s_0\)-chart pivot strata: ranks \(11,10,9\) on
   \(A\ne0\), \(A=0,B\ne0\), and \(A=B=0\), respectively.
10. a four-parameter family in the common boundary \(A=B=\mu_2=0\),
    strictly containing the earlier two-parameter plane, is excluded
    exactly by moments \(3,\ldots,7\).

Fifteen further exact results sharpen the remaining SIC problem:

11. on the normalized anti-Weyl locus, all odd moments vanish by symmetry
    and the even moments through order fourteen generate the unit ideal
    over \(\mathbb Q\);
12. every mixed coefficient support of size at most five is classified
    exactly, and all \(11{,}200\) seven-entry coefficient tori are excluded
    over \(\mathbb Q\) by moments through order twelve;
13. among the \(7{,}588\) six-entry coefficient tori, exactly two
    normalized complex points survive, as proved by degree-one RURs; they
    are Weyl/torus copies of the full-rank Rodrigues form;
14. the all-order Rodrigues and integration-by-parts identities prove that
    orbit SIC-safe with the explicit multiplier bound above; and
15. on its normalized null-quadratic \(s_6\)-chart, it is an isolated
    length-five local component through \(\mu_{11}\), while its natural
    five-variable tangent slice has no other reduced point.
16. the rank-two nine-moment fiber has an explicitly isolated semistable
    real point on the anti-Weyl chart; its projective quotient point is
    reduced, and it is excluded from both corrected systems already by
    the certified sign \(\mu_{10}>0\).
17. all \(12{,}780\) mixed eight-entry coefficient tori are classified:
    \(12{,}765\) are unit ideals through order twelve, the unique timed-out
    odd-parity system is a unit ideal through order fourteen, and every
    complex point on the other fourteen systems is explicitly one-sided.
18. on the generic rank-two Hurwitz chart with \(B(0)=1\),
    \(\Delta M_{01}\ne0\), exact characteristic-zero elimination excludes
    the complete \(\lambda=0\) fibre and the doubly exceptional
    \(\mu_2\)-pivot branch \(P_1=P_2=0\), both already through
    \(\mu_8\).
19. all sixteen nine-entry \(3\times3\) rectangle supports are SIC-safe.
    Six transpose/reversal orbits cover them.  On each dense coefficient
    torus the exact finite scheme through \(\mu_{14}\) has every
    \(2\times2\) minor in its radical, as proved by 54 characteristic-zero
    Rabinowitsch unit ideals; hence its reduced points have rank one.
    Every coordinate boundary has support at most eight and is covered by
    result 17.
20. all 96 complete-two-row/column fringe supports are SIC-safe.
    They form 24 transpose/reversal orbits; on every representative
    the dense coefficient-torus ideal through \(\mu_{14}\) is the unit
    ideal over \(\mathbb Q\), while every boundary has support at most
    eight.
21. the complete exact-rank-two reversal-parity factor family is SIC-safe.
    The reversal centralizer moves every nonzero even/odd \(U\)-pair into
    the normalized chart \(U_{0,*}=(1,1)\).  On both projective
    semistable \(W\)-charts the exact-rank-two opens are unit
    ideals through \(\mu_6\).  The invariant-zero slice has exactly two
    one-dimensional components through \(\mu_6\); both are fixed-flag
    one-sided families over \(\mathbb Q(q)\), with recurrence
    \(\nu_{m+1}=0\), a nonzero degree-two mixed value, and mixed tail
    \(2m>e\).
22. all 576 cross-plus-two supports are SIC-safe.  They consist of one
    complete row, one complete column, and any two further entries, and
    form 156 transpose/reversal orbits.  Every dense coefficient-torus
    ideal is the unit ideal over \(\mathbb Q\) through \(\mu_{10}\).
23. all 480 regular three-line supports are SIC-safe.  They have three
    occupied rows with three entries each, or are transposes of such a
    support, with the sixteen rectangles removed.  Of their 120 symmetry
    orbits, 114 have unit dense-torus ideals through \(\mu_{10}\).  Each
    of the other six has one rational point; all six have coefficient
    rank two and an exact flag change into the fixed chamber \(i>j\), so
    their recurrence is \(\nu_{m+1}=0\) and their mixed cutoff is \(m>e\).
24. all 1,148 three-line supports with line counts \(4+3+2\), or their
    transposes, are SIC-safe.  Their 287 symmetry representatives have
    exact dense coefficient-torus unit ideals through \(\mu_{10}\).
    Together with the rectangle, fringe, and regular-three-line classes,
    this closes every nine-entry support with an empty row or column.
25. all 1,244 full-line supports with line-count partition
    \(4+3+1+1\), no complete line on the other axis, or their transposes,
    are SIC-safe.  Every one of their 311 dense symmetry representatives
    has unit moment ideal over \(\mathbb Q\) through \(\mu_{10}\).

The fourth result excludes every SIC(2) counterexample lying in a single
irreducible summand.  The sixth gives dimension-sized moment coordinates,
not a zero-fiber theorem.  Results 11--25 remove the first semistable
all-order survivor from the SIC search and force any counterexample to
have at least nine nonzero entries in the displayed coefficient basis.
They do not classify the remaining full mixed
\(\operatorname{Sym}^6\oplus\operatorname{Sym}^4\oplus
\operatorname{Sym}^2\) locus.  Consequently bidegree-\((3,3)\) SIC safety,
and hence minimality of the known bidegree-\((4,4)\) counterexample,
remains open.

Use contraction pairs \((W,Z),(V,Y)\) and the basis

\[
 M_{ij}=W^{3-i}V^iZ^{3-j}Y^j,\qquad 0\le i,j\le3.          \tag{1.1}
\]

## 2. The full one-sided nullcone

For a fixed flag, the positive-weight positions are the six pairs \(i>j\).
On the finite flag chart use

\[
 W'=W,\qquad V'=V-qW,\qquad Z'=Z+qY,\qquad Y'=Y.           \tag{2.1}
\]

This preserves the contraction pairing:

\[
 W'Z'+V'Y'=WZ+VY.
\]

The general one-sided form on the chart is

\[
 f_0=\sum_{i>j}a_{ij}
 (W')^{3-i}(V')^i(Z')^{3-j}(Y')^j.                        \tag{2.2}
\]

There are six fiber parameters \(a_{ij}\) and one flag parameter \(q\).
Eliminating these seven parameters from the sixteen expanded coefficients
over \(\mathbb Q\) gives a prime coefficient ideal \(J_{33}\).  Exact
Singular computation gives

\[
 \dim V(J_{33})=7,\qquad
 \deg \mathbb P V(J_{33})=20,                              \tag{2.3}
\]

and the reduced degree-compatible Gröbner basis has \(148\) elements.
The finite \(q\)-chart is dense in the projective flag line, so its image
closure is the complete pair-linear one-sided nullcone.

This calculation supplies the exact target ideal for the full
bidegree-\((3,3)\) moment comparison.  It does not assert that the moment
zero set already equals this nullcone.

### 2.1 The first minimal set fails, and a corrected candidate

Write the general coefficient matrix as \(C=(c_{ij})_{0\leq i,j\leq3}\)
and introduce

\[
 F_C(x,y)=\sum_{i,j=0}^3c_{ij}x^iy^j.
\]

The full scalar moments can be evaluated from the diagonal of \(F_C^m\):

\[
 \mu_m(C)=
 \sum_{I=0}^{3m}(3m-I)!\,I!\,
 [x^Iy^I]F_C(x,y)^m.                                     \tag{2.4}
\]

Consequently

\[
 \frac{\partial\mu_m}{\partial c_{ab}}(C)=
 m\sum_{I=0}^{3m}(3m-I)!\,I!\,
 [x^{I-a}y^{I-b}]F_C(x,y)^{m-1}.                         \tag{2.5}
\]

At the integral point

\[
 C_*=
 \begin{pmatrix}
 2&2&-4&0\\
 4&3&2&0\\
 3&1&-1&4\\
 -2&0&-2&-3
 \end{pmatrix},                                          \tag{2.6}
\]

take the Jacobian of \(\mu_1,\ldots,\mu_{13}\) and the thirteen columns

\[
 c_{03},c_{10},c_{11},c_{12},c_{13},c_{20},c_{21},
 c_{22},c_{23},c_{30},c_{31},c_{32},c_{33}.
\]

Exact integer evaluation of (2.5) gives the nonzero determinant

\[
\begin{aligned}
 -{}&2^{256}3^{107}5^{48}7^{29}11^{13}13^7 17^3
       \cdot19\cdot139\cdot4493\cdot886069\\
    &{}\cdot651443921434147\\
    &{}\cdot108355984865758174686774716693198468303195999878781902931.
\end{aligned}                                             \tag{2.7}
\]

> **Proposition 2.1.** The first thirteen full moments
> \(\mu_1,\ldots,\mu_{13}\) are algebraically independent over
> \(\mathbb Q\).

Indeed, (2.7) gives Jacobian rank thirteen at \(C_*\).  Since the generic
\(\mathrm{SL}_2\)-stabilizer on this sixteen-dimensional representation is
finite, the invariant quotient has dimension \(16-3=13\).  Thus these
moments form a transcendence basis of its invariant function field and
attain the Krull-height lower bound for any moment set defining the
nullcone.

They are, however, **not** a homogeneous system of parameters. This can
be decided before attempting the zero fiber. The
\(\mathrm{SL}_2\)-weights of

\[
 V_3=\operatorname{Sym}^0\oplus\operatorname{Sym}^2
          \oplus\operatorname{Sym}^4\oplus\operatorname{Sym}^6
\]

are

\[
 \mathcal W=\{0^4,2^3,(-2)^3,4^2,(-4)^2,6,-6\}.
 \tag{2.7a}
\]

If \(H(t)=\sum_{n\geq0}h_nt^n\) is the Hilbert series of
\(\mathbb Q[V_3]^{\mathrm{SL}_2}\), the elementary
\(\mathrm{SL}_2\) weight-multiplicity formula gives

\[
 h_n=[t^nq^0]\prod_{w\in\mathcal W}(1-tq^w)^{-1}
     -[t^nq^2]\prod_{w\in\mathcal W}(1-tq^w)^{-1}.        \tag{2.7b}
\]

Exact expansion begins

\[
 1,1,4,8,26,53,146,305,704,1417,2920,5533,10500,18825.
\]

Now form the numerator required by parameter degrees \(1,\ldots,13\):

\[
 P_{13}(t)=H(t)\prod_{m=1}^{13}(1-t^m).
\]

The exact coefficient calculation gives

\[
 \boxed{[t^{63}]P_{13}(t)=-2186}.                         \tag{2.7c}
\]

The invariant ring is Cohen--Macaulay in characteristic zero. Therefore
for any homogeneous system of parameters the corresponding product is
the Hilbert series of the Artinian quotient and has nonnegative
coefficients. Equation (2.7c) proves:

> **Proposition 2.2.** No homogeneous system of parameters of
> \(\mathbb Q[V_3]^{\mathrm{SL}_2}\) has degrees
> \(1,2,\ldots,13\). In particular,
> \(\mu_1,\ldots,\mu_{13}\) do not define the nullcone.

Since these thirteen moments are algebraically independent, this is
stronger than a warning about an unfinished computation: their common
zero set in the invariant quotient necessarily contains a non-origin
point. Equivalently, the first thirteen moment equations have a
semistable component outside the one-sided nullcone.

The least-total-degree repair permitted by the same Hilbert test is

\[
 \boxed{\mu_1,\mu_2,\ldots,\mu_{12},\mu_{14}}.             \tag{2.7d}
\]

At the same integral point \(C_*\), the \(13\times13\) minor on columns

\[
 c_{00},c_{01},c_{02},c_{03},c_{10},c_{11},c_{12},
 c_{13},c_{20},c_{21},c_{22},c_{23},c_{30}
\]

is a nonzero \(303\)-digit integer. Hence:

> **Proposition 2.3.** The moments
> \(\mu_1,\ldots,\mu_{12},\mu_{14}\) are algebraically independent.

For these corrected degrees, the exact Hilbert numerator computation
through degree \(100\) is nonnegative, has last observed nonzero
coefficient \(1\) in degree \(76\), and vanishes in degrees
\(77,\ldots,100\). Its checked coefficient sum is \(9\,226\,602\).
This is strong compatibility evidence, not yet a proof that (2.7d) is a
homogeneous system of parameters. The corrected zero-fiber equality is
the next geometric target.

The completed invariant algebra gives a smaller corrected candidate.
In the conventions below, its first missing quadratic Casimir is
\[
 q_2=80\Delta_2.
\tag{2.7e}
\]
The system
\[
 \boxed{\mu_1,\ldots,\mu_{12},q_2}
\tag{2.7f}
\]
also has exact Jacobian rank thirteen.  Its Hilbert numerator is
nonnegative through degree \(120\), has predicted top degree \(64\), no
observed tail, and coefficient sum \(1\,318\,086\), compared with
\(9\,226\,602\) for (2.7d).  Two degree-sum-\(92\) mixed
moment/Casimir controls likewise pass the same necessary test and have
full rank.

This does not yet identify the unavoidable semistable component of
\((\mu_1,\ldots,\mu_{13})\).  It does show that (2.7f) removes every
point with \(\Delta_2\ne0\), leaving only the normalized null-quadratic
stratum \(F_2=L^2\ne0\) and the boundary \(F_2=0\).  A complete
weight-\(14\) evaluation test additionally proves that
\(\mu_{14}\) is independent, modulo the lower-moment span, from every
degree-\(14\) polynomial in \(q_2,q_4\) inside the generated completed
algebra.  Thus \(\mu_{14}\) is not simply a polynomial Casimir
replacement, although equality of their zero divisors in the full
invariant quotient remains open.  The exact systems, ranks, and
nonrelation certificate are canonical in
[`COMPLETED_MOMENT_ALGEBRA_RESEARCH.md`](COMPLETED_MOMENT_ALGEBRA_RESEARCH.md#51-the-cubic-corrected-momentcasimir-comparison).

The first synchronization step is now generically established on the
nonzero null-quadratic chart.  After \(F_2=X^2\), the exact linear normal
symbols of \(\mu_2,\mu_3,\mu_4\) have generic rank three.  Their nonzero
maximal minors share one explicit irreducible cubic pivot \(P(s_2,t_1)\).
At the integral allowed base \((20,27,36,47,60)\), the full seven-normal
fiber of \(\mu_2,\ldots,\mu_{12}\) is zero-dimensional modulo \(32003\),
of length \(195\); good reduction therefore proves characteristic-zero
transverse isolation on a nonempty open subset.  This replaces the
generic part of the proposed synchronization calculation by a
seven-variable normal fiber.  The residual rank locus splits exactly
into two three-dimensional quadratic-field components and one
two-dimensional rational locus, all disjoint from \(P=0\).  Full
seven-normal fibers at good reductions of exact algebraic points on
\(P=0\) and both quadratic components again have dimension zero and
length \(195\); the lower rational locus has dimension zero and length
\(197\).  Thus transverse isolation holds on a nonempty open subset of
every linear-rank stratum.  Proper closed subsets inside those strata
and the separate \(F_2=0\) chart remain open; see
[`COMPLETED_MOMENT_ALGEBRA_RESEARCH.md`](COMPLETED_MOMENT_ALGEBRA_RESEARCH.md#52-generic-synchronization-on-the-null-quadratic-chart).

### 2.2 The global quadratic anchor and the torus-fixed test

The complete Clebsch--Gordan change of basis can be inverted exactly. In
the divided-power convention

\[
 F_2(X,T)=r_0X^2+2r_1XT+r_2T^2,
\]

the quadratic coordinates of the general matrix \(C=(c_{ij})\) are

\[
\begin{aligned}
 r_0&=\frac{3c_{10}+2c_{21}+3c_{32}}{10},\\
 r_1&=\frac{-9c_{00}-c_{11}+c_{22}+9c_{33}}{20},\\
 r_2&=-\frac{3c_{01}+2c_{12}+3c_{23}}{10}.
\end{aligned}                                             \tag{2.8}
\]

Thus the global quadratic-anchor invariant is explicitly

\[
 \boxed{\Delta_2=r_1^2-r_0r_2.}                           \tag{2.9}
\]

This is the invariant that must be placed in the radical of the full
moment ideal.

There is a first exact zero-fiber test. Restrict to the maximal-torus fixed
subspace

\[
 C=\operatorname{diag}(c_{00},c_{11},c_{22},c_{33}).       \tag{2.10}
\]

Equivalently, put

\[
 H(q)=c_{00}+c_{11}q+c_{22}q^2+c_{33}q^3.
\]

Formula (2.4) becomes

\[
 \mu_m=\sum_{I=0}^{3m}(3m-I)!\,I!\,[q^I]H(q)^m.           \tag{2.11}
\]

Let \(I_T=(\mu_1,\mu_2,\mu_3,\mu_4)\) in the four diagonal
coefficients. Exact Gröbner reduction gives

\[
 c_{00}^7,c_{11}^7,c_{22}^7,c_{33}^7\in I_T.             \tag{2.12}
\]

Hence

\[
 \sqrt{I_T}=(c_{00},c_{11},c_{22},c_{33}).                \tag{2.13}
\]

> **Proposition 2.2.** The first four moments have no nonzero common zero
> on the torus-fixed diagonal slice.

Every nonzero diagonal endomorphism is semistable: if its orbit closure
contained zero, all characteristic-polynomial invariants would vanish, so
the diagonal endomorphism would be nilpotent and hence zero. Proposition
2.2 therefore excludes the entire maximal-torus fixed semistable slice,
including its \(\Delta_2\ne0\) locus.

This is a meaningful anchor test but not the global anchor theorem.
Generic closed \(\mathrm{SL}_2\)-orbits have finite stabilizer and need not
meet the torus-fixed subspace.

## 3. The binary-sextic summand

Let

\[
 D=Y\partial_Z-W\partial_V.                                \tag{3.1}
\]

Starting from the highest-weight vector \(V^3Z^3\), the seven vectors

\[
 e_k=\frac{D^k}{k!}(V^3Z^3),\qquad 0\le k\le6,             \tag{3.2}
\]

span the highest irreducible summand \(\operatorname{Sym}^6\).  Write

\[
 f_s=\sum_{k=0}^6s_ke_k.
\]

In the row-major basis (1.1), its coefficient matrix is

\[
 \begin{pmatrix}
 -s_3&-3s_4&-3s_5&-s_6\\
 3s_2&9s_3&9s_4&3s_5\\
 -3s_1&-9s_2&-9s_3&-3s_4\\
 s_0&3s_1&3s_2&s_3
 \end{pmatrix}.                                           \tag{3.3}
\]

Associate to \(s\) the binary sextic

\[
 H_s(X,T)=\sum_{k=0}^6\binom6k s_kX^{6-k}T^k.             \tag{3.4}
\]

The nullcone in this summand consists exactly of sextics with a root of
multiplicity at least four:

\[
 H_s=L^4Q,                                                 \tag{3.5}
\]

where \(L\) is linear and \(Q\) is quadratic.  On the finite root chart,

\[
 H_s=(X+qT)^4(uX^2+vXT+wT^2).                             \tag{3.6}
\]

Eliminating \(u,v,w,q\) gives a prime ideal

\[
 J_6\subset\mathbb Q[s_0,\ldots,s_6]
\]

with ten Gröbner generators, affine dimension four, and projective degree
twelve.

## 4. Exact moment radical

Put

\[
 \mu_m(s)=\mathcal E_2(f_s^m)
\]

and form the ideal

\[
 I_6=(\mu_2,\mu_4,\mu_6,\mu_{10})
 \subset\mathbb Q[s_0,\ldots,s_6].                        \tag{4.1}
\]

The checker constructs every \(\mu_m\) directly from the multinomial
contraction formula, without using the nullcone parametrization.  Exact
Gröbner reduction gives

\[
 \dim V(I_6)=4,\qquad |G(I_6)|=65,\qquad
 I_6\subseteq J_6.                                        \tag{4.2}
\]

For the ten reduced Gröbner generators \(j_1,\ldots,j_{10}\) of \(J_6\),
the exact power certificate is

\[
 j_1\in I_6,\qquad j_r^5\in I_6\quad(2\le r\le10).         \tag{4.3}
\]

Therefore

\[
 \sqrt{I_6}=J_6.                                          \tag{4.4}
\]

> **Theorem 4.1.** Let \(f\) lie in the binary-sextic
> \(\operatorname{Sym}^6\) summand of the two-pair bidegree-\((3,3)\)
> coefficient space over a characteristic-zero field.  If all pure
> contractions \(\mathcal E_2(f^m)\) vanish, then after scalar extension to
> an algebraic closure, \(f\) is pair-linearly one-sided.  Consequently
> \(\mathcal E_2(gf^m)=0\) for every fixed \(g\) and all sufficiently large
> \(m\).  This summand contains no SIC(2) counterexample.

Only the four displayed moment orders are needed for the implication.

### The lower pure summands

The same lowering-chain construction identifies the
\(\operatorname{Sym}^4\) summand with binary quartics

\[
 K_t(X,T)=\sum_{k=0}^4\binom4k t_kX^{4-k}T^k.
\]

Its nullcone is \(K_t=L^3R\).  Exact contraction gives

\[
\begin{aligned}
\mu_2&=t_0t_4-4t_1t_3+3t_2^2,\\
\mu_3&=-t_0t_2t_4+t_0t_3^2+t_1^2t_4
       -2t_1t_2t_3+t_2^3.
\end{aligned}                                             \tag{4.5}
\]

Elimination of \((X+qT)^3(uX+vT)\) proves directly that
\((\mu_2,\mu_3)\) is the prime quartic nullcone ideal.  It has affine
dimension three and projective degree six.

Likewise, on the \(\operatorname{Sym}^2\) summand write

\[
 R_r(X,T)=r_0X^2+2r_1XT+r_2T^2.
\]

The second moment is, up to a nonzero scalar,

\[
 \mu_2=-r_0r_2+r_1^2,                                     \tag{4.6}
\]

which generates the prime \(L^2\) nullcone ideal.  It has affine dimension
two and projective degree two.

> **Corollary 4.2.** None of the three pure irreducible summands
> \(\operatorname{Sym}^6\), \(\operatorname{Sym}^4\), or
> \(\operatorname{Sym}^2\) contains an SIC(2) counterexample.

## 5. The normalized non-null quadratic branch

Under the diagonal \(\mathrm{SL}_2\)-action,

\[
 \operatorname{End}(\operatorname{Sym}^3)
 \cong
 \operatorname{Sym}^6\oplus
 \operatorname{Sym}^4\oplus
 \operatorname{Sym}^2\oplus
 \operatorname{Sym}^0.                                   \tag{5.1}
\]

This is specifically a decomposition for the contraction-preserving
diagonal group.  Before restriction, the coefficient space is the single
external-tensor representation
\[
 \operatorname{Sym}^3(k^2)^*\boxtimes\operatorname{Sym}^3(k^2)
\]
of \(\mathrm{SL}_2\times\mathrm{SL}_2\).  Independent changes in the two
factors do not preserve \(WZ+VY\), hence do not preserve the SIC moments.
They may organize left/right binary-cubic rank or root data, but they
cannot be used as the orbit quotient for (5.1), the quadratic
null/non-null split, or the moment equations.

The first moment removes the scalar summand.  The unresolved locus is
therefore the mixed fifteen-dimensional representation

\[
 \operatorname{Sym}^6\oplus
 \operatorname{Sym}^4\oplus
 \operatorname{Sym}^2.                                   \tag{5.2}
\]

A non-null binary quadratic can be moved by \(\mathrm{SL}_2\) to

\[
 R=2cXT,\qquad c\ne0.                                     \tag{5.3}
\]

Keeping \(c\) as a variable preserves homogeneity.  On the
\(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) branch, let

\[
 I_{42}=(\mu_2,\mu_3,\mu_4,\mu_5,\mu_6).
\]

Exact characteristic-zero Gröbner reduction gives

\[
 c^6\in I_{42},\qquad c^5\notin I_{42}.                   \tag{5.4}
\]

> **Theorem 5.1.** The mixed
> \(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) branch with non-null
> quadratic component contains no SIC(2) counterexample.

Indeed, (5.4) contradicts \(c\ne0\) whenever the five displayed moments
vanish.

The analogous \(\operatorname{Sym}^6\oplus\operatorname{Sym}^2\) branch
has also been sharply reduced, but not yet promoted to characteristic
zero.  Over \(\mathbb F_{32003}\), the homogeneous ideal generated by the
even moments

\[
 \mu_2,\mu_4,\mu_6,\mu_8,\mu_{10},\mu_{12},\mu_{14}
\]

has dimension four and a Gröbner basis of size \(7576\), with

\[
 c^{25}\in I,\qquad c^{24}\notin I.                       \tag{5.5}
\]

This is finite-field computation, not a theorem over \(\mathbb Q\).
Ordinary rational Gröbner reduction and exact modular reconstruction did
not finish within the recorded bounded runs.  Equation (5.5) identifies a
specific characteristic-zero membership certificate to reconstruct; it
does not justify silently transferring the modular result.

### 5.1 The five full mixed anchor charts

Set \(c=1\) in (5.3), which is legitimate on the moment-zero locus by
overall scalar homogeneity.  The stabilizer of \(XT\) contains the diagonal
torus.  In the lowering-chain coordinates its weights are

\[
\begin{array}{c|rrrrrrr}
 &s_0&s_1&s_2&s_3&s_4&s_5&s_6\\ \hline
\operatorname{wt}&3&2&1&0&-1&-2&-3
\end{array},
\qquad
\begin{array}{c|rrrrr}
 &t_0&t_1&t_2&t_3&t_4\\ \hline
\operatorname{wt}&2&1&0&-1&-2.
\end{array}                                               \tag{5.6}
\]

If every nonzero-weight coordinate vanishes, only \(s_3,t_2\), and the
normalized quadratic remain; this is contained in the torus-fixed diagonal
slice closed by Proposition 2.2.  Otherwise, after scalar extension to an
algebraic closure, the residual torus sets one nonzero-weight coordinate
to one.  Weyl reflection exchanges positive and negative weights, leaving
the five representative charts

\[
 s_0=1,\quad s_1=1,\quad s_2=1,\quad t_0=1,\quad t_1=1.    \tag{5.7}
\]

Each chart has eleven free coordinates.  Restrict the eleven moments

\[
 \mu_2,\mu_3,\ldots,\mu_{12}.                              \tag{5.8}
\]

At one displayed integral point on each chart, exact evaluation of the
eleven-by-eleven Jacobian gives a nonzero integer determinant.

> **Proposition 5.2.** On every chart orbit in (5.7), the restricted
> moments (5.8) are algebraically independent.

The exact points and determinants are stored in
[`two_pair_sic_bidegree33_anchor_jacobians.json`](../artifacts/generated-results/two_pair_sic_bidegree33_anchor_jacobians.json).
This proves that (5.8) is dimension-sized on every full mixed anchor chart.
It does **not** prove that its zero fiber is empty.  The remaining full
anchor problem is now a finite collection of eleven-equation affine
zero-fiber tests rather than an unspecified fifteen-variable locus.

There is also an exact triangular reduction relevant to the corrected
moment set. In the normalized coordinates,

\[
 \frac{\mu_2}{24}=
 -3s_0s_6+18s_1s_5-45s_2s_4+30s_3^2
 +14t_0t_4-56t_1t_3+42t_2^2+70.                         \tag{5.9}
\]

Consequently \(\mu_2=0\) solves one variable without saturation on every
representative chart:

\[
\begin{array}{c|ccccc}
\text{chart}&s_0=1&s_1=1&s_2=1&t_0=1&t_1=1\\ \hline
\text{variable}&s_6&s_5&s_4&t_4&t_3\\
\partial\mu_2&-72&432&-1080&336&-1344.
\end{array}                                               \tag{5.10}
\]

Thus \(\mu_2,\ldots,\mu_{12},\mu_{14}\) reduces exactly on every
chart to eleven equations in ten effective variables.

The next step is explicit on \(s_0=1\). After solving (5.9) for \(s_6\),
the third moment has the form

\[
 \mu_3=-103680\,A\,s_5-17280\,B\,t_4+C,                  \tag{5.11}
\]

where \(C\) is independent of \(s_5,t_4\), and

\[
\begin{aligned}
A={}&6s_1^2t_1-3s_1s_2t_0-3s_1t_2-3s_2t_1
      +2s_3t_0-3t_0+t_3,\\
B={}&12s_1s_3+28s_1t_0t_1-18s_1-9s_2^2-14s_2t_0^2\\
   &\quad-3s_4-2t_0t_2-12t_1^2.
\end{aligned}                                             \tag{5.12}
\]

Hence this chart splits into \(A\ne0\), where \(\mu_3\) eliminates
\(s_5\); \(A=0,B\ne0\), where it eliminates \(t_4\); and the common
boundary \(A=B=0\). This is an exact two-step triangularization.

There is a further constant-pivot layer which is useful computationally
and does not require a saturation.  Namely,

\[
 \frac{\partial A}{\partial t_3}=1,\qquad
 \frac{\partial A}{\partial s_4}=0,\qquad
 \frac{\partial B}{\partial s_4}=-3,\qquad
 \frac{\partial B}{\partial t_3}=0.                       \tag{5.12a}
\]

Thus \(A=0\) eliminates \(t_3\) globally.  On \(A=0,B\ne0\), the
\(\mu_3\) pivot then eliminates \(t_4\), leaving nine effective variables
after the inverse for \(B\) is included.  On \(A=B=0\), equations
\(A,B\) eliminate \(t_3,s_4\) globally, leaving eight effective variables.
In particular, the common-boundary system should not be submitted as ten
variables together with the two equations \(A,B\); the substitutions in
(5.12a) are exact affine coordinate eliminations.

The resulting eight-variable boundary has a useful base--fiber form.
Take the rational point

\[
 (s_1,s_2,s_3,t_0,t_1,t_2)
 =\left(-3,-3,\frac{45973}{1026},-3,3,-3\right).          \tag{5.12b}
\]

After using \(A=B=\mu_2=0\) to recover \(t_3,s_4,s_6\), the third
moment vanishes identically in the two remaining fiber variables
\((s_5,t_4)\).  The next two moments are

The exact polynomials \(\mu_4,\mu_5\in\mathbb Q[s_5,t_4]\) have total
degrees two and three, respectively, and are coprime.  Their
degree-reverse-lexicographic
initial monomials are \(t_4^4,s_5t_4^2,s_5^2\), so the quotient has basis

\[
 1,t_4,t_4^2,t_4^3,s_5,s_5t_4
\]

and length six.  Since finiteness and preservation of this standard
monomial basis hold on a Zariski-open neighborhood, this proves:

> **Proposition 5.2a (generic rank-six boundary fiber).** On a nonempty
> characteristic-zero open subset of the six-variable \(\mu_3=0\) base,
> the algebra \(\mathbb Q[s_5,t_4]/(\mu_4,\mu_5)\) is finite locally free
> of rank six, with basis
> \(1,t_4,t_4^2,t_4^3,s_5,s_5t_4\).

This does not exclude the common moment zero fiber.  It replaces the
eight-variable Gröbner problem on this open by six coefficient equations
on the base for every later moment.

The generic fiber calculation can be performed exactly over the rational
function field
\(K=\mathbb Q(s_1,s_2,s_3,t_0,t_1,t_2)\).  A sparse integral
construction of the moments followed by a characteristic-zero Gröbner
calculation gives a three-element basis with leading exponents
\((2,0),(1,2),(0,4)\) in \((s_5,t_4)\), quotient length six, and
six-term normal forms for each of \(\mu_6,\mu_7\).  Put

\[
 L=s_1t_0-t_1,\qquad
 Q=s_1^2-s_2-\frac{13}{3}t_0^2,                           \tag{5.12d}
\]

and

\[
\begin{aligned}
J={}&9801s_1^4-19602s_1^2s_2-23832s_1^2t_0^2
       -60840s_1t_0t_1\\
  &\quad+9801s_2^2+54252s_2t_0^2+75076t_0^4+30420t_1^2.
\end{aligned}                                             \tag{5.12e}
\]

The three Gröbner leading coefficients are nonzero rational multiples
of \(Q,L,J\), respectively.  The denominators in the normal form of
\(\mu_6\) have support \(LQ\), while those for \(\mu_7\) also contain
\(J\).  This corrects the two-prime reconstruction, which saw only the
\(\mu_6\) denominators \(LQ,LQ^2\): the later-moment principal open is
controlled by \(LQJ\), and \(J=0\) is an additional special divisor.

Those divisor bases can also be computed without submitting the full
eight-variable ideal.  Exact characteristic-zero calculations over the
appropriate rational-function fields of the remaining base variables
give:

\[
\begin{array}{c|c|c|c}
\text{stratum}&\operatorname{in}(\mu_4,\mu_5)&
 \text{standard basis}&\text{length}\\ \hline
L=0&(s_5^2,t_4^3)&
1,t_4,t_4^2,s_5,s_5t_4,s_5t_4^2&6\\
Q=0&(s_5t_4,t_4^3,s_5^3)&
1,t_4,t_4^2,s_5,s_5^2&5\\
L=Q=0&(t_4^2,s_5^2t_4,s_5^3)&
1,t_4,s_5,s_5t_4,s_5^2&5.
\end{array}                                               \tag{5.12f}
\]

The normal forms of both \(\mu_6\) and \(\mu_7\) occupy all the displayed
fiber coordinates.  Hence the continuation requires six base
coefficients on the \(L\)-divisor but only five on the generic
\(Q\)-divisor and on the intersection.

There is also a new exact constant pivot among those base coefficients.  The
\(t_4^3\)-coefficient of
\(\operatorname{NF}_{(\mu_4,\mu_5)}(\mu_6)\) is \(N/(39LQ)\), where
\(N\) has 42 terms and is affine-linear in both \(t_2\) and \(s_3\),
with no \(s_3t_2\) term.  Exact differentiation gives

\[
 \frac{\partial N}{\partial t_2}=-3903051350016000LQ,
 \qquad
 \frac{\partial}{\partial t_2}\left(\frac{N}{39LQ}\right)
 =-100078239744000.                                      \tag{5.12g}
\]

Thus on the generic \(LQJ\)-open, vanishing of this one \(\mu_6\)
coefficient eliminates \(t_2\) globally, with no further split.  This
reduces the subsequent base calculation from six variables to five
before the remaining normal-form coefficients are imposed.

For completeness, the alternate \(s_3\)-coefficient is, up to sign,

\[
\begin{aligned}
H={}&430353s_1^4-860706s_1^2s_2-1591461s_1^2t_0^2
      -1946880s_1t_0t_1\\
  &\quad+430353s_2^2+2564901s_2t_0^2
      +3802298t_0^4+973440t_1^2.                         \tag{5.12h}
\end{aligned}
\]

On \(H\ne0\), the same equation may instead eliminate \(s_3\), but
\(H=0\) is not a separate obstruction because the \(t_2\) pivot remains
constant.

The apparent pair of new quartic branches simplifies in the adapted
coordinates.  If

\[
 W=99Q+155t_0^2,
\]

then exact expansion gives

\[
 J=W^2+30420L^2,\qquad H=32J+1179QW.                    \tag{5.12i}
\]

Consequently \(J=H=0\) has no point on \(LQ\ne0\): after inverting \(Q\),
the second identity puts \(W\) in \((J,H)\), and the first then puts
\(L^2\) there.  Equivalently,

\[
 (J,H):(LQ)^\infty=(1).                                  \tag{5.12j}
\]

Thus the \(J\)- and \(H\)-special loci are disjoint on the principal
\(LQ\)-open.  Over \(\mathbb Q(\sqrt{-30420})\), the \(J\)-divisor
splits into the conjugate linear conditions
\(W=\pm\sqrt{-30420}\,L\), suggesting a two-variable
quadratic-extension continuation.

That continuation has an exact finite point even though a direct
transcendental coefficient-field calculation is expensive.  Let

\[
\begin{aligned}
m(\beta)={}&441554190069069\beta^4
 +15795130399581456\beta^3\\
&+193851580108553334\beta^2
 +319468919863825776\beta\\
&+1067521643767708429.
\end{aligned}
\]

This quartic is irreducible over \(\mathbb Q\).  In
\(\mathbb Q[\beta]/(m)\), put

\[
\alpha=
\frac{3(1026265600730531007\beta^3
+41799868694363859156\beta^2
+506411570533205547441\beta
+545569851002913527492)}
{18525795986003750110}.
\]

Exact reduction gives \(\alpha^2=-30420\).  The base point

\[
 s_1=t_0=L=1,\qquad s_3=\beta,\qquad t_2=0,\qquad
 Q=\frac{\alpha-155}{99}                                \tag{5.12j'}
\]

satisfies \(J=\mu_3=0\).  At this point,
\((\mu_4,\mu_5)\) has leading ideal

\[
 (s_5^2,t_4^3,s_5t_4^2)
\]

and hence length five with standard basis

\[
 1,t_4,t_4^2,s_5,s_5t_4.                               \tag{5.12j''}
\]

The special fiber can be promoted generically without adjoining \(J\)
as a third polynomial variable.  Work directly over

\[
 K_J=\mathbb Q(\alpha)(s_1,s_3,t_0,L,t_2),\qquad
 \alpha^2=-30420,
\]

and put \(99Q+155t_0^2=\alpha L\).  The specialized \(\mu_4,\mu_5\)
have only six and eight fiber monomials.  Representing every coefficient
as \(u+v\alpha\), fraction-free pseudo-reduction constructs three
relations with respective supports \(6,7,6\) and leading monomials

\[
 s_5^2,\qquad s_5t_4^2,\qquad t_4^3.                   \tag{5.12j'''}
\]

The first two construction reductions take two and four steps.  Of the
three Buchberger pairs, the pair \(s_5^2,t_4^3\) is removed by the
product criterion, while the only remaining pair reduces exactly to
zero in five steps over \(K_J\).  Thus (5.12j''') is the generic
initial ideal on either conjugate component of \(J=0\).  Consequently
the generic \(J\)-divisor quotient has length five with basis
(5.12j'').  The quartic point (5.12j') additionally proves that this
rank-five locus meets \(\mu_3=0\).

The same quadratic-pair arithmetic pseudo-reduces \(\mu_6\) and
\(\mu_7\) without dividing by a base polynomial.  Their normal forms
occupy all five standard monomials in respectively five and ten reduction
steps.  Thus every one of the five \(L,Q,J\) strata now has exact
\(\mu_6,\mu_7\) normal forms; this statement does not compute the
successive common-root quotients.

Generic rational-function-field calculations on both split components
in characteristics \(47\) and \(101\) independently reproduce the same
basis.  The earlier direct computation that kept \(J\) as an additional
polynomial variable did not finish within the recorded \(180\)-second
Singular bound; the fraction-free quadratic-pair calculation avoids that
software bottleneck.

> **Proposition 5.2b (generic boundary quotient and divisor bases).**
> Over \(K\), the quotient
> \(K[s_5,t_4]/(\mu_4,\mu_5)\) has length six with standard basis
> \(1,t_4,t_4^2,t_4^3,s_5,s_5t_4\).  At the generic points of
> \(L=0\), \(Q=0\), and \(L=Q=0\), its characteristic-zero lengths and
> bases are those in (5.12f).  The exact normal forms have denominator
> support and the constant \(t_2\)-pivot described in
> (5.12d)--(5.12g).
> Moreover the localized exceptional intersection is empty as in
> (5.12j).  At the generic point of \(J=0\), the quotient has the
> length-five basis (5.12j'') and initial ideal (5.12j'''); the exact
> point (5.12j') shows that this rank-five locus meets \(\mu_3=0\).

Independent reductions in characteristics \(47\) and \(101\) reproduce
all four quotient shapes.  Proposition 5.2b is a finite-quotient theorem,
not a zero-fiber exclusion.

There is one further exact finite algebra on the deepest divisor.  On
\(L=Q=0\), use \(t_2,s_5,t_4\) as fiber variables and
\(s_1,s_3,t_0\) as base variables.  The third moment is cubic in \(t_2\)
with constant leading coefficient \(933120\).  Exact characteristic-zero
reduction gives

\[
 \dim_{\mathbb Q(s_1,s_3,t_0)}
 \frac{\mathbb Q(s_1,s_3,t_0)[t_2,s_5,t_4]}
      {(\mu_3,\mu_4,\mu_5)}=15.                         \tag{5.12p}
\]

The reduced basis has eleven elements; its leading exponents are recorded
in the corrected-boundary artifact.  This promotes the deepest
\((\mu_4,\mu_5)\) rank-five algebra to a rank-fifteen algebra after the
\(\mu_3\) pivot.  Content-preserving pseudo-reduction in this algebra gives
primitive normal forms for \(\mu_6\) and \(\mu_7\) in respectively \(19\)
and \(35\) reduction steps.  Both normal forms occupy all fifteen standard
monomials and have maximum fiber degree three.

There is also an exact divisor exclusion across all five branches.  At
\(t_0=0\), the generic, \(L\), \(Q\), \(J\), and \(L=Q\) strata use the
respective principal opens

\[
 LQ(9801Q^2+30420L^2),\qquad Q,\qquad L,\qquad LQ,\qquad 1.
\]

These strata partition the adapted \((L,Q)\)-plane over an algebraic
closure.  Exact sparse linear algebra gives the unit ideal through
\(\mu_{10}\) on every stratum.  In particular, on the deepest stratum,

\[
 (\mu_3,\mu_4,\ldots,\mu_{10})
 =\mathbb Q[t_2,s_5,t_4,s_1,s_3]
 \qquad\text{on }L=Q=t_0=0.                            \tag{5.12q}
\]

Thus the entire branchwise common boundary at \(t_0=0\) is excluded; any
surviving point on any of the five strata lies on the \(t_0\)-open.

The remaining open has a smaller exact presentation than the direct
eight-variable system suggests.  Normalize \(t_0=1\), put \(u=s_0^{-1}\),
and use \(\mu_2,A,B\) to eliminate \(t_4,t_3,s_4\), respectively.  In
the resulting variables,

\[
 x=s_1^2u-s_2,\qquad \ell=s_1u-t_1,
\]

the transformed \(\mu_3\) is a 73-term base equation independent of
\(s_5,s_6\).  Over
\(\mathbb Q(s_1,s_2,s_3,t_1,t_2,u)\), exact reduction gives

\[
 \dim
 \frac{\mathbb Q(s_1,s_2,s_3,t_1,t_2,u)[s_6,s_5]}
      {(\mu_4,\mu_5)}=6,                              \tag{5.12r}
\]

with leading monomials \(s_6^2,s_6s_5^2,s_5^4\) and standard basis

\[
 1,\ s_6,\ s_5,\ s_6s_5,\ s_5^2,\ s_5^3.
\]

The exact \(\mu_6\) normal form occupies all six coordinates.  With the
reverse fiber order \((s_6,s_5)\), the three leading coefficients factor,
up to the nonzero rational constants \(311040,324,122472\), as

\[
 K,\qquad H,\qquad Q_*KJ_*H,
\]

where

\[
\begin{aligned}
 Q_*&=3x-13u,\\
 A_*&=99x-274u,\\
 J_*&=(99x-274u)^2+30420\ell^2,\\
 K&=351x-901u,\\
 H&=(99x-274u)K+121680\ell^2.
\end{aligned}
\]

These factors satisfy the simpler identities

\[
 K=4A_*-15Q_*,\qquad H=4J_*-15A_*Q_* .              \tag{5.12s}
\]

The apparent \(K\)-exception can itself be resolved into finite algebras.
On \(K=0\), equivalently
\(s_2=s_1^2u-(901/351)u\), the generic leading ideal changes to

\[
 (s_6s_5,s_6^3,s_5^4),
\]

the quotient still has length six, and the leading coefficients have
support \(\ell,\ell^3,\ell^2J_*\).  The exact \(\mu_6\) normal form again
occupies all six coordinates.  Since \(H=121680\ell^2\) on \(K=0\), the
reduced \(K=H=0\) intersection is the linear slice

\[
 s_2=s_1^2u-\frac{901}{351}u,\qquad t_1=s_1u .
\]

There \((\mu_4,\mu_5)\) has leading ideal \((s_5^2,s_6^3)\), length six,
and basis

\[
 1,\ s_6,\ s_6^2,\ s_5,\ s_6s_5,\ s_6^2s_5;
\]

\(\mu_6\) occupies all six coordinates yet again.  Together with the
separate exact \(J_*=0\) algebra, this shows that \(K\) introduces no new
positive-dimensional fiber degeneration.

The residual \(H=0\) divisor is a rational conic.  A dense
parametrization, with parameter \(r\) and
\(D=1287r^2+40560\), is

\[
 \ell=-\frac{775ru}{D},\qquad
 Q_*=-\frac{155}{33}u+r\ell,\qquad
 x=\frac{Q_*+13u}{3}.                                \tag{5.12t}
\]

Substitute \(s_2=s_1^2u-x\) and \(t_1=s_1u-\ell\).  Over the resulting
rational-function field, \((\mu_4,\mu_5)\) has leading ideal
\((s_6^2,s_5^3)\), length six, and basis

\[
 1,\ s_6,\ s_5,\ s_6s_5,\ s_5^2,\ s_6s_5^2.
\]

The omitted parametrization point has \(A_*=\ell=0\) and hence lies on
the separately certified \(J_*=0\) branch.  Thus \(H\), like \(K\), does
not cause a generic positive-dimensional fiber; lower-dimensional
specializations and its later-moment normal forms are still open.

Thus (5.12r) reduces the \(t_0\)-open continuation to one base equation,
a rank-six fiber, with generic rank-six presentations on every new
\(K,H\) leading-factor divisor.  This does not yet resolve the
lower-dimensional coefficient specializations, the inherited
\(Q_*,J_*\) later-moment radicals, or the common radical.

The first common-root step can be certified without expanding the
six-parameter norm.  At the exact rational base

\[
 (s_1,s_2,s_3,t_1,t_2,u)
 =
 \left(-3,-3,\frac{33467}{26028},0,1,3\right),       \tag{5.12u}
\]

the base equation \(\mu_3=0\), while
\[
 (Q_*,J_*,K,H)=(51,7077924,7827,26668476).
\]
Over \(\mathbb Q[s_6,s_5]\), the ideal
\((\mu_4,\mu_5)\) has the same leading monomials as in (5.12r) and
quotient length six, but
\[
 (\mu_4,\mu_5,\mu_6)=(1).
\]
Consequently the multiplication norm of \(\mu_6\) is not identically
zero on the local base component through (5.12u).  The common-root locus
is contained in a proper norm divisor there, and a nonempty
characteristic-zero base neighborhood is excluded.  This fixed-fiber
certificate is exact and completes in under a second; it does not
compute the norm divisor or cover its specializations.

The point (5.12u) lies on an exact rational curve in the same base:

\[
\begin{gathered}
 s_1=r,\quad s_2=-3,\quad t_1=0,\quad t_2=1,\quad u=3,\\
 s_3=
 \frac{2916r^3+810r^2-2511r-3025}
      {54(36r^3-r+5)}.                               \tag{5.12v}
\end{gathered}
\]

Substitution proves \(\mu_3=0\) identically.  In the rank-six algebra
defined by \((\mu_4,\mu_5)\), exact multiplication matrices give

\[
 \det M_{\mu_6}=\frac{N_{198}(r)}{D_{144}(r)},        \tag{5.12w}
\]

where the primitive numerator \(N_{198}\) is irreducible over
\(\mathbb Q\).  The normal form of \(\mu_7\) occupies all six standard
monomials.  More decisively, the coefficient of \(z\) in

\[
 \det(M_{\mu_6}+zM_{\mu_7})
\]

has numerator degree \(209\), denominator degree \(153\), and numerator
coprime to \(N_{198}\).  Hence away from the curve and border
denominators, either \(\mu_6\) is a unit or some linear combination of
\(\mu_6,\mu_7\) is a unit.

The denominator factors are exactly
\[
 (3r^2-10)^3(36r^3-r+5)^{42}
 (9801r^4-4230r^2+30625)^3
\]
for (5.12w), with multiplicities \(4,43,4\) for the next pencil
coefficient.  The quadratic and quartic factors are respectively the
specialized \(Q_*=0\) and \(J_*=0\) divisors.  Direct exact calculations
over their degree-two and degree-four number fields give quotient length
five for \((\mu_4,\mu_5)\) and
\[
 (\mu_4,\mu_5,\mu_6,\mu_7)=(1)
\]
on both.  The cubic factor is the denominator in (5.12v); its numerator
is coprime, so its roots are poles rather than affine points of the
parametrized curve.  Every defined point of the rational curve (5.12v)
is therefore excluded.  The full calculation takes under three seconds
and is the first successive Fitting exclusion on the \(t_0\)-open.

The global interpolation size can now be measured without expanding the
five-variable determinant.  Write
\[
 \mu_3=a_2s_3^2+a_1s_3+a_0
\]
over the remaining base variables.  At a specialization where this
quadratic has two distinct roots, two rank-six evaluations recover the
\(1,s_3\) residue coordinates of \(c_0=\det M_{\mu_6}\) and
\(c_1=[z]\det(M_{\mu_6}+zM_{\mu_7})\) modulo \(\mu_3\).
Fifteen exact directional reconstructions, using two unrelated base
points and the primes \(1019,2039\), each fit on 400 paired samples and
verify on 50 unused pairs.  Every line gives
\[
 \operatorname{den}(c_0)=a_2^{41}Q_*^3J_*^3,\qquad
 \operatorname{den}(c_1)=a_2^{42}Q_*^4J_*^4.         \tag{5.12x}
\]
The largest observed numerator/denominator profiles occur in the
\(u\)-direction: \(177/91,176/91\) for the two \(c_0\) coordinates and
\(185/96,184/96\) for \(c_1\).  All directional numerator slices are
nearly dense.  This is stable modular interpolation evidence and a
degree/denominator budget, not a reconstruction of the multivariate
numerators or a common-root theorem.  In particular, a tensor grid would
recreate the original blow-up and is not the next algorithm.

The same paired-root engine now evaluates all seven coefficients
\(c_0,\ldots,c_6\) of the determinant pencil.  Forty-four deterministic
random shards at \(p=43,47,59,71\), each protected by a \(20\)-second
internal guard, tried \(46475\) base points and accepted \(19800\) bases
where \(\mu_3(s_3)\) has two distinct roots.  Thus the bounded scout
tested
\[
 39600\ \text{\(\mu_3\)-roots},\qquad
 20\ \text{points with }c_0=\cdots=c_6=0.            \tag{5.12y}
\]
Direct Gröbner replay of
\((\mu_4,\mu_5,\mu_6,\mu_7)\) at all twenty points gives a reduced
length-one quotient.  The block matrices satisfy
\(\operatorname{rank}[M_6\ M_7]=5\) and
\(\operatorname{rank}[M_6\ M_7\ M_8]=6\) at every point.  Nineteen use
the leading five columns of \(M_6\) as a nonzero \(5\)-by-\(5\) pivot;
one replaces the fifth column by the first column of \(M_7\).  The first
column of \(M_8\) restores full rank in both charts.

Direct specialization-safe scouts at \(p=43\) also test \(900\)
\(\mu_3\)-roots on each of \(Q_*=0,J_*=0,K_*=0,H_*=0\).  The generic
\((\mu_4,\mu_5)\) lengths are five on \(Q_*,J_*\) and six on
\(K_*,H_*\), with lower length drops recorded on the first two.
\(Q_*\) and \(J_*\) each contain one sampled reduced common point through
\(\mu_7\); both are excluded by \(\mu_8\), while no such point occurs
in the \(K_*,H_*\) samples.  Reconstructing every full normalized
parameter point and evaluating the original corrected moment formula
gives \(\mu_8\ne0\).  Fixed-point evaluation stops at the first nonzero
later moment and never expands a global \(\mu_{14}\).
This is exact bounded finite-field evidence, not an exhaustion: the
common-\(\mu_6,\mu_7\) incidence is expected to have positive dimension,
so the absence of a \(\mu_8\)-survivor in these samples is not a global
or divisor exclusion.

The lower-dimensional continuation can nevertheless be made explicit
without a large elimination.  Additional direct \(p=43\) scouts cover
\(K=H=0\), the three relevant \(Q,J,K,H\) intersections, \(a_2=0\), and
the repeated-root divisor of \(\mu_3(s_3)\).  Together with enlarged
\(Q=0\) and \(J=0\) shards they retain every specialized fiber at
\(31\) deterministic artifacts.  Seventeen reduced common points through
\(\mu_7\) occur across these strata; direct evaluation of the original
corrected moment formula gives \(\mu_8\ne0\) at all seventeen.

A separate rank-complement scout evaluates \(6300\) further
\(\mu_3\)-roots.  Four points miss both selected \(5\)-by-\(5\) pivots,
including one with \(\operatorname{rank}M_6=3\), but all four satisfy
\[
 \operatorname{rank}[M_6\ M_7]=6.
\]
No sampled point has joint rank at most four.  Four additional full-pencil
zeros through \(\mu_7\) are again excluded by \(\mu_8\).  Thus the two
named pivots are not a global chart cover, while the full joint-rank
condition remains the correct invariant.

More importantly, `liftstd` exposes the specialization border hidden by
the reduced rational-function-field bases.  The least common multiple of
the lifted leading coefficients is a single irreducible polynomial modulo
\(43\) on every displayed stratum:

\[
\begin{array}{c|c|c|c}
\text{stratum}&\text{generic length}&(\deg,\text{terms})&
 \text{sampled drops / border zeros}\\ \hline
Q&5&(36,578)&132/132\\
J&5&(36,1245)&136/136\\
K&6&(22,16)&0/0\\
H&6&(36,15)&0/0\\
K\cap H&6&(4,1)&0/0\\
Q\cap J\cap H&4&(59,2997)&171/263\\
J\cap H&4&(57,1399)&151/226\\
J\cap K&5&(29,166)&77/77\\
a_2=0&6&(33,187)&0/0 .
\end{array}                                             \tag{5.12z}
\]

Thus every sampled length drop lies on the exact modular border.
Conversely the border is exact on the \(Q,J,J\cap K\) samples; on the two
deeper \(H\)-intersections some border points use another basis of the
same length.  Enlarging the \(Q,J\) point clouds gives \(132,136\) drop
points.  In both cases the degree-at-most-four evaluation matrix has
\(126\) columns and rank \(120\); its six-dimensional kernel is exactly
the ambient cubic component equation and its five linear multiples.
There is no extra drop equation through degree four.

Since both \(\mu_3\) and the borders are quadratic in \(s_3\), their
projection is a small resultant rather than a multivariate Gröbner
problem.  Factoring these resultants modulo \(43\), after removing
\(u=0\) and inherited linear intersections, gives the residual components

\[
\begin{array}{c|c}
\text{stratum}&(\deg,\text{terms},\text{multiplicity})\\ \hline
Q&(20,195,2)\\
J&(24,612,2)\\
J\cap K&(24,161,2)\\
Q\cap J\cap H&(20,101,1),(24,161,1),(24,163,3)\\
J\cap H&(20,54,3),(24,143,1),(24,143,1).
\end{array}                                             \tag{5.12aa}
\]

Writing the two quadratics as
\(as_3^2+bs_3+c\) and \(ds_3^2+es_3+f\), their first
pseudo-remainder is
\[
 A s_3+B,\qquad A=db-ae,\quad B=dc-af.                \tag{5.12ab}
\]
On every residual factor in (5.12aa), direct modular gcd computation gives
\(\gcd(A,R)=\gcd(B,R)=1\).  Hence every such component has the dense
linear pivot \(s_3=-B/A\); the failures of this pivot are precisely the
recorded \(u=0\) and inherited linear factors.

The \(Q=0\) row admits a direct characteristic-zero promotion.  Over
\(\mathbb Q(s_1,s_3,t_1,t_2,u)\), the lifted leading border is irreducible
of degree \(36\) with \(588\) terms.  Its exact resultant with \(\mu_3\)
has \(5563\) terms, degree \(76\), and factors as
\[
 c\,u^{20}\,J_Q^4\,R_{20}^2,                          \tag{5.12ac}
\]
where \(J_Q\) is the irreducible four-term quartic cutting out the
inherited \(Q\cap J\) locus and \(R_{20}\) is irreducible over
\(\mathbb Q\), with degree \(20\) and \(200\) terms.  The exact
pseudo-remainder coefficients \(A,B\) have term/degree profiles
\((262,33)\), \((535,38)\), and both are coprime to \(R_{20}\).
Reduction modulo \(43\) gives the corresponding \(Q\)-row data above
up to nonzero scalars.

There is now one exact characteristic-zero component exclusion after the
dense pivot.  On the slice
\[
 s_1=0,\qquad \ell=s_1u-t_1=0,\qquad T=t_2/u^2,
 \qquad u=s_0^{-1},                                    \tag{5.12ad}
\]
the residual degree-five equation gives a quartic Kummer extension and
\((\mu _4,\mu _5)\) gives the length-four fiber basis
\[
 1,\ s_5,\ s_6,\ s_5^2.
\]
The multiplication norm of \(\mu _6\) has two irreducible factors in
\(\mathbb Q[T]\): one of degree \(100\), and
\[
 N_3(T)=909T^3+25521T^2+6189560T-22223500.             \tag{5.12ae}
\]
The degree-\(100\) factor is coprime to a coefficient of
\(\det(M_{\mu _6}+zM_{\mu _7})\), so the \((\mu _6,\mu _7)\)-ideal is
the unit ideal there.  Over \(\mathbb Q[T]/(N_3)\), the exact Kummer
modulus becomes
\[
 u^4.                                                   \tag{5.12af}
\]
Thus the complete scheme-theoretic cubic fiber is supported at \(u=0\)
and disappears after the localization \(u=s_0^{-1}\ne0\).  Consequently:

> **Proposition 5.7.**  The defined slice (5.12ad) contains no common
> zero of \(\mu _3,\ldots,\mu _7\).  In particular, the apparent cubic
> norm survivor is a localization-boundary component, not a
> characteristic-zero moment-zero point.

This closes the whole one-parameter slice, not the unspecialized
\(Q\)-residual component.

The residual continuation must use
\[
 \ell=s_1u-t_1,\qquad t_1=s_1u-\ell.                  \tag{5.12ag}
\]
An earlier exploratory driver instead substituted \(t_1=u(s_1-\ell)\);
that formula is valid only when its displayed \(\ell\) is reinterpreted
as the ratio \(\lambda=\ell/u\).  The corrected exact closed-fibre
calculations at
\[
 (s_1,\ell,u)=(5,7,2),\ (7,4,3),\ (11,9,5)           \tag{5.12ah}
\]
all give an irreducible degree-five specialization of \(R_{20}\), verify
that the dense pivot coefficient is a unit, annihilate both \(\mu _3\)
and the leading border, and produce a length-four
\((\mu _4,\mu _5)\)-algebra.  In all three cases adjoining
\(\mu _6,\mu _7\) gives the unit ideal.  These exact computations show
that the common-zero incidence does not dominate the residual base, so
its possible image is contained in a proper Fitting-closed exceptional
locus.  They do not show that this exceptional locus is empty.

There is now a bounded black-box realization of the remaining Fitting
problem.  At a closed base point it forms the full rank-twenty algebra
over the ground field, constructs the multiplication matrices
\(M_{\mu _6},M_{\mu _7}\), and recovers all twenty-one coefficients of
\[
 \det(M_{\mu _6}+zM_{\mu _7}).                        \tag{5.12ai}
\]
At \((5,7,2)\) modulo \(43\), the pivot, border, and \(\mu _3\) checks
pass, the joint block matrix has rank \(20\), and every coefficient in
(5.12ai) is nonzero; the standard-basis and matrix part takes less than
one tenth of a second.  This is a modular oracle for interpolation, not
a characteristic-zero component exclusion.

There is also an exact projective explanation for the length drop from
six to four.  Put \(h\) for the fibre homogenizing variable.  The
highest fibre-degree parts factor as
\[
\begin{aligned}
 \mu _4^{\rm top}
   &=c_4u(6s_1us_5-s_6)
     (1092\ell s_5+930s_1us_5-155s_6),\\
 \mu _5^{\rm top}
   &=c_5(6s_1us_5-s_6)
     (180s_1^2u^2s_5^2-60s_1us_5s_6
       +196u^2s_5^2+5s_6^2),
\end{aligned}                                             \tag{5.12aj}
\]
with \(c_4,c_5\ne0\).  The resultant of the two remaining factors on
the infinity line is
\[
 980(6084\ell^2+4805u^2)=980J_Q.             \tag{5.12ak}
\]
Moreover, if \(D\) denotes the residual border factor, exact symbolic
calculation gives
\[
 \operatorname{border}|_Q=c\,u^{10}J_Q^2D^2,\qquad
 \det(\text{two tangent rows at }[6s_1u:1:0])
   =c'u^5D.                                    \tag{5.12al}
\]
Consequently, on the actual residual open \(uJ_Q\ne0,D=0\), the
projective \((\mu _4,\mu _5)\)-intersection has one forced double point
at infinity,
\[
 [s_6:s_5:h]=[6s_1u:1:0],                     \tag{5.12am}
\]
and a residual affine scheme of length four.  The exact closed probe at
\((5,7,2)\) verifies coprimality of the full homogenizations, total
intersection number \(2\cdot3=6\), and infinity length two.  This is the
complete-intersection geometry behind (5.12ai), not merely a pattern in
the three closed fibres.

This structure supplies two smaller elimination formulations.  First,
the ordinary ternary resultant
\[
 \operatorname{Res}_{\mathbb P^2}
 \bigl(\mu _4^h,\mu _5^h,h\mu _6^h+z\mu _7^h\bigr)          \tag{5.12an}
\]
has a forced \(z^2\) infinity factor on \(R_{20}\); its quotient is the
degree-four pencil over the residual quintic whose norm is (5.12ai).
This is the residual-resultant situation studied by
[Busé--Elkadi--Mourrain](https://doi.org/10.1016/S0022-4049(00)00144-4),
while the explicit-factor analysis of iterated eliminations in
[Busé--Mourrain](https://arxiv.org/abs/cs/0612050) explains how to
separate the infinity/projection factors.  Second, division by the
quadratic \(\mu _4\) makes each of \(\mu _5,\mu _6,\mu _7\) linear in
\(s_6\).  At the exact closed probe their bidegrees in
\((s_6,s_5)\) are respectively
\[
 (1,3),\qquad(1,3),\qquad(1,4),               \tag{5.12ao}
\]
and these three remainders together with \(\mu _4\) generate the unit
ideal.  Thus one may eliminate \(s_6\) explicitly and replace the
rank-twenty pencil by univariate equations in \(s_5\) of degrees at most
\(6,5,6\), with the separate exceptional case \(A_5=B_5=0\).  This is
the next exact route; it has substantially smaller subresultant matrices
than either generic Gröbner basis attempted below.

The order of these two eliminations matters computationally.  If the dense
pivot \(s_3=-B/(6A)\) is substituted first, even the modular
\(\mu _5\)-remainder takes about \(212\) seconds and has coefficients with
up to \(624\,000\) terms.  Pseudo-dividing by \(\mu _4\) first leaves,
modulo \(43\), only \(25,29,41\) monomials in the raw linear remainders of
\(\mu _5,\mu _6,\mu _7\), with largest coefficient supports
\(433,714,1130\).  Forming
\[
\begin{aligned}
 E_4&=aQ_5^2-f_1P_5Q_5+f_0P_5^2,\\
 E_6&=P_5Q_6-P_6Q_5,\qquad
 E_7=P_5Q_7-P_7Q_5
\end{aligned}                                             \tag{5.12ap}
\]
before the pivot takes less than two seconds.  Their
\((\deg_{s_5},\deg_{s_3},\text{terms})\) profiles are
\[
 E_4:(5,12,46),\qquad E_6:(5,10,40),\qquad
 E_7:(6,12,48).                                           \tag{5.12aq}
\]
The coefficient content of each is a unit modulo \(43\).  These are exact
finite-field identities over
\(\mathbb F_{43}(s_1,\lambda,v)[T]/(R_{20})\), not a
characteristic-zero exclusion.

By contrast, explicitly substituting the pivot into \(E_6\) does not
finish within the recorded \(180\)-second bound, although it stays below
\(300\) MB.  The bounded formulation is therefore the four-equation
system
\[
 6A s_3+B=E_4=E_6=E_7=0                                  \tag{5.12ar}
\]
in \((s_5,s_3)\), together with the exceptional ideal
\((P_5,Q_5)\), rather than the expanded univariate coefficients.  A
residual resultant or Bezoutian/Macaulay matrix can consume (5.12ar)
without expanding the linear pivot.  This is also compatible with the
trace-matrix construction of
[Janovitz-Freireich--Mourrain--Rónyai--Szántó](https://arxiv.org/abs/0901.2778),
which gives Macaulay/Bezoutian formulations for finite projective roots
and discusses components at infinity.

The batched pencil oracle also quantifies why dense three-variable
interpolation is not attractive.  Modulo \(1009\), 499 of the first 500
points on the \(s_1\)-line \((\ell,u)=(7,2)\) are good, and rational
reconstruction from 400 points validates on all 99 held-out points.  A
representative coefficient has numerator/denominator degrees \(209/91\)
in \(s_1\); transverse full-field scans give \(420/212\) in \(\ell\),
and \(270/100\) in \(v=u^2\) on the natural \(\lambda=0\) chart.  These
are modular degree scouts, not characteristic-zero bounds, but a dense
box at those degrees would contain tens of millions of monomials.  The
univariate-remainder construction (5.12ao), rather than dense
interpolation of (5.12ai), is therefore the bounded continuation.

The two remaining direct symbolic implementations have also been
localized as resource dead ends.  Function-field standard-basis
computations fail in coefficient normalization under the recorded
memory caps.  In the corrected ratio chart
\(\lambda=(s_1u-t_1)/u,\ v=u^2\), sparse FLINT arithmetic computes
\(R_{20}\), the dense pivot, and all four input moments quickly, and
substitutes the pivot into \(\mu _4,\mu _5\) modulo \(43\) in about
thirty seconds.  The subsequent unspecialized triangular-basis step
still exceeds four minutes under a \(3\)-GB cap and was stopped.
Equations (5.12ap)--(5.12ar) remove that memory bottleneck: all three
pre-pivot remainders and all three raw elimination equations are formed
in about five seconds and below \(300\) MB.  The remaining viable routes
are the projective resultant (5.12an) and a resultant/subresultant
calculation on the compact system (5.12ar), followed by recursion only
on their explicit exceptional factors.

The remaining rows of (5.12z)--(5.12ab) are exact finite-field
calculations and sampled incidence checks, not characteristic-zero
factorizations or component exclusions.  Even on the exact
\(Q\)-component, general-purpose Gröbner recomputation exceeds the
\(20\)-second guard after the linear pivot.  The next bounded
implementation should therefore perform unspecialized custom quotient
arithmetic in the degree-five and degree-six extensions defined by the
residual factors, reduce \(\mu_4,\ldots,\mu_8\) there, and recurse only
on the exceptional subresultant loci.  Proposition 5.7 supplies one exact
specialized calibration for that arithmetic.

An attempted branchwise continuation through corrected \(\mu_8\) exposed
an important content issue.  Singular's `cleardenom` removes polynomial
content, so applying it to a fiber-constant base equation can replace that
equation by \(1\).  Such an export gives a spurious unit ideal and is not a
Nullstellensatz certificate.  The corrected exporter preserves base
content explicitly.  With that correction, the reordered \(600\)-second
mod-\(47\) deepest solve still does not finish.  Exact and mod-\(47\)
fraction-free reductions of \(\mu_8\) also exceed the recorded
\(600\)-second bound.  Thus \(\mu_6,\mu_7\) are reduced exactly on the
deepest branch and (5.12q), together with the other four saturated unit
computations, closes the whole \(t_0=0\) divisor.  On the \(t_0\)-open,
(5.12r) reduces \(\mu_6\) in a generic rank-six algebra and (5.12u)
proves that its first norm is nonzero on one exact local base component,
while (5.12v)--(5.12w) close one dense rational curve through its
irreducible norm divisor.  The bounded scout (5.12y) confirms that the
full determinant-pencil test and direct higher-moment replay can be run
in short shards, but it does not replace elimination.  The next exact
step is to construct the common-pencil incidence itself and intersect
its components with \(\mu_8\), before continuing to later moments.  The
global norm divisor, its successive common-root/Fitting equations, the
lower-dimensional
coefficient strata and inherited \(Q_*,J_*\) branch radicals, and
corrected orders \(7,8,9,10,11,12,14\) remain unresolved.  Neither a
semistable moment-zero point on this normalized non-null branch nor a
radical equality follows from the present certificates.  The Rodrigues
point found below lies instead on the null-quadratic branch.

### 5.15 The anti-Weyl locus

There is one exact involutive subbranch of the normalized non-null
quadratic chart.  Put
\[
\omega_a(X,T)=(aT,-a^{-1}X),\qquad q=a^2.
\]
The condition \(\omega_a(F)=-F\) is
\[
\begin{gathered}
s_4=-qs_2,\qquad s_5=q^2s_1,\qquad s_6=-q^3s_0,\\
t_2=0,\qquad t_3=qt_1,\qquad t_4=-q^2t_0.              \tag{5.15a}
\end{gathered}
\]
The residual torus conjugates \(a\) to one.  The retained coordinates are
\[
(s_0,s_1,s_2,s_3,t_0,t_1),
\]
and (5.15a) becomes
\[
s_4=-s_2,\quad s_5=s_1,\quad s_6=-s_0,\quad
t_2=0,\quad t_3=t_1,\quad t_4=-t_0.                   \tag{5.15b}
\]
Haar invariance gives
\[
\mu_m(F)=\mu_m(\omega_1F)=(-1)^m\mu_m(F),
\]
so every odd moment vanishes identically.  Up to the common nonzero
factor in the raw contraction, the quadratic equation is
\[
3s_0^2+18s_1^2+45s_2^2+30s_3^2
-14t_0^2-56t_1^2+70.                                  \tag{5.15c}
\]

Exact characteristic-zero `msolve` computation gives
\[
(\mu_2,\mu_4,\mu_6,\mu_8,\mu_{10},\mu_{12},\mu_{14})
=\mathbb Q[s_0,s_1,s_2,s_3,t_0,t_1].                  \tag{5.15d}
\]

> **Proposition 5.6 (anti-Weyl exclusion).** The anti-Weyl locus
> (5.15a) in the normalized non-null quadratic branch contains no
> all-order pure-moment point and hence no SIC counterexample.

The exact unit output is stored in
[`two_pair_sic_bidegree33_anti_weyl_normalized_msolve14_char0.json`](../artifacts/generated-results/two_pair_sic_bidegree33_anti_weyl_normalized_msolve14_char0.json).

The full even-moment unit ideal does not say where a finite-prefix component
first fails.  The
[isolated rank-two finite-prefix theorem](TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_FINITE_PREFIX.md)
supplies that missing local calculation.  It gives a rational isolating box
for a unique exact-rank-two semistable zero of
\(\mu_1,\ldots,\mu_9\), proves that the projective quotient point is reduced
and isolated, and certifies \(\mu_{10}>0,\mu_{12}<0,\mu_{14}>0\).
Thus this explicit realization of the Hilbert-series component does not
survive either corrected system.  The same checker derives the five-variable
square-invariant quotient and proves over \(\mathbb Q\) that its two rank
quartics together with \(\mu_4,\mu_6,\mu_8,\mu_{12}\) generate the unit
ideal.  Hence the corrected rank-two system has no anti-Weyl point at all.

### Trace/norm reconnaissance on the \(L\)-open

There is a second useful finite algebra after weighted normalization
\(L=1\).  Regard \(s_3,s_5,t_4\) as fiber variables and
\(s_1,t_0,Q,t_2\) as base variables.  Modulo \(47\), the quotient

\[
 \mathbb F_{47}(s_1,t_0,Q,t_2)[s_3,s_5,t_4]/
 (\mu_3,\mu_4,\mu_5)
\]

has length twelve and leading ideal

\[
 (s_3^2,s_5^2,s_5t_4^2,t_4^4).                       \tag{5.12b'}
\]

Its standard basis is

\[
\begin{split}
 1,t_4,t_4^2,t_4^3,s_5,s_5t_4,\;&s_3,s_3t_4,
 s_3t_4^2,s_3t_4^3,\\
 &s_3s_5,s_3s_5t_4.
\end{split}
\]

Exact generator multiplication matrices and the normal form of \(\mu_6\)
have been exported over the rational-function field.  At specialized
base points, the checker constructs the multiplication matrices for
\(\mu_6,\mu_7\), their traces, norms and characteristic polynomials, and
the rank of the joint map

\[
 [M_{\mu_6}\ M_{\mu_7}]\colon A^{\oplus2}\longrightarrow A.
\]

Overflow-safe specialization is essential: base values are substituted
inside the finite-field ring before passage to the three-variable fiber
ring.

In the recorded deterministic reconnaissance, all \(1200\) accepted
points modulo \(47\) and all \(250\) accepted points modulo \(101\) have
the same leading ideal (5.12b').  Modulo \(47\), three points have joint
rank eleven and a reduced common quotient of length one.  Their base
coordinates \((s_1,t_0,Q,t_2)\) are

\[
 (37,42,6,10),\qquad(36,21,17,21),\qquad(40,19,7,17).
\]

Direct replay with the corrected moment set
\(\mu_1,\ldots,\mu_{12},\mu_{14}\) gives respective \(\mu_8\)-values
\(-9,9,7\), so all three sampled trace candidates are excluded already
by \(\mu_8\).  These calculations are finite-field reconnaissance, not
an exhaustion of the rank-drop locus and not a characteristic-zero
nullcone certificate.

There is also one complete characteristic-zero slice certificate.  Fix

\[
 L=s_1=t_0=1.
\]

The remaining variables are \(Q,t_2,s_3,s_5,t_4\).  Modulo both \(47\)
and \(101\), the ideal
\((\mu_3,\ldots,\mu_7)\) is zero-dimensional of length \(1128\).
Adjoining corrected \(\mu_8\) gives the unit ideal in both
characteristics.  More importantly, Singular's verified modular
reconstruction (`modStd` with exactness one) gives

\[
 (\mu_3,\mu_4,\mu_5,\mu_6,\mu_7,\mu_8)
 =\mathbb Q[Q,t_2,s_3,s_5,t_4]                         \tag{5.12b''}
\]

on this exact slice.  Thus the whole five-variable slice contains no
corrected zero-fiber point.  Equation (5.12b'') is a
characteristic-zero slice exclusion, not a statement about the full
four-dimensional base.

Two six-variable hyperslices admit the same exact treatment.  Keeping
\(s_1\) free and fixing \(t_0=1\), or keeping \(t_0\) free and fixing
\(s_1=1\), `msolve` computes

\[
 (\mu_3,\ldots,\mu_8)=(1)
\]

modulo both \(47\) and \(101\), and also directly over \(\mathbb Q\).
Thus

\[
 (\mu_3,\ldots,\mu_8)
 =\mathbb Q[s_1,Q,t_2,s_3,s_5,t_4]
 \quad\text{when }t_0=1,                              \tag{5.12b'''}
\]

and the analogous equality holds in
\(\mathbb Q[t_0,Q,t_2,s_3,s_5,t_4]\) when \(s_1=1\).
This excludes both complete hyperslices in characteristic zero.
Singular's verified modular reconstruction of the first hyperslice did
not finish within the recorded \(600\)-second bound, but the independent
exact `msolve` computations returned the one-element Gröbner basis
\([1]\) using deterministic exact sparse linear algebra.  The unfixed
seven-variable \(L=1\) system did not finish
within its recorded \(600\)-second exact `msolve` run, so these two
hyperslice exclusions are not promoted to a full \(L\)-open theorem.

The constant \(t_2\)-pivot permits one more exact elimination without a
five-variable Gröbner basis.  Use the adapted base coordinates
\((s_1,s_3,t_0,L,Q)\), solve (5.12g) for \(t_2\), and let
\(P(s_3)\) be the numerator of the resulting third moment.  Then \(P\)
has 642 terms, total degree 21, and degree three in \(s_3\).  Since the
rank-six \((s_5,t_4)\)-basis is uniform after localizing the three
leading coefficients \(LQJ\), adjoining \(P\) gives a generic
rank-\(18\) iterated quotient over
\(\mathbb Q(s_1,t_0,L,Q)\).

After the \(t_2\)-pivot, the leading remaining coefficient of
\(\operatorname{NF}(\mu_6)\) is the coefficient \(C(s_3)\) of
\(s_5t_4\).  It is also cubic in \(s_3\), and exact elimination gives

\[
 \operatorname{Res}_{s_3}(P,C)=L^6Q^6\mathcal R_{63},    \tag{5.12k}
\]

where \(\mathcal R_{63}\in
\mathbb Q[s_1,t_0,L,Q]\) has total degree 63 and 6702 terms.  Its
degree-preserving reduction modulo \(47\) is irreducible (with 6565
surviving terms), so \(\mathcal R_{63}\) is irreducible over
\(\mathbb Q\).  Reduction modulo \(101\) independently gives an
irreducible degree-63 polynomial with 6633 surviving terms.

The next \(t_4^2\)-coefficient \(C_2(s_3)\) is again cubic.  Its exact
resultant has the form

\[
 \operatorname{Res}_{s_3}(P,C_2)=L^9Q^6\mathcal T_{66},  \tag{5.12l}
\]

where \(\mathcal T_{66}\) has degree 66.  Its reductions modulo \(47\)
and \(101\) are irreducible of degree 66, with respectively 6951 and
7038 surviving terms, and in both characteristics

\[
 \gcd(\mathcal R_{63},\mathcal T_{66})=1.
\]

One good reduction already proves coprimality over \(\mathbb Q\).
Therefore simultaneous vanishing of the first two remaining
\(\mu_6\)-coefficients has codimension at least two in the
four-parameter principal base; the degree-63 divisor is not itself a
component of the full \(\mu_6\)-zero locus.

> **Proposition 5.2c (first residual base divisor).** On the
> \(LQJ\)-open, the constant \(t_2\)-pivot and
> \(\mu_3,\mu_4,\mu_5\) give a generic rank-18 finite quotient.
> Vanishing of the next \(s_5t_4\)-coefficient of \(\mu_6\) forces the
> base onto the irreducible divisor \(\mathcal R_{63}=0\).  Imposing
> also the \(t_4^2\)-coefficient confines it to the codimension-at-least
> two intersection
> \(\mathcal R_{63}=\mathcal T_{66}=0\).

There is a further exact simplification on the degree-63 divisor.  Write

\[
\begin{aligned}
 P&=as_3^3+bs_3^2+cs_3+d,\\
 C&=es_3^3+fs_3^2+gs_3+h,
\end{aligned}
\]

and set

\[
 A=eb-af,\qquad B=ec-ag,\qquad C_0=ed-ah.
\]

Twice taking the pseudo-remainder, without division, gives the linear
polynomial

\[
 V_1s_3+V_0,\qquad
 \begin{cases}
 V_1=A(Ac-aC_0)-(Ab-aB)B,\\
 V_0=A^2d-(Ab-aB)C_0.
 \end{cases}                                             \tag{5.12m}
\]

It vanishes at every common zero of \(P\) and \(C\).  Exact expansion
gives

\[
\deg V_1=60,\quad |V_1|=2105,\qquad
\deg V_0=63,\quad |V_0|=5170,                            \tag{5.12n}
\]

where \(|\cdot|\) denotes the number of rational terms.  In
characteristics \(47\) and \(101\), direct gcd calculations give

\[
 \gcd(\mathcal R_{63},V_1)
 =\gcd(\mathcal R_{63},V_0)=1.
\]

Either good reduction proves the corresponding characteristic-zero
coprimality.  Hence \(V_1\) is nonzero on a dense open of the irreducible
divisor \(\mathcal R_{63}=0\), and there

\[
 s_3=-V_0/V_1.                                          \tag{5.12o}
\]

> **Proposition 5.2d (dense linear pivot on the residual divisor).**
> On a dense open of the irreducible degree-63 residual divisor, the
> common cubic equations \(P=C=0\) have the rational linear pivot
> (5.12o).  Thus the divisorial incidence cover is birational to its
> base rather than a generically cubic extension.

Thus the generic continuation is no longer an unrestricted
five-variable coefficient ideal or even a divisorial branch: it is
confined to a codimension-at-least-two locus in four adapted parameters.
The pivot (5.12o) simplifies the generic degree-63 branch, but it does
not prove that \(V_1\) stays nonzero on every component of
\(\mathcal R_{63}=\mathcal T_{66}=0\).  That exceptional subresultant
intersection, the remaining four \(\mu_6\)-coefficients, and the separate
lower-dimensional coefficient ideal on \(J=0\) remain open.  The latter
now sits in a generic finite rank-five quotient; the remaining
obstruction is its moment-coordinate ideal, not generic finiteness.

One natural plane in the common boundary can already be closed. Put

\[
\begin{gathered}
s_0=1,\qquad
s_1=s_2=s_3=s_4=s_5=t_1=t_2=0,\\
t_0=a,\qquad t_3=3a,\qquad t_4=b,\qquad
s_6=\frac{14ab+70}{3}.                                   \tag{5.13}
\end{gathered}
\]

Then \(A=B=\mu_2=0\) identically, while exact contraction gives

\[
\mu_3=1866240\,a^3,\qquad
\mu_4=138240(11249-8776ab-901a^2b^2).                    \tag{5.14}
\]

Writing \(u=ab\) and \(q=11249-8776u-901u^2\), the polynomial

\[
r=\frac1{11249}
  +\frac{8776}{11249^2}u
  +\frac{8776^2+901\cdot11249}{11249^3}u^2
\]

satisfies \(rq\equiv1\pmod{a^3}\). Therefore (5.14) gives an explicit
certificate

\[
 1\in(\mu_3,\mu_4)\subset\mathbb Q[a,b].
\]

> **Proposition 5.3.** The two-parameter common-boundary slice (5.13)
> contains no moment-zero point.

The plane sits in a larger exact family.  Introduce \(h,q\) and put

\[
\begin{gathered}
s_0=1,\qquad s_1=s_2=s_3=t_2=0,\\
s_4=-4q^2,\qquad s_5=h,\qquad
t_0=a,\qquad t_1=q,\qquad t_3=3a,\qquad t_4=b,\\
s_6=\frac{14ab-168aq+70}{3}.
\end{gathered}                                            \tag{5.14a}
\]

Then \(A=B=\mu_2=0\) identically.  After removing nonzero rational
contents, the third moment is

\[
 \mu_3^{\mathrm{prim}}=3a^3+4aq^4+8q^3.
\tag{5.14b}
\]

Let \(I_{3\ldots r}\subset\mathbb Q[a,b,h,q]\) be generated by the
restricted primitive moments of orders \(3,\ldots,r\).  Exact
characteristic-zero Gröbner computation gives

\[
\begin{array}{c|cccc}
 &I_{3,4}&I_{3,4,5}&I_{3,4,5,6}&I_{3,4,5,6,7}\\ \hline
\dim&2&1&0&-1\\
\text{quotient length}&-&-&372&0.
\end{array}                                               \tag{5.14c}
\]

In particular,

\[
 I_{3,4,5,6,7}=(1).
\tag{5.14d}
\]

> **Proposition 5.4.** The four-parameter common-boundary family (5.14a)
> contains no moment-zero point.

The earlier plane (5.13) is the specialization \(q=h=0\).  Proposition
5.4 still closes only a structured family, not the full \(A=B=0\)
boundary.

The three pivot strata are not intrinsically rank-degenerate. In the
coordinate order

\[
(s_0,s_1,s_2,s_3,s_4,s_5,s_6,t_0,t_1,t_2,t_3,t_4),
\]

take

\[
\begin{aligned}
P_A&=(1,0,2,0,0,2,3,-2,-1,3,-2,-2),\\
P_B&=(1,0,0,-3,3,2,3,-2,1,2,-18,2),\\
P_0&=(1,0,2,0,-148/3,2,3,-2,-1,3,-12,-2).
\end{aligned}                                             \tag{5.15}
\]

They satisfy

\[
(A,B)(P_A)=(10,-148),\qquad
(A,B)(P_B)=(0,-13),\qquad
(A,B)(P_0)=(0,0).
\]

At every point in (5.15), the full \(11\times11\) chart Jacobian of
\(\mu_2,\ldots,\mu_{12}\) is invertible over \(\mathbb Q\). Moreover,
\(\partial A/\partial t_3=1\) and
\(\partial B/\partial s_4=-3\), so the imposed boundary equations are
independent.

> **Proposition 5.5.** The restrictions of
> \((\mu_2,\ldots,\mu_{12})\) to \(A\ne0\), to
> \(A=0,B\ne0\), and to \(A=B=0\) have generic differential ranks
> \(11,10,9\), respectively.

Thus the branch difficulty is not a forced differential-rank defect. The
appropriate target is a finite quotient or resultant on each stratum,
with \(\mu_{14}\) tested afterward.

Direct modular reconnaissance has not yet closed the first chart. The
corrected system timed out within the recorded \(180\)-second runs over
\(\mathbb F_{101}\) and \(\mathbb F_{43}\), with both Singular and
`msolve`, after direct submission and the first pivot.  After the
additional constant substitutions (5.12a), the eight-variable
\(A=B=0\) boundary still timed out in recorded \(180\)-second Singular
and `msolve` runs over \(\mathbb F_{47}\) and \(\mathbb F_{101}\).
Export of the fully substituted \(A\ne0\) branch produces
eleven equations in ten variables, but the recorded `msolve` run
terminated inside the solver. Sparse encodings retain \(s_5\) or \(t_4\)
and add the inverse-pivot relation, avoiding the expansion; the full
corrected \(A\)- and \(B\)-open systems nevertheless exceed the recorded
bounds. These are not mathematical results. The next computation should
construct the finite special quotient of \(\mu_2,\ldots,\mu_{12}\) on
each stratum and test \(\mu_{14}\) there, while decomposing the common
boundary beyond (5.13).

## 6. Degree-three continuation after the Rodrigues survivor

The former global quadratic-anchor target is false.  The
[Rodrigues survivor](TWO_PAIR_SIC_BIDEGREE33_RODRIGUES_SURVIVOR.md)
has null but nonzero quadratic component, invertible coefficient matrix,
and every pure moment zero.  Its all-order beta--Rodrigues identity
nevertheless proves eventual mixed vanishing for every fixed multiplier.
Thus the degree-three task is no longer to prove moment--nullcone
equality.  It is to classify the remaining semistable pure-moment orbits
and prove each one SIC-safe, or to find one with a persistent regular
multiplier.

The
[moment--nullcone program](TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md)
still places (5.4)--(5.5) inside the decomposition

\[
 \operatorname{End}(\operatorname{Sym}^d)
 \cong\bigoplus_{r=0}^d\operatorname{Sym}^{2r}
\]

for arbitrary \(d\), but the relevant degree-three split now has three
parts:

1. the normalized non-null quadratic branch of Section 5, including its
   explicit residual \(Q\)- and \(J\)-divisors;
2. the nonzero null-quadratic charts, where the Rodrigues orbit is locally
   isolated to length five through \(\mu_{11}\); and
3. the boundary \(F_2=0\).

The complete sparse census proves that a new SIC counterexample, if one
exists, has at least nine nonzero entries in the displayed coefficient
basis.  For size six, full complex degree-one RURs show that the two
nonunit systems contain only the Rodrigues copies.  At size eight the
sole hard odd-parity support is an exact characteristic-zero unit through
\(\mu_{14}\).  The remaining fourteen systems contain twenty normalized
complex points in total.  Twelve are the rectangular families

\[
\begin{aligned}
 &(W^{3-a}V^a-3^{b-a}W^{3-b}V^b)(Z+Y/3)^3,\\
 &(W-rV)^3(Z^{3-c}Y^c+Z^{3-d}Y^d),
 \qquad r^{d-c}=(-1)^{d-c+1},
\end{aligned}
\]

and are visibly one-sided.  One exceptional system consists of the three
points \(q^3=1\), each obtained from the finite-flag one-sided normal form
with parameters

\[
(-3q^2,-3q,0,-1,3q,-3q^2)
\]

in the ordered positions \(i>j\).  The last point factors, after the
pair-preserving change
\(A=W-V,B=W+V,P=Z+Y,Q=Z-Y\), as

\[
 -AP(AQ-BP)(AQ+BP)/4
\]

and has only the two negative-weight positions \((0,1),(2,3)\) in the
primed pair coordinates.  The RUR eliminant degrees are respectively
\(1\), \(d-c\), \(3\), and \(1\), so these displayed points exhaust the
complex zero sets, not merely their real loci.

This support statement is not invariant under the diagonal
\(\mathrm{SL}_2\)-action, but it removes every support of size at most
eight from later coefficient-torus calculations.

The first structured size-nine class is now closed as well.  Delete one
row and one column from the \(4\times4\) coefficient matrix.  Transpose
and simultaneous row/column reversal preserve all moments and split the
sixteen resulting \(3\times3\) rectangles into six orbits.  On a
representative dense coefficient torus, normalization leaves seven
coordinates.  Exact rational-univariate calculation of
\(\mu_1,\ldots,\mu_{14}\) gives a finite scheme.  Localizing that scheme
at each of the nine \(2\times2\) minors gives the unit ideal over
\(\mathbb Q\), in all \(6\mathbin{\cdot}9=54\) cases.  Therefore every
reduced dense point has coefficient rank one and is SIC-safe by the
balanced cubic rank-one theorem.  A rectangle boundary has at most eight
entries and is already safe by the complete smaller-support census.
Consequently:

> **Proposition 6.1.** Every bidegree-\((3,3)\) form supported in a
> \(3\times3\) coefficient rectangle is SIC-safe.  In particular none of
> the sixteen dense nine-entry rectangle tori contains a counterexample.

This is an exact characteristic-zero component calculation, not a
collection of modular fibres.  It closes sixteen of the \(11{,}420\)
mixed size-nine supports.  Support remains basis-dependent, so this does
not replace the diagonal-\(\mathrm{SL}_2\) chart classification.  The
second structured class closes 96 further supports.  Its row-fringe
members contain any two complete rows and one entry in a third row; the
column-fringe members are their transposes.  The 96 supports form 24
four-element transpose/reversal orbits.  Exact coefficient-torus saturation gives the
unit ideal through \(\mu_{14}\) on every representative, and every
boundary is covered by the support-eight theorem.  Hence:

> **Proposition 6.2.** Every bidegree-\((3,3)\) form supported on a
> complete-two-row/column fringe is SIC-safe.

The third class consists of one complete row, one complete column, and
any two entries outside their seven-entry cross.  There are 576 such
supports: 288 have the extra entries aligned in a secondary row or column,
and 288 have them in distinct secondary rows and columns.  They form 156
transpose/reversal orbits, 24 of size two and 132 of size four.  Exact
coefficient-torus saturation through \(\mu_{10}\) gives the unit ideal on
every representative.  Thus:

> **Proposition 6.3.** Every bidegree-\((3,3)\) form supported on a
> complete-row/complete-column cross plus two further entries is
> SIC-safe.

Six of the 576 cross supports lie entirely in one triangular half of the
coefficient matrix.  Thus Propositions 6.1--6.3 close 688 coordinate
subspaces, of which 682 belong to the mixed size-nine census.

The fourth class has three occupied rows, exactly three entries in each,
or the transposed column pattern; the sixteen \(3\times3\) rectangles are
removed.  There are 480 such mixed supports in 120 four-element symmetry
orbits.  Exact coefficient-torus elimination through \(\mu_{10}\) gives
the unit ideal on 114 representatives.  The remaining six systems each
have a degree-one rational RUR.  In the normalized coefficient order,
every residual point has matrix rank two.  For a rational flag parameter
\(q\), the contraction-preserving change

\[
 W'=W,\qquad V'=V-qW,\qquad Z'=Z+qY,\qquad Y'=Y
\]

sends its dehomogenized polynomial to a sum supported only at \(i>j\).
The relative-period integrand therefore has strictly positive
\(u\)-weight, and

\[
 \operatorname{CT}_u P(u,t)^m=0,\qquad
 \nu_{m+1}=0,\qquad \nu_1=0.
\]

For a balanced multiplier of bidegree \((e,e)\), its weight is at least
\(-e\), whereas every monomial in the \(m\)-th power has weight at least
\(m\).  Hence every mixed contraction vanishes for \(m>e\).  Exact checks
at multiplier degrees one and two also find nonzero low-order values, so
the cutoff is substantive rather than an identically zero mixed family.
Consequently:

> **Proposition 6.4.** Every bidegree-\((3,3)\) form supported on a
> regular three-row/column nine-entry support is SIC-safe.  The only
> nonempty dense moment systems through order ten are six unique rational
> rank-two points, and every one is fixed-flag one-sided.

The remaining three-line degree pattern is \(4+3+2\), together with its
transpose.  It contains 1,148 mixed supports in 287 four-element symmetry
orbits.  Exact coefficient-torus saturation through \(\mu_{10}\) gives
the unit ideal on every representative.  Its boundaries again have
support at most eight.  Hence:

> **Proposition 6.5.** Every bidegree-\((3,3)\) form supported on a
> nine-entry \(4+3+2\) three-row/column support is SIC-safe.

The row-count partitions \(4+4+1\), \(4+3+2\), and \(3+3+3\) exhaust
nine entries on three occupied rows.  Propositions 6.1, 6.2, 6.4, and
6.5 and their transposes therefore close all 1,740 nine-entry supports
having an empty row or column.  Including Proposition 6.3, the five
structured theorems close 2,310 mixed size-nine supports, leaving
\(9{,}110\) at this basis-dependent level.  The full discrete census has
\(2{,}924\) orbits under transpose and simultaneous reversal: 138 of size
two and 2,786 of size four.  The five propositions close 591 of these
mixed orbits and leave 2,333.  This is a support orbit census, not a
classification of continuous diagonal-\(\mathrm{SL}_2\) orbits.

There is one further complete-line family with row-count partition
\(4+3+1+1\) and no complete column, together with its transpose.  It has
1,244 mixed supports in 311 four-element symmetry orbits.  Exact dense
coefficient-torus elimination gives the unit ideal through \(\mu_{10}\)
on every representative, and support-eight handles every boundary.
Therefore:

> **Proposition 6.6.** Every bidegree-\((3,3)\) form supported on a
> full-line \(4+3+1+1\) nine-entry support, with no complete line on the
> other axis, is SIC-safe.

The six structured support theorems now close 3,554 mixed size-nine
supports in 902 discrete transpose/reversal orbits, leaving 7,866
supports in 2,022 orbits.  Every remaining support uses all four rows and
columns and belongs to one of five unordered row/column degree-partition
types.

The next global computation should saturate each normalized chart by the
nullcone and the explicit Rodrigues orbit before eliminating later
moments.  Reconstructing a specialized norm or a single finite-field
quotient remains useful
evidence, but cannot replace this residual component calculation.

On the exact-rank-two locus, the
[Hurwitz-chart calculation](TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_FINITE_PREFIX.md#4-two-exact-exclusions-on-the-generic-hurwitz-chart)
now removes the full \(\lambda=0\) fibre and the branch on which both
successive \(\mu_2\) pivots vanish.  The two principal-open branches
\(P_1\ne0\) and \(P_1=0,\ P_2\ne0\), the localization boundaries, other
channel charts, and exceptional cubic pencils remain.  The original
seven-equation system has mixed volume \(74\,144\); this is a complexity
measurement, not a root count or an exclusion theorem.

A complementary direct factor family uses reversal-even and reversal-odd
channel lines.  Before normalizing the internal gauge, write its two dual
channels as
\[
 \begin{aligned}
 A_+&=a(W^3+V^3)+b(W^2V+WV^2),\\
 A_-&=c(W^3-V^3)+d(W^2V-WV^2).
 \end{aligned}                                          \tag{6.0a}
\]
For \(\lambda\ne0\), put
\[
 \alpha=\frac{\lambda+\lambda^{-1}}2,\qquad
 \beta=\frac{\lambda-\lambda^{-1}}2,\qquad
 g_\lambda=\begin{pmatrix}\alpha&\beta\\\beta&\alpha\end{pmatrix}.
                                                               \tag{6.0b}
\]
Then \(\det g_\lambda=1\), and \(g_\lambda\) commutes with reversal.
In the eigen-coordinates \(X=W+V,Y=W-V\), it acts by
\(X\mapsto\lambda X,Y\mapsto\lambda^{-1}Y\).  The transformed endpoint
coefficients of (6.0a) are
\[
 \widetilde a=
 \frac{(a+b)\lambda^4+3a-b}{4\lambda},\qquad
 \widetilde c=
 \frac{(3c+d)\lambda^4+c-d}{4\lambda^3}.                \tag{6.0c}
\]
The coefficient matrices
\[
 \begin{pmatrix}1&1\\3&-1\end{pmatrix},\qquad
 \begin{pmatrix}3&1\\1&-1\end{pmatrix}
\]
both have determinant \(-4\).  Hence for every nonzero even/odd channel
pair the two bad endpoint conditions are proper affine-linear equations
in \(z=\lambda^4\).  Over the algebraic closure, a nonzero \(z\) avoids
both.  The inverse target action also commutes with reversal, and internal
diagonal gauge then normalizes \(U_{0,*}=(1,1)\).  Thus the normalized
chart below meets every exact-rank-two orbit in the parity family:
\[
 U=\begin{pmatrix}
 1&1\\ b&d\\ b&-d\\1&-1
 \end{pmatrix},\qquad
 W=\begin{pmatrix}
 p&q&q&p\\r&s&-s&-r
 \end{pmatrix}.                                         \tag{6.1}
\]
The parity decomposition fixes the two channel lines, and their internal
scalings are removed by the first row of \(U\).  On \(r=1\), the first
moment eliminates
\[
 p=-\frac{bq+ds+3}{3}.                                  \tag{6.2}
\]
The rows zero and three of \(U\) have determinant \(-2\), so \(U\) is
always rank two.  On this chart the odd row of \(W\) is nonzero, while
the even row is nonzero on the complete two-open cover
\[
 q\ne0\qquad\text{or}\qquad bq+ds+3\ne0.                \tag{6.3}
\]
Moreover
\[
 \operatorname {tr}\!\left((C\operatorname {diag}(6,2,2,6))^2\right)
 =32(ds+3)^2,                                           \tag{6.4}
\]
so inverting \(ds+3\) places the chart off the diagonal
\(\mathrm{SL}_2\)-nullcone.  On each open in (6.3), exact
characteristic-zero calculations in the variable orders
\((b,d,q,s)\) and \((s,q,d,b)\) agree:
\[
 \begin{array}{c|c}
  \mu_2,\ldots,\mu_5&\text{positive-dimensional},\\
  \mu_2,\ldots,\mu_6&\text{unit ideal}.
 \end{array}                                             \tag{6.5}
\]
The complementary projective chart has \(r=0,s=1\) and
\(p=-(bq+d)/3\).  Its exact-rank-two locus is likewise covered by
\(q\ne0\) and \(bq+d\ne0\).  On both opens the scheme is
zero-dimensional through \(\mu_5\) and is the unit ideal through
\(\mu_6\), again in two variable orders.  If \(r=s=0\), then \(W\) has
rank at most one.

It remains to classify the boundary \(ds+3=0\) on the \(r=1\) chart.
Put
\[
 d=-\frac3s,\qquad p=-\frac{bq}{3}.
\]
Here \(s\ne0\), and exact rank two is exactly \(q\ne0\).  After adjoining
\(hsq-1\), an exact characteristic-zero minimal-prime decomposition of
the ideal through \(\mu_6\) is
\[
 \begin{aligned}
  \mathfrak p_+&=(s-1,b+1,hq-1),\\
  \mathfrak p_-&=(s+3,b-3,3hq+1).
 \end{aligned}                                          \tag{6.6}
\]
Thus these are complete components over \(\mathbb Q(q)\), rather than
isolated modular fibres.

Both components are one-sided.  On \(\mathfrak p_+\), use
\((W',V',Z',Y')=(W,V-W,Z+Y,Y)\); on \(\mathfrak p_-\), use
\((W',V',Z',Y')=(W,V+W,Z-Y,Y)\).  The two forms become
\[
 \begin{aligned}
 F_+={}&(V')^2(Z')^2
 \left(\frac q3(2W'+V')Z'-V'(Z'-2Y')\right),\\
 F_-={}&(V')^2(Z')^2
 \left(-qV'(Z'+2Y')+(2W'-V')Z'\right).
 \end{aligned}                                          \tag{6.7}
\]
Every displayed monomial has dual-minus-target weight at least two.
Writing
\[
 P(u,t)=\sum_{i,j}c_{ij}u^{j-i}t^j(1-t)^{3-j},
\]
the relative period is
\[
 \frac{\mu_m}{(3m+1)!}
 =\operatorname {CT}_u\int_0^1P(u,t)^m\,dt.             \tag{6.8}
\]
Every term of \(P\) has \(u\)-degree at most \(-2\).  Hence the
function-field creative-telescoping recurrence degenerates to the exact
valuation telescoper
\[
 \nu_{m+1}=0\quad(m\ge0),\qquad
 \nu_m=\frac{\mu_m}{(3m+1)!},                            \tag{6.9}
\]
with forward coefficient one, no singular step, and initial vanishing
\(\mu_1=0\).

The degree-two multiplier \(M_{2,0,2}\) has relative-period numerator
\(u^2t^2\).  At \(m=1\) its mixed contractions are respectively
\[
 8(q+3),\qquad 24(1-q),                                  \tag{6.10}
\]
which are nonzero in \(\mathbb Q(q)\).  For an arbitrary degree-\(e\)
multiplier its numerator has \(u\)-degree at most \(e\), so the mixed
sequence vanishes whenever \(2m>e\).  The two components are therefore
SIC-safe.  Equations (6.0a)--(6.10), together with the \(r=0\) chart,
classify the complete exact-rank-two reversal-parity factor family up to
pair-preserving \(\mathrm{SL}_2\)-orbit and internal gauge.  The full
rank-two factor space remains open.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_sextic_slice.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_anchor_jacobians.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_boundary_family.py
```

The checker evaluates (2.5), verifies the exact nonzero minor (2.7), proves
the Hilbert obstruction (2.7c), and verifies the corrected rank-thirteen
minor for (2.7d). It then inverts the full Clebsch--Gordan basis to verify
(2.8)--(2.9), proves the diagonal radical and power certificates
(2.12)--(2.13), reconstructs all three lowering-chain embeddings, performs
the exact eliminations over \(\mathbb Q\), generates the restricted
moments, verifies the containments and power certificates (4.2)--(4.3),
verifies the direct ideal equalities (4.5)--(4.6), proves (5.4),
reproduces the finite-field calculation (5.5), and writes
`artifacts/generated-results/two_pair_sic_bidegree33_frontier.json`.
The second checker proves Proposition 5.2 and (5.9)--(5.12), then writes
the five exact Jacobian and pivot certificates. Modular reconnaissance of
the full chart ideals, which is evidence only, is run separately with

```bash
.venv/bin/python scripts/explore_two_pair_sic_bidegree33_full_anchor.py \
  --prime 43 \
  --orders 2,3,4,5,6,7,8,9,10,11,12,14 \
  --timeout 180 --backend msolve --charts s0
```

The exact Rodrigues, anti-Weyl, sparse, and null-chart certificates are
replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rodrigues_survivor.py

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_anti_weyl.py \
  --prime 0 --through 14 --backend msolve \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_anti_weyl_normalized_msolve14_char0.json

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_null_quadratic_s6.py \
  --orders 2,3,4,5,6,7,8,9,10,11 --skip-solver \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_null_quadratic_s6_local.json

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_survivor_rur.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_support8.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_rectangle9.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_two_row_fringe9.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_cross_two9.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_three_line9.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_three_line4329.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_full_line43119.py

.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_parity_channels.py
```

The two characteristic-zero Hurwitz exclusions are replayed by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz.py \
  --characteristic-zero --backend msolve --minor 01 \
  --lambda-value 0 --orders 2,3,4,5,6,7,8 \
  --timeout 60 --memory-gb 3 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_lambda0_char0.json

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz.py \
  --characteristic-zero --backend msolve --minor 01 \
  --orders 2,3,4,5,6,7,8 \
  --mu2-pivot-boundary-reduced secondary \
  --timeout 120 --memory-gb 3 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_secondary_boundary_char0.json
```

The sharded size-six and size-seven census commands are recorded in
[`REPRODUCE.md`](../REPRODUCE.md).
