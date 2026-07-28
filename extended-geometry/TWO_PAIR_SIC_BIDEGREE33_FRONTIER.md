# The two-pair bidegree-\((3,3)\) frontier

## 1. Status

This note starts the first genuinely new two-pair bidegree after the
complete bidegree-\((2,2)\) theorem.  Unequal bidegrees are already
one-sided for the total dual-minus-coordinate grading, so the next balanced
case is \((3,3)\).

Ten exact characteristic-zero results are proved here:

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

The fourth result excludes every SIC(2) counterexample lying in a single
irreducible summand.  The sixth gives dimension-sized moment coordinates,
not the required zero-fiber theorem.  The full mixed
\(\operatorname{Sym}^6\oplus\operatorname{Sym}^4\oplus
\operatorname{Sym}^2\) problem remains open.

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
`msolve`, after direct submission, the first pivot, and the boundary
\(A=B=0\). Export of the fully substituted \(A\ne0\) branch produces
eleven equations in ten variables, but the recorded `msolve` run
terminated inside the solver. Sparse encodings retain \(s_5\) or \(t_4\)
and add the inverse-pivot relation, avoiding the expansion; the full
corrected \(A\)- and \(B\)-open systems nevertheless exceed the recorded
bounds. These are not mathematical results. The next computation should
construct the finite special quotient of \(\mu_2,\ldots,\mu_{12}\) on
each stratum and test \(\mu_{14}\) there, while decomposing the common
boundary beyond (5.13).

## 6. Generalized continuation

Further full-coefficient Gröbner elimination is no longer the primary
attack. The
[moment--nullcone program](TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md)
places (5.4)--(5.5) inside the decomposition

\[
 \operatorname{End}(\operatorname{Sym}^d)
 \cong\bigoplus_{r=0}^d\operatorname{Sym}^{2r}
\]

for arbitrary \(d\). Its first target is a global invariant
quadratic-anchor certificate: all moments should force the discriminant of
the \(\operatorname{Sym}^2\) component to vanish. The next target is a
common-root synchronization identity forcing the higher binary-form
components into the same destabilizing flag.

Reconstructing \(c^{25}\) over \(\mathbb Q\) remains a legitimate
specialized check, but by itself it does not supply either general
identity.

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
