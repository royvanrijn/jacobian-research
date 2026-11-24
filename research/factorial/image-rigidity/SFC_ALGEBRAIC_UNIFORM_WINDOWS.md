# Uniform algebraic windows and rational block-average SFC

Internal research proof, 21 November 2026. Inspected base commit:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

**What is proved.** Part A gives a fixed finite moment window at every
starting order for every monomial-times-linear factorial phase, including
zero exponents. Part B sharpens the bound when the auxiliary primitive is
polynomial, and proves a general trace-centered polynomial-observable
estimate. Part C uses that estimate to prove the actual support-size SFC
window for an unbounded family of rational block-average exponents.

**What remains open.** The general fixed windows are larger than the
support size. Neither the sharp contact bound for every monomial-linear
phase nor unrestricted SFC is proved. No literature-priority, external
refereeing, or formal-verification claim is made. All coefficient-field
statements cover arbitrary characteristic-zero fields. Written arguments
supply universal quantifiers; exact computations check their algebra.

## A. A fixed window for every monomial-linear factorial frame

### A.1. Statement and normalization

Let `s>=1`, let `r_i>=0` be integers, and let `c_i` be nonzero elements of
a characteristic-zero field. Equal coefficient values are allowed. Set

\[
 F=\left(\prod_{i=1}^s X_i^{r_i}\right)\sum_{i=1}^s c_iX_i,
 \quad R=\sum_i r_i,\quad d=R+1,
\]

\[
 P(z)=\prod_i(1-c_i z),\quad A(z)=\prod_i(1-c_i z)^{r_i},
 \quad a(z)=zA(z),\quad
 D(z)=P(z)+zP(z)A'(z)/A(z).
 \tag{A.1}
\]

The polynomial `a` has degree `d`, and `D` has degree `s`, constant
term one, and leading coefficient `(R+1) lc(P)`. Define the normalized
moments and their generating series by

\[
 g_m=\frac{\Lambda(F^m)}{m!\prod_i(mr_i)!}
     =[z^m]A(z)^{-m}P(z)^{-1},
 \qquad \rho(t)=\sum_{m\ge0}g_mt^m.
 \tag{A.2}
\]

Let `H(0)=0,H'=A/P` be a formal primitive, and let `B=H composed with
a^{-1}`. Formal Lagrange inversion gives

\[
 g_m=(m+1)[t^{m+1}]B(t),\qquad
 \rho=B',\qquad \rho(a(z))=\frac1{D(z)}.
 \tag{A.3}
\]

No polynomiality of `H` is assumed. In particular this includes mixed
zero and positive common exponents.

**Theorem A.** Put

\[
 e=s-\left\lceil\frac{s}{d}\right\rceil,
 \qquad
 C_{s,R}=d e+\frac{s d(d+1)}2.
 \tag{A.4}
\]

For every `n>=1`, at least one of

\[
                 g_n,g_{n+1},\ldots,g_{n+C_{s,R}-1}
 \tag{A.5}
\]

is nonzero. Equivalently, for every polynomial `b` of degree at most
`n>=1`,

\[
                 \operatorname{ord}_0(H-b\circ a)
                =\operatorname{ord}_0(B-b)
                \le n+C_{s,R}.
 \tag{A.6}
\]

All discarded factorials are nonzero integers, so (A.5) holds for the
ordinary factorial moments as well. If `R=0`, the sharper window given
by the number of distinct values among the `c_i` is valid.

The proof below supplies a more precise bound using the actual singular
values and the actual order of a scalar differential equation. Formula
(A.4) is a coefficient-independent upper bound for that quantity.

### A.2. Polynomial representation of the rational observable

First work over the complex numbers. The final argument descends to
every characteristic-zero field. Put

\[
 \Delta(t)=\operatorname{Res}_z(a(z)-t,D(z)).
 \tag{A.7}
\]

This is a nonzero polynomial of degree exactly `s`: up to a nonzero
constant it is the product of `a(alpha)-t` over the `s` roots of `D`,
counted with multiplicity. Its leading coefficient cannot cancel.

The polynomial `D(z)` divides `Delta(a(z))`. Indeed, if `alpha` is a
root of `D` of multiplicity `m`, the corresponding factor
`(a(alpha)-a(z))^m` already vanishes at least `m` times at `alpha`.
Thus

\[
                  Q(z)=\frac{\Delta(a(z))}{D(z)}
 \tag{A.8}
\]

is polynomial of degree exactly `s(d-1)`. Reducing it modulo
`a(z)-t` gives a unique polynomial vector `U(t)` such that

\[
 \frac1{D(z)}=\frac1{\Delta(t)}
               \sum_{i=0}^{d-1}U_i(t)z^i
 \quad\text{in }\mathbb C(t)[z]/(a(z)-t).
 \tag{A.9}
\]

Assign weights `wt(z)=1, wt(t)=d`. Polynomial reduction does not
increase weighted degree, since the relation replacing `z^d` has
terms of weights at most `d`. Consequently

\[
 \deg U_i\le\left\lfloor\frac{s(d-1)-i}{d}\right\rfloor,
 \quad
 u:=\max_{U_i\ne0}\deg U_i
       \le e=s-\left\lceil\frac{s}{d}\right\rceil.
 \tag{A.10}
\]

This is a polynomial division argument. It needs neither an inverse
Vandermonde asymptotic nor a genericity assumption.

### A.3. The connection and the actual possible singular values

For `d>=2`, let `sigma(t)` be the squarefree polynomial whose roots are
the distinct finite critical values of `a`, and put `c=deg sigma`.
Normalization of its nonzero leading coefficient is immaterial. At a
critical point of local degree `v`, the polynomial `sigma(a(z))`
vanishes at least `v` times, while `a'` vanishes `v-1` times. Therefore

\[
                        \frac{\sigma(a(z))}{a'(z)}
 \tag{A.11}
\]

is polynomial, of degree `(c-1)d+1`. In the basis
`1,z,...,z^(d-1)`, the derivation with respect to `t` is

\[
               \partial_t=\partial_t^{\mathrm{coeff}}
                                  +\frac{M(t)}{\sigma(t)},
               \qquad \deg M\le c-1.
 \tag{A.12}
\]

To see the degree directly, differentiate a basis monomial: its
numerator is `i z^(i-1) sigma(a)/a'`, of degree at most `cd-1`.
Reduction modulo `a-t` gives coefficient degrees at most `c-1`.
For `d=1`, take `sigma=1` and `M=0`; there is just the basis vector 1.

The actual frame relation `a'=(A/P)D` gives the useful support bound.
Every critical value other than zero is the image of a root of `D`.
Indeed, at a critical point where `D` is nonzero, `A/P` must have a
zero, so `A` and then `a` vanish. Hence

\[
               \sigma\mid\operatorname{rad}(t\Delta).
 \tag{A.13}
\]

Set

\[
              \Omega=\operatorname{rad}(\Delta\sigma),
              \qquad C=\deg\Omega.
 \tag{A.14}
\]

It follows that `C<=s+1`. If `Omega(0)!=0`, then `C<=s`. Both
`Omega Delta'/Delta` and `(Omega/sigma)M` are polynomial, of degree at
most `C-1`. The multiplicities of the poles of the observable are
handled by `Delta`; only their locations enter `Omega`.

### A.4. A degree-controlled scalar differential equation

Define polynomial vectors `N_k` by

\[
 \partial_t^k\rho=\frac{N_k}{\Delta\Omega^k},
 \qquad N_0=U.
 \tag{A.15}
\]

Here equality means equality in the function field with its displayed
basis, and then equality on the distinguished inverse branch. Exact
differentiation gives

\[
 N_{k+1}=\Omega N_k'
 +\left(\frac{\Omega}{\sigma}M
        -\Omega\frac{\Delta'}{\Delta}I-k\Omega'I\right)N_k.
 \tag{A.16}
\]

Every entry is polynomial and

\[
                       \deg N_k\le u+k(C-1).
 \tag{A.17}
\]

Let `q` be the smallest order of a nonzero rational-coefficient linear
differential equation satisfied by `rho`. It exists with `1<=q<=d`.
The columns `N_0,...,N_(q-1)` are independent over `C(t)`, and adjoining
`N_q` gives rank `q`. Choose `q` rows on which the first `q` columns
are independent. Their signed maximal minors `m_k`, omitting column
`k`, give the full vector relation `sum m_k N_k=0`, with `m_q!=0`.
All other rows satisfy the same relation because the full column rank
is `q`. The degree estimate is

\[
 \deg m_k\le q u+(C-1)\left(\frac{q(q+1)}2-k\right).
 \tag{A.18}
\]

Cancelling the common `Delta`, the scalar operator

\[
 \mathcal L=\sum_{k=0}^q p_k(t)\partial_t^k,
 \qquad p_k=m_k\Omega^k,
 \qquad
 M_*=q u+(C-1)\frac{q(q+1)}2
 \tag{A.19}
\]

has order `q`, annihilates `rho`, and satisfies

\[
                \deg p_k-k\le M_*.
 \tag{A.20}
\]

If `Omega(0)=0`, it also satisfies `ord_0 p_k>=k`.

### A.5. Why this minimal operator has no nonzero polynomial solution

The solution space of the minimal rational differential equation of an
algebraic function is the complex linear span of its algebraic
conjugates. Here is the precise argument used. Choose a maximal set
of conjugates linearly independent over the constants. Their Wronskian
is nonzero. The monic Wronskian operator has exactly their span as its
local solution space; its coefficients are fixed by every automorphism
of the finite normal closure of `C(t)`, because such an automorphism
permutes the span. Its coefficients therefore lie in `C(t)`. Conversely,
any rational operator annihilating one branch annihilates every
conjugate and has order at least the dimension of their constant span.
Thus that dimension is exactly the minimal order `q`. The cofactor
operator in (A.19), having the same minimal order, has the same solution
space. Multiplying an operator by a nonzero rational function does not
change that space away from its finitely many singular points.

In the present problem every root of `a(z)=t` tends to infinity at
order `t^(1/d)` along a sector at infinity. Since `deg D=s>=1`, every
conjugate of `1/D(z)` is `O(t^(-s/d))`. Every constant linear combination
therefore tends to zero. No nonzero polynomial can be such a solution.
This excludes constants as well as nonconstant polynomials; no
positivity premise is used.

### A.6. Contact, zero windows, and the uniform bound

Let `b` have degree at most `n>=1`. If it is constant, its contact with
`B=t+O(t^2)` is at most one. Otherwise `b'` is a nonzero polynomial,
so Section A.5 gives `L(b')!=0`. Write

\[
                         K=\operatorname{ord}_0(B-b).
\]

Since `L rho=0`, one has
`L(b')=-L((B-b)')`. Formula (A.20) gives

\[
                         \deg\mathcal L(b')\le n-1+M_*.
 \tag{A.21}
\]

If `Omega(0)=0`, each term on the right has order at least `K-1`,
because differentiation loses at most its order and `p_k` restores
that order. If `Omega(0)!=0`, the loss is at most `q`, so the lower
bound is `K-1-q`. Comparing the degree of a nonzero polynomial with
its vanishing order yields

\[
 \boxed{K\le n+L_*,\qquad
 L_*=q u+(C-1)\frac{q(q+1)}2
              +q\,\boldsymbol1_{\Omega(0)\ne0}.}
 \tag{A.22}
\]

The derivative-order inequalities are valid even if `K` is smaller
than some differentiation order; negative lower bounds are harmless.
They also show `B` is not a polynomial, since its derivative could
otherwise be a nonzero polynomial annihilated by `L`.

If `Omega(0)=0`, use `C-1<=s`. Otherwise use `C-1<=s-1` and
`q<=q(q+1)/2`. In either case,

\[
 L_*\le q e+\frac{s q(q+1)}2
      \le d e+\frac{s d(d+1)}2=C_{s,R}.
 \tag{A.23}
\]

A vanishing block of `C_(s,R)` moments starting at `n` would make the
degree-`n` truncation of `B` have contact at least `n+C_(s,R)+1`,
contradicting (A.6). The same argument uses the smaller actual bound
`L_*` when that data is available.

When `R=0`, `a=z`, `rho=1/P`, `u=0`, `q=1`, and `Omega=rad(P)` has
degree equal to the number of distinct coefficient values. Its
constant term is nonzero. Equation (A.22) therefore gives exactly that
number as a valid window. This also follows directly from the
first-order equation `Omega rho' + (Omega P'/P)rho=0`.

Finally, all coefficients involved in any proposed counterexample lie
in a finitely generated characteristic-zero field, which embeds in the
complex numbers. The integer exponent profile and ordinary factorial
identities survive that embedding. A zero block in the original field
would remain a zero block over the complex numbers. This proves the
uniform statement over every characteristic-zero field.

### A.7. Status and comparison

The result is uniform in both support size and starting order, but its
window is exponent dependent. For example, the primitive profile
`(1,1,2)` has `s=3,d=5,e=2`, giving the coarse fixed window 55. The
polynomial-observable trace-centered argument gives the stronger fixed
window 30 for this example; the reduced-divisor derivative calculation
gives window `3n`. These can be combined by taking the smallest proved
bound. None proves the desired window 3 for every starting order.

The actual alignment enters at (A.13). Using only a polynomial quotient
`a'/H'` without the zero-fiber condition would not justify the
support-based estimate for `C`. The proof does not assume a scalar
recurrence of order `s`, nonsingular period observability, or stability
under arbitrary row merging.

The formal identities in (A.1)--(A.3) are those of the retained
[contact note](SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md).
Admissible common-factor descent from the retained
[pullback note](SFC_RAMIFIED_PULLBACK_CLOSURE.md)
preserves `B` and `rho`, so one may first apply that exact descent and
then use this bound on the resulting seed.

## B. Trace centering and polynomial observables

### B.1. A trace-centered polynomial-observable bound

Work first over the complex numbers. Let

\[
a(z)=z+O(z^2)\in\mathbf C[z],\qquad d=\deg a\ge2,
\]

and let \(\sigma(t)\) be the monic squarefree polynomial whose roots are
the finite critical values of \(a\). Write

\[
c=\deg\sigma,\qquad
\epsilon=\begin{cases}0,&\sigma(0)=0,\\1,&\sigma(0)\ne0.\end{cases}
\]

There is at least one finite critical value, so \(c\ge1\).
Let \(H\in\mathbf C[z]\) have degree \(h\), put

\[
m=\lfloor h/d\rfloor,\qquad
\tau(t)=\frac1d\operatorname{Tr}_{\mathbf C(z)/\mathbf C(a)}H(z),
\qquad Y(t)=H(a^{-1}(t))-\tau(t),
\]

and suppose \(H\notin\mathbf C[a]\), so \(Y\ne0\).
The trace \(\tau\) is polynomial of degree at most \(m\).
Let \(q\) denote the minimal order of a nonzero linear differential
operator over \(\mathbf C(t)\) annihilating \(Y\). Then

\[
1\le q\le d-1.
\]

**Theorem B.** If \(Q\in\mathbf C[t]\), \(\deg Q\le N\),
\(N\ge m\), and \(Q\ne\tau\), then

\[
\boxed{\operatorname{ord}_0\bigl(H(a^{-1}(t))-Q(t)\bigr)
\le N+qm+(c-1)\frac{q(q+1)}2+\epsilon q.}
\tag{B.1}
\]

If \(Q=\tau\), the contact is at most \(h\). Two applications below
make this exceptional case harmless:

* if \(h<d\), \(H'(0)=1\), and \(N\ge1\), then its contact is at most
  one, since \(\tau\) is constant;
* if \(H=z^b\) and \(N\ge b\), its contact is at most \(b\le N\).

#### B.1.1 A connection with controlled polynomial degree

The extension \(\mathbf C(z)/\mathbf C(t)\), \(t=a(z)\), has basis
\(1,z,\ldots,z^{d-1}\). Reduction modulo \(a(z)-t\) respects the
weighted degree with weights \(\mathrm{wt}(z)=1\) and
\(\mathrm{wt}(t)=d\). Thus a polynomial of degree \(v\) in \(z\)
reduces to coordinates of \(t\)-degree at most \(\lfloor v/d\rfloor\).

If \(\alpha\) is a critical point of local degree \(e\), then
\(a'(z)\) has order \(e-1\) at \(\alpha\), whereas
\(\sigma(a(z))\) has order at least \(e\). Consequently

\[
J_a(z)=\frac{\sigma(a(z))}{a'(z)}\in\mathbf C[z],
\qquad \deg J_a=(c-1)d+1.
\]

The derivative of the basis element \(z^j\) is

\[
\frac{\partial z^j}{\partial t}
=\frac{jz^{j-1}J_a(z)}{\sigma(t)}.
\]

For \(1\le j<d\), its numerator has degree at most \(cd-1\).
Reducing these numerators gives a matrix \(C(t)\) satisfying

\[
\partial_t v=v'+\frac{C(t)}{\sigma(t)}v,
\qquad \deg C\le c-1,
\tag{B.2}
\]

where \(v\) is a column of coordinates in the displayed basis.
The zeroth column of \(C\) is zero.

The reduction of \(H\) has coordinate degrees at most \(m\).
The traces of \(1,z,\ldots,z^{d-1}\) are constants, by the Newton
identities for \(a(z)-t\), whose nonconstant coefficient is its constant
term. This also proves \(\deg\tau\le m\). The reduced coordinate
vector \(N_0\) of \(H-\tau(a)\) therefore has degree at most \(m\).

Define polynomial vectors recursively by

\[
N_{k+1}=\sigma N_k'-k\sigma' N_k+C N_k.
\tag{B.3}
\]

They represent \(\partial_t^kY=N_k/\sigma^k\), and induction gives

\[
\deg N_k\le m+k(c-1).
\tag{B.4}
\]

#### B.1.2 Why the minimal operator has no rational solutions

Let \(z_1(t),\ldots,z_d(t)\) be the roots of \(a(z)-t\) in its
splitting field \(E\) over \(\mathbf C(t)\). The centered conjugates

\[
Y_i=H(z_i)-\tau(t)
\]

sum to zero. Their span over the constant field \(\mathbf C\) has
dimension \(r\le d-1\), and is invariant under the Galois group.
For any constant-linearly independent basis of that span, its Wronskian
is nonzero. The Wronskian criterion here follows, for instance, by
induction after dividing all functions by a nonzero first function and
differentiating the remaining ratios. The constant field of \(E\) is
\(\mathbf C\).

The monic order-\(r\) Wronskian operator is unchanged by a constant
change of basis. Its coefficients are therefore Galois invariant and
lie in \(\mathbf C(t)\). Hence \(q\le r\). Conversely, every rational
operator annihilating \(Y\) annihilates all its conjugates. On a small
ordinary disk an order-\(q\) equation has a \(q\)-dimensional solution
space, so \(r\le q\). Thus \(q=r\), and the solutions of the minimal
operator are precisely the constant span of these conjugates.

There is no nonzero rational function in that span. Indeed, averaging
over the finite Galois group fixes any rational function and kills each
\(Y_i\), because the action on the roots is transitive and their sum
is zero. Therefore the minimal operator has no nonzero rational
solution. This argument excludes polynomial solutions as well; it does
not assume that algebraic degree equals differential order.

#### B.1.3 Exact polynomial coefficients and contact

The first \(q\) coordinate vectors \(N_0,\ldots,N_{q-1}\) are
independent over \(\mathbf C(t)\), and \(N_q\) is dependent on them.
Choose \(q\) coordinate rows with a nonzero minor in the first \(q\)
columns. The signed maximal minors \(u_k\), \(0\le k\le q\), of
the resulting \(q\)-by-\((q+1)\) matrix give the relation

\[
\sum_{k=0}^q u_kN_k=0.
\]

The relation holds in every coordinate row, since the full matrix has
rank \(q\). Moreover \(u_q\ne0\), and (B.4) gives

\[
\deg u_k\le qm+(c-1)\left(\frac{q(q+1)}2-k\right).
\]

Put \(p_k=u_k\sigma^k\). The operator

\[
L=\sum_{k=0}^q p_k\partial_t^k
\]

annihilates \(Y\), has minimal order, and satisfies

\[
\deg p_k-k\le M:=qm+(c-1)\frac{q(q+1)}2.
\tag{B.5}
\]

If \(\sigma(0)=0\), it also satisfies
\(\operatorname{ord}_0p_k\ge k\). Keep these unreduced coefficients:
dividing out their common factor without tracking its order would
invalidate this particular valuation argument.

Let \(\widetilde Q=Q-\tau\ne0\). Its degree is at most \(N\), and
the previous subsection gives \(L\widetilde Q\ne0\). By (B.5),

\[
\deg(L\widetilde Q)\le N+M.
\]

Writing \(K=\operatorname{ord}_0(Y-\widetilde Q)\), the elementary
valuation bound on each summand of \(L(Y-\widetilde Q)\) gives

\[
\operatorname{ord}_0(L\widetilde Q)\ge K-\epsilon q.
\]

A nonzero polynomial has order at most its degree, proving (B.1).
If \(Q=\tau\), substitute \(t=a(z)\): the nonzero polynomial
\(H(z)-\tau(a(z))\) has degree at most \(h\), giving the stated
exceptional bound.

### B.2. Application when the factorial primitive is polynomial

Group equal nonzero coefficients into distinct values \(c_j\), of
multiplicity \(p_j\), and let \(R_j\) be the sum of the corresponding
nonnegative exponents. Set

\[
S=\sum_jp_j,\quad R=\sum_jR_j,\quad d=R+1,
\]
\[
P(z)=\prod_j(1-c_jz)^{p_j},\quad
A(z)=\prod_j(1-c_jz)^{R_j},\quad a=zA,
\]
\[
H(0)=0,\quad H'=A/P,\qquad
D=P+zP A'/A.
\]

Suppose \(R_j\ge p_j\) for every group. Then \(H\) is polynomial,
\(H'(0)=1\), \(\deg H=d-S<d\), and

\[
a'=H'D,\qquad \deg D=S.
\]

Every zero of \(H'\) maps to zero under \(a\). Thus the number of
nonzero critical values of \(a\) is at most \(S\). In Theorem B,
\(m=0\) and \(q\le d-1\). If zero is a critical value, \(c-1\le S\);
otherwise \(c\le S\), and

\[
(c-1)\frac{q(q+1)}2+q\le S\frac{q(q+1)}2.
\]

It follows that, for every polynomial \(b\) of degree at most \(n\ge1\),

\[
\boxed{\operatorname{ord}_0(H-b\circ a)
\le n+W_{\mathrm{poly}},\qquad
W_{\mathrm{poly}}=\frac{S\,d(d-1)}2.}
\tag{B.6}
\]

The exact coefficient identity from the monomial-linear contact proof is

\[
g_m:=\frac{\Lambda(F^m)}{m!\prod_i(mr_i)!}
=(m+1)[t^{m+1}]H(a^{-1}(t)).
\]

If \(W_{\mathrm{poly}}\) consecutive moments starting at \(n\) vanished,
the degree-\(n\) truncation would have contact at least
\(n+W_{\mathrm{poly}}+1\), contradicting (B.6). Hence every such window
contains a nonzero moment, uniformly in its starting order and in all
active coefficients. A smaller certificate is available by retaining
the actual \(q,c,\epsilon\) in Theorem B.

For any active coefficients and exponents \((1,1,2)\), this yields
the fixed window \(30\). The reduced-divisor estimate yields \(3n\):
only the group containing the exponent two can have nonzero exponent
excess, so its reduced logarithmic degree is at most one and its
critical-quotient radical has degree at most three. The proofs combine to give \(\min(3n,30)\). The desired SFC
window is \(3\), and is still unproved by these arguments for arbitrary
starting order. The number \(30\) is a guaranteed upper bound, not a
claim that a zero run of that length occurs.

### B.3. Arbitrary characteristic-zero coefficient fields

A hypothetical violation of any displayed bound uses only finitely many
coefficients, hence is defined over a finitely generated extension of
\(\mathbf Q\). Such a field embeds in \(\mathbf C\): choose an
algebraically independent complex image of a transcendence basis and
extend through its finite algebraic extension. The embedding preserves
all vanishing and nonvanishing statements. Applying the complex proof
contradicts the hypothetical violation. This transfers the conclusions
to every characteristic-zero field; no ordering or positivity on that
field is used.


## C. Rational block-average SFC under an explicit size condition

### C.1. An explicit fixed window for powers of a polynomial inverse

Let `A(z)=z+O(z^2)` be polynomial of degree `D>=2` over `C`, let
`V=A^{-1}`, and let `b>=1` be an integer. Put

\[
 m=\lfloor b/D\rfloor,\qquad
 \tau(t)=\frac1D\operatorname{Tr}_{\mathbb C(z)/\mathbb C(A)}z^b,
 \qquad Y(t)=V(t)^b-\tau(t).
 \tag{C.1}
\]

Let `sigma` be the monic squarefree polynomial of all finite critical
values of `A`, including zero when it is a critical value. Write

\[
 c=\deg\sigma,\qquad
 \epsilon=\begin{cases}0,&\sigma(0)=0,\\1,&\sigma(0)\ne0.\end{cases}
 \tag{C.2}
\]

Let `q` be the minimal differential order of `Y` over `C(t)`. Then
`1<=q<=D-1`, and `deg tau<=m`. Define the integer

\[
 C_{A,b}=qm+(c-1)\frac{q(q+1)}2+\epsilon q.
 \tag{C.3}
\]

**Theorem C1.** If `N>=b` and `deg Q<=N`, then

\[
             \operatorname{ord}_0(V^b-Q)\le N+C_{A,b}.
 \tag{C.4}
\]

For `0<=N<b`, the elementary bound is instead
`ord_0(V^b-Q)<=b`. Thus the assertion does not silently use the
differential estimate below the degree of the trace polynomial.

#### Proof, including the exceptional cases

Reduction modulo `A(z)-t` in the basis `1,z,...,z^(D-1)` turns a
degree-`b` polynomial into coordinates of degree at most `m` in `t`.
The traces of these basis elements are constants, by Newton identities,
so `tau` is polynomial of degree at most `m`.

The centered observation is nonzero. If `z^b=tau(A(z))`, evaluating at
zero gives `tau(0)=0`. Since `A` has a simple zero at zero and degree
at least two, it has a nonzero root `alpha`; evaluating there would give
`alpha^b=0`, a contradiction. Therefore `q=0` cannot occur in the stated
domain. For a general polynomial observation belonging to `C[A]`,
centering would give zero and this argument would not apply. Likewise
`D=1`, where `V^b=t^b`, is explicitly excluded.

For `N>=b`, Part B applies with polynomial observation `H=z^b`:
`N>=b>=m` ensures the centered approximant retains degree at most `N`.
If `Q!=tau`, its no-rational-solution argument gives (C.4). If `Q=tau`,
the nonzero polynomial `z^b-tau(A(z))` has degree at most `b`, so its
order is at most `b<=N`. Thus the exceptional case also satisfies (C.4).

If `N<b`, any nonzero `Q` has a first nonzero term below degree `b`,
while `Q=0` has contact exactly `b`. This proves the small-degree clause,
including `N<m`. The connection, trace, cofactor degrees and complete
minimal-operator proof are given in Part B. QED.

### C.2. The ordinary factorial family

Let `a,b,v_1,...,v_h` be positive integers. Write

\[
 \nu=\sum_{j=1}^h v_j,\qquad D=1+a\nu,\qquad S=b\nu.
 \tag{C.8}
\]

For each group `j`, use `p_j=bv_j` independent ordinary variables
`Y_(j,l)`, and choose any nonnegative integer exponents satisfying

\[
             \sum_{l=1}^{p_j} r_{j,l}=av_j.
 \tag{C.9}
\]

Let `c_1,...,c_h` be distinct nonzero field elements and set

\[
 F=\left(\prod_{j,l}Y_{j,l}^{r_{j,l}}\right)
                 \sum_{j=1}^h c_j\sum_{l=1}^{p_j}Y_{j,l}.
 \tag{C.10}
\]

Its actual monomial support size is `S`: the linear terms involve
different variables, so multiplying by the common monomial does not
identify them. Put

\[
 R(z)=\prod_j(1-c_jz)^{v_j},\qquad A_*(z)=zR(z^a).
 \tag{C.11}
\]

**Theorem C2.** If

\[
 \boxed{\quad 2\bigl(b-\lfloor b/D\rfloor\bigr)\ge ahD,\quad}
 \tag{C.12}
\]

then every window of `S` consecutive ordinary factorial moments of
`F`, beginning at any positive order, contains a nonzero moment.
The coefficients, allocations (C.9), support sizes and starting orders
are otherwise unrestricted within the statement.

A simpler sufficient condition is

\[
             b\ge\left\lceil\frac{ahD^2}{2(D-1)}\right\rceil
              =\left\lceil\frac{hD^2}{2\nu}\right\rceil.
 \tag{C.13}
\]

If every `v_j=1`, the sharper sufficient test

\[
           2\bigl(b-\lfloor b/D\rfloor\bigr)
                      \ge(ah-1)D+2
 \tag{C.14}
\]

is available. A single active coefficient group already has every
moment nonzero by the binomial formula, without these inequalities.

#### Proof

Direct expansion, or convolution of the equal-coefficient variables,
gives the normalized ordinary moment

\[
 g_m=\frac{\Lambda(F^m)}{m!\prod_{j,l}(mr_{j,l})!}
                =[z^m]R(z)^{-(am+b)}.
 \tag{C.15}
\]

Lagrange inversion for `V_*=A_*^{-1}` gives

\[
        [t^{am+b}]V_*(t)^b=\frac{b}{am+b}g_m.
 \tag{C.16}
\]

Only degrees congruent to `b` modulo `a` occur in this inverse power.
All displayed scalar factors are nonzero in characteristic zero.

The nonzero critical values of `A_*` number at most `ah`: after
discarding derivative roots over zero, its derivative equation has
degree at most `h` in `z^a`. If some `v_j>1`, zero is a critical value,
so `c-1<=ah` and `epsilon=0`. If all `v_j=1`, zero is not a critical
value, so `c<=ah` and `epsilon=1`. With `q<=D-1` in Theorem C1, both
cases give the uniform bound

\[
 C_{A_*,b}\le C_*:=(D-1)\lfloor b/D\rfloor
                         +\frac{ahD(D-1)}2.
 \tag{C.17}
\]

For the zero-free critical fiber case, the sharper bound is

\[
 C_{A_*,b}\le C_*^{(0)}:=(D-1)\lfloor b/D\rfloor
                 +\frac{(ah-1)D(D-1)}2+(D-1),
 \tag{C.18}
\]

which is at most (C.17) because `D>=2`.

Suppose that `g_n,...,g_(n+S-1)` all vanished. Truncate `V_*^b` at

\[
                        N=a(n-1)+b\ge b.
 \tag{C.19}
\]

The next possible nonzero degree would be at least
`N+a(S+1)`. Condition (C.12) is exactly the sufficient inequality

\[
                  C_*\le b(D-1)=aS.
 \tag{C.20}
\]

Thus (C.16) would force contact strictly greater than `N+C_*`,
contradicting Theorem C1. The same argument with (C.18) proves (C.14).
Finally, `b-floor(b/D)>=b(D-1)/D` proves (C.13). QED.

The argument also gives the explicit, sometimes smaller window
`floor(C_*/a)` for this family, or `floor(C_*^(0)/a)` when all
`v_j=1`. The theorem uses the support-size consequence, not an
assumption that this bound is sharp.

### C.3. Coefficient boundaries, examples and retained-scope comparison

The theorem is stated on the active coefficient groups. Equal values
can be regrouped by adding their `v_j`; their exponent totals and
variable counts still satisfy (C.9). This lowers `h` without changing
`nu`, so (C.12) remains sufficient. Zero coefficient groups may be
discarded; their common-monomial variables contribute only nonzero
factorial scalars. The criterion then uses the remaining active data.

For a family asserted simultaneously on all coefficient faces, (C.13)
is a convenient uniform sufficient hypothesis: the function
`D^2/(D-1)` is nondecreasing for integer `D>=2`, and deleting groups
decreases both `h` and `D`. Thus its threshold does not increase.
The completely zero polynomial is excluded, as in SFC.

For a small example take `a=1,b=4,v_1=v_2=1`. Then `D=3,S=8`, and
`2(4-floor(4/3))=6=ahD`. Consequently

\[
 F=X_1Y_1\bigl(c(X_1+X_2+X_3+X_4)
                   +d(Y_1+Y_2+Y_3+Y_4)\bigr)
 \tag{C.21}
\]

satisfies SFC at every starting order for arbitrary coefficients
`c,d` not both zero. If `c=d` or one is zero, the single-group binomial
formula handles that boundary directly. In fact (C.18) gives a window
of seven for the two distinct active groups.

For arbitrary `h>=2`, take `a=1`, every `v_j=1`, and
`b=(h+1)^2+1`. Condition (C.12) holds, and this gives a family with
`h b` monomials and arbitrary active group coefficients. Its total
common exponent is `h`; moreover

\[
                 \gcd(hb,h+1)=1.
 \tag{C.22}
\]

Thus it cannot be obtained by a nontrivial support-multiplying
polynomial pullback of a smaller canonical factorial frame: the
retained degree obstruction requires the pullback degree to divide
`gcd(S,1+sum r_i)`. This is an infinite primitive family in that
precise sense; it is not the still-open positive-exponent primitive
core.

These examples have common group-average exponent `a/b`, generally
nonintegral. The retained weighted-radial theorem requires the
integer ratio `R_j/p_j=q`; it does not cover them. Independent linear
padding does not directly produce (C.21), since all added terms carry
the factor `X_1Y_1`, which uses variables present in the linear form.
This comparison is with the inspected retained transfer theorems,
not a claim of exhaustive novelty in the literature.

The new closure is for coefficient groups constrained to have equal
values within each group. It is not a theorem for fully independent
coefficients on all `S` monomials. The group-average and large-
denominator conditions must not be dropped.

### C.4. Fields and the sharp target still open

All coefficient formulas are algebraic in the finitely many parameters.
For a characteristic-zero coefficient field, embed its finitely
generated coefficient subfield into `C`. A putative zero block would
be preserved, contradicting the complex proof. This proves the same
ordinary factorial statement over every characteristic-zero field.
The analytic language in the minimal-operator argument therefore
imposes no ordering or positivity assumption on the original field.

The earlier target

\[
 \operatorname{ord}(V^b-Q)
             \le \deg Q+b\,\#\operatorname{CV}^{\times}(A)
 \tag{C.23}
\]

remains unproved in general. Theorem C1 supplies a larger explicit
constant. Theorem C2 succeeds where the actual support size is large
enough to absorb that constant; it does not establish (C.23) elsewhere.
The arbitrary unequal positive-exponent monomial-linear sharp window,
and the further transfer from that family to general monomial
supports, remain separate obligations.

The companion standard-library checker reconstructs the connection,
trace, numerator recurrence and cofactor annihilator for stated exact
examples, as well as (C.15)--(C.16) from ordinary factorials. It tests the
zero-center and small-approximation-degree controls. These computations
check the displayed identities; the preceding proofs carry all
unrestricted quantifiers.

## D. Literature, remaining target, and replay boundary

### D.1. Relation to algebraic gap methods

This is a controlled-degree differential-operator argument. It does not
use a point-counting or finite-field Stepanov estimate. Dutter's
*Gaps in Taylor series of algebraic functions*, Bull. Austral. Math.
Soc. 91 (2015), 412–418, introduction on p. 413, describes earlier
Chudnovsky–Chudnovsky methods giving effectively computable, potentially
large bounds on consecutive gaps. That context is consistent with the
role of the present calculation; no priority claim is made. Dutter's
Theorem 2.1 itself bounds the position of the successive nonzero terms
and is not imported as a sharp support-size window theorem.

Primary source inspected:
<https://www.cambridge.org/core/services/aop-cambridge-core/content/view/01445BE33E476E479BE182522E79FBFA/S0004972715000040a.pdf/div-class-title-gaps-in-taylor-series-of-algebraic-functions-div.pdf>.

The proof above supplies its own connection, annihilator-degree bounds,
and exclusion of rational solutions. No RH or Alpöge claim is used.

### D.2. Scope still requiring a sharper argument

The support-size target remains
\(\operatorname{ord}_0(H-b\circ a)\le\deg b+S\).
The coefficient minors in (B.5) establish a fixed bound but do not reduce
it to \(S\). Common polynomial pullback factors may transfer a smaller
seed bound, but primitive frames remain. Even a sharp theorem for all
monomial-times-linear polynomials would still require an additional
valid transfer to arbitrary monomial supports to prove full SFC.

### D.3. Exact replay

Run from the repository root:

```sh
python research/factorial/image-rigidity/verify_sfc_uniform_rational_window.py
python research/factorial/image-rigidity/verify_sfc_trace_centered_window.py
python research/factorial/image-rigidity/verify_sfc_rational_block_sfc.py
```

The three standard-library programs pin this complete proof and their
retained arithmetic dependencies. The rational-observable checker verifies
six actual frames, their complete resultant identities, weighted reductions,
connections, minimal cofactor equations, and 54 ordinary factorial moments.
The trace-centered checker verifies four exact connections and minimal
operators, 11 determinant values establishing complete critical-value
resultant identities, 48 truncation contacts, and 39 normalized factorial
moments. The rational-block checker records its exact inverse-power and
ordinary-factorial controls in its receipt. None treats its finite cases as
proof of the unrestricted quantifiers established above.

Companion results:
[SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md),
[SFC_REDUCED_DIVISOR_CONTACT_BOUND.md](SFC_REDUCED_DIVISOR_CONTACT_BOUND.md),
and [SFC_INVERSE_POWER_CONTACT_AUDIT.md](SFC_INVERSE_POWER_CONTACT_AUDIT.md).
The last note records the sharp inverse-power target and an exact failed
critical-value budget; Part C here does not assume that target.
