# Generic Bernstein equations and actual formal-symbol gaps

Internal written proof, 16 November 2026. The mathematical authority is
[MATH_STATUS.json](../../MATH_STATUS.json). Claims
**GVCFORMALORDERPOLY20261116** and **GVCFORMALALLLINEAR20261116**.
Independent replay, formal verification and external review are all **false**.
The exact checker exercises transcription controls; the universal proof is written below.


## Target

Let $K$ be any characteristic-zero field, $P\in K[x_1,\ldots,x_n]$,
$\lambda\in (\xi)K[[\xi]]$, and $T=\lambda(D_x)$.
Assume $T^mP^m=0$ for every $m\geq1$, and that for every polynomial $Q$
the sequence $T^m(QP^m)$ is eventually zero. The target is an integer
$C$ such that $T^m(QP^m)=0$ whenever $m>C\deg Q$.
Zero $P$ or zero $T$ are immediate; constants are also immediate because
$T$ has positive order. Below $P$ and $\lambda$ are nonzero.

## Primary import, with its formal hypotheses

Mebkhout–Narváez-Macarro, Ann. sci. ENS (4)24 (1991),227--256,
[DOI 10.24033/asens.1627](https://doi.org/10.24033/asens.1627), Sections 1.1.2, 1.2, 2 and 3:
for a regular equicodimensional characteristic-zero algebra whose closed
residue fields are algebraic over the ground field and whose derivations
have coordinate bases, minimal-dimension differential modules have finite
length. Theorem3.1.1 supplies the Bernstein equation for an element of
such a module and an arbitrary coefficient-ring element. Formal power
series rings satisfy these hypotheses. Remark2.2.5 also gives preservation
of minimal dimension under local coefficient completion. The proof of
Theorem3.1.1 supplies the generic-power minimal lattice used below.
No analytic-only specializability theorem is imported.

### The actual Fourier module

Set $R=K[[\xi_1,\ldots,\xi_n,\eta]]$ and
$\mathcal D=R\langle\partial_{\xi_1},\ldots,\partial_{\xi_n},\partial_\eta\rangle$.
On $E=K[x][[\eta]]$ let
$$
 \xi_i=D_{x_i}+\eta P_i(x),\quad \partial_{\xi_i}=-x_i,
 \partial_\eta=D_\eta+P(x).
$$
These actions define a $\mathcal D$-module. Its actual cyclic image
$M=\mathcal D1\subset E$ is holonomic and $\eta$-torsion-free.
For $v_k=\lambda(\xi)^k1$,
$$
 v_k=e^{-\eta P}T^ke^{\eta P},\qquad
 \nu_\eta(v_k)=\tau_k:=\min\{j\geq0:T^kP^j\ne 0\}.
$$
For $a\in R$, its action on the leading time coefficient is $a(D_x,0)$.
If $a(0,0)\ne 0$, this last operator is invertible on $K[x]$.

**Proof.**
All Weyl relations follow by direct commutators. With $D=\max(1,\deg P)$,
a word of length $N$ in the $\xi_i$ applied to a degree-$q$ polynomial
has zero coefficient at time degree $j$ if $N>q+Dj$:
a term using $j$ multiplication factors has degree at most $q+Dj-N$.
This proves coefficientwise well-definedness for every formal coefficient
and every input of $E$, without completing the multiplier space in an
eventuality assumption. Formal coefficients preserve $\eta$-order;
$\partial_\eta$ preserves $E$, though it can lower positive order.

Start with the algebraic rank-one exponential connection
$K[x,\eta]e^{\eta P}$, a holonomic algebraic Weyl module. Partial Fourier
transform in $x$ preserves holonomicity. Extend its coefficients from
$K[\xi,\eta]$ to $R$. A good order filtration and flat coefficient completion
show that the resulting module is holonomic or zero. Sending its generator
to $1\in E$ is well-defined: every algebraic relation holds under the displayed
actions, and the tensor extension maps by $a\otimes m\mapsto a m$.
Its image is exactly $M$, not a larger space of solutions. Holonomic modules
are closed under quotients. The embedding in $E$ gives time torsion-freeness.
The conjugation identity follows first on polynomials by the product rule
and then coefficientwise by the word-length bound. The exponential is a
time-adic unit, proving the order formula. Finally a zero-constant formal
differential operator lowers polynomial degree, so the formal inverse of
$a(D_x,0)$ acts finitely on each polynomial.
∎

### Generic power and one further Bernstein parameter - GVCFORMALORDERPOLY20261116

There exist $b(s,z)\in K[s,z]\setminus\{0\}$ and an integer $k_0$ such that
for every integer $k\geq k_0$ with $v_k\ne 0$,
$$
                         b(k,-\tau_k-1)=0.             \tag{1}
$$

**Proof.**
First consider $\mathcal D(s)=K(s)\otimes_K\mathcal D$ and the symbol
$w=\lambda^s1$ in $M_\lambda(s)\lambda^s$, with
$\partial_i(\lambda^s m)=\lambda^s(\partial_i m+s\lambda_i m/\lambda)$;
the time derivative has no correction. The proof of the primary
Theorem3.1.1, using Proposition1.2.7 and Theorems2.1.4,2.2.3, provides a
minimal-dimension $\mathcal D(s)$-lattice $L$ whose localization is this
twisted module. Thus $\lambda^\ell w\in L$ for some $\ell\geq0$.
The cyclic module $\mathcal D(s)\lambda^\ell w$ is minimal dimensional.
The semilinear parameter translation
$$
 a(s)m\lambda^s\longmapsto a(s-\ell)\lambda^{-\ell}m\lambda^s
$$
intertwines $\mathcal D(s)$ with its scalar automorphism $s\mapsto s-\ell$
and sends $\lambda^\ell w$ to $w$. It preserves the order filtration and
characteristic dimension. Consequently $N=\mathcal D(s)w$ itself has
minimal dimension. This argument does not assume that generic twisting
commutes with an unbounded parameter specialization.

Put $B=(R[s])_{(\xi,\eta)}$, localizing at the prime generated by the
coordinate variables, and $\widehat B=K(s)[[\xi,\eta]]$.
This is the maximal-ideal completion of the regular local ring $B$;
its residue field is $K(s)$, and $B\to\widehat B$ is faithfully flat.
The localized cyclic module $N_B=\mathcal D_Bw$ is holonomic: a localized
good order filtration has dimension at most $n+1$, and the Bernstein
inequality for $B$ gives the reverse inequality if it is nonzero.
Remark2.2.5 now gives holonomicity of
$\widehat N=\mathcal D_{\widehat B}w$ over ground field $K(s)$.
This is coefficient extension of a finite cyclic module, not replacement
of $M$ by its time completion $E$.

Apply Theorem3.1.1 to $\widehat N$, the coefficient function $\eta$,
and the new exponent $z$. It gives
$$
 b(s,z)\eta^z w=Q(s,z)\eta^{z+1}w,                    \tag{2}
$$
where $0\ne b\in K(s)[z]$ and
$Q\in\mathcal D_{\widehat B}[z]$, of finite differential order and finite
$z$-degree. Clear the scalar $s$ denominators of $b$.

Equation(2) descends to $B$. Indeed use the ambient
$B$-module $M_B[\lambda^{-1},\eta^{-1}][z]\lambda^s\eta^z$.
The image of $\mathcal D_B[z]\eta w\eta^z$ extends to the corresponding
image over $\widehat B$ by PBW and flatness. Faithful flatness applied
to the quotient by that image descends the membership of $b(s,z)w\eta^z$.
It supplies an operator in $\mathcal D_B[z]$, not a coefficientwise
infinite series in the exponent. Left multiply by one common denominator
$g(s,\xi,\eta)\in R[s]$ with $g(s,0,0)\ne 0$. The result has an operator
$Q_0\in\mathcal D_R[s,z]$ and equality after coefficient localization:
$$
 g(s)b(s,z)\lambda^s\eta^z1
                  =Q_0(s,z)\lambda^s\eta^{z+1}1.       \tag{3}
$$
The displayed coefficient arguments of $g$ are suppressed from now on.

Here is the uniform torsion argument required before specializing (3).
Expand the finite operator $Q_0$ by Leibniz. For some uniform integers
$a,c\geq0$, the difference of the two sides is
$\lambda^{s-a}\eta^{z-c}u(s,z)$, with $u(s,z)\in M[s,z]$.
Its equality to zero in the localized ambient means that for one
$h(s)\in R[s]$, $h(s,0,0)\ne 0$, and some $A,C\geq0$,
$h(s)\lambda^A\eta^Cu(s,z)=0$ in $M[s,z]$.
There are only finitely many coefficients to clear; time torsion-freeness
removes $\eta^C$. For every integer $k\geq a+A$ and every integer $j$,
$h(k)\lambda^{k-a}u(k,j)=0$ in $M$. Thus, in $M[\eta^{-1}]$,
$$
 h(k)g(k)b(k,j)\eta^jv_k
                  =h(k)Q_0(k,j)\eta^{j+1}v_k.          \tag{4}
$$
There is no fixed shift of $v_k$: the already present large positive
power $\lambda^{k-a}$ kills the finite localization torsion.

Exclude the finitely many integer roots of the nonzero polynomials
$h(s,0,0)g(s,0,0)$. Set $j=-\tau_k-1$ in (4).
Then $\eta^{j+1}v_k=\eta^{-\tau_k}v_k\in E$;
every operator on the right preserves $E$. The left has a potential
time pole of order one. If $A_k(x)\ne 0$ is the leading coefficient of
$v_k$, that pole has coefficient
$$
 h(k,D_x,0)g(k,D_x,0)b(k,-\tau_k-1)A_k(x).
$$
The first two factors are invertible on polynomials, so the scalar
$b(k,-\tau_k-1)$ must be zero. This proves(1) in the actual module.
∎

For the actual first-survival indices, put $p(s,t)=b(s,-t-1)$.
The result is a **single nonzero polynomial** $p\in K[s,t]$ with
$p(k,\tau_k)=0$ at **every sufficiently large integer** $k$ for which
$\tau_k$ is finite. No purity or eventuality was used to produce this equation.

### Integer root constraint forces a gap

Let $\tau_k$ be finite nondecreasing integers with $\tau_k>k$ and
$\tau_k-k\to\infty$. If a nonzero $p(s,t)\in K[s,t]$ vanishes at
$(k,\tau_k)$ for all large $k$, then $\liminf\tau_k/k>1$.

**Proof.**
Resolve the finitely many coefficients of $p$ in a basis of their
$\mathbb Q$-span. Evaluation at integer pairs shows that every rational
component vanishes; select a nonzero component and assume $p\in\mathbb Q[s,t]$.
Let $p_d$ be its highest total homogeneous part and $a(x)=p_d(1,x)$,
a nonzero univariate polynomial. On every compact interval,
$k^{-d}p(k,kx)\to a(x)$ uniformly.

Suppose $\liminf\tau_k/k=1$. Then $a(1)=0$. Choose $\delta>0$ so
$a$ has no roots in $(1,1+3\delta]$. For large $k$, no root ratio lies
in $[1+\delta,1+2\delta]$. Thus every ratio is low, below $1+\delta$,
or high, above $1+2\delta$. A high-to-low transition at consecutive
large indices contradicts monotonicity, since
$(1+2\delta)k>(1+\delta)(k+1)$.
Arbitrarily late low indices therefore imply that all late indices
are low. For each $0<\epsilon<\delta$, uniform convergence excludes
ratios in $[1+\epsilon,1+\delta]$ eventually. Hence
$g(k)=\tau_k-k=o(k)$ for all $k$, not just along a subsequence.

The polynomial $q(s,z)=p(s,s+z)$ is nonzero. There are only finitely
many integers $z$ for which $q(s,z)$ is identically zero in $s$:
such $z$ are common roots of its finitely many coefficient polynomials.
The condition $g(k)\to\infty$ avoids these exceptions eventually.
If $d_s=\deg_s q$, every remaining integer $z$ can occur as $g(k)$
for at most $d_s$ indices. But for any fixed $\epsilon>0$ and all
$K\leq k\leq N$ with $K$ sufficiently large, $0<g(k)\leq\epsilon N$.
There are at most $\epsilon N+1$ possible values. Thus
$N-K+1\leq d_s(\epsilon N+1)$, impossible for $\epsilon d_s<1$.
If $d_s=0$, there are no nonexceptional roots and the contradiction is immediate.
∎

### Full field implication - GVCFORMALALLLINEAR20261116

The target implication holds for every formal symbol and every polynomial
phase in every dimension and degree over every characteristic-zero field.

**Proof.**
The actual orders $\tau_k$ are nondecreasing, since $T^{k+1}P^j\ne 0$
implies $T^kP^j\ne 0$. Purity gives $\tau_k>k$.
If some $\tau_N=\infty$, then every height
$h_m=\max\{k:T^kP^m\ne 0\}$ is less than $N$, and purity supplies a
positive global deficit fraction. Otherwise all $\tau_k$ are finite.

The retained formal commutator lemma identifies individual eventuality
with eventual vanishing of $T^kP^{k+d}$ for each fixed $d\geq0$.
Together with purity, applying this to the finitely many shifts
$0\leq d\leq D_0$ proves $\tau_k-k\to\infty$.
The preceding two lemmas, with $p(s,t)=b(s,-t-1)$, give
$\liminf\tau_k/k>1$. Choose $a>1$ and $K$ so $\tau_k\geq ak$ for $k\geq K$.
If $h_m\geq K$, then $\tau_{h_m}\leq m$, hence $h_m\leq m/a$.
If $h_m<K$ and $m\geq2K$, then $h_m/m<1/2$.
The finitely many earlier $m$ satisfy $h_m\leq m-1$ by purity.
Therefore $\epsilon=\inf_{m\geq1}(m-h_m)/m>0$.
For a marker of degree $q$, the formal commutator lemma gives
$T^m(QP^m)=0$ as soon as $m-q>h_m$.
Any integer $C\geq1/\epsilon$ consequently works for $m>Cq$.
∎

## Boundary and controls

This is a field, fixed-input existence statement. It does not prove
pure-to-safe GVC, compute $C$, or give nonreduced-family uniformity.
It does not assert a raw or reversed scalar-time ODE, a finitely generated
order/time product Rees module, or the proposed quasi-superadditivity of
$\tau$. All polynomial-marker assumptions stay in the original amplitude
space; completed coefficient extension is used only with faithful descent.
The old safe Hermite Rees obstruction and lacunary no-time-ODE control remain valid.
For $P=x^D$ and $\lambda=\xi^r$ times an arbitrary formal unit,
$\tau_k=\lceil rk/D\rceil$ and
$\prod_{j=0}^{D-1}(D\tau_k-rk-j)=0$.
The unsafe seed has $\tau_k=k+1$, which satisfies a polynomial relation
but fails $\tau_k-k\to\infty$. Abstract $\tau_k=k+\lceil\sqrt k\rceil$
is monotone and has an unbounded sublinear excess, but has no such
all-index polynomial relation; it is not an actual counterexample.
The root constraint alone is insufficient without monotonicity: choose
$\tau_k=k+\sqrt k$ at square $k$, and $\tau_k=2k$ otherwise.
Then $((\tau_k-k)^2-k)(\tau_k-2k)=0$ and the excess diverges,
but high-to-low transitions at squares violate monotonicity.


## Evidence and retained boundaries

- [Exact checker](verify_generic_bernstein_gap.py) and [receipt](generic_bernstein_gap_verification.json).
- [Primary-source and proof audit](generic_bernstein_gap_source_audit.json).
- [Pre-update status and protected records](generic_bernstein_gap_preupdate_20261116.json).
- [Primary 1991 source](sources/mebkhout_narvaez_1991_generic_twist_20261116.pdf), physical pages 5–12 and 16 (printed pages 230–237 and 241), read in text and rendered form.
- [Formal commutator and exact gap equivalence](FORMAL_SHIFTED_POWERS_AND_APOLAR_NILPOTENCE_GAPS.md).
- [Raw-source rigidity](FORMAL_SYMBOL_SOURCE_RIGIDITY.md), [reversed-time obstruction](UNIVERSAL_REVERSED_TIME_EQUATION_COUNTERCONTROL.md), and [actual product obstruction](IDEAL_TOPOLOGY_BRIDGE_AND_ACTUAL_APOLAR_CONTRACTIONS.md) remain valid.

The earlier formal-lattice, order/time Rees and separable-shift drafts in `research/.scratch/`
remain sealed exploratory evidence. This proof supplies an estimate on actual powers
through a new finite polynomial constraint. It does not revive their failed estimates.
Both complete exceptional Haar MZ questions and the broader all-degree programme remain open.
