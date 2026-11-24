# Uniform upper-response bounds from the actual scalar leaves

## Scope

Let `q=p^s`, where `p` is prime and `s>=1`. Let `n>=1`, let `gamma>0`
be rational with `v_p(gamma)>=0`, and let `1<=h<=q`. All conclusions
concern the actual pair

\[
 T_{n,\gamma,q}(t)=\sum_{k=0}^{n}\binom nk(\gamma)_{qk}t^k,
 \qquad T_{n+1,\gamma+h,q}(t).
\]

The imported complete transfer is
[SFC_PRIME_POWER_SCALAR_TRANSFER.md](SFC_PRIME_POWER_SCALAR_TRANSFER.md),
especially Sections 5–8. Its actual seeds, closed-disk convergence,
operator convention, and upper functional are retained. This note
proves two extensions:

* a uniform stronger coprimality criterion at every prime-power gap;
* `deg gcd <= h-1` at every **prime** gap and every allowed offset.

It also proves a denominator-depth bound at all prime powers. It does
not assert the `h-1` bound at every proper prime power, or unrestricted
Strong Factorial.

Write

\[
 b=\frac{q-1}{p-1},\quad B=b-s,\quad
 C=\gamma+qn+q,\quad\kappa=\frac{q(n+1)}{(C)_h},
 \quad k=v_p(\kappa),\quad\delta=(-1)^q.
\]

Every lower root yields `u=(-1)^(q+1)/t` with `v_p(u)>=b`. Put
`v=u/p^b`. The retained transfer gives a single monic polynomial family

\[
 Q(X;v)=Q_0(X)+\sum_{j\ne j_*}b_j(v)G_j(X),
 \quad b_j\in qv\mathbb Z_p\langle v\rangle,
 \tag{1}
\]

where `c=gamma+q-1`, `j_*=c mod q`, `a_j=c-j`,
`Q_0=prod_(j!=j_*)(X-a_j)` and `G_j=Q_0/(X-a_j)`. For a polynomial
`P` of degree at most `q-1`, define the actual integral upper functional
by the triangular expansion

\[
 P(X)-P(-qn)=\sum_{j=0}^{q-2}d_j\ell_j(X),\qquad
 \ell_j(X)=(X+qn)\prod_{r=0}^{j-1}(c-X-r),
\]
\[
 \mathcal T_h(P)=\sum_{j=\max(0,q-h)}^{q-2}
               d_j(C+q-j)_{h-q+j}.
\]

An empty sum is zero. Set

\[
 T_0=\mathcal T_h(Q_0),\quad
 D=\min_j v_p(\mathcal T_h(G_j)),\quad
 A_0=Q_0(-qn),\qquad v_p(A_0)=B.
\]

For finite `D`, the retained estimates are

\[
 Q(-qn;v)=A_0+p^{B+1}v A(v),\qquad
 \mathcal T_h(Q(X;v))=T_0+p^{s+D}v R(v),
 \quad A,R\in\mathbb Z_p\langle v\rangle.
 \tag{2}
\]

Every actual common root must satisfy

\[
 F(v):=Q(-qn;v)-\delta\kappa p^b v\mathcal T_h(Q(X;v))=0.
 \tag{3}
\]

This is a necessary condition retaining actual seeds; the omitted
resonant leaf is not asserted to follow from it.

## 1. A stronger uniform coprimality criterion

**Theorem 1.** If `D` is finite and

\[
                           k+2s+D>0,
 \tag{4}
\]

the mixed pair is coprime over `Q[t]` and over every characteristic-zero
field. In particular the parameter-only sufficient condition
`v_p(kappa)>-2s` holds at every `1<=h<=q`. If all `T_h(G_j)=0`, the
pair is also coprime, for arbitrary `kappa`.

**Proof.** The constant coefficient of `F` has valuation exactly `B`.
By (2), every coefficient of degree at least two has valuation at least

\[
              \min\{B+1,\ k+b+s+D\}>B.
\]

Only the linear coefficient can have valuation at most `B`. If it
vanishes, or has valuation greater than `B`, the normalized restricted
series has unit constant and all other coefficients in `pZ_p`, so it
has no integral root. Otherwise divide by an actual coefficient of
minimal valuation, the minimum of the constant and linear valuations.
The normalized series is integral and restricted, its linear coefficient
is a unit, and every coefficient of degree at least two is in `pZ_p`.
Its divided difference at any two integral arguments is a unit, so it
has at most one integral root. This includes cancellation in the linear
coefficient: the comparison uses its actual resulting value.

All common roots of the rational pair are represented in this disk.
Consequently the pair has at most one distinct common root, which must
be rational by separability and conjugacy over `Q`. No squarefreeness
assumption is needed for this uniqueness-to-rationality step.

A rational common root is impossible by the retained rational Gamma
Gram factorization: one of the two powers is even, and the corresponding
positive rational shape gives a strictly positive rationally weighted sum of
squares. Explicitly, `M_beta(X^j)=(beta)_j` has

\[
 (\beta)_{i+j}=\sum_r L_{ir}\,r!(\beta)_r\,L_{jr},\qquad
 L_{ir}=\binom ir(\beta+r)_{i-r},
\]

with unit triangular `L` and positive rational diagonal weights.
This applies to `beta=gamma` or `gamma+h`. Rational Bezout descent
then gives the stated field scope.

If `D=+infinity`, the upper functional is the constant `T_0` on the
entire family. All coefficients of `F` of degree at least two again
have valuation greater than `B`, for arbitrary `kappa`. The identical
constant-or-linear comparison and rational exclusion apply. ∎

For ordinary binomials with `gamma=an+1`, `h=a<=q`, `N=n+1`, this gives
the single support-independent condition

\[
 v_p(q+a)+\sum_{j=1}^{a-1}v_p((q+a)N-j)<3s+D.
 \tag{5}
\]

Dropping `D>=0` gives the simpler threshold `<3s`. Thus the earlier
`<2s` condition improves uniformly, without an index or support scan.
For the entire power family `a=p^w<q`, it suffices that
`(a-1)/(p-1)<3s`.

## 2. A stronger coefficient filtration

The final actual leaf seed in the retained construction has the form

\[
 H_j(y)=\sum_{m\ge0}h_{jm}y^m,\qquad
                  v_p(h_{jm})\ge m,
 \tag{6}
\]

and its constant is a unit. This stronger statement survives absorption
of the entire finite rising product.

Indeed the seed immediately after the last residue extraction has
coefficients `g_(pm+r)(p/b_unit)^m`, giving (6). The final multiplication
part is

\[
 F_*(y)=y+\sum_{r\ge1}a_r^*y^{p^r},\qquad
                        v_p(a_r^*)=p^r-1.
\]

Each absorbed factor is
`e/p^rho+p^(s-rho)(L+ell)`, where `rho=v_p(e)<=s-1`.
Its diagonal Euler part preserves the bound at degree `m`.
Its linear multiplication term raises degree by one and contributes
valuation at least one. Its term of degree `p^r` contributes valuation
at least `1+(p^r-1)=p^r`. Thus each factor preserves (6), including
constant-unit reduction. The inverse seed has the same coefficient
bound, by formal inversion after the rescaling `w=py`.

Let `T=L-L(0)=-yD_y+F_*(y)`. If a series has coefficient valuations
at least `m-r`, applying `T` gives valuations at least `m-r-1`.
All operators are also integral. It follows that the exact leaf terms

\[
 \tau_{jk}(v)=\left.\frac{T^kH_j}{H_j}\right|_{y=\nu v}
\]

satisfy the coefficientwise bound

\[
              v_p([v^m]\tau_{jk})\ge\max\{0,m-k\}.
 \tag{7}
\]

Using the retained divided-Taylor denominator bounds, the coefficient
of `v^m`, `m>=1`, in the normalized interpolation matrix has valuation
at least `m`; the forcing coefficient has valuation at least `s+m-1`.
The inverse matrix has the same slope-one filtration, since after
`w=pv` it is an integral formal matrix with identity constant term.
Consequently

\[
 v_p([v^m]b_j(v))\ge s+m-1\quad(m\ge1).
 \tag{8}
\]

This is a coefficient estimate on the actual restricted family. It does
not evaluate an unbounded formal series on a unit disk.

Since `v_p(G_j(-qn))>=B-s+1`, (8) gives

\[
 v_p([v^m]Q(-qn;v))\ge B+m\quad(m\ge1),
\]
\[
 v_p([v^m]\mathcal T_h(Q(X;v)))\ge s+m-1+D
                                                 \quad(m\ge1).
 \tag{9}
\]

**Theorem 2 (denominator-depth bound).** For finite `D`, let
`e=k+2s+D` and `R=max(1,2-e)`. The actual mixed pair has at most `R`
distinct common roots. Its gcd has even degree at most
`min(n,2 floor(R/2))`.

**Proof.** The constant coefficient of `F` has valuation `B`. For
`r>=2`, (9) gives

\[
 v_p([v^r]F)\ge\min\{B+r,\ B+e+r-2\}.
\]

Both bounds exceed `B` when `r>R`. Thus the minimal coefficient
valuation of the restricted series is attained at an index at most
`R`. After normalization, its largest unit coefficient has index at
most `R`; the restricted-series root-count lemma bounds its integral
roots by that index. All actual common roots lie in that disk.
The retained Theorem `thm:prime-power-lower-one` in
[ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex](ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex) proves that the lower
polynomial is squarefree for the present positive rational shapes.
This last input converts a distinct-root count to the asserted ordinary
gcd-degree bound `min(n,R)`. The parity refinement follows from the
real-root exclusion immediately below. ∎

**Real-root exclusion and parity.** Neither adjacent pair has a common
real root. If `n` is even, the lower polynomial is the Gamma moment of
the nonzero real polynomial `(1+tX^q)^n`, an even power; if `n` is odd,
the upper polynomial has even power `n+1`. Apply the displayed Gram
factorization to the square of `(1+tX^q)^(n/2)` or
`(1+tX^q)^((n+1)/2)`, respectively. The positive real diagonal weights
and invertible triangular matrix make that moment strictly positive.
Thus a real common root is impossible. The rational gcd has real
coefficients, and its nonreal roots occur in conjugate pairs with the
same multiplicities. Its degree is consequently even. This argument
uses a specific positive definite moment functional; it does not assume
positivity for arbitrary binomial coefficients or coefficient fields.

The bound depends on a valuation depth, not on a finite list of
supports or moment indices. Section 3 provides an independent offset
bound at every prime gap.

## 3. Every offset at a prime gap

**Theorem 3.** Let `p` be any prime, let `n>=1`, and let `gamma>0` be
rational with `v_p(gamma)>=0`. For every `1<=h<=p`,

\[
 \deg\gcd(T_{n,\gamma,p},T_{n+1,\gamma+h,p})
               \le 2\left\lfloor\frac{h-1}{2}\right\rfloor.
 \tag{10}
\]

For `h=1,2` the pair is coprime. For larger `h`, this is a bounded
degree reduction, not a coprimality assertion.

**Proof.** The cases `h=1,2` are in the retained transfer. Suppose
`3<=h<=p`; thus `p>=3`. Here `s=1`, `b=1`, and `B=0`.
Put `w=p nu v`. By (8),

\[
          \widehat Q(X;w)=Q(X;w/(p\nu))\in\mathbb Z_p[[w]][X].
\]

Work formally in the inner variable `W=py`, keeping the outer `w`
and all coefficients of `Q` constant under Euler operations. The final
inner operator satisfies

\[
 p(L-L(0))=-pW D_W+pF_*(W/p)
                  \equiv G(W)\pmod p,
\]
\[
             G(W)=\sum_{r\ge0}W^{p^r},\qquad G-G^p=W.
 \tag{11}
\]

For the last identity of operators, the final coefficient formula gives
`a_r^*/p^(p^r-1)=1 modp`; it holds for odd `p` and for `p=2`.
The derivative part is divisible by `p`. The scaled seed
`H_j(W/p)` is an integral formal series with unit constant by (6), so
it can be canceled formally after reduction. No evaluation of `G` at a
unit argument is made: this formal calculation controls coefficients,
and the actual parameter `w=p nu v` has positive valuation.

The exact scalar leaf equations reduce to

\[
                \overline{\widehat Q}(a_j+G(w);w)=0
                         \qquad(j\ne j_*).
\]

There are `p-1` distinct node residues and their pairwise differences
are units. The polynomial is monic of degree `p-1`, so these equations
give the exact identity in `F_p[[w]][X]`

\[
              \overline{\widehat Q}(X;w)=\bar Q_0(X-G(w)).
 \tag{12}
\]

Put `d=h-2`. The functional `T_h` annihilates polynomials of degree
less than `p-h+1`. On `X^(p-h+1)` its value is `(-1)^(p-h)`, because
the first contributing basis polynomial is `ell_(p-h)` and its
rising-basis upper coefficient is one. Therefore

\[
 Z\longmapsto\mathcal T_h(Q_0(X-Z))
\]

is a polynomial of exact degree `d`; its leading coefficient is

\[
                    (-1)^p\binom{p-1}{d},
 \tag{13}
\]

a `p`-adic unit. Lower terms of `Q_0` cannot contribute to that
coefficient after applying the degree cutoff of the functional.

Since `d<=p-2` and `G(w)=w+O(w^p)`, the coefficient of `w^d` in
the reduced upper response is the unit (13). Every coefficient of
`T_h(Qhat)` is integral. Returning to `v`, we obtain

\[
 v_p([v^m]\mathcal T_h(Q(X;v)))\ge m\quad(m\ge0),
 \qquad v_p([v^d]\mathcal T_h(Q(X;v)))=d.
 \tag{14}
\]

Also `v_p([v^m]Q(-pn;v))>=m` for `m>=1` by (9).
If `k>=0`, Theorem 1 already gives coprimality. Suppose `k<0`.
In the mixed auxiliary (3), its coefficient at index `h-1=d+1`
has valuation **exactly** `k+h-1`: the upper term has that valuation
by (14), whereas the lower evaluation term has valuation at least
`h-1`, strictly larger. Every coefficient at index `r>h-1` has
valuation strictly greater than `k+h-1`, since its two contributions
have valuations at least `r` and `k+r`.

Thus the minimal coefficient valuation is attained among the first
`h` coefficients. Normalize by one actual coefficient of least valuation
in that finite set. The resulting integral restricted series has largest
unit index at most `h-1`, and hence at most `h-1` integral roots, with
multiplicities. All actual common roots are among these roots.
The lower moment polynomial is squarefree by the retained theorem,
which gives gcd degree at most `h-1`. The even-degree observation in
Section 2 proves (10). ∎

The root-count lemma used above remains valid for ramified roots.
If an integral restricted series has largest unit index `r`, division
by any integral root `alpha` gives a restricted quotient with largest
unit index `r-1`, using coefficients
`b_i=sum_(k>=i+1)a_k alpha^(k-i-1)`. Iteration gives the bound. This
requires neither unramifiedness nor positivity of the coefficients.

## 4. Why the same translation proof does not automatically cover p^s

For `s>1`, after the rescaling `w=p nu v`, (8) gives

\[
 \widehat Q(X;w)-Q_0(X)\in p^{s-1}w\mathbb Z_p[[w]][X].
\]

Thus its reduction is `Q_0`, not the translation in (12). The naive
extension `Qhat=Q0(X-G(w)) modp` is false: the coefficient of `w` on
its right side is `-Q0'`, whose leading coefficient is `-(q-1)`, a
unit, while the actual coefficient vanishes modulo `p`.

The nodes form clusters modulo `p^(s-1)`. A correct higher-power
extension must retain the cluster equations and their factorial
coefficient carries. The slope-one bound and Theorems 1–2 remain
valid, but (13) cannot simply be imported with `p` replaced by `q`.
No all-prime-power `h-1` bound is claimed here.

The actual guard `p=q=5,h=5,n=1,gamma=11` has
`T0=-32400`, `D=0`, `T5(Q0')=4820`, and `v5(kappa)=-1`.
Its constant and first normalized responses are degenerate at the
precision used in the earlier criterion. The present proof retains
them and finds the guaranteed unit at translation degree `h-2=3`;
Theorem 1 also excludes a common root because `k+2s+D=1>0`.
This guard prevents replacement of the actual upper polynomial by an
assumed nonzero linear response.

## 5. Exact controls

The companion [verify_sfc_prime_power_upper_response.py](verify_sfc_prime_power_upper_response.py) constructs actual complete
leaf operators and seeds, checks the stronger slope filtration, and
reconstructs the response coefficients through explicit finite orders.
It tests positive rational shapes, the `p5,gamma11` guard, a deeper
prime denominator, and a proper-prime-power countercontrol to the
unscaled translation shortcut. These checks verify exact algebra in
their recorded finite domains. The all-parameter claims follow from
the proofs above and the retained transfer, not from a parameter scan.


Replay from the repository root:

```sh
python research/factorial/image-rigidity/verify_sfc_prime_power_upper_response.py
```

The exact output is recorded in
[sfc_prime_power_upper_response_verification.json](sfc_prime_power_upper_response_verification.json).
The verifier pins this note, the complete scalar-transfer proof and checker,
and that transfer's three retained proof dependencies. The mathematical
inputs are the actual mixed upper functional, closed-disk leaf transfer,
and prime-power squarefreeness; the latter is Theorem
`thm:prime-power-lower-one` of the retained arithmetic-frontier proof.

Internal mathematical review: the root agent independently reconstructed
all three theorems; the factorial-review agent separately checked seed
slope preservation, loss under the inner operator, interpolation inversion,
the denominator-depth root count, formal prime translation with the outer
parameter fixed, and the top functional coefficient. These are internal
mathematical reviews. The exact replay is a separate finite algebra check;
neither claim implies external peer review or formal verification.
