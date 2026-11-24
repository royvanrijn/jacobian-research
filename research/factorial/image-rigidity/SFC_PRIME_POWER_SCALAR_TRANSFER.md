# A uniform prime-power collision transfer

## Scope and retained inputs

Fix a prime `p`, an integer `s>=1`, and put

\[
 q=p^s,\qquad b_s=\frac{q-1}{p-1},\qquad B_s=b_s-s.
\]

For `n>=1` and a positive rational `gamma` with `v_p(gamma)>=0`, write

\[
 T_{n,\gamma,q}(t)=\sum_{k=0}^n\binom nk(\gamma)_{qk}t^k.
\]

This note proves a uniform normalization of every admissible collision
polynomial, an all-index mixed-offset criterion, and a closed-disk
scalar-leaf transfer. Its statements apply simultaneously to all the
displayed parameters; no support list is searched.

The indispensable retained inputs are in
`research/factorial/image-rigidity/`:

* [ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex](ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex), Lemma
  `lem:artin-hasse-closure`, and Lemma `lem:frobenius-noncollision`.
  The latter forbids a primitive zero cyclic column when its polynomial
  has reduced degree strictly below `p^s-1`.
* [LOWER_TWO_FOURTH_POWER_GAP_SFC.tex](LOWER_TWO_FOURTH_POWER_GAP_SFC.tex), Section 1 and the proof of
  Theorem `thm:boundary`: explicit nonresonant unit bases and their
  divided-Taylor consequences. The extra condition on `v_p(n)` in that
  note's application is not assumed here.
* [DIAGONAL_PRIME_GAP_SFC.tex](DIAGONAL_PRIME_GAP_SFC.tex), Theorem `thm:offset-rank`: the actual
  mixed-shape matrix for every `1<=h<=q`, including multiple lower roots.

These sources were inspected in the worktree based at
`19eb9912fb60142f8b620dac24896931d4ff3e46`. The argument does not import
an unproved arbitrary-seed rank assertion.

## 1. The exact stage class

Work in a complete algebraically closed extension `K` of `Q_p`, with
valuation `v_p(p)=1`, ring `O`, maximal ideal `m`, and residue field `k`.
The value group may be nondiscrete. At stage `r>=1` use exactly

\[
 A=c-\theta_z+F(z),\quad c\in\mathbb Z_p+p^r\mathcal O,
 \quad F(z)=z+\sum_{i\ge1}a_i z^{p^i},\quad a_i\in p\mathbb Z_p,
\]
\[
 h\in\mathcal O[[z]],\qquad \bar h=\overline{h(0)}\in k^\times,
 \quad\Pi_{n,p^r}(X)=\prod_{\ell=0}^{n-1}(X+p^r\ell).
 \tag{1}
\]

The cyclic column `C_(p^r)(f;u)` consists of the coefficients of the
remainder at `z^(p^r)=u`, evaluated after all inner Euler operations.
At `v_p(u)>0`, its vanishing means that `f` vanishes at every root of
`z^(p^r)=u`. This includes the convergent series in (1).

Assume `v_p(u)>=b_r=(p^r-1)/(p-1)` for the following stage assertion.
In the exact Artin--Hasse split from `p^r` to `p^(r-1)`, the nonresonant
blocks have operators `e+pL`, `e=1,...,p-1`. After absorbing the unit
rising product, their actual seeds still reduce to constants units.
The first `p^(r-1)` powers of `L` applied to such a seed form an
integral unit basis after cyclic evaluation. For a primitive integral
polynomial `Q`, a zero original cyclic column therefore implies

\[
 Q^{[j]}(e)\in p^{p^{r-1}-j}\mathcal O
 \quad(0\le j<p^{r-1},\ 1\le e<p).
 \tag{2}
\]

The superscript is the divided Taylor coefficient. These are actual
ideal containments, not merely statements of positive valuation.

## 2. Forced maximal normalization at every resonant stage

**Theorem 1.** Suppose (1) holds at stage `s`, `n>=1`, `u!=0`, and
`v_p(u)>=b_s`. Let `Q in O[X]` be primitive of ordinary degree at most
`q-1`. If

\[
        \mathcal C_q(\Pi_{n,q}(A)Q(A)h;u)=0,
 \tag{3}
\]

then

\[
        \widetilde Q(X)=\frac{Q(p^{s-1}X)}{p^{B_s}}
             \in\mathcal O[X]
\]

is primitive, its constant coefficient is a unit, the coefficient of `X^(q-1)`
in the original `Q` is a unit, and

\[
        v_p(Q(p^s a))=B_s\quad\text{for every }a\in\mathcal O.
 \tag{4}
\]

In particular `v_p(Q(0))=B_s`. This theorem allows ramified coefficients.

**Proof.** At `s=1`, reduction of (3) is
`Abar^n Qbar(Abar) hbar=0` on `k[z]/(z^p)`. The operator is cyclic with
the distinct eigenvalues `F_p`, and the seed has a nonzero component in
every eigenspace. Hence `Qbar` vanishes at all `p-1` nonzero eigenvalues.
Since its degree is at most `p-1`, primitivity forces

\[
                 \bar Q=a(X^{p-1}-1),\qquad a\ne0.
\]

Its constant coefficient is a unit. Evaluation at `pa`, `a in O`,
therefore remains a unit. This proves the base case `B_1=0`.

Suppose `s>=2`, and put `q'=p^(s-1)`. Equation (2) implies

\[
 \bar Q=(X^{p-1}-1)^{q'}R(X),\qquad \deg R\le q'-1.
\]

The first unit coefficient of `Q` consequently occurs at an index
`r_0<=q'-1`. Writing `Q=sum c_jX^j`, choose an actual nonzero
coefficient `p^j c_j` of minimal valuation

\[
             \mu=\min_j\{v_p(c_j)+j\}\le r_0\le q'-1.
\]

The resonant block carries `Q(pX)` and the exact factor
`Pi_(n,q)(pL)=p^n Pi_(n,q')(L)`. Dividing by that actual minimal
coefficient gives an integral primitive polynomial with reduced degree
at most `floor(mu)`. If `mu<q'-1`, the retained reduced-degree
noncollision lemma at stage `s-1` contradicts this zero cyclic column.
Thus `mu=q'-1`. Its first unit coefficient is at index `q'-1`.
The displayed reduced factorization then forces
`R=a X^(q'-1)`, `a!=0`: its lower coefficients vanish because the
prefactor has unit constant. Hence the coefficient of `X^(q-1)`
in `Qbar` is a unit. This also holds in the prime-stage base case.
The integer normalization

\[
                    Q_1(X)=Q(pX)/p^{q'-1}
\]

is integral and primitive. All its coefficients of index at least
`q'` have positive valuation. Its actual ordinary degree may exceed
`q'-1`, but its reduced degree is at most `q'-1`.

For the subsequent induction use the same proof with a polynomial of
unrestricted ordinary degree and reduced degree at most the stage
boundary. The nonresonant divided-Taylor expansion is finite and all
its higher terms are integral. The first-unit argument and resonant
normalization are unchanged. At each stage the minimum must equal
its maximal permitted value `p^(r-1)-1`; otherwise the retained
noncollision lemma applies, which itself permits unrestricted ordinary
degree. At the final prime stage the reduced polynomial is again a
nonzero multiple of `X^(p-1)-1`.

The total normalization exponent is

\[
 \sum_{r=1}^{s-1}(p^r-1)=\frac{p^s-1}{p-1}-s=B_s.
\]

Thus the displayed `Qtilde` is integral with unit constant. Evaluating
it at `pa` proves (4). Every intermediate division was by an actual
coefficient before its valuation was proved to equal an integer; no
fractional power of `p` was introduced. ∎

**Degree qualification.** The same proof gives Theorem 1 whenever the
initial reduced degree is at most `q-1`; ordinary degree may be larger,
provided the initial polynomial is integral. At each later stage the
high-degree coefficients gain the strictly positive valuation required
by the reduced-degree argument. The actual mixed relation below has
ordinary degree at most `q-1`.

## 3. A uniform all-index mixed-offset criterion

Put

\[
 C=\gamma+qn+q,\quad \epsilon=(-1)^{q+1},\quad
 \kappa_h=\frac{q(n+1)}{(C)_h},\qquad1\le h\le q.
\]

**Theorem 2.** For every positive rational `gamma` integral at `p`,
and every `n>=1`,
if

\[
                          v_p(\kappa_h)>-s,
 \tag{5}
\]

then `T_(n,gamma,q)` and `T_(n+1,gamma+h,q)` are coprime over `Q[t]`
and over every characteristic-zero field.

**Proof.** A common root gives `u=epsilon/t`. Among `qk` consecutive translates of `gamma in Z_p`, exactly
`qk/p^r` are divisible by `p^r` for `1<=r<=s`. The lower monic reciprocal therefore implies
`v_p(u)>=b_s`.

Use the actual fixed-shape operator
`A=gamma+q-1-zD_z+z`, seed 1, and polynomials

\[
 \ell_j(X)=(X+qn)\prod_{r=0}^{j-1}(\gamma+q-1-X-r),
                         \quad0\le j\le q-2.
\]

The mixed theorem supplies a relation, normalized with integral
coefficients `a_j`, one a unit, of the cyclic columns
`Pi_(n,q)(A)[ell_j(A)-epsilon*kappa_h*u*v_(q-1-j)]1`, where

\[
 v_k=(C+k+1)_{h-1-k}\ (k<h),\qquad v_k=0\ (k\ge h).
\]

Let `T=sum_j a_j v_(q-1-j)` and

\[
 Q(X)=\sum_j a_j\ell_j(X)-\epsilon\kappa_huT.
 \tag{6}
\]

Its nonconstant coefficients are integral and one is a unit, because
the `ell_j` have successive degrees and leading coefficients `(-1)^j`.
The correction has valuation strictly above `b_s-s=B_s>=0`, so `Q`
is integral and primitive. Theorem 1 gives `v_p(Q(-qn))=B_s`.
But (6) gives `Q(-qn)=-epsilon*kappa_h*u*T`, of valuation greater than
`B_s`. This contradiction proves the theorem. ∎

Thus every actual common root, without condition (5), necessarily
satisfies

\[
 v_p(\kappa_h)\le-s,\qquad
 v_p(u)+v_p(T)=B_s-v_p(\kappa_h).
 \tag{7}
\]

For this necessity, (6) is still integral: if its constant term had
negative valuation, dividing by that term would create a primitive
polynomial with nonzero constant reduction and all higher coefficients
in the maximal ideal, contradicting reduced-degree noncollision.

## 4. An ordinary-binomial condition without support enumeration

For ordinary lower exponent `a`, take `gamma=an+1` and `h=a<=q`.
With `N=n+1`, cancellation gives

\[
 \kappa_a=\frac{q}{(q+a)\prod_{j=1}^{a-1}((q+a)N-j)}.
 \tag{8}
\]

Therefore all such supports satisfy the following single sufficient
condition at every terminal index:

\[
 v_p(q+a)+\sum_{j=1}^{a-1}v_p((q+a)N-j)<2s.
 \tag{9}
\]

No coefficient choices enter this test, and its assertion is over every
characteristic-zero field. Empty products at `a=1` have their usual value.

**Uniform power corollary.** If `a=p^w<q=p^s`, `w>=0`, then

\[
 v_p(\kappa_a)=s-\frac{a-1}{p-1}
\]

at every terminal index. Hence the entire support family is closed
whenever `(a-1)/(p-1)<2s`. Indeed `v_p(q+a)=w` and, for `1<=j<a`,
`v_p((q+a)N-j)=v_p(j)`, since the first term is divisible by `p^w`
and the second has smaller valuation. Summing and using Legendre's
formula gives the displayed constant. Some of this parameter range is
already covered by other retained theorems; no independent novelty
claim is attached to overlapping cases.

More generally, write `a=p^w A<q` with `p` not dividing `A`. Then (8)
is equivalent to the valuation identity

\[
 v_p(\kappa_a)=s-v_p(a!)-
 v_p\binom{(A+p^{s-w})N-1}{A-1}.
 \tag{10}
\]

To verify it, factor `(a-1)!` out of the consecutive product in (8).
Kummer's digit-sum formula, or the factorial valuation formula directly,
gives
`v_p binom((q+a)N-1,a-1)=v_p binom((A+p^(s-w))N-1,A-1)`:
both top and bottom end in the same `w` digits `p-1`, which contribute
no carries. The remaining constant is `v_p(a!)`.

Thus the criterion fails only when that one binomial coefficient has
valuation at least `K=2s-v_p(a!)`. When `A<p` and `K>=1`, the exceptional
terminal set is exactly the union of at most `A-1` residue classes

\[
              (A+p^{s-w})N\equiv j\pmod {p^K},
                         \quad1\le j<A.
 \tag{11}
\]

The classes are distinct. This is a uniform arithmetic description,
not a finite bound on `N`: any nonempty class contains arbitrarily
large terminal indices. If `K<=0`, this particular inequality provides
no exclusion, and no closure is inferred from it.

For arbitrary `A>1` and `K>0`, put `m=A-1` and `L=floor(log_p(m))`.
The exceptional condition in (10) depends only on `N mod p^(K+L)`.
Indeed, Vandermonde's identity and
`v_p binom(p^(K+L),j)>=K+L-v_p(j)>=K` for `1<=j<=m` show that
shifting the top argument by `p^(K+L)` changes the binomial coefficient
by a multiple of `p^K`. The multiplier `A+p^(s-w)` is a unit, so the
same period applies to `N`. This is a precision bound, not a finite
terminal cutoff.

If `K>L`, there is a sharper exact description. Set

\[
 \delta_j=v_p\!\left(m\binom{m-1}{j-1}\right),\qquad1\le j\le m.
\]

The exceptional set is the disjoint union

\[
 (A+p^{s-w})N\equiv j\pmod {p^{K+\delta_j}},
                         \qquad1\le j\le m.
 \tag{11a}
\]

For proof, let `x=(A+p^(s-w))N` and choose `j` maximizing
`v_p(x-j)` among `1<=j<=m`. The remaining factor valuations are at
most the valuations of `j-i`, so
`v_p binom(x-1,m)<=v_p(x-j)-delta_j`. If the left side is at least
`K>L`, the closest factor has valuation greater than every nonzero
difference `j-i`. Consequently all the other factor valuations equal
those differences, and equality holds in the displayed bound.
Conversely every congruence in (11a) has that property. Distinct `j`
cannot share a residue modulo `p^K`, because `m<p^K`; hence the union
is disjoint. Its exact `p`-adic density is
`sum_(j=1)^m p^(-K-delta_j)`. This bounds the failure set by one formula
across all supports covered by the criterion.

## 5. Complete scalar leaves on the closed disk

This section specializes to rational `c=gamma+q-1 in Z_p`, original multiplication
part `F(z)=z`, and original seed 1. All auxiliary coefficients are in
`Z_p`, while the evaluated roots and relation coefficients may lie in
ramified extensions.

At a stage with multiplication part `F(z)=z+sum a_r z^(p^r)`, the
Artin--Hasse gauge and residue `j in {0,...,p-1}` give

\[
 e\equiv c-j\pmod p,\quad0\le e<p,\quad
 b=a_1-1\in\mathbb Z_p^\times,\quad y=bz^p/p,
\]
\[
 A=e+pL,\quad
 L=\frac{c-j-e}{p}-yD_y+y+
 \sum_{r\ge1}\frac{(a_{r+1}-1)p^{p^r-1}}{b^{p^r}}y^{p^r}.
 \tag{12}
\]

Writing `E_p(z)^(-1)h(z)=sum g_m z^m`, the actual new seed is

\[
                h'_j(y)=\sum_{\ell\ge0}g_{p\ell+j}(py/b)^\ell.
 \tag{13}
\]

Its constant is a unit. Every coefficient at degree `ell>=1` has
valuation at least `ell`, so it is a **restricted** series on the
closed unit disk. The high coefficients of the multiplication part
in (12) have valuation exactly `p^r-1`; that series is restricted too.

Perform all `s` residue splits, including the final size-`p` to
size-one split. Even when `v_p(u)=b_s`, before this last split the inner
root has valuation `1/p>0`; hence its Artin--Hasse gauge converges.
More generally, after `r<s` splits the cyclic parameter has valuation
at least `b_(s-r)`, and its cyclic roots have valuation at least
`b_(s-r)/p^(s-r)>0`. Thus every gauge used before a split is evaluated
strictly inside its convergence disk.
The last seed and multiplication series are restricted by (12)--(13),
so evaluation on the resulting unit disk is legitimate. No extra
Artin--Hasse gauge is evaluated at a unit argument.

The shared final variable is

\[
                    y=\nu\,u/p^{b_s}=\nu v,
                    \qquad\nu\in\mathbb Z_p^\times.
 \tag{14}
\]

All scaling units are independent of the residue choices: (12) changes
the multiplication part independently of `c` and `j`. A leaf has index
`j=j_0+p j_1+...+p^(s-1)j_(s-1)`, `0<=j<q`, and the original operator
on it is `e_j+qL_j`, with

\[
 e_j\equiv c-j\pmod q,\quad0\le e_j<q,\qquad
 L_j(0)=(c-j-e_j)/q.
\]

There is one resonant leaf `j_*=c mod q`. Discard that equation and
retain the other `q-1`. If `r_j=v_p(e_j)<s`, absorb

\[
 p^{-nr_j}\Pi_{n,q}(e_j+qL_j)
  =\prod_{\ell=0}^{n-1}
       (e_j/p^{r_j}+p^{s-r_j}(L_j+\ell))
\]

into its actual seed. Each factor is an integral Euler operator reducing
to a scalar unit. The absorbed seed `H_j` remains restricted and reduces
to a constant unit. Finite Euler operations preserve restrictedness:
Euler differentiation multiplies coefficients by integers and multiplication
of restricted series is closed. Consequently `H_j` has a restricted
inverse on the entire closed unit disk.

Explicitly, Euler differentiation and multiplication by the final
restricted `F` have Gauss operator norm at most one. Each factor in the
last product is therefore bounded and integral; because `s-r_j>=1`,
its reduction is the scalar `e_j/p^r_j`. Applying finitely many such
factors preserves both coefficient decay and constant-unit reduction.
Hence `H_j=H_j(0)+pK_j` with `K_j` restricted and `H_j(0)` a unit.
The inverse is
`H_j(0)^(-1) sum_(m>=0)(-pK_j/H_j(0))^m`, convergent in Gauss norm.
This proves invertibility on the entire closed disk, rather than only
formal invertibility from a unit constant coefficient.

Let `a_j=c-j` and put

\[
 \tau_{jk}(v)=
 \left.\frac{(L_j-L_j(0))^kH_j}{H_j}\right|_{y=\nu v}.
\]

These belong to `v Z_p< v >`. Their first member satisfies
`tau_(j,1)=nu*v mod p`, because the seed reduces to a constant and
`L_j-L_j(0)` reduces to `y-yD_y`. Every actual collision polynomial
satisfies the exact nonresonant leaf equations

\[
 \mathcal F_j(Q)=Q(a_j)+
       \sum_{k=1}^{q-1}q^kQ^{[k]}(a_j)\tau_{jk}(v)=0.
 \tag{15}
\]

All outer parameters and coefficients of `Q(X;v)` remain constant under
the inner Euler operators. These operators act on entire inner series
before evaluation; no characteristic-zero quotient action is assumed.

## 6. Uniform interpolation and its first variation

Set

\[
 Q_0(X)=\prod_{j\ne j_*}(X-a_j),\quad
 G_j(X)=Q_0(X)/(X-a_j),\quad w_j=G_j(a_j).
\]

All nonzero node differences have valuation at most `s-1`. Thus, for
every `k>=1`,

\[
 v_p\!\left(q^kG_\ell^{[k]}(a_j)/w_j\right)\ge k,
 \qquad
 v_p\!\left(q^kQ_0^{[k]}(a_j)/w_j\right)\ge s+k-1.
 \tag{16}
\]

For the first bound, a divided derivative removes exactly `k` factors
from the product of node differences, including the factor `a_j-a_ell`
when `ell!=j`. For the second it removes `k-1`, since `Q_0(a_j)=0`.
Each reciprocal difference costs at most `s-1`; multiplication by
`q^k=p^(sk)` proves (16).

Write a monic degree-`q-1` polynomial as

\[
                   Q=Q_0+\sum_j b_j(v)G_j.
\]

The leaf equations, divided by `w_j`, have matrix
`I+p v U(v)` with `U in Mat_(q-1)(Z_p< v >)` and right-hand side in
`q v Z_p< v >^(q-1)`. The inverse exists in the same restricted ring,
by the `p`-adically convergent geometric series. Hence

\[
 b_j(v)\in qv\mathbb Z_p\langle v\rangle,
 \qquad b_j(v)/(qv)\equiv-\nu\pmod p.
 \tag{17}
\]

The second statement follows because after division of the forcing by
`qv`, only its `k=1` term survives: `Q_0'(a_j)/w_j=1` and
`tau_(j,1)/v=nu mod p`. The matrix is the identity modulo `p`.
Consequently

\[
           Q(X;v)\equiv Q_0(X)-\nu qv Q_0'(X)
                  \pmod {pqv\mathbb Z_p\langle v\rangle[X]}.
 \tag{18}
\]

Theorem 1 and reduced-degree noncollision ensure that every primitive
actual collision polynomial has unit leading coefficient. Normalize it
to be monic. The invertible leaf system shows that it equals this one
specific restricted family `Q(X;v)` at its actual parameter.

The numbers `-qn-a_j` run through the nonzero residue classes modulo
`q`. Therefore

\[
 v_p(Q_0(-qn))=B_s,\qquad
 v_p(G_j(-qn))\ge B_s-s+1,
\]

and the transfer gives the uniform evaluation estimate

\[
 Q(-qn;v)=Q_0(-qn)+p^{B_s+1}v A(v),
                     \quad A\in\mathbb Z_p\langle v\rangle.
 \tag{19}
\]

This is an all-parameter estimate on the closed disk. It is stronger
than merely asserting that the actual kernel is one-dimensional.

## 7. The remaining scalar functional

For any polynomial `P` of degree at most `q-1`, expand its nonconstant
coefficients in the integral triangular basis `ell_0,...,ell_(q-2)`:

\[
 P(X)-P(-qn)=\sum_j d_j\ell_j(X),\qquad
 \mathcal T_h(P)=\sum_j d_j(C+q-j)_{h-q+j},
 \tag{20}
\]

where a summand is zero when `q-1-j>=h`. This is the integral linear
functional `sum d_j v_(q-1-j)` from the actual mixed equation. Put

\[
 T_0=\mathcal T_h(Q_0),\quad
 D=\min_j v_p(\mathcal T_h(G_j)),\quad
 t_0=v_p(T_0),
 \tag{21}
\]

with `D=+infinity` if all terms vanish.
Every actual common root satisfies

\[
 Q(-qn;v)=(-1)^q\kappa_h p^{b_s}v\,
                       \mathcal T_h(Q(X;v)).
 \tag{22}
\]

The transfer gives

\[
 \mathcal T_h(Q)=T_0+p^{s+D}v R(v),
                  \quad R\in\mathbb Z_p\langle v\rangle.
 \tag{23}
\]

When `D` is finite, its first reduction is

\[
 R(v)\equiv-\nu\,p^{-D}\mathcal T_h(Q_0')\pmod p.
 \tag{24}
\]

**Theorem 3 (scalar nondegeneracy).** If `T_0!=0` and

\[
                          t_0<s+D,
 \tag{25}
\]

then the mixed pair is coprime. The same conclusion holds when
`t_0=s+D` and `v_p(T_h(Q_0'))>D`.
If `D=+infinity`, `T_h(Q)` is the constant `T_0`, and the conclusion
holds whenever `T_0!=0`. If that constant is zero, (19) already forbids
the mixed equation.

**Proof.** In the finite-`D` case, (23) and the stated condition make
`T_h(Q)` have valuation exactly `t_0` on the whole integral disk,
and its difference from `T_0` is in `p^(t_0+1) v Z_p< v >`.
By (19), any common root must have

\[
 v_p(v)=-s-v_p(\kappa_h)-t_0=:r.
\]

If `r<0`, there is no integral solution. Otherwise `r` is an integer.
Substitute `v=p^r w` in (22) and divide by `p^B_s`. The constant and
linear coefficients are units; all remaining coefficients are divisible
by `p`. The resulting restricted series is injective on the closed
integral disk: its divided difference is a unit. Thus there is at most
one distinct common root. The constant-functional cases have the same
argument or an immediate contradiction.

Since both polynomials have rational coefficients, every conjugate of
a common root would be another common root in `C_p`. Uniqueness forces
that root to be rational. A real common root is impossible: one of the
two moment powers is even. For every positive rational shape `beta`, the functional
`M_beta(X^j)=(beta)_j` has the exact Gram factorization

\[
 (\beta)_{i+j}=\sum_{k\ge0}L_{ik}\,k!(\beta)_k\,L_{jk},
 \quad L_{ik}=\binom ik\frac{(\beta)_i}{(\beta)_k}
              =\binom ik(\beta+k)_{i-k}\quad(k\le i).
\]

Vandermonde's identity proves the formula (or prove it at every positive
integer beta from the factorial identity, then use polynomial identity
in beta). The triangular matrix has unit diagonal and its diagonal
weights are positive rationals. Thus every nonzero rational polynomial
has a strictly positive `M_beta` square. The exact identities
`T_(n,gamma,q)(t)=M_gamma((1+tX^q)^n)` and the upper counterpart at
`beta=gamma+h` give the contradiction. Rational Bezout descent supplies
the characteristic-zero field statement. ∎

**Uniform offset-two corollary.** For every prime power `q>=2`, every
positive rational `gamma` integral at `p`, and every `n>=1`, the pair
with `h=2<=q` is
coprime. Indeed `T_2` is, up to sign, the leading coefficient of `P`:
it is a constant unit on monic `Q`, and is zero on every `G_j`.
The offset-one functional is identically zero and gives the same
conclusion. These conclusions concern all positive rational shapes integral at `p`; they are
not merely a list of ordinary supports.

## 8. Uniform offset-three quadratic reduction

For `q>=3`, the actual offset-three functional satisfies

\[
 \mathcal T_3(G_j)=-(-1)^q,\qquad
 \mathcal T_3(Q_0)=(-1)^q(2(\gamma+qn)+4+j_*).
 \tag{26}
\]

These identities follow by taking the top two coefficients in (20).
For example, the coefficient of `X^(q-2)` in `Q_0` is
`-(q-1)(gamma+q-1)+q(q-1)/2-j_*`; substituting the top two coefficients
of the `ell_j` basis gives (26).

If `v_p(2(gamma+qn)+4+j_*)<s`, Theorem 3 proves coprimality.
In every remaining case, the mixed pair has **at most two distinct
common roots**. To see this, (18) and (26) give

\[
 \mathcal T_3(Q)=(-1)^q q
      [z+\nu(q-1)v+p vR(v)],
 \quad z\in\mathbb Z_p,
 \quad R\in\mathbb Z_p\langle v\rangle.
\]

Normalize (22) by `p^B_s`, and set `e=v_p(kappa_3)+2s`.
If `e>0`, its reduction is a nonzero constant and there is no common
root. If `e<=0`, clear its scalar denominator by multiplying by
`p^(-e)`. The result is an integral restricted series whose coefficient
of `v^2` is a unit and every coefficient of degree greater than two is
divisible by `p`. Hence it has at most two zeros on the integral disk.

For completeness, the elementary restricted-series root bound used
here needs no unramified-root assumption. If an integral restricted
series has largest unit coefficient at index `d`, division by any
integral root `alpha` gives a restricted quotient whose largest unit
index is `d-1`: its coefficients are
`b_i=sum_(k>=i+1) a_k alpha^(k-i-1)`. Iteration bounds the number of
zeros, with multiplicities, by `d`.

Thus every irreducible common factor has degree at most two. The
retained Theorem `thm:prime-power-lower-one` additionally proves
squarefreeness of `T_(n,gamma,q)` for positive rational `gamma` integral
at `p`; under that inspected input the ordinary gcd itself has degree
at most two. A degree-two possibility is not excluded by this theorem.
No quadratic obstruction is promoted to full coprimality without an
additional argument.

## 9. Remaining structural obstruction

The scalar-leaf transfer is complete at every prime-power gap, including
the closed root radius. Its remaining uncertainty is the actual upper
functional, not generic cyclicity. For general `h`, conditions (21),
(23), and (24) give a finite algebraic description of the leading
response. When `t_0>s+D` and `T_h(Q_0')` is also divisible by `p^(D+1)`,
both the constant and first response vanish at this precision. The
leaf rank theorem alone supplies no nonvanishing next response.

This can happen for actual parameters even when the coarse valuation
criterion also fails. Take `p=q=5`, `s=1`, `h=5`, `n=1`, `gamma=11`.
The exact upper functional has

\[
 \kappa_5=\frac{10}{21\cdot22\cdot23\cdot24\cdot25},
 \qquad T_0=-32400,
\]
\[
 (T_5(G_j))_j=(1146,1184,1224,1266),\qquad
 T_5(Q_0')=4820.
\]

Thus `v_5(kappa_5)=-1`, `D=0`, `t_0=2>s+D=1`, and
`v_5(T_5(Q_0'))=1>D`. The scalar nondegeneracy conditions cannot be
declared automatic. This is a failure of the criterion, not an SFC
counterexample: the lower polynomial is linear and direct exact
substitution in the upper polynomial gives a nonzero value.

Nor does a degree-two bound settle SFC: irreducible quadratics may have
nonreal conjugate roots, so the rational sum-of-squares exclusion does
not apply merely from a two-root count. The adjacent-powers-of-three
trace identity is an additional arithmetic mechanism, not a consequence
of full leaf rank for arbitrary offsets.

Offsets `h>q` also fall outside the mixed-rank theorem used here; their
upper operators require a further exact reduction. This is material
for ordinary binomials with lower exponent greater than their gap.

At a gap with more than one prime factor, splitting only its `p`-power
part leaves a nontrivial cyclic block. The full scalar proof above does
not transfer by replacing that block's actual seed with a constant.
The explicit gap-six rank failures recorded in the previous investigation
are negative controls for such a shortcut.

An unrestricted SFC closure therefore still needs an independent theorem
for the remaining upper response or the remaining non-prime-power
blocks. This note supplies unconditional partial theorems and exact
interfaces; it does not assume that missing assertion under a new name.

## 10. Exact replay, overlap, and review scope

Run:

```sh
python3 research/factorial/image-rigidity/verify_sfc_prime_power_scalar_transfer.py
```

The verifier uses the Python standard library and exact rational
arithmetic. It refuses optimized Python that disables assertions and
checks SHA256 pins for this note and all three retained proof sources.
Its [receipt](sfc_prime_power_scalar_transfer_verification.json) records
finite controls: actual residue splits, absorbed seeds, interpolation
matrices and their first variations, sharp stage examples, rational
Gamma Gram identities, the ordinary valuation formulas, and the
explicit difficult-response example. These controls check algebra and
implementation; the unrestricted quantifiers follow from the proofs.

The final-stage convergence and inverse argument has been independently
reconstructed internally by both the root reviewer and a separate
factorial reviewer. The latter requested the explicit closed stage
radius in Section 1; it is now stated. Root review requested the
leading-unit argument and the rational-shape scope; both are explicit.
These are internal checks, not external review or formal verification.

The retained sources already prove prime-power squarefreeness and the
offset-one conclusion on their stated domains. The present note uses
those exact inputs and records a uniform collision normalization,
closed-disk transfer and upper-response criteria. It makes no priority
claim and does not assert unrestricted SFC. In particular, the remaining
higher-response and non-prime-power issues in Section 9 are still
mathematical obligations.
