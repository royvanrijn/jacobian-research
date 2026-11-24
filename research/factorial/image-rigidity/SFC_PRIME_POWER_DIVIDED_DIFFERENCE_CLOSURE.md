# Divided differences and uniform prime-power support closure

## Result and scope

Let `p` be any prime and let `0<=w<=s`, with `s>=1`. Set

\[
 a=p^w,\qquad q=p^s.
\]

For every integer Gamma shape `1<=beta<=a`, every nonzero binomial

\[
                  X^a(A+BX^q)
\]

has a nonzero Gamma-`beta` moment in every consecutive pair of positive
orders. In particular, ordinary factorial moments (`beta=1`) satisfy
this statement for **all supports `{p^w,p^w+p^s}` with `w<=s`**.
The phase coefficients may belong to any characteristic-zero field.
The theorem is uniform in `p,w,s`, the starting order, and the displayed
Gamma shapes; no support or order enumeration is used in its proof.

A second result supplies a nonvacuous arithmetic coprimality criterion
for every ordinary lower exponent `1<=a<q`. It removes the large
factorial valuation from the earlier threshold by proving that the
actual upper functional has corresponding divisibility.

The imported analytical object is the **actual** scalar family of
[SFC_PRIME_POWER_SCALAR_TRANSFER.md](SFC_PRIME_POWER_SCALAR_TRANSFER.md),
including its seeds and closed-disk estimates. The stronger criterion
proved in [SFC_PRIME_POWER_UPPER_RESPONSE.md](SFC_PRIME_POWER_UPPER_RESPONSE.md),
Theorem 1, is used with its precise hypotheses. This note proves the
functional identities and valuation bounds needed to apply it.
No arbitrary-seed assertion, complex positivity, or RH input is used.

## 1. The functional and its divided differences

Fix `q=p^s`, `n>=1`, positive rational `gamma` integral at `p`, and
`1<=h<=q`. Put

\[
 C=\gamma+q(n+1),\quad c=\gamma+q-1,\quad
 r=c\bmod q\in\{0,\ldots,q-1\},\quad a_j=c-j.
\]

Write

\[
 Q_r(X)=\prod_{0\le j<q,\ j\ne r}(X-a_j),\qquad
 G_j(X)=\frac{Q_r(X)}{X-a_j}\quad(j\ne r).
\]

Thus `Q_r` is the retained `Q_0`; the subscript here records the omitted
node. Let

\[
 \ell_k(X)=(X+qn)\prod_{i=0}^{k-1}(c-X-i),\quad0\le k\le q-2.
\]

For a polynomial `P` of degree at most `q-1`, define the actual integral
functional by

\[
 P(X)-P(-qn)=\sum_{k=0}^{q-2}d_k\ell_k(X),\qquad
 \mathcal T_h(P)=\sum_{k=q-h}^{q-2}
             d_k(C+q-k)_{h-q+k}.
 \tag{1}
\]

An empty sum is zero. For a polynomial `f`, the notation
`f[x_0,...,x_m]` denotes its order-`m` divided difference, including its
polynomial continuation when nodes coincide. Put

\[
             f_h(Y)=(Y)_h,\qquad A_j=j-q+1.
\]

**Lemma 1 (exact functional formula).** For every displayed parameter,

\[
 \boxed{\mathcal T_h(Q_r)=(-1)^q f_h[C,C-q,A_r]},
 \qquad
 \boxed{\mathcal T_h(G_j)=(-1)^{q+1}
                   f_h[C,C-q,A_r,A_j]}.
 \tag{2}
\]

These are characteristic-zero polynomial identities, not merely
congruences or asymptotic formulas.

**Proof.** Substitute `Y=c-X` and set `M=c+qn=C-1`. Use the falling
Newton basis `Y^{\underline i}`. If

\[
             P(c-Y)=\sum_{i=0}^{q-1}e_iY^{\underline i},
\]

then

\[
 (M-Y)Y^{\underline k}=(M-k)Y^{\underline k}
                                  -Y^{\underline{k+1}}
\]

gives, by triangular elimination,

\[
 d_k=-\sum_{i=k+1}^{q-1}e_i
                  (C-k-2)^{\underline{i-k-1}}.
\]

Substitution in (1), followed by telescoping consecutive rising
products, yields

\[
 \mathcal T_h(P)=-\frac1q\sum_{i=q-h+1}^{q-1}e_i
 \bigl[(C+q-i)_{h-q+i}-(C-i)_{h-q+i}\bigr].
 \tag{3}
\]

One direct verification of this telescoping is to put `l=q-1-i`.
The coefficient of `-e_i` before telescoping is

\[
 \sum_{m=l+1}^{h-1}
 (C+m+1)_{h-1-m}
 (C-q+m-1)^{\underline{m-l-1}}.
\]

It equals
`[(C+l+1)_(h-l-1)-(C-q+l+1)_(h-l-1)]/q`.
This follows by subtracting adjacent reciprocal rising products, or
as a polynomial identity after clearing the displayed denominators.

The quotient of `Y^{\underline q}` by its root `Y-r` has Newton
coefficients

\[
 \frac{Y^{\underline q}}{Y-r}
 =\sum_{i=0}^{q-1}(r-i-1)^{\underline{q-i-1}}
                                      Y^{\underline i}.
\]

Consequently, for `P=Q_r`, putting `l=q-1-i` gives
`e_i=(-1)^(q-1)(A_r)_l`. Equation (3) becomes

\[
 \frac{(-1)^q}{q}\sum_{l=0}^{h-2}(A_r)_l
       \bigl[(C+l+1)_{h-1-l}-(C-q+l+1)_{h-1-l}\bigr].
\]

The elementary product-difference identity

\[
 \sum_{l=0}^{h-1}(A)_l(B+l+1)_{h-1-l}
                         =\frac{(B)_h-(A)_h}{B-A}
\]

turns this into `(-1)^q f_h[C,C-q,A_r]`: the omitted term `l=h-1`
cancels between the two brackets. Finally,

\[
             Q_r-Q_j=(j-r)G_j.
\]

Dividing the difference of the first formula in (2) by `j-r`, while
`A_r-A_j=r-j`, proves the second formula with its stated sign. ∎

For orientation, (2) gives `T_1=0`, gives the constant monic value
`T_2(Q_r)=(-1)^q`, and at offset three gives

\[
 \mathcal T_3(Q_r)=(-1)^q(2C-2q+r+4),\qquad
 \mathcal T_3(G_j)=(-1)^{q+1},
\]

recovering the retained formulas exactly.

## 2. A factorial divisibility bound for arbitrary integral nodes

For `m>=0` define

\[
 B_{p,m}(h)=\sum_{\ell\ge1}
          \max\{\lfloor h/p^\ell\rfloor-m,0\}.
 \tag{4}
\]

Only finitely many summands are nonzero.

**Lemma 2.** For arbitrary nodes `x_0,...,x_m in Z_p`,

\[
           v_p\bigl(f_h[x_0,\ldots,x_m]\bigr)
                                      \ge B_{p,m}(h).
 \tag{5}
\]

If `h<m`, the divided difference is zero and has infinite valuation.

**Proof.** Apply the divided-difference product rule to the ordered
product `f_h(X)=prod_(i=0)^(h-1)(X+i)`. Each nonzero order-`m` term
chooses `m` distinct linear factors on which to take one difference.
Those factors contribute one. The remaining factors form `m+1`
consecutive blocks, evaluated successively at `x_0,...,x_m`.
If their lengths are `L_0,...,L_m`, then

\[
                     \sum_i L_i=h-m.
\]

Every consecutive block of `L` translates of an element of `Z_p`
has valuation at least `v_p(L!)`: among any `L` consecutive translates
at least `floor(L/p^ell)` are zero modulo `p^ell`. Hence a term has
valuation at least `sum_i v_p(L_i!)`.

For each `P=p^ell`, write `L_i=P b_i+r_i`, `0<=r_i<P`. Then

\[
 \left\lfloor\frac hP\right\rfloor-\sum_i b_i
     =\left\lfloor\frac{m+\sum_i r_i}{P}\right\rfloor\le m,
\]

because `m+sum r_i <=(m+1)P-1`. Thus
`sum_i floor(L_i/P)>=max(floor(h/P)-m,0)`.
Sum over `ell` and then over product-rule terms using the nonarchimedean
valuation inequality. This proves (5), including coincident nodes by
the polynomial form of the product rule. ∎

In particular, with

\[
 D=\min_{j\ne r}v_p\bigl(\mathcal T_h(G_j)\bigr),\qquad
 B_p(h):=B_{p,3}(h),
\]

Lemma 1 gives the uniform structural estimate

\[
                         \boxed{D\ge B_p(h)}.
 \tag{6}
\]

This estimates the actual functional on all interpolation directions.
It does not replace those directions or their seeds by a generic
polynomial.

## 3. Consecutive moment criterion

The reviewed stronger scalar criterion states that the mixed pair is
coprime if

\[
 k+2s+D>0,\qquad k=v_p\left(\frac{q(n+1)}{(C)_h}\right).
 \tag{7}
\]

If all functional values on `G_j` vanish, the pair is coprime directly.
For clarity, the reason for (7) is recalled here. The retained scalar
equation has constant coefficient of valuation `B=(q-1)/(p-1)-s`.
Every coefficient of degree at least two has valuation greater than
`B` under (7). Only its actual linear coefficient can compete with
the constant. The resulting restricted series has at most one zero
on the entire integral root disk, by its unit divided difference.
Rational conjugacy forces a common root to be rational. The exact
Gamma Gram factorization excludes it because one of the two moment
powers is even. That use of positive rational squares is after
rational descent, and is not a positivity premise over the eventual
coefficient field.

Combining (6) and (7) proves the following criterion for all mixed
pairs in the scope of Section 1:

\[
                 \boxed{k+2s+B_p(h)>0\ \Longrightarrow\ \gcd=1.}
 \tag{8}
\]

## 4. Every lower power of the same prime is closed

**Theorem 3.** Let `p` be prime, `0<=w<=s`, `s>=1`,
`a=p^w`, `q=p^s`, and `1<=beta<=a` an integer. For every `n>=1`,

\[
  \gcd\bigl(T_{n,an+\beta,q},\
                  T_{n+1,a(n+1)+\beta,q}\bigr)=1
                         \quad\text{over }\mathbb Q[t].
 \tag{9}
\]

**Proof.** The case `w=0` is the retained offset-one theorem, so assume
`w>=1`. Set `N=n+1`. For `h=a` and `gamma=an+beta`,

\[
 C=(a+q)N+\beta-a.
\]

The denominator `(C)_a` contains exactly the factor `(a+q)N` at
offset `a-beta`. Cancel it with the numerator `qN`. Every other
factor is `(a+q)N+j`, with

\[
             -(a-\beta)\le j\le\beta-1,\qquad j\ne0.
\]

Since `0<|j|<a=p^w`, one has `v_p(j)<w<=v_p((a+q)N)`.
Thus its valuation is exactly `v_p(j)`. Consequently

\[
 k=s-v_p(a+q)-v_p((a-\beta)!)-v_p((\beta-1)!).
 \tag{10}
\]

Every binomial coefficient `binom(p^w-1,j)` is a `p`-adic unit:
the factors `p^w-i` and `i` have equal valuations for `1<=i<p^w`.
It follows that the last two factorial valuations sum to
`v_p((a-1)!)`. If

\[
 \delta=v_p(a+q)-w,
\]

then `delta=0` except when `p=2` and `w=s`, when `delta=1`. Hence

\[
              k=s-v_p(a!)-\delta.
 \tag{11}
\]

Let `b_w=(p^w-1)/(p-1)=v_p(a!)`. From (4),

\[
 B_p(a)=b_w-\sum_{\ell=1}^{w}
                             \min\{p^{w-\ell},3\}
                         \ge b_w-3w+2.
 \tag{12}
\]

The last inequality uses the final summand `1` and bounds each of
the other `w-1` summands by `3`. Therefore

\[
 k+2s+B_p(a)\ge3(s-w)+2-\delta\ge1.
\]

Criterion (8) proves (9). ∎

For the moment interpretation, the exact identity is

\[
 \Lambda_\beta\bigl((X^a(A+BX^q))^m\bigr)
       =A^m(\beta)_{am}T_{m,am+\beta,q}(B/A),
 \qquad\Lambda_\beta(X^j)=(\beta)_j.
\]

The displayed factors are nonzero for `AB!=0`. If either coefficient
vanishes, every positive moment of the resulting nonzero monomial is
nonzero. A rational Bezout identity for (9) remains valid after any
characteristic-zero field extension, proving the stated coefficient
scope.

For odd `p` and `w<s`, the ordinary case also falls within the retained
relative-gap-three theorem. Dyadic `w=s` is covered by the retained
dyadic unit theorem. These overlaps do not make the new proof depend
on those results. The divided-difference argument treats all primes,
including the diagonal `w=s`, by the same algebraic estimate.
For `p=2` and `w<=s-2`, the ordinary case is also within the
relative-gap-three theorem. The directions not covered by those two
overlap arguments are the odd-prime diagonal `w=s` and the dyadic
ratio `q=2a`; some initial instances already have other retained
proofs, including the diagonal prime case `s=w=1`. The integer Gamma
range and the uniform functional estimate are additional aspects of
the present theorem. No priority claim is made.

## 5. A nonvacuous threshold for every ordinary lower exponent

For integers `1<=a<q=p^s`, write `a=p^w A` with `p` not dividing `A`.
Put

\[
 R_p(a)=\sum_{\ell\ge1}\min\{\lfloor a/p^\ell\rfloor,3\},
 \qquad K_p(a,s)=3s-R_p(a).
 \tag{13}
\]

The retained exact ordinary cancellation and factorial valuation give

\[
 k=s-v_p(a!)-V_N,\qquad
 V_N=v_p\binom{(A+p^{s-w})N-1}{A-1}.
 \tag{14}
\]

Since `B_p(a)=v_p(a!)-R_p(a)`, criterion (8) proves

\[
               \boxed{V_N<K_p(a,s)\ \Longrightarrow\ \gcd=1.}
 \tag{15}
\]

Moreover `floor(a/p^ell)=0` for `ell>=s`, so
`R_p(a)<=3(s-1)` and **`K_p(a,s)>=3` for every `a<q`**, irrespective of the size of `v_p(a!)`. Thus this criterion
has a positive threshold throughout its full parameter range. It
does not assert that `V_N` is uniformly bounded when `A>1`.

The existing binomial-valuation period and residue formulas transfer
with this new `K`. For example, if `1<A<p`, failure of (15) is exactly
the union of the `A-1` distinct classes

\[
       (A+p^{s-w})N\equiv j\pmod {p^{K_p(a,s)}},
                                  \qquad1\le j<A.
 \tag{16}
\]

Indeed the denominator `(A-1)!` is a unit and at most one numerator
factor is divisible by `p`. If `A=1`, `V_N=0` for every order and the
criterion has no exceptional class, as in Theorem 3.

For general `A>1`, with `L=floor(log_p(A-1))`, failure depends only on
`N mod p^(K+L)`. Vandermonde and
`v_p binom(p^(K+L),j)>=K+L-v_p(j)>=K` for `1<=j<=A-1` prove that
period exactly as in the retained arithmetic discussion. A surviving
class is an unresolved class for this criterion, not an SFC
counterexample or a finite cutoff.

## 6. Checks and remaining boundary

Run from the repository root:

```sh
python3 research/factorial/image-rigidity/verify_sfc_prime_power_divided_difference.py
```

The [exact standard-library checker](verify_sfc_prime_power_divided_difference.py)
refuses optimized Python and verifies fixed SHA256 pins for this note and
both imported scalar-transfer proof notes. Its full
[receipt](sfc_prime_power_divided_difference_verification.json) records all
finite controls. It independently builds
the original triangular functional, compares it with (2), verifies
the contiguous-block divisibility and cancellation formulas, and
checks bounded actual moment pairs as algebra controls. The proofs
above, not those bounded tests, establish the unrestricted parameters.

This theorem does not close arbitrary `a<q` in the surviving classes
of (15), lower exponents `a>q`, gaps with several prime factors, or
unrestricted SFC. In particular an offset-three quadratic factor at
an arbitrary positive rational shape remains a separate issue.

The strictness of (8) is retained: replacing `>0` by `>=0` without
an additional scalar argument is not justified. General prime-power
root-degree bounds and residual split behavior do not by themselves
exclude an irreducible imaginary quadratic over `Q`.


The functional identities were independently reconstructed and replayed
internally by a separate reviewer, using complete homogeneous symmetric
polynomials for divided differences rather than division by node
differences. That reviewer also checked the four-block count and exact
Gamma cancellation; the root reviewer reconstructed the valuation
argument. These are internal checks, not external or formal verification.
