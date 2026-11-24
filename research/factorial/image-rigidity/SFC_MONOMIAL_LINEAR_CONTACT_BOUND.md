# Ordinary monomial–linear observables: uniform bounds, residue modules, and exact obstructions

Internal research note, 21 November 2026. The central proved result is the
all-exponent, all-starting-index moment window of length
`(2s-1)n-s+1`, obtained from the contact bound `(2n-1)s+1`.
Its initial window has the desired length `s`. The note also proves
a residue-module construction. It does **not** prove the desired
support-size window for arbitrary starting order or unrestricted SFC.
No existing status registry is changed by this note.

## 1. Setup and exact identities

Let `s>=1`, let `r_i>=0` be integers, and let all `c_i` be nonzero elements
of a characteristic-zero field. Repeated coefficient values are allowed
except in the residue-module section, which expressly assumes distinct
values. Inactive variables can be removed: they multiply the `m`th moment
by the nonzero scalar `(r_i m)!`.

Put

\[
 F=\prod_{i=1}^sY_i^{r_i}\sum_i c_iY_i,
 \quad P(z)=\prod_i(1-c_i z),
 \quad A(z)=\prod_i(1-c_i z)^{r_i},
 \quad a(z)=zA(z),\quad d=\deg a=1+\sum_i r_i.
\]

Define the formal primitive `H(0)=0`, `H'=A/P`, and
`B=H\circ a^{-1}`. When all `r_i>=1`, `H` is polynomial. Allowing zero
exponents requires only a rational derivative and formal integration at
zero. In all cases the formal series are normalized at zero.

Set

\[
 E=P\frac{A'}A=-\sum_i r_i c_i\prod_{j\ne i}(1-c_jz),
 \qquad D=P+zE.
\tag{1}
\]

Then

\[
 a'=H'D,\qquad B'(a)=\frac1D,\qquad
 \deg D=s,\qquad [z^s]D=(-1)^s\left(\prod_i c_i\right)d\ne0.
\tag{2}
\]

Direct moment expansion and formal Lagrange inversion give

\[
 g_m:=\frac{\mathcal L(F^m)}{m!\prod_i(r_i m)!}
 =(m+1)[t^{m+1}]B(t)
 =[z^m]A(z)^{-m}P(z)^{-1}.
\tag{3}
\]

Equivalently, the entire sequence has the exponent-free differential
representation

\[
 \boxed{g_m=\operatorname{Res}_{z=0}a(z)^{-m}\omega(z),
 \qquad \omega=\frac{dH}a=\frac{dz}{zP(z)}.}
\tag{4}
\]

The powers of `a` retain the exponents; only the differential `omega` is
independent of them. This is not a moment-preserving substitution.

## 2. A proved initial window, uniform in all exponents

**Theorem 1.** The first `s` positive ordinary moments of `F` cannot all
vanish. This holds for every `r_i>=0` and arbitrary nonzero `c_i` over a
characteristic-zero field.

**Proof.** If `g_1=...=g_s=0`, (3) gives
`B(t)=t+O(t^(s+2))`, hence `H-a=O(z^(s+2))`. But

\[
 (H-a)'=H'(1-D).
\]

The first factor is a unit at zero, and `1-D` is a nonzero polynomial of
degree `s`. Its vanishing order is at most `s`, contradicting the required
order at least `s+1`. QED.

This is a uniform initial-window theorem, not the all-starting-orders
Strong Factorial statement. No new literature-priority claim is made.

## 3. A proved contact bound independent of the exponents

**Theorem 2.** For every polynomial `b` of degree at most `n>=1`,

\[
 \boxed{\operatorname{ord}_0(B-b)\le(2n-1)s+1.}
\tag{5}
\]

Equivalently the same bound holds for `ord_0(H-b\circ a)`. Consequently
every moment window beginning at `n>=1` with length

\[
 W(n,s)=(2s-1)n-s+1
\tag{6}
\]

contains a nonzero moment. The desired SFC length is `s`; (6) supplies it
at `n=1` and is generally weaker for later starts.

**Proof.** For `k>=2`, define polynomials `M_k` by

\[
 M_2=-D',\qquad
 M_{k+1}=PD M_k'
 +\bigl(P'D-(k-1)ED-(2k-1)PD'\bigr)M_k.
\tag{7}
\]

Since `PH''/H'=E-P'`, differentiation of (2) proves inductively

\[
 B^{(k)}(a(z))=
 \frac{M_k(z)}{H'(z)^{k-1}D(z)^{2k-1}P(z)^{k-2}}.
\tag{8}
\]

All denominators in (8) are units at zero. Equation (7) shows

\[
 \deg M_k\le(s-1)+(k-2)(2s-1).
\tag{9}
\]

No `M_k` is identically zero. Otherwise `B^(k)=0`, so `B` is a polynomial
in characteristic zero, and (2) would give the polynomial identity
`D(z)B'(a(z))=1`. This is impossible because `D` is nonconstant.

If `K=ord_0(B-b)>=n+1`, differentiation `n+1` times and (8)–(9) give

\[
 K-n-1=\operatorname{ord}_0 M_{n+1}
 \le(2n-1)s-n.
\]

This is (5); smaller `K` already satisfies it. A zero moment window of
length `W(n,s)` gives `K>=n+W(n,s)+1=(2n-1)s+2`, a contradiction. QED.

This bound is substantive but does not use enough of the special root
profile to reach `n+s`. A higher-order scalar recurrence is not being
assumed. Complex coefficients and repeated coefficient values cause no
problem for this proof.

## 4. The logarithmic approximation problem

In this section assume the `c_i` are distinct, and put `alpha_i=1/c_i`.
The partial-fraction residues of `omega` are

\[
 \operatorname{Res}_0\omega=1,\qquad
 \beta_i=\operatorname{Res}_{\alpha_i}\omega
 =\frac1{\alpha_iP'(\alpha_i)}
 =-\frac{c_i^{s-1}}{\prod_{j\ne i}(c_i-c_j)}.
\tag{10}
\]

In particular `sum beta_i=-1`. Let
`b(u)=u+sum_(k=2)^n b_k u^k`, and define

\[
 \Phi(u)=\sum_{k=2}^n\frac{k}{k-1}b_k u^{k-1},\qquad
 \Lambda(z)=\sum_i(\beta_i-r_i)\log(1-c_i z).
\]

Then the exact differential identity is

\[
 \boxed{\frac{d(H-b(a))}{a}
 =\omega-\frac{da}{a}-d\Phi(a)
 =d\bigl(\Lambda-\Phi(a)\bigr).}
\tag{11}
\]

The logarithmic residues are `beta_i-r_i`, whose sum is `-d`. They are
not arbitrary free parameters. If `K=ord_0(H-b(a))>=2`, (11) gives

\[
 \operatorname{ord}_0(\Lambda-\Phi(a))=K-1.
\tag{12}
\]

Thus the desired contact bound becomes a very specific logarithmic
approximation bound. Replacing these prescribed residues by arbitrary
weights would strengthen the problem and is not justified.

## 5. An exact rank-s residue module and invertible shifts

Here assume **distinct nonzero** `c_i` and **positive** integers `r_i`.
Put `Q=zP`. For every integer `m>=1`, define

\[
 L_m R=QR'-mDR,
 \qquad
 V_m=\mathbb C[z]_{\le md+s-1}ig/
      L_m\bigl(\mathbb C[z]_{\le md-1}\bigr).
\tag{13}
\]

**Proposition 3.** `V_m` has dimension `s`, with basis represented by
`1,z,...,z^(s-1)`. Its residue map identifies it with

\[
 W=\{(v_0,...,v_s)\in\mathbb C^{s+1}:\sum v_i=0\}.
\tag{14}
\]

The map is

\[
 [f]\longmapsto
 \left(\operatorname{Res}_{\alpha_i}
       \frac{f(z)\,dz}{Q(z)a(z)^m}\right)_{i=0}^s,
 \qquad\alpha_0=0.
\tag{15}
\]

Moreover multiplication by `a` gives an isomorphism
`V_m -> V_(m+1)` and preserves the residue vector.

**Proof.** For `0<=j<=md-1`, the leading term of `L_m(z^j)` has degree
`j+s` and coefficient `[z^(s+1)]Q (j-md)`, which is nonzero. Triangular
reduction proves the asserted basis and dimension.

The rational differential in (15) has poles only at the listed points
and is `O(z^(-2)) dz` at infinity. Its residues sum to zero, and

\[
 \frac{L_mR}{Q a^m}\,dz=d(Ra^{-m}),
\]

so (15) is well defined. Conversely, a rational differential with every
residue zero has a rational primitive, by partial fractions. At zero its
primitive has pole order at most `m`; at `alpha_i` the order is at most
`m r_i`. Choose the primitive to vanish at infinity. It is therefore
`R/a^m` with `deg R<=md-1`, proving that the kernel is exactly the
denominator subspace in (13). Since both dimensions are `s`, this is an
isomorphism onto (14).

Finally `L_(m+1)(aR)=a L_m(R)`. The multiplication map is therefore well
defined, and its residue vector is unchanged. Equation (15) makes it an
isomorphism. QED.

These are vector spaces, **not quotient algebras**: `im L_m` need not be
an ideal. For example `a=Q=z(1-z)`, `D=1-2z`, `m=1` give
`L_1(1)=-1+2z`, but `z L_1(1)` reduces to the nonzero constant `-1/2`.
Multiplication by `a` works between the different spaces for the stated
differential identity; arbitrary polynomial multiplication does not
descend to one fixed space.

## 6. Exact obstruction to consecutive scalar observability

The module does **not** prove that every `s` consecutive scalar moments
determine its residue state. The following countercontrol lies inside
the exact positive-exponent, distinct-coefficient family.

Let

\[
 a=z(1-z)(1-tz)^2,\quad H'=1-tz,
 \quad \omega=\frac{dz}{z(1-z)(1-tz)},
\]

and use the first two independent residue coordinates

\[
 v_m=\left(\operatorname{Res}_0a^{-m}\omega,
            \operatorname{Res}_1a^{-m}\omega\right).
\]

Exact residues are

| `m` | First coordinate | Second coordinate |
|---|---|---|
| 0 | `1` | `1/(t-1)` |
| 1 | `2+3t` | `(5t-2)/(t-1)^4` |
| 2 | `6+15t+15t^2` | `3(12t^2-9t+2)/(t-1)^7` |
| 3 | `20+70t+112t^2+84t^3` | `2(143t^3-156t^2+65t-10)/(t-1)^10` |

Consequently

\[
 \det(v_1,v_2)=
 \frac{-15t^3(t-2)(5t^2-2t+2)}{(t-1)^7},
\tag{16}
\]

\[
 \det(v_2,v_3)=
 \frac{-252t^5(t-2)(12t^2-5t+5)}{(t-1)^{10}}.
\tag{17}
\]

At `t=(1+3i)/5`, the quadratic in (16) vanishes, whereas
`12t^2-5t+5=(1-t)/5` is nonzero. Thus `v_1,v_2` are dependent while
`v_2,v_3` are independent. This is not merely a permanent reduction in
the dimension of the observable subspace.

There is even a fixed nonzero linear combination of the two residue
cycles whose sequence has two consecutive zeros followed by a nonzero
term. Write `v_m=(x_m,y_m)` and set
`ell(v)=y_1 x-x_1 y`. Then

\[
 \ell(v_1)=\ell(v_2)=0,\qquad
 \ell(v_0)\ne0,\qquad\ell(v_3)\ne0.
\tag{18}
\]

For the specified complex embedding, `ell(v_0)=(144-42i)/625`.
After normalization this is a period sequence with zeroth term one and
a forbidden two-zero window. It is **not** the ordinary factorial
sequence: the ordinary sequence uses the distinguished first coordinate,
which has `x_1=2+3t!=0` here. No SFC counterexample is claimed.

What fails is the argument that rank `s`, invertible vector shifts, and
`s` consecutive scalar observations *automatically* imply scalar
zero-window propagation. This example does not refute all rank-s systems,
all scalar recurrences, or a proof that uses the distinguished cycle.
It shows that nonsingular consecutive observability must be proved and
is false in this natural residue realization at some complex parameters.

## 7. Correspondence countercontrols

The coordinate relation `u=H(z), v=b(a(z))` is generally an algebraic
correspondence, not a single-valued rational self-map of `u`.

Already `a=z(1-z)^2`, `H=z-z^2/2`, and `b(z)=z` have this issue. At
`z=1+w`,

\[
 u=\tfrac12-\tfrac12w^2,\qquad v=w^2+w^3.
\]

Thus near `u=1/2`, the two branches are
`v=-2(u-1/2) +/- (-2(u-1/2))^(3/2)`. Equal orders of ramification
make the first derivative finite but do not remove the branch.

More seriously, `H^{-1}(0)={0,2}` and `a(2)=2`. Therefore the
correspondence sends zero to both zero and two. Making the distinguished
zero branch repelling does not show that zero has a fixed forward orbit
under the whole correspondence. The polynomial proof excluding critical
value zero cannot be copied without additional control of branch choices.

Counting attracting fixed points by distinct critical values also fails
for unrestricted polynomial-type correspondences. For an integer `e>=2`
take the parameterization `u=z^e+epsilon^e`, `v=z^(e+1)`, with small
positive `epsilon`. The fixed equation becomes, on writing `z=epsilon w`,

\[
 \epsilon w^{e+1}-w^e-1=0.
\]

The implicit function theorem gives `e` solutions near the distinct roots
of `w^e=-1`. They determine distinct fixed input values, and their local
multipliers are `(e+1)z/e`, hence are attracting for small `epsilon`.
The forward projection has just one finite critical value, zero. Thus
different attracting branches can use the same critical value. This
example is a countercontrol for a general correspondence counting rule,
not a counterexample to the special SFC contact conjecture.

## 8. What remains, and verification boundary

The desired improvement is still
`ord_0(H-b\circ a)<=deg b+s`, uniformly for the special root profile.
Theorems 1–2 and Proposition 3 are proved progress toward that structure;
the observability and correspondence countercontrols identify what they
do not supply. Neither finite checks nor a dimension count is promoted
to the missing theorem.

The companion standard-library checker verifies the derivative recurrence
and moment normalization on exact rational cases; verifies the two
determinant identities as complete polynomial identities; performs the
specialization in `Q[t]/(5t^2-2t+2)` exactly; checks finite module reductions
and shift determinants; and checks the correspondence identities.

Run the standard-library checker from the repository root:

```sh
python -B research/factorial/image-rigidity/verify_sfc_monomial_linear_contact.py \
  --output /tmp/sfc-monomial-linear-contact.json
```

The [retained receipt](sfc_monomial_linear_contact_verification.json) records
the exact finite controls and SHA-256 hashes of this note and its checker.
The checker supports the written proof; finite test ranges do not establish
its unrestricted quantifiers by extrapolation.

The primary-literature check did not supply an applicable correspondence
Fatou bound. In particular Dinh–Kaufmann–Wu,
[arXiv:1808.10130](https://arxiv.org/abs/1808.10130), explicitly treats
equal forward and backward topological degrees, whereas this proposed
correspondence has unequal degrees. Siqueira,
[arXiv:1711.02365v2](https://arxiv.org/abs/1711.02365v2), studies the special
family `(w-c)^q=z^p`; no theorem from that paper is imported here. These
scope observations are not a claim that the desired theorem is absent
from the entire literature.
