# The fixed column-6 word has finitely many Selmer specializations

For the specific coefficientwise column-6 continuation defined below,

\[
 \boxed{\#\{t\in\mathbf Q:\ E_t\text{ smooth},\ \beta(t)\text{ defined},\
                   [\beta(t)]\in\operatorname{Sel}_2(E_t)\}<\infty.}
\]

Consequently its labelled 2-cover is rationally soluble for only finitely
many eligible rational parameters. This strengthens the
[earlier good-divisor obstruction and genus bound](CARRIER_CLOSURE_AND_DEPENDENCY_CONTINUATION_2026-09-12.md)
by adding the uniform finite-prime argument proposed in the conversation.
It closes this fixed-output continuation branch. It does not constrain a
constructor that chooses new principal dependencies at different parameters.

The proof is unconditional, using the retained exact polynomial identities
and Faltings's theorem. It computes neither the finite parameter set nor a
useful numerical bound. No resultant is factored, no genus-five twist is
searched, and column 7 remains untouched.

## Fixed inputs and the specialization domain

The [expanded definition](../../artifacts/generated-results/elliptic-curves/frozen_dependency_continuation_v1/definition.json)
contains the entire 1,676-atom word, all sixteen generic section maps, the
field equation, normalization constants and source hashes. The
[original Bézout packet](../../artifacts/generated-results/elliptic-curves/frozen_dependency_continuation_v1/result.json)
certifies its polynomial coprimalities. These are also included in the
[complete proof packet](../../artifacts/generated-results/elliptic-curves/fixed_word_selmer_finiteness_v1/proof-packet.zip),
with the [member manifest](../../artifacts/generated-results/elliptic-curves/fixed_word_selmer_finiteness_v1/proof-packet-manifest.json).

For the frozen compact MW16-05 equation `Y^2=X^3+A(t)*X+B(t)`, put `u=289/2`
and

```
c1(t) = (u^4*A(t)+27)/81,
c0(t) = (u^6*B(t)+3*u^4*A(t)+27)/729,
f_t(Z) = Z^3+Z^2+c1(t)*Z+c0(t),
z = (u^2*X-3)/9,       w = u^3*Y/27.
```

Thus `w^2=f_t(z)`. Write `L_t=Q[theta]/(f_t)` on each smooth rational fibre;
this is an etale cubic algebra, possibly a product of fields. No irreducibility
assumption on every specialized cubic is needed.

For every index in the fixed set `I6`, let

```
alpha_i(t) = a_i+b_i*theta,
n_i(t) = a_i^3-a_i^2*b_i+c1(t)*a_i*b_i^2-c0(t)*b_i^3,
pi_i(t) = n_i(t)*alpha_i(t).
```

The generic factors are `gamma_j(t)=s_j*(z_j(t)-theta)`, where each `s_j`
is a positive rational square and `w_j(t)^2=f_t(z_j(t))`. The continued word is

```
beta(t) = product_(i in I6) pi_i(t) * product_(j in J6) gamma_j(t),
J6 = [3,5,6,7,9,10].
```

All `n_i` and `w_j` are nonzero polynomials. The eligible affine domain is

```
disc_Z(f_t) != 0,
n_i(t) != 0 for every i in I6,
w_j(t) != 0 for j in J6.
```

The word is then invertible in `L_t`, and its norm is a rational square
because `Norm(pi_i)=n_i^4` and `Norm(gamma_j)=s_j^3*w_j^2`. This supplies its
usual 2-descent class and cover. The omitted polynomial zeros are a finite
set, so allowing a separately justified extension of the class at those
parameters would not change the finiteness conclusion. The theorem concerns
`t in Q`; adding the projective parameter at infinity adds at most one fibre.

The distinguished factor is

```
alpha_2(t) = -6324860115183736256 + 286*theta,
n_2(t) = lambda*N(t),       lambda = -10793861/54.
```

Here `N` is the fixed primitive integer polynomial of degree 12 in the
definition packet. It is squarefree. Its full coefficients are preserved;
neither interpolation nor a change of generator is performed.

## Exact binding of the finite exceptional-prime set

Let the rational polynomials `q` range over

```
N',             disc_Z(f_t),
n_i  (i in I6, i != 2),
w_j  (j=0,...,15).
```

There are 1,693 of them. For each let `D_q` be the positive least common
multiple of its coefficient denominators, and set `Q_q=D_q*q in Z[t]`.
No primitive-content normalization is applied. Define the integer

\[
 R_q=\operatorname{Res}_t(N,Q_q)=\det\operatorname{Sylvester}(N,Q_q).
\]

The [finite-prime binding](../../artifacts/generated-results/elliptic-curves/fixed_word_selmer_finiteness_v1/sigma-binding.json)
specifies every `D_q`, every expanded integer polynomial `Q_q`, the exact
Sylvester-determinant definition of `R_q`, and its nonzero residue modulo
`1000003`. Resultant values and their product remain unexpanded integer
expressions; this still specifies exact integers, not unspecified witnesses.
The existing complete-product Bézout certificate is retained and checked
after denominator clearing. An independent SymPy replay recomputes every
resultant residue from the integer coefficients.

Let `C` be the product of the distinct positive integer guard values listed
in that same packet. There are **105 integers, not 105 asserted primes**.
The product is explicitly stored as a 591-digit integer. Its factors include:

- `2`, `3` and the nonzero leading coefficient of `N`;
- numerators and denominators of `lambda`, `b_2` and the coordinate scale `u`;
- coefficient denominators of the cubic, all atoms and all generic maps;
- numerators and denominators of the generic scales `s_j`;
- every denominator-clearing multiplier `D_q`.

Define, without factoring,

\[
 M=C\prod_q R_q\in\mathbf Z\setminus\{0\},\qquad
 \Sigma=\{p\text{ rational prime}:p\mid |M|\}.
\]

The checks give **`M = 700418 modulo 1000003`**, so this exact determinant
expression is nonzero. Thus `Sigma` is a specific finite set defined by
complete integer data. Its prime divisors are not enumerated. Including
the derivative and all sixteen ordinate resultants enlarges `Sigma`
conservatively and keeps the binding aligned with the prior certificate.

This is neither the control fibre's bad-prime set nor the diagnostic primes
through 1009. Fresh fibres can have bad primes outside `Sigma`; the resultants
ensure those primes cannot coincide with the relevant divisors of `N(t)`
when the parameter is integral there.

## Finite-prime lemma

Let `t=m/n` with coprime integers `m,n` and `n>0`, in the eligible domain.
Set

\[
 \mathcal N(m,n)=n^{12}N(m/n)
     =\sum_{j=0}^{12}N_jm^jn^{12-j}\in\mathbf Z\setminus\{0\}.
\]

**Lemma.** If `p` is outside `Sigma` and `k=v_p(mathcal N(m,n))` is odd,
then the specialization of `beta(t)` fails the local point-Kummer condition
at `p`. In particular its 2-cover has no `Q_p`-point.

**Proof.** Since `mathcal N` is an integer, odd `k` here is positive. If `p`
divided `n`, coprimality would give

```
N_hom(m,n) = lc(N)*m^12 != 0 modulo p,
```

contradicting `p | mathcal N`. Hence `p` does not divide `n`, the parameter
is `p`-integral, and `v_p(N(t))=k`.

Every scalar required for the chosen model and elements is integral at `p`;
the designated normalization and slope scalars are units. Because `p` does
not divide `Res(N,Q_q)`, no `Q_q` can vanish at the same residue as `N`.
Indeed, a common residue root would make the reduced Sylvester determinant
zero. Since `D_q` is a unit, every `q(t)` in the list is also a unit.

In particular `disc_Z(f_t)` is a unit. The monic integral cubic gives good
elliptic reduction at this odd prime, and its cubic algebra is unramified.
Also all other `n_i(t)` are units. The corresponding `alpha_i(t)` are
integral elements with unit norm, hence units in every component. The
generic `s_j`, `z_j(t)` and `w_j(t)` are integral, with `s_j` and `w_j(t)`
units, so their norm identities make every `gamma_j(t)` a unit as well.

Since `n_2(t)=lambda*N(t)` and `lambda` is a unit, its valuation is `k`.
The slope `b_2` is a unit, so `r=-a_2/b_2` is integral and is a simple root
of the cubic modulo `p`. After an unramified splitting extension, precisely
one of the three factors `a_2+b_2*theta_l` is nonunit. The product has
valuation `k`, so their valuations are `(k,0,0)`. Multiplication by the
base scalar `n_2(t)` gives

```
valuations of the complete beta(t): (2k,k,k).
```

This notation uses the three geometric sheets. An irreducible residual
quadratic contributes one degree-two unramified prime of valuation `k`.
No split-cubic assumption was used.

Local point-Kummer classes at good odd reduction are unramified, hence have
even valuations in each cubic component. This is the standard descent
criterion; see [Poonen–Schaefer, Propositions 12.2–12.5][PS]. The elementary
root-factor proof is also given in the earlier note. For odd `k`, the
pattern `(0,1,1)` violates that condition. This proves the lemma.

Squares change valuations by even numbers. Generic point-Kummer corrections
also cannot remove the obstruction. More generally, at each eligible fibre,
multiplication by any rational-point Kummer class preserves membership in
the Selmer subgroup, even if the choice of correction varies with `t`.

## Finite twist containment and Faltings

If `[beta(t)]` is in `Sel_2(E_t)`, the lemma implies

\[
 v_p(\mathcal N(m,n))\equiv0\pmod2\qquad(p\notin\Sigma).
\]

Since `n^12` is a rational square, there is a signed squarefree integer
`delta`, supported on `Sigma`, for which

\[
 N(t)=\delta y^2,\qquad y\in\mathbf Q.
\]

There are only finitely many such `delta`: both signs and all subsets of
the finite prime set are allowed. This is a necessary containment, not a
converse criterion for Selmer membership.

For every such `delta`, the curve `delta*y^2=N(t)` is geometrically integral.
Its smooth projective normalization has twelve distinct finite geometric
branch points and no branch point at infinity. Riemann–Hurwitz gives
`2g-2=-4+12`, hence genus five. This includes twists with no rational point.
The projective normalization accounts for points at infinity explicitly.

[Faltings's theorem][F] gives finitely many rational points on each of
these curves. A finite union of their images on the parameter line is
finite. It contains every eligible Selmer specialization of this fixed
word and therefore every rationally soluble specialization of its labelled
cover. No assumption about GRH, BSD or finiteness of Sha is involved.

The earlier genus bound alone did not prove this statement: it concerned
each fixed parameter cover. The finite set `Sigma` now confines all locally
soluble specializations, across all denominators, to finitely many twists.

## Verification and research boundary

The independent replay checks the full 1,676-atom input and specializations,
all sixteen generic section identities, all 1,693 denominator clearings and
modular resultants, the scalar guards, the scaled prior Bézout identity, and
the homogeneous denominator identity. It passes in 1.242 CPU seconds; the
binding producer used 0.335 CPU seconds. These are component times. Both
use one process and peak below 159 MB RSS. The written deduction invokes
the cited Faltings theorem; it is not a proof-assistant formalization.

```sh
sage -python research/elliptic-curves/rank-jump/verify_fixed_word_selmer_finiteness.py
```

The replay does not overwrite its receipt, expand or factor resultants,
enumerate `Sigma` or twists, search parameters, or recompute number-field
ideal factorizations. The complete proof packet includes the previously
missing expanded definition and Bézout payload, both new scripts and their
allowlisted replay inputs. Restore its research-relative paths only in an
empty replay checkout.

The fixed output is unsuitable as an infinite-family template. The reusable
object remains the successful adaptive map from an equation and its generic
subgroup to a newly constructed principal dependency. Its atom word can vary
with the parameter; this theorem puts no finiteness restriction on that
constructor's possible outputs. Rational lifting, rank independence,
strictness and ordinary ideal-class independence remain separate endpoints.
The certified soluble control at `3/17` remains one member of the finite set;
the proof does not identify the other members.

The carrier bank remains closed, the 182 exact misses and interrupted
candidate remain preserved, and no new collector is launched. The theorem
does not bound high-rank fibres of MW16-05, rule out other continuations of
the marked class, or exclude further isolated soluble outputs of this word.

[PS]: https://math.mit.edu/~poonen/papers/descent.pdf
[F]: https://link.springer.com/article/10.1007/BF01388432
