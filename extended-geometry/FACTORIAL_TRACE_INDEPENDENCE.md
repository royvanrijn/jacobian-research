# Factorial trace independence and gamma-affine classification

## 1. Status and scope

This note upgrades the scaled-factorial rigidity theorem in
[`BINARY_GVC_FIRST_GHOST_SOURCE_COLLAPSE_AND_RAY_RIGIDITY.md`](BINARY_GVC_FIRST_GHOST_SOURCE_COLLAPSE_AND_RAY_RIGIDITY.md).
That theorem reconstructs one positive multiplicity partition from one
all-scale factorial product.  The result below first separates an arbitrary
finite sum of **signed** integer-affine factorial rays even when every ray is
multiplied by an exponential-polynomial sequence.  It then gives a complete
classifier for positive-integer-slope gamma products with arbitrary complex
offsets.  In both statements the coefficients may contain rational functions
of the scale.

These are characteristic-zero statements about exact hypergeometric
sequences.  It does not assert that a Hall--jet filtration exposes one fixed
scale-compatible packet or distinguish two marked states with the same
canonical gamma signature, classify factorial units or carries in positive
characteristic, or decompose a general holonomic sequence into finitely many
hypergeometric terms.  Section 7.3 records an exact reason that raw
characteristic-`p` valuations cannot inherit the characteristic-zero theorem.

## 2. Exponential-rational coefficients and canonical rays

Work with sequences modulo eventual equality.  Let

\[
 \mathcal E_{\mathrm{rat}}
 =\left\{
   \sum_{\lambda\in\Lambda}r_\lambda(n)\lambda^n:
   \Lambda\subset\mathbb C^\times\text{ finite},\quad
   r_\lambda(n)\in\mathbb C(n)
  \right\}.
\tag{2.1}
\]

This contains the usual exponential-polynomial sequences, for which every
`r_lambda` is a polynomial.  Roots of unity make the sequence ring have zero
divisors, so the conclusion below is stated as a separation property rather
than relying on abstract freeness terminology.

For a finite-support vector

\[
 v=(v_1,v_2,\ldots)\in\mathbb Z^{(\mathbb N_{>0})},
\]

put

\[
 \Phi_v(n)=\prod_{a\geq1}\Gamma(an+1)^{v_a}
          =\prod_{a\geq1}((an)!)^{v_a}.
\tag{2.2}
\]

Negative entries encode factorials in the denominator.  Thus multinomial,
binomial, radial, and denominator-cleared moment rays all lie in this class
after rational prefactors and ordinary exponentials have been removed.

> **Theorem 2.1 (factorial-trace separation).**
> Let `v^(1),...,v^(s)` be pairwise distinct finite-support integer vectors.
> If
> \[
>  \sum_{j=1}^s c_j(n)\Phi_{v^{(j)}}(n)=0
>  \qquad(n\gg0),
> \tag{2.3}
> \]
> with `c_j` in `E_rat`, then every `c_j` is the zero sequence.

The term *trace* is justified by the following immediate case.  Every entry
of `A^n`, for a fixed finite matrix `A`, is an exponential-polynomial sequence
by Cayley--Hamilton (equivalently, by Jordan form after extending scalars).
Consequently a finite identity

\[
 \sum_j\Phi_{v^{(j)}}(n)
 \operatorname {Tr}(B_jA_j^n)=0
\tag{2.4}
\]

splits canonically by the distinct vectors `v^(j)`.

## 3. Rational independence of dissimilar hypergeometric terms

A nonzero sequence `h` is hypergeometric when

\[
 q_h(n)=\frac{h(n+1)}{h(n)}\in\mathbb C(n).
\]

Two hypergeometric terms are **similar** when their quotient is rational.
The following standard lemma is Theorem 5.1 of
[Petkovšek's 1992 paper on hypergeometric recurrence solutions](https://doi.org/10.1016/0747-7171(92)90038-6).
Its short proof is included so that Theorem 2.1 has no external proof gap.

> **Lemma 3.1 (dissimilar hypergeometric terms separate).**
> Pairwise dissimilar hypergeometric terms are linearly independent over
> `C(n)`.

### Proof

Suppose a dependence exists, and choose one with the fewest terms.  Absorb
its nonzero rational coefficients into the terms to write

\[
 b_1+\cdots+b_m=0.
\tag{3.1}
\]

Put `s_i(n)=b_i(n+1)/b_i(n)`.  Shift (3.1), subtract `s_m` times (3.1), and
obtain

\[
 \sum_{i<m}(s_i-s_m)b_i=0.
\tag{3.2}
\]

Minimality forces every `s_i=s_m`.  Hence `b_i/b_m` is shift-invariant and
therefore constant.  Undoing the absorbed rational coefficients makes
`h_i/h_m` rational, contradicting dissimilarity.  \(\square\)

## 4. The shift-orbit classifier for factorial rays

For `lambda` nonzero define

\[
 h_{v,\lambda}(n)=\lambda^n\Phi_v(n).
\]

Its consecutive quotient is

\[
 \frac{h_{v,\lambda}(n+1)}{h_{v,\lambda}(n)}
 =\lambda\prod_{a\geq1}\prod_{u=1}^a(an+u)^{v_a}.
\tag{4.1}
\]

> **Lemma 4.1 (factorial-ray dissimilarity).**
> The terms `h_(v,lambda)` and `h_(w,mu)` are similar exactly when
> `v=w` and `lambda=mu`.

### Proof

Assume their quotient is `r(n)` in `C(n)` and set `delta=v-w`.  Then

\[
 \frac{r(n+1)}{r(n)}
 =\frac\lambda\mu
  \prod_{a\geq1}\prod_{u=1}^a(an+u)^{\delta_a}.
\tag{4.2}
\]

For a rational function `r`, the divisor of `r(n+1)/r(n)` has total
multiplicity zero in every translation orbit in `C/Z`: a zero or pole of
`r(n+1)` is the translate of one of `r(n)`.  The nonconstant factors on the
right of (4.2) have roots

\[
 -\frac ua\in[-1,0).
\tag{4.3}
\]

Let `A` be the largest index with `delta_A` nonzero.  In the translation
orbit of `-1/A`, all roots in (4.3) are in the same half-open interval, so
they can be translates only when they are equal.  The equality
`u/a=1/A` gives `a=uA`; maximality of `A` leaves only `a=A,u=1`.
The orbit multiplicity in (4.2) is therefore `delta_A`, contradicting the
zero orbit sum.  Descending induction gives `delta=0`.  Equation (4.2) now
says `r(n+1)/r(n)=lambda/mu`.  Every nonzero rational function has
`r(n+1)/r(n)` tending to one at infinity, so `lambda=mu`.  The converse is
immediate.  \(\square\)

### Proof of Theorem 2.1

Expand every coefficient in (2.3) and combine repeated pairs `(v,lambda)`.
The result is a rational-function dependence among terms
`h_(v,lambda)`.  Lemma 4.1 makes them pairwise dissimilar, and Lemma 3.1
makes every rational coefficient zero.  Grouping by `v` gives `c_j=0` for
every `j`.  \(\square\)

## 5. Gamma-affine rays and all explicit symmetries

### 5.1 Integer offsets

Consider an integer-affine factorial expression

\[
 H(n)=\alpha^n
 \prod_i(a_in+b_i)!^{\varepsilon_i},
 \qquad
 a_i>0,\quad b_i,\varepsilon_i\in\mathbb Z,
\tag{5.1}
\]

defined for all sufficiently large `n`.  Integer translation of a factorial
argument changes it by a rational function:

\[
 \frac{(an+b)!}{(an)!}\in\mathbb C(n)^\times.
\tag{5.2}
\]

After cancelling numerator and denominator factors, (5.1) has the canonical
signed slope vector

\[
 v_a=\sum_{i:a_i=a}\varepsilon_i
\tag{5.3}
\]

and can be written

\[
 H(n)=\alpha^nr(n)\Phi_v(n),
 \qquad r(n)\in\mathbb C(n)^\times.
\tag{5.4}
\]

Thus the complete equivalence relation over exponential-rational
coefficients is equality of (5.3).  It explicitly quotients:

1. permutation of factorial factors;
2. numerator--denominator cancellation;
3. integer shifts of affine arguments;
4. exponential rescaling; and
5. any geometric or marked-side redistribution which preserves the combined
   signed slope vector.

The fifth item is important for the Hall application.  The scalar factorial
functional cannot recover a marking which has already been forgotten.

### 5.2 Arbitrary offsets and Gauss refinement

Let `A` be a finite list of triples `(a_i,b_i,epsilon_i)` with
`a_i` a positive integer, `b_i` in `C`, and `epsilon_i` in `Z`, and put

\[
 G_A(n)=\alpha_A^n
 \prod_i\Gamma(a_in+b_i)^{\varepsilon_i}.
\tag{5.5}
\]

Positive slopes make this an eventually finite, nonzero sequence.  Its
**gamma-orbit signature** is the finite signed divisor

\[
 \Sigma_A
 =\sum_i\varepsilon_i\sum_{u=0}^{a_i-1}
  \left[\frac{b_i+u}{a_i}\right]
 \in\mathbb Z[\mathbb C/\mathbb Z].
\tag{5.6}
\]

> **Theorem 5.1 (complete gamma-affine classifier and separation).**
> For two positive-integer-slope gamma-affine rays `G_A,G_B`, the following
> are equivalent:
>
> 1. `Sigma_A=Sigma_B`;
> 2. there are \(C,\lambda\in\mathbb C^\times\) and
>    \(r\in\mathbb C(n)^\times\) such that
>    \[
>      \frac{G_A(n)}{G_B(n)}=C\lambda^nr(n)\qquad(n\gg0).
>    \tag{5.7}
>    \]
>
> Consequently rays with pairwise distinct gamma-orbit signatures have the
> separation property of Theorem 2.1 over `E_rat`.

### Proof

The consecutive quotient is

\[
 q_A(n)=\frac{G_A(n+1)}{G_A(n)}
 =\alpha_A\prod_i\prod_{u=0}^{a_i-1}
  (a_in+b_i+u)^{\varepsilon_i}.
\tag{5.8}
\]

Its finite zeros and poles are the negatives of the representatives in
(5.6).  If (5.7) holds, then

\[
 \frac{q_A(n)}{q_B(n)}
 =\lambda\frac{r(n+1)}{r(n)}.
\tag{5.9}
\]

The divisor of the rational coboundary on the right has total multiplicity
zero on every translation orbit, so `Sigma_A=Sigma_B`.

Conversely, [Gauss's multiplication formula](https://dlmf.nist.gov/5.5.E6)
gives the exact refinement

\[
 \Gamma(an+b)
 =(2\pi)^{(1-a)/2}a^{an+b-1/2}
  \prod_{u=0}^{a-1}\Gamma\!\left(n+\frac{b+u}{a}\right).
\tag{5.10}
\]

Thus (5.6) is precisely the signed multiset of unit-slope gamma factors after
all slopes are refined.  Equality modulo `Z`, followed by
`Gamma(z+1)=z Gamma(z)`, changes paired unit-slope factors only by rational
functions of `n`.  The constants and powers `a^(an)` in (5.10) supply `C`
and `lambda` in (5.7).  This proves the equivalence.

For the separation statement, expand every `E_rat` coefficient into
rational multiples of exponentials.  Terms belonging to distinct signatures
cannot be similar by (5.9); different exponential bases multiplying one fixed
ray cannot be similar unless the bases agree.  Lemma 3.1 then kills every
rational coefficient.  \(\square\)

For rational offsets, reduce each `(b_i+u)/a_i` to `[0,1)` and collect signed
multiplicities.  This is a directly computable canonical form.  The theorem
also gives a complete generator description of the symmetries in this class:

1. permutation and numerator--denominator cancellation;
2. Gauss multiplication or refinement;
3. integer translation of unit-slope offsets via the gamma recurrence; and
4. multiplication by constants, exponentials, and rational functions.

For example, duplication gives

\[
 \Gamma(2n+1)
 =\pi^{-1/2}4^n\Gamma(n+\tfrac12)\Gamma(n+1),
\tag{5.11}
\]

and both sides have signature `[0]+[1/2]`.  In contrast,
`Gamma(2n+1)/Gamma(n+1)^2` has nonzero quotient signature
`[1/2]-[0]`.  Its leading asymptotic `4^n/sqrt(pi*n)` is therefore not an
exact exponential-rational simplification.

### 5.3 Recurrence-facing form

The same classifier does not require a gamma presentation.  For
\(q(n)\in\mathbb C(n)^\times\), define

\[
 \sigma(q)
 =\sum_{[\rho]\in\mathbb C/\mathbb Z}
   \left(\sum_{\rho'\in[\rho]}\operatorname {ord}_{\rho'}q\right)[\rho].
\tag{5.12}
\]

A finite divisor on one translation orbit is of the form
`div(r(n+1)/r(n))` exactly when its orbit sum is zero: order the finitely many
points on that orbit and take cumulative sums.  Hence

\[
 \frac{q_i(n)}{q_j(n)}
 =c\frac{r(n+1)}{r(n)}
 \quad\Longleftrightarrow\quad
 \sigma(q_i/q_j)=0.
\tag{5.13}
\]

This is the divisor version of the multiplicative rational normal forms of
[Abramov and Petkovšek](https://doi.org/10.1006/jsco.2002.0522).  With
rational linear factors, (5.12) and the cumulative-sum certificate use exact
rational arithmetic.  Over a general coefficient field, polynomial
dispersion or a Gosper--Petkovšek normal-form implementation detects the
integer translates without explicitly adjoining all roots.  Thus a
recurrence search can group certified first-order factors by (5.12) before
connection constants or initial conditions are computed.

Lemma 3.1 therefore gives the same separation conclusion for arbitrary
hypergeometric terms whose shift quotients have distinct signatures (5.12),
even when no gamma presentation has been chosen.

An implementation must remove rational coboundaries in both numerator and
denominator.  A one-sided polynomial Gosper decomposition can leave a
nonconstant residual even when (5.13) holds with rational `r`.  The safe test
is the full multiplicative rational normal form, or the orbit-divisor
integration used by the checker followed by direct verification of (5.13).

### 5.4 `m`-fold hypergeometric terms

A nonzero sequence `H` is `m`-fold hypergeometric when

\[
 q_H^{[m]}(n)=\frac{H(n+m)}{H(n)}\in\mathbb C(n).
\tag{5.14}
\]

Replace translation by `1` in (5.12) with translation by `m`; that is, push
the divisor of \(q_H^{[m]}\) to \(\mathbb C/(m\mathbb Z)\).  The
cumulative-sum proof of (5.13) works verbatim on these larger orbits.

> **Corollary 5.2 (`m`-fold separation).**
> Pairwise distinct `m`-shift signatures separate `m`-fold hypergeometric
> terms over `E_rat`.  Equal signatures are the complete symmetry classes.

Indeed, restrict a relation to `n=mk+r`.  Each coefficient remains
exponential-rational in `k`, while the divisor on `C/(mZ)` becomes, after an
affine rescaling, the ordinary divisor on `C/Z`.  Lemma 3.1 applies on every
residue.  Conversely, equal signatures give

\[
 \frac{q_i^{[m]}(n)}{q_j^{[m]}(n)}
 =c\frac{r(n+m)}{r(n)}.
\tag{5.15}
\]

Choose \(\eta\) with \(\eta^m=c\).  Then
\(H_i(n)/(\eta^n r(n)H_j(n))\) is `m`-periodic, and every complex periodic
sequence is an exponential-polynomial by the finite Fourier transform on
`Z/mZ`.  It is therefore already an allowed coefficient symmetry.

This periodic factor cannot be omitted.  The sequences `H_1(n)=1` and
`H_2(2k)=1,H_2(2k+1)=2` have the same two-step quotient `1`, and

\[
H_2(n)=\frac{3-(-1)^n}{2}H_1(n).
\tag{5.16}
\]

Positive rational slopes are therefore included as well.  Choose `m` divisible
by every slope denominator.  Then shifting `n` by `m` increments every gamma
argument by a positive integer, so the ray is `m`-fold hypergeometric and
Corollary 5.2 gives its complete signature on `C/(mZ)`.  Equivalently, each
residue subsequence is a positive-integer-slope ray to which Theorem 5.1
applies.  For example, duplication classifies `Gamma(n+1)` together with
`Gamma(n/2+1/2)Gamma(n/2+1)` at step two, with exponential base `4` per
two-step shift.

Negative slopes which create infinitely many poles still require a different
choice of eventual domain and are not included here.

## 6. Stirling transseries and inverse moments

The exact difference-algebra proof has an asymptotic audit.  Stirling's
series, in the normalization recorded in
[DLMF 5.11](https://dlmf.nist.gov/5.11), gives

\[
\begin{aligned}
 \log\Phi_v(n)
 ={}&S(v)n\log n+\bigl(E(v)-S(v)\bigr)n
      +\frac{L(v)}2\log n+C(v)\\
 &+\sum_{m\geq1}
   \frac{B_{2m}}{2m(2m-1)n^{2m-1}}
   M_{2m-1}(v),
\end{aligned}
\tag{6.1}
\]

where

\[
\begin{aligned}
 S(v)&=\sum_aav_a,&
 E(v)&=\sum_aav_a\log a,\\
 L(v)&=\sum_av_a,&
 C(v)&=\frac12\sum_av_a\log(2\pi a),\\
 M_{2m-1}(v)&=\sum_a\frac{v_a}{a^{2m-1}}.
\end{aligned}
\tag{6.2}
\]

For a positive partition, `E(v)` is the proposed entropy term
`sum_i a_i log(a_i)`.  It is not complete.  The partitions

\[
 (12,6,4,4,4,1),\qquad(9,8,8,2,2,2)
\tag{6.3}
\]

have the same total `31`, length `6`, product `4608`, and entropy product
`product_i a_i^(a_i)`, so every term of (6.1) through `C(v)` agrees.  Their
first inverse moments are respectively

\[
 2,\qquad\frac{67}{36}.
\]

The negative odd moments are complete in the zero-shift class.  If
`delta` is nonzero and `a_0` is its smallest active slope, then

\[
 M_{2m-1}(\delta)
 =a_0^{-(2m-1)}\left(\delta_{a_0}+o(1)\right)
 \quad(m\longrightarrow\infty).
\tag{6.4}
\]

Hence equality of all the inverse moments forces `delta=0`.  The Stirling
hierarchy therefore gives the right diagnostic.  Difference algebra is still
needed for the exact theorem because an allowed exponential absorbs an
entropy difference and a rational prefactor has its own infinite inverse
power expansion.

## 7. Consequences and limits

### 7.1 Fixed factorial packets

Suppose a fixed finite packet has an exact all-scale expression

\[
 \sum_x P_x(N)\Lambda_x^N\Phi_{v_x}(N)=0,
 \qquad P_x\in\mathbb C(N).
\tag{7.1}
\]

Theorem 2.1 splits it into the identities

\[
 \sum_{x:v_x=v}P_x(N)\Lambda_x^N=0
\tag{7.2}
\]

for every canonical vector `v`.  Thus distinct factorial rays cannot sustain
an all-order cancellation, even in an arbitrary finite mixed packet.  Only
states with the same canonical factorial vector and a cancelling
exponential-rational trace remain.

This does not prove unrestricted `GVC(2)`.  The open Hall--jet step is still
to expose a packet whose states, markings, and congruence classes are fixed
and scale-compatible.  Moreover, the all-span collision family

\[
 R_{s+6}B_aB_{a+1}=R_sB_{a+3}B_{a+4}
\]

in [`BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md`](BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md)
has the same canonical factorial vector on both sides.  It is an explicit
infinite projected family left untouched by Theorem 2.1; the strengthened
span-seven census finds no other fully decorated `C2,C3` survivor, and its
first separator on the primitive pair is the marked `C4` character.  The
complete blind scaled fibre is nevertheless terminal: its Franel rows at the
first two periods of any fixed finite character force support loss.  Hence the
family survives factorial-trace independence only before fixed-packet
exposure; it supplies no nonzero scale-compatible packet afterward.

### 7.2 SIC and holonomic recurrences

If a certified recurrence solution is known to be a finite sum of
gamma-hypergeometric rays, Theorem 5.1 makes the decomposition into
gamma-orbit signature classes unique.  Petkovšek's `Hyper` algorithm, or
first-order factorization in the shift Ore ring, can certify all such
components.  Formula (5.12) then groups the factors without choosing a
possibly noncanonical gamma presentation.  Only factors inside one signature
class need connection-constant or initial-value cancellation tests.

There is an immediate exact application to the SIC2C4 moment family in
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md):

\[
 M_{d,r}(m)
 =4^{rm}\frac{\Gamma(dm+3)\Gamma(rm+1)^2}{\Gamma(2rm+2)},
 \qquad d=4r+k,\quad k\geq0.
\]

Its signature is `U_(d,3)+2U_(r,1)-U_(2r,2)`, where

\[
 U_{a,b}=\sum_{u=0}^{a-1}\left[\frac{b+u}{a}\right].
\]

> **Corollary 7.1 (SIC radial-family separation).**
> Distinct pairs `(d,r)` with `d>=4r` give distinct gamma-orbit signatures.
> Hence no finite exponential-rational identity can mix distinct SIC radial
> moment families.

To prove this, compare two signatures.  If their `d` values differ, the
primitive residue `1/D` for the larger value `D` occurs in its `U_(D,3)` and
cannot occur in any smaller slope.  Thus the `d` values agree and those terms
cancel.  If the `r` values differ, the primitive residue `1/(2R)` for the
larger `R` occurs with coefficient `-1` in `U_(2R,2)` and in none of the
remaining smaller slopes.  Thus the `r` values also agree.  The checker
replays all 276 pairs `(d,r)` with `d<=48`; the argument is unbounded.

A proper-hypergeometric **multisum** is generally only P-recursive.  It need
not be hypergeometric or a finite sum of hypergeometric terms.  Consequently
Theorems 2.1 and 5.1 are exact prefilters for the SIC recurrence and holonomic
counterexample searches, not a completeness theorem for them.  Formal
Birkhoff--Trjitzinsky branches likewise require an actual connection theorem
and nonzero Stokes data before they become an exact decomposition of a
specific sequence.

### 7.3 Positive characteristic

The proof uses characteristic zero essentially.  There is already an exact
collapse at the level of raw valuations.  Put

\[
 F_a^{(p)}(n)=\nu_p((an)!).
\]

Legendre's formula gives

\[
 F_{pa}^{(p)}(n)
 =\sum_{j\geq1}\left\lfloor\frac{pan}{p^j}\right\rfloor
 =an+F_a^{(p)}(n),
\]

and hence, by iteration,

\[
 F_{p^ku}^{(p)}(n)
 =F_u^{(p)}(n)+u\frac{p^k-1}{p-1}n.
\tag{7.3}
\]

Thus all slopes in one `p`-power orbit become equal modulo a linear sequence.
For any positive `a,b`, there is the exact four-ray relation

\[
 b\bigl(F_{pa}^{(p)}-F_a^{(p)}\bigr)
 -a\bigl(F_{pb}^{(p)}-F_b^{(p)}\bigr)=0.
\tag{7.4}
\]

For example, at `p=2`,
`3F_2-3F_1-F_6+F_3=0` for every scale.  More generally a signed slope vector
reduces canonically to

\[
 \sum_{u:p\nmid u}
   \left(\sum_{k\geq0}v_{p^ku}\right)F_u^{(p)}(n)
 +\frac{n}{p-1}
  \sum_{u:p\nmid u}\sum_{k\geq0}u(p^k-1)v_{p^ku}.
\tag{7.5}
\]

This is the correct first quotient for a characteristic-`p` valuation phase:
aggregate each primitive slope tower and retain its linear drift.  It is not
complete.  Raw factorial rays eventually vanish modulo `p`, the
shift-constant field changes, exponential bases can coalesce, and normalized
`p`-free units depend on base-`p` digits and carries.  Theorem 2.1 may still
separate characteristic-zero components before integral reduction, but the
primitive recurrence factor ledger and, in higher order, the integral
companion-lattice/Smith data remain necessary.

### 7.4 Practical search integration

The theorem is most useful as an exact early partition, not as a late
asymptotic heuristic:

1. **Certify the finite packet.**  In a Hall search this means one packet with
   fixed markings and scale dependence; in a recurrence search it means
   certified hypergeometric components or first-order shift factors.
2. **Compute the signature.**  Use (5.6) for gamma data, (5.12) for an
   ordinary rational shift quotient, and the `C/(mZ)` version for an
   `m`-step quotient.
3. **Reject cross-signature cancellation.**  No connection constant, matrix
   trace, or exponential-polynomial coefficient can cancel two different
   classes.
4. **Normalize equal classes.**  Integrate the zero-sum orbit divisor to
   obtain `r`, remove the exponential base, and combine all terms into one
   rational/exponential trace coefficient.
5. **Spend expensive computation only on the residual class.**  Test its
   markings, initial values, Stokes/connection data, or finite-character
   traces.  If the object is only P-recursive and has no certified finite
   hypergeometric decomposition, stop: factorial-trace independence is not a
   completeness argument.
6. **Before reduction modulo `p`, retain the characteristic-zero classes.**
   After taking valuations, first apply the primitive-slope reduction (7.5),
   then compute the unit, digit, carry, and integral-lattice data.

## 8. Exact regression

Run

```bash
python3 scripts/verify_factorial_trace_independence.py
```

For a direct exact comparison, encode each factor as
`slope:offset[:multiplicity]`:

```bash
python3 scripts/verify_factorial_trace_independence.py \
  --compare '2:1' '1:1/2,1:1'
python3 scripts/verify_factorial_trace_independence.py \
  --compare '2:1' '1:1:2'
```

The first command returns the duplication symmetry, exponential base `4`,
and trivial rational coboundary.  The second returns the nonzero signature
`[1/2]-[0]`.  Equal signatures with integer-shifted offsets also print the
factored rational coboundary in (5.13).

An optional independent CAS replay is

```bash
make verify-factorial-trace-independence-sympy
```

It asks SymPy to simplify all 322 configured Gauss/shift certificates and all
1,384 census collision certificates as rational-function identities.  The
default checker remains dependency-free.

The dependency-free checker:

1. reconstructs every signed vector on slopes `1,...,7` with entries in
   `[-2,2]` from its exact translation-orbit divisor;
2. integrates all 8,134 nonzero zero-sum divisors with coefficients in
   `[-2,2]` on a seven-point rational translation orbit;
3. classifies 67,524 products of at most three rational-offset gamma atoms
   into 66,140 signatures, certifying all 1,384 symmetry collisions;
4. verifies 161 Gauss refinements before and after nontrivial integer shifts,
   plus 1,000 seeded signed transformation cases;
5. verifies representative positive and negative integer-affine shift
   reductions to the canonical vector;
6. replays a one-scale factorial collision and the entropy collision (6.3);
7. checks the periodic exponential-polynomial symmetry, a two-step
   coboundary certificate, and rational-slope duplication for `m`-fold terms;
   and
8. verifies the `p`-primitive valuation reduction (7.5) for 2,187 signed
   profiles at `p=2,3,5`, together with exact four-ray kernels; and
9. separates all 276 SIC radial-moment families `(d,r)` with `d<=48`.

These finite checks are regressions for the exact formulas.  The all-scale
separation results are Theorems 2.1 and 5.1 and Corollary 5.2; the valuation
collapse follows directly from Legendre's formula.
