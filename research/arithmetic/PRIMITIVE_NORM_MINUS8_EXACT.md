# Exact sums of two primitive discriminant -8 norms

Let
\[
 Q(x,y)=x^2+2y^2,\qquad
 P=\{Q(x,y)>0:\gcd(x,y)=1\},
\]
and let \(C\) contain 1 and the positive integers whose prime factors
are all \(1\) or \(3\pmod8\). Zero coordinates are allowed; summands
are positive.

**Theorem M8.** The following exact identities hold:
\[
 \begin{split}
 C+C&=\{N\ge2:N\text{ even},\ 8\nmid N\},\\
 C+2C&=\{N\ge3:N\text{ odd}\},\\
 P+P&=\{N\ge2:16\nmid N\}.
 \end{split}                                               \tag{M8.1}
\]
There are no exceptional admissible integers. For every \(N\ge2^{78}\)
with \(16\nmid N\), the existing explicit count also gives
\[
 \#\{(a,b)\in P^2:a+b=N\}\ge N^{2/3}/4800.                  \tag{M8.2}
\]

The proof combines the previously established quaternionic Hecke bounds,
sharper lattice/divisor estimates, exact rational interval certificates
and complete finite coverage through \(2^{27}\). The interval bounds and
descent make the finite computation complete for the unrestricted claim.
This is an internally reviewed deduction from classical inputs, without
formal verification, external review or an established priority claim.

## The primitive criterion and the two strong cores

The exact primitive value criterion is
\[
                              P=C\cup2C.                  \tag{M8.3}
\]
Here is its arithmetic justification. The ring \(\mathbb Z[\sqrt{-2}]\)
is norm Euclidean: rounding both coordinates leaves error norm at most
\(1/4+2/4=3/4<1\). An inert prime \(5,7\pmod8\) dividing a norm divides
both coordinates. The ramified prime of norm 2 may occur to exponent at
most one in a primitive value. At a split prime, selecting its entire
power on one of the two principal prime ideals produces a primitive
generator. These observations give both directions of (M8.3); the full
argument is also E3 of the [exact arithmetic source](PRIMITIVE_NORM_EXACT_CLASSIFICATIONS.md).
In particular every element of \(C\) is \(1\) or \(3\pmod8\).

It follows that
\[
                  P+P=(C+C)\cup2(C+C)\cup(C+2C).            \tag{M8.4}
\]
The notation \(2(C+C)\) denotes dilation. The first set has only
residues \(2,4,6\pmod8\); the second has only \(4,8,12\pmod{16}\);
the last consists of odd integers at least 3. Thus multiples of 16 are
excluded, and the first two identities in (M8.1) imply the last one.

For clarity, the two strong core identities are also necessary for the
full primitive identity. A target of valuation one at 2 must be in
\(C+C\). If a target \(N\) has valuation two, apply the primitive
identity to \(2N\): its valuation is three, so its summands must both
belong to \(2C\); dividing them by 2 represents \(N\) in \(C+C\).
Odd targets are exactly the mixed core. No unrestricted closure under
doubling a primitive coordinate representation is assumed.

Multiplication by \(\varpi=\sqrt{-2}\) acts on a coordinate pair as
\[
 (x,y)\longmapsto(-2y,x),\qquad Q(\varpi v)=2Q(v).
\]
It bijects the primitive odd-norm vectors with the primitive even-norm
vectors; the latter have even first coordinate and odd second coordinate.
The inverse is \((X,Y)\mapsto(Y,-X/2)\). At odd primes its determinant
is a unit. This verifies every dyadic division used below directly.

## Previously proved shell and Hecke inputs

We use precisely the following results from the
[analytic bridge](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md), equations (29)--(40a).
Write \(N=2^e m\), \(m\) odd, \(r(N)\) for the full signed ordered
\(Q\oplus Q\) shell count, and \(G_2(N)\) for its subset with neither
binary coordinate pair divisible by 2. Then
\[
 r(2^e m)=4c_e\sigma(m),\qquad
                   (c_0,c_1,c_2,c_3)=(1,2,6,6),
\]
\[
 (G_2(m),G_2(2m),G_2(4m),G_2(8m))
       =(2\sigma(m)-2a_m,4\sigma(m),20\sigma(m),16\sigma(m)),
                                                               \tag{M8.5}
\]
where \(a_m\) is the coefficient of
\(\eta(4\tau)^2\eta(8\tau)^2\) and
\(|a_m|\le\tau(m)\sqrt m\). Let
\[
 \mathcal I=\{p\text{ prime}:p\equiv5,7\pmod8\},\qquad
 E(m)=\frac{\tau(m)\sqrt m}{\sigma(m)},\qquad
 M(X)=\max_{1\le t\le X}
                   \prod_{\substack{q\mid t\\q\equiv1,3(8)}}(v_q(t)+1).
\]
For every inert prime \(p\), let \(b_p(N)\) count full-shell points
with at least one field component divisible by \(p\). The all-ideal
Hecke estimate gives, for all target valuations,
\[
 b_p(N)/r(N)\le\mu_p+E(m),\qquad \mu_p=2/[p(p-1)].          \tag{M8.6}
\]
For \(p\nmid N\), this is normalizer averaging on the finite
projective Hecke quotient. For \(p\mid N\), anisotropy forces both
components divisible by \(p\), and the exact ratio is zero or
\(\sigma(m/p^2)/\sigma(m)<p^{-2}\le\mu_p\). Thus (M8.6) does not
apply an operator at its own level prime.

The underlying order is the index-four Hurwitz suborder of norm
\(Q\oplus Q\); its locally principal right-ideal class number is one.
The cited proof checks this by unit transitivity on six local cosets,
normalizes by its four units, transfers the noncharacter spectrum by
Jacquet--Langlands, and bounds the all-ideal Hecke roots by the weight-two
Ramanujan theorem. Its fixed dyadic operators have degrees \(c_e\);
no unramified estimate at 2 is used. The source ledger pins that proof
and the exact primary eta/Sturm imports used for (M8.5).

Two further already-proved conclusions are needed:

* Every prime at least 2 belongs to \(P+P\). The proof has an analytic
  prime threshold \(2^{22}\) and an exhaustive certificate for all
  295947 primes through that threshold.
* For every odd inert core \(m\ge2^{60}\) and \(0\le e\le3\),
  \(2^e m\in P+P\).

Both have complete proofs and retained certificates in the cited source;
they do not depend on the exact identity now being proved. The second
implies all three large-core conclusions needed here: at \(m\) the
representation is mixed; at \(2m\) it is in \(C+C\); and at \(8m\)
both values lie in \(2C\), so division by 2 gives \(4m\in C+C\).
The quantitative bound (M8.2) is the explicit surviving-shell/multiplicity
calculation in that same source, with its numerical onset unchanged.

## Sharper point counts and the selected shells

For every real \(X\ge0\),
\[
 S_Q(X):=\#\{(x,y):0<Q(x,y)\le X\}\le\frac83 X,\qquad
 S_o(X):=\#\{(x,y):Q(x,y)\le X,\ Q(x,y)\text{ odd}\}\le2X.
                                                               \tag{M8.7}
\]
Here are complete elementary bounds. Row integration gives
\[
 S_Q(X)\le\frac\pi{\sqrt2}X+(2+\sqrt2)\sqrt X.
\]
For integer \(X\ge64\), use \(\pi<22/7\), \(\sqrt2>24/17\), and
\(\sqrt2<10/7\). The ratio is less than
\(187/84+3/7=223/84<8/3\). Exact lattice enumeration checks integers
1 through 63. An odd norm has odd first coordinate; counting these
points row by row gives
\[
 S_o(X)\le\frac\pi{2\sqrt2}X+(1+\sqrt2)\sqrt X+1.
\]
For integer \(X\ge16\) the same rational estimates give a ratio less
than \(187/168+17/28+1/16<2\); integers 1 through 15 are checked
exactly. Applying the integer estimates to \(\lfloor X\rfloor\)
extends both statements to all real \(X\). At \(X<1\) both counts vanish.

We sieve the following three subsets of full shells, for odd \(m\):
\[
 \begin{array}{c|c|c|c}
 N&\text{selected subset }S&|S|/r(N)\text{ at least}&L\\\hline
 m&G_2(m)&1/2-E(m)/2&8/3\\
 2m&\text{both field norms odd}&1/2&2\\
 4m&\text{both field norms odd}&2/3&4/3.
 \end{array}                                               \tag{M8.8}
\]
The even subset counts are exact. Multiplying both field components by
\(\varpi\) bijects the full shell of norm \(N/2\) with the part of
the norm-\(N\) shell having both field norms even. Hence its odd-odd
complement has size \(r(N)-r(N/2)\): these sizes are respectively
\(4\sigma(m)\) and \(16\sigma(m)\). The same count at \(4m\) is
\(G_2(4m)-G_2(2m)=16\sigma(m)\).

All selected components are nonzero. Unique ideal factorization gives
\(r_Q(t)\le2M(N)\) for \(1\le t\le N\). For the first row, count a
\(p\)-divisible nonzero component using \(S_Q(N/p^2)\), the other
using \(2M(N)\), and both orientations. Its bad count is at most
\((32/3)M(N)N/p^2\); divided by \(r(m)\ge4m\) this gives \(L=8/3\).
For the even rows use \(S_o(N/p^2)\le2N/p^2\) instead. The numerator
is \(8M(N)N/p^2\). Dividing by \(r(2m)\ge8m\) or
\(r(4m)\ge24m\) gives the other two constants.

Put
\[
 B=\sum_{p\in\mathcal I}\mu_p<3/16,\quad
 \Pi(P)=\#\{p\in\mathcal I:p\le P\},\quad
 T(P)=\sum_{\substack{p\in\mathcal I\\p>P}}p^{-2}.
\]
The stated bound for \(B\) is the existing complete finite-plus-tail
bound from the analytic bridge. With \((\lambda,L,\delta)\) equal to
\((1/2,8/3,1/2)\), \((1/2,2,0)\), or \((2/3,4/3,0)\), a sufficient
criterion for the corresponding row is
\[
 B+E(m)(\Pi(P)+\delta)+LM(N)T(P)<\lambda.                  \tag{M8.9}
\]
Indeed small bad-prime counts may be bounded in the full shell using
(M8.6); large bad-prime counts use the selected subset and (M8.7).
Subtract their union from \(|S|\). The selected subset need not itself
be stable under the full-shell Hecke operator. On survival no inert
prime divides either norm, and the dyadic conditions put the values in
\(C+2C\) or \(C+C\), respectively, by (M8.3).

## Exact interval bounds for inert cores

If every prime factor of the odd integer \(m\) is inert, then
\[
 E(m)m^{1/3}
       =\prod_{p^a\parallel m}\frac{(a+1)p^{5a/6}}{\sigma(p^a)}
                             <18/5.                       \tag{M8.10}
\]
The finite certificate makes this an unrestricted inequality. Raise
each factor to the sixth power. For \(p\ge64\), it is at most
\((a+1)^6/p^a\le1\). For inert \(p<64\), only exponents 0 through 3
need be checked: for \(a\ge3\) the next/current ratio is strictly less
than \((5/4)^6/p\le(5/4)^6/5<1\). The exact product of these maxima
is below \((18/5)^6\). Every prime and exponent has been accounted for.

The divisor envelope \(M(X)\) is exactly computable in the bounded
range needed below. Remove nonsplit factors, sort the remaining
exponents in decreasing order, and assign them to the smallest split
primes. This preserves the divisor product and decreases the integer:
exchanging \(a<b\) at \(p<q\) multiplies it by \((p/q)^{b-a}<1\).
Thus initial split-prime prefixes with nonincreasing exponents suffice.
The checker enumerates all such profiles through \(2^{63}\), with an
explicit prime-prefix and exponent-cap completeness check, and retains
the resulting exact maximum frontier.

Let \(R=65536\) and \(D=2^{64}\). For each integer \(P<R\), exactly
sieve the inert primes through \(R\) and use the rational upper bound
\[
 V(P)=\frac1D\sum_{\substack{p\in\mathcal I\\P<p\le R}}
           \left\lceil\frac D{p^2}\right\rceil
                            +\frac1{4R}+\frac2{R^2}\ge T(P).
                                                               \tag{M8.11}
\]
The last two terms majorize all remaining primes by the two progressions
\(5,7\pmod8\), integrating \(x^{-2}\) separately in each progression.
The prime count \(\Pi(P)\) is exact within this finite sieve.

For each row of (M8.8), put \(N=fm\), where \(f=1,2,4\). In a
complete dyadic interval \(2^k\le m<2^{k+1}\), set
\[
 d_{f,k}=M(f(2^{k+1}-1)),\qquad
 e_k=\frac{18}{5\lfloor\sqrt[3]{2^k}\rfloor}.
\]
Choose the integer threshold recorded in the receipt (each is less than
\(R\)). The complete exact rational certificate verifies
\[
 \frac3{16}+e_k(\Pi(P)+\delta)+L d_{f,k}V(P)
                            <\lambda-\frac1{256}           \tag{M8.12}
\]
for every interval in the following table:

| Target family | Every interval index checked | Number of intervals |
|---|---|---:|
| Odd \(m\), mixed core | \(27\le k\le59\) | 33 |
| \(2m\), even strong core | \(26\le k\le59\) | 34 |
| \(4m\), even strong core | \(20\le k\le59\) | 40 |

All 107 inequalities and their exact inputs are retained. The existing
\(m\ge2^{60}\) theorem supplies the strong-core conclusions above the
last interval, as explained earlier. Thus the only finite work left is
mixed targets below \(2^{27}\), \(2m\) targets below \(2^{27}\), and
\(4m\) targets below \(2^{22}\).

## Complete finite coverage and descent

The checker covers a larger finite range than required: both strong
identities through \(B_0=2^{27}\). It sieves every forbidden prime to
construct the exact flags of \(C\cap[1,B_0]\), with no prime beyond
\(B_0\) able to divide a value in that interval. Independent primitive
coordinate enumeration agrees with these flags through 10000.

Pack the flags into an integer \(B_C\), with bit \(j\) indicating
\(j\in C\); likewise form \(B_{2C}\). Exact bit shifts and unions give
\[
 U_e=\bigvee_{\substack{a\le5000\\a\in C}}(B_C\ll a),
\]
\[
 U_o=\bigvee_{\substack{a\le5000\\a\in C}}(B_{2C}\ll a)
       \ \vee\!
       \bigvee_{\substack{b\le5000\\b\in2C}}(B_C\ll b).     \tag{M8.13}
\]
After truncation at \(B_0\), the first union equals exactly the periodic
mask with residues \(2,4,6\pmod8\), and the second equals exactly the
odd-target mask with 1 removed. There are 699 short \(C\)-summands
and 365 short \(2C\)-summands. This certifies all
\[
 50331648\text{ even strong targets},\qquad
 67108863\text{ mixed odd targets},
\]
with zero omissions and zero forbidden positions. No assumption that a
short summand must exist enters the proof: the computed union actually
covers every required bit. The receipt retains the exact ranges,
counts, membership and coverage hashes. The algorithm reconstructs all
bits; a hash alone is not the mathematical certificate.

Suppose first that a least eligible exception to \(C+C\) exists.
It has valuation one or two at 2. Any split prime factor \(q\) can be
removed: the quotient remains an eligible even integer at least 2,
and a representation of that smaller target lifts by multiplication
by \(q\). Thus the least exception has the form \(2m\) or \(4m\)
with odd inert core \(m\). The finite certificate, its corresponding
interval family, and the large-core theorem cover every such \(m\),
a contradiction.

For the mixed identity, a least missing odd target \(N\ge3\) either
has no split prime factor, or is itself a single split prime. A composite
target with a split prime factor has quotient at least 3, so the same
multiplicative descent applies. The prime boundary is supplied by the
already-proved all-prime theorem: on an odd target, every \(P+P\)
representation has one value in \(C\) and one in \(2C\). An inert-core
exception is excluded by the finite certificate, the odd intervals and
the large-core theorem. This proves both strong identities. Their
congruence exclusions were proved before (M8.4), which now gives the
full primitive identity (M8.1). \(\square\)

## Evidence and scope

Run the [checker](verify_primitive_norm_minus8.py), standard library only,
from the repository root:

```sh
python3 research/arithmetic/verify_primitive_norm_minus8.py \
  --output research/arithmetic/primitive_norm_minus8_checks.json
```

The [receipt](primitive_norm_minus8_checks.json) retains all interval
and finite-coverage data. The [source ledger](primitive_norm_minus8_source_ledger.json)
pins the prior shell/Hecke/prime-family proofs, their primary imports,
and the old finite prime certificate; its current replay matches the
retained record. The new checker certifies finite arithmetic, complete
Euler/profile searches and the stated target masks.

Certificate verification requires active Python assertions. The checker rejects optimized execution (-O, -OO, or optimization enabled through PYTHONOPTIMIZE) before generating a receipt.

The separate [coordinate membership verifier](verify_primitive_norm_coordinate_membership.py)
authenticates every membership flag through the actual bound, rather than
only the coordinate regression prefix. Independent sign changes allow
\(x\ge1\) odd and \(y\ge0\); the exact inequality
\(x^2+2y^2\le2^{27}\), with \(\gcd(x,y)=1\), enumerates all of \(C\)
in that range. Its complete bytearray must equal the original prime-support
sieve's bytearray before either is hashed. The
[membership receipt](primitive_norm_coordinate_membership_checks.json)
also checks the companion \(-11\) ellipse and rejects inserted inert
prime powers and a deleted unit. Run both domains with

```sh
python3 research/arithmetic/verify_primitive_norm_coordinate_membership.py \
  --output research/arithmetic/primitive_norm_coordinate_membership_checks.json
```

Target-mask equality alone does not authenticate summand membership;
the full coordinate comparison supplies this separate finite check. The infinite
spectral theorem and the completion arguments remain written
mathematics, rather than consequences of the checker merely passing.
Several internal agents reconstructed the local, spectral, geometric
and descent steps separately. External review, formal verification
and historical priority have not been established. The independent
\(C_{12}+C_{12}\) barrier and general unweighted leading constants
remain open.
