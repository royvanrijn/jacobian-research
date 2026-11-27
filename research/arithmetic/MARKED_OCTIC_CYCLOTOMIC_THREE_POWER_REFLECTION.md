# Higher three-power quotients on the marked octic's cyclotomic mirror

The fourth marked norm can give a cyclic quotient of
order \(3^n\) in the class group of a cyclotomic layer
of the imaginary mirror. The exponent at eleven in
the parameter scale makes the same norm a local
\(3^n\)-th power at three. The other three markings
supply independent cubic quotients. On actual prime
ratios, a fixed rational residue test controls the
numerator prime, and the denominator prime splits
completely.

## The higher power scale

**Proposition YR1.** Fix \(m=3^n\) for \(n\ge1\), choose
\(0\le a<m\) with \(4a\equiv3\bmod m\), and put
\[
 c_m=\frac{3}{7}11^a,\qquad t=c_mu^m,\qquad
 \gamma=\frac{\beta^2}{11^3},\qquad
 \delta_i=\frac{\alpha_i}{\bar\alpha_i}.
\]
Here \(F,\beta,\alpha_i\) are the data of
[MO1–MO3](QUADRATIC_MIRROR_CLASS_RANK_FOUR_FROM_MARKED_OCTIC.md).
There is a nonempty fixed local slice for \(u\) on
which \(K=\mathbb Q(\sqrt{F(t)})\) is real quadratic,
\(\gamma\) is an \(m\)-virtual element and a local
\(m\)-th power above three, and all \(\delta_i\) are
virtual cubic elements and cubes above three.
The four elements are independent modulo cubes.
The primes \(2,3,5,7,11\) all split in \(K\).

Use the trace \(T=49t^2/9+125/3+196/t^2\)
of [BS1](MARKED_OCTIC_BETA_CUBIC_SCALE_CLASSIFICATION.md).
Choose sufficiently large positive \(N\). At
\(2,3,7,11\) prescribe \(v_p(u)=-N\), and at
\(5\) prescribe \(v_5(u)=N\). At these centers
\(F(t)\) is a square locally, since \(t\) is close
to zero or infinity and the corresponding leading
coefficient is a nonzero square. The negative
trace valuations at \(2,5,7,11,3\) are
\[
       -2mN,\quad-2mN,\quad-2mN,\quad
                   2(a-mN),\quad-2mN ,
\]
respectively. At \(p\ne11\) the two valuations of
\(\gamma\) are plus and minus twice the trace
valuation. At eleven they are
\(4(a-mN)-3\) and its negative. These are all
multiples of \(m\). Outside these primes,
when \(t\) is a unit both values are units.
When it is a nonunit the trace valuation is
\(-2|v_p(t)|\) and \(v_p(t)\equiv0\bmod m\).

At the three-adic center let \(w=3^{-v_3(t)}t\).
For the root \(\beta\) of smaller valuation,
\(\beta/T\) tends to one as \(N\) grows, since the
constant term \(11^3\) is a unit. Removing the
\(3\)-power from \(\gamma\), its leading unit is
\[
         (7w)^4/11^3=11^{\,4a-3} v^{\,4m},
         \qquad v=3^{-v_3(u)}u\in\mathbb Z_3^\times .
\]
This is an \(m\)-th power. The error tends to zero
with \(N\) and lies in \(1+3^{n+1}\mathbb Z_3\)
for large enough \(N\), a subgroup of the \(m\)-th
powers. The other root gives an inverse because
\(\gamma\) has norm one.

For the cubic marks, the leading terms of
\(\alpha_1,\alpha_2,\alpha_3\) at the four boundary
points have respective orders \(-3,-3,0\), with
the following cube roots of their leading
coefficients. At infinity the parameter is \(1/t\).
\[
\begin{array}{c|rrr}
 &\alpha_1&\alpha_2&\alpha_3\\ \hline
 0^+&-7/11&21/11&15/11\\
 0^-&7/11&-21/11&-5/11\\
 \infty^+&7/66&7/22&-5/11\\
 \infty^-&-7/66&-7/22&15/11
\end{array}
\]
These follow by substituting the boundary series
into the elliptic maps of MO1. Thus at all five
fixed primes the \(\alpha_i\) specialize to cubes
on both branches when \(N\) is sufficiently large.
At the other primes the good-prime argument of
MO2 gives virtuality. Passing to the quotient
by the conjugate preserves those properties and
doubles each cubic character.

Here are four independent residue rows, written for
\(\alpha_1,\alpha_2,\alpha_3,\beta\). Set
\(t=3u_0^3/7\bmod q\). The \(R\) column is the chosen
square root of \(F(t)\).
\[
\begin{array}{c|rrr|c}
 q&u_0&R&\zeta&\text{exponents}\\ \hline
13&7&1&3&(0,2,2,1)\\
31&1&5&25&(0,1,2,2)\\
43&10&14&36&(1,1,0,0)\\
61&4&17&47&(0,2,1,0)
\end{array}
\]
At each of these primes \(v_3(q-1)=1\), and
\(a\equiv0\bmod3\). Thus every listed \(t/c_m\)
is in the image of the \(m\)-th power map, for
every \(n\). For \(\delta_i,\gamma\) the rows are
doubled, still giving rank four. Weak approximation
makes the tests simultaneous. The squarefree
part of the degree-\(8m\) form provides quadratic
specializations. Positivity of \(F\) gives the real
sign. The [checker](verify_marked_octic_cyclotomic_three_power_reflection.py)
recomputes the boundary coefficients, character
rows and higher-power local data at \(m=3,9,27\).
The valuation and leading-unit arguments apply at
every exponent.

## Descent to the cyclotomic layer

**Theorem YR2.** For a quadratic field in YR1 with
fundamental discriminant \(D\), put
\[
 M=\mathbb Q(\sqrt{-3D}),\qquad
 B=M\mathbb Q(\zeta_m)^+ .
\]
There is an explicit abelian extension \(W/B\),
unramified at every finite prime and split above
three, with
\[
 \operatorname{Gal}(W/B)\simeq
           \mathbb Z/m\times(\mathbb Z/3)^3 .
 \tag{1}
\]
Over \(M\) the extension is Galois and
\[
 \operatorname{Gal}(W/M)\simeq
  \bigl((\mathbb Z/m)\rtimes_4(\mathbb Z/(m/3))\bigr)
                     \times(\mathbb Z/3)^3 ,
 \tag{2}
\]
where a generator acts by multiplication by four
on the order-\(m\) factor. Moreover
\[
               |\operatorname{Disc}(B)|
                       =3^{(2n-1)m/3}D^{m/3}.
 \tag{3}
\]

Put \(L=K\mathbb Q(\zeta_m)\).
The quadratic field \(K\) is real, whereas the only
quadratic subfield of the cyclotomic field is
imaginary. Thus \(L/K\) is abelian and disjoint
from the cubic extension given by the root of
any noncubic element of \(K\). Indeed a cubic
extension obtained from such an element would
require \(\zeta_3\) in order to be Galois, and has
degree three over the real field \(K\).
Consequently the four independent cube classes of
YR1 remain independent modulo cubes in \(L\).

Adjoin the \(m\)-th root of \(\gamma\) and third
roots of the three \(\delta_i\), giving \(L'\).
The classes \(\gamma,\delta_i^{m/3}\) generate
the direct product of one order-\(m\) and three
order-three factors in \(L^\times/L^{\times m}\).
In a putative relation divide the exponents by
their greatest common power of three, say \(3^r\).
Here \(r<n\), and the possible factor in
\(\mu_{3^r}\) is itself a cube in \(L\).
A nonzero cube relation contradicts YR1.
Kummer duality gives the stated group over \(L\).
The valuations and local powers give an unramified
extension \(L'/L\) split above three.

The involution \(s\) acting simultaneously by
quadratic conjugation on \(K\) and by inversion on
\(\zeta_m\) fixes \(B\) inside \(L\). It lifts to
\(L'\) by inverting the four radicals, since their
quadratic norms are one. This order-two lift
centralizes the Kummer group. Take \(W\) to be its
fixed field. Unramifiedness and splitting descend
through the extension \(L/B\) of degree two,
since the Kummer extension has odd degree.
This proves (1). The order-three-power cyclotomic
generator \(\zeta_m\mapsto\zeta_m^4\) fixes \(K\),
lifts by fixing the radicals, and commutes with
\(s\). Its action on the Kummer group is exactly
the one in (2). In particular, for \(n>1\), the
order-\(m\) factor contributes only order three
to the abelianization over \(M\).

Since \(3\) splits in \(K\), the local field of
\(M\) at three is \(\mathbb Q_3(\sqrt{-3})\).
Together with the real cyclotomic local subfield
it gives the full \(\mathbb Q_3(\zeta_m)\). Its
different has discriminant exponent
\((2n-1)3^{n-1}\). At rational primes other than
three the cyclotomic extension is unramified,
so the discriminant exponents of \(M\) are
multiplied by \(3^{n-1}\). Since \(-3D\) is
fundamental, this gives (3).

Hubbard–Washington study closely related Kummer
generators in cyclotomic layers. Their
[Theorem 2 and Proposition 6](sources/hubbard-washington-kummer-lambda-0810.1691.pdf)
use a base imaginary radicand coprime to three
and fundamental units under specified class-number
and \(\lambda\) hypotheses. The construction here
has an imaginary radicand ramified at three and
uses explicitly marked virtual elements. The
quadratic base field varies with the exponent.

## Moving prime types and the field count

**Corollary YR3.** For each fixed \(m=3^n\), the
fields of YR1 have
\(\gg_m X^{1/(4m)}/(\log X)^2\) distinct \(D\le X\).
The fields \(B\) in YR2 therefore number
\(\gg_m Y^{3/(4m^2)}/(\log Y)^2\) with absolute
discriminant at most \(Y\).
The three elliptic twists of MO1 have positive rank
and conductors \(1540D^2,13860D^2,6930D^2\).

For the same fixed \(m\), there is a rational scale
\(s_m\) in the class \(c_m\) modulo \(m\)-th powers
such that the constructions also hold on genuine
ratios \(t=s_m(P/Q)^m\), with distinct rational
primes \(P,Q\equiv1\bmod m\), at scale \(X^{1/(8m)}\).
There are \(\gg_m X^{1/(4m)}/(\log X)^{B_m}\)
distinct resulting \(D\le X\). In \(W/B\),
all primes above \(Q\) split; the Frobenius order
above \(P\) can be prescribed as any divisor of
\(m\) without changing the exponent. The cubic
factors from the three \(\delta_i\) split at both
moving primes.

For the first count apply
[Kulkarni–Levin, Theorem 1.2](sources/kulkarni-levin-hilbert-irreducibility-2111.15582.pdf)
to the degree-two map of the squarefree degree-\(8m\)
curve \(R^2=F(c_mu^m)\) and the fixed neighborhoods
of YR1. Rational parameter height \(H\) gives
\(D\ll_m H^{8m}\). The unique quadratic subfield
of \(B\) is \(M\), so \(B\) determines \(D\).
Formula (3) gives the exponent in \(Y\).
The independent marks give positive rank on the
three fixed elliptic twists as in MO3, and the
splitting at all their bad primes gives the
conductors as in BS4.

Choose \(u_*=5^N/(2^N3^N7^N11^N)\) for \(N\)
sufficiently large, and put \(s_m=c_m u_*^m\).
Then \(u=u_*P/Q\) has the same local valuations
at \(2,3,5,7,11\) as required in YR1, and no unit
restriction is needed there. At the four residue
primes one imposes fixed ratio classes on \(P/Q\).
Those classes exist since \(u_*\) is a unit there
and the necessary ratio \(m\)-th powers were
checked in YR1.

Both moving primes split in \(K\) by the boundary
values, and they split in \(L\) because they are
\(1\bmod m\). The leading coefficients in YR1
show that the three cubic marks split there.
For \(\gamma\), choose its smaller-valuation branch.
At the numerator prime \(P\), \(t\) is close to
zero and the leading trace term is \(196/t^2\).
On removing the local power of \(P\), the unit
class modulo \(m\)-th powers is that of
\[
             A=\frac{(98/3)^4}{11^6}.
\]
Indeed the coefficient is \(14^4/s_m^4\cdot 11^{-3}\),
and \(4a\equiv3\bmod m\). The unit from \(Q\) is
already an \(m\)-th power. At the denominator prime
\(Q\) the leading trace term is \(49t^2/9\), giving
the coefficient
\[
          (7/3)^4s_m^4/11^3
                        \equiv1\pmod{\mathbb Q^{\times m}}.
\]
The other branch has reciprocal unit class.
This proves the claimed moving residue tests.

At primes over two in \(\mathbb Q(\zeta_m)\),
the valuation of \(A\) is four. Thus
\(\mathbb Q(\zeta_m,\sqrt[m]{A})\) has a cyclic
subgroup of order \(m\) over \(\mathbb Q(\zeta_m)\),
and every order in that subgroup occurs at positive
density by Chebotarev. Its ramified primes are
among \(2,3,7,11\). The fixed residue moduli
\(13,31,43,61\) are disjoint from these primes.
The intersection of the two Galois fields is
therefore unramified at every finite prime of
\(\mathbb Q\), so is \(\mathbb Q\).
All prescribed Frobenius orders are compatible
with the fixed ratio residues. For each fixed
choice, Chebotarev for \(P\) and the prime number
theorem in progressions for \(Q\) give
\(\gg_m H^2/(\log H)^2\) prime pairs in a
rectangle of distinct prime sizes. The fixed-form
prime factor truncation and field multiplicity
argument of
[MP3](MIXED_THREE_POWER_CLASS_GROUPS_FROM_TWO_PRIME_FIBRES.md#two-high-order-split-prime-classes-and-an-extra-cubic-class)
work with degree \(8m\) and genus \(4m-1\).
This gives the claimed field count after a
logarithmic loss. The explicit residue tests
are local to each moving prime, including the
case \(m=3\).
