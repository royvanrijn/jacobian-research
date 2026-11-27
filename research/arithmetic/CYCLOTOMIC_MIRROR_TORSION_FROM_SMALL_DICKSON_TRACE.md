# Cyclotomic mirror torsion from a degree-\(m+1\) trace family

A norm-one trace centered at \(2\) reduces the squarefree
discriminant form to degree \(m+1\) for \(m=3^n\).
The norm is a local \(m\)-th power above three and
is noncubic at a split denominator prime. This gives
a cyclic order-\(m\) ordinary class-group quotient
over the cyclotomic layer of the quadratic mirror.
A genuine prime denominator can generate that quotient.
For the imaginary base sign it also gives a quadratic
prime ideal class of exact \(3\)-primary order \(m\).
The same base fields give positive-rank twists of
\(Y^2=X^3+4\).

## The trace, squarefree form, and local conditions

Fix \(m=3^n\) with \(n\ge2\), choose \(N=\lceil n/3\rceil\),
and put
\[
 a=3^{2N+1},\qquad A=2a^3=2\cdot3^{6N+3},\qquad
 T=2+A t^m,\qquad
 H(t)=A t(4+A t^m).
\]
Set \(K=\mathbb Q(\sqrt{H(t)})\) and
\[
                         \gamma=\frac{T+\sqrt{T^2-4}}2.
\]
Here \(T^2-4=t^{m-1}H(t)\), so \(\gamma\) belongs to \(K\)
and has norm one. We impose
\[
                   t\in 2(3+8\mathbb Z_2),
                   \qquad t\equiv1\pmod3.
                   \tag{1}
\]
Exclude the rational specializations and the single quadratic
field \(\mathbb Q(\zeta_3)\). Use a rational parameter
interval \(t>0\) for the real base sign or
\(-(4/A)^{1/m}<t<0\) for the imaginary base sign.

**Proposition SD1.** On (1), \(\gamma\) is \(m\)-virtual,
\(2\) splits in \(K\), and \(K_3=\mathbb Q_3(\zeta_3)\).
At three the norm is a local \(m\)-th power.
It can be made noncubic in \(K\) on either sign by
imposing a pole of \(t\) at a fixed prime \(q\equiv1\bmod m\)
where \(A\) is a noncube modulo \(q\).

Write \(\delta=A t^m=T-2\). At two it has even valuation
\(m+1\) and unit congruent to one modulo eight, by (1).
The factor \(4+\delta\) is a square since
\(\delta/4\in 8\mathbb Z_2\).
Thus \(T^2-4=\delta(4+\delta)\) is a square there.
At three the valuation of \(\delta\) is \(6N+3\)
and its unit is \(-1\bmod3\). As \(4+\delta\) is a
square unit, the local quadratic field is
\(\mathbb Q_3(\sqrt{-3})\). In its uniformizer
normalization, \((\gamma-1)^2/\gamma=T-2\) gives
\[
                         v_{K_3}(\gamma-1)=6N+3.
\]
By
[LD2](LOCAL_DICKSON_NORM_ONE_UNIT_TOWERS_AT_THREE.md),
this is on the pure principal norm-one subgroup with
free exponent of valuation \(3N\ge n\); thus it
is an \(m\)-th power.

At \(p\ne2,3\), whenever \(v_p(t)\ge0\), the two roots
of \(X^2-TX+1\) are units. A pole of \(t\) gives
negative trace valuation \(m v_p(t)\) and split roots
of valuations \(m v_p(t),-m v_p(t)\).
On the pole branch,
\[
                       \gamma/(A t^m)\equiv1\pmod p.
                       \tag{2}
\]
At two and three the roots are units. This proves
virtuality and the fixed noncube test. A suitable
\(q\) exists by Chebotarev: \(A\) has valuation one at
each prime above two in \(\mathbb Q(\zeta_m)\), so
the radical field for \(A\) has a cyclic translation
subgroup of order \(m\) over that field. Frobenius
generators in the subgroup give \(q\equiv1\bmod m\)
with \(A\) a noncube modulo \(q\).

## Ordinary mirror quotients and a fixed elliptic twist

**Theorem SD2.** On the slice of SD1 with \(\gamma\)
noncubic, write \(D=\operatorname{disc}K\), let \(M\)
be the quadratic field of fundamental discriminant
\(-D/3\), and put \(B=M\mathbb Q(\zeta_m)^+\).
The root field \(L\) of \(\mathcal D_m(X,1)-T\) has
degree \(m\) and full affine normal closure \(W/\mathbb Q\).
The translation-fixed field is \(B\), and \(W/B\) is
cyclic of order \(m\), everywhere unramified and
completely split above three. The following discriminants
are exact:
\[
 \begin{split}
 |\operatorname{disc} B|&=
      3^{(2n-1)m/3-1}|D/3|^{m/3},\\
 |\operatorname{disc} L|&=
      3^{(n-1)(m-1)}|D/3|^{(m-1)/2}.
 \end{split}
 \tag{3}
\]
The \(D\)-twist of \(Y^2=X^3+4\) has positive rank and
conductor \(12D^2\). For either sign, the slice with a
fixed denominator test in SD1 has at least
\[
 \begin{array}{c|c}
 \text{object and bound}&\text{number of distinct objects}\\ \hline
 |D|\le X & \gg_m X^{2/(m+1)}/(\log X)^2\\
 |\operatorname{disc}B|\le Y & \gg_m Y^{6/(m(m+1))}/(\log Y)^2\\
 |\operatorname{disc}L|\le Y & \gg_m Y^{4/(m^2-1)}/(\log Y)^2\\
 12D^2\le Z & \gg_m Z^{1/(m+1)}/(\log Z)^2 .
 \end{array}
 \tag{4}
\]

The noncube over \(K\) stays noncubic in the abelian
extension \(K(\zeta_m)/K\). A degree-three radical
over \(K\) would be nonnormal since \(\zeta_3\notin K\),
whereas every subfield of that cyclotomic extension
over \(K\) is normal. As \(K\) is disjoint from
\(\mathbb Q(\zeta_m)\), Kummer degree \(m\) and simultaneous
inversion give the full affine group, as in
[WD2](GENUS_THREE_PELL_WILD_DICKSON_FIELDS.md).
The fixed field of translations is \(M\) composed
with the real cyclotomic field. The discriminant
\(-D/3\) is fundamental because \(D\) is odd with
\(v_3(D)=1\). It is coprime to three, which splits in \(M\).

At each prime above three, \(\gamma\) becomes an \(m\)-th
power over \(K(\zeta_m)\) locally. The simultaneous
inversion descent then makes the relative cyclic
extension split over \(B\), since the auxiliary extension
has degree at most two and the cyclic degree is odd.
At primes away from three, virtuality makes the radical
unramified over \(K(\zeta_m)\); all residual inertia
in \(W/\mathbb Q\) maps injectively to the affine
multipliers. Consequently \(W/B\) is unramified.
At real places the odd cyclic degree also implies
complete splitting. Thus \(W/B\) is a quotient of the
ordinary class group. Since the two discriminants
of \(M\) and the real cyclotomic layer are coprime,
their discriminant product gives the \(B\) formula.

At primes dividing \(D/3\), the order-two inertia in
the affine normal closure gives exponent \((m-1)/2\)
for the root-field discriminant. Two is split in \(K\)
and contributes zero. At three the parameter is in
the pure, \(m\)-divisible case of
[LD3](LOCAL_DICKSON_NORM_ONE_UNIT_TOWERS_AT_THREE.md).
There is one rational root and two local factors of
degree \(3^j\) for each \(0\le j\le n-1\).
The two factors of degree \(3^j\) jointly contribute
\((2j+1)3^j-1\) to the discriminant exponent.
Their sum is \((n-1)(m-1)\), giving the \(L\) formula.

For the elliptic connection put \(z=a t^{m/3}\) and
\(\omega=\sqrt{z(z^3+2)}\in K\), with the branch chosen
so that \(\gamma=1+z^3+z\omega\). Then
\[
       X=\frac2z,\qquad Y=\frac{2\omega}{z^2},\qquad
       Y^2=X^3+4,\qquad
       \gamma=\frac{(Y+2)^2}{X^3}.
       \tag{5}
\]
Hence the standard cubic descent value \(Y+2\)
is noncubic. The point is anti-invariant over \(K\).
The cube-class evaluation vanishes on torsion of
order prime to three, and there is no nonzero
anti-invariant \(3\)-torsion over a quadratic field:
the only rational abscissa of a nonidentity
\(3\)-torsion point is \(X=0\), since the residual
division factor \(X^3+16\) is irreducible, and
both points there are rational. Thus the point has
infinite order on the twist. The base conductor is
\(108\), and its twist by \(-3\) retains exponent
three at \(3\). The exact checker verifies these
two local conductor inputs. At two the twist
is locally trivial, and at every other prime
dividing \(D\) the base curve has good reduction.
This gives \(108(D/3)^2=12D^2\).

The polynomial \(H(t)\) is squarefree of degree
\(m+1\). Its factors have disjoint roots and the
nonzero roots of \(4+A t^m\) are simple. On a fixed
translated local neighborhood realizing (1), the
fixed pole test, and either real sign, use
[Kulkarni–Levin, Theorem 1.2](sources/kulkarni-levin-hilbert-irreducibility-2111.15582.pdf)
for the degree-two map of \(H\). Rational parameters
of height \(U\) give \(\gg_m U^2/(\log U)^2\) distinct
quadratic fields with \(|D|\ll_m U^{m+1}\). The
discriminant formulas distinguish all corresponding
objects in each sign and give (4).

## Genuine numerator and denominator primes

**Corollary SD3.** For either base sign one may instead
take \(t=\pm2 P/Q\) with genuine distinct primes.
The numerator prime \(P\) ramifies in \(K\) and \(M\)
and splits completely in the relative extension \(W/B\).
The denominator prime \(Q\) splits completely in \(B\),
and every prime of \(B\) above it has Frobenius order
\(m\) in \(W/B\). The distinct quadratic discriminants
with \(|D|\le X\) number
\(\gg_m X^{2/(m+1)}/(\log X)^{C_m}\)
on these prime ratios. If \(K\) is imaginary, either
prime ideal above \(Q\) has class of exact \(3\)-primary
order \(m\). One can force its full order to be \(2m\)
on a fixed further slice, with the same count exponent.

Choose \(Q\equiv1\bmod m\) with prescribed generator
Frobenius for the fixed radical field of \(A\) over
\(\mathbb Q(\zeta_m)\), as in SD1.
A fixed class of such primes modulo eight can be
selected with positive density by Chebotarev in the
compositum with \(\mathbb Q(\zeta_8)\).
Place \(P\) in the fixed nonzero residues modulo
eight and three making (1) hold. Use the plus sign
and a positive parameter interval for real \(K\),
and the minus sign and an interval inside
\((-(4/A)^{1/m},0)\) for imaginary \(K\).
Both can be realized by positive primes in rectangles
with side lengths comparable to their heights.
Chebotarev and the prime number theorem in progressions
give \(\gg_m U^2/(\log U)^2\) prime pairs at height \(U\).

At \(P\), the radicand has odd valuation \(m\) before
square removal. The norm-one unit is congruent to one
in the ramified local field, since its polynomial
reduces to \((X-1)^2\). It is an \(m\)-th power there
by Hensel, so the prime splits over \(B\) by the same
odd-degree descent as at three. At \(Q\) equation (2)
gives the residue \(A\) after removing the \(m\)-th
power and valuation. Both \(K\) and the full cyclotomic
field split at \(Q\), and the two local branches
give inverse residues of Frobenius order \(m\).
In particular the residue enforces the global noncube
without any auxiliary fixed pole condition.

At all other primes \(\gamma\) is a unit. Thus, for
one orientation,
\[
                         (\gamma)
                         =(\mathfrak q/\bar{\mathfrak q})^m.
\]
The full order of \([\mathfrak q]\) divides \(2m\).
In imaginary \(K\ne\mathbb Q(\zeta_3)\) all units
are cubes. If the \(3\)-primary order of
\([\mathfrak q]\) were smaller than \(m\), the divided
ideal \((\gamma)^{1/3}\) would be principal, and
\(\gamma\) would be a cube, a contradiction.

For full order \(2m\), choose a fixed prime
\(\ell\equiv1\bmod4\), not over \(6\), where
\(-4/A\) has an \(m\)-th root and \(\ell\nmid m\).
Such a prime is provided by complete splitting in
the fixed cyclotomic and radical fields. Set the
parameter in a fixed class modulo \(\ell^2\) where
\(4+A t^m\) has valuation one. This is possible
by its nonzero derivative at the chosen root.
Make \(Q\) a quadratic nonresidue modulo \(\ell\).
The latter condition can be imposed by fixing a
unit residue of \(Q\bmod\ell^2\), independently of
the previous \(Q\)-conditions: the governing fields
for those conditions ramify only at \(2,3\) and
are disjoint from \(\mathbb Q(\zeta_{\ell^2})\).
The congruence on \(t\) is then a fixed residue for
\(P\) modulo \(\ell^2\).
Since \(\ell\) is a prime factor of the fundamental
discriminant of \(K\), \(K(\sqrt\ell)/K\) is an
unramified genus extension, since the two coprime\nfundamental discriminants \(\ell,D/\ell\) give absolute\ncompositum discriminant \(D^2\). Its character at
\(\mathfrak q\) is \((Q/\ell)=-1\), giving the even
part of order two.

For both prime counts, the separable binary form
for \(H\) has degree \(m+1\) after fixing the scale
\(2\). The prime-factor truncation of
[PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#quantitative-count-of-genuine-fields)
keeps all but \(O_m(U^2/(\log U)^4)\) pairs with
\(\omega(D)=O_m(\log\log U)\).
For the fixed form, its genus \((m-1)/2\ge4\) and the
same fixed-field \(2\)-descent argument as in PP2
give rank \(O_m(1+\omega(D))\) for each quadratic-twist
Jacobian. Uniform Mordell–Lang then bounds the
multiplicity of each field by a power of \(\log U\).
Taking \(U\asymp_m X^{1/(m+1)}\) proves the count.

## Scope and evidence

The classical rational norm construction and the
cyclotomic Kummer arguments are the inputs. The comparison
with
[Hubbard–Washington, Theorem 2 and Proposition 6](sources/hubbard-washington-kummer-lambda-0810.1691.pdf)
is as in YR2. Their fixed-field unit and Iwasawa
hypotheses differ from this varying quadratic trace
family, where the discriminant counts run separately
for each fixed \(m\). Relative to YR3, this family uses
the smaller degree \(m+1\) and gives one cyclic quotient
without the three extra marked cubic factors.
The construction makes no fixed-field assertion as \(n\)
grows. In real \(K\) the absolute class order of the
prime can still be affected by the unit kernel.
The [exact checker](verify_cyclotomic_mirror_small_dickson_trace.py)
replays the polynomial reduction, elliptic map, finite
local coefficients and representative Frobenius tests.

## Tame inertia gives the exact generic torsion

**Theorem SD4.** Let \(m>1\) be odd, \(A\in\mathbb Q^\times\),
and let \(J\) be the Jacobian of
\[
                    C:\ y^2=A t(4+A t^m).
\]
If \(\gcd(m,v_2(A))=1\), its rational torsion of odd
order is cyclic of order \(m\), generated by the
difference of the two infinities. If additionally
\(\gcd(m,2-v_2(A))=1\), it is the full rational torsion
subgroup. In particular this holds for the scales
in SD1. On those scales every rational norm-one
function on \(C\) with divisor divisible by \(m\)
has class modulo \(m\)-th powers in the
cyclic subgroup generated by \(\gamma\).

Write \(u=1/t\) and
\(\rho=y/(A t^{(m+1)/2})\), so that the curve is
\[
                         \rho^2=1+(4/A)u^m .
\]
The points \(Q_\pm=(0,\pm1)\) and the single point
\(Q_\infty\) at infinity are rational. The quotient
by \(\tau:u\mapsto\zeta_m u\) is the \(\rho\)-line,
totally ramified at \(Q_+,Q_-,Q_\infty\). The last
assertion uses \(m\) odd, since the orders in
\(\operatorname{div}(u)\) are \(1,1,-2\).

The group of geometric fixed points \(J^\tau\)
has a direct divisor description. Every fixed line
bundle can be linearized for the cyclic action by
scaling a lift of its generator, and Hilbert 90 on
its generic fiber gives an invariant divisor
representative. An invariant divisor is a pullback
plus multiples of the three ramification points.
Modulo \(m\), its three coefficients sum to zero in
degree zero. A principal invariant divisor comes
from an eigenfunction \(u^j f(\rho)\); its residue
vector is therefore a multiple of \((1,1,-2)\).
If the residues vanish, the divisor is a pullback
of degree zero from the projective line and is
principal. Every residue vector of sum zero occurs.
Hence
\[
 J^\tau\simeq
 \{(b_+,b_-,b_\infty)\in(\mathbb Z/m)^3:
                      b_++b_-+b_\infty=0\}
                   /\langle(1,1,-2)\rangle
             \simeq\mathbb Z/m .
\]
The class of \(Q_+-Q_-\) generates it, since it
represents twice the generator after the displayed
relation and \(m\) is odd.

To bound rational odd torsion, work at two. The local
field \(\mathbb Q_2(\zeta_m)\) is unramified.
Adjoining \(\alpha=A^{1/m}\) gives a cyclic totally
ramified extension of degree \(m\): a root has
valuation \(v_2(A)/m\) of denominator \(m\).
Over this extension use \(u=\alpha X\),
\(\rho=1+2Y\). The good-reduction model becomes
\[
                              Y^2+Y=X^m .
\]
The affine special fiber is smooth by the partial
derivative in \(Y\). At infinity use \(v=1/X\) and
\(Z=v^{(m+1)/2}Y\), with equation
\(Z^2+v^{(m+1)/2}Z=v\). The derivative in \(v\)
is nonzero at the unique point \(v=Z=0\).
These two charts give a smooth proper model, finite
over the two charts of the projective line. Odd
torsion on the Jacobian of this model is unramified.
A generator of local inertia on \(\alpha\) acts
on \(C\), under this comparison of models, as
\(X\mapsto\zeta_m X\) or its inverse.
Thus a rational odd-order point of \(J\) must
belong to \(J^\tau\). The rational generator
\(Q_+-Q_-\) proves equality.

If \(\gcd(m,2-v_2(A))=1\), the nonzero branch
polynomial \(t^m+4/A\) is irreducible over
\(\mathbb Q_2\), by the denominator of the valuation
of its roots. Together with the rational branch
point \(t=0\), the two branch orbits have sizes
\(m,1\). A geometric Jacobian \(2\)-torsion point
is an even subset of the branch points modulo
complement. Taking the representative omitting
the rational branch point, a rational class can
only use the empty set, since the other orbit
has odd size. This eliminates rational \(2\)-torsion
and hence all rational even torsion.

For the last assertion take any norm-one rational
function \(\lambda\) with \(m\)-divisible divisor.
The divided divisor defines a rational point of
\(J[m]\), which is a multiple of the class of the
divided divisor of \(\gamma\). The latter follows
from the two branches at \(t=\infty\), with
valuations \(m,-m\), and no other zero or pole.
The difference from the appropriate power of
\(\gamma\) has the form \(c f^m\) with
\(c\in\mathbb Q^\times,\ f\in\mathbb Q(C)^\times\),
by descent of a rational principal divisor.
Its norm is one, so \(c^2\) is a rational
\(m\)-th power. Hence \(c\) is an \(m\)-th power,
as \(m\) is odd.

Jędrzejak’s Theorem 1 already gives the order-\(m\)
lower subgroup and general prime-divisor restrictions
for \(y^2=x(x^m+c)\), which is this curve after
dividing \(y\) by \(A\). The exact upper bound on
this subfamily comes from the full tame inertia
action at two. The ramification-divisor subgroup
also fits the structure theorem of
Wawrów, Proposition 1 for coprime superelliptic
degree and exponent. These comparisons distinguish
the lower subgroup from the rational upper bound.
[Jędrzejak, Theorem 1](https://www.impan.pl/shop/en/publication/transaction/download/product/91624)
and [Wawrów, Proposition 1](https://arxiv.org/pdf/1911.05431).

## The rational fibers and pole support

**Proposition SD5.** For the scales in SD1 the rational
points on \(C\) are exactly its two infinities and the
branch point \(t=0\). Consequently any rational norm-one
function whose trace has all poles at rational
parameters is a signed integral power of \(\gamma\).
More generally, the rational norm-one functions with
\(m\)-divisible divisors have exactly the one direction
in SD4. A fixed finite list of other rational norm-one
functions can yield \(m\)-virtual values on only
\(O_m(U^{3/2}\log U)\) parameters of height at most \(U\),
with the constant also depending on the list. This bound
uses the same fixed-parameter height as (4).

By (5), any rational point with \(t\) finite and nonzero
would give a rational point of \(E:Y^2=X^3+4\) with
\(X\) finite and nonzero. Here
\(E(\mathbb Q)=\{O,(0,2),(0,-2)\}\).
For a direct rank proof apply
[Cohen–Pazuki, Definition 1.3 and Proposition 2.2](sources/cohen-pazuki-2009-three-isogeny.pdf).
The image of the first descent on \(E\) has order three.
The norm identity \(Y^2-4=X^3\) makes all
valuations of a descent value outside two divisible
by three, and a marked torsion point gives a nontrivial
class of \(2\). Its dual curve is \(E':Y^2=X^3-108\),
and the second descent uses \(Y+6\sqrt{-3}\) up to
conjugation. Away from \(2,3\), the pair of conjugate
elements has difference a local unit; at any pole
the rational ordinate has valuation divisible by three.
Thus every valuation of the second value is divisible
by three. At \(2,3\) the same follows from the cube
norm and the unique prime of \(\mathbb Q(\sqrt{-3})\)
over each. As this ring of integers is Euclidean, its
virtual classes come from units, so the only remaining
possibility is the class of a cube root of unity.

In fact no rational point on \(E'\) has \(v_3(X)>0\).
For \(v_3(X)\ge2\), \(X^3-108\) has odd valuation
three. For \(v_3(X)=1\), write \(X=3x\). When
\(x\equiv-1\bmod3\), \(27(x^3-4)\) has odd valuation
three. When \(x\equiv1\bmod3\), \(x^3\equiv1\bmod9\),
and \((x^3-4)/3\equiv-1\bmod3\) is a nonsquare
unit, again impossible. For \(v_3(X)=0\), the descent
value is congruent to a rational unit modulo
\((\sqrt{-3})^2\). For \(v_3(X)<0\) the rational ordinate
has valuation a multiple of three, so division by
its \(3\)-power is division by a cube and gives
the same unit congruence. A nontrivial cube root of
unity modulo cubes cannot have that congruence,
because the cube of a local unit reduces to a
rational unit modulo \((\sqrt{-3})^2\).
Thus the second descent image is trivial and the
rank is zero. Reduction at \(5,7\) gives point
counts \(6,3\), bounding rational torsion by
three. This proves the rational point claim.

For the pole assertion, if one conjugate of a
norm-one function has a pole at a point, its
conjugate has a zero there and its trace has a pole.
Over a rational nonsplit fiber there is only one
closed point, fixed by conjugation, so the
norm-one divisor has zero coefficient there.
The same applies at the ramified fiber \(t=0\).
Thus rational trace poles can occur only at
infinity. The full subgroup of principal divisors
supported there is generated by \(m(Q_+-Q_-)\)
by SD4, and a rational constant of norm one is
\(\pm1\).

Finally any fixed norm-one function outside the
generic span of SD4 has some divisor coefficient
not divisible by \(m\). Apply
[VB1–VB2](QUADRATIC_VIRTUALITY_BINARY_FORM_POWER_OBSTRUCTIONS.md)
to each such function on this quadratic map.
The finite union gives the height bound. It leaves
open field-dependent classes, thin subfamilies and
representations from different curves.
