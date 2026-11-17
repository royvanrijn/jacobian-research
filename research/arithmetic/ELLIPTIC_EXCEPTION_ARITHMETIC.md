# Four quadratic fields and infinitely many field fibres on an elliptic exception

The exceptional base pair `(pi,ell)=(128/25,9/10)` in the
[resolvent fibration](FIXED_MAP_RESOLVENT_FIBRATION.md) has more arithmetic
structure than infinitude of matching-squareclass targets alone. Its
quadratic factor is always irreducible at admissible rational points,
and its field belongs to an exact list of four imaginary quadratic
fields. The cubic factor is irreducible except at finitely many rational
points of this elliptic slice. Each of the four quadratic fields occurs
in infinitely many distinct complete `2+3` field-fibre targets.
In fact each of these four quadratic classes contains infinitely many
isomorphism classes of cubic fields and of complete fibre algebras.

These statements do not assert infinitely many Hasse failures: the
finite-prime local conditions remain to be imposed.

## 1. The fixed elliptic slice

Throughout this note fix `pi=128/25`, `ell=9/10`. The normalized inverse
polynomial and its factor chart are those of equations (1.1)--(1.5) of
the resolvent note. In particular, put

\[
 r=-\frac{12}{5},\quad s=s_0-v,\quad
 s_0=\frac{1816347}{819200},\quad
 t=\frac{10}{9}v^2-\frac{10s_0}{9}v+\frac{125}{1024},
\]
\[
 q(W)=W^2+\frac9{10}W+v,\qquad
 h(W)=W^3+rW^2+sW+t.                                 \tag{1.1}
\]

The target is the fixed-map target

\[
 \left(-\pi,-\frac{32\pi^2}{9}
      \left(t+\ell s+vr+\frac5{4\pi}\right),
      -\frac{32vt}{9}\right).                       \tag{1.2}
\]

Define

\[
 L(v)=819200v-243483,\quad v_L=\frac{243483}{819200},
\]
\[
 Q(v)=20480000v^2-110592000v+143467857,\quad
 J(v)=(400v-81)Q(v),\quad K=41231686041600000000.
\]

The cubic discriminant and normalized matching-resolvent curve are

\[
 \Delta(v)=-\frac{L(v)^2Q(v)}{412316860416000000},
 \qquad E:\ Z^2=\frac{J(v)}K.                       \tag{1.3}
\]

The original resolvent variable is `Y=L(v)Z`. Retain the admissible open

\[
 D(v)\Delta(v)\operatorname{Res}_W(q,h)\ne0,
 \qquad D(v)=\frac{81}{100}-4v.                     \tag{1.4}
\]

As proved in Propositions 3.7--3.8 of the resolvent note, this elliptic
curve is isomorphic to

\[
 T^2=X^3-10758382X+6600759744,
 \quad X=\frac{6400}{3}v-3984,\quad
 T=\frac{20971520000}{3}Z,                           \tag{1.5}
\]

and contains the nontorsion point
`P=(199200/49,59459616/343)` in `(X,T)` coordinates. Translating the
rational point of order two to the origin, put `x=X+3552`, `y=T`.
Then

\[
 E_0:\ y^2=x^3-10656x^2+27091730x,
 \qquad D=-\frac{3x}{1600}.                         \tag{1.6}
\]

Write `p_0=2709173`; this is prime, and the coefficient of `x` is
`b=10p_0`. The point of order two is `T_0=(0,0)` in this model.

## 2. Exact quadratic-field classification

**Theorem 2.1 (four quadratic fields, with no reducible quadratic
exception).** At every admissible rational point of `E`, the quadratic
factor `q` is irreducible and its quotient field is one of

\[
 \boxed{\mathbb Q(\sqrt{-3}),\quad
 \mathbb Q(\sqrt{-6}),\quad
 \mathbb Q(\sqrt{-40637595}),\quad
 \mathbb Q(\sqrt{-81275190}).}                       \tag{2.1}
\]

Each field in this list occurs at infinitely many distinct admissible
rational matching-resolvent target images (1.2).

**Proof.** Use the usual two-isogeny squareclass map

\[
 \alpha:E_0(\mathbb Q)\longrightarrow\mathbb Q^*/\mathbb Q^{*2},
 \quad \alpha(O)=1,\quad\alpha(T_0)=b,\quad
 \alpha(x,y)=x\pmod{\mathbb Q^{*2}}.                 \tag{2.2}
\]

The standard descent interface is given in A. Dujella, *Number Theory*
(2021), section 15.5, pp. 507--509, equations (15.29)--(15.31), in
[the retained primary section](sources/dujella-number-theory-2021-section15-5.pdf)
and on the [author's book page](https://web.math.pmf.unizg.hr/~duje/numbertheorybook.html).
Here its required facts also have short elementary proofs. A line
`y=mx+n` meets the cubic in three `x` coordinates whose product is
`n^2`, proving the squareclass homomorphism law when none is zero.
For a line through `T_0`, the other two coordinates have product `b`.
Vertical lines and the tangent cases give the same law by their
intersection multiplicities. These identities include all cases in
(2.2).

The image is represented by signed squarefree divisors of `b`. Indeed,
if `v_p(x)<0`, the term `x^3` has strictly smallest valuation in (1.6),
forcing `v_p(x)` even. If `v_p(x)>0` and `p` does not divide `b`, the
valuation of the right side equals `v_p(x)`, again forcing it even.
Only primes dividing `b` can therefore occur to odd exponent.

No negative class occurs. For real `x<0`, the factor
`x^2-10656x+27091730` is positive, making the right side of (1.6)
negative. Thus every finite rational point has `x>=0`, and every
admissible point has `x>0`.

In (1.6), the nontorsion point `P` has

\[
 x(P)=\frac{373248}{49}=2\left(\frac{432}{7}\right)^2.
\]

Together with `alpha(T_0)=10p_0`, this shows that the image contains

\[
 H=\{1,2,5p_0,10p_0\}.                             \tag{2.3}
\]

The other four positive divisor classes are exactly the coset `5H`.
To exclude this coset, it suffices by the homomorphism law to exclude
the class `5`. If that class occurred, write `x=5(M/e)^2` with
coprime integers `M,e`. Substitution in (1.6) gives the necessary
integral quartic equation

\[
 N^2=5M^4-10656M^2e^2+5418346e^4.                   \tag{2.4}
\]

The rational square root here is integral because its square is an
integer. Modulo sixteen the coefficients of (2.4) are `5,0,10`.
Since `M,e` are not both even, its right side is respectively `5`,
`10`, or `15` according to their parity. None is a square modulo
sixteen. Consequently `alpha(E_0(Q))=H` exactly.

By (1.6), the squareclass of the quadratic discriminant is `-3alpha`.
For `alpha=1,2,5p_0,10p_0`, this is precisely (2.1). It is negative
at every admissible point, so `q` is always irreducible.

Finally each class in `H` has a rational representative among
`O,P,T_0,P+T_0`. For any such representative `R`, the points
`R+2nP` are infinitely many distinct rational points and have the same
`alpha` class. The function `v` has degree two, and the conditions
(1.4) remove only finitely many values of `v`; hence infinitely many
admissible parameters remain in each class. The target map has fibres
of size at most two because its `b(v)` coordinate is quadratic. Thus
each field occurs in infinitely many distinct target images. \(\square\)

**Theorem 2.2 (rank two, without a basis assertion).** The elliptic
curve `E` has Mordell--Weil rank exactly two over `Q`.

**Proof.** The curve two-isogenous to (1.6) is

\[
 E_0':\quad y^2=x^3+21312x^2+5183416x,
 \qquad 5183416=2^3 7^3\cdot1889.                   \tag{2.5}
\]

Its descent image is represented by the sixteen signed squarefree
divisors of `2*7*1889`. For a candidate `d`, the necessary primitive
quartic is

\[
 N^2=dM^4+21312M^2e^2+(5183416/d)e^4.               \tag{2.6}
\]

The following complete residue exclusions leave only four candidates.

| Modulus | Classes excluded by primitive insolubility |
| --- | --- |
| `16` | `2,7,3778,13223,-1,-14,-1889,-26446` |
| `25` | `-2,-7,-3778,-13223` among those still surviving |

For completeness, these are projective primitive checks: modulo `p^k`,
normalize `e=1` if `e` is a unit; otherwise normalize `M=1` and let `e`
range over multiples of `p`. Scaling the homogeneous quartic by a
fourth power preserves whether its value is a square. The checker
enumerates every residue in these two charts, so the table excludes
every primitive integer solution in the indicated classes.

The surviving classes `1,14,1889,26446` all occur. The identity gives
`1`, the rational point `(0,0)` gives `26446`, and

\[
 Q'=\left(\frac{3584}{25},\frac{4301248}{125}\right)\in E_0'(\mathbb Q),
 \qquad \frac{3584}{25}=14\left(\frac{16}{5}\right)^2
\]

gives `14`. Their product gives `1889`. Hence both descent images
have size four, using Theorem 2.1 for `E_0`. The two-isogeny rank
formula in the same retained Dujella section, p. 508, is

\[
 2^{r+2}=|\operatorname{im}\alpha|
         |\operatorname{im}\alpha'|.
\]

It gives `2^{r+2}=16`, so `r=2`. This imports the stated classical
rank formula and does not claim that particular points form a
saturated basis. \(\square\)

**Corollary 2.3 (exact torsion subgroup).** One has
`E(Q)_tors={O,T_0}`, isomorphic to `Z/2Z`.

**Proof.** The primes `3,11` are odd primes of good reduction. Direct
finite-field point counts give `#E(F_3)=4` and `#E(F_11)=10`.
The rational torsion subgroup injects into either group, by
Proposition 15.7, p. 490, of Dujella's *Number Theory*, in
[the retained primary section 15.3](sources/dujella-number-theory-2021-section15-3.pdf).
Its order therefore divides `gcd(4,10)=2`. The known nonidentity
point `T_0` has order two, so the subgroup is exactly as stated.
\(\square\)

## 3. The cubic-root exception has genus four

**Theorem 3.1 (finitely many reducible cubic factors).** The smooth
projective normalization `C` of the cubic-root incidence over `E`,

\[
 Z^2=J(v)/K,\qquad h(W,v)=0,                         \tag{3.1}
\]

is geometrically integral of genus four. In particular only finitely
many admissible rational points of `E` have a reducible cubic `h`.

**Proof.** First consider the degree-three root cover over the
projective `v`-line. It is geometrically irreducible. At `v=infinity`,
the coefficient valuations of `h` in powers `W^0,W^1,W^2,W^3` are
`-2,-1,0,0`, respectively. The lower Newton edge joins `(0,-2)` to
`(3,0)`, and each root has valuation `-2/3`. Equivalently, comparing
the valuations of the four terms directly shows that the smallest
valuation can cancel only at this value. There is no root in the
Laurent series field `Qbar((1/v))`; a cubic without a root is
irreducible. The cover is therefore geometrically integral, with
ramification index three at infinity.

Its finite branch values can occur only at zeros of its discriminant
in (1.3). The two distinct roots of `Q` are simple discriminant zeros,
so each has inertia a transposition. At the double discriminant value
`v_L`, the exact expansion is

\[
 h\left(\frac45+w,v_L+z\right)
   =w^3-wz+\frac{10}{9}z^2-\frac{42649}{16384}z.     \tag{3.2}
\]

The nonzero coefficient of `z` makes the root curve smooth there and
expresses `z` as a unit times `w^3` to leading order. Thus this value
has ramification index three.

The elliptic double cover `E -> P^1_v` branches exactly at the two
roots of `Q`, at `v_0=81/400`, and at infinity. These branch values
are distinct, and `v_L` is none of them: explicitly
`Q(v_0)=121912785` and `Q(v_L)=3683347832025/32768` are nonzero,
with `v_L!=v_0`. Base-changing the irreducible cubic extension to this
quadratic extension leaves it irreducible, since a degree-three root
cannot lie in an extension of degree two. This argument applies over
`Qbar` as well, proving geometric integrality of `C`.

The ramification of the degree-three map `C -> E` is now explicit.
At a root of `Q`, the index-two ramification is cancelled by the
index-two base change: locally normalizing `sigma^2=tau^2` gives
unramified branches. At `v_L`, the elliptic cover is unramified and
has two geometric preimages, each with index three in `C`. At
infinity the indices three and two are coprime; normalizing
`sigma^3=tau^2` gives index three over the elliptic curve. The extra
branch `v_0` of `E` introduces no ramification in `C -> E`, because
the root cover is unramified there. There is no other ramification.
Thus the total ramification degree is `2+2+2=6`. Riemann--Hurwitz gives

\[
 2g(C)-2=3(2g(E)-2)+6=6,\qquad g(C)=4.              \tag{3.3}
\]

Apply Faltings's finiteness theorem to this smooth projective
geometrically integral curve. The primary import is Satz 7, p. 365,
of G. Faltings, *Endlichkeitssatze fur abelsche Varietaten uber
Zahlkorpern*, Inventiones Mathematicae **73** (1983), 349--366:
[primary scan](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf),
[retained copy](sources/faltings-1983-endlichkeitssaetze.pdf).
It gives finiteness, without an effective list or uniform bound.

On the admissible open (1.4), the cubic-root incidence is finite etale
over the smooth curve `E`; a rational cubic root therefore gives a
rational point of its normalization `C`. Every reducible cubic over
`Q` has a rational root. Such exceptions are consequently finite.
\(\square\)

**Corollary 3.2 (infinitely many genuine `2+3` field fibres in each
quadratic class).** For each of the four fields `K_2` in (2.1),
infinitely many distinct rational targets of the fixed map `Phi` on
this elliptic slice have complete fibre algebra

\[
 K_2\times L_3,
\]

where `L_3` is a non-Galois cubic field with quadratic resolvent `K_2`.

**Proof.** In each quadratic class, Theorem 2.1 gives infinitely many
admissible target images. Theorem 3.1 removes only finitely many
parameters with reducible `h`. The remaining quadratic and cubic
factors are irreducible, squarefree, coprime, and have matching
nonsquare discriminant classes. The exact fixed-map inverse chart
identifies their product as the complete fibre. \(\square\)

**Corollary 3.3 (infinitely many distinct cubic-field and fibre-algebra
classes in each quadratic class).** For each `K_2` in (2.1), the
cubic fields `L_3` and complete algebras `K_2 times L_3` in Corollary
3.2 represent infinitely many isomorphism classes over `Q`.

**Proof.** Fix one cubic number field `L`. The genus-four curve `C`
has only finitely many `L`-rational points by the same imported
Faltings theorem, which applies over every number field. If an
admissible cubic factor has quotient field isomorphic to `L`, an
isomorphism sends the class of `W` to a root of `h` in `L`. Together
with its rational elliptic point, this gives an `L`-point of the smooth
incidence normalization `C`. Hence only finitely many admissible target
parameters can yield this fixed field isomorphism class. Corollary 3.2
gives infinitely many target images in each quadratic class, so
infinitely many cubic-field classes are necessary. The degree-two and
degree-three field components of a product algebra are intrinsic;
therefore distinct cubic classes also give distinct complete
fibre-algebra classes for fixed `K_2`. No uniform multiplicity bound
or quantitative distinct-class asymptotic is asserted. \(\square\)

The local Hasse-failure test from the ambient classification still
requires checking the finite ramified primes. Neither Theorem 3.1 nor
Corollary 3.2 proves that infinitely many of these targets pass those
tests.

An [explicit Hasse-failing target on this same elliptic slice](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md)
has quadratic field
`Q(sqrt(-81275190))` and comes with exact certificates at every required
finite prime. It establishes existence with a different quadratic field
from the two earlier lines, while leaving Hasse infinitude on this
elliptic slice open.

## 4. Verification scope

Run `.venv/bin/python research/arithmetic/verify_elliptic_exception_arithmetic.py`.
The [checker](verify_elliptic_exception_arithmetic.py) verifies the curve
translations, point and discriminant identities, the full divisor-class
partition, the primitive quartic obstruction modulo sixteen,
dual-descent primitive residue tables modulo sixteen and twenty-five,
the dual-curve witness and rank formula arithmetic, the cubic
discriminant and exact triple-root expansion, the branch separations and
valuation arithmetic, the genus-four Riemann--Hurwitz calculation, and
the two finite-field point counts determining torsion.
The [receipt](elliptic_exception_arithmetic_checks.json) records source
hashes and exact controls. The written valuation and ramification
arguments establish the global geometric claims; the retained Faltings
theorem is a stated external import. No formal verification,
independent whole-proof replay or external review is claimed.
