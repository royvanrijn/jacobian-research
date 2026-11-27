# Arithmetic intermediate fields of a two-prime mirror extension

Two independent prime-path norms give an affine
Galois extension with a rank-two translation
group. The intermediate fields disjoint from the
cyclotomic mirror have an exact classification.
Their discriminant identifies the type of the
translation quotient, and the two actual moving
primes give a sharper classification within
each type. The fixed-resolvent count is a subgroup
count in a local-power virtual group.

Fix an odd prime \(\ell\), let \(m=\ell^n\), and
take one of the quadratic fields in
[OPR1–OPR2](ODD_PRIME_POWER_RECIPROCAL_NORM_CURVES.md)
on the basis prescription for the moving primes
\(P,Q\). Write \(D=\operatorname{disc}K\),
\(E=(K\mathbb Q(\zeta_m))^{(\sigma_K,-1)}\), and
\(W/E\) for the resulting unramified relative
extension with group \(V=(\mathbb Z/m)^2\).
The two actual primes give a basis \(e_P,e_Q\)
of translations, with an orientation at each
prime. The extension \(W/\mathbb Q\) is Galois and
\(\operatorname{Gal}(E/\mathbb Q)=U_m
=(\mathbb Z/m)^\times\) acts on \(V\) by scalars.
The action follows from Kummer: on the
radicands the cyclotomic coordinate is trivial,
the quadratic coordinate acts by inversion,
and on roots of unity the two actions give
the indicated scalar on the diagonal quotient.

## Discriminants identify the translation type

**Theorem AFM1.** The embedded intermediate
fields \(L\subseteq W\) with \(L\cap E=\mathbb Q\),
taken up to \(\mathbb Q\)-isomorphism, correspond
bijectively to subgroups \(T\le V\).
Their degrees are \(|G_T|\), where \(G_T=V/T\).
For each \(j\ge1\) let \(h_j(G_T)\) be the number
of elements of order \(\ell^j\) in \(G_T\).
Their local algebra at \(\ell\) is
\[
 \mathbb Q_\ell\ \times\
 \prod_{j\ge1}\mathbb Q_\ell(\zeta_{\ell^j})^{
                                 h_j(G_T)/\varphi(\ell^j)}.
\]
Their absolute discriminant is
\[
 |D|^{(|G_T|-1)/2}\ell^{\,f(G_T)},\qquad
 f(G_T)=\sum_{j\ge1} jh_j(G_T)
                       -\frac{|G_T|-1}{\ell-1}.
 \tag{1}
\]
If \(K\) is imaginary then \(L\) is totally
real. If \(K\) is real then \(L\) has one
real place and \((|G_T|-1)/2\) pairs of
complex places. At any fixed degree the
\(\ell\)-part of the discriminant determines
the isomorphism type of the translation
quotient.

For the classification, an extension of
\(U_m\) by an odd abelian translation
group with scalar action has conjugate
complements. In fact a lift \(t\) of
\(-1\) acts by inversion and has
\(t^2=1\), since \(t\) fixes \(t^2\).
Every other lift is conjugate to \(t\)
by a unique translation because squaring
is an automorphism of the odd group.
It follows that the centralizer of \(t\)
maps isomorphically onto \(U_m\) and
is a complement. All complements have
this form. Thus \(W/\mathbb Q\) has the
affine group \(V\rtimes U_m\).

Write \(H=\operatorname{Gal}(W/L)\).
The condition \(L\cap E=\mathbb Q\)
means \(H\) projects onto \(U_m\).
Then \(T=H\cap V\) is a subgroup
invariant under scalar action, hence
normal in the full group. Modulo \(T\)
the subgroup \(H\) is a complement and
its conjugacy class is determined by
\(T\). This proves existence and
uniqueness up to conjugacy. Distinct
subgroups \(T\) give distinct classes:
any isomorphism of number fields inside
the normal field \(W\) extends to an
automorphism of \(W\) over \(\mathbb Q\),
and \(T\) is normal. Explicitly the
permutation representation is the
affine action on \(G_T\).

At \(\ell\) the relative extension
splits over \(E\), which has one prime
with the full local cyclotomic completion.
Its decomposition group in \(W\) is a
complement. So the local action is
multiplication by the cyclotomic units
on \(G_T\). A point of exact order
\(\ell^j\) has orbit size \(\varphi(\ell^j)\)
with the local cyclotomic field of that
degree. Its discriminant exponent is
\((j(\ell-1)-1)\ell^{j-1}\). Summing
gives \(f(G_T)\). Away from \(\ell\)
only the primes dividing \(D\) ramify.
Their tame inertia has order two and
acts by a translate of inversion on
\(G_T\). It fixes one point and pairs
the rest, giving (1). At infinity the
mirror is real on the imaginary
quadratic side and totally complex
on the real quadratic side. In the
latter case complex conjugation is
also a translate of inversion. This
gives the signatures.

Write \(G_T=C_{\ell^a}\times C_{\ell^b}\)
with \(a\ge b\ge0\).
For fixed \(a+b\), increasing \(a\)
decreases
\(|G_T[\ell^j]|=\ell^{\min(a,j)+\min(b,j)}\)
at every threshold where it changes,
with a strict decrease for at least
one threshold. Hence the sum of
logarithmic element orders, the
first term of \(f(G_T)\), strictly
increases. Since \(\ell\nmid D\),
the valuation at \(\ell\) in (1)
recovers the type. The signature and
the remaining discriminant then
recover the signed \(D\) for \(a\ge1\).
This extends the one-direction
cyclotomic orbit calculation of
[YR5](MARKED_OCTIC_CYCLOTOMIC_THREE_POWER_REFLECTION.md).

## Exact moving-prime order profiles

**Theorem AFM2.** For a fixed type \((a,b)\)
with \(n\ge a\ge b\ge0\) and \(a\ge1\), let
the ordered local degrees at \(P,Q\) be
\((\ell^i,\ell^j)\). A field of AFM1 with
this profile exists exactly when
\(\max(i,j)=a,\ \min(i,j)\ge b\).
The number of fields at the indicated
profile from the same \(W\) is one
when \(a=b\) or when \(\min(i,j)=b\).
In the remaining cases it is
\((\ell-1)\ell^{\min(i,j)-b-1}\).
At a moving prime with local degree
\(\ell^i\), every prime factor is
unramified with that residue degree
and there are \(\ell^{a+b-i}\) factors.
For each permitted profile on a fixed OPR1
slice, distinct fields of this type of absolute
discriminant at most \(Y\) number
\(\gg Y^{2/(m(\ell^{a+b}-1))}/(\log Y)^C\).

The Frobenius at \(P,Q\) acts on \(G_T\)
by translations by the images of the
basis. Thus the local degrees are
their orders. Count surjections
\(V\to G=C_{\ell^a}\times C_{\ell^b}\)
with the prescribed image orders,
and divide by \(|\operatorname{Aut}G|\)
to count kernels. When \(a=b\) this
requires both orders \(\ell^a\), and
all such surjections have kernel
\(\ell^a V\). For \(a>b\), orient the
profile so that \(i=a,\ j=c\le a\).
The first image has unit first
coordinate, with
\(\varphi(\ell^a)\ell^b\) choices.
If \(c<a\), the first coordinate
of the second image is divisible
by \(\ell\); its second coordinate
must generate \(C_{\ell^b}\). If
\(c=b\), the first coordinate can
be any of the \(\ell^b\) elements
of order at most \(\ell^b\).
If \(b<c<a\), it must have order
\(\ell^c\), giving
\(\varphi(\ell^c)\) choices.
If \(c=a\), the first coordinate
is a unit and the second coordinate
must avoid one residue modulo
\(\ell\), giving the same number
\(\varphi(\ell^b)\) of choices.
For \(b=0\) the second coordinate
is absent and contributes one
choice. Orders \(c<b\) cannot
generate with the first image.
The number for the profile \(a,b\)
before dividing is exactly
\(|\operatorname{Aut}G|\), since
the two images then define an
automorphism of \(G\). This
gives the asserted ratios.

For example, degree \(m\) with
cyclic translation has
\(m+m/\ell\) fields of the same
discriminant in one \(W\).
At least one actual moving prime
is inert in each field. Counting
all degree \(m\) fields here
requires specifying the translation
type: degree alone also permits
noncyclic quotients with smaller
local degrees at both primes. The
count by \(Y\) follows from OPR1 and
(1), transporting its exact type
and signature on each fixed slice.

## The exact fixed-resolvent mass

**Theorem AFM3.** Let \(K\) be any nontrivial
quadratic field with \(\ell\) split and
with odd fundamental discriminant \(D\).
Put \(m=\ell^n\) and form \(E\) as in OPR2.
Let \(S_m(K)\) be the \(m\)-virtual classes
in \(K^\times/K^{\times m}\) whose two
local classes above \(\ell\) are \(m\)-th
powers. For a prescribed type
\(G=C_{\ell^a}\times C_{\ell^b}\) with
\(n\ge a\ge b\ge0\), \(a\ge1\),
the everywhere unramified extensions
\(W'/E\) of group \(G\) that split above
\(\ell\), are Galois over \(\mathbb Q\),
and have scalar \(U_m\)-action on \(G\),
are in bijection with the subgroups
of type \(G\) in \(S_m(K)\).
Each extension supplies one
\(\mathbb Q\)-isomorphism class of the
fixed fields of its complements, with
discriminant and signature (1).
Writing \(s_j=|S_m(K)[\ell^j]|\) and
\(s_0=1\), the number is
\[
  \begin{cases}
     (s_a-s_{a-1})/\varphi(\ell^a),& b=0,\\[2mm]
     (s_a-s_{a-1})(s_b-\ell s_{b-1})/
               |\operatorname{Aut}G|,& b\ge1.
  \end{cases}
 \tag{2}
\]
For fixed \(\ell,a,b\) and a class of
signed resolvent discriminants, summing
(2) for \(|D|\le
(Y/\ell^{f(G)})^{2/(|G|-1)}\) counts
exactly the associated degree-\(|G|\)
fields of discriminant at most \(Y\).
This mass uses the local-power
virtual group at varying \(K\).

For the descent, writing \(C=\mathbb
Q(\zeta_m)\), any such \(W'\) gives a
Kummer \(G\)-extension \(W'C/KC\):
adjoining \(C\) to \(E\) has degree
two and \(W'/E\) has odd degree.
The Kummer subgroup is fixed by
\(\operatorname{Gal}(KC/K)=U_m\)
and conjugation on \(K\) acts by
inversion. Indeed the action on
the radicals is dual to the scalar
action on the translations, multiplied
by the cyclotomic action on roots of
unity. For \((\sigma_K^i,a)\) their
quotient is \(a/((-1)^ia)=(-1)^i\).
This is the higher-power analogue
of the Kummer eigensigns in
[Cohen–Rubinstein-Salzedo–Thorne, §3](https://thornef.github.io/ohno-nakagawa3.pdf).

The restriction map from \(K^\times/
K^{\times m}\) onto the invariants in
\((KC)^\times/(KC)^{\times m}\) is an
isomorphism even when \(n>1\).
First take invariants under the
cyclotomic involution \(-1\).
They come from the fixed field \(KC^+\)
by the norm, since the norm times
restriction is multiplication by two
on the invariant group. The fixed
field contains no nontrivial
\(m\)-th roots of unity. For an
invariant class there, conjugate
ratios under \(\operatorname{Gal}(KC^+/K)\)
have unique \(m\)-th roots. Those
roots form a cocycle. Hilbert's
Theorem 90 divides the representative
by an \(m\)-th power to make it
rational over \(K\). Injectivity
follows by the same uniqueness
at this step and by the degree-two
norm at the previous step.
The argument holds over the
completions at \(\ell\) as well:
the relevant cyclotomic group
there is the full \(U_m\) because
\(\ell\) splits in \(K\).
At all other places \(KC/K\)
is unramified, so unramifiedness of
the Kummer extension gives
\(m\)-divisible valuations in \(K\).
The local splitting of \(W'/E\)
at \(\ell\) gives local powers in
\(K\) by the injectivity just proved.
This puts the descended subgroup
in \(S_m(K)\).

Conversely any virtual class in
\(K\) is anti-invariant: its rational
norm has all valuations divisible
by \(m\), hence is an \(m\)-th power
since \(m\) is odd. A subgroup of
\(S_m(K)\) gives a Kummer extension
over \(KC\), unramified away
from \(\ell\) and split there.
The diagonal involution of \(KC\)
acts trivially by conjugation on
its odd translation group. Its
unique order-two lift is central
in the normal extension over
\(\mathbb Q\). To get the order-two
lift adjust the square of any lift
by the unique inverse square root
in the odd central kernel. Its
fixed field gives precisely the
required scalar extension over
\(E\). The prime-to-two descent
preserves the unramified and split
conditions. The Kummer construction
and restriction are inverse.

The complements give one class of
fields per extension by the
conjugacy argument in AFM1.
When \(a<n\), the action on \(G\)
factors through \(U_{\ell^a}\)
and the core of a complement
is the kernel of this scalar
action. Thus the normal closure
is the fixed field of that core,
and recovers its maximal abelian
subfield \(E_{\ell^a}\), the
mirror at exponent \(\ell^a\).
From \(E_{\ell^a}\) one recovers
\(K\) by adjoining \(\zeta_{\ell^a}\)
and using split \(\ell\), as in
OPR2. Within fixed \(K\), the normal
closure also recovers \(W'\)
by compositing with \(E\): the
core lies inside a complement,
with zero intersection with the
translation group, and the two
fixed fields generate \(W'\).
Hence distinct extensions and
distinct \(K\)'s count separately.

For (2), choose an element \(x\)
of exact order \(\ell^a\), which
has \(s_a-s_{a-1}\) possibilities.
For \(b\ge1\), the second element
\(y\) must lie in \(S_m(K)[\ell^b]\)
with \(\ell^{b-1}y\) outside the
line generated by \(\ell^{a-1}x\).
This ensures exact order and
independence of the two cyclic
subgroups. The line is in the
image of multiplication by
\(\ell^{b-1}\) on \(S_m(K)[\ell^b]\).
Each of its \(\ell\) elements
has \(s_{b-1}\) preimages.
Divide the count of injections
by the automorphism order of \(G\).
This gives (2) and its summation
with (1). For imaginary \(K\),
\(S_m(K)\) embeds in
\(\operatorname{Cl}(K)[m]\),
since the units are \(m\)-divisible.
Its count is bounded by the
corresponding class group subgroup
count. In particular for \(a=b=1\)
it is controlled by the ordinary
\(\ell\)-torsion second moment.
For higher exponents, using only
\(\operatorname{Cl}(K)[\ell]\)
bounds each summand by a constant
times \(|\operatorname{Cl}(K)[\ell]|^{a+b}\).
An upper estimate for the ordinary
second moment alone does not supply
the needed higher moments.

## An equivalence for split quadratic torsion moments

**Theorem AFM4.** Fix \(\ell\) and
\(1\le a\le n,\ 0\le b\le a\).
For either sign let \(\mathcal K(X)\)
be the quadratic fields of that
sign with odd fundamental discriminant
of absolute value at most \(X\) and
split \(\ell\). Let \(N_{e,f}(X)\) count
the affine fields of AFM3 with
translation type
\(C_{\ell^e}\times C_{\ell^f}\)
and with resolvent in \(\mathcal K(X)\).
Then
\[
 \sum_{K\in\mathcal K(X)}
       |\operatorname{Cl}(K)[\ell^a]|\,
       |\operatorname{Cl}(K)[\ell^b]|
 \ \asymp_{\ell,n}\
 |\mathcal K(X)|+
   \sum_{\substack{1\le e\le a\\0\le f\le\min(e,b)}} N_{e,f}(X).
 \tag{3}
\]
The resolvent bound for each type is
equivalent to a discriminant bound
\(\ell^{f(G)}X^{(|G|-1)/2}\) as in
(1). In particular, for \(a=b=1\),
\[
 \sum_{K\in\mathcal K(X)}
           |\operatorname{Cl}(K)[\ell]|^2
       \ \asymp_\ell\
        |\mathcal K(X)|+N_{1,1}(X).
 \tag{4}
\]
For \(\ell=3\), the split odd
quadratic second moment with exponent
\(1+\epsilon\) is thus equivalent to
the natural \(Y^{1/4+\epsilon}\)
upper count for this split-at-three
degree-nine affine family, for
either sign separately.

Let \(V_m(K)\) be the full group of
virtual classes. The map to the
valuation-zero local quotient at one
prime above \(\ell\) has kernel
\(S_m(K)\). The other local class
is its inverse. The local quotient
consists of units modulo \(m\)-th
powers in \(\mathbb Q_\ell\) and
has order \(m\). Thus
\([V_m(K):S_m(K)]\le m\).
The virtual exact sequence to
\(\operatorname{Cl}(K)[m]\) has
kernel the units modulo powers,
of size at most \(m\) on both
signs. On taking \(\ell^j\)-torsion,
the kernel and cokernel of the
resulting map both have size at
most \(m\), the latter because
it injects into the quotient of
the unit kernel by \(\ell^j\).
These facts give uniform
comparability of the sizes of
\(S_m(K)[\ell^j]\) and
\(\operatorname{Cl}(K)[\ell^j]\)
for fixed \(m\).

Now count ordered pairs in
\(S_m(K)[\ell^a]\times S_m(K)[\ell^b]\)
according to the subgroup they
generate. The zero subgroup gives
one pair. Every other generated
subgroup has type \(C_{\ell^e}
\times C_{\ell^f}\), with
\(e\le a\) and \(f\le b\).
For the latter inequality, the
image of the first generator in
the quotient by the second is
cyclic, so the second generator
must have order at least
\(\ell^f\). Equivalently, if
\(f>b\), every element killed
by \(\ell^b\) is zero in the
quotient by \(\ell\) times the
whole subgroup and cannot serve
as one of its two generators.
Conversely any subgroup of a
listed type admits a generating
pair obeying the two bounds, by
taking its two invariant-factor
generators. Each subgroup
contributes between one and
\(\ell^{4n}\) ordered pairs, a
constant independent of \(K\).
Use (2) to sum their numbers.
This proves (3).

For the rank-one terms in the
special case \(a=b=1\), when
\(\dim S_\ell(K)\ge2\) the
number of one-dimensional
subspaces is bounded by a
constant times the number of
two-dimensional subspaces.
When the dimension is at most
one it is bounded by one.
This gives (4). For
\(\ell=3\) the exponent of
\(|D|\) for type \((1,1)\)
is four, yielding the stated
change of discriminant cutoff.
The equivalence concerns the
split odd slice; extending the
upper bound to other dyadic and
three-adic types requires
their corresponding local counts
or an independent comparison.

## The full prime-level second moment

**Theorem AFM5.** Fix an odd prime \(\ell\) and
a sign. For each quadratic \(K\) of that sign
outside the unique quadratic cyclotomic
exception, form the prime-level mirror
\(E=(K\mathbb Q(\zeta_\ell))^{(\sigma_K,-1)}\).
Count degree-\(\ell^2\) affine fields from
extensions \(W'/E\) with translation
group \(C_\ell^2\), scalar action,
everywhere unramified and completely
split above \(\ell\), as in AFM3.
There is no local restriction on \(K\)
at two or at \(\ell\).
For \(D=\operatorname{disc}K\) put
\(e_\ell=e(E_\ell/\mathbb Q_\ell)\).
The affine field has absolute
discriminant
\[
  \ell^{(\ell^2-1)(1-1/e_\ell)}
          (|D|/\ell^{v_\ell(D)})^{(\ell^2-1)/2},
 \qquad
 e_\ell=
 \begin{cases}
  \ell-1 & \ell\nmid D\ \text{or }\ell\equiv1\bmod4,\\
  (\ell-1)/2 & \ell\mid D,\ \ell\equiv3\bmod4 .
 \end{cases}
 \tag{5}
\]
It determines \(K\) by its normal
closure and scalar action.

Write \(N_\ell(Y)\) for the number of
these fields of discriminant at most
\(Y\), on the chosen quadratic sign.
The full quadratic \(\ell\)-torsion
second-moment bound \(O_{\ell,\epsilon}
(X^{1+\epsilon})\) on that sign is
equivalent to
\[
          N_\ell(Y)
             =O_{\ell,\epsilon}
                      (Y^{2/(\ell^2-1)+\epsilon})
\]
for every \(\epsilon>0\).
For \(\ell=3\), this gives the
degree-nine upper count at exponent
\(1/4+\epsilon\) as an equivalent
form of the full second moment.
The retained
[prime-path second moment](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md)
gives \(N_3(Y)\ll_\epsilon
Y^{17/56+\epsilon}\) from this
equivalence.

At prime level the Kummer descent
of AFM3 has prime-to-\(\ell\)
degree, globally and locally:
\([K(\zeta_\ell):K]=\ell-1\) and local
degrees divide \(\ell-1\). Thus its
bijection with two-dimensional
subgroups of \(S_\ell(K)\) still
holds. At both local places above
\(\ell\), or at the single place,
the condition defining \(S_\ell(K)\)
has uniformly bounded index in
\(V_\ell(K)\). Indeed its image in
the product of local unit
quotients modulo \(\ell\)-th powers
is bounded by
the corresponding orders for the
finite set of degree-one and
quadratic extensions of \(\mathbb
Q_\ell\). Taking one split unit
coordinate is sufficient in the
split case. The global unit
kernel has order at most \(\ell\)
on these quadratic fields.
Consequently
\[
           |S_\ell(K)|\asymp_\ell
                       |\operatorname{Cl}(K)[\ell]|
\]
on both signs, with constants
independent of \(D\).
Equation (2) for \(a=b=1\), with
the bounded ranks zero and one
absorbed by the number of quadratic
fields, gives the second-moment
comparison as in AFM4.

For (5), \(e_\ell\) is given by
[Cohen–Rubinstein-Salzedo–Thorne,
Lemma 2.11](https://thornef.github.io/ohno-nakagawa3.pdf).
The relative extension splits at
\(\ell\), so the inertia action
on the translation set is a tame
scalar subgroup of order \(e_\ell\):
any lift of an order-prime-to-\(\ell\)
subgroup of \(U_\ell\) is conjugate
to the scalar subgroup by averaging
the translation cocycle. Each
nonzero translation point has
orbit size \(e_\ell\), giving the
tame conductor
\((\ell^2-1)(1-1/e_\ell)\).
For a rational prime \(p\ne\ell\)
dividing \(D\), the inertia action
is a translate of inversion.
The nontrivial quadratic
character on inertia occurs
\((\ell^2-1)/2\) times in the
permutation representation.
Since \(\mathbb Q(\zeta_\ell)\)
is unramified at \(p\), that
quadratic character has the same
local conductor as the character
of \(K_p/\mathbb Q_p\), including
the wild dyadic cases. The relative
extension is unramified. This
gives the second factor in (5)
also when \(D\) is even.

The normal closure is \(W'\) since
the scalar action is faithful.
Its commutator subgroup is the
translation group, so one recovers
\(E\). Compositing with the
cyclotomic field recovers
\(K\mathbb Q(\zeta_\ell)\).
There are at most two possible
quadratic fields \(K,K'\) outside
the cyclotomic field with that
compositum. They differ by the
quadratic cyclotomic character.
If \(\ell\equiv3\bmod4\), their
diagonal involutions are different
and their mirrors differ. If
\(\ell\equiv1\bmod4\), their
mirrors coincide, but the second
choice changes the scalar
character \(u\) on \(U_\ell\)
to \(\chi_\ell(u)u\), where
\(\chi_\ell\) is the quadratic
character of \(\mathbb F_\ell^\times\).
The actions disagree on nonsquares
and cannot be conjugate by an
automorphism of \(C_\ell^2\),
since scalar actions commute with
all such automorphisms. Thus the
affine field determines the
resolvent in all cases.

Formula (5) changes the cutoff
\(X\) to a fixed multiple of
\(X^{(\ell^2-1)/2}\), with finitely
many local profiles and constants
depending only on \(\ell\).
The number of quadratic fields
contributes \(O(X)\). The one
excluded resolvent contributes a
fixed term. This proves both
directions of the upper-count
equivalence.
