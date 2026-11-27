# Mixed three-power class groups from two prime fibers

For every fixed \(m=3^n\), \(n\ge1\), a cubic marking
on one elliptic curve survives beside two independent
order-\(m\) boundary markings. The same local centers
and character rows work for all \(n\).
In imaginary specializations the three classes give
\((\mathbb Z/m)^2\times\mathbb Z/3\) in the ideal class
group. The first two factors can be represented by
the three-primary parts of classes over two actual
moving rational primes. In real specializations the
unit loss still leaves \(\mathbb Z/m\times\mathbb Z/3\).
The fields of either sign yield positive-rank twists
of one fixed elliptic curve, with exact conductor
\(27161D^2\).

## An elliptic curve with three compatible marks

**Proposition MP1.** Let
\[
 b=\frac{813074625}{2768896},\qquad c=-1664,\qquad
 E:y^2+(x+b)y=x^2(x+b).
\]
It has a point \(T=(0,0)\) of order four, with
\(S=2T=(-b,0)\), and the point
\[
 R=\left(-\frac{4492125}{173056},
         -\frac{166770140625}{287965184}\right)
\]
of order three. For \(t=y/[x(x+b)]\) the involution
of the quadratic parameter map is
\(\iota(P)=T-P\). Put
\[
 F(t)=(bt^2+t-1)^2-4bt^3,\qquad
 L=\frac{y^2}{x+b},\quad
 \gamma_0=\frac{L^2}{b^3},\quad
 \gamma_\infty(P)=\gamma_0(P+S).
\]
Let \(l_R\) be the tangent line function at \(R\)
with \(y\)-coefficient one, and set
\(\eta(P)=l_R(P+S)/l_R(\iota(P)+S)\).
The three functions have norm one to \(\mathbb Q(t)\).
On \(C_m:W^2=F(cu^m)\), their pullbacks give rational
Jacobian torsion \((\mathbb Z/m)^2\times\mathbb Z/3\)
via the divided divisors for \(m,m,3\), respectively.

The tangent at \(T\) is \(y=0\), meeting \(E\) in
\(2T+S\). The divisor of \(x+b\) is \(2S-2O\).
Thus
\[
  \operatorname{div}(t)=O+T-S-(-T),\qquad
  \operatorname{div}(L)=4(T-O).
\]
The two degree-two fibers have sum \(T\), giving
\(\iota\). Substitution gives
\(t^2x^2+(bt^2+t-1)x+bt=0\), so the quadratic field
of the map is defined by \(F\). The two roots have
product \(b/t\), and after adding \(b\) their product
is \(b/t^2\). Since \(L=t^2x^2(x+b)\), its norm is
\(b^3\). Translation by \(S\) gives \(t\mapsto1/(bt)\)
and commutes with \(\iota\). Thus the two norm-one
functions have respective divisors \(8(T-O)\) and
\(8((-T)-S)\). The tangent at \(R\) is a flex line
of divisor \(3R-3O\). Consequently \(\eta\) has
divisor three times a divisor of class \(2R\), since
\(\iota\) acts as minus one on \(\operatorname{Pic}^0(E)\).

The polynomial \(F\) is squarefree and has \(F(0)=1\).
Over \(\overline{\mathbb Q}\) adjoining \(u\) with
\(t=cu^m\) gives a cyclic extension of degree \(m\),
totally ramified over the four displayed boundary
points. The hyperelliptic curve \(C_m\) has genus
\(2m-1\). If a product of powers of the three marks
becomes a cube there, its class on \(E\) modulo cubes
is in the line of \(t\), by the unique cubic
subextension. At the two points above \(t=0\), the
orders of \(\gamma_0\) are opposite and nonzero
modulo three, whereas the orders of \(t\) agree.
The orders of \(\eta\) are zero modulo three.
So the exponent of \(\gamma_0\) must be zero
modulo three. The poles give the same conclusion for
\(\gamma_\infty\). Finally \(\eta\) has cube-divisible
divisor, so it cannot belong to a nonzero power of
the line of \(t\), and its divided divisor of order
three shows that it is itself noncubic. This proves
geometric independence modulo cubes.

The divided divisors upstairs exist by ramification
for the two boundary marks and by the flex identity
for \(\eta\). A principal relation on them would make
a product of \(\gamma_0,\gamma_\infty,\eta^{m/3}\)
an \(m\)-th power over the algebraic closure. In a
putative nontrivial relation divide the exponents
by their common maximal power of three. The remaining
relation modulo cubes contradicts the independence.
This proves the asserted mixed torsion group.

## The same finite places work for every exponent

**Proposition MP2.** For every \(m=3^n\), there is a
nonempty open slice at finitely many finite places
for \(u\), with \(u\) a unit at each of them, on which
the three marks in MP1 are ordinary virtual elements
of exponents \(m,m,3\) in
\(K=\mathbb Q(\sqrt{F(cu^m)})\). Their reductions
modulo cubes are independent. The primes
\(2,3,5,11,13,181\) ramify in \(K\), and \(157,173\)
split.

At the ramified primes use the following exact
centers in the \(t\)-coordinate. The third column
is the valuation of \(F(t_p)\).
\[
\begin{array}{c|r|r|r}
p&t_p&v_p(F(t_p))&
 \text{residue of the unit }t_p/c\\ \hline
2&30336&9&1\\
3&-1664&3&1\\
5&1&3&1\\
11&1&3&7\\
13&10686&1&8\\
181&1&1&31
\end{array}
\]
At \(3\) the ratio is exactly one. At the other primes
the residue of the ratio has order prime to three,
so the ratio is a \(3^n\)-th power of a local unit for
every \(n\). In each ramified quadratic completion
every norm-one element is a unit, giving all three
valuation conditions at once.

At \(157,173\), take \(t=cu_0^3\) modulo the prime,
with \(u_0=1,2\), respectively. Use the square roots
\(W=46,37\) modulo those primes. The reduced
\((\gamma_0,\gamma_\infty,\eta)\) values are
\((75,39,141)\) and \((149,10,7)\).
All are units on either branch by the norm-one
identities. Each residue of \(u_0^3\) belongs to the
image of the \(3^n\)-power map for all \(n\), since
\(v_3(p-1)\le1\). These tests give the two split primes.

For the independence tests use the following rows.
Here \(t=cu_0^3\bmod q\), \(W^2=F(t)\), and the exponents
are taken relative to \(\zeta\) after raising each mark
to \((q-1)/3\).
\[
 \begin{array}{c|rrr|c}
 q&u_0&W&\zeta&\text{exponents}\\ \hline
 7&3&1&2&(0,2,2)\\
 31&1&9&25&(2,0,2)\\
 43&4&11&36&(1,1,1)
 \end{array}
\]
The rows are independent, and in all three cases
\(v_3(q-1)=1\). Again the same residues of \(t/c\)
can be reached for every exponent. These checks
and the ramified centers are replayed by the
[exact checker](verify_mixed_three_power_boundary_torsion.py).
For clarity, \(x=(1-t-bt^2+W)/(2t^2)\) is used in
the evaluations.

At every other prime, \(\eta\) has cube-divisible
valuations. Indeed \(E\) has good reduction away
from the eight listed bad primes. Its involution
extends on the smooth proper model, with irreducible
special fiber. The norm-one identity forces the
vertical coefficient of \(\operatorname{div}\eta\)
to vanish there. The horizontal coefficients are
multiples of three, which gives the valuation
statement on pulling back to a valuation ring.
For the boundary marks, \(b\) is a unit outside
the six ramified primes. If \(t\) is also a unit,
both \(x\)'s are integral in the quadratic equation
above. Then both values of \(L=t^2x^2(x+b)\) are
integral, with unit norm. The same holds for its
translate since \(1/(bt)\) is a unit.

When \(v_p(t)=k>0\) at such a prime, the completion
splits, and the \(x\)-valuations are \(k\) and \(-2k\).
Thus \(\gamma_0\) has valuations \(8k,-8k\).
The translate is
\[
          L(P+S)=b^2(1+tx)^2/(x+b),
\]
a unit on both branches from the same two valuations.
For \(v_p(t)<0\) translation exchanges the roles of
the two functions. Since \(c=-2^7\cdot13\) and these
primes ramify on the slice, every remaining \(v_p(t)\)
is divisible by \(m\). This proves all valuation
assertions. Weak approximation gives the simultaneous
open slice.

## Two high-order split-prime classes and an extra cubic class

**Theorem MP3.** Fix \(m=3^n\), \(n\ge1\).
There are
\[
            \gg_m X^{1/(2m)}/(\log X)^2
\]
distinct fields of each sign with absolute
fundamental discriminant at most \(X\) from the slice
of MP2. The imaginary class groups contain
\((\mathbb Z/m)^2\times\mathbb Z/3\), and the real
class groups contain \(\mathbb Z/m\times\mathbb Z/3\).

There are also \(\gg_m X^{1/(2m)}/(\log X)^{B_m}\)
distinct fields of each sign on this slice with
\(u=\pm P/Q\), for two distinct positive rational
primes at scale \(X^{1/(4m)}\). Both primes split in
\(K\). In the imaginary case the three-primary
classes of a prime above \(P\) and one above \(Q\)
have exact order \(m\) and give the two order-\(m\)
factors, independent of the divided ideal of
\(\eta\). In the real case at least one of those
two prime classes has three-primary order \(m\).

For either family the three classes
\(\gamma_0,\gamma_\infty,\eta^{m/3}\) generate
\((\mathbb Z/m)^2\times\mathbb Z/3\) modulo
\(m\)-th powers in \(K\). In a putative nonzero
relation, divide the exponents by their common
highest power of three. This leaves a cube relation
because \(K\) has no nontrivial three-power roots
of unity, by the fixed ramification conditions.
The character rows give a contradiction.
The virtual-element exact sequence maps the subgroup
to \(\operatorname{Cl}(K)[m]\) with kernel contained
in the units modulo \(m\)-th powers. That kernel is
trivial in the imaginary case and at most cyclic
in the real case. Quotienting this subgroup by a
cyclic subgroup leaves rank at least two modulo
three. For \(m>3\), the image of its subgroup
\(3^{n-1}(\mathbb Z/m)^2\) still has positive rank,
so an element of exact order \(m\) survives.
This gives the real mixed subgroup as well.

For the unrestricted count, take a real open
interval near \(t=1/17\) or \(t=-1/17\), where
\(F(t)<0\) or \(F(t)>0\), respectively.
Such intervals give real open sets in \(u\).
Choose a rational center by weak approximation
and translate the parameter in
[Kulkarni–Levin, Theorem 1.2](sources/kulkarni-levin-hilbert-irreducibility-2111.15582.pdf).
The fixed neighborhoods allow a sufficiently small
ball at every specified place. The degree-two
map is the one on the smooth curve \(C_m\), and
the finitely many undefined evaluation parameters
may be deleted. That theorem gives
\(\gg_m H^2/(\log H)^2\) distinct quadratic fields
with parameters of height at most \(H\).
Since the defining binary form has degree \(4m\),
its squarefree part has absolute size \(O_m(H^{4m})\),
proving the discriminant lower bound.

For the prime count, every local condition can be
expressed by \(u\) in unit congruence classes modulo
fixed powers of the relevant primes. For each real
sign fix a rectangle with \(P,Q\asymp H\) and
\(P<Q\), with \(\pm P/Q\) in the chosen real interval
in \(u\). Impose \(Q\equiv1\) and the appropriate
unit residue for \(P\) modulo the fixed product.
The prime number theorem for arithmetic progressions
gives \(\gg_m H^2/(\log H)^2\) pairs there.

Clear denominators by a fixed square and take the
primitive part of the degree-\(4m\) binary form.
The prime-factor truncation in
[PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#quantitative-count-of-genuine-fields)
applies with \(4m\) as a bound for its projective roots.
In detail, for \(k=A_m\log\log H\) and
\(Q_0=H^{1/(2k)}\), the number of prime factors
larger than \(Q_0\) is \(O_m(k)\).
The number of pairs with at least \(k\) distinct
prime factors at most \(Q_0\) is
\[
    O_m\left(H^2\frac{(4m\log\log H+O_m(1))^k}{k!}\right),
\]
using the projective root bound and moduli at most
\(\sqrt H\). Taking \(A_m\) sufficiently large makes
this \(O_m(H^2/(\log H)^4)\). Thus most of the prime
pairs have \(\omega(D)=O_m(\log\log H)\).

The proof of PP2 applies to the fixed degree-\(4m\)
form with genus \(2m-1\). Its fixed splitting field
makes all the two-torsion of every quadratic twist
Jacobian rational. Kummer theory and good reduction
away from a fixed set and the primes dividing \(D\)
give rank \(O_m(1+\omega(D))\). The Abel–Jacobi
embedding has fixed degree \(2m-1\) for the theta
polarization. By Gao–Ge–Kühne, Theorem 1.1, the number
of rational points is at most \(C_m^{1+\omega(D)}\).
A genus at least two curve in its Jacobian contains
no positive-dimensional abelian translate.
This bounds the number of kept parameters over a
fixed field by a power of \(\log H\), giving the
prime-pair count at \(H\asymp_m X^{1/(4m)}\).

The only nonzero valuations of \(\gamma_0\) now
occur at the two places above \(P\), where they
are \(8m,-8m\). Those of \(\gamma_\infty\) occur at
the two places above \(Q\) with the same values.
Their divided ideal classes are the \(16\)-th
powers of the corresponding oriented prime classes,
because the conjugate classes are inverses.
Since \(16\) is prime to \(3\), this proves the
claimed prime descriptions. In the real case a
cyclic unit kernel cannot kill both independent
order-three classes obtained by multiplying the
two order-\(m\) directions by \(m/3\).

Kulkarni–Levin, Corollary 1.4 and Theorem 1.5,
give faster field counts for two unrestricted
order-\(m\) directions and for elementary imaginary
three-rank at least three. MP3 concerns their
combination at mixed exponents and identifies the
two explicit moving primes. No historical priority
claim is used.

## The elliptic twist and its conductor

**Corollary MP4.** Apart from finitely many parameters,
the fields in MP3 give a positive-rank twist \(E^D\)
of conductor \(27161D^2\), for both signs.
The distinct conductor counts inherit the field
count exponent \(1/(4m)\) and the same respective
logarithmic losses with conductor ordered by size.

The point \(2P-T\), for the point \(P\) on \(E\)
over the quadratic parameter \(t\), is negated
under \(\iota\). It therefore gives a rational point
on \(E^D\). Only finitely many points of \(E\) of
degree at most two can be torsion, by the bounded
height of elliptic torsion and Northcott's theorem.
The nonconstant map \(P\mapsto2P-T\) thus permits
only finitely many exceptional parameters.

The base conductor is
\(2\cdot3\cdot5\cdot11\cdot13\cdot157\cdot173\cdot181\).
At \(157,173\) the twist is locally trivial by the
split tests of MP2. At the six fixed ramified primes,
the local twist models and conductor exponents are
\[
\begin{array}{c|rrrrrr}
 p&2&3&5&11&13&181\\ \hline
 \text{local twist}&6&6&10&22&13&181\\
 v_p(N(E^D))&6&2&2&2&2&2.
\end{array}
\]
The exact checker verifies the squareclasses and
conductor exponents. At other ramified primes the
original curve has good reduction, which gives
conductor exponent two. These exponents are
\(2v_p(D)\) away from \(157,173\), and both these
primes have exponent one in the conductor. Their
product is \(27161\).

## The two-primary relation between the actual prime classes

**Proposition MP5.** Take the prime-ratio fields of
MP3 and orient \(\mathfrak P,\mathfrak Q\) by the
zeros of \(L,L\circ(+S)\) at the respective boundary
fibers. Let \(r\) be the class of the product of the
four ramified primes over \(3,5,11,181\).
In additive class group notation,
\[
 8m[\mathfrak P]=8m[\mathfrak Q]=r,\qquad
 2m([\mathfrak P]-[\mathfrak Q])=r .
 \tag{5}
\]
If the field is imaginary and \(|D|>4\cdot29865\),
the two classes individually have exact order
\(16m\), and together generate
\(\mathbb Z/(16m)\times\mathbb Z/(2m)\).
The inherited cubic class gives a further independent
\(\mathbb Z/3\) factor. The prime-ratio field count
of MP3 holds with these conclusions.

The ratio of the two tangent functions is a fourth
power up to a rational constant:
\[
     \frac{L}{L\circ(+S)}
            = \frac{(y/x)^4}{b^2},\qquad
             \operatorname N(y/x)=b .
 \tag{6}
\]
For the first identity use
\(L\circ(+S)=b^2(1+tx)^2/(x+b)\) and
\(y(y+x+b)=x^2(x+b)\). For the norm use
\(\operatorname N(x)=b/t\) and
\(\operatorname N(y)=b^2/t\).
The same formulas also give
\(\operatorname N(L)=\operatorname N(L\circ(+S))=b^3\).
The exact checker verifies (6) in the quadratic
function field.

At the two moving primes the ideals of \(L\) and
\(L\circ(+S)\) have opposing valuations \(4m,-4m\)
on the indicated primes and their conjugates.
At the six fixed ramified primes their normalized
valuations are \(3v_p(b)\), by the norm \(b^3\).
At all other primes their valuations vanish by the
unit argument of MP2. The odd valuations of \(b\)
occur at \(3,5,11,181\) only. Each corresponding
ramified ideal class has order at most two.
This gives the first two equalities in (5).
Formula (6) gives \(m,-m\) for \(y/x\) above
\(P\) and \(-m,m\) above \(Q\) in the stated
orientations. At the ramified primes its valuations
are \(v_p(b)\). This gives the last equality in (5).

For an imaginary field, a nonrational integral
element has norm at least \(|D|/4\). The norm of
the ideal defining \(r\) is the squarefree number
\(3\cdot5\cdot11\cdot181=29865\).
For \(|D|>4\cdot29865\) it cannot be principal.
Thus \(r\) has exact order two, and the two-primary
parts \(P_2,Q_2\) of the prime classes have
exact order sixteen. By (5),
\[
             Q_2=5P_2+e,\qquad 2e=0 .
\]
Indeed \(m\) is odd and \(2(P_2-Q_2)=8P_2\).

To see that \(e\) is independent of \(P_2\), use
the unramified quadratic extension \(K(\sqrt{13})/K\).
The prime \(13\) occurs to exponent one in \(D\).
At \(13\) the two ramified quadratic local fields
have unramified or trivial compositum over \(K\),
and away from \(13\) the additional quadratic
field is unramified. Real places, if present,
split in the extension. Its Artin character on
the classes above \(P,Q\) is
\((P/13),(Q/13)\), respectively, by quadratic
reciprocity. The local center at \(13\) in MP2
has \(t/c\equiv8\bmod13\) with valuation zero.
Since \(u^m=t/c\) and \(m\) is odd, \(u\bmod13\)
is a nonsquare. The sign in \(u=\pm P/Q\) has
no effect since \(-1\) is a square modulo \(13\).
Thus the character on \(e\) is nontrivial, while
the character on \(8P_2\) is trivial. This proves
that the two-primary subgroup is
\(\mathbb Z/16\times\mathbb Z/2\).
Combining with the three-primary independence
from MP3 proves the result. Deleting the finitely
many fields with small \(|D|\) preserves the count.

Equivalently in the large imaginary fields, the
full principal-relation lattice on the two
oriented primes consists exactly of exponent
pairs \((a,b')\in\mathbb Z^2\) satisfying
\[
        m\mid a,\quad m\mid b',\quad
        a+5b'\equiv0\bmod16,\quad b'\equiv0\bmod2.
\]
For either quadratic sign, the genus character
at \(13\) separates the two prime classes for
every choice of orientations. In particular
the two actual prime supports cannot be put in
one coset modulo the odd \(m\)-torsion, despite
the independent three-primary markings in the
imaginary case.

## Infinite towers with prescribed moving-prime decomposition

**Lemma MP6 (cuts preserving the abelianization).** Let \(K\) be an
imaginary quadratic field other than \(\mathbb Q(\sqrt{-3})\), put
\(A=\operatorname{Cl}(K)_{(3)}\), \(d=\dim_{\mathbb F_3}A/3A\), and let
\(G=\operatorname{Gal}(K^{\mathrm{ur},3}/K)\) be its maximal unramified
pro-\(3\) group. For finitely many primes \(\mathfrak l_i\) of \(K\),
choose Frobenius representatives \(x_i\in G\) and powers
\(e_i=3^{a_i}\), \(a_i\ge0\), killing their images in \(A\). Set
\[
 b_i=\begin{cases}e_i&e_i>1,\\2&e_i=1.\end{cases}
\]
If, for some \(0<t<1\),
\[
                 1-dt+dt^3+\sum_i t^{b_i}<0,
                 \tag{7}
\]
then the quotient
\[
 \Gamma=G/\overline{\langle\!\langle x_i^{e_i}:i\rangle\!\rangle}
                 \tag{8}
\]
is infinite and its natural abelianization is still \(A\).
If the image of \(x_i\) in \(A\) has exact order \(e_i\), the
decomposition group at \(\mathfrak l_i\) in the corresponding extension
is exactly \(C_{e_i}\). In particular, \(e_i=1\) gives complete splitting.
If the marked prime/exponent pairs are stable under quadratic
conjugation, that extension is also Galois over \(\mathbb Q\).

Proof. The arithmetic input in
[Koch–Venkov, §1, pp. 57–58](sources/koch-venkov-imaginary-quadratic-p-towers-1975.pdf)
gives a minimal presentation \(G=F/R\) with \(d\) generators and
\(d\) defining relations. Their Satz 2 and proof, pp. 60–62,
place the relations in odd Zassenhaus depths, hence at least three.
Thus their contribution to the Golod–Shafarevich polynomial is at most
\(dt^3\) for \(0<t<1\).

For \(e_i>1\), any lift \(f_i\in F\) of \(x_i\) has
\(f_i^{3^{a_i}}\in F_{3^{a_i}}\). This also follows directly from
\((f_i-1)^{3^{a_i}}=f_i^{3^{a_i}}-1\) in the completed group algebra
over \(\mathbb F_3\). For \(e_i=1\), the class hypothesis puts
\(x_i\) in \(\overline{[G,G]}\subseteq G_2\), so it has a lift in
\(F_2\). These are precisely the power-depth and quotient-filtration
properties in [Hajir–Maire–Ramakrishna, Propositions 1.2(i), 1.6 and
1.9, pp. 4–5](sources/hajir-maire-ramakrishna-cutting-towers-1901.04354v1.pdf).
Adjoining these lifted relators leaves the generator rank unchanged
and bounds the new polynomial by (7). Their Theorem 1.4, p. 4,
therefore proves infinitude; no equality of relation ranks is asserted
for the quotient.

Class field theory identifies \(G^{\mathrm{ab}}\) with \(A\).
Every \(x_i^{e_i}\) is in \(\overline{[G,G]}\), so their closed normal
closure is contained there. Thus (8) preserves the whole
abelianization, a stronger assertion than preservation of generator
rank alone. In \(\Gamma\), the order of \(x_i\) divides \(e_i\) and
is at least its order in \(A\). Unramified decomposition groups are
procyclic and generated by Frobenius, proving the exact-order claim.
Finally quadratic conjugation takes Frobenius to a conjugate of the
Frobenius at the conjugate prime. Stability of the marked pairs
makes the closed normal subgroup in (8) invariant, proving the last
assertion. The field in question is a subextension of
\(K^{\mathrm{ur},3}\), so it remains everywhere unramified over \(K\).

**Corollary MP7 (two moving local degrees).** Fix \(m=3^n\ge9\).
Every imaginary prime-ratio field \(K\) of MP3 admits an infinite
everywhere unramified pro-\(3\) extension \(L/K\), Galois over
\(\mathbb Q\), such that
\[
 \operatorname{Gal}(L/K)^{\mathrm{ab}}\simeq\operatorname{Cl}(K)_{(3)},
 \qquad D_{\mathfrak l}(L/K)\simeq C_m
 \quad(\mathfrak l\mid PQ).
\]
Consequently there are
\[
       \gg_m X^{1/(2m)}/(\log X)^{B_m}
\]
distinct imaginary quadratic fields of absolute fundamental
discriminant at most \(X\) with these properties. For \(m=9\) the
exponent is \(1/18\) and all four marked local degrees are nine.

Proof. MP3 gives \(d\ge3\) and exact order \(m\) for the images in
\(A\) of a prime above each of the distinct split rational primes
\(P,Q\). Their conjugate classes are inverses and have the same
orders. Apply MP6 to all four primes with \(e_i=m\). At \(t=1/2\),
\[
 1-\frac d2+\frac d8+\frac4{2^m}
       =1-\frac{3d}{8}+2^{2-m}
       \le-\frac18+2^{2-m}
       \le-\frac{15}{128}<0.
\]
The four cuts are stable under conjugation. MP6 gives both infinitude
and exact Frobenius orders, using the unchanged full abelianization.
It applies to every field in the imaginary prime-ratio family, so
the distinct-field count transfers without another counting loss.
No assertion for \(m=3\) follows from this polynomial estimate.

The imported tower mechanism is classical; MP7 is its application to
the exact moving-prime markings of MP3. No historical novelty is
claimed. This is a written derived corollary with internal source
inspection, not independent review or a computation of the infinite
Galois group. The finite local checker above verifies the underlying
marked-field data, not the tower argument.

Verification scope (2026-11-27): the retained Koch–Venkov §1 and
Satz 2 proof, and the retained arXiv:1901.04354v1 §1.1–1.4, were
inspected for the presentation, power-depth, quotient and infinitude
inputs. The mixed-field and octic local checkers were internally
replayed with their existing \(\texttt{sage -python}\) commands.
Exact rational arithmetic confirms the two extremal polynomial
values \(-15/128\) and \(-3/125\) in MP7 and MO4a. These checks do
not supply independent mathematical review of the underlying
families or of these corollaries.
