# Splitting at the prime marks of a reflected cubic

The reflected cubic field of a coherent norm triangle also records
the decomposition of its three marked primes. In the prime
progression construction these moving splitting tests reduce to
Frobenius in two fixed pure cubic splitting fields. One can
prescribe them jointly with Frobenius at finitely many fixed
primes while retaining the distinct-field bound of the progression
construction.

## The moving primes

**Proposition PS1.** Let \(M\) be the reflected cubic field of an
irreducible triangle in [RF1](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#the-exact-reflected-kummer-class).
For any vertex \(p_i\ne3\), \(M\) is unramified at \(p_i\).
If \(p_i\equiv2\bmod3\), then its residue degrees are \(1,2\).
If \(p_i\equiv1\bmod3\), then it splits completely exactly when
the listed unit is a cube modulo \(p_i\), and is inert otherwise:
\[
 \begin{array}{c|c}
  \text{prime}&\text{unit}\\ \hline
  p_1 & u_{12}x/z\\
  p_2 & E xy\\
  p_3 & u_{23}y/z
 \end{array}                                                     \tag{1}
\]
All quotients in the table have unit denominators at their
designated primes. For a progression \(p_1=p-h,\ p_2=p,\
p_3=p+h\) with normalized coefficients \((1,1,1)\), the cubic
residue tests are those for \(7,2,7\), respectively.

Each vertex is split in \(K\) and coprime to \(3\Delta_K\) in
the stated case. The quadratic character of the reflected
field at that prime is therefore the character of \(-3\).
The unramified Frobenius is a transposition when \(p_i\equiv2\bmod3\).
Otherwise the prime splits in \(K(\zeta_3)\), and (2) of RF1
gives a split-local test by evaluating
\(\gamma=\overline\alpha_{12}\alpha_{23}\) at a prime above
\(p_i\) where it is a unit.

At \(\overline{\mathfrak p}_1\), \(\alpha_{12}=0\) modulo that
prime and \(\overline\alpha_{12}=u_{12}\). With
\[
 s=\frac{u_{12}v_{23}-u_{23}v_{12}}2=\frac{\delta(A-C)}{z}
\]
as in [RT1](ABC_COHERENT_REFLECTION_AT_THREE.md#the-residue-test),
one has \(\alpha_{23}=-s/v_{12}=(C-A)/(xz)\) at the same
prime. The residue of \(\gamma\) is \(p_3^3 u_{12}x/z\).
At \(\mathfrak p_2\), both \(\alpha_{12}\) and
\(\overline\alpha_{23}\) vanish, so the residue of \(\gamma\)
is \(u_{12}u_{23}=Exy\bmod p_2\) by coherence.
At \(\overline{\mathfrak p}_3\) the residue is
\(u_{23}s/v_{23}=p_1^3 u_{23}y/z\). This proves (1), including
when a vertex divides its opposite normalized coefficient.
In a progression, \(u_{12}\equiv-7h^3\bmod p_1\),
\(E\equiv-4h^6\bmod p_2\), and
\(u_{23}\equiv7h^3\bmod p_3\). Signs and squaring do not
alter the test for a cube.

For \((p_1,p_2,p_3)=(7,37,67)\) and \((x,y,z)=(2,7,5)\),
the cubic residue powers in (1) give \(1,26,1\) in the
respective prime fields. The reflected cubic field has
discriminant \(595009\) and is split, inert, split at the three
marks, by PARI/GP `idealprimedec`. Here the first vertex
divides the opposite normalized coefficient \(y\).

## Independent local choices on the progression family

**Theorem PS2.** Prescribe split or inert at each of the three
ordered prime marks. Also fix a finite set of rational primes
\(q>169\), with one of the three unramified cubic splitting
types prescribed at each \(q\). There is a constant \(B>0\)
such that at least \(\gg X^{1/3}/(\log X)^B\) distinct
totally real cubic fields of fundamental discriminant at most
\(X\) occur in [RF4](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#a-family-of-cubic-points-and-the-fixed-field-count)
with all these prescriptions. The marks are primes congruent
to \(1\bmod3\). The fields still have \(3\) completely split,
and their imaginary reflections have the two independent
three-torsion classes given by the progression. Constants may
depend on the finite prescriptions. Before passing to distinct
fields, the relative asymptotic frequency of each moving
profile in the constructed progression subfamily is the
product of \(1/3\) for each split mark and \(2/3\) for each
inert mark.

First consider the fixed-prime conditions. The pencil of lines
with coefficients \((p-h,p,p+h)\) on the Fermat cubic has
center \((1:-2:1)\), off the curve modulo \(q>3\). Projection
from the center gives a connected separable degree-three cover
of the \(p/h\)-line over \(\overline{\mathbb F}_q\).
The branch polynomial is \(F(t,1)\) of PP2, with six simple
roots for \(q>7\), and infinity is unramified because
\(F(1,0)=3\). To see that no extra branch comes from the choice
of coordinates in the line restriction, its cubic discriminant
on \(p_3\ne0\) is \(27p_3^6 F(p,h)\) by (4) of RF3. The charts
with \(p_1\) or \(p_2\) nonzero give the same square class and
cover the other directions. Each simple branch has a
transposition. The transitive geometric group is therefore
\(S_3\), as is the arithmetic group, and the Galois closure
is geometrically connected. It has genus four by
Riemann–Hurwitz: each of six simple branch places contributes
three in the degree-six closure.

For any desired permutation \(g\in S_3\), twist the Galois
cover by \(g\) as the image of finite-field Frobenius.
The twist has the same genus and its rational points over
unramified directions specify Frobenius in the conjugacy
class of \(g\): a point on the twist corresponds to a
geometric point on the cover on which Frobenius acts by
\(g\) in its deck orbit, with the inverse convention giving
the same conjugacy class. Hasse–Weil gives at least
\(q+1-8\sqrt q\) points on the twist. Excluding the six
branch values and \(t=0,1,-1,\infty\) removes at most
\(60\) points. The remaining count is positive for \(q>169\).
Choose one remaining \(t_q\in\mathbb F_q\) and impose
\(p\equiv t_q,\ h\equiv1\bmod q\).
Then the primes are units at \(q\), the cubic discriminant is
a unit there, and the prescribed factorization lifts over
\(\mathbb Q_q\).

Next impose independent-ideal tests as in PP1, available
away from any of the fixed primes. In the control progression
\((17,23,29)\) of PP1 let \(\alpha,\beta\) denote the first
two generators. Their cube classes are independent over
\(K_*=\mathbb Q(\sqrt{-185684403})\) by its two exact residue
tests, and remain independent over \(B_*=K_*(\zeta_3)\) by
the quadratic norm. Adjoining the class of \(\zeta_3\) gives
three independent cube classes over \(B_*\). Indeed
\([\zeta_3]\ne0\) by degree, and conjugation of \(\zeta_3\)
fixing \(K_*\) acts by inversion on this class and fixes
\(\alpha,\beta\). The Kummer extension is normal over
\(\mathbb Q\), since conjugating either generator of \(K_*\)
replaces its class by its inverse, by the rational cube
norm. Chebotarev in this fixed extension gives arbitrarily
large auxiliary rational primes with Frobenius trivial on
\(B_*(\zeta_9)\), and nontrivial on exactly one of
\(\sqrt[3]\alpha,\sqrt[3]\beta\). The two separate choices
give primes \(r_1,r_2\) outside the fixed set and outside
\(2,3,7,17,23,29\). The nontrivial tests are unchanged if the
chosen Frobenius is conjugated, since conjugation preserves
the respective Kummer lines.

Impose \((p,h)\equiv(23,6)\bmod r_1r_2\), together with \(p\)
odd and \(p\equiv1\bmod3,\ h\equiv0\bmod6\). The same leading
residues and square-root choices as in PP1 imply that the
two consecutive ideal classes are independent in every
valid specialization. The auxiliary primes are \(1\bmod9\),
so an imaginary quadratic unit can never upset their cube
tests. All these conditions and the fixed-prime conditions
are simultaneous by CRT. The positivity window of PP3 and
the \(p_i>7\) condition still give the coherent triangle, and
because the three \(p_i\) agree modulo \(3\), RF3 applies.

Write the resulting fixed congruences on \((p,h)\) as two
residues modulo a fixed modulus. The three prime forms then
have nonparallel linear parts, as in PP3, with no fixed
prime divisor. For primes dividing the modulus this follows
from the unit choices. At the other primes one may take
\(p=1,h=0\) locally. Apply
[Frei–Tinková, Theorem 1.4](https://arxiv.org/pdf/2408.17441)
to these forms and the PP3 rectangle, with Galois fields
\(K_7,K_2,K_7\) where \(K_a=\mathbb Q(\zeta_3,\sqrt[3]a)\).
The desired classes are the identity or the class of
three-cycles. Their restrictions to the common abelian
subfield \(\mathbb Q(\zeta_3)\) agree, and \(p_i\equiv1\bmod3\)
makes the density there positive. The Euler factors away
from the abelian conductor are the usual local densities.
Conditioning on the congruences, the factors accompanying
the identity and the pair of three-cycles are thus in ratio
\(1:2\) at each mark, and the main term has positive constant
for every profile. The theorem gives \(\gg Z^2/(\log Z)^3\)
progressions in each case, and its product formula gives
the stated frequencies. Prime-power contributions are
\(o(Z^2)\) in the weighted count, and the logs of the
prime forms are uniformly \(\log Z+O(1)\) in this rectangle.

Finally the same squarefree-part multiplicity bound PP2 and
exceptional \(\omega(d)\) estimate of PP3 apply to each
restricted subfamily, giving the asserted logarithmic-loss
distinct-field count. The fundamental reflected discriminant
is \(d/3\) by RF3, so distinct \(d\) give distinct cubic
fields. Delete the finitely many fields of degree at most
three carrying a torsion line-section point, as in RF4.
