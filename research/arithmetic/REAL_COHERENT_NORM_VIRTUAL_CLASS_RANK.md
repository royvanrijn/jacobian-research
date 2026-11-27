# Real coherent norms and the unit direction

In a real quadratic field an ideal class can disappear from
an edge label because of a unit. A different obstruction gives
control of the virtual-unit labels: the irrational coefficient
of an odd power of a nonrational integer has a lower bound.
For cubic triangles this combines with the Fermat line section
to force two independent virtual-unit classes. A squarefree
family realizes this with all three primes in one dyadic interval.

## Short coefficients distinguish virtual-unit labels

**Proposition RV1.** Let \(d>1\) be squarefree and real,
and \(\ell\) an odd prime. A nonrational quadratic integer
\(\alpha=(u+v\sqrt d)/2\) with
\[
                   0<|v|<2^{1-\ell}d^{(\ell-1)/2}
\]
is not an \(\ell\)-th power in \(\mathbb Q(\sqrt d)\).
In a real coherent prime part with all coefficients below this
bound, the relative virtual-unit labels at the vertices are
distinct. Each ordinary ideal-class label occurs at most
\(\ell\) times.

If \(\alpha=\gamma^\ell\), integrality of \(\alpha\) implies
integrality of \(\gamma\). Write
\(\gamma=(t+s\sqrt d)/2\) with \(t,s\) integers and \(s\ne0\).
Expansion gives
\[
       v=2^{1-\ell}
           \sum_{j=0}^{(\ell-1)/2}
               \binom{\ell}{2j+1}t^{\ell-2j-1}s^{2j+1}d^j .
\]
The power of \(t\) is even and the summands have the same
sign when nonzero. The last term alone has the required
lower bound. For the labels, use power classes modulo
\(\ell\)-th powers of \(\alpha_{0i}\) for a chosen anchor.
Coherence makes the difference of two labels the power
class of \(\alpha_{ij}\), because the rational norm factor
at the anchor is an \(\ell\)-th power. Thus the labels are
distinct. The map from virtual units modulo \(\ell\)-th powers
to \(\mathrm{Cl}(\mathbb Q(\sqrt d))[\ell]\) has kernel the
units modulo \(\ell\)-th powers, of size \(\ell\) in a real
quadratic field. This gives the multiplicity bound.

In the [bounded-unit-logarithm selection](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#coherent-witnesses-with-a-bounded-unit-logarithm),
the coefficients are \(O_\ell(Z^\ell/\sqrt d)\).
Thus whenever \(Z=o(\sqrt d)\), keeping one prime per
ordinary label loses at most a fixed factor in each part.
In parts of growing size the number of affine-collinear
triples of distinct retained labels is \(O_\ell(m^2)\),
since a fixed pair determines a line with \(\ell\) points.
One can therefore select independent edge-class differences
at the same order of magnitude as the original number of
coherent triangles. This is a real counterpart to the
nonrational norm lower bound used in the imaginary witness
selection.

## A same-colour triangle gives two independent virtual units

**Proposition RV2.** Let \(K=\mathbb Q(\sqrt d)\) be real
quadratic, and take a coherent cubic prime triangle at
distinct odd split primes. Suppose the three nonzero
coefficients have absolute value below \(d/4\) and have
the same cube class \(w\). The two consecutive virtual
units are independent modulo cubes. In particular
\(\mathrm{Cl}(K)[3]\ne0\), and the rank of
\(y^2=x^3+d w^2/4\) is at least two.

The three edge classes in the virtual-unit quotient are
nonzero by RV1. The sum of the two consecutive classes
is the class on the third edge, also nonzero.
For their difference, use the binary cubic (1) of
[RF1](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#the-exact-reflected-kummer-class)
with the signed parameter \(E=-d\delta^2<0\), where
\(\delta\) is the common coefficient gcd. The trace and
Cardano identities are polynomial identities and still
give the Kummer class \(\overline\alpha_{12}\alpha_{23}\)
over \(K(\zeta_3)\). The same cube class on the three
edges makes the normalized coefficients signed cubes.
Thus the associated line section is on the Fermat cubic,
and its three line coefficients have distinct prime
square classes. No rational Fermat point lies on the
line. A zero Kummer class would put the splitting field
inside the biquadratic field \(K(\zeta_3)\), which cannot
contain a root of an irreducible cubic. This excludes
equality of the two consecutive virtual-unit classes.
Together the sum and difference tests give independence
over \(\mathbb F_3\).

The unit kernel of the class map has dimension one, so
its image on this span is nonzero. For the rank assertion,
the evaluation in
[CR1](QUADRATIC_NORM_COLOUR_CLASS_RANK.md#from-one-colour-to-one-class-subspace)
already takes the points to the independent virtual-unit
classes before projecting to ideal classes.
The Jacobian has no rational 3-primary torsion since the
real quadratic field is disjoint from \(\mathbb Q(\zeta_3)\).
An integral dependence modulo torsion would then give a
dependence of those classes with a primitive coefficient
vector, which the independence excludes.

One exact example showing the unit contribution is
\(d=86957347778285=5\cdot17\cdot83\cdot107\cdot823\cdot139967\).
The primes \((577,457,47)\) have normalized coefficients
\((8,1,1)\) and doubled traces on edges \(12,23\) equal
to \(280899354,-11251421\). The signed Heron parameter
is \(-d\). The real quadratic class group has invariant
factors \(12,2,2,2\), verified by
`bnfinit(x^2-d,1)` and `bnfcertify` in PARI/GP.
Thus its three-torsion has dimension one, although
the two virtual-unit classes are independent.

## A dyadic squarefree family over real fields

**Theorem RV3.** For every sufficiently large \(X\), at least
\(\gg X^{1/3}/(\log X)^3\) distinct squarefree \(d\in[X,2X)\)
with \(d\equiv1\bmod4\) occur in genuine coherent real
prime triangles with all three primes in a common dyadic
interval at scale \(X^{1/6}\). The absolute edge coefficients
are \((10^3,1,481^3)\), all cubes. The triangles give rank
at least two on \(y^2=x^3+d/4\) and nontrivial
\(\mathrm{Cl}(\mathbb Q(\sqrt d))[3]\). Their reflected
cubic fields have negative fundamental discriminant \(-3d\).

Use \(x=10^3,y=1,z=481^3\), where \(481=13\cdot37\).
Write the scaled prime directions as
\[
   (P_1,P_2,P_3)=(p_1,481^2p_2,100p_3),\qquad
   P_1+P_2=\lambda P_3,\qquad \lambda=\frac{166667}{50}.
\]
Equivalently \(p_1+481^2p_2=333334p_3\).
The direction
\[
       (p_1,p_2,p_3)=
         \left(\frac65,\frac32,\frac{3470427}{3333340}\right)
\]
has strictly negative \(H_0(P_1,P_2,P_3)\), with \(H_0\)
the Heron form of SW1. All three coordinates are in
\((1,2)\). The line is a secant of the Fermat dual, from
the conjugate points with \(X_1+X_2=1,\ X_3=-1/100\).
Substituting
\((P_1,P_2,P_3)=(166667U,166667V,50(U+V))\) gives,
up to a nonzero scalar and on \(V=1\), the factors
\[
 \begin{split}
 &111110888889\,U^2+111110222222\,U+111110888889,\\
 &192903935191550904321\,U^4-192902777796296441358\,U^3\\
 &\qquad -2314846065037037\,U^2-192902777796296441358\,U
                      +192903935191550904321 .
 \end{split}
\]
The quadratic is irreducible over \(\mathbb Q\), and the
quartic is irreducible modulo \(23\). They are separable
and coprime.

For the full norm congruences at \(2\), take locally
\(p_3=1\), \(P_1=166667+4,\ P_2=166667-4\).
Then \(v_2(P_1^3-P_2^3)=3,\ v_2(P_3^3)=6\), so
\(E=H_0/(xyz)^2\equiv3\bmod4\).
At \(5\), take the scaled \(P_1=P_2=166667,\ p_3=1\);
then \(H_0=P_3^3(4P_1^3-P_3^3)\) and the normalized
\(E\) is a unit. At each \(q=13,37\),
\(v_q(\lambda^3-1)=2\). Let \(\zeta\) be the 3rd root
of unity in \(\mathbb Z_q\) reducing to
\(\lambda\bmod q\). Choose \(P_3=100,\ P_1=\zeta P_3,\
P_2=(\lambda-\zeta)P_3\). This gives primitive prime
residues, \(P_1^3=P_3^3\), and \(E\) a unit since
\(v_q(P_2)=2\). At \(3\), all weights are units and
\(P_1,P_2\) agree modulo \(3\) with \(P_3\) of opposite
residue, so \(E\equiv1\bmod3\).
All choices give open local sets in the integer relation.

Use the squarefree sieve on this relation exactly as in
[SW4](ABC_COHERENT_PRIME_TRIANGLES_SQUAREFREE_WINDOW.md#one-colour-with-all-three-primes-in-the-same-dyadic-interval).
Outside the weight support and \(2,3\), every remaining
prime is at least \(7\); for \(q\ge11\) there is a unit
direction avoiding the six roots and the three vanishing
prime forms. At \(7\) take \(p_3=1,\ P_1=P_2=166667\);
the three prime values are units and \(H_0\) is nonzero
modulo \(7\). Thus all finite bad primes can be treated
without keeping any square factor of \(E\).
The binary form has a real linear factor over \(\mathbb R\),
since it is negative at the displayed direction and
positive on the line at \(P_1=P_2\).
The quartic large-square estimate and the distinct-value
bound therefore give the stated \(\gg X^{1/3}/(\log X)^3\)
distinct \(d=-E\) in a chosen dyadic interval. The
incident-prime unit condition follows as in SW4 from
\(\lambda^3-1\) and the opposite cube residues at \(p_3\).

The exact reconstruction identities use the signed
\(E\) here. They yield integral real quadratic
generators, now with norm equation
\(4(p_i p_j)^3=u_{ij}^2-dv_{ij}^2\). The incident
prime units give their oriented prime-cube ideals.
The signs can be chosen totally positive: the sign
of each generator agrees in the two embeddings and
the signs multiply along a coherent path. None of
this changes the cube colour, so RV2 applies.

Finally the signed cubic reflection formula gives a
reflected splitting field with quadratic subfield
\(\mathbb Q(\sqrt{-3d})\), unramified over that field
outside \(3\) by its virtual-unit radicand.
At \(3\) the mixed signs of the scaled marks allow the
same Hensel lift as in RF3, giving a local rational root
and a quadratic splitting field. This proves that the
reflected discriminant is the fundamental value
\(-3d\). The examples also give rank growth on the fixed
Fermat elliptic curve over these reflected cubic fields
apart from finitely many fields with a bounded-degree
torsion line-section point, as in RF4.

An exact squarefree control in this dyadic family is
\((p_1,p_2,p_3)=(130330379471,160333237511,111284742913)\).
The value of \(-E\) is the product of the seven distinct primes
\[
 29,\ 31,\ 43,\ 44129,\ 976489,\ 26065050478471,\
 4845445037788068598690132208784197839929308433179 .
\]
Primality of all the factors and the reflected discriminant
\(-3d\) were checked with PARI/GP `isprime` and `nfdisc`
using the displayed relation and the form in RF1.
