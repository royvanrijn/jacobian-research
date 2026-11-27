# Small partial-degree gcds from standard planar corners

## Planar pairs

**Theorem PlanarSmallPartialGcdCorners.**
Let \(P,Q\) be a characteristic-zero planar Keller pair with positive
partial degrees \(au,bu\) in the same source variable, where \(a,b>1\)
are coprime. Let \(h\) be the common monic leading factor in the other
variable, so that the leading coefficients are scalar multiples of
\(h^a,h^b\). If the pair is noninvertible, then \(u\ge4\).
For the next three gcds the geometric multiplicities of \(h\) satisfy
the following restrictions.

| Partial-degree gcd | Geometric root multiplicities of \(h\) |
|---:|---|
| 4 | One root of multiplicity 12 |
| 5 | One root of multiplicity 20 |
| 6 | One root of multiplicity 30, or multiplicities 3 and 15 with at least one 15 |

In the last alternative, if there is more than one root, then \(6\mid\deg h\)
and the numbers \(s,r\) of roots of multiplicity three and fifteen satisfy
\(s\le3r\).
In particular, if \(u=6\) and \(\deg h=30\), the possible multiplicity
lists are \((30)\) and \((15,15)\).

Work over an algebraic closure. Extension of the ground field detects
invertibility by faithful flatness. The
[simultaneous triangular normalization](PARTIAL_DEGREE_FOURTEEN_AND_WEYL_TRANSFER.md#simultaneous-triangular-normalization)
preserves the two leading rows and gives total degrees \(a(u+d),b(u+d)\),
where \(d=\deg h\). For a noninvertible pair we have \(u+d\ge16\) by
Guccione–Guccione–Valqui, Corollary 6.6. Thus \(d>u\) when \(u\le6\).
The [axial root lemma](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#one-translation-for-the-two-leading-coefficients)
gives a root of \(h\) of multiplicity \(\mu>u\). Translate it to zero
and swap the axes. Translation preserves both total degrees.
The resulting pair is standard in GGV Definition 4.3 with parameters
\(a,b\). Indeed, the total-degree and right-face degree ratios are
both \(a:b\), and the normalized start and end of the right face
are \((u,\mu)\) and \((u,d)\). Both have negative \((1,-1)\)-weight.
The same standardization works at every root of multiplicity \(>u\).
[GGV, primary text](https://arxiv.org/pdf/1401.1784v3)

Use GGV's primitive corner \(A_0=(i,j)=a^{-1}\operatorname{en}_{\rho,\sigma}P\)
in these swapped axes. By GGV Proposition 6.1 and Definition 6.2,
the corner is of type II.b and \(\sigma<0\). Its first coordinate
cannot exceed \(u\), since it comes from the support of \(P\).
Propositions 6.4, 6.5 and 6.7 give
\[
  4\le i\le u,\qquad i<j\le i(i-1),\qquad i+j\ge16,\qquad
  \gcd(i,j)>1.
\]
There is a bracket mate with corner \((f_1,f_2)=(p_0/q_0)A_0\),
where \(p_0/q_0\) is reduced and \(2\le f_1<i\). Proposition 7.8(4)
excludes \(\gcd(i,j)=2\) and \(q_0=2\).

We use two consequences of the regular-corner sequence in GGV
Theorem 7.6 and Remarks 5.10–5.12. If \(A_0\ne (u,d)\),
the last regular corner is \((u,d)\) at \((1,0)\), with index
\(k>0\). By Theorem 7.6(7)–(8), its power exponent \(d_k\)
satisfies \(q_0\mid d_k\mid\gcd(u,d)\), so \(q_0\mid u\).
When \(i=u\), the negative second weight shows that \(A_0\)
is the bottom of the right face, namely \((u,\mu)\).
If also \(q_0=u\), GGV Corollary 7.10 forces that face
to have a single endpoint. Hence \(d=\mu=j\).

For \(u<4\) the first-coordinate bound is impossible.
For \(i=4\), the sum and upper bounds give \(j=12\).
Theorem 7.11 gives \((f_1,f_2)=(3,9)\), so \(q_0=4\).
This gives the asserted result for \(u=4\). For \(u=5,6\)
this primitive corner is impossible since \(4\nmid u\).
For \(i=5\) the same theorem forces \(i+j>20\).
Now \(15<j\le20\) and \(\gcd(i,j)>1\) give \(j=20\).
As \(2\le f_1<5\), we have \(q_0=5\), which proves the
claim for \(u=5\) and excludes \(i=5\) for \(u=6\).

For \(u=6\) it follows that \(i=6\). The inequalities now give
\(14<j\le30\) and \(3\mid j\). In Proposition 6.4, \(f_1=3\)
would give \(q_0=2\), leaving \(f_1=2,4,5\).
Keep only integral \(f_2=j f_1/6\) and directions
\[
 (\rho,\sigma)=
 \left(\frac{f_2-1}{e},\frac{1-f_1}{e}\right),
 \qquad e=\gcd(f_2-1,f_1-1),
 \qquad \rho\le6.
\]
The table lists those choices and the possible first coordinate \(i'\)
of the other endpoint \(A'_0=(i',j')\), which must satisfy
\(0\le j'<i'<6\) and have the same weight as \(A_0\).
All the candidate directions have \(\sigma=-1\), so this last
test is \(W/\rho\le i'< W/(\rho-1)\), with \(W=6\rho-j\).

| \(j\) | \((f_1,f_2)\) | \(\rho\) | \(W\) | \(i'\) |
|---:|---:|---:|---:|---:|
| 15 | \((2,5)\) | 4 | 9 | none |
| 15 | \((4,10)\) | 3 | 3 | 1 |
| 18 | \((2,6)\) | 5 | 12 | none |
| 21 | \((2,7)\) | 6 | 15 | none |
| 24 | \((4,16)\) | 5 | 6 | none |
| 30 | \((5,25)\) | 6 | 6 | 1 |

For \(j=30\) we have \(q_0=6\), so \(h\) has just the
translated root of multiplicity thirty. For \(j=15\) we have
\(q_0=3\). If the right face has one endpoint then \(h\) has only
a root of multiplicity fifteen. Otherwise the last power exponent
\(d_k\) is divisible by three, so the right initial form
\(x^{6a}h(y)^a\), up to a scalar, is an \(a d_k\)-th power in
the allowed Laurent polynomial ring. Unique factorization forces
each root multiplicity of \(h\) to be divisible by three.
The root lemma excludes multiplicity six. Standardizing again
at any root of multiplicity \(>6\) forces its multiplicity
to be fifteen or thirty, and the latter requires a single root.

For the remaining multiple-root lists let \(B\) be the polynomial
in the axial residue equation \(B'-(h'/6h)B=1\) from the root lemma.
Its simple zeros at the roots of \(h\) give \(\deg B\ge r+s\ge2\), so cancellation
of the leading term of the residue equation requires
\(\deg B=d/6=(5r+s)/2\). This gives \(6\mid d\) and \(s\le3r\). When \(d=30\), one root of multiplicity
fifteen and five roots of multiplicity three would violate that
inequality, leaving \((15,15)\) or the single root of multiplicity
thirty.

The [two-corner exclusion](PLANAR_PARTIAL_DEGREE_TWENTY.md)
settles the reduced profile in the next corollary and gives inversion
through common partial degree twenty.

**Corollary PlanarCommonPartialDegreeSeventeen.**
Every characteristic-zero planar Keller pair with both partial
degrees at most seventeen in one source variable is invertible.
Under a common partial-degree bound of twenty, every noninvertible
pair reduces by target power subtraction to partial degrees \((12,18)\).
In this reduced pair \(\deg h\ge30\), the leading coefficient
of the degree-eighteen output has degree at least ninety, and
if \(\deg h=30\) the normalized total degrees are \((72,108)\).

With positive partial degrees, use the leading-coefficient
power relation to subtract powers whenever one degree divides
the other, as in the [previous reduction](PARTIAL_DEGREE_FOURTEEN_AND_WEYL_TRANSFER.md#the-planar-frontier).
Each subtraction reduces one partial degree and preserves the
Keller property and invertibility. A zero partial degree gives
a triangular automorphism. A noninvertible pair thus ends at
nondividing degrees \(au,bu\), with \(2\le a<b\) and \(u\ge4\)
by the theorem. Under the common cap twenty we have \(u\le6\)
since \(bu\ge3u\). For \(u=4\), the normalization gives larger
total degree \(16b\le80\). For \(u=5\), it gives
\(25b\le100\). Moh's inclusive degree-100 theorem excludes
both cases. Its complex result applies over any
characteristic-zero field by embedding the finitely generated
coefficient field into \(\mathbb C\) and descending invertibility.
For \(u=6\) only \(a=2,b=3\) fits the cap. Its larger
total degree \(3(6+d)\) must exceed one hundred. Thus \(d\ge28\).
The root restrictions give \(3\mid d\) and hence \(d\ge30\).
The leading coefficient in question has degree \(3d\).
[Moh's statement of the bound](https://www.math.purdue.edu/~ttm/jacobian)

If the original pair has maximum partial degree eighteen, every
degree-eighteen member has leading-coefficient degree below ninety,
and the reduced maximum remains eighteen, then an original member
of degree eighteen survives unchanged up to scaling. Indeed a
dividing subtraction either lowers the original degree-eighteen
member or leaves one of the two original members of that degree
fixed. Once a member drops below eighteen it stays below.
The reduced pair cannot have the required leading coefficient.
This proves a coefficient-sensitive criterion at partial degree
eighteen.

## Transfer to polynomial Weyl pairs

**Corollary WeylCommonDifferentialOrderTwenty.**
Over characteristic zero, a scalar-commutator pair with nonzero
commutator in the first Weyl algebra generates if both differential
orders are at most twenty. More generally a nondividing order
pair with gcd four or five can occur only with common leading
factor degree twelve or twenty, respectively. In particular,
a proper pair containing a member of order fifteen reduces,
with that member fixed, to mate order twenty-five or eighty.

For a nondividing Weyl pair of positive orders \(au,bu\), the
[central reduction preserving the axial profile](PARTIAL_DEGREE_FOURTEEN_AND_WEYL_TRANSFER.md#the-weyl-to-plane-transfer-with-partial-degrees)
gives a complex planar Keller pair of the same partial orders
and the same common leading-factor degree. After the simultaneous
planar shear its total degrees are \(a(u+d),b(u+d)\), so by
degree divisibility for plane automorphisms the planar pair is
noninvertible. Theorem PlanarSmallPartialGcdCorners applies,
including its degree conclusions for gcd four and five.

Apply the same target power subtractions to the Weyl pair as
for a planar pair. They preserve the scalar commutator and
the generated subalgebra. If one order is zero, its member
is \(f(x)\) and the leading commutator with an order-\(n>0\)
mate has order \(n-1\) with leading coefficient, up to sign,
\(nf'(x)\) times the mate's leading coefficient. Hence \(n=1\)
and \(f'\) and the mate's leading coefficient are nonzero
constants. Then the members generate \(x,D\). Otherwise
the positive nondividing orders through twenty give just
\((12,18)\) by the planar corollary. The
[Weyl order-twelve reduction](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-remaining-differential-order-pair)
excludes that pair. Finally the
[fixed order-fifteen classification](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#the-order-fifteen-case)
gives mate orders \(10,20,25,80\). The first two pairs are
within the new cap.

## Applications to Hessian potentials and sextic maps

**Corollary HessianDirectionalEighteenthCoefficientBound.**
Assume Tang's retained complex structural alternative.
A four-variable complex polynomial \(H\) of nonzero constant
Hessian determinant has invertible gradient if there are
three independent constant directions \(w\) such that
\(D_w^{19}H\) is constant and \(\deg D_w^{18}H\le90\).
In particular it suffices that \(D_w^{18}H\) is constant
in each of the three directions.

In the cotangent case, as in the
[previous directional criterion](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#four-variable-hessian-potentials-with-small-directional-degrees),
one of the three directions projects to a nonzero direction
of the base plane. The coefficients of the cotangent variables
in \(D_w^{19}H\) and \(D_w^{18}H\) respectively bound both
planar partial degrees by eighteen and each possible
degree-eighteen leading coefficient by degree eighty-nine.
Use the planar corollary and invert the gradient in the same
way. If \(D_w^{18}H\) is constant then \(D_w^{19}H=0\).

**Corollary SexticCubicGenericFibreChart.**
Let \((P,Q,R)\) be a characteristic-zero sextic Keller triple.
Suppose there is a polynomial map \(Z(u,v)\) over \(k(t)\)
into \(R=t\) of geometric rank two everywhere. If
\(\deg_v Z_i\le2\) for every \(i\), the triple is invertible.
The same holds if \(\deg_v Z_i\le3\) and the coefficient
vector \(V_3(u)\) of \(v^3\) satisfies \(\deg_u V_3\le14\).
More generally in the cubic case it suffices that each nonzero
\(P_6(V_3),Q_6(V_3)\) has degree less than ninety in \(u\),
where \(P_6,Q_6\) are the ordinary degree-six homogeneous parts.

By the [rank-two polynomial generic-fibre argument](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#a-three-variable-generic-fibre-criterion),
\(P(Z),Q(Z)\) is a planar Keller pair over \(k(t)\), and its
inversion gives denominator-free inversion of the original triple.
Its common partial degree in \(v\) is at most twelve in the
quadratic case and eighteen in the cubic case. In the latter
case the only terms of \(v\)-degree eighteen come from
\(P_6(V_3),Q_6(V_3)\). The degree bound of fourteen gives
degree at most \(6\cdot14=84\) for both coefficients.
The coefficient-sensitive planar criterion applies. The
generic-fibre chart remains a hypothesis of this sextic result.
