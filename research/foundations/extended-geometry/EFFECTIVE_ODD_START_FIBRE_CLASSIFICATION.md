# An effective classification of every fixed algebraic odd-height fibre

Fix an odd positive integer \(a\) and \(u\in K(x)\) with
\(u(0)=0,\ u'(0)=1\), where \(K\) is a presented number field.
Write \(N=\deg u\). The fibre consists of the normalized polynomial pairs
\(C(0)=-D(0)=1\) for the two-height sphere phase
\[
 p=t^a(C+t^2D)/x,\qquad t^2+2xy=1,
\]
with every positive pure moment zero and selected endpoint coordinate
exactly \(u\). No coefficient-degree bound is imposed.

The retained
[fixed-coordinate theorem](ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md)
classifies this fibre using reciprocal units. Its finite representatives
were obtained from Dickson's lemma without an effective bound. This
packet closes that particular gap. It does not classify all endpoint
coordinates at once, arbitrary sphere phases, or unrestricted primitive
degrees five and seven.

## 1. A terminating classification theorem

**Theorem EffectiveAlgebraicOddFibre.** There is a terminating procedure
which, from \(K,a,u\), returns the complete fibre over \(K\) as a finite
list of explicit isolated pairs and explicit one-parameter families
\[
 \eta=\eta_iR^{j_i+d k},\quad k\ge0,\qquad
 X=\frac{u}{xu'}\frac{1+\eta}{2},\quad
 C=\frac{xX^a}{u},\quad D=-\frac{xX^{a+2}}u.                 \tag{E1}
\]
Here \(R\in K[x]\), \(R(0)=1\), is the unique minimal nonconstant
polynomial reciprocal unit, when one exists. All surviving families
have the same positive modulus \(d\). In particular emptiness,
finiteness and infinitude are decidable, and a finite fibre can be
listed in full.

The same procedure returns the geometric fibre over
\(\overline{\mathbb Q}\) after one effectively constructible finite
extension of the input number field. Its outputs contain every
geometric pair; extensions do not introduce continuous coefficient
parameters.

The procedure consists of bounded polynomial-factor enumeration, a
finite integer box when \(R\) does not exist, an explicit finite box of
pole twists when it does, and finitely many rational-matrix Orbit
Problem instances. Only the last decision subproblem is imported as a
polynomial-time result. No polynomial-time bound, practical small
cutoff, or general-purpose implementation of the whole classification
is asserted here.

## 2. The exact retained input

Let \(\mathcal C_u\) be the marked irreducible component of
\(u(x)+u(y)=0\) through \((0,0)\), on its smooth projective model.
The marked return need not be rational on the original \(x\)-line.
The retained theorem gives a bijection (E1) with rational functions
\(\eta\) satisfying all of the following:

* \(\eta(0)=1\) and \(\eta(x)\eta(y)=1\) on \(\mathcal C_u\).
* Finite poles occur only at poles of \(u\). Above a pole of order \(m\),
  their order is at most \(1+\lfloor m/(a+2)\rfloor\).
* \(\eta\equiv-1\pmod J\), where \(J=J_{a,u}\) is monic, has no root
  at zero, and has multiplicity \(\lceil s/a\rceil-1\) at a nonzero
  zero of \(u\) of order \(s\), multiplicity
  \(\operatorname{ord}_\alpha u'\) at a point with finite nonzero
  \(u(\alpha)\), and multiplicity zero otherwise.

Only positive multiplicities contribute to \(J\). Its roots are
disjoint from the poles of \(u\), so all admissible denominators are
invertible modulo \(J\).

Set \(c=u(\infty)\), with \(-\infty=\infty\), and
\[
 S=u^{-1}(\infty)\cup u^{-1}(-c)\cup\{\infty\}.             \tag{E2}
\]
If its finite nonzero geometric points are
\(\alpha_1,\ldots,\alpha_s\), every pre-congruence candidate is uniquely
\[
 \eta=\prod_{i=1}^s(1-x/\alpha_i)^{n_i},\quad
 Mn=0,\quad n_i\ge\ell_i,                                 \tag{E3}
\]
where \(M\) records the full pullback divisors of
\((1-x/\alpha_i)(1-y/\alpha_i)\) on the marked curve, and
\(\ell_i=-1-\lfloor m_i/(a+2)\rfloor\) at poles of \(u\),
zero otherwise. Normalization makes a zero divisor equivalent to
reciprocity. Ramification multiplicities must be retained.

The other retained structural input is that polynomial reciprocal
units are either only \(1\), or precisely \(R^k,\ k\ge0\), for one
unique normalized \(R\). In the latter case every rational reciprocal
unit and \(u\) belong to \(K(R)\), and
\[
 u=U(R),\quad U(1/z)=-U(z),\quad \deg R\mid N.              \tag{E4}
\]
The geometric canonical \(R\) descends to \(K\) by uniqueness.

The component and all identity tests are effective: factor the
numerator of \(u(x)+u(y)\) and select its unique factor through
\((0,0)\). That point is smooth since both partial derivatives are
one. Thus its geometric component is unique and defined over \(K\).
After clearing denominators, polynomial division by this factor
tests reciprocity exactly. The full reducible correspondence would
give an incorrect test.

## 3. A bounded search for the canonical polynomial

**Lemma BoundedCanonicalPolynomial.** If \(R\) exists, every finite
zero of \(R\) lies in \(u^{-1}(-c)\setminus\{0\}\), and
\[
                  1\le\deg R\le N,\qquad \deg R\mid N.    \tag{E5}
\]
Consequently \(R\), or its nonexistence, can be determined by a finite
factor search of degree at most \(N\).

Indeed, above a zero \(x=\alpha\) of a polynomial reciprocal unit,
reciprocity forces the other coordinate \(y\) to be its pole
\(\infty\). The identity \(u(x)=-u(y)\) therefore gives
\(u(\alpha)=-c\). Both projections are finite nonconstant maps on the
smooth curve, so a point above \(\alpha\) exists. The degree assertion
is (E4).

Factor the polynomial defining the finite fibre \(u=-c\), discard a
factor supported at zero, and normalize every remaining irreducible
factor \(f\) to \(f/f(0)\). Enumerate products of their nonnegative
integer powers having positive total degree dividing \(N\).
Their degree is at most \(N\), so this list is finite. Test reciprocity
on the marked component. A successful product of least degree is
the unique \(R\); if none succeeds, no polynomial reciprocal unit
exists even geometrically. The last assertion uses the retained
Galois descent of the unique geometric generator, so factoring over
the original \(K\) suffices for this step.

This does not search merely for divisors with multiplicities bounded
by the original fibre polynomial: the bound is total degree \(N\).

## 4. A numerical finite box when no polynomial unit exists

**Lemma EffectiveFiniteReciprocalBox.** Suppose that the preceding
search finds no \(R\). Put
\[
 B=1+\lfloor N/(a+2)\rfloor,\qquad
 H=s!\,B(2N)^{s-1}\quad(s\ge1).                           \tag{E6}
\]
Every candidate (E3) satisfies
\[
                       \ell_i\le n_i\le H.               \tag{E7}
\]
For \(s=0\), the only pre-congruence candidate is \(\eta=1\).
Thus at most \((H+B+1)^s\) vectors need be tested when \(s\ge1\).
Also \(s\le2N\), so this is a bound determined by the input degree
and height.

Each projection of \(\mathcal C_u\) has degree at most \(N\).
Every matrix entry in (E3), being the sum of two pullback orders of
a linear factor, has absolute value at most \(2N\).
Also \(-B\le\ell_i\le0\).

Consider the nonempty rational polyhedron
\[
 \mathcal P=\{v\in\mathbb R^s:Mv=0,\ v_i\ge\ell_i\}.
\]
Its recession cone is \(\{v:Mv=0,\ v_i\ge0\}\).
A nonzero vector in that rational cone yields a nonzero rational
one: keep its positive support and approximate in the rational
kernel on that support. Clearing denominators gives an integer
nonnegative kernel vector, hence a nonconstant polynomial reciprocal
unit. This is excluded. Thus the recession cone is zero and
\(\mathcal P\) is bounded. Explicitly, an unbounded sequence,
normalized by its norm and restricted to a convergent subsequence,
would produce a nonzero vector in that cone.

Every coordinate extremum is attained at a vertex. One can see this
without an optimization oracle: at an extremum, if the active rows
of \(M\) and the equalities \(v_i=\ell_i\) do not span the full dual
space, move in their common nullspace until another bound becomes
active. Both local directions are feasible, so the extremal
coordinate is unchanged. Boundedness forces a new independent active
bound. Repeating gives a vertex.

A vertex is determined by \(s\) independent integer rows taken from
\(M\) and coordinate rows. Their entries have magnitude at most
\(2N\), their right sides magnitude at most \(B\), and their
nonzero integer determinant has magnitude at least one.
Cramer's rule bounds each coordinate by
\(s!B(2N)^{s-1}\). This proves (E7).

For the algorithm, \(M\) need not actually be computed. Enumerate
the box (E7), use the marked-component identity test for reciprocity,
and then test the exact congruence modulo \(J\). Retain those
functions whose coefficients lie in \(K\). Formula (E1) returns the
complete finite fibre.

## 5. Explicit pole twists when the canonical polynomial exists

Write \(z=R(x)\). Compute \(U\) in (E4) by solving
\(A(R)=uB(R)\) for numerator and denominator degrees at most
\(N/\deg R\), then reduce \(A/B\). This is finite linear algebra
over \(K\), with a nonzero denominator.

Take one representative \(\alpha\) from each pair
\(\{\alpha,\alpha^{-1}\}\) of nonzero finite poles of \(U\), excluding
the fixed values \(1,-1\). The value \(1\) is not a pole since
\(U(1)=0\). The pole set is invariant under inversion. Set
\[
 t_\alpha(z)=\frac{z-\alpha}{1-\alpha z};
 \qquad t_\alpha(1)=1,\quad t_\alpha(1/z)=t_\alpha(z)^{-1}.
                                                               \tag{E8}
\]
For a finite pole value \(\alpha\) of \(U\), let \(m_\alpha\) be its
pole order and let
\[
 e_\alpha=\max_{R(\beta)=\alpha}\operatorname{ord}_\beta(R-\alpha),
 \qquad
 q_\alpha=\left\lfloor\frac{m_\alpha}{a+2}
                                  +\frac1{e_\alpha}\right\rfloor.
                                                               \tag{E9}
\]
All these preimages are finite, since \(R\) is polynomial.
Define
\[
 k_{\min}=
 \begin{cases}
 -q_0,&U\text{ has a pole at }0,\\
 0,&U\text{ is regular at }0.
 \end{cases}                                                   \tag{E10}
\]

**Lemma ExplicitReciprocalTwistBox.** The complete set of normalized
reciprocal units satisfying the pole bounds, before the \(J\) test,
is the disjoint union
\[
 \eta=\eta_{\mathbf n}R^j,\quad j\ge0,\qquad
 \eta_{\mathbf n}
 =R^{k_{\min}}\prod_\alpha t_\alpha(R)^{n_\alpha},
 \qquad -q_\alpha\le n_\alpha\le q_{\alpha^{-1}}.           \tag{E11}
\]
This is an explicitly bounded finite set of forward orbits.

To prove completeness, all reciprocal units lie in \(K(R)\), or in
the corresponding finite splitting-field extension. For the
involution \(z\mapsto1/z\), their divisors are anti-invariant.
Hence orders at \(1,-1\) are zero, the orders at \(0,\infty\) are
opposites, and the orders at \(\alpha,\alpha^{-1}\) are opposites.
A nonzero finite pole of the descended unit must be a pole of \(U\),
because every finite \(z\)-value has a finite preimage under \(R\).
The same holds for its paired zero. This yields the factorization
in (E11), with integer \(k=k_{\min}+j\); normalization fixes the
remaining scalar to one.

A pole of order \(h\) in the \(z\)-line pulls back at a point of
ramification index \(e\) to order \(eh\), while \(u\) has pole
order \(em_\alpha\). The exact retained bound is
\[
 eh\le1+\left\lfloor\frac{em_\alpha}{a+2}\right\rfloor
 \quad\Longleftrightarrow\quad
 h\le\frac{m_\alpha}{a+2}+\frac1e .
\]
The strongest condition is the largest \(e\), giving (E9).
The sign of \(n_\alpha\) specifies which member of the pair is its
pole and gives precisely the interval in (E11). At \(z=0\) only
negative \(k\) give finite poles; this proves (E10).
No upper bound arises from infinity.

Distinct twist vectors give distinct divisors away from \(0,\infty\),
so the forward orbits are disjoint. Compute them in a splitting
field and retain exactly the bases \(\eta_{\mathbf n}\) defined over
\(K\). This discards an entire orbit or none of it: multiplication
by \(R^j\in K(x)^\times\) cannot change whether the coefficients
descend to \(K\).

## 6. The exact Orbit Problem bridge

For each base \(\eta_0\), \(J=1\) accepts all \(j\ge0\).
Otherwise work in the finite \(K\)-algebra
\(\mathcal A=K[x]/J\). Inverses of admissible denominators are
computed by the extended Euclidean algorithm.
If \([\eta_0]\) is not a unit, no index can work. If \([R]\) is not
a unit, no positive index can work; test \(j=0\) separately.
In the remaining case the exact condition is
\[
                  [R]^j=-[\eta_0]^{-1}.                 \tag{E12}
\]
This algebra can be nonreduced: the full powers in \(J\) are retained.

**Lemma RationalOrbitReduction.** Equation (E12) is an instance of
the rational point-to-point Orbit Problem, including all higher jets.

Choose a rational basis of \(K\), then the induced rational basis of
\(\mathcal A\). Its dimension is \([K:\mathbb Q]\deg J\).
Multiplication by \([R]\) is a rational matrix \(M_R\). If \(v\)
denotes the coordinate vector, (E12) is exactly
\[
            M_R^jv(1)=v(-[\eta_0]^{-1}),\qquad j\ge0.     \tag{E13}
\]
There is equality of complete vectors, not the vanishing of a single
linear observation of an orbit.

For a directly inspected primary decision theorem, use Chonev,
Ouaknine and Worrell, Theorem 2.2, its one-dimensional-target clause.
It applies to rational matrices of arbitrary dimension and a rational
target line. Reduce (E13) to that clause by
\[
 \widetilde M=\operatorname{diag}(M_R,1),\quad
 \widetilde v=(v(1),1),\quad
 V=\mathbb Q\,(v(-[\eta_0]^{-1}),1).
                                                               \tag{E14}
\]
The last coordinate fixes the line multiplier to one. Thus a
nonnegative orbit hits \(V\) exactly when (E13) holds.
The primary proof in Sections 3--4 includes repeated eigenvalues;
no diagonalizability or reduced-algebra hypothesis is inserted.

After the decision returns yes, exact sequential evaluation finds a
least witness \(j_0\) and terminates. To describe all witnesses,
compute the minimal polynomial of \(M_R\). It has finite
multiplicative order exactly when that polynomial is squarefree
and is a product of cyclotomic polynomials. In that case the order
\(d\) is the least common multiple of their indices.
This test is finite: an index \(h\) with \(\varphi(h)\le d_0\),
where \(d_0\) is the matrix dimension, satisfies \(h\le2d_0^2\).
Indeed the product formula gives \(\varphi(h)^2/h\ge1/2\):
every odd-prime factor contributes at least one, and the factor
from two contributes at least \(1/2\).

The map \(j\mapsto[R]^j\) is a homomorphism from \(\mathbb Z\).
If \([R]\) has infinite order, a nonempty solution set is the single
integer \(j_0\). If it has finite order \(d\), the nonnegative
solutions are exactly \(j_0+d\mathbb Z_{\ge0}\), where \(j_0\) is
the least nonnegative solution. This proves all output and
termination assertions in Section 1.

For the geometric classification, the finitely many roots in (E2)
and in the pole factorization of \(U\) belong to an effectively
constructed finite number field. Normalization in (E3) and (E11)
puts every candidate in that one field, so running the same
procedure there returns the entire geometric fibre.

## 7. Controls, scope and what becomes possible

The accompanying exact checker exercises bounded generator extraction
on both rational and nonrational original returns, the finite
no-generator box on \(u=x/(1-x^2)\), ramification-sensitive twist
bounds, nontrivial inverse-paired pole twists, and rational orbit
matrices with empty, singleton and periodic solutions. Repeated-root
controls reject the incorrect shortcut of checking only values.
It is finite support for the written proof, not an implementation
or independent proof of the imported general Orbit algorithm.

The practical improvement is a complete finite prescription for a
specified algebraic endpoint coordinate. A failure to find a generator
now has a proved stopping bound; finite residual fibres have a
numerical box; and the exponent condition is decidable. The boxes can
be large. This theorem does not claim a global bound on endpoint
degree or solve the remaining global small-degree moment problem.

## 8. Primary source and evidence

V. Chonev, J. Ouaknine and J. Worrell, *On the Complexity of the Orbit
Problem*, Journal of the ACM 63(3), Article 23 (2016).
The inspected 40-page author manuscript is retained at
[sources/chonev-ouaknine-worrell-orbit-2016.pdf](sources/chonev-ouaknine-worrell-orbit-2016.pdf);
the exact import is Theorem 2.2, one-dimensional clause, page 3;
Sections 3--4, pages 5--9 supply the reduction and decision proof.
Primary author-hosted source:
https://people.mpi-sws.org/~joel/publications/orbit_journal_14.pdf .

The source ledger records exact source bytes, theorem scope, internal
proof review and finite-check limitations. No external human review,
formal verification or claim of historical priority is inferred.
