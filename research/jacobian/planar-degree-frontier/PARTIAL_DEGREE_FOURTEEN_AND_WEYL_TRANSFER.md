# Partial degree fourteen and preservation of axial data in Weyl reduction

## Simultaneous triangular normalization

**Theorem AxisDegreeToTotalDegreeNormalization.**
Let \(k\) have characteristic zero. Consider a planar Keller pair or
a scalar-commutator pair in the first Weyl algebra, with positive
partial degrees or differential orders \(m=au,n=bu\) in the
second variable, \(\gcd(a,b)=1\). Write the top coefficients as
\(\alpha h(x)^a,\beta h(x)^b\), with \(h\) monic and \(v=\deg h\).
A single polynomial shear in the second variable, fixing \(x\),
preserves the top rows and makes the total or Bernstein degrees
exactly \(a(u+v),b(u+v)\).

The common-factor identity is the
[axial leading-factor calculation](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#one-translation-for-the-two-leading-coefficients).
In the Weyl setting use PBW coordinates \(x^iD^j\), where
\([D,x]=1\). We use the
[positive-weight classification of rational slices](../../dixmier/dixmier-formal-lifting/WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#positive-weights-and-compatibility-of-exposed-edges)
for the leading symbols in both settings. The shift in the Weyl
case provides a short alternative to the existing
[weighted normalization](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree)
when the pair has positive orders.

The two top corners are \((av,m)\) and \((bv,n)\). Let \(L\) be
the maximum of \(1\) and of all ratios \((i-av)/(m-j)\) from
supported monomials of the first member with \(j<m\), and all
ratios \((i-bv)/(n-j)\) from the second with \(j<n\). At \(L=1\)
the desired total degrees follow from the defining inequalities
and the actual corner coefficients. Otherwise write \(L=s/r>1\)
in lowest terms. Both initial forms at \((r,s)\) contain their
top corner, and at least one has another term. Their summed
weights are strictly greater than \(r+s\), since \(m+n\ge2\)
and \(s>r\). Therefore the leading Poisson bracket vanishes.
For weighted forms of positive weights \(a(r v+s u)\) and
\(b(r v+s u)\), contraction of \(dp\wedge dq=0\) with the
weighted Euler field gives \(a p\,dq=b q\,dp\). Their powers
\(p^b,q^a\) are proportional.

The positive-weight classification for \(s>r\) forces
\(r=1\) and the nonmonomial form to be a scalar times
\(x^i(y-\lambda x^s)^j\). Its top \(y\)-corner determines
\(i=av,j=m\) or \(i=bv,j=n\), respectively. By the proportionality
of powers the other form uses the same \(\lambda\). The coefficient
of the term one step below the top \(y\)-degree determines
\(\lambda\) in \(k\), by characteristic zero. The simultaneous
shear \(y\mapsto y+\lambda x^s\) removes that edge from both
forms while preserving both \(y\)-degrees and both top rows.
In the Weyl algebra \(D\mapsto D+\lambda x^s\) has the same
associated graded substitution, since every contraction
lowers weight by \(1+s\). In either case no term formerly of
strictly smaller weight can reach the top weight. Thus the
new \(L\) is strictly smaller. Each value of \(L>1\) is an
integer by the classification, so the procedure terminates.
The composition still fixes \(x\) and adds a polynomial in
\(x\) to the second variable.

## The planar frontier

**Theorem PlanarPartialDegreeFourteenResidual.**
Over a characteristic-zero field, planar Keller pairs of
positive partial-degree gcd at most two are invertible.
For a noninvertible planar Keller pair with both partial
degrees at most fourteen in one source variable, reduce
dividing partial degrees by subtracting target powers. The
result has one of the following two unordered profiles. In
the table \(h\) is the common monic leading factor in that
reduced pair.

| Reduced partial degrees | Gcd \(u\) | Possible \(\deg h=v\) | Degree of leading coefficient of degree-twelve output |
|---|---:|---|---:|
| \((8,12)\) | 4 | \(4\mid v,\quad v\ge32\) | \(3v\ge96\) |
| \((9,12)\) | 3 | \(3\mid v,\quad v\ge24\) | \(4v\ge96\) |

It suffices to work over an algebraic closure, with inverse
descent by faithful flatness. When \(m,n>0\), put \(u=\gcd(m,n)\)
and \(v=\deg h\). The normalization gives total-degree gcd
\(u+v\). By Guccione–Guccione–Valqui, Corollary 6.6, this sum
is at least sixteen in a counterexample. Their \(B\) is the
minimum over all counterexamples, so a lower bound on \(B\)
is a lower bound for each pair. [Primary text](https://arxiv.org/pdf/1401.1784v3)

There are two further tests for the axial data. First, if
\(\gcd(u,v)=1\), a sufficiently large shear \(y\mapsto y+x^N\)
gives total degrees \(a(v+Nu),b(v+Nu)\), since the joint
top row dominates all lower rows. By Dirichlet we can choose
\(N\) with \(v+Nu\) prime. Apply the
[Nagata prime-gcd theorem](EXTERNAL_JC2_AND_HC4_INPUTS.md#nagata-prime-total-degree-gcd-interface).
Second, \(v=u\) is impossible because a sufficiently large
positive weight on \(y\), relative to the \(x\)-weight, exposes
\(x^{au}y^{au}\), which has no rational slice. If \(v>u\)
and \(a,b>1\), the
[standard endpoint obstruction](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#the-endpoint-interface-and-the-planar-theorem)
requires \(\gcd(u,v)>2\). The total degrees required by
that interface can now be taken to be \(a(u+v),b(u+v)\).

When one positive partial degree divides the other, keeping
the smaller component fixed and subtracting powers of it
lowers the larger partial degree until it is zero or no
longer a multiple. For an original pair of gcd at most two,
if the smaller degree divides the larger then the smaller
is itself at most two. The final gcd is still at most two
or the zero case has been reached. A zero partial degree
component gives a triangular automorphism directly from the
Jacobian condition. For nondividing degrees with \(u=1\)
the prime-gcd test applies. For \(u=2\), \(v=0,1,2\) is
excluded for a counterexample by the three tests above,
and \(v>2\) contradicts \(\gcd(u,v)>2\).
In particular, the degree pair \((10,12)\) is settled
without a separate low-scale calculation.

Under the cap fourteen, after elimination of divisibility
and gcd at most two, the only pairs are \((6,9),(8,12),(9,12)\).
Indeed they are the only products \(u(a,b)\) in the range
with \(u\ge3\) and coprime \(2\le a<b\). The
[previous max-11 theorem](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#the-endpoint-interface-and-the-planar-theorem)
supplies \((6,9)\). In either degree-twelve case,
\(u+v\ge16\) implies \(v>u\). Thus the endpoint obstruction
forces \(u\mid v\). Under the simultaneous normalization
the larger total degree is \(3(4+v)\) or \(4(3+v)\),
respectively. Moh's inclusive degree-100 bound gives
\(3(4+v)>100\) or \(4(3+v)>100\), yielding the table.
Moh's complex theorem applies over any characteristic-zero
field by embedding the finitely generated coefficient field
into \(\mathbb C\), using faithful flatness there, and
extending scalars back to the original field.
[Author's statement of the degree bound](https://www.math.purdue.edu/~ttm/jacobian)

In particular, for a Keller pair of maximum partial degree
at most twelve, if every degree-twelve member has leading
coefficient of degree less than \(96\) in the other variable,
the pair is invertible. A dividing reduction either drops
the maximum below twelve, or retains one of the original
degree-twelve members unchanged up to nonzero scaling.
In the latter case the unchanged coefficient violates the
necessary lower bound in the table. This suffices to improve
the fibre and Hessian tests below without making any
exhaustiveness claim about the possible ordinary total
degrees above Moh's bound.

## The Weyl-to-plane transfer with partial degrees

**Proposition WeylAxisDegreeProfileTransfersToPlane.**
A characteristic-zero scalar-commutator Weyl pair of positive
differential orders \(m,n\) and leading-coefficient degrees
\(i,j\) yields a complex planar Keller pair of partial degrees
\(m,n\) in one variable and leading-coefficient degrees \(i,j\)
in the other. In particular if the orders are nondividing,
with common factor degree \(v\) and order gcd \(u\), then
\[
       \gcd(u,v)>1,\qquad
       v>u\ \Longrightarrow\ \gcd(u,v)>2.
\]
These constraints use the axis data of the Weyl pair itself.

Scale to commutator one, ordering \(P,Q\) as the images of
\(x,D\). First specialize the finitely many coefficients
to a number field, preserving the commutator, the support
bounds and the two top corner coefficients. This uses
the finite polynomial equations for the PBW coefficients
and inverse variables for the nonvanishing conditions.
Spread out the specialization over a localization of the
ring of integers, including inverses of the two corner
coefficients. Reduce at maximal ideals of arbitrarily large
prime residue characteristic \(p\). The images of the
central generators \(x^p,D^p\) are \(P^p,Q^p\). By
Lauritzen–Thomsen, Theorem 1.7 and Propositions 1.10–1.12,
they belong to the center and the induced polynomial pair
has nonzero constant Jacobian. [Primary text](https://arxiv.org/pdf/1508.02499v1)

In characteristic \(p\) the differential-order filtration
has associated graded ring \(k[x,\xi]\). Thus the principal
differential term of \(P^p\) is \(\operatorname{lc}(P)^p
\xi^{pm}\), with degree \(pi\) in \(x\) of its coefficient.
The analogous term for \(Q^p\) has order \(pn\) and coefficient
degree \(pj\). Dividing the monomial exponents by \(p\) when
writing each central polynomial in \(x^p,D^p\) gives the
claimed four degrees. Independently the Bernstein degree
of each central polynomial in those variables is at most
the original Bernstein degree of the corresponding operator.
These uniform finite support bounds, the exact axial degrees,
and the unit Jacobian are polynomial equations with
nonvanishing conditions. A complex solution exists: otherwise
the Nullstellensatz, with inverse variables for the
nonvanishing conditions, gives an identity over a number
field certifying inconsistency, and after clearing denominators
the finite integer system has no solution at almost every
prime. Equivalently, use faithful flatness to obtain the
inconsistency identity over \(\mathbb Q\) for the integer
equations. This is the same arithmetic mechanism as in the
[total-degree transfer](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FORTY_NINE_FROM_F7.md#two-inputs-from-polynomial-degree-pairs)
with the differential-order filtration retaining the extra data.

For nondividing \(m,n\), the complex pair cannot be an
automorphism: after the simultaneous shear its total degrees
are \(a(u+v),b(u+v)\) with coprime \(a,b>1\), which fail degree
divisibility for plane polynomial automorphisms. For
\(\gcd(u,v)=1\) this contradicts the prime-gcd test. For \(v>u\)
the planar endpoint theorem gives the second inequality.
In particular a single-root leading factor of multiplicity
\(v>u\) with \(\gcd(u,v)=2\) fails this test irrespective of
the quotient orders. It is the planar standard endpoint
obstruction that provides the inequality on these Weyl data.

## Hessian and three-variable consequences

**Corollary HessianDirectionalTwelfthCoefficientBound.**
Assume Tang's complex structural alternative. Let \(H\) have
four variables and nonzero constant Hessian determinant.
If there are three linearly independent constant directions
\(w\) with \(D_w^{13}H\) constant and \(\deg D_w^{12}H\le96\),
then \(\nabla H\) is a polynomial automorphism.

In the cotangent alternative of
[Tang's structural theorem](../../hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md),
one of the three directions projects nontrivially to the
source plane. The coefficients of the two fibre variables
in \(D_w^{13}H\) show that both plane partial degrees in this
projected direction are at most twelve. The same coefficients
in \(D_w^{12}H\) bound each leading coefficient of degree
twelve in that plane frame by total degree \(95\).
The planar criterion applies, and the gradient is inverted
as in the [previous directional consequence](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#four-variable-hessian-potentials-with-small-directional-degrees).

**Corollary SexticQuadraticGenericFibreChart.**
Let \((P,Q,R)\) be a characteristic-zero sextic Keller
triple. Suppose the generic fibre \(R=t\) has a polynomial
parametrization \(Z(u,v)\) over \(k(t)\) of geometric rank
two everywhere, of degree at most two in \(v\). If the
coefficient vector of \(v^2\) has degree at most fifteen in
\(u\), then the triple is invertible. The other coefficients
of the parametrization are unrestricted.

The proof of the
[generic-fibre chart criterion](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#a-three-variable-generic-fibre-criterion)
gives a planar Keller pair \(P(Z),Q(Z)\). Both have
\(v\)-degree at most twelve. A \(v^{12}\) coefficient
is obtained by evaluating the ordinary degree-six part
on the coefficient vector of \(v^2\). Its \(u\)-degree
is at most \(6\cdot15=90\), below the bound \(96\).
Apply the planar result and use the same denominator-free
descent from the generic fibre. The full-rank polynomial
hypothesis is retained.
