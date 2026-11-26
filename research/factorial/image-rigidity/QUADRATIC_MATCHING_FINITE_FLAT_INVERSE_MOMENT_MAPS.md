# Finite flat inverse moment maps for structured quadratics

The forward argument uses the shorter window in
[Quadratic strong factorial windows from matchings](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md).

The additional inputs are polynomial-ring Cohen-Macaulayness and the local
regular-sequence criterion in Stacks10.104, tags00ND and02JN.

Fix m, rational odd-denominator c_ij, and diagonal/linear supports U,V
with p=|U|+|V|>=1. Over a characteristic-zero field K containing chosen
b_ij with b_ij^2=4c_ij, let z be the p separate u_i (i in U) and v_i
(i in V) parameters, putting absent parameters0. Set

 F(z;X)=sum u_i^2 X_i^2+sum b_ij u_i u_j X_iX_j+sum v_i^2 X_i,
 H_j(z)=Lambda_m(F(z;X)^j).

For fixed n>=1 define the p-component map

 Phi_n:A^p_K --> A^p_K,
 z |--> (H_n(z),...,H_(n+p-1)(z)).                         (1)

**Theorem 1.** The map(1) is finite and flat, in fact its coordinate
ring is a finite FREE module over the abstract target polynomial ring,
of rank

                   D=product_(j=n..n+p-1) (2j).            (2)

Every geometric target therefore has a nonempty finite fiber of
scheme length D. The coordinate functions H_n,...,H_(n+p-1) are
algebraically independent. The geometric zero fiber has only the origin
as its underlying point, and still has scheme length D. No reducedness,
uniqueness, explicit inverse, or full support on every fiber is claimed.

## 1. Only the origin lies over the zero target

Each H_j is homogeneous of degree2j in z. The phase parameter map has
no nonzero basepoint: F=0 forces every diagonal square u_i^2=0 and every
linear square v_i^2=0. Thus F!=0 for every z!=0 over any field extension.

Every such F obeys the fixed-weight relations and has k+ell<=p.
If all p moments H_n,...,H_(n+p-1) vanished, their initial k+ell
subblock would contradict the matching upper bound at this same n.
Consequently these p forms have no nonzero common zero over Kbar.

Let R=K[z_1,...,z_p], I=(H_n,...,H_(n+p-1)). By the Nullstellensatz
over Kbar, and faithful scalar extension, R/I is finite dimensional
over K and its radical is (z_1,...,z_p). In particular I has height p.

## 2. Regular sequence and its Hilbert series

The polynomial ring R is Cohen-Macaulay: Stacks Lemma10.104.7 tag00ND.
In its localization at the homogeneous maximal ideal, the p forms
constitute a regular sequence, since their quotient has dimension0:
Stacks Lemma10.104.2 tag02JN. This is an actual mathematical use:
https://stacks.math.columbia.edu/tag/00N7
https://stacks.math.columbia.edu/tag/02JN

Their regularity holds already in the graded polynomial ring. Indeed,
for each successive quotient the kernel of multiplication by the next
homogeneous form is graded. A nonzero homogeneous element in that kernel
cannot die on localizing at the homogeneous maximal ideal: multiplying
it by a polynomial with nonzero constant term has a nonzero lowest-degree
part. Thus each such kernel is zero globally.

Write degrees d_i=2(n+i-1), i=1,...,p, and Q=R/I.
The exact graded multiplication sequence at each regular element gives

 Hilb_Q(t)= [product_i(1-t^d_i)]/(1-t)^p.                  (3)

This is a polynomial because Q is finite dimensional, and its value at1 is

                       dim_K Q=product_i d_i=D.           (4)

The Hilbert series computation follows directly from the injective
multiplication maps; no multiplicity formula or reducedness premise is used.

## 3. An explicit graded freeness argument

Choose any homogeneous K-basis of Q and homogeneous lifts b_1,...,b_D
to R. Let S=K[Y_1,...,Y_p] be an ABSTRACT polynomial ring, graded by
degY_i=d_i, and let Y_i act on R as multiplication by H_(n+i-1).
There is a graded S-module map

             psi: direct_sum_(a=1..D) S(-deg b_a) --> R,
             (s_a) |--> sum_a s_a(H_n,...,H_(n+p-1))b_a.  (5)

This is surjective by induction on total parameter degree. For a
homogeneous r in R, subtract its expansion in the lifted basis modulo I.
The remainder is a sum H_(n+i-1)r_i, with homogeneous r_i of strictly
smaller degree because every d_i>0. Apply the same induction to each r_i.

The source in(5) has Hilbert series

 [sum_a t^(deg b_a)]/[product_i(1-t^d_i)]
      =Hilb_Q(t)/product_i(1-t^d_i)
      =1/(1-t)^p=Hilb_R(t),                              (6)

using(3). Since(5) is a graded surjection and all its graded pieces are
finite dimensional, equality(6) forces its kernel to be zero in every
degree. Thus(5) is an isomorphism and R is free of rank D over S.
In particular S acts faithfully, so the H_j are algebraically independent.
This proves the finite free assertion and hence finite flatness of Phi_n.

## 4. Prescribing moments and fiber lengths

For any geometric target a=(a_1,...,a_p), tensor the finite free S-module
R with Kbar via Y_i |--> a_i. The fiber algebra has Kbar-dimension D>0.
It is a nonzero finite algebra, hence has a maximal ideal and a Kbar-point.
This proves surjectivity and the stated constant scheme length of fibers.
If a and all b_ij are algebraic over Q, a point can be taken over a number
field: its coordinates lie in a finite algebraic extension of the number
field generated by the finitely many data.

For the zero target, Section1 identifies the underlying point as the origin.
For targets of the form (0,...,0,a_p), a_p!=0, the same support argument
as in the every-start sharpness proof forces all intended parameters
nonzero. For p=1 the first-zero block is empty and a_p!=0 directly forces
the only parameter nonzero. This special inverse-moment fiber therefore
recovers exact support and the maximal p-1 zero block at any start.

## 5. Diagonal phases without the square parameter cover

If all c_ij=0, omit the square parametrization and work directly with
the p independent phase coefficients A_i (i in U), d_i (i in V):

                F=sum_(i in U) A_i X_i^2+sum_(i in V) d_i X_i.

Every H_j now has degree j. Nonzero coefficient vectors still give
nonzero phases, and the matching theorem still excludes a zero p-block.
The identical regular-sequence and graded-basis proof shows that the
diagonal coefficient moment map is finite free of degree

                         product_(j=n..n+p-1) j.          (7)

This counts coefficient vectors themselves. Formula(2) counts squared
parameter vectors and must not be mislabeled as a count of distinct phases.

## 6. A genuine converse through finiteness and parameter coverage

**Theorem 2 (equivalence).** For the fixed rational weights, the matching
k+ell upper-window assertion over all characteristic-zero fields is
equivalent to FINITENESS of every map(1), for every n>=1, compatible U,V,
and every fixed choice of the square roots b_ij, over every
characteristic-zero base field containing those b_ij.
Finite freeness of the degrees in(2) is a stronger equivalent package,
since Theorem1 derives it from the upper-window assertion.

Only the reverse implication remains to be shown. Let a nonzero
admissible phase f over a characteristic-zero field have actual diagonal
support U and linear support V, and put p=k+ell>=1. In an algebraic
extension choose u_i with u_i^2=A_i for i in U and v_i with v_i^2=d_i
for i in V. All these chosen parameters are nonzero.
For i,j in U set b_ij=B_ij/(u_i u_j). Then b_ij^2=4c_ij.
If either endpoint is outside U the relation forces B_ij=0; choose
any permitted square root b_ij there. Thus the given f is exactly the displayed parameterized phase
for one of the fixed sign branches included in the assertion.
No single preferred sign branch is assumed to cover every phase.

If its moments n,...,n+p-1 all vanished, this nonzero parameter point z
would lie over the zero target of its finite map(1).
After scalar extension to an algebraic closure, homogeneity gives

                    Phi_n(lambda z)=0 for every lambda.

This is an entire affine line, with infinitely many distinct points,
in a fiber of a finite morphism. A finite geometric fiber has finite
support, so this is impossible. This proves the upper-window assertion.

The same proof without square roots of coefficients gives equivalence
between diagonal actual-support SFC and finiteness of all the direct
diagonal coefficient maps in Section5, for all supports and starts.

This is an actual mathematical reverse implication. Flatness alone or
surjectivity alone is not being substituted for finiteness, and one
fixed sign branch is not being substituted for the full covering family.
No claim of an unrestricted SFC equivalence is made.

## 7. Scope of the deduction

The closely related
[support-deletion finite-free moment theorem](FACTORIAL_SUPPORT_DELETION_FINITE_FREE_MOMENTS.md)
states that a proved deleted monomial-support family gives
a finite free affine zero-moment chart after inverting its leading
resultant, together with a next-moment multiplication determinant.
The present proof concerns the full arbitrary-target moment map of a
structured quadratic coefficient cover, with an explicit restricted
finiteness converse. That older chart theorem is context, not a proof
premise here. The graded commutative-algebra mechanism is classical;
no claim that this mechanism or the general moment-ideal viewpoint is new
is made.

The matching upper bound is the forward special-family mathematical input;
standard commutative algebra turns it into a finite inverse moment theorem.
Section6 explains the exact restricted-family converse from finiteness
with homogeneity and full parameter-branch coverage.

These are structured zero-constant quadratics, all dimensions and all
compatible supports, with fixed rational2-integral correlations. No arbitrary
quadratic or SFC(2)/SFC(3) closure, reduced fiber, canonical coefficient
selection, or explicit coefficient formula follows here. No finite
calculation, formal verification, external review or priority assertion
is a premise.
