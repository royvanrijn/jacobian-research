# Root-polynomial feedback: exact acyclicity and polynomial factors

15 November 2026. Three written all-degree claims with selected exact controls. No independent whole-proof replay, formal verification, external review or novelty certificate. OP-KDESC remains open. This closes a class of non-descended feedback left outside the [higher-return theorem](../universal-keller-higher-return/PROOF.md), without classifying all higher-rank selectors.

Work over a characteristic-zero field k. Identity parameters u are allowed. For N>=1 independent source triples P_i=(x_i,y_i,z_i), put

    t_i=1+x_i*y_i, q_i=t_i^2*z_i+y_i^2*(1+3*t_i),
    sigma_i=t_i*q_i, X_i=x_i*q_i,
    B0_i=-y_i/2-3*x_i*q_i/2,
    C0_i=x_i*(5-3*t_i)-x_i^3*z_i.

Consider the ENTIRE polynomial map, keeping u unchanged,

    T_i=(sigma_i, B0_i+F_i(sigma,X,u), C0_i+H_i(sigma,X,u)),       (0.1)

where F_i,H_i are arbitrary polynomials in ALL signals sigma and ALL root expressions X. These X_i are polynomial source functions, not old output coordinates. No degree, number-of-monomials, graph or native-profile ansatz is imposed. The restrictions are that every first output sigma_i is retained and that the two corrections lie in the named polynomial subring k[u,sigma,X]. General source corrections, changed first outputs, nonreduced coefficient rings and positive characteristic are outside the theorem.

## 1. UAROOTDAG1: Keller is equivalent to a self equation and an acyclic graph

Define the effective dependency graph on the N blocks by i->j, i!=j, exactly when F_i or H_i depends on X_j as a polynomial. Dependence only on sigma_j does not give an edge. Collect actual polynomials first: apparent coefficient-list cycles can cancel.

Then the following are equivalent:

1. The whole map (0.1) has constant nonzero Jacobian.
2. Its Jacobian is one.
3. For every i,

       sigma_i^2 * partial_(X_i) H_i
          =2*X_i^2 * partial_(X_i) F_i,                         (1.1)

   and the effective dependency graph is acyclic.

The exact self-equation normal form is

    F_i=f_i0+sigma_i^2 * sum_(n>=1) b_in * X_i^n,
    H_i=h_i0+sum_(n>=1) (2*n/(n+2))*b_in * X_i^(n+2),          (1.2)

with finite sums, and f_i0,h_i0,b_in in k[u,sigma,X_others]. This characterizes ALL root-polynomial corrections in (0.1), not only the earlier native higher profiles.

### The independent diagonal variables

UAM1 gives birational coordinates sigma_i,S_i,r_i on each source field:

    S_i=x_i/t_i, r_i=1/t_i, X_i=sigma_i*S_i,
    x_i=S_i/r_i, y_i=(1-r_i)/S_i,
    B0_i=(r_i-1-3*sigma_i*S_i^2)/(2*S_i),
    C0_i=(r_i+1)*S_i-sigma_i*S_i^3.                           (1.3)

In particular all sigma_i,X_i,r_i,u are algebraically independent and generate the full source function field. This is a proof chart only; (0.1) remains defined on all affine space, including t_i=0 and sigma_i=0. The [fixed-signal leaf calculation](../universal-keller-higher-return/PROOF.md) explains why such charts are not polynomial affine-plane coordinates.

Work first at fixed sigma,u, in variables r,S. Both corrected outputs are affine in r. Eliminate the r columns using partial B_i/partial r_i=1/(2*S_i) and partial C_i/partial r_i=S_i. Comparing with the determinant-one uncorrected product gives the exact rational identity

    det DT = det(diag(r_1,...,r_N)+A)/product_i r_i,
    A_ij=(sigma_j/2)*(partial_(X_j) H_i
                       -2*X_i^2/sigma_i^2 * partial_(X_j) F_i).       (1.4)

The matrix A is independent of ALL r_i. If det DT=c is constant, comparison of the coefficient of product r_i gives c=1. Expanding the determinant in the independent r_i gives

    det A_J=0 for EVERY nonempty principal index set J.         (1.5)

This is stronger than ordinary nilpotence, which only sets sums of principal minors to zero. Conversely (1.5) makes (1.4) equal one.

### The classical matrix lemma and the actual dependency graph

A matrix over a field, or domain, has all nonempty principal minors zero exactly when its directed nonzero pattern is acyclic. This is classical: see Berson, Lemma5.7 and Remark5.8, explicitly referring to van den Essen's Corollary6.3.9. We include the short proof. If the pattern had a shortest directed cycle, its vertices could have no chord: a chord gives a shorter cycle. On those vertices the determinant therefore has exactly one nonzero cycle term, the product of the cycle entries up to sign, contrary to (1.5). An acyclic pattern can be ordered strictly triangular, which gives the converse. The domain assumption prevents a product of nonzero entries from vanishing.

The diagonal part A_ii=0 is (1.1). For i!=j put

    K_ij=partial_(X_j) H_i-2*X_i^2/sigma_i^2*partial_(X_j) F_i.

Differentiating (1.1) gives

    partial_(X_i) K_ij=-4*X_i/sigma_i^2 * partial_(X_j) F_i.    (1.6)

Thus K_ij=0 implies that both partial_(X_j) F_i and partial_(X_j) H_i vanish. Characteristic zero makes this equivalent to actual independence of X_j. The converse is immediate. Hence the matrix pattern in (1.5) is exactly the effective dependency graph of the POLYNOMIAL outputs, not a chosen presentation of their coefficients. This proves the equivalence.

Finally expand F_i in X_i. Since sigma_i and X_i are coprime independent polynomial variables, (1.1) forces sigma_i^2 to divide every positive-X_i-degree coefficient of F_i. Integrating H_i coefficientwise, using characteristic zero, gives exactly (1.2). Conversely (1.2) satisfies (1.1).

All identities hold in the function field of the entire affine source; equality there proves the polynomial Jacobian identity globally. No inverse-chart restriction is made to obtain a Keller map.

## 2. UAROOTFACT1: every allowed graph gives actual atomic stages

Use the graph convention above: i->j means block i requires the original root X_j. Choose a vertex with NO INCOMING edges, process it first, and delete it. Such a vertex exists in an acyclic graph. Its own output may require roots of still unprocessed blocks, which are still available. No remaining output needs the now-lost original root of the processed block. Repeat.

At every stage, an unprocessed signal sigma_j is computed from its original source triple; a processed signal is its first current coordinate. Every required other root belongs to an unprocessed block. Formula (1.2) therefore defines an actual polynomial whole-affine stage, retaining every other current coordinate. It has determinant one by the one-block case of UAROOTDAG1. Its f_i0,h_i0 terms can be removed by polynomial target translations in the retained signal and identity parameters.

Each stage has generic degree d_i>=3 and arithmetic/geometric monodromy S_(d_i). To see this, view all its other current coordinates as identity parameters and use the source reconstruction (1.3). After the target translations, its inverse equation is

    h_i(S)+b*S^2-c/2=0,
    h_i(S)=S+sigma_i*S^3
       -sum_(n>=1) (2/(n+2))*b_in*sigma_i^(n+2)*S^(n+2).      (2.1)

Its cubic coefficient is sigma_i*(1-(2/3)*b_i1*sigma_i^2), which cannot vanish identically because b_i1 is POLYNOMIAL in sigma_i. Hence d_i=deg_S h_i>=3 and h_i'(0)=1. UAQTILT1 applies over the field of identity parameters and sigma_i. The primitive-monodromy/birational-etale argument of UAALLATOM1 then proves absolute stable self-atomicity of this WHOLE polynomial stage. Allowing n=1 extends the native coefficient list to index3; the general quadratic-tilt theorem, rather than a frozen cubic coefficient assumption, is the needed input.

Consequently

    T=stage_N composed ... composed stage_1,
    deg T=product_i d_i,

with a displayed chain of N absolute stable atoms. In particular N>=2 is never absolutely stably self-atomic. We do not classify every possible maximal chain or claim a product monodromy group for source-root couplings. No uniqueness or classification of other chains is asserted.

More generally, every subset J with NO edge from outside J into J has a polynomial partial-processing model and a polynomial outer factor. These subsets form the order-ideal lattice of the dependency partial order. This supplies intermediate models; it is not an assertion that every intermediate field is in that lattice. The edgeless case recovers the Boolean partial-processing mechanism of UAPARALLEL1. Nonzero acyclic root dependencies are additional source-dependent, generally non-descended examples, but they remain factorized.

Multiplying these corrections by polynomials in retained signals, including arbitrary powers of the marked divisors in UAMARKEDDEF1, does not escape (0.1). Whenever the resulting map is Keller, the same theorem applies. Agreement to arbitrarily high finite order at marked planes therefore cannot rescue atomicity INSIDE this class. Arbitrary vanishing corrections are not asserted Keller.

### A genuinely non-descended rank-two example is still factorized

Take two blocks P,U. Give the P block constant native profile a_4=1, and the U block native profile a_4=X_P. The graph has the single edge U->P. Relative to auxiliary G=K(U), its last-two-output return selector is beta=K(U)-Phi_4(U;X_P), with zero first component. It is not in k[P,K(U)]: hold P=(1,0,1), so X_P=1, and compare U=(0,0,-1) and U=(-1/2,3,26). Their K(U) values both equal(-1,0,0), while their profiled auxiliary outputs differ in both last coordinates. Thus this selector really escapes the descended-time hypothesis UABASICTIME1. It still factors: process the U profile while P is retained, then process P. The two displayed absolute stable atomic degrees are4 and4. This theorem closes more than the basic time-map class, without claiming to classify every non-descended selector.

### Native root profiles are a proper included subclass

Allow a_ij to be any polynomial in u,sigma,X in the old native corrections

    F_i=-(sigma_i^2/2)*sum_(j>=4) j*a_ij*X_i^(j-2),
    H_i=-sum_(j>=4) (j-2)*a_ij*X_i^j.                         (2.2)

They satisfy the self equation exactly when

    sum_(j>=4) X_i^j*partial_(X_i) a_ij=0.                   (2.3)

Together with effective acyclicity, this is the full Keller criterion in this subclass. Apparent self-dependent coefficient lists can cancel; for example a_i4=X_i^2, a_i5=-2*X_i, a_i6=1 give ZERO corrections. Use the collected outputs, or the equivalent effective potential, before declaring a cycle.

## 3. UAROOTSYM1: equivariant orbit factors and the transitive fibre obstruction

Suppose a finite group Gamma permutes the N source/target triples, acts compatibly on the identity parameters u, and (0.1) is Gamma-equivariant. It then preserves the effective dependency graph. Define a vertex's height as the length of its longest outgoing directed path. Heights are Gamma-invariant and strictly decrease along edges. Thus no edge joins vertices in one Gamma-orbit, and the quotient graph on block orbits is acyclic.

Choose a topological order of these orbits, again processing an orbit having no incoming edge from remaining orbits. Within that orbit there are no cross-root dependencies. Process its blocks simultaneously, retaining all other coordinates and original signals as above. Each orbit stage is a Gamma-equivariant polynomial Keller self-map. This gives actual equivariant orbit factors, not merely non-equivariant factors exposed after forgetting symmetry. Different choices between incomparable orbits are allowed; a naturally chosen ordering or a universal arbitrary-base trivialization is not asserted.

If Gamma acts TRANSITIVELY on ALL the nonidentity blocks, the graph has no edges at all. After self normalization (1.2), only the retained signals and parameters can enter another block. For a Gamma-fixed target, identity parameters and signals are fixed. All block equations are then copies of one three-variable equation. Therefore the COMPLETE scheme fibre is

    T0^N,                                                     (3.1)

with the permutation action. Over an algebraically closed field its rank is m^N, where m is the finite rank of T0, and every nonempty fibre contains a Gamma-fixed diagonal point. For N>1 it cannot be the natural transitive N-label Gamma-set: that set has no global fixed label, and N is not an N-th power. This is a complete-fibre statement, not just a selected subconfiguration count.

This extends the UAPARALLEL1 diagonal obstruction from retained-output coefficients to ALL root-polynomial corrections (0.1). The actual canonical trace receiver UACANFIB1 has an EXTRA auxiliary orbit and shifts its first output; it does NOT meet the hypotheses of (3.1). With extra orbits we claim the equivariant orbit factorization above, not a universal diagonal fixed point in every nonempty fibre. No contradiction with the constructed canonical receiver follows.

## 4. Exact boundaries and the next gate

Ordinary nilpotence alone cannot replace (1.5). For example

    A=((p*q,p^2),(-q^2,-p*q)), A^2=0,
    det(diag(r1,r2)+A)-r1*r2=p*q*(r2-r1).

The two independent diagonal variables detect the cycle that a scalar characteristic polynomial misses. This is why the cyclic nilpotent time families from the [higher-return audit](../universal-keller-higher-return/LITERATURE_AND_GAPS.md) do not refute this theorem: their ordinary time variables and descent hypotheses are different.

The characteristic-zero field/domain scope is essential. In odd characteristic p, native coefficients a_14=X_2^p and a_24=X_1^p give a Jacobian-one polynomial map with a genuine two-way polynomial dependency: the cross derivatives vanish. Over the nonreduced characteristic-zero base k[epsilon]/(epsilon^2), coefficients a_14=epsilon*X_2 and a_24=epsilon*X_1 also give determinant one, because the two-cycle product is zero. These are boundaries of the proof, not atomicity or canonical-fibre certificates. The earlier universal marked-fibre deformation identities remain valid over nonreduced bases in their stated scope.

A remaining characteristic-zero atomic canonical candidate must leave both the descended-time class UABASICTIME1 and this retained-signal root-polynomial class. It may change a first output (as the existing trace filter already does), involve further source channels beyond k[u,sigma,X], or change the compiler/auxiliary geometry. Such departures still require a whole-source determinant identity, the exact complete marked scheme, effective gluing and a separate absolute-atomicity proof. No bounded-degree sweep is left to perform inside (0.1).
