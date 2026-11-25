# A finite presentation-to-witness algorithm for compact rational homogeneous products

The principal-face and algebraic density reconstruction theorems imply an effective construction once their exact cell and moment data are known. This proof computes those data from the promised rational Hopf or matrix presentations below, then produces the actual rational all-order counterexample. The result concerns existence and termination of an algorithm; no software implementation or complexity bound is asserted.

## 1. Effective input and output

An effective affine-group presentation means a finitely presented rational coordinate Hopf algebra A=Q[z_1,...,z_s]/I, with explicit polynomial formulas for comultiplication, counit and antipode, and a finitely generated Hopf ideal J specifying a closed subgroup H. Equivalently one may start with defining ideals for actual matrix groups H subset G subset GL_m; determinant inverse is included as an affine coordinate. Ideals and structure maps, not a list of generators of an abstract group of rational points, are the input.

PROMISE: G is connected semisimple, G(R) is compact, H0 is reductive, and h=Lie(H) is not an ideal of g=Lie(G). The procedure is required to halt only on this promised class; it does not decide that classification or promise validity.

There is a terminating exact procedure producing integers L>=1 and e>=1 and an explicitly represented rational right-H-invariant Laurent polynomial

    F in O(G)[t,t^-1],       U_marker=t^e,

which defines an actual regular function on (G/H) times Gm. Every pure moment of F is zero and EVERY positive-order U_marker-marked moment is a nonzero positive rational number under the native normalized trace/constant term. The same rational pair remains pure/unsafe over every characteristic-zero scalar extension. F is absolutely nonpower even up to a unit and has primitive literal Laurent support. Over every formally real extension both factors are individually Mathieu and their product is not.

The output format is an ambient right-H-invariant polynomial representing the quotient function. Computing a generating presentation of O(G/H) is unnecessary. If a quotient coordinate presentation and pullback are also supplied, one can recover an expression in those generators by rational polynomial enumeration and ideal equality, which terminates by the fact that the supplied quotient generators generate O(G/H).

The proof uses the established NIP geometry and SPT/PDA moment arguments, and only classical effective real quantifier elimination (QE), polynomial-ideal normal forms, and exact algebraic-number operations as algorithmic primitives.

## 2. From Hopf data to a matrix embedding

Here is a direct finite construction if a matrix model is not initially supplied. Compute a Groebner basis of I. For each algebra generator f, reduce its coproduct and write

    Delta(f)=sum_i a_i tensor b_i

where the b_i are distinct standard monomials, hence linearly independent in A. Coassociativity shows that V_f=span_Q{a_i} is a finite coideal with Delta(V_f) subset V_f tensor A: apply linear functionals selecting the b_i in the third tensor factor to (Delta tensor id)Delta(f)=(id tensor Delta)Delta(f). The counit gives f=sum_i a_i epsilon(b_i), so V_f contains f. Include f=1 and take V=sum_f V_f. All its elements are explicit rational polynomials and its rational basis is computed by finite normal-form linear algebra.

Express Delta(v_i)=sum_j v_j tensor a_ji in that basis. These coefficients give a rational representation of G by right translation on V, with inverse matrix supplied by the antipode. The coefficient functions a_ji generate A as an algebra: applying the counit to the first tensor factor recovers every v_i and in particular every original algebra generator. Thus the representation map G->GL(V) is a CLOSED immersion, not just a pointwise faithful map. Its defining ideal and the image ideal of H are computed by elimination from their explicit coefficient formulas. This reduces the general effective Hopf input to the matrix input by finite rational computations.

For completeness the ideal operations used throughout are effective in arbitrary dimension. Polynomial reduction terminates for a fixed admissible monomial order, meaning a multiplicative well-order. Repeatedly add nonzero reduced S-polynomials; each addition strictly enlarges the leading-monomial ideal. Dickson's lemma prevents infinitely many such additions: every infinite sequence in N^m has a coordinatewise increasing pair, by taking a nondecreasing subsequence in one coordinate and inducting on m. Hence completion terminates. Buchberger's criterion then gives a Groebner basis and membership is equivalent to zero normal form. The necessary criterion and normal-form facts, not a zero-dimensional restriction, are the scope of [B] below. The source writes degree-lex order; its criterion proof uses well-foundedness and multiplicativity and therefore applies also to the admissible block-elimination orders used here.

## 3. Rational Lie input and an actual nonzero moment projection

Differentiate the defining ideals at the identity to compute rational bases of g and h. Matrix commutators give their rational brackets. B=-Killing_g is rational and positive definite on the actual compact real Lie algebra; compute it from the adjoint matrices and let p=h-perp. The nonideal promise implies [h,p]!=0. Indeed [h,h] subset h, while B-invariance gives [h,p] subset p, and if the latter were0 then h would be an ideal.

Choose basis vectors a in h and x in p with b=[x,a]!=0. This is a finite exact search through basis pairs. For z(s)=Ad(exp(sb))x,

    d/ds B(a,projection_h z(s)) at0 = B(a,[b,x])=B(b,b)>0.

Consequently x is nonzero rational, projection_h x=0, and its compact adjoint orbit M has a nonconstant, nonzero H-moment projection. Compute 2D=rank(ad_x) by rational linear algebra. All of this is the explicit positive-input direction of the NIP theorem.

Every algebraic formula below is now effective over the real algebraic numbers ALG. One can use standard ALG arithmetic and QE. Alternatively encode each ALG number by a rational polynomial and an isolating rational interval, with uniqueness decided by QE. Membership, equality and order then reduce to rational QE sentences. Enumerating these representations enumerates ALG. A nonempty semialgebraic set defined over ALG contains an ALG point, by QE/transfer between the real closed fields ALG and R; enumeration with exact membership therefore finds a point and terminates. This inefficient alternative is enough for every sampling step below.

## 4. A rational Cartan and the full geometric finite action

For a variable a in h, compute the maximum rank of ad_a using its polynomial minors. A rational a attaining that maximum is found by enumerating rational tuples until the rank is attained. It exists because the maximum-rank locus is a nonempty rational Zariski-open set. Compactness makes every real element semisimple; maximal rank therefore means regular semisimple. Its Lie centralizer t=ker(ad_a|h) is a rational Cartan, the Lie algebra of an actual maximal Q-torus T_H. The abelian case has maximum rank0 and t=h. No torus classification or rational diagonalizing group element is required.

Simultaneously diagonalize the commuting adjoint t-matrices over the algebraic numbers to obtain their finite root list, or enumerate and verify simultaneous algebraic eigen-data. Compact root values become real algebraic linear forms after division by i. The rational regular a chooses a positive root system by its root signs. The closed chamber C is given by those finite algebraic linear inequalities. Connected Weyl reflections and their finite closure compute W0, acting trivially on the centre.

The full normalizer of t in H is described by polynomial equations in the matrix group. Its action on t gives the semialgebraic set

    W_real={R : exists h in H(R), Ad(h)t=t and Ad(h)|t=R}.

QE computes this set. It is FINITE and equals the full geometric Weyl matrix image: the PDA compact-component polar lemma supplies real representatives for all geometric normalizer components. The matrix image, rather than an unjustified rational representative, is the input required by the circuit. Enumerate its algebraic points, removing each point found and testing whether any remain; finiteness gives termination. Equivalently enumerate finite ALG lists and test equality with W_real by QE. The set is Galois-stable because the algebraic normalizer and rational Cartan are defined over Q.

Compute Omega={w in W:wC=C} by finite QE containment/equality tests. W0 and Omega have the chamber decomposition established in PDA. B restricts to a positive rational W-invariant quadratic form on the rational coordinate space t.

## 5. The principal face and its actual algebraic linear support

The actual radial image is exactly the following semialgebraic set:

    R={u in C : exists g in G(R), projection_h(Ad(g^-1)x)=u}.       (1)

Indeed the H0-action makes the full moment image conjugation-invariant, and each compact h-orbit meets C. No H0-parametrization is needed in (1). Inversion is a regular polynomial in matrix coordinates and determinant inverse, so this is an effective ALG formula. The set is compact and contains0.

Enumerate the finite root faces of C. QE tests which face sigma has R subset closure(sigma) and R intersect sigma nonempty; the NIP principal-face theorem guarantees the unique such principal face. Equivalently choose the smallest root-face closure containing R and verify its principal intersection. These are finite linear-root tests.

Compute A=span_R R: start with the zero span; ask by QE whether (1) contains a point outside the current span; if yes choose an ALG such point and enlarge the span. At most dim(t) such enlargements occur. The NIP linearity argument proves that this is exactly the effective centre moment flat. Its ALG basis is thus computed, n=dim A>=1, and R subset A subset span(sigma). Neither Q-definition nor full-W-stability of A is asserted. Omega preserves sigma and A; this also follows from the genuine actual symmetry proved in PDA.

Choose an ALG eta in sigma and compute l=C_h(eta), which is constant on this face, together with its centre s. These are finite bracket-kernel computations. The effective torus kernel may be recovered as the B-orthogonal complement of A inside s, but no defining equations for that quotient torus are needed in the procedure.

## 6. Exact critical-value formula: the clean tangent uses l

Write mu(z)=projection_h z and Y={z in M:mu(z) in sigma}. The principal cross-section is locally the preimage of the open l*-slice. KKS transversality across h/l therefore gives

    T_zY={ [a,z] : a in g, projection_(l-perp in h)[a,z]=0 }.       (2)

One must NOT obtain (2) by naively differentiating all constraints into the centre s. Although the actual Y-moment image lies in sigma subset s, those extra central restrictions are consequences of the principal-face geometry. Clean transversality is to l, whose orthogonal complement consists of the noncentral KKS root directions.

Let C_z:g->l-perp be a->projection_lperp[a,z], and E_z:g->A be a->projection_A[a,z]. The first projection may be composed with projection_h without changing it. These matrices have ALG entries polynomial in z, and

    rank(C_z)=dim h-dim l,
    rank(dmu_Y at z)=rank([C_z;E_z])-rank(C_z).                   (3)

The second identity is the elementary formula for the rank of E_z restricted to ker(C_z); ker(ad_z) is already killed by E_z. Formula (2), rather than a possibly nonreduced equation tangent, identifies this linear space with the actual smooth tangent.

The critical-value set inside sigma intersect A is consequently the effective semialgebraic set

    Crit={u in A intersect sigma : exists g in G(R),
          z=Ad(g^-1)x, mu(z)=u,
          rank([C_z;E_z]) < dim(h)-dim(l)+n }.                (4)

Use the vanishing of all maximal minors to express the last inequality. This is an exact existential formula over ALG. The target may be lower dimensional in t and the subgroup moment may be singular as an h-valued map; neither invalidates (2)-(4).

## 7. Finding a finite wall cover by decidable enumeration

The NIP finite-weight/centralizer argument proves that (4) is contained in finitely many REALALG affine hyperplanes which are genuinely nonconstant on A. We need not compute their minimal set, nor take a Zariski closure.

Enumerate finite families of proper affine hyperplanes in the computed ALG coordinate space A. For each family test by QE

    for every u, u in Crit implies u lies in their union.       (5)

Each test terminates. Some finite ALG family passes by the NIP existence theorem, so the enumeration terminates. Extra walls are harmless. Include the proper restrictions to A of the root walls bounding sigma; discard equations identically zero on A.

Find an integer R_box with R strictly inside the coordinate box |u_i|<R_box by QE and increasing integer search. Compactness guarantees success. Add its bounding hyperplanes. Close this finite arrangement under Omega and intersect the resulting translated boxes; the support remains strictly inside the resulting bounded Omega-invariant polytope, because the actual R is Omega-invariant.

Every nonempty full-dimensional arrangement sign cell inside sigma and the box is an explicitly given bounded open convex polyhedron with ALG data. Their finite list is found by sign enumeration and QE feasibility. The local proper finite-stabilizer DH theorem in NIP makes the centre density polynomial on each such regular connected chamber. The root-Pfaffian induction factor has degree q=(dim h-dim l)/2, so the actual radial density has degree at most

    b=D-n.                                                     (6)

No CDKW derivative/regular-orbit-support hypothesis enters this argument. Compactness and the NIP full-measure argument exclude additional wall masses; the radial probability is A-absolutely continuous. Each chamber has either empty fibre and density identically0 or positive reduced volume throughout. If necessary, keeping every cell with merely nonnegative polynomial density is also allowed at this stage.

Thus we have COMPUTED the finite ALG cell and degree inputs that PDA formerly accepted as supplied data. Faces can be appended with zero actual measures to form a compact cell complex. Omega degree-orbit maxima are finite; the common bound (6) already supplies one compatible bound.

## 8. Exact full-H Chevalley lifting

Let P in Q[t]^W be a polynomial of degree d, treating homogeneous parts separately. Write a general homogeneous polynomial Q_d on h with unknown rational coefficients. Impose the following finite rational linear conditions:

1. Q_d restricted to the rational Cartan t equals the corresponding part of P;
2. Q_d(Ad(h)Z)=Q_d(Z) as a polynomial identity in O(H)[Z].

The second identity uses the FULL H coaction, including every disconnected component. Expand it and reduce by the defining ideal of H; the remaining coefficients are rational linear expressions in the unknowns. Groebner normal forms therefore give a finite rational linear system. The established full-H rational Chevalley restriction theorem guarantees a solution in the same homogeneous degree. Rational Gaussian elimination produces one. No averaging over hypothetical rational component representatives is involved.

Compose the lift with projection_h Ad(g^-1)x. This gives an explicitly represented regular rational function F_P on G, right-H-invariant. Its unique descended function on G/H is theta(P). This construction evaluates every native invariant test required by PDA without computing O(G/H) generators.

## 9. The native rational moment oracle is finite linear algebra

Given F_P in O(G), write it as a polynomial in the matrix entries and d=(det g)^-1. Choose bounds e,s for its entry degree and its power of d, and include1. The finite span of all entry monomials of degree at most e times d^j, 0<=j<=s, maps to a finite-dimensional LEFT-G-stable rational subspace V of O(G): left translation replaces entries linearly and multiplies d by a scalar character. Groebner reduction and rational linear algebra compute a basis of this image, even if normal forms themselves have larger degree.

For a rational basis a of g, compute the left infinitesimal derivatives D_a on this finite basis. They preserve V and have rational matrices. Solve for the linear functional ell on V satisfying

    ell(1)=1,       ell(D_a v)=0 for every a and v in the bases.    (7)

There is exactly one solution. Characteristic-zero semisimplicity of the rational G-module gives V/gV isomorphic to V^G, and V^G=Q because left translation is transitive on G and1 belongs to V. The native normalized group trace satisfies (7), hence is this rational solution. Its restriction to right-H-invariant functions is precisely the native G/H trace.

This is an exact algorithm for the rational value lambda(theta(P)). The finite module is constructed from the input polynomial, not assumed to be supplied externally. It also gives a moment oracle at ANY requested order for the eventual output witness. It is not a numerical Haar quadrature or a finite-prefix test.

## 10. Recovering the actual densities

The real-component polar lemma in PDA supplies genuine Omega symmetry of the actual radial law; it is valid for this full H and the W computed in Section4. Native moments from Sections8-9 are rational. The computed ALG cells and bound (6) therefore meet the PDA finite-dimensional reconstruction theorem.

Construct an ALG basis for the Omega-invariant signed polynomial cell-measure space. Exact polynomial integrals over the bounded ALG polytopes are computable: enumerate vertices by linear systems and use a pulling triangulation, then the usual simplex monomial integral and ALG Gram determinant for relative volume. Enumerate rational monomials and apply the finite W-Reynolds operator to get rational invariant tests. Galois stability makes their coefficients rational; exact algebraic arithmetic recognizes them as such. Invertible moment rows are found by retaining independent rows, and existence/injectivity in PDA guarantees finite termination. Query their moments using (7), solve the ALG system and recover every density coefficient exactly.

Identically zero cells are discarded by exact coefficient comparison. The NIP regular-fibre argument makes every retained density positive on its open chamber. The actual0 support, noatom0 and positive-dimensional Origin properties are those already proved in NIP for the selected x. No transcendental volume constant needs to be supplied: the native normalization and moment equations recover the actual probability coefficients.

## 11. An entirely finite choice of the rational circuit

The supporting relative hyperplanes of the retained cells are now ALG. Through-origin facets need no zero because their Euler flux vanishes. For all other facets, the rational invariant face-polynomial existence theorem in SPT gives an f_base in Q[t]^W with f_base(0)=1 vanishing there. It is COMPUTABLE without any extra number-field algorithm: enumerate rational polynomials; check invariance under the finite ALG W matrices, value1 at0, and identically zero restriction to every listed nonradial ALG facet. Each check is finite algebraic linear substitution, and existence guarantees termination.

Use the rational positive quadratic q=B restricted to t. QE finds a positive rational c with1-c*q positive on the compact cell closures. The exact algebraic density coefficients now allow the finite damping/simplex/tail-bound choices proved in PDA. Choose M so the Away bound is below1/4; empty Away gives M=0. Then find an ALG simplex in a positive Origin cone with |f|>1/2 and positive first effective density; compute A_plus>0. Compute an ALG upper bound C for the finite Away terms, express the densities in the linear A-coordinates centred at0 and determine d_star from their exact first homogeneous pieces, and choose

    K>=max(B_star,2d_star,1),
    C(1+K)^d_star4^-K < A_plus,
    L>=max(B_star,K,1).                                      (8)

All comparisons are exact ALG comparisons; every search terminates for the reasons proved in PDA. This single inequality proves the ENTIRE subsequent positive tail by its monotone ratio bound. It is not an enumeration of moment prefixes.

Apply the homogeneous lifts of Section8 to f((1+t)u) coefficientwise. This gives an actual rational right-H-invariant a(T), T=1+t. Output

    p=T/t * a(T)^2,
    F=p^L(1+p),       U_marker=t^d_star.                     (9)

SPT proves that (9) is pure at every order and the SAME unit marks EVERY positive order nontrivially. The native values are rational, and their proven real positivity makes them nonzero over every characteristic-zero extension. Every operation producing (9) terminates from the promised presentations; no cell description, torus action matrix, or moment oracle remains an external input.

## 12. Direct geometric certificates and the precise conclusion

At the actual identity coset, projection_h x=0 and theta(P)(eH)=P(0). Hence a(T)(eH)=1, p(eH,t)=1+t^-1, and at t=-1/2 the modifier1+p has derivative-4 while p=-1 is a unit. Smooth geometric integrality of the actual affine G/H therefore gives an actual prime divisor of order1 in F. Thus F is not a unit times any proper power after algebraic closure. The same identity specialization is

    F(eH,t)=(1+t^-1)^L(2+t^-1),

whose adjacent coefficients2 and2L+1 prove primitive literal Laurent support. These certificates can be attached to the output independently of a numerical moment computation. Native trace basechange and the rational Gram/scalar Laurent theorems give the field statements in Section1.

This proves a terminating algorithm on every promised effective rational Hopf/matrix input. It supplies no complexity estimate or implemented software and makes no claim outside the promise, in particular no arbitrary finite/ideal-H counterexample classification. The exact INPUT condition h nonideal is used to select x, not asserted as a converse criterion for Mathieu failure. The [principal-face existence theorem](COMPACT_RATIONAL_NONIDEAL_HOMOGENEOUS_PRODUCT_COUNTEREXAMPLES.md) and [conditional density reconstruction theorem](ALGEBRAIC_RECONSTRUCTION_AND_FINITE_AMPLIFICATION_OF_PRINCIPAL_FACE_DENSITIES.md) are separate mathematical inputs; this effective refinement does not change their scopes.

## 13. Classical algorithms and scope

[Q] S. Basu, R. Pollack and M.-F. Roy, [On the combinatorial and algebraic complexity of quantifier elimination](https://doi.org/10.1145/235809.235813), Journal of the ACM43(6)(1996),1002-1045. The [retained author PostScript](sources/basu-pollack-roy-quantifier-elimination-1996-author.ps) is the exact source downloaded from https://www.math.purdue.edu/~sbasu/jacm95.ps. Theorems1.3.1-1.3.3, printed1004-1005, supply effective QE and sentence decision over real closed fields. Section6, printed1041-1044, assembles the algorithm from its earlier sampling, sign-determination and general-decision subroutines. Those classical subroutines remain imported; only existence and termination, not complexity bounds, are used here.

[B] B. Buchberger, A theoretical basis for the reduction of polynomials to canonical forms, ACM SIGSAM Bulletin10(3)(1976),19-29; [retained author PDF](sources/buchberger-canonical-forms-1976-author.pdf), downloaded from https://www3.risc.jku.at/people/buchberger/papers/1976-00-00-A.pdf. Theorem3.3 and Proposition3.7 provide the general-field S-polynomial criterion and zero-normal-form membership equivalence. The source uses degree-lex order; the admissible-order variant and completion termination needed here are explained in Section2. There is no restriction to zero-dimensional ideals.

The normalized quotient trace, full-H rational Chevalley restriction, compact principal-face geometry and its finite critical-wall cover are the scoped inputs from the [actual nonideal homogeneous theorem](COMPACT_RATIONAL_NONIDEAL_HOMOGENEOUS_PRODUCT_COUNTEREXAMPLES.md). The [positive stratified invariant-test criterion](POSITIVE_POLYHEDRAL_STRATIFIED_INVARIANT_TEST_COUNTEREXAMPLES.md) and [algebraic reconstruction theorem](ALGEBRAIC_RECONSTRUCTION_AND_FINITE_AMPLIFICATION_OF_PRINCIPAL_FACE_DENSITIES.md) provide the moment/circuit and finite-tail stages. These scopes include their stated classical inputs and field hypotheses.

Classical computation of DH measures predates this composition. Christandl-Doran-Kousidis-Walter, [arXiv1204.0741v2](https://arxiv.org/abs/1204.0741v2), Algorithm4.4 and Remark4.6 give coadjoint-projection computational methods, while Algorithm4.17 treats projective torus measures. This is a context comparison with its own fixed-point and derivative-principle hypotheses; the present principal-face construction has no CDKW premise. No first-computability or new complex-existence assertion is made.
