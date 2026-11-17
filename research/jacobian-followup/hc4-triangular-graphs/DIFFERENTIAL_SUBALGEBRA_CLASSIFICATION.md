# Differential subalgebras of an actual plane Keller source

Claims HC4MYKV0/1. The full target remains HC4/JC2. This note completes
an all-degree classification INSIDE AN ACTUAL polynomial plane Keller
source. It does not construct such a source for a candidate cover or prove
the last source-to-trace-regularity implication.

## 1. Statement: a finite Boolean family over each intermediate field

Let B=C[P,Q] subset R=C[x,y] be an ACTUAL complex plane Keller inclusion,
F=Frac(B), L=Frac(R), and W_B the target Weyl algebra acting on R.
For ANY intermediate field F subset E subset L, let S_E be the full finite
normalization of B in E and put

    Lambda_E={q in E:Tr_(E/F)(q*S_E) subset B},
    A_E=B[Lambda_E], U_E=Spec(A_E).

KU0/1 proves that U_E is the entire smooth affine etale locus of S_E/B
and that A_E subset R. The actual source map X=Spec(R)->U_E is etale.
Let V_E be its OPEN IMAGE. Let C_(E,1),...,C_(E,b_E) be the distinct
codimension-one components of U_E minus V_E. These are prime curves
in the SMOOTH surface U_E, not ramified primes of S_E/B. There may also
be missing isolated points. No disjointness, smoothness or residue-degree
one assertion for these intermediate curves is needed below.

For J subset {1,...,b_E}, define

    U_(E,J)=U_E minus union_(j in J) C_(E,j),
    C_(E,J)=Gamma(U_(E,J),O).

Then the COMPLETE classification is

    {C : B subset C subset R, W_B*C subset C, Frac(C)=E}
          = {C_(E,J) : J subset {1,...,b_E}}.          (1)

No finite-generation, normality or smoothness hypothesis on C is imposed.
Every such C is AUTOMATICALLY normal, finitely generated, smooth affine
and etale over B. The family is a Boolean lattice with exactly2^b_E
members. Its minimum is A_E, and its maximum is R intersect E. Increasing
J ALLOWS poles on more missing curves and makes the ring larger.

There are only finitely many intermediate fields of the finite separable
extension L/F. Thus there are only FINITELY MANY target-differential
B-subalgebras of the fixed actual R, even allowing different fraction
fields. This is conditional on the actual source and its KQ0 layer theorem;
it is not a classification of arbitrary differential algebras in rational
function fields or of subalgebras without target-derivation stability.

## 2. The maximal ring in a field is an affine etale open

First

    R intersect E=Gamma(V_E,O).                       (2)

One inclusion is pullback. Conversely, if a rational function in E pulls
back to a regular function on X, it is regular on V_E: X->V_E is an
etale surjection, so regularity descends faithfully flat locally. Equally,
every height-one local ring of the normal V_E has an etale source lift,
and its valuation cannot be negative. This argument identifies functions
in the SPECIFIED embedded field E; it is not invariant-theoretic descent
for an arbitrary source action.

The complement U_E minus V_E has finitely many irreducible components.
After removing its prime curves, only finitely many isolated points can
remain missing. Normality of the smooth surface extends functions over
these codimension-two points. Consequently

    R_E:=R intersect E=C_(E,{1,...,b_E}).              (3)

Every finite union of the curves C_(E,j) is an effective Cartier divisor
on U_E, since U_E is regular. Its complement is affine. To see this
without assuming the divisor principal, cover the affine U_E by principal
affines on which its ideal is invertible with one generator. On each the
complement is a principal affine. The open immersion is therefore an
affine morphism, since affineness is local on the base; its source is
affine. The open is of finite type over C and etale over B.

Thus all the displayed C_(E,J) are finitely generated normal etale
B-subalgebras of R_E, stable under the uniquely extended target
derivations, with fraction field E. This also proves these properties
for R_E. This geometry does NOT assume that E has an affine-plane model.

## 3. The minimal differential algebra is the canonical trace-dual algebra

Let I_L=Soc_W(R) and T_L=R/I_L. KQ0 proves that both I_L and T_L are
SEMISIMPLE, I_L consists of full-support IC factors, and T_L consists
only of curve-support IC factors with NO point-support factor. Let I_E
be the IC extension of the generic permutation subrepresentation E
inside I_L. KU2 proves

    A_E intersect I_L=I_E, I_E*I_E=A_E.              (4)

Now let C be ANY algebra on the left of (1). Its generic B-localization
is E: B-localizing an algebraic domain over the field F gives a field,
so C tensor_B F=Frac(C). The W-submodule C intersect I_L has this same
generic subspace E. Indeed R/I_L is B-torsion, so intersecting with I_L
does not change that generic subspace. Semisimplicity of I_L therefore
gives

    C intersect I_L=I_E.

As C is a RING, (4) forces A_E subset C. Hence

    A_E subset C subset R_E.                        (5)

Moreover R_E intersect I_L=I_E by the same argument. Therefore
R_E/I_E embeds in T_L, and

    Q_E:=R_E/A_E

is SEMISIMPLE with only curve IC composition factors. This is the crucial
transfer through the ACTUAL R. We do not apply the plane boundary theorem
directly to the intermediate normalization, and we do not silently discard
point-supported extension terms for an unrelated affine surface.

## 4. Why multiplication forbids mixed or partial boundary choices

The W-submodule M=C/A_E of Q_E is also an A_E-submodule because C is
a ring containing A_E. Analyze it at the generic point of any target
divisor D. Complete in a normal parameter u and split the finite etale
residue extensions. Equivalently keep each finite residue field K_j
over K=C(D). On the unramified factors of U_E the completed coefficient
algebra is a product of DVRs K_j[[u]]. On a missing curve j, its quotient
in Q_E is the principal-parts module

    H_j=K_j((u))/K_j[[u]].                          (6)

Other factors give zero: retained unramified curves already have regular
functions in both rings, and ramified factors were already removed in
U_E and contribute no difference in (5). Formula(3) supplies EXACTLY
the principal parts (6), with no unidentified generic lattice.

The completed A_E action contains the projectors to its distinct prime
factors and the ENTIRE coefficient field K_j on each. Thus even if two
curves have the SAME target image and isomorphic W-character modules,
a diagonal or a nontrivial slope between their principal parts cannot be
an A_E-submodule. Likewise a proper subrepresentation of a residue-field
factor cannot be chosen just because its target local system decomposes:
multiplication by K_j has no nonzero proper K_j-submodule of its own
one-dimensional vector space.

Finally H_j is SIMPLE under K_j[[u]] and the normal derivative. Any
nonzero principal part has a leading pole u^-m times a unit. Multiplying
by an appropriate regular unit and u^(m-1) recovers u^-1. Successive
normal derivatives give

    (d/du)^k(u^-1)=(-1)^k*k!*u^(-k-1),             (7)

so all principal parts are generated in characteristic zero. Completion
and residue splitting are faithfully flat, and the normal derivative is
a rational target-linear combination of the original derivatives locally
at the generic divisor. This proves that M selects either ALL or NONE
of each prime-curve module (6). The distinction is a prime of U_E,
not just a target image or an abstract character summand.

We obtain a UNIQUE subset J of the missing curves. The already constructed
C_(E,J)/A_E has exactly these same generic submodules. Both it and M are
submodules of the semisimple Q_E. Their sum modulo their intersection has
zero generic restriction on every curve, so could have only point support.
But every subquotient of Q_E is semisimple with curve IC composition
factors, and hence has no nonzero point-supported subquotient. Therefore
the two submodules are EQUAL. Their inverse images in R_E give
C=C_(E,J), proving (1) globally, INCLUDING all special target points.

Distinct J give distinct rings by their generic allowed poles. Intersections
correspond to intersection of J; the algebra generated by two members is
W-stable and permits exactly the union of their allowed poles, so it is
the union member. This completes the Boolean-lattice assertion.

## 5. Derivative algebras and one generic linear source generator

For ANY finite family Sigma of actual source polynomials, put E=F(Sigma)
and let D_B(Sigma) be the B-algebra generated by Sigma and ALL their target
partial derivatives. The field E is stable under these derivatives because
it is algebraic and separable over F. Its fraction field is E. Applying
(1) shows the EXACT formula

    D_B(Sigma)=C_(E,J(Sigma)),
    J(Sigma)={j : some s in Sigma has v_(C_(E,j))(s)<0}. (8)

Indeed a source polynomial lies in C_(E,J) exactly when none of its
unramified finite poles occurs outside J. There are no further poles on
U_E: R_E is its specified open function ring. Each C_(E,J) is already
differentially stable. Thus the smallest such algebra containing Sigma
has exactly the union of those allowed poles. Normalization, finite
generation and smoothness of D_B(Sigma) are CONSEQUENCES, not hypotheses.

In particular, if t is ANY integral primitive element of S_E/B, then

    D_B(t)=A_E.                                    (9)

The primitive ORDER B[t] may be nonnormal and singular. Equation(9)
does not require its trace-dual to agree with that of S_E: it differentiates
inside the ACTUAL source and uses (1). This contrasts with KU0's order
recipe, which correctly returns the ORDER'S own etale locus. No general
order-to-normalization theorem outside the present source is claimed.

For E=L, write b_0=b_L for the number of omitted unramified primes. By
KU1 these are exactly the curves of U_L minus X. Localization from the
actual X=A2 gives

    Pic(C_(L,J))=Z^(b_0-|J|), units=C*,
    C_(L,J) is factorial iff J is ALL iff C_(L,J)=R. (10)

These are actual-source completion statements. They do not characterize
arbitrary factorial etale surfaces, as the credited Kulikov control shows.

Choose t_lambda=x+lambda*y. For ALL BUT FINITELY MANY lambda in C,

    F(t_lambda)=L and D_B(t_lambda)=R.              (11)

For primitivity, take the d embeddings of L into its normal closure.
For each pair, equality of their values on x+lambda*y excludes at most
one lambda: equality for two lambdas would make them agree on x,y and
therefore on L. At most d*(d-1)/2 constants are excluded. At each of the
b_0 missing unramified primes, at least one of x,y has a pole; otherwise
their regular values would give a centre in A2_source for its DVR.
The induced map to the separated finite normalization agrees generically
with its existing centre, so that centre would be in the source open,
contrary to the curve being omitted. Cancellation of a pole in
x+lambda*y excludes at most one additional lambda per prime. The remaining
lambda are primitive and polar on every such curve, so (8) gives (11).
Thus d*(d-1)/2+b_0 is a valid upper bound for this sufficient exceptional
set; no sharpness is claimed. No coefficient-depth search is required.

## 6. Literature audit and the remaining full theorem

Makar-Limanov's [MPIM22-11, A conjecture on D algebras](https://pure.mpg.de/rest/items/item_3658498_1/component/file_3658499/content)
defines precisely the algebra generated by an algebraic function and all
its target derivatives. Its unrestricted nonembedding statement is a
CONJECTURE, not a theorem available to finish HC4. Its Section2 proves
the nonsingular primitive-hypersurface case by recovering the reciprocal
of P_t. It separately conjectures that a suitable nonsingular algebraic
element exists in D(t). Our automatic SMOOTHNESS OF THE ENTIRE ALGEBRA
in (8) does NOT supply one globally nonsingular primitive hypersurface
element. Smooth affine surfaces need not have that monogenic presentation.
The retained KN1/HY2 obstruction already rules out global CI/monogenic
finite normalizations for a nontrivial actual embedded intermediate.
The distinction remains essential after (9).

For clarity, the nonsingular-order argument is elementary: on P=0,
P_i+P_t*delta_i(t)=0. If P and all its partials generate1 in the ambient
polynomial ring, dividing that identity by P_t expresses 1/P_t in D(t).
If P_t is nonconstant as an element, this produces a nonconstant unit,
impossible inside R. In the finite etale case P_t can instead be a scalar
unit; complex affine-plane finite-cover triviality handles that case.
One cannot assume the nonsingular-order premise merely because D(t) is
smooth, etale and finitely generated. No unproved statement from the
preprint is imported here; no claim of priority is made for its definition
or its proved special cases. Its all-dimensional conjecture is also not
a substitute for the already retained higher-dimensional countercontrols.

Adamus--Crespo--Hajto's [SIGMA2019 Theorem3.1](https://doi.org/10.3842/SIGMA.2019.034)
uses SURJECTIVITY of a particular differential tensor-constant map to
obtain strong normality, then the Galois case. Finite generation in (8)
does not prove that surjectivity or Galoisness. This is context, not an
additional imported theorem or an automatic descent of source actions.

FULL THEOREM TARGET: prove full HC4/JC2 or construct an ACTUAL counterexample.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.          OPEN

Classification (1) and finite generation (8) are complete; this implication
is not. It would give Lambda_L=S_L, finite etale normalization and inversion,
as KU0/1 explains. The alternative full marked canonical-energy integrality
gate also remains open. The useful next input must constrain the actual
source embedding or marked multiplication, beyond classifying its available
differential subalgebras. External review and outreach remain DEFERRED;
no assurance flag is upgraded by these internal proofs and exact controls.
