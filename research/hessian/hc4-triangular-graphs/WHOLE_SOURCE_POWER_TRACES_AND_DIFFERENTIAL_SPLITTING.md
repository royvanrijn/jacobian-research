# Whole-source power traces and differential splitting

Theorem claims: HC4MYKO0/1/2. Characteristic zero, over C throughout.
These are all-degree criteria and obstruction descriptions, not a proof of
HC4. In particular, regularity of traces on a finite normalization must not
be substituted for regularity of traces on its larger source-open ring.

## 1. Setup and the retained actual-sheet input

Let K:A^n_X -> A^n_Y be polynomial with constant nonzero Jacobian. Write
R=C[x_1,...,x_n], B=C[K_1,...,K_n], F=Frac(B), L=Frac(R), d=[L:F].
Let S be the finite normalization of B in L. The retained normality and
Zariski-Main argument gives S subset R and an open immersion X -> Z=Spec(S).
This is the FULL source field, not an arbitrary intermediate model.

For every target prime divisor P, there is an actual source prime C over P
with ramification index one. Indeed, an irreducible target equation p has
nonconstant polynomial pullback p(K); each irreducible factor defines a
source divisor mapping dominantly to P by quasi-finiteness. Etaleness gives
index one. This is retained KMFI1/HC4MYHW0, not a new surjectivity theorem.
Every r in R is integral at this source DVR. Its image therefore has
nonnegative valuation in at least one conjugate above P. The residue degree
can exceed one; the argument only needs one integral conjugate.

The retained inputs, proof/checker/receipt hashes and the inspected primary
context are bound in [the source context](sources/whole-source-trace-context.json).
No old registered checker is rerun.

## 2. All but the last power trace detect integrality: HC4MYKO0

For r in R, put m=[F(r):F]. Then

    r belongs to S
      iff Tr_(L/F)(r^k) belongs to B for 1 <= k <= m-1.

In particular the uniform cutoff d-1 suffices, even if r is not primitive.
For m=1 the trace list is empty and the assertion follows from the same
actual-sheet argument below (equivalently retained R intersect F = B).

Proof. Set E=F(r). Divide the listed L/F traces by [L:E], which is an
invertible complex integer, to obtain E/F power sums. Newton identities
give all coefficients except the constant term of the monic minimal
polynomial

    f(T)=T^m+c_1 T^(m-1)+...+c_(m-1) T+c_m

in B. At every P choose the actual source DVR from Section 1. All terms
of f(r) except c_m have nonnegative valuation there. Hence c_m also has
nonnegative valuation there. Since c_m is in F, its target P-valuation is
nonnegative. Intersecting the height-one DVRs of the polynomial ring B gives
c_m in B. Thus r is integral over B. This also handles m=1 directly. The
reverse implication is ordinary integrality of conjugates and normality of
B. It uses no finite-free basis for S and no monogenicity of S.

There is a useful local refinement. Fix P and extend its valuation to a
splitting field of L/F. Count the negative-valuation conjugates of r with
embedding multiplicity, say q. At least one conjugate is integral, so
q <= d-1. If q>0, some trace Tr(r^k), 1<=k<=q, has a pole at P. Otherwise
subtract the integral conjugates' power sums from the integral total sums.
Newton identities make all elementary symmetric functions of the negative
block integral. Its product has strictly negative valuation, a contradiction.
This proof allows unequal pole orders and ramification and needs no choice
of a primitive element. Characteristic zero makes every Newton denominator
a valuation unit.

The local cutoff is sharp for this valuation statement: the d=q+1 roots
0 and t^(-1) times all q-th roots of unity have power sums zero through
q-1 and q*t^(-q) at q. This is a split root model, NOT an irreducible Keller
field or a sharpness example among actual Keller maps.

## 3. Finite source-generator and generic-linear criteria: HC4MYKO1

If the first d-1 power traces of EVERY source coordinate are in B, then all
coordinates belong to S by Section 2. Hence R=S is finite over B. The
Keller map is finite etale, so retained complex affine simple connectedness
makes its degree one and gives its polynomial inverse. The same conclusion
holds for any finite algebra-generating set of R, with each generator's
individual cutoff m-1. Conversely an inverse has d=1 and all these trace
conditions hold. These are exact criteria, not automatic trace regularity.

Here is why one GENERIC linear form is enough. List the finitely many
omitted boundary prime divisors D_i in Z. At each D_i some x_j has a pole.
If all coordinates were in O_(Z,D_i), then R would map into that DVR.
The resulting morphism Spec O_(Z,D_i) -> X -> Z agrees with the canonical
one at the generic point; separatedness of Z makes the two agree everywhere,
forcing the omitted generic point of D_i into X, a contradiction.

For each D_i let a_i=min_j v_(D_i)(x_j)<0. The constant coefficients c for
which ell_c=sum c_j x_j has v_(D_i)(ell_c)>a_i form a proper C-linear
subspace: reduce the terms of minimal valuation modulo the maximal ideal
after multiplying by a uniformizer to the power -a_i. Avoid these finitely
many subspaces. The resulting nonempty Zariski open U in C^n consists of
linear forms with a pole on EVERY omitted D_i. If ell_c for c in U has its
first d-1 power traces in B, Section 2 removes every D_i. With no omitted
prime divisor, normal Hartogs gives R=S; since both are affine the open
immersion is an isomorphism. Finite etaleness then gives the inverse.
If there are no D_i, the same Hartogs argument already gives the inverse.

Equivalently the set

    V = {c in C^n : Tr(ell_c^k) in B for 1<=k<=d-1}
      = {c : ell_c in S}

is a C-linear subspace. It is proper under nonautomorphy; any spanning set
of good directions, or any Zariski-dense set of them, therefore forces an
inverse. The boundary argument supplies a specific generic open on which
every direction detects all omitted divisors. It does NOT supply the missing
regularity of any such direction. No uniform degree or primitive-element
assumption is needed.

More precisely, a target height-one prime P occurs as a pole of one of
Tr(x_j^k), 1<=k<=d-1, if and only if at least one omitted boundary divisor
lies over P. If such a divisor exists, some coordinate has a pole there
and the local refinement detects it. Conversely a trace pole means that
some conjugate of that coordinate is not integral at a prime of S over P;
such a prime cannot lie in X. This describes the divisorial support exactly;
it makes no unsupported claim about higher-codimension boundary components.

## 4. Differential retractions are exactly the missing global trace: HC4MYKO2

Let delta_i be the commuting polynomial Nambu derivations dual to K:
delta_i(K_j)=delta_ij. They act on R because det(JK) is a scalar unit and
span its derivations over R. Let D_B be the target Weyl algebra, acting
on R by delta_i and on B by target partial derivatives.

The following are equivalent:

1. K has a polynomial inverse.
2. Tr_(L/F)(R) is contained in B.
3. There is a D_B-linear retraction p:R -> B of B -> R.

For 2 -> 1, apply Section 2 to the powers of every coordinate. For 1 -> 3,
R=B and the identity works. To prove 3 -> 2 and uniqueness, localize a
retraction over F. The localization R tensor_B F equals L: it is a finitely
generated algebraic domain over F and hence a field. Nondegeneracy of the
separable field trace pairing gives a unique a in L with

    p(r)=Tr_(L/F)(a*r).

Field trace commutes with each uniquely extended delta_i. Commuting p with
delta_i consequently gives Tr((delta_i a)*r)=0 for every r in L, whence
delta_i a=0. Since the delta_i span the ordinary source partial derivations,
a is in C. Finally p(1)=1 gives a=1/d. Thus every such retraction is exactly
the NORMALIZED FIELD TRACE, and its existence says that this field map is
regular on the WHOLE source R. It is not enough to have it on S, where it
is always regular, or only after a target localization.

Commutation of algebraic field trace with differential operators is already
retained in PF2TRACE1 and QZM1; the argument above proves the precise global
retraction assertion directly. No semisimplicity, D-module splitting theorem
or automatic consequence of de Rham acyclicity is imported.

## 5. Actual higher-dimensional Keller countercontrol

Credit Alpoge's announced map and Gallagher's inspected
[weighted-lift derivation](https://github.com/algal/jacobianfun/blob/main/RESEARCH.md),
Sections 1-2. We rederive the needed identities here. Put u=1+xy and

    A=u^3*z+y^2*u*(4+3*x*y),
    B0=y+3*x*u^2*z+3*x*y^2*(4+3*x*y),
    C0=2*x-3*x^2*y-x^3*z.

The exact Jacobian is -2. The three points (0,0,-1/4),
(1,-3/2,13/2),(-1,3/2,13/2) all map to (-1/4,0,0).
This is a genuine whole-source Keller control in THREE dimensions.
It is not a plane map or a four-variable unit-Hessian gradient.

Put gamma=2-3xy-x^2z, w=u*gamma, P=B0*C0, Q=A*C0^2. Direct expansion gives

    P=2*gamma+4*w-3*w^2,
    Q=w^2+w*gamma-w^3,
    w^3-2*w^2+P*w-2*Q=0.

This cubic is irreducible over C(P,Q,C0): its quotient over C[P,Q,C0]
is C[P,w,C0], a domain, and its ideal meets C[P,Q,C0] only at zero.
Its leading coefficient as a polynomial in w is one, so Gauss's lemma
gives field irreducibility. Conversely a generic root w reconstructs
gamma=(P-4*w+3*w^2)/2, x=C0/gamma,
y=(w/gamma-1)/x and z=(5-3*w/gamma-gamma)/x^2. Thus L=F(w) and d=3.
The exact companion-matrix checks below verify the field computations.

The source polynomial w is integral, with Tr(w)=2 and Tr(w^2)=4-2P;
an arbitrary integral source polynomial therefore need not give an inverse.
For the different source polynomial u=w/gamma one obtains

    c*u^3+b*u-4Q=0,
    c=P^3-P^2-18*P*Q+27*Q^2+16*Q, b=P^2-12*Q,
    Tr(u)=0, Tr(u^2)=-2*b/c.

The latter trace has a genuine pole: on P=0,Q=-16/27, c=0 but b=64/9.
This proves that trace cancellation at the first power is possible for an
ACTUAL Keller source, while the second power detects nonintegrality. It
also refutes automatic whole-source trace regularity, and hence automatic
D_B splitting, from the Keller condition in arbitrary dimension. It does
NOT refute a possible additional plane-specific or gradient-specific
argument. No external theorem about all higher dimensions is needed for this
exact control, and no novelty priority is claimed for the known map.

## 6. Relation to the finite trace-lattice programme and the remaining gate

The inspected [hostile-matrix primary note](https://github.com/alok/jacobian-two/blob/main/docs/trace-hostile-matrices.md)
already distinguishes a symmetric trace Gram matrix from an actual unital,
associative regular-representation algebra. Its Section 4 adds the cubic
trace tensor and divisibility needed to recover multiplication. Those
constraints concern the FINITE algebra S. Its field trace S -> B is regular
regardless of whether the larger open ring R has trace poles. Even all
multiplication identities for S do not assert Tr(R) subset B. The two
particular hostile curves were subsequently excluded by topology in that
source; this does not exclude every degree-six packet or every finite cover.
We import no unproved external trace-lattice exclusion.

The new criterion applies DIRECTLY to grad(h) in HC4 with n=4, independently
of Tang or a native source alignment. It also applies to its plane-Keller
branch through the retained HC4/JC2 equivalence. The unresolved step is to
force regular traces of actual source generators, one boundary-detecting
linear form, or a global differential retraction using genuine extra
plane/gradient data. The canonical nontorsion condition from KN remains
compatible with trace poles. Neither a Frobenius pairing for S nor a free
target-module pushforward supplies this missing regularity.

[The checker](verify_whole_source_traces.py) verifies exact Newton controls,
the local sharpness model, the actual three-dimensional map and trace
formulas, and differential trace pairing identities. Universal statements
are proved above, not inferred from its bounded examples. The receipt binds
immutable prior artifacts without replaying their checkers. Independent
replay, formal verification and external review are all absent. Full
HC4/JC2 and the user's original proof-or-counterexample objective remain open.
