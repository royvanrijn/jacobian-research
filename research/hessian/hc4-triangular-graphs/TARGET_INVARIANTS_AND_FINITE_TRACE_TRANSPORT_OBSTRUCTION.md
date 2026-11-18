# Target invariants and finite trace-transport obstruction

18 November 2026. Actual complex polynomial plane Keller maps. The
finite normalization, its source base change, and its canonical etale
completion are distinct throughout.

Let B=C[P,Q] subset R=C[x,y], F=Frac B, L=Frac R, with
Jac(P,Q)=1. Let S be the full finite normalization of B in L and
Lambda=Hom_B(S,B), identified with the ORDINARY trace-dual in L.
Fix ANY polynomial source coordinates (x,y), and put delta=partial_x.
No generic or source-finite assumption on this frame is made.

Three consequences of the actual finite target are proved:

1. If the map is not invertible, then B intersect C[y]=C.
2. Put T=R tensor_B S. If the map is not invertible, NO nonzero
   replica q(y)*delta admits a regular derivation lift to T.
3. For a noninverse, R tensor_B (Lambda/S) admits NO connection
   relative to any such replica, even as an abstract R-module.

The second statement concerns a finite cover of the SOURCE, not
preservation of the original ring S inside L. It explains an exact
obstruction to transporting the common translation through finite
trace pairings. No automatic codifferent membership or fold exclusion
is asserted.

## 1. A polynomial target invariant forces inversion

Suppose B intersect C[y] contains a nonconstant polynomial. Choose

    h(y)=H(P,Q)

of the least positive y-degree N. For any c in C, factor H-c in the
UFD C[P,Q]. Every nonconstant target factor has nonconstant pullback,
since P,Q are algebraically independent. That pullback divides
h(y)-c in C[x,y], hence belongs to C[y]: the x-degree of a product
is the sum of the x-degrees. Each such pullback has positive degree
at least N by minimality.

The product has degree N. Thus H-c has precisely one irreducible
factor, with exponent ONE. In particular H-c is reduced for EVERY c.
This argument includes the possibility of a repeated single factor;
its multiplicity would also contradict minimality.

The morphism X=Spec R -> Y=Spec B is etale. Base change to the reduced
curve H=c is etale, so its source is reduced. Equivalently,

    R/(h(y)-c) is reduced for every c.                    (1)

A reduced C-algebra remains reduced under an etale algebra extension;
this assertion does not require the curve H=c to be smooth. Since the
whole source is C[x,y], (1) says that h(y)-c is squarefree in C[y].
If N>1, h' has a root a in C, and y=a is a repeated root of
h(y)-h(a), contradiction. Hence N=1 and y belongs to B.

Write y=H_0(P,Q). The Jacobian identity gives

    P_x=H_0,Q(P,Q),       Q_x=-H_0,P(P,Q).                (2)

Thus delta preserves B and its field F. Apply the
[direct replica index theorem](DIRECT_REPLICA_INDEX_FROM_A_DESCENDED_DIFFERENTIAL.md)
to the GIVEN replica delta=1*partial_x and intermediate field E=F.
It yields [L:F] dividing deg(1)+1=1. The map is birational, and the
[full-normalization inverse criterion](CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md)
then makes it invertible. This endpoint uses no source-finite frame,
coordinate-descent classification, or Rentschler theorem.

Consequently every noninverse actual map satisfies

    B intersect C[y]=C                                  (3)

in EVERY polynomial source frame.

## 2. The actual finite source-basechanged algebra

The algebra

    T=R tensor_B S

is finite locally free over R. Indeed a normal surface S is Cohen--Macaulay.
At every closed target point, a regular base parameter pair is a system
of parameters upstairs, hence is a regular sequence on S. Its B-depth
is two; the regular-base projective-dimension formula gives B-projectivity,
as in the retained normal-surface duality proof. It is also normal,
possibly a PRODUCT of normal domains, because S->T is the etale base
change of B->R.
It need not be a domain or the original normalization S.

Flat base change and finite projectivity identify its ordinary trace-dual:

    T^#=R tensor_B Lambda,
    T^#/T=R tensor_B (Lambda/S).                         (4)

Here the trace is the trace of the finite generic R-algebra
L tensor_F L. It is not a new use of Tr_(L/F) on just one factor.

There is a diagonal R factor of T. Indeed the diagonal source X in
X times_Y Z is a closed section of a finite morphism to X; it is
also open because X->Z is open and X->Y is etale. Hence

    T=R times T_other.                                  (5)

The first factor has zero trace defect. Nonzero trace defect in (4)
can survive on the other factors. Thus the diagonal source sheet alone
does not annihilate the original trace defect.

Let Delta in B generate the discriminant ideal of the finite projective
B-algebra S. Such a generator exists because the determinant line over
the factorial polynomial ring B is free. It is nonzero by separability.
The discriminant ideal of T/R is its base change:

    disc(T/R)=Delta(P,Q)*R.                              (6)

In (6), Delta(P,Q) means the image in R of the target polynomial.
For a noninverse map Delta is NONCONSTANT: a unit discriminant would
make S/B finite etale, and the retained inverse criterion applies.
The target embedding is injective, so its image is also nonconstant.

## 3. Exact trace transport does not assert regularity

There is an exact identity for the common translation, even without a
regular lift to T. Put a=P_x,b=Q_x and let D_P,D_Q be the uniquely
extended target coordinate derivations on L. Then

    delta=a*D_P+b*D_Q,
    D_P(a)+D_Q(b)=partial_x Jac(P,Q)=0.

Trace commutes with the target derivations, so for every f in L,

    Tr_(L/F)(delta f)
      =D_P Tr_(L/F)(a*f)+D_Q Tr_(L/F)(b*f).              (7)

The traces on the right are rational functions in F. Their regularity
is NOT supplied when a*f or b*f is outside Lambda. Also delta itself
need not preserve F, so (7) must not be replaced by the incorrect
identity Tr(delta f)=delta Tr(f). The divergence identity is compatible
with poles and supplies no endomorphism of Lambda/S by itself.

## 4. A regular trace connection preserves the discriminant ideal

Let D_R be any derivation of R, and suppose it has a regular derivation
lift D_T:T->T. The lift on the finite generic algebra is unique because
that algebra is a product of finite separable extensions of L.
Differentiation commutes with its field traces. This follows either
from the characteristic polynomial or from the unique extensions of
derivations to a common normal closure and the sum over embeddings.

On a local free chart choose an R-basis s_1,...,s_d and write

    D_T(s_i)=sum_j M_ij*s_j,    G_ij=Tr(s_i*s_j).

The coefficients of M are regular on that chart. Trace differentiation
gives the exact matrix identity

    D_R(G)=M*G+G*M^t,
    D_R(det G)=2*Tr(M)*det G.                            (8)

The determinant formula holds first generically, where G is invertible,
and then in the domain. A change of basis multiplies det G by the square
of a local unit, so the assertion that its principal ideal is D_R-stable
is independent of the chart. Thus globally

    D_R(Delta(P,Q)) belongs to Delta(P,Q)*R.             (9)

This is a necessary consequence of a regular lift. It is not a claim
that an arbitrary derivation preserving a discriminant ideal lifts.

## 5. No common replica lifts in the noninverse case

Set D_R=q(y)*partial_x with q in C[y], q!=0. This is locally nilpotent
on R. A nonzero principal D_R-stable ideal in R has an invariant
generator. In the present frame one can see this without any theorem:
if D_R(f)=a*f is nonzero, the left side has x-degree deg_x(f)-1,
whereas the right side has x-degree deg_x(a)+deg_x(f), impossible.
Hence D_R(f)=0, and f belongs to C[y].

If a regular lift to T existed, (9) would therefore imply

    Delta(P,Q) in C[y] intersect B.                    (10)

For a noninverse map (3) makes this intersection C, contradicting the
nonconstant discriminant in Section2. This proves the second assertion.
The lift itself was NOT assumed locally nilpotent; local nilpotence
was needed only on the base ring R.

In particular, when y is the proved common source coordinate of the
canonical etale completion A=B[Lambda], the complete translation on A
has NO regular lift to the finite source-basechanged normalization T.
Its induced generic derivation must have a pole on at least one finite
ramification divisor: otherwise the height-one intersection criterion
for each normal factor would make the derivatives of finite algebra
generators regular everywhere. On the etale locus a base derivation
has a unique regular lift, so such poles are necessarily ramified.
One cannot infer a regular action on (4), or a
trace-defect commutator, from completeness on A.

This does NOT contradict
[the invariant replica-extension criterion](PLINTH_IDEALS_AND_INVARIANT_REPLICA_EXTENSION.md).
That criterion concerns q(y)*partial_x preserving the ORIGINAL S inside L,
and is equivalent to y in S. A lift to T instead extends a derivation of
the LEFT source factor and determines its action on the right algebraic
factors through the target relations. The two requests are different.

## 6. Even an abstract trace-defect connection is obstructed

Put M_0=Lambda/S and M=R tensor_B M_0. The exact
[reduced-branch annihilator theorem](REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md)
gives Ann_B(M_0)=H*B, where H is the reduced equation of the ENTIRE
branch divisor. This is the full reduced support used in
[the logarithmic defect calculation, Sections1--2](LOGARITHMIC_TRACE_DEFECT_AND_NON_DENSE_DERIVATIVE_DOMAINS.md),
not merely an equality away from finitely many points.

Flatness of R over B gives the EXACT annihilator identity

    Ann_R(M)=H(P,Q)*R.                                  (11)

To see this without a base-change assumption on Hom, choose finite
B-generators m_1,...,m_N of M_0. The annihilator is the kernel of
B->M_0^N sending b to (b*m_i)_i. Tensor this kernel sequence with
the flat B-algebra R. The images of the m_i still generate M, so
the resulting kernel is exactly Ann_R(M). This proves (11).

Let D_R=q(y)*partial_x, q!=0. A D_R-CONNECTION on M means an
additive C-linear operator nabla satisfying

    nabla(r*m)=D_R(r)*m+r*nabla(m).

It need not be locally nilpotent. If r annihilates M, applying this
identity to every m shows that D_R(r) also annihilates M. Thus the
principal ideal (11) must be D_R-stable. The degree argument of
Section5 now forces H(P,Q) into B intersect C[y]. For a noninverse,
H is nonconstant and (3) forbids this. Hence

    M admits NO D_R-connection for a noninverse map.      (12)

This is stronger than the obstruction to regular algebra lifts.
Indeed a regular lift to T preserves T^# by trace differentiation:
for lambda in T^# and s in T,

    Tr(D_T(lambda)*s)=D_R Tr(lambda*s)-Tr(lambda*D_T(s)) in R.

It therefore induces precisely such a connection on T^#/T=M.
The exact logarithmic target operator on Lambda/S causes no conflict:
its base derivation is H*D_target, not the source replica D_R.
Equation(12) identifies the missing transport structure rather than
asserting that the trace defect vanishes.

## 7. Exact controls and proof boundary

The checker controls the least-degree/reducedness distinction, the
critical-root obstruction, cyclic trace Gram matrices and their connection
matrices, arbitrary triangular changes of basis, discriminant derivatives,
and positive/negative cyclic trace-defect connection relations,
and positive/negative regular-lift examples. The finite cover
u^e=h(x,y) has

    D(u)=D(h)*u/(e*h),
    Tr(connection matrix)=(e-1)*D(h)/(2*h).

For D=partial_x, h=y gives an invariant-branch regular lift, whereas
h=x+y^2 gives a polar lift. The nonlocally-nilpotent derivation
x*partial_x lifts to u^e=x, illustrating why local nilpotence on the
base is essential to the conclusion of Section5.

Two controls guard the source hypotheses. For (P,Q)=(x,y^m), the target
coordinate Q pulls back with a repeated critical fibre; the Jacobian
vanishes. For the Laurent Keller map P=u^e,Q=v/u^(e-1), the critical
root u=0 is missing, so the whole-plane step (1) cannot be used to
force u^e to have degree one. Neither control is an actual plane
Keller counterexample. The branched cyclic covers are finite-algebra
controls, not normalizations of hypothetical noninverse plane Keller maps.

No automatic ordinary-codifferent membership, density of a derivative
domain, preservation of the original S by the common translation, or
exclusion of either surviving fold family is proved. The remaining
landing problem still requires an additional actual-source implication.
The evidence is an internal written proof with exact controls, not
formal verification or external human review.
