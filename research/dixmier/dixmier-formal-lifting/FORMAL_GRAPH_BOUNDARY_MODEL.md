# Formal graph freeness can hide failure of finite generation

14 November 2026. **QGM1.** There is an explicit cyclic coherent algebraic
Weyl bimodule, flat over its parameter, whose formal completion is free
of left rank one and whose special fiber is the diagonal classical graph,
but whose generic fiber is not finitely generated from the left. The
generic fiber is nevertheless holonomic as a bimodule.

Thus formal left freeness, algebraic coherence, holonomicity and smooth
special-fiber graph support do not force a chosen algebraization to
remain finite from the left. The example locates the missing boundary
condition in a graph-module approach to [QFL2](BOUNDED_REALIZATION.md).
It does not exclude the existence of a better algebraization: this
example has one after a formal gauge.

## The exact rational canonical transformation

Let K have characteristic zero and put
\[
 R=K[t],\qquad A=R\langle x,p\rangle/([p,x]-t),\qquad
 q=1-t^2x.
\]
The powers of q form an Ore denominator set. Explicitly, q is a position
polynomial and [p,q]=-t^3, so moving differential operators past a
denominator introduces only further powers of q. Write
\[
 S=A[q^{-1}]=K[t,x,q^{-1}][p;t\partial_x].
\]
Define a homomorphism phi:A->S by
\[
             X=\phi(x)=x/q,\qquad
             P=\phi(p)=q^2p-t^3q.                   \tag{1}
\]
The scalar correction is the symmetric-order correction for q^2p:
(t/2) partial_x(q^2)=-t^3q. Directly,
\[
                    [P,X]=t q^2\partial_x(x/q)=t.
\]
Thus the assignment preserves the Weyl relation exactly.

Let M=S with its natural left A action and right A action
m*a=m phi(a). These actions commute, so M is an A-bimodule.

## Cyclic and coherent as a bimodule

Every element of S lies in A q^(-k) for some k>=0. Equation (1) gives
\[
                         \phi(1+t^2x)=q^{-1}.
\]
Consequently the right orbit of 1 contains q^(-k) for every k, and
left multiplication then gives all of S. Thus M is cyclic over
A^e=A tensor_R A^op.

The algebra A^e is a noetherian iterated Ore extension, with the right
pair having the opposite commutator sign and the two pairs commuting.
Its cyclic module M is therefore finitely presented. Also M is flat
over R: its position coefficient ring is a localization of R[x], and
its PBW expansion is free over that coefficient ring.

## Not finite from the left, even generically

The ascending union
\[
                    M=\bigcup_{k\ge0} A q^{-k}      \tag{2}
\]
is strictly ascending. If q^(-(k+1)) belonged to A q^(-k), multiplying
on the right by q^k would give q^(-1) in A. This is impossible: q is
not a unit in the finite PBW algebra, as its positive spatial degree
shows.

Every finite set of proposed left generators is contained in one term
of (2), so it cannot generate M. This argument still holds after
extending scalars to K(t): q is a nonconstant position polynomial
there. Hence M_K(t) is not finitely generated as a left A_K(t)-module.
Localizing only the parameter ring cannot remove this generic defect.

## Free of left rank one after formal completion

For every N, q is a unit modulo t^N, with inverse the finite geometric
sum in t^2x. Therefore
\[
                  S/t^NS \cong A/t^NA,\qquad
                         \widehat M\cong\widehat A
                         \quad\text{as left modules}.       \tag{3}
\]
The right action is the formal canonical transformation (1). Modulo t,
it is the identity, so
\[
 M/tM\cong K[x,p]
\]
as the diagonal graph bimodule over the two copies of K[x,p]. Its
support is the smooth diagonal Lagrangian in the product with opposite
symplectic form.

The divisor q=0 is absent from every formal t neighborhood, because
q=1 modulo t. At any nonzero parameter it is the affine divisor
x=t^(-2), which moves to infinity as t approaches zero. That is why
the formal completion does not see the failure in (2).

## The generic bimodule is holonomic

Set L=K(t), so A_L is a Weyl algebra and A_L^e is a rank-two Weyl
algebra. Filter the cyclic bimodule M_L by words of length at most m
in the four left/right generators applied to 1.

Write an element in rational PBW form sum_b c_b(x)p^b. If a term
arises from a word of length <=m, let k be a bound on its q pole
order and d=deg_infinity(c_b). The following two bounds hold:
\[
                              k+b\le m,\qquad d+b\le3m.
 \tag{4}
\]
They follow directly from the four generator actions:

- Left x raises d by one.
- Left p either raises b by one or differentiates c_b, raising k by
  at most one and lowering its degree at infinity.
- Right x multiplies by x/q. Its r-th derivative has q pole order
  at most r+1, while the corresponding Ore term lowers b by r.
- Right p multiplies by q^2p-t^3q. Ore differentiation of q^2 or q
  introduces no new q pole, and raises d+b by at most three.

After choosing the common denominator q^m, every such term belongs to
the span of
\[
                    x^a q^{-m}p^b,\qquad
                           0\le a\le4m,\quad0\le b\le m.
\]
The filtered dimension is at most (4m+1)(m+1). Thus GKdim(M_L)<=2.
Since M_L is nonzero and finitely generated over A_L^e isomorphic
to A_2(L), Bernstein's inequality gives GKdim(M_L)>=2.
Equality proves holonomicity.

For the role of finite generation from a specified side in graph
bimodules, compare
[Lauritzen--Thomsen, The graph of a Weyl algebra endomorphism](https://arxiv.org/html/1904.03128).
The example here satisfies bimodule coherence and generic holonomicity
while deliberately failing the required left finite generation.

## This is a bad algebraization, not an intrinsic impossibility example

The formal graph in (3) is isomorphic to the regular formal bimodule.
Indeed let
\[
                 G=x^2p+t x,\qquad U=\exp(tG)\in\widehat A^\times.
\]
Then phi=Ad(U) on the completion. For x, the iterated commutators
sum to x/(1-t^2x); for p the commutator series terminates at the
second term and gives q^2p-t^3q. Multiplication m->mU identifies the
completed graph bimodule with the regular one.

The regular algebraic bimodule A is therefore another algebraization
of the same formal isomorphism class, and it is free of left rank one.
The intertwining unit U is formal, with unbounded spatial support.

This distinction is essential. QGM1 disproves the inference that an
arbitrary coherent algebraization with formal left freeness will be
left finite generically. It does not disprove the existence of some
left-finite algebraization after a suitable formal gauge.

## The same failure in every Weyl rank

Adjoin any number of unchanged spectator Weyl pairs and let phi act as
the identity on them. Equivalently, use A_n over K[t] and the same
q and first-pair formulas (1). The strict union (2), cyclic generation
and formal left freeness are unchanged. The generic growth bound is
O(m^2) for the first pair times O(m^(2n-2)) for the spectators.
Bernstein's inequality therefore gives the holonomic dimension 2n
over A_n^e, while generic left finite generation still fails.
In particular the phenomenon occurs at Weyl rank two.

## Consequence for the research direction

A module approach must control a left lattice or its behavior at
infinity, rather than only producing an algebraic coherent bimodule
with the correct formal graph. The good algebraization must be chosen
with the necessary left finiteness and rank; formal triviality alone
does not certify those properties for a proposed choice.

For the Dixmier construction problem, the missing statement remains
existence of such an effective choice for a noninvertible classical
map. [QFL2](BOUNDED_REALIZATION.md) shows that a common spatial bound
would suffice. This example explains why coherence and holonomicity
alone cannot supply that bound by inspecting an arbitrary algebraic
model.

[Exact symbolic controls](check_curve_and_module.py) check (1), the
generator identity, and the two nonzero commutators in the formal
conjugation calculation. The strict infinite chain, completion,
noetherian coherence and generic holonomicity use the written proof.
No formal verification, external review, independent theorem replay,
novelty priority or Dixmier counterexample is claimed.
