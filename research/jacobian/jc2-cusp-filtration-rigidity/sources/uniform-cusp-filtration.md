# Uniform rigidity for coprime cusp filtrations

Producer: swarmHQ ROOT coordinator (configured gpt-6-astra; hosted identity unverified).
Independent hostile reviewer: configured gpt-5.6-sol / ultra; hosted identity unverified.
Date: 2026-11-07 UTC. Basis: f8c3ebfc4fe785cba1414e8ec3e1b3b9ef923886.
Claim: COPRIME-CUSP-FILTERED-EULER-RIGIDITY-1.
Evidence tier: MANUAL. Lifecycle: PROVISIONAL.
Independent verdict: CONFIRM_UNIFORM_CUSP_THEOREM. Literature novelty unassessed.

This extends the [earlier exact (2,5) case](filtered-euler-lift-rigidity-swarmHQ-root-sol56-20261107T1755Z.md)
to every coprime binomial cusp. The general Poisson/divergence lemma replaces
small-piece enumeration. The proof below is self-contained polynomial calculus
and commutative algebra; the earlier case is a comparison, not a logical premise.
[Exact producer and independent review](../box/uniform-cusp-filtration-rigidity-swarmHQ-20261107T1825Z/README.md).

## Statement and scope

Let p,q>=2 be coprime integers, s=pq-p-q, B=C[x,y], and w=x^p-y^q.
Then s>=1. Every f has unique normal form sum_{i=0}^{p-1} x^i f_i(y,w).
Define delta by the maximal normal-monomial weight, with delta(x,y,w)=(q,p,s).
Suppose a polynomial Keller map F^*:x->P,y->Q, Jac(P,Q)=J in C*, satisfies

    delta(F^*f)=e delta(f) for EVERY nonzero f in B, integer e>=1.    (1)

Then e=1, P=alpha x, Q=beta y, alpha^p=beta^q, alpha beta!=0.
For J=1 also alpha beta=1 and alpha^(p+q)=1. Conversely these diagonal maps
satisfy (1) with e=1. No reduction of an arbitrary Keller map to (1) is asserted.

## General degree-zero Poisson rigidity lemma

Let delta be a nonnegative integer semidegree on C[x,y], zero exactly on
nonzero constants. Suppose A=gr_delta(B) is a domain, the ordinary Jacobian
bracket satisfies delta({f,g})<=delta(f)+delta(g), and its induced degree-zero
bracket on A has fraction-field Poisson center C.
Then divergence is an injective C-linear map

    {polynomial derivations L: delta(Lf)<=delta(f) for all f} -> C.  (2)

Proof. Since the induced bracket is nonzero, choose f,g with
 delta({f,g})=delta(f)+delta(g). For L=a partial_x+b partial_y, direct expansion gives

    L{f,g}={Lf,g}+{f,Lg}-(div L){f,g},  div L=a_x+b_y.

Each of the first three terms has degree at most delta(f)+delta(g), because
L preserves the filtration. Multiplicativity therefore gives delta(div L)<=0
when div L is nonzero. Thus div L is a constant (possibly zero).

If div L=0, polynomial integration over C gives H in B with H_x=b,H_y=-a,
so L={H,-}. This uses only a_x+b_y=0 and termwise antiderivatives, not a
rational primitive or a boundary extension. If L!=0, H is nonconstant and
delta(H)>0. Its leading form Hbar is not a scalar. Since the fraction-field
Poisson center is C, there is a leading form fbar in A with {Hbar,fbar}!=0:
otherwise Hbar commutes with A and, by quotient rule, its entire fraction field.
Consequently delta(Lf)=delta(H)+delta(f)>delta(f), a contradiction.
The kernel of divergence is zero, proving (2); this vector space has dimension
at most one. In particular, if it contains an actual polynomial E with
nonzero constant divergence, every member is a constant multiple of E.
This lemma does NOT assert that such an E exists for an arbitrary filtration.

## Verifying every cusp hypothesis

The monic relation x^p=y^q+w gives the stated normal basis. Its highest-weight
relation is X^p=Y^q; comparison of the same normal basis shows

    A=C[X,Y,W]/(X^p-Y^q), deg(X,Y,W)=(q,p,s).

Coprimality makes A a domain. Thus delta is a semidegree, with positive degree
on every nonconstant polynomial and finite-dimensional filtration pieces.
Its fraction field is C(z,W), with X=z^q and Y=z^p. Bezout for p,q shows that
z belongs to Frac(A), so this is the actual same fraction field.

The Jacobian brackets on B are

    {x,y}=1, {y,w}=-p x^(p-1), {x,w}=-q y^(q-1).

The latter two have degrees p+s and q+s exactly. The normal basis and
Leibniz rule therefore prove delta({f,g})<=delta(f)+delta(g) for ALL f,g.
The leading degree-zero bracket has

    {X,Y}=0, {Y,W}=-p X^(p-1), {X,W}=-q Y^(q-1).

On the fraction field this is {z,W}=-z^(s+1). A rational function commuting
with z has zero W-derivative, hence lies in C(z); if it also commutes with W
its z-derivative vanishes, hence it is in C. The Poisson center is exactly C.
No normality or everywhere symplectic property of A is asserted or required.

The actual derivation E=q x partial_x+p y partial_y preserves delta, since
E(w)=pq w and every normal monomial is an E-eigenvector. Its divergence is
p+q. The lemma now proves that every delta-preserving polynomial derivation
is a scalar multiple of E. This replaces all case-by-case generator enumeration.

## Transport through the entire polynomial Keller pair

Define polynomial D_P=(Q_y partial_x-Q_x partial_y)/J and
D_Q=(-P_y partial_x+P_x partial_y)/J. Set D=q P D_P+p Q D_Q.
Then D(F^*f)=F^*(Ef). Condition (1) induces an injective graded map
phi:A -> A multiplying degrees by e.

Every polynomial derivation L has finite filtered order bounded, and if L!=0
attained on a generator, by

    r=max(delta(Lx)-q, delta(Ly)-p, delta(Lw)-s).

Indeed the normal basis and Leibniz rule bound each differentiated monomial.
If r>0 for D, its leading part is a nonzero degree-r derivation of A.
But delta(D(F^*f))=delta(F^*(Ef))<=e delta(f), so it annihilates phi(A).
The field Frac(A) is finite algebraic over Frac(phi(A)): both are finitely
generated fields of transcendence degree two and the latter is a subfield.
In characteristic zero this extension is separable. Differentiating the minimal
polynomial of any element shows that a derivation zero on the base field is zero
on the extension. This contradicts the nonzero leading derivation.
Therefore D preserves delta and is mu E by (2).

Both D_P,D_Q have divergence zero. Hence div(D)=q D_P(P)+p D_Q(Q)=p+q,
which forces mu=1. It follows that E(P)=qP,E(Q)=pQ. Coprimality p,q>=2 implies
that the only monomial of E-weight q is x, and the only one of weight p is y.
Thus P=alpha x,Q=beta y; J!=0 implies alpha beta!=0. Applying (1) to y gives
e=1. Applying it to w then gives alpha^p=beta^q: otherwise its degree is pq>s.
This proves the theorem. No integrality of B over F^*B is used.

## Full homogeneous polynomial lift consumer

Let S=C[X,Y,W,T]/(X^p-Y^q-WT^(p+q)), weights(q,p,s,1), with
{X,Y}=T^(p+q), {Y,W}=-pX^(p-1), {X,W}=-qY^(q-1), T central.
This is the actual Rees ring via X=x t^q,Y=y t^p,W=w t^s,T=t.
Suppose U,V,H are homogeneous polynomial elements of weights(qe,pe,se),
U^p-V^q=T^((p+q)e)H, and their reductions give an injective Poisson map A->A.
Bracket the full relation with V:

    -p U^(p-1){U,V}=T^((p+q)e){V,H}.

T is prime and does not divide U, so T^((p+q)e) divides {U,V}. Homogeneity
makes the quotient a scalar c. Cancellation gives {V,H}=-pc U^(p-1).
Reduction modulo T, the prescribed Poisson map, and Ubar!=0 force c=1.
At T=1 the pair is Keller with Jacobian one. To see full scaling, the relation
defines a graded endomorphism Theta of S by X->U,Y->V,W->H,T->T^e. For f of
delta-degree d, its homogeneous Rees lift maps to a homogeneous element of
degree ed. Its reduction is the injective central image of fbar, so it is not
divisible by T. Its specialization is F^*f, of degree exactly ed. Thus (1) holds
for every f and the theorem gives e=1. This excludes every e>1
full polynomial lift in the entire coprime binomial-cusp family.
Rational/completed-local lifts and finite jets are not excluded by this statement.

## Controls, review and remaining gap

For ordinary positive weighted degree the leading degree-zero Jacobian bracket
is zero, so the center hypothesis fails. This proof cannot be applied there.
For arbitrary infinitely-near filtrations neither scalar graded Poisson center
nor a nonzero filtration-preserving Liouville derivation has been established.
The whole source being A2 does not establish either input. The theorem above
verifies every required cusp hypothesis directly in the actual source ring.

The independent hostile review confirmed the general lemma, polynomial
Hamiltonian integration, every coprime-cusp hypothesis, finite-order transport,
and the full Rees consumer. It explicitly checked the homogeneous-Rees-lift
argument giving full scaling, expanded above. ROOT read the entire sealed
204-LF review and checked its complete execution and reading record. A help
subtotal reported as29LF was30LF, still within its100LF cap; this does not
change the mathematical scope. Exact original report bytes are preserved.

Replay is desk-only. The evidence package supplies canonical seals and a
SHA256SUMS manifest; no CAS, formal proof or new external source audit is claimed.
The proof uses characteristic-zero polynomial integration and separability;
these standard algebraic foundations are not newly formalized.

JC2 remains unresolved. Applying this result to arbitrary hypothetical Keller
counterexamples still requires an actual-map reduction to suitable filtrations
or another argument handling the remaining valuations. No generic Liouville
existence, properness, surjectivity, outside-boundary descent, or general
finite-normalization conclusion follows. Further special-case extensions need
a demonstrated connection to that missing reduction.

## OPENS RAISED

None. Existing actual-map realization and finite-normalization gaps remain open.

## COLLISIONS

status: EMPTY

- NONE — the report contains no explicitly raised `OPEN[...]` entries.

Corpus basis: `f8c3ebfc4fe785cba1414e8ec3e1b3b9ef923886` (Git blobs only).

<!-- BODY-END -->

## Seal

- Body definition: every byte through the unique standalone `<!-- BODY-END -->` line,
  including its terminating newline; this seal is outside the body.
- Body bytes: `9567`.
- Body SHA-256:
  `2e8ff8d38a4452fa6bd4ec6f6405eda06fd2b80506a45a3d64e3cf97dcbba4ce`.
- Frozen basis: `f8c3ebfc4fe785cba1414e8ec3e1b3b9ef923886`.
