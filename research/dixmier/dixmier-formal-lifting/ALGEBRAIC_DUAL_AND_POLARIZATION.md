# Algebraic dual sections and the missing polarization

14 November 2026. **QEF1.** For a torsion-free rank-one left module M
over a left noetherian Ore domain A,
\[
       M\text{ is finitely generated over }A
              \quad\Longleftrightarrow\quad
                   \operatorname{Hom}_A(M,A)\ne0.   \tag{1}
\]
Thus a nonzero algebraic left-dual section is a precise test for the
left-finiteness part of a proposed rank-one graph realization.

This criterion does not descend automatically from the formal graph.
For the QGM1 boundary model the algebraic left dual is zero, while
the completed left dual is the whole completed Weyl algebra. Bounded
polynomial quantum assignments also show that classical coordinate
agreement does not force actual local nilpotence of a graph constraint,
and that raw differential orders can jump arbitrarily under finite
quantum gauges. These are explicit limits on extending QSO1 to every
bounded full lift.

## The rank-one dual criterion

Let Q(A) be the Ore division ring. Assume M is nonzero and left
torsion-free, and
\[
                    \dim_{Q(A)}(Q(A)\otimes_A M)=1.
\]
If 0!=f in Hom_A(M,A), Ore localization makes f a nonzero map between
one-dimensional left Q(A)-vector spaces. It is injective after
localization. Its original kernel is torsion and therefore zero.
Hence f embeds M into A. Left noetherianity of A makes M finitely
generated.

Conversely suppose M is finitely generated. Embed M into Q(A) by
choosing a basis after localization. Its finitely many generators
d_i have a common nonzero right denominator c: d_i c is in A.
Right multiplication by c gives a nonzero left A-linear map M->A.
This proves (1).

No holonomicity is needed for this lemma. For Weyl bimodules the
algebraic left dual in (1) is the dual used in
[Lauritzen--Thomsen, Theorem 2.5](https://arxiv.org/html/1904.03128).
That theorem proves right finite generation of the dual of a holonomic
bimodule. It does not assert that the dual is nonzero. The model below
shows why nonvanishing is substantive.

The rank-one and torsion-free hypotheses must be established
algebraically for the model in question. They cannot simply be inferred
from a formal completion. Also (1) supplies left finiteness, not a free
left basis; the latter remains necessary to read a right action as an
endomorphism of the fixed Weyl algebra rather than of another
endomorphism ring.

## The algebraic dual of the boundary model vanishes

Use [QGM1](FORMAL_GRAPH_BOUNDARY_MODEL.md):
\[
 A=K[t]\langle x,p\rangle/([p,x]-t),\quad
 q=1-t^2x,\quad M=A[q^{-1}],
\]
with right action x->x/q, p->q^2p-t^3q. Its left module is the
strict union of free left modules A q^(-k).

If f:M->A is left linear, then
\[
              f(1)=q^k f(q^{-k})\in q^k A
                         \quad\text{for every }k.
\]
But
\[
                          \bigcap_{k\ge0}q^kA=0,
\]
because a nonzero element divisible on the left by q^k has spatial
degree at least k. Thus f(1)=0, and torsion-freeness gives
f(q^(-k))=0 for every k. Therefore
\[
                        \operatorname{Hom}_A(M,A)=0. \tag{2}
\]
The same proof works generically over K(t).

In contrast, q is a unit modulo every t^N, and QGM1 gives
M/t^NM isomorphic to A/t^NA as a left module. Consequently
\[
 \operatorname{Hom}_{A/t^N A}(M/t^NM,A/t^NA)
                   \cong (A/t^NA)^{op}
\]
as endomorphism rings, or simply A/t^NA as right modules.
Their compatible inverse limit is nonzero. In particular
\[
 \operatorname{Hom}_{\widehat A}(\widehat M,\widehat A)
                        \cong\widehat A
\]
as a right module, whereas the completion of (2) is zero.
Left duality and formal completion do not commute for this model.

This is failure of effectivity for the chosen algebraization. The
regular bimodule A is another algebraization of the same formal
isomorphism class and has nonzero dual. There is no intrinsic
nonalgebraizability conclusion.

## Even generic simplicity does not remove this failure

The generic boundary bimodule from QGM1 is simple over its enveloping
Weyl algebra, in addition to being holonomic. Here is a direct proof.

Work over K(t) and write its right images as X=x/q and
P=q^2p-t^3q. The right image algebra contains q^(-1)=1+t^2X and also
\[
                          p=q^{-2}P+t^3q^{-1}.
\]
Let N be a nonzero sub-bimodule. Choose a nonzero element n of N of
minimal differential order b in p. The operation
\[
                         n\longmapsto q\,nX-x\,n
                                   =q[n,X]
\]
preserves N. If b>0 and n has leading coefficient a(x), this operation
has exact order b-1 and leading coefficient b*t*a(x)/q, since
partial_x X=q^(-2). This contradicts minimality. Hence N contains
a nonzero position function f in K(t)[x,q^(-1)].

Clear its q denominator using left multiplication to get a nonzero
polynomial r(x) in N. Since p belongs to both the left and right
acting algebras, [p,r]=t r' is again in N. Repeated differentiation
gives a nonzero scalar, so 1 belongs to N. Cyclic generation from QGM1
then gives N=M. This proves generic simplicity.

Unchanged spectator pairs preserve this conclusion in higher rank:
their left/right commutators differentiate the finite spectator PBW
polynomial until a nonzero first-pair coefficient remains, after which
the preceding argument applies. Thus even a simple holonomic
rank-one generic bimodule can have zero algebraic left dual and a
formally free completion.

## A concrete derived boundary class

The preceding failure also has an exact Ext description. Identify
A q^(-k) with A. Its inclusion into A q^(-(k+1)) is right
multiplication by q. The countable direct-limit resolution by free
left modules gives
\[
 \operatorname{Ext}^1_A(M,A)
   \cong
 \operatorname{coker}\left[
   \prod_{k\ge0} A\longrightarrow\prod_{k\ge0} A,\quad
       (y_k)\longmapsto(y_k-qy_{k+1})
                         \right].                  \tag{3}
\]
The sequence (1,1,...) represents a nonzero class.

Indeed a preimage would satisfy
\[
 y_0=1+q+\cdots+q^{N-1}+q^Ny_N
                         \quad\text{for every }N.
\]
Multiplying by 1-q on the left gives
\[
             (1-q)y_0-1\in q^N A
                         \quad\text{for every }N.
\]
The intersection above is zero, so (1-q)y_0=1. This is impossible:
1-q=t^2x is not right invertible in A. This proves nonvanishing
without computing finite support boxes.

In the completed algebra q is a unit, and the recurrence
y_k-qy_(k+1)=1 has the solution obtained recursively from y_0=0.
Thus this inverse-system obstruction disappears after completion.
The calculation gives a concrete boundary class invisible in every
formal parameter jet.

## A bounded canonical lift need not give local nilpotence

In the same parameter Weyl algebra, consider the exact polynomial
automorphism
\[
                  \phi(x)=x+t^2p,\qquad \phi(p)=p.   \tag{4}
\]
It has spatial degree one, a polynomial inverse, and identity classical
reduction. Let A^phi be its graph bimodule and put
\[
                      G=R_x-L_x.
\]
On normal PBW symbols m(x,p),
\[
 G(m)=m\,(x+t^2p)-x\,m
     =t\,\partial_p m+t^2p\,m.                      \tag{5}
\]
The products in the middle are Weyl products; the expression on the
right is in normal PBW symbols.

Equation (5) implies G(A) is contained in tA, so G is nilpotent modulo
every t^N and topologically nilpotent on the completion. Nevertheless
the highest p term of G^n(1) is
\[
                            t^{2n}p^n
\]
for every n. Thus G is not locally nilpotent on the algebraic generic
fiber.

A fuller expression is
\[
 G^n(1)=\sum_{j=0}^{\lfloor n/2\rfloor}
       \frac{n!}{2^j j!(n-2j)!}
                         t^{2n-j}p^{n-2j}.
\]
This exact degree-one example shows that even boundedness together
with coordinate agreement at t=0 does not imply the algebraic local
nilpotence needed for a finite Taylor or Kashiwara-type reduction.

In
[Kanel-Belov--Fedorov, Lemmas 4.4--4.6](https://arxiv.org/html/2610.33919v1),
actual local nilpotence of the graph constraint follows from the
prescribed arithmetic center relation and the support argument. It
does not follow from classical h=0 agreement alone. The example above
explains the distinction.

The inverse automorphism in (4) restores the identity graph. Hence
the example refutes a statement about arbitrary representatives,
not a possible normalization theorem asserting existence of a good
representative.

The [birational endomorphism theorem of Lauritzen--Thomsen](https://arxiv.org/html/1508.02499)
also requires an actual endomorphism of the finite Weyl algebra.
The boundary model already has a birational division-ring change:
x=X/(1+t^2X) and p=(1+t^2X)^2 P+t^3(1+t^2X).
But X=x/q is not an element of the finite Weyl algebra. That theorem
does not supply the missing polynomial frame or left-finite model.

## Raw differential orders are not fixed by the classical map

If [P,Q]=t, then for every N>=1
\[
                     Q\longmapsto Q+t^2P^N         \tag{6}
\]
is an exact finite quantum shear, with inverse obtained by subtraction.
It preserves the classical t=0 map and a common finite spatial bound
for each fixed N. If ord(P)=m>0 and Nm>ord(Q), its new generic
differential order is Nm.

For example, with [p,x]=t,
\[
 P=p^2+x,\qquad Q=-p,\qquad Q_N=-p+t^2P^N
\]
has [P,Q_N]=t and order profile (2,2N) for N>=1. Its classical
reduction remains (p^2+x,-p), of orders (2,1). Adding unchanged
spectator pairs embeds both examples into rank two and higher.

Accordingly, QSO1 cannot be extended to every bounded lift by assuming
that the displayed classical orders persist. A useful all-degree
normalization theorem would need to guarantee the existence of a
finite algebraic polarization and a representative with the desired
profile. Bounds or order growth attached to one raw representative
are insufficient.

## Applicable research and the remaining positive problem

[Ben-Zvi--Nevins, D-bundles and integrable hierarchies](https://ems.press/content/serial-article-files/31805)
study framed D-bundles on curves. Their framing controls the
associated graded module at infinity, and their finite-colength
invariant records the local second Chern class. This suggests a
concrete direction: seek an algebraic lattice with an effective framing
at infinity, so that a nonzero algebraic dual section in (1) can be
constructed and bounded.

This is a proposed application. Their theory already works with
D-bundles and specified framing data; it does not automatically turn
our formal four-dimensional graph into such an object or bound its
charge. A new comparison and an actual bound would be required.

The precise positive questions are now:

1. Can a nonzero formal dual section be made algebraic for some
   appropriate rank-one model, with controlled behavior at infinity?
2. Can that model be given a free left frame compatible with the right
   Weyl action?
3. If a scalar polarization is used, can actual local nilpotence and a
   useful minimal order profile be achieved by finite algebraic changes?

None is answered merely by formal existence, holonomicity, or bounded
coefficients of an arbitrary representative. A positive construction
for a noninvertible classical map would meet the
[QFL2 realization criterion](BOUNDED_REALIZATION.md); the common
spatial bound and those realization conditions remain unknown for
Long's map.

[Exact controls](check_effectivity_controls.py) verify the bounded
canonical shear, its graph-constraint action, the Hermite formula
through five finite controls, and the order-jumping examples.
The rank-one criterion, infinite dual intersection, Ext class and
all-N conclusions use the written proofs. No bounded lift of Long,
Dixmier counterexample, formal verification, external review,
independent theorem replay or novelty priority is claimed.
