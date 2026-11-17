# Quantum radical towers cannot realize the Long division extension

**QRG1.** A bounded full quantum lift of Long cannot have its ambient
quantum division ring generated from its image division ring by a finite
internal tower of radical adjunctions. This excludes the construction
orbit of finite monomial cotangent coverings, including their finite
compositions and arbitrary actual quantum division-ring changes of
coordinates. Parameter poles in those coordinate changes do not avoid
the obstruction.

The proof uses the non-Galois cubic classical extension and elementary
noncommutative Hensel lifting. It applies to the actual generic division
extension, without requiring a coordinate change to preserve the
\(h\)-valuation. It is not a classification of all finite Weyl pairs.
No positive bounded Long lift or all-representative nonlifting is proved.

## The valuation setting and the meaning of an internal radical tower

Let \(E\subset Q\) be division rings with a central discretely valued
coefficient field, common central uniformizer \(\pi\), and the same
value group \(\mathbb Z\). Their completions embed as
\(\widehat E\subset\widehat Q\). Assume the residue division rings are
commutative characteristic-zero fields \(P\subset H\).

An internal radical tower is a finite chain of division subrings
\[
 E=E_0\subset E_1\subset\cdots\subset E_m=Q,\qquad
 E_i=\langle E_{i-1},\theta_i\rangle_{\rm div},\quad
                   \theta_i^{k_i}\in E_{i-1},\quad k_i\ge2. \tag{1}
\]
Every \(\theta_i\) belongs to \(Q\). There is no assumption that it
commutes with \(E_{i-1}\). In particular (1) denotes a generating
condition; it does not assert that adjoining a noncentral root always
has degree \(k_i\), or that noncommutative extension degrees multiply.

Suppose the residue extension has no new radicals:
\[
            b\in H,\quad b^k\in P,\ k\ge2
                              \quad\Longrightarrow\quad b\in P.
                                                               \tag{2}
\]
Then every division ring obtained from \(E\) by a tower (1) inside \(Q\)
lies in \(\widehat E\). If \(H\ne P\), it therefore cannot be all of \(Q\).

## A root with prescribed nonzero residue is unique

Let \(A\) be a valuation-zero element of \(\widehat E\), and let
\(b\in P^\times\) satisfy \(b^k=\bar A\). There is exactly one
valuation-zero root \(X\in\widehat Q\) of \(X^k=A\) whose residue is \(b\);
that root belongs to \(\widehat E\).

Choose a lift \(X_1\in\widehat E\) of \(b\).
If \(A-X_n^k\) has valuation at least \(n\), choose
\(c_n\in\widehat E\), of nonnegative valuation, with residue
\[
       \bar c_n=
       \frac{\overline{\pi^{-n}(A-X_n^k)}}{k b^{k-1}}\in P,
 \qquad X_{n+1}=X_n+\pi^n c_n .
\]
Because the residue field is commutative,
\[
 (X_n+\pi^n c_n)^k-X_n^k
                  \equiv \pi^n k b^{k-1}\bar c_n
                                      \pmod{\pi^{n+1}}.
\]
The error improves by one. Completeness gives the required root
in \(\widehat E\).

For uniqueness, the exact identity
\[
 X^k-Y^k=\sum_{j=0}^{k-1}X^j(X-Y)Y^{k-1-j}
\]
holds without commutation. If two roots have the same residue \(b\)
and their difference first occurs at valuation \(n\ge1\), dividing
by \(\pi^n\) gives the nonzero residue
\(k b^{k-1}\overline{\pi^{-n}(X-Y)}\), a contradiction.
No centrality of \(A\), \(X\), or \(Y\) was used.

For (1), proceed inductively. If \(E_{i-1}\subset\widehat E\),
multiply \(\theta_i\) by a central integral power of \(\pi\) to make
its valuation zero; this is possible because the value group is
\(\mathbb Z\). Its \(k_i\)-th power remains in \(E_{i-1}\).
Condition (2) puts its nonzero residue in \(P\). The preceding
uniqueness identifies it with the root in \(\widehat E\).
Thus \(E_i\subset\widehat E\).

If the last \(E_i\) were \(Q\), taking closures would give
\(\widehat Q=\widehat E\), since the complete subring
\(\widehat E\) is closed. Their residue fields would agree.
This proves the obstruction.

## The Long residue extension satisfies the no-radical condition

Use [QDI1's rational Long chart](QUANTUM_DIVISION_AND_COMPLETION.md).
Its source function field has birational coordinates \(r,v,a,d\);
the image field has coordinates \(r,s,t,D\), with the fourth image
independent of \(a\) after the displayed chart correction. Eliminating
\(v\) gives
\[
                  r a^3-2a^2+s a+4t=0,\qquad
                   v=(s-4a+3ra^2)/2.                     \tag{3}
\]
The extension \(H/P\) of full classical function fields has degree
three. For example over \(K(r,s,D)\), equation (3) is irreducible:
it is linear in \(t\) with unit coefficient and its quotient
polynomial ring is \(K(r,s,D)[a]\); Gauss's lemma then gives
irreducibility over \(K(r,s,D)(t)\).

The discriminant in \(a\) is
\[
 \Delta=4s^2-4rs^3+128t-144rst-432r^2t^2 .
\]
As a quadratic polynomial in \(t\), it has discriminant
\[
                           256(4-3rs)^3\ne0.
\]
Thus \(\Delta\) is nonconstant and squarefree in \(t\) over \(K(r,s,D)\),
and is not a square in \(P\). The degree-three extension is not Galois;
its normal closure has group \(S_3\). These assertions survive any
extension of the characteristic-zero constant field.

We may first extend the constants to contain all roots of unity.
If \(b\in H\) has \(b^k\in P\), then \(P(b)/P\) is Galois:
every root of the binomial lies in \(P(b)\), even if that binomial
is reducible, and its minimal polynomial is separable.
A nontrivial subextension of a degree-three extension would be
the entire \(H\), contradicting its non-Galois nature.
Therefore (2) holds.

This is an assertion about radicals lying inside \(H\).
An \(S_3\) cubic is solvable by radicals after adjoining elements
outside that cubic field. No failure of solvability by radicals
is claimed here.

## Application to actual bounded full lifts

For a bounded full formal Long lift, its generic generator images lie
in the finite PBW rank-two Weyl algebra over \(K((h))\).
Let \(Q\) be its ambient Weyl division ring, and \(E\) the image of
the source Weyl division ring. The \(h\)-Gauss valuation has group
\(\mathbb Z\), central uniformizer \(h\), commutative source residue \(H\)
and image residue \(P\), as proved in QDI1. In particular
\[
                  [\widehat Q:\widehat E]_\ell
                      =[\widehat Q:\widehat E]_r=3.
\]
The preceding theorem rules out (1) for this \(Q/E\).

Over constants not containing all roots of unity, extend the constants
first. In the ambient Weyl division ring over the enlarged scalar field,
take the division closure of \(E\) and the new scalars.
The original Weyl generators are finite division expressions in \(E\)
and the \(	heta_i\), so the same generating presentation and radical
relations persist in that enlarged ambient Weyl division ring.
No identification with a tensor product of division rings is needed.
The geometric cubic discriminant remains nonsquare, giving the same
contradiction. A ramified scalar parameter may be
treated using its common central uniformizer instead of \(h\), provided
the finite lift is regular there and the same classical Long residue
is retained.

The argument concerns the completed residue degree. It does not
identify the possibly larger uncompleted division indices, and it does
not prove \(Q=E\) when the residue extension is trivial.

## The actual monomial coverings and their complete construction orbit

In the Laurent Weyl algebra \(K[v,v^{-1}][b;h\partial_v]\),
the finite pair
\[
               P=v^k,\qquad Q_0=\frac1k v^{1-k}b,\qquad k\ge2
\]
satisfies \([Q_0,P]=h\). Its ambient division ring is generated over
the pair's image division ring by \(v\):
\[
                         v^k=P,\qquad b=k v^{k-1}Q_0.      \tag{4}
\]
Thus it is an internal radical adjunction. Identity spectator Weyl
pairs do not change that fact.

All coordinate automorphisms here are over the generic scalar field
\(K((h))\), and in particular fix \(h\).
Actual source division-ring automorphisms leave the image division
ring unchanged; target automorphisms transport (4) and its generators.
No valuation preservation or parameter regularity of those
automorphisms is needed. Finite compositions of these coverings and
actual division-ring automorphisms also give finite internal radical
towers: for two embeddings use the intermediate image division ring
and transport the inner tower along the outer embedding.

Consequently no member of this entire construction orbit can be a
bounded full Long lift. This includes actual coordinate changes with
arbitrary finite parameter poles and birational quantum changes
defined on the division rings. An arbitrary formal Hamiltonian gauge
need not preserve the actual division ring, and is outside this
transport statement.

## Controls that delimit the claim

A proper actual radical extension can have identical completions.
In the commutative control
\[
       E=K((h))(u),\quad Q=K((h))(x),\quad
       u=x+hx^2,\quad \theta=1+2hx,
\]
one has \(\theta^2=1+4hu\) and \(\bar\theta=1\).
Its binomial expansion belongs to \(\widehat E\), while
\(\theta\notin E\). Thus Hensel lifting puts roots into the completion,
not automatically into the actual image division ring.
This control has residue degree one and does not meet the strict
residue-extension hypothesis of the obstruction.

Commutativity of the residue is essential for the uniqueness proof.
Over matrices let \(X=\operatorname{diag}(1,-1)\) and
\(Y=X+h e_{12}\). Then \(X^2=Y^2=I\) and \(X,Y\) have the same
reduction but are distinct. Their residue algebra is noncommutative.
Accordingly no unrestricted matrix or projective-carrier statement
is being transferred from the scalar valued-division setting.

[Exact controls](check_quantum_radical_towers.py) check the cubic
discriminants, the finite monomial cover pair, noncommutative
first-correction and uniqueness identities in a Weyl algebra,
the binomial completion control, and the matrix residue boundary.
The unrestricted tower exclusion and all allowed coordinate changes
use the written induction and valuation proof.

No actual mixed Long pair, bounded lift, DC1/DC2 resolution,
formal verification, external review or classification of all finite
embeddings into the Laurent Weyl algebra is claimed.

