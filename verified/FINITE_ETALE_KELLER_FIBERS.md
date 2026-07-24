# Finite étale algebras as Keller fibers

## 1. Definition

Let \(K\) be a field of characteristic zero. A nonzero finite étale
\(K\)-algebra \(A\) is a **Keller fiber** if there are a polynomial Keller map

\[
F:\mathbb A^m_K\longrightarrow\mathbb A^m_K
\]

and a point \(y\in\mathbb A^m(K)\) such that

\[
F^{-1}(y)\simeq\operatorname{Spec}A,\qquad
\dim_KA=\operatorname{gdeg}(F).
\]

The equality means that the fiber is full: every generic inverse sheet is
present, with none lost at the finite-normalization boundary.

## 2. Rank-classification theorem

The possible ranks of nonzero Keller fibers over a characteristic-zero field
are exactly

\[
oxed{1,3,4,5,\ldots}.
\]

More precisely, let \(P\in K[T]\) be squarefree of degree \(N\ge3\). Choose
\(a\in K\) with

\[
P'(a)P'''(a)
e0
\]

and put

\[
G(S)=P(a+S)-P(a).
\]

The root-engineered quadratic gauge attached to \(G\) is an explicit map

\[
F_G:\mathbb A^3_K\longrightarrow\mathbb A^3_K,\qquad
\det DF_G=-2,\qquad \operatorname{gdeg}(F_G)=N,
\]

with full fiber

\[
oxed{
F_G^{-1}\left(1,0,-rac{2P(a)}{P'(a)}ight)\simeq
\operatorname{Spec}K[T]/(P).
}
\]

Every finite étale algebra over an infinite field is monogenic, so this
realizes every finite étale \(K\)-algebra of rank at least three.

Rank one is realized by the identity. Rank two is impossible by the exact
Galois case below.

## 3. Exact generality of the degree-two exclusion

Let \(F:\mathbb A^m_K	o\mathbb A^m_K\) be Keller with geometric degree two.
Then

\[
K(x_1,\ldots,x_m)/K(F_1,\ldots,F_m)
\]

is a separable quadratic extension and hence Galois. Bass--Connell--Wright,
Chapter I, Theorem 2.1, states the Galois implication over an arbitrary ground
field: in characteristic zero, a polynomial map with invertible Jacobian and
Galois function-field extension is a polynomial automorphism. Thus no
algebraic-closedness assumption is required.

Reference:

- H. Bass, E. Connell, D. Wright,
  [*The Jacobian conjecture: reduction of degree and formal expansion of the inverse*](https://doi.org/10.1090/S0273-0979-1982-15032-7),
  Bull. Amer. Math. Soc. (N.S.) 7 (1982), 287--330.

Campbell proved the complex case; Razar and Wright gave later algebraic
accounts.

## 4. Scheme-theoretic reconstruction

For a quadratic-gauge target \((\pi,b,c)\) with \(\pi
e0\), write

\[
E(S)=G_\pi(S)-rac{g_1}{2}(bS^2+c),\qquad
R=K[S]/(E),\qquad s=Smod E.
\]

If \(E\) is squarefree, Bézout gives \(U,V\in K[S]\) with

\[
UE+VE'=1.
\]

Hence \(E'(s)\) is a unit in \(R\); its inverse is explicitly \(V(s)\). Put

\[
d=rac{E'(s)}{g_1},\qquad
eta(\pi,s)=rac{G_\pi'(s)/g_1-1-\pi s^2}{s},\qquad
Q=b-eta(\pi,s),
\]

and reconstruct

\[
 t=d^{-1},\quad x=sd^{-1},\quad y=Q-\pi s,\quad q=\pi d,
\]

\[
 z=d^2\left(q-rac{g_1}{g_3}y^2(1+3t)ight).
\]

These are elements of \(R\), not rational functions on geometric points. They
satisfy the source equations and map to \((\pi,b,c)\), defining one ring map
from the fiber algebra to \(R\).

Conversely, on the entire fiber scheme,

\[
tq=\pi\in K^*.
\]

Therefore \(t\) and \(q\) are units globally, and

\[
S=x/t,\qquad Q=y+xq
\]

are global fiber-ring elements satisfying \(E(S)=0\). The two constructions
are inverse on coordinate rings. Thus

\[
F_G^{-1}(\pi,b,c)\simeq\operatorname{Spec}K[S]/(E)
\]

scheme-theoretically, not merely as a set of geometric points.

## 5. Explicit optimal Hasse fiber and scaling ledger

Let

\[
P_5(T)=(T^3-19)(T^2+T+1)
\]

and

\[
G(T)=P_5(T)-P_5(0)=T^5+T^4+T^3-19T^2-19T.
\]

The normalized quadratic gauge has coefficient vector

\[
(g_1,g_2,g_3,g_4,g_5)=(-19,-19,1,1,1)
\]

and determinant \(-2\). If its target coordinates are
\((\Pi,B_0,C_0)\), the denominator-free displayed map is exactly

\[
(\Pi,19B_0,19C_0).
\]

Therefore its determinant is

\[
19^2(-2)=-722,
\]

and normalized target \((1,0,-2)\) becomes displayed target
\((1,0,-38)\). At the normalized target,

\[
G(T)-rac{-19}{2}(-2)=G(T)-19=P_5(T).
\]

Hence the complete fiber is

\[
\operatorname{Spec}\mathbb Q[T]/((T^3-19)(T^2+T+1)).
\]

It has points over every completion of \(\mathbb Q\) but no rational point.
A Chebotarev argument excludes such finite étale schemes in total degree at
most four, so the minimum Hasse-failing Keller-fiber rank is exactly five.

## 6. Consequences and chain update

The theorem immediately preserves the whole input algebra, including:

- connectedness and field decomposition;
- signatures at real places;
- splitting fields and the Galois action on geometric points;
- all local algebras \(A\otimes_KK_v\);
- good-prime factorization types;
- rational points, local points, and intersectivity.

In particular every finite separable field extension of degree at least three
is a connected Keller fiber.

This changes the earlier arithmetic chain:

1. **Absolute occurrence is settled.** Every finite étale algebra of rank at
   least three occurs directly, without tangent normalization, Hilbert
   irreducibility, or auxiliary factors.
2. The weighted tangent-admissibility theorem remains useful for occurrence
   inside the specific linear-tilt weighted family.
3. The quadratic-stabilized intersective transfer theorem remains valid as a
   weighted-presentation result, but its degree-two overhead is no longer
   needed for Keller-fiber existence.
4. The open arithmetic questions now concern realization inside one fixed map,
   bounded coordinate complexity, stable-equivalence multiplicity, and
   additional geometric restrictions on the ambient map—not abstract
   occurrence.

## 7. Verification and novelty audit

Run

```bash
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

The checker verifies:

- determinant \(-2\) for translated examples in degrees three, four, and five;
- explicit Bézout inversion and quotient-ring reconstruction;
- both compositions of the reconstruction identities;
- the exact target scaling of the Berend--Bilu map;
- determinant \(-722\), target \(-38\), and inverse polynomial \(P_5\);
- the fixed-map infinite Hasse-family identity.

The dated literature search supporting the qualified novelty statement is
recorded in
[`papers/common-arithmetic-fibers/LITERATURE_AUDIT.md`](../papers/common-arithmetic-fibers/LITERATURE_AUDIT.md).
