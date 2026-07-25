# Finite étale algebras as Keller fibers

## 1. Definition

Let \(K\) be a field of characteristic zero. A nonzero finite étale
\(K\)-algebra \(A\) is a **Keller fiber** if there are a polynomial Keller map

\[
F:\mathbb A^m_K\longrightarrow\mathbb A^m_K
\]

and a point \(y\in\mathbb A^m(K)\) such that

\[
F^{-1}(y)\simeq\operatorname{Spec}A,
\qquad
\dim_KA=\operatorname{gdeg}(F).
\]

The equality says that the fiber is **full**: every generic inverse sheet is
present. A shorter special fiber can occur when sheets escape through the
nonproperness boundary; fullness excludes precisely that loss.

## 2. Effective Jacobian-one realization

Let \(P\in K[T]\) be squarefree of degree \(N\ge3\). Choose \(a\in K\) with

\[
P'(a)P'''(a)\ne0
\]

and put

\[
G(S)=P(a+S)-P(a)=g_1S+\cdots+g_NS^N.
\]

Then \(g_1=P'(a)\), \(g_3=P'''(a)/6\), and \(g_N\ne0\), so the
root-engineered quadratic gauge applies. Its normalized map

\[
F_G=(\Pi,B,C):\mathbb A^3_K\longrightarrow\mathbb A^3_K
\]

has determinant \(-2\). Compose with the target automorphism

\[
L(\Pi,B,C)=\left(\Pi,-\frac{B}{2},C\right)
\]

and write \(\widetilde F_G=L\circ F_G\). Then

\[
\boxed{\det D\widetilde F_G=1},
\qquad
\operatorname{gdeg}(\widetilde F_G)=N.
\]

The distinguished target has second coordinate zero, so it is fixed by this
normalization:

\[
y_{P,a}=\left(1,0,-\frac{2P(a)}{P'(a)}\right).
\]

At this target the full fiber is

\[
\boxed{
\widetilde F_G^{-1}(y_{P,a})\simeq
\operatorname{Spec}K[T]/(P).
}
\]

The construction is effective. With total degree,

\[
\deg\Pi\le7,
\qquad
\deg B\le6N+2,
\qquad
\deg C\le6N,
\]

so

\[
\boxed{\max_i\deg(\widetilde F_G)_i\le6N+2.}
\]

The estimates use \(\deg t=2\), \(\deg q\le5\), and

\[
\deg\bigl(t^2x^{k-2}q^k\bigr)\le6k+2,
\qquad
\deg\bigl((xq)^k\bigr)\le6k.
\]

## 3. Complete rank classification

Every finite étale algebra over an infinite field is monogenic. After
extension to a separable closure, the algebra is a product of \(N\) copies of
the field; the primitive elements form the nonempty open set on which the
Vandermonde discriminant

\[
\prod_{i<j}(X_i-X_j)^2
\]

is nonzero. Thus every finite étale \(K\)-algebra of rank \(N\) is
\(K[T]/(P)\) for a squarefree degree-\(N\) polynomial.

Consequently every rank at least three occurs by the construction above.
Rank one is realized by the identity. Rank two is impossible: a degree-two
Keller map has a separable quadratic, hence Galois, function-field extension;
the Campbell--Razar--Wright Galois case makes it an automorphism. For an
arbitrary characteristic-zero ground field, descend to a finitely generated
subfield, embed it in \(\mathbb C\), apply the complex theorem, and descend the
unique formal inverse.

Hence the possible ranks of nonzero Keller fibers are exactly

\[
\boxed{1,3,4,5,\ldots}.
\]

## 4. Scheme-theoretic reconstruction

For a quadratic-gauge target \((\pi,b,c)\) with \(\pi\ne0\), set

\[
E(S)=G_\pi(S)-\frac{g_1}{2}(bS^2+c),
\qquad
R=K[S]/(E),
\qquad
s=S\bmod E.
\]

If \(E\) is squarefree, Bézout gives \(U,V\in K[S]\) with

\[
UE+VE'=1.
\]

Therefore \(E'(s)\) is a unit in \(R\), with explicit inverse \(V(s)\). Put

\[
d=\frac{E'(s)}{g_1},
\qquad
Q=b-\beta(\pi,s),
\qquad
\beta(\pi,S)=\frac{G_\pi'(S)/g_1-1-\pi S^2}{S},
\]

and reconstruct in \(R\)

\[
t=d^{-1},
\qquad
x=sd^{-1},
\qquad
y=Q-\pi s,
\qquad
q=\pi d,
\]

\[
z=d^2\left(q-\frac{g_1}{g_3}y^2(1+3t)\right).
\]

These are elements of the quotient ring, not pointwise rational functions.
Conversely, on the entire fiber ring,

\[
tq=\pi\in K^*,
\]

so \(t\) and \(q\) are units globally and

\[
S=x/t,
\qquad
Q=y+xq
\]

are global fiber-ring elements satisfying \(E(S)=0\). The two maps are inverse
on coordinate rings. Thus

\[
F_G^{-1}(\pi,b,c)\simeq\operatorname{Spec}K[S]/(E)
\]

scheme-theoretically.

## 5. Scalar-extension compatibility

For every characteristic-zero field extension \(K\hookrightarrow K'\), the
construction commutes with coefficientwise scalar extension:

\[
\widetilde F_{G,K'}=\widetilde F_G\otimes_KK'.
\]

The distinguished fiber becomes

\[
(K[T]/(P))\otimes_KK'\simeq K'[T]/(P).
\]

Hence connectedness, field decomposition, real signatures, splitting fields,
Galois actions, and all local algebras are transported functorially.

In particular, every finite separable field extension \(L/K\) of degree at
least three occurs as a connected full fiber of a Jacobian-one map of
\(\mathbb A^3_K\), with coordinate degree at most \(6[L:K]+2\).

## 6. Explicit optimal Hasse fiber

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

and determinant \(-2\). If its target coordinates are \((\Pi,B_0,C_0)\), the
denominator-free displayed map is exactly

\[
(\Pi,19B_0,19C_0).
\]

Therefore its determinant is

\[
19^2(-2)=-722,
\]

and normalized target \((1,0,-2)\) becomes displayed target \((1,0,-38)\).
At the normalized target,

\[
G(T)-\frac{-19}{2}(-2)=G(T)-19=P_5(T).
\]

Hence

\[
F^{-1}(1,0,-38)\simeq
\operatorname{Spec}\mathbb Q[T]/((T^3-19)(T^2+T+1)).
\]

This fiber has points over every completion of \(\mathbb Q\) but no rational
point. A Chebotarev argument excludes such finite étale schemes in total
degree at most four, so the minimum Hasse-failing Keller-fiber rank is exactly
five.

## 7. Consequences for the earlier chain

1. **Absolute occurrence is settled.** Every finite étale algebra of rank at
   least three occurs directly, with determinant one and an effective degree
   bound.
2. Weighted tangent-admissibility remains useful for occurrence inside the
   specific weighted linear-pencil family, but it is not an absolute
   existence condition.
3. Quadratic-stabilized intersective transfer remains valid as a
   weighted-presentation theorem, but its degree-two overhead is unnecessary
   for general Keller-fiber existence.
4. The remaining arithmetic questions concern one fixed map, minimal
   coordinate complexity, stable-equivalence multiplicity, and additional
   restrictions on the ambient map.

The geometric-degree spectrum theorem now also has a second direct existence
proof: apply this realization to any squarefree polynomial of degree
\(N\ge3\).

## 8. Verification and formalization

Run

```bash
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

The exact checker verifies translated examples in degrees three, four, and
five, determinant-one normalization, the \(6N+2\) degree bound, Bézout
inversion, both quotient-ring reconstruction compositions, the explicit
quintic scaling, determinant \(-722\), target \(-38\), and the fixed-map
infinite Hasse-family identity.

A staged Lean project is stored in
[`formal/finite-etale-keller`](../formal/finite-etale-keller). Stage one
formalizes the explicit quintic Jacobian, output normalizations,
inverse-polynomial identity, and a constructive Bézout inverse in the quotient
algebra. The universal marked-line identities and full scheme reconstruction
are the next formal modules. The historical degree-two theorem remains a
separate external input until formalized.

The dated qualified novelty search is recorded in
[`papers/common-arithmetic-fibers/LITERATURE_AUDIT.md`](../papers/common-arithmetic-fibers/LITERATURE_AUDIT.md).
