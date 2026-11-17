# Nonlinear base-chart parameters: the exact omitted-point endpoint

14 November 2026. The native five-chart theorem does not transfer a
Hessian determinant through nonlinear source changes. This note handles
an entire remaining class instead: every base-only graph x=h(y) or
y=h(x), with arbitrary polynomial automorphisms of its three retained
dual variables over the free-base polynomial ring.

If h(0)=0, no such reparameterization has nonzero constant Hessian.
If h(0) is nonzero, existence of any unit-Hessian reparameterization
is exactly existence of a complex plane Keller pair generating the
unit ideal. One fixed physical graph x=1 already realizes every such
pair, with genuine polynomial retained-dual and embedding inverses.
The existence of that pair is left unresolved.

This is an omitted-point condition, not the
[singleton-zero endpoint](P_GRAPH_SINGLETON_ENDPOINT.md): the two
coefficient ideals are (1) and (s,U), respectively.

## 1. A composite trace on a three-dimensional hyperplane forbids unit Hessian

**HC4MYGU0.** Let f(s,z1,z2,z3) be polynomial in characteristic zero.
If

\[
 f(0,z)=H(a(z)),\qquad \deg H\geq2,              \tag{1}
\]

then det Hess f is not a nonzero constant. No coordinate or embedding
hypothesis on a is required, and no Tang theorem is used.

Over the algebraic closure, if a is nonconstant choose a root xi of
H' and a point with a=xi. Such a point exists by restricting a to a
line where it is nonconstant. At that point the three-by-three dual
Hessian block is

\[
 H''(\xi)\,\nabla a\nabla a^T,                  \tag{2}
\]

of rank at most one. Adding the one normal row and column raises rank
by at most two, so the full four-by-four Hessian is singular. If a is
constant, the dual block is zero everywhere on s=0 and the same
conclusion is immediate. This proves the claim. Algebraically, for
nonconstant a the boundary determinant is divisible by H'(a): a
three-by-three block H'(a)S+H''(a)vv^T has its adjugate divisible by
H'(a), and its determinant divisible by H'(a)^2.

Apply this to either graph eliminating a base variable when its
omitted base coordinate vanishes identically on the free-base-zero
slice. For example x=R(y,p,q,r) with R(0,p,q,r)=0, or y=R(x,p,q,r)
with R(0,p,q,r)=0, has boundary potential Lp0^2+Mp0. Arbitrary
relative polynomial retained-dual automorphisms preserve that identity.
In fact arbitrary polynomial retained-dual substitutions still satisfy
(1). Thus the whole zero boundary-trace class has no unit-Hessian
member, in every degree. Merely having one point over the old base
origin is not the asserted hypothesis; the whole parameter slice must
map to that base point.

## 2. Every base-only graph avoiding the origin is a relative coordinate

Let a and b be the three-entry coefficient rows of the original
linear forms A and B in (p,q,r). Put u=1+xy. Their ordered two-by-two
minors are exactly

\[
 m_1=yu^3,\qquad m_2=xu(u+1),\qquad
 m_3=2x^2(2u+1).                               \tag{3}
\]

Their ideal is (x,y). The following literal identities certify this:

\[
 K=(2u+1)m_3/2-4xm_2=x^2,
\]
\[
 x=m_2/2-(3y+xy^2)K/2,
\]
\[
 y=m_1-(3y^2+3xy^3+x^2y^4)x.                   \tag{4}
\]

Reverse containment is evident. Now substitute a polynomial base
curve (x,y)=(X(s),Y(s)) with (X,Y)=(1) in K[s]. Equations (3)-(4)
give a certificate sum_i c_i*m_i=1 over K[s]. The third row
w=(c3,-c2,c1) completes a,b to a three-by-three matrix of determinant
one, since their cross product is (m3,-m2,m1).

Thus (A,B,W) is a relative linear automorphism over K[s]. Replacing
B by H(A)+NB, with N nonzero, is a polynomial shear and rescaling.
Consequently Psi is a coordinate of K[s][p,q,r], with polynomial
relative companions and inverse. No division by X, Y or u is used.

For the actual graph curves (s,h(s)) or (h(s),s), the ideal condition
is exactly h(0) nonzero. Write h(s)=c+s*k(s), c nonzero; the Bezout
identity is h/c-s*k/c=1. These curves are genuine closed graphs,
with parameter s recovered by the remaining base projection.

**HC4MYGU1.** On every such origin-avoiding base-only graph, arbitrary
relative retained-dual automorphisms give precisely all polynomial
coordinates of K[s][z1,z2,z3]. Indeed a fixed coordinate Psi can be
carried to any other coordinate, with chosen polynomial companions,
by composing their actual relative automorphisms. If h(0)=0,
HC4MYGU0 instead forbids unit Hessian after every such automorphism.

## 3. Classify unit-Hessian relative coordinates

Work over C. A relative coordinate f in C[s][z1,z2,z3] has no zero
of its three dual partial derivatives at any source point: its
relative inverse specializes at every s, so its coordinate Jacobian
is invertible. In particular grad f misses the zero target and cannot
be a polynomial automorphism.

Suppose now det Hess f is a nonzero constant. The reproduced
[Tang alternative HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
supplies a jointly affine constant two-plane W. Project to the single
free base s. Only ranks zero and one occur. All adaptations below are
constant affine changes preserving the free-base fibration, so both
the unit Hessian and relative-coordinate property are preserved.

### Rank one cannot occur

Put the normal form into

\[
 f=D(U,Z)+sA(U,Z)+VC(U,Z),\qquad
 \operatorname{Jac}(A,C)\in C^\times.            \tag{5}
\]

C is nonconstant, so choose (U0,Z0) with C=0. The two equations

\[
 \nabla D+s\nabla A+V\nabla C=0                 \tag{6}
\]

have a unique solution (s,V), because Jac(A,C) is a nonzero constant.
At the resulting point all three dual derivatives f_U,f_Z,f_V vanish.
This contradicts the relative-coordinate property. No restriction to
a generic base value is substituted for the property at every s.

### Rank zero is exactly a unit-ideal plane Keller pair

Constant linear dual coordinates give

\[
 f=D(s,U)+P F(s,U)+Q G(s,U),\qquad
 \det\operatorname{Hess}f=\operatorname{Jac}_{s,U}(F,G)^2. \tag{7}
\]

If F=G=0 at some (s0,U0), their U derivatives cannot both vanish by
the unit Jacobian. Choose P,Q so that D_U+P F_U+Q G_U=0. Then again
all three dual derivatives vanish. Hence (F,G) has no common zero,
equivalently

\[
 (F,G)=(1)\quad\hbox{in }C[s,U].                \tag{8}
\]

Conversely take any F,G satisfying (8), and choose a,b with aF+bG=1.
For arbitrary polynomial D(s,U), the literal relative coordinate map

\[
 (U,P,Q)\longmapsto
 (U,\ D+PF+QG,\ -bP+aQ)                         \tag{9}
\]

has a polynomial inverse. If V=D+PF+QG and W=-bP+aQ, then

\[
 P=a(V-D)-GW,\qquad Q=b(V-D)+FW.                \tag{10}
\]

The two-by-two matrix determinant is aF+bG=1; no inverse of the
plane map (F,G) is assumed. When (F,G) is Keller, (7) is unit Hessian
and grad f omits zero because its P,Q outputs are F,G.

This proves **HC4MYGU2**: existence of a unit-Hessian relative
coordinate in four variables is exactly existence of a unit-ideal
plane Keller pair. Every such relative coordinate has the cotangent
form (7) after a constant affine adaptation. The forward implication
uses Tang; the converse and explicit inverses are elementary. This
is not an existence proof for the plane pair.

## 4. One fixed graph gives the whole endpoint by genuine embeddings

For completeness the fixed graph x=1 has especially short formulas.
Write s=y, u=1+s, and define

\[
 A=u^3p+3u^2q-r,\quad
 T=u(u+1)p+2(2u+1)q,\quad W=-(2u+1)p/2-4q.
\]

The frame (A,T,W) has constant determinant -1, and B=T-(2-3s)A.
For any relative coordinate f(s,U,P,Q) with companions U,W from
(9), set

\[
 T=(f-H(U))/N+(2-3s)U,
\]
\[
 p=-4T-2(2u+1)W,\quad
 q=(2u+1)T/2+u(u+1)W,\quad
 r=u^3p+3u^2q-U.                               \tag{11}
\]

These are polynomial and give A=U, Psi=f. Invert the physical frame,
recover f=Psi and U=A, and then use (10) to recover P,Q. Thus (11)
is a polynomial automorphism of the three retained dual coordinates
over C[s], with an explicit inverse. Its graph map into C^5 is a
genuine closed embedding onto the single hyperplane x=1. This is not
a merely Keller parameterization or a construction away from a divisor.

**HC4MYGU3.** For every base-only graph x=h(y) or y=h(x) with h(0)
nonzero, existence of a unit-Hessian pullback after relative retained-dual
automorphisms is equivalent to existence of a unit-ideal plane Keller
pair. One fixed graph x=1 suffices to realize the entire endpoint.
For h(0)=0 unit Hessian is impossible by HC4MYGU0.

By translating a missing target and the Nullstellensatz, the absence
of unit-ideal plane Keller pairs is exactly the assertion that every
complex plane Keller map is surjective. This is the precise plane
assertion left unresolved here. We do not identify it with all JC2;
in particular no argument here rules out a surjective nonautomorphic
plane Keller map. It is also distinct from the singleton-zero
properness problem, which has a nonempty affine zero fibre of length one.

The fixed x=1 graph cannot contain both original v2 marked base points
whose x coordinates are opposite. No automatic transfer of that fixed
collision is asserted. Any unit-Hessian member constructed from a
hypothetical unit-ideal Keller pair would already have nonautomorphic
gradient by its omitted zero target.

## Controls, reproduction and remaining scope

The row F=s, G=sU-1 is unimodular, with a=U,b=-1, but its Jacobian
is s. It gives a genuine relative-coordinate embedded graph through
(11), with actual Hessian determinant s^2, not a counterexample.
Conversely F=s,G=U has unit Jacobian but a common zero. Its cotangent
potential is not a relative coordinate over C[s]. Both controls are
retained to distinguish the coordinate and Keller hypotheses.

The [checker](verify_relative_coordinate.py) verifies the bordered
trace factor, all minors and ideal certificates, frame and graph
inverses, both plane normal determinants, critical-point reconstruction,
and the two opposite hypothesis controls. The all-rank necessity,
Nullstellensatz steps and endpoint quantifiers are written proofs.
[Actual receipt](receipt-relative-coordinate.json).

```sh
timeout 60s .venv/bin/python research/hessian/hc4-triangular-graphs/verify_relative_coordinate.py
```

General dual-dependent base graphs after nonlinear retained-dual
reparameterization, nonlinear changes mixing the free base with dual
parameters, and uncovered embeddings remain. No unit-ideal Keller
pair, surjectivity theorem, singleton-zero properness theorem, formal
verification, independent external review or unrestricted HC4/JC2
solution is asserted.
