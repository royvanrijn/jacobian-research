# All-degree base graphs in one arbitrary dual direction

14 November 2026. Let `Psi=L*A^2+M*A+N*B`, with `L*N != 0`, be the
scaled Meng--Yang potential. Put

\[
 w=a(s)p+b(s)q+c(s)r,\qquad R\in K[s,w],\quad\deg_w R>0,
\]

where a,b,c are polynomials, not all zero. A base-dependent affine offset
in w may be absorbed into R. Coefficients need not form a unimodular row.
Allow every actual polynomial automorphism of the three retained duals
over K[s], with polynomial inverse and unrestricted word length.

**HC4MYHD1. Every y=R(x,w) graph has an invertible gradient under actual
unit Hessian, in every direction (a,b,c).** This includes genuinely old-r
dependent graphs without any boundary-trace restriction.

**HC4MYHD2. The same conclusion holds for x=R(y,w) whenever (a,b) is
not identically (0,0).** Thus all mixed directions except the pure old-r
direction are closed, in arbitrary graph and coefficient degree.

**HC4MYHD3. In the remaining pure-old-r opposite family, unit Hessian
forces its boundary trace to be a nonzero constant or a nonconstant
linear polynomial in r.** A nonautomorphic member would then supply,
respectively, the existing unit-ideal or singleton-zero plane Keller
endpoint. This is a necessary reduction, not a converse realization or
a properness proof for that subfamily.

The proof counts whole fixed-base value fibres using ranks of two linear
rows. A fixed cubic/quartic calculation controls every degree of R at
once. No coefficient search, branch expansion or boundary enumeration
is used. The retained Hodge-trace theorem, Tang and the plane-line inputs
remain explicit; no new classification is imported.

## 1. A general rank-change fibre formula (HC4MYHD0)

Fix a scalar s where the free base is nonzero and the vector
ell=(a(s),b(s),c(s)) is nonzero. For this section these are field
constants. At old omitted-base value z, write A=A_z*d and B=B_z*d,
where d=(p,q,r) and the two coefficient rows A_z,B_z are independent
for every z. Suppose their entries are polynomial in z. On the graph
z=R(s,w), with w=ell*d, restrict those rows to the two-dimensional
direction plane ker ell. Their rank is two except at the roots of

\[
             K_\ell(z)=\ell\mathbin{\cdot}(A_z\times B_z).   \tag{1}
\]

Assume K_ell is a nonconstant polynomial. A rank-zero restriction is
impossible: both independent ambient rows would then be multiples of
ell. At a rank-one restriction there are two possibilities:

- **Quadratic:** A_z is not proportional to ell. On the affine plane
  ell*d=w, use T=A as a coordinate. B is affine in T. A generic potential
  value has two T roots, each with a free affine-line coordinate.
- **Linear:** A_z is proportional to ell. It is constant on the affine
  plane, while B is nonconstant. Every value fibre is one affine line.

At a rank-two restriction, (A,B) are affine coordinates and
`L*A^2+M*A+N*B=lambda` is one affine line.

Let m be the number of distinct w values with K_ell(R(s,w))=0, and k
the number of those values of quadratic type. All are counted once,
regardless of multiplicity. Off them, solving the affine two-row system
gives an actual product `(A1 minus m points) times A1`. At each quadratic
value there are two affine lines; at each linear value one. Thus, for
generic lambda,

\[
 E(f_s^{-1}(\lambda))=(t-m)t+2kt+(m-k)t=t^2+kt,\quad t=uv. \tag{2}
\]

All products are explicit stratum isomorphisms. This is not multiplicativity
for finite covers. A generic lambda simultaneously avoids the finitely
many quadratic discriminant values. If K_ell has even one quadratic-type
root z0, nonconstant R(s,w) attains z0, so k>=1. The
[HC4MYHC0 single-trace theorem](HODGE_TRACE_CLOSURE.md) then forces
gradient invertibility under actual unit Hessian.

Completing ell to linear coordinates is done only on this fixed scalar
fibre. No global polynomial inverse for the row ell over K[s] is assumed.
Actual retained-dual automorphisms preserve the whole fixed-base fibre.

## 2. Every direction in the y graph (HC4MYHD1)

Here old x=s, old y=z, and u=1+sz. The actual coefficient rows are

\[
 A_z=(u^3,3su^2,-s^3),
\]
\[
 B_z=(z^2u(4+3sz),\ z+3sz^2(4+3sz),\ 2s-3s^2z).
\]

Their cross product is exactly

\[
 (2s^2(2sz+3),\ -su(u+1),\ zu^3).                         \tag{3}
\]

It never vanishes for s nonzero: the zero of its first component has
sz=-3/2 and its middle component is then nonzero. Thus the ambient
rows are always independent. Formula (1) becomes

\[
 K_\ell=2as^2(2sz+3)-bsu(u+1)+czu^3.                       \tag{4}
\]

If c=0, this polynomial has degree two when b is nonzero and degree
one otherwise. A_z cannot be proportional to ell, since its r entry
is -s^3. Every root is therefore of quadratic type.

If c is nonzero, K_ell has degree four. There is at most one value z0
where A_z is proportional to ell: after normalizing its nonzero third
entry, its other entries determine both u^2 and u^3, hence u uniquely.
If no such value exists, every root is quadratic. Otherwise put
v=sz and tau=sz0. After rescaling ell to A_z0 and dividing by s^2,
the polynomial is

\[
 P_\tau(v)=2(1+\tau)^3(2v+3)
       -3(1+\tau)^2(v+1)(v+2)-v(v+1)^3.                    \tag{5}
\]

It vanishes at v=tau and

\[
                         P_\tau'(\tau)=-6(1+\tau)^3.
\]

For tau!=-1, the single possible linear-type root is simple, so the
quartic has another root, necessarily quadratic. For tau=-1 the
polynomial is `-v*(v+1)^3`: v=0 is quadratic. This exhausts all
directions and gives k>=1 in (2).

Choose generic s so the original coefficient row is nonzero and R(s,w)
is nonconstant. This proves the all-degree y-chart result, including
the pure old-r direction and coefficients with arbitrary common factors.

## 3. Every non-r direction in the opposite graph (HC4MYHD2)

Now old y=s, old x=z, u=1+sz. The rows are

\[
 A_z=(u^3,3zu^2,-z^3),
\]
\[
 B_z=(s^2u(4+3sz),\ s+3zs^2(4+3sz),\ 2z-3z^2s),
\]

with cross product

\[
               (2z^2(2sz+3),\ -zu(u+1),\ su^3).             \tag{6}
\]

This is never zero for s nonzero: its last entry can vanish only at
u=0, where its first entry is nonzero. Scale v=sz and set
alpha=a, beta=sb, gamma=s^3*c. Multiplying (1) by s^2 gives

\[
 K(v)=2\alpha v^2(2v+3)-\beta v(v+1)(v+2)+\gamma(v+1)^3.    \tag{7}
\]

This is nonconstant for every nonzero row. Indeed its cubic, quadratic
and linear coefficients are the matrix

\[
 \begin{pmatrix}4&-1&1\\6&-3&3\\0&-2&3\end{pmatrix}
 (\alpha,\beta,\gamma)^T,
\]

whose determinant is -6.

If alpha=0 and beta is nonzero, no A_z can be proportional to ell:
its p entry could vanish only at u=0, where its q entry also vanishes.
Thus every root is quadratic. If alpha is nonzero, at most one finite
value z0 can have A_z proportional to ell: the ratio of its q and p
entries determines `v/(1+v)` uniquely. With tau=sz0, rescale to that
row. Equation (7) factors as

\[
 K(v)=(v-\tau)Q_\tau(v),
\]
\[
 Q_\tau(v)=(6\tau^2+9\tau+4)v^2
                  +(9\tau^2+13\tau+6)v+\tau^2.             \tag{8}
\]

The two leading coefficients cannot both vanish: twice the second
minus three times the first is -tau, and at tau=0 the first is 4.
So Q_tau is nonconstant. For tau!=0,-1,
`Q_tau(tau)=6*tau*(tau+1)^3` is nonzero, giving another, quadratic-type
root. Tau=-1 would be the excluded pure-r row. At tau=0 the row is
the pure-p row and K=2*v^2*(2v+3); v=-3/2 is quadratic.

Consequently k>=1 in every direction except alpha=beta=0. In the
exceptional pure-r direction, (7) is gamma*(v+1)^3, all of linear
type, and (2) gives exactly E=t^2. No positive coefficient is invented
for that exception.

As before, generic s avoids only finitely many coefficient degenerations;
no field-unit entry in the polynomial direction row is required. This
proves HC4MYHD2 under every actual retained-dual automorphism over K[y].

## 4. The remaining pure-r opposite route reaches the known endpoints (HC4MYHD3)

Absorb any nonzero r coefficient of the direction into the graph and
write x=R(s,r), y=s. Put R0(r)=R(0,r), g(r)=r*R0(r).

There is no relative critical point at s nonzero. On u=1+sR nonzero,
T=A eliminates p and the literal potential is

\[
 LT^2+\left(M+\frac{Ns^2(4+3sR)}{u^2}\right)T
       +Nsq+\frac{N R(2+sR)r}{u^2}.                        \tag{9}
\]

The q derivative is Ns. At u=0 the original q derivative is -2Ns,
again nonzero. At s=0 the polynomial trace coordinates
`T=p+3R0*q-R0^3*r,q,r` give

\[
                          H(T)+2Ng(r),\quad H=LT^2+MT.
\]

Thus the whole critical set is the union of the literal closed lines

\[
            s=0,\quad T=-M/(2L),\quad g'(r)=0,
                              \quad q\text{ free}.          \tag{10}
\]

The full native Hessian on such a line is exactly

\[
                              -4LN^3g''(r).                 \tag{11}
\]

To check its independence of all off-slice terms, use (9) near s=0,
where its relative coordinate Jacobian is u^3=1 on the line. The q row
has only the s cross entry N. The remaining (T,r) Hessian is diagonal
with entries 2L and 2N*g''. The bordered determinant is (11).
Nonlinear coordinate corrections vanish on the relative critical scheme,
as in [HC4MYGK0](RELATIVE_CRITICAL_LOCUS.md).

If R0=0, (10) contains a critical two-plane and unit Hessian is impossible.
If deg R0=d>=2, g' has degree d. Under any unit-Hessian retained-dual
pullback, (11) must have the same nonzero constant value on all critical
lines. A multiple root of g' is immediately excluded. At its d distinct
roots, g'' must have one common nonzero value C. But g''-C has degree
less than d, so it vanishes identically. This would make g' linear,
contradicting d>=2. Hence **all such nonlinear boundary traces forbid
unit Hessian**, in every off-slice degree.

Only two cases remain under unit Hessian:

- R0 is a nonzero constant. The whole relative critical scheme is empty.
- R0=a*r+b, a nonzero. The whole scheme is one closed A1 at s=0,
  with full native Hessian -8*L*N^3*a.

For a hypothetical nonautomorphic member, Tang rank one is impossible:
its whole critical scheme is a nonempty plane Keller fibre, and an A1
component would make its pair automorphic by HC4MYGR0. Tang rank zero
therefore gives a plane Keller pair (F,G). Every common zero contributes
one closed critical A1, since F_U,G_U cannot both vanish when J(F,G)
is a field unit, and the remaining equation is affine linear in the
two cotangent variables. Distinct zeros give disjoint components.

In the empty case (F,G) has unit ideal. In the single-line case it has
exactly one zero; its etale zero scheme is reduced, so its ideal is the
maximal ideal of that point. After a source translation this is exactly
the retained singleton-zero endpoint. These are the distinct unresolved
[unit-ideal](RELATIVE_COORDINATE_ENDPOINT.md) and
[singleton-zero](P_GRAPH_SINGLETON_ENDPOINT.md) plane problems.
No converse realization inside this pure-r subfamily or properness theorem
is supplied. A passed constant critical determinant is not sufficient.

These residual trace signatures also occur in the conditional cotangent
models of the two known endpoints. A unit-ideal plane Keller pair has
coprime restricted rows on every base line, so its cotangent trace has
generic E=t^2. For a singleton-zero pair, translate the zero to the base
origin. The restricted row gcd is linear, since the zero is unique and
the unit Jacobian makes at least one U derivative nonzero. Its cotangent
trace is polynomially equivalent to U*V with one free variable, with
E=t^2-t. This identifies the limitations of the zero/negative trace
test; it does not construct either nonautomorphic pair or realize all
such pairs inside the pure-r graph orbit.

## Scope, descent and controls

All nonexistence and invertibility arguments first hold over C. Embed
the finitely generated coefficient field of a characteristic-zero candidate
and its actual retained inverse into C; a polynomial gradient inverse
descends uniquely. The necessary residual plane endpoint is a complex
one after that passage, as in the earlier reductions.

The entire one-affine-dual-direction y family and every non-r opposite
direction are retired. More general dual dependence, other
boundary and mixing classes, uncovered embeddings and the two existing
plane endpoints remain open. No unrestricted HC4/JC2 theorem, formal
verification, independent external review or novelty priority is claimed.

[verify_directional_base.py](verify_directional_base.py) and its
[receipt](receipt-directional-base.json) check the native rows, both
cross products, the fixed-degree root exhaustion, rank-stratum identities,
pure-r critical scheme and arbitrary off-slice Hessian formula. The
unrestricted-degree and whole-fibre proofs are written. External inputs
are exactly the retained Hodge, Tang, critical-line and endpoint statements.
