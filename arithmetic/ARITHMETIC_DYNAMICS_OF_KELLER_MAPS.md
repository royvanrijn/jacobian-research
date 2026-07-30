# Programme 9: arithmetic dynamics of Keller maps

This note starts the iteration programme with the foundational map.  It keeps
three levels separate:

- the exact degree, generic-preimage, and exceptional-orbit results proved
  below;
- exact first computations that define the next problems; and
- claims in the recent external record that have not been imported as
  repository theorems.

Write

\[
\begin{aligned}
F_1&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
F_2&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
F_3&=2x-3x^2y-x^3z.
\end{aligned}
\]

The foundational results used here are the [determinant and
collision](../verified/FOUNDATIONAL_GEOMETRY.md), the
[marked-root inverse equation](../verified/MARKED_ROOT_MODEL.md), and the
[exact image and nonproperness set](../verified/IMAGE_AND_NONPROPERNESS.md).

## 1. Exact degree growth

Let

\[
(\alpha_n,\beta_n,\gamma_n)
=\bigl(\deg (F^n)_1,\deg (F^n)_2,\deg (F^n)_3\bigr),
\qquad
(\alpha_0,\beta_0,\gamma_0)=(1,1,1).
\]

### Theorem 1

For every \(n\geq 0\),

\[
\begin{pmatrix}\alpha_{n+1}\\ \beta_{n+1}\\ \gamma_{n+1}\end{pmatrix}
=
\begin{pmatrix}
3&3&1\\
3&2&1\\
3&0&1
\end{pmatrix}
\begin{pmatrix}\alpha_n\\ \beta_n\\ \gamma_n\end{pmatrix}.
\tag{1}
\]

In particular, if \(d_n=\deg F^n=\alpha_n\), then

\[
d_0=1,\qquad d_1=7,\qquad
d_{n+2}=6d_{n+1}+d_n,
\tag{2}
\]

and

\[
d_n=
\frac{(4+\sqrt {10})(3+\sqrt {10})^n
(\sqrt {10}-4)(3-\sqrt {10})^n}{2\sqrt {10}}.
\tag{3}
\]

Consequently the first dynamical degree is

\[
\lambda_1(F)=3+\sqrt {10}.
\]

Thus the extension of \(F\) to the standard projective compactification is
not algebraically stable: \(\deg F^n\neq 7^n\).

### Proof

Suppose a polynomial triple \(G=(A,B,C)\) has component degrees
\((\alpha,\beta,\gamma)\) and put
\(\Delta=\alpha-\beta+\gamma\).  In \(F_1\circ G\), the candidate of largest
degree from \((1+AB)^3C\) has degree
\(3\alpha+3\beta+\gamma\); its only competing top candidate, from
\(3A^2B^4\), has degree \(2\alpha+4\beta\).  Their difference is
\(\Delta\).  The corresponding comparisons in the other coordinates are

\[
\begin{array}{c|c|c|c}
&\text{selected degree}&\text{competitor}&\text{difference}\\ \hline
F_1\circ G&3\alpha+3\beta+\gamma&2\alpha+4\beta&\Delta\\
F_2\circ G&3\alpha+2\beta+\gamma&2\alpha+3\beta&\Delta\\
F_3\circ G&3\alpha+\gamma&2\alpha+\beta&\Delta.
\end{array}
\]

At \(n=0\), \(\Delta_0=1\).  If the selected terms win at stage \(n\), then

\[
\Delta_{n+1}
=\alpha_{n+1}-\beta_{n+1}+\gamma_{n+1}
=3\alpha_n+\beta_n+\gamma_n>0.
\]

They therefore win uniquely at every stage.  Their leading homogeneous
forms are nonzero products, so cancellation is impossible.  This proves
(1).  The characteristic polynomial of the matrix in (1) is
\(T(T^2-6T-1)\), which gives (2), (3), and the claimed limit. \(\square\)

The first six multidegrees are

\[
(7,6,4),\ (43,37,25),\ (265,228,154),\
(1633,1405,949),\ (10063,8658,5848),\
(62011,53353,36037).
\]

The proof is all-order; these values are checks, not the evidence on which
the theorem rests.

## 2. Generic preimage trees and field multiplicity

Let \(K=k(x,y,z)\) in characteristic zero and let
\(\sigma=F^\ast:K\hookrightarrow K\).  The marked-root cubic gives
\([K:\sigma K]=3\).  Since \(\sigma^i\) identifies the inclusion
\(\sigma K\subset K\) with
\(\sigma^{i+1}K\subset\sigma^iK\),

\[
[K:\sigma^nK]
=\prod_{i=0}^{n-1}[\sigma^iK:\sigma^{i+1}K]
=3^n.
\tag{4}
\]

### Corollary 2

The generic degree of \(F^n\) is exactly \(3^n\).  Moreover

\[
\det D(F^n)=(-2)^n,
\]

so its geometric generic fiber is reduced and has \(3^n\) points.  The
geometric generic preimage tree through level \(n\) is therefore the regular
rooted ternary tree.

Let \(H_n\) be the Galois group of the normal closure of the \(n\)-th generic
fiber.  Parent maps in the tree give

\[
H_n\hookrightarrow
\operatorname{Aut}(T_{3,n})
\cong S_3\wr S_3\wr\cdots\wr S_3,
\tag{5}
\]

and the level-one quotient is \(H_1=S_3\).  Equality in (5) is **open**.
This is the first precise form of “variation of fiber Galois groups under
iteration”: prove or disprove maximal iterated monodromy, then determine the
special-target drop loci.

## 3. An exceptional arithmetic and escaping orbit

Let

\[
L_x=V(y,z),\qquad L_z=V(x,y).
\]

Direct substitution gives

\[
F(q,0,0)=(0,0,2q),\qquad
F(0,0,q)=(q,0,0).
\tag{6}
\]

Hence \(L_x\cup L_z\) is invariant, the two components are exchanged, and
\(F^2\) acts as multiplication by \(2\) on each component.  The origin is
fixed.  Every other complex point on this cross escapes in the usual
archimedean norm.  For \(q\in\mathbb Q^\times\), the multiplicative and
logarithmic heights satisfy

\[
H(F^n(q,0,0))=2^{n/2+O(1)},\qquad
h(F^n(q,0,0))=\frac n2\log 2+O(1),
\tag{7}
\]

and the same holds on \(L_z\).  In the standard arithmetic-dynamics
convention
\(\alpha_F(P)=\lim h^+(F^n(P))^{1/n}\), these nonzero rational points have
arithmetic degree \(1\).  This is an explicit exceptional set on which height
growth is far smaller than the ambient dynamical degree \(3+\sqrt {10}\).

At the fixed point, the derivative has characteristic polynomial

\[
(\lambda-1)(\lambda^2-2),
\]

consistent with the invariant cross and its \(F^2\)-multiplier \(2\).

## 4. Boundary dynamics without affine critical points

The affine critical locus is empty for every iterate, but the
nonproperness hypersurface

\[
D=27a^2c^2-18abc+16a+b^3c-b^2=0
\]

still controls the loss of sheets at infinity.  Its first pullback is the
degree-eight polynomial

\[
\begin{aligned}
D\circ F={}&-9x^4y^2z^2-54x^3y^3z-18x^3yz^2-81x^2y^4\\
&-72x^2y^2z-9x^2z^2-54xy^3+6xyz+63y^2+16z.
\end{aligned}
\tag{8}
\]

This exact degree \(8\), far below the naive bound \(28\), is the first
signal that boundary iteration has its own cancellation dynamics.  No
factorization, invariance, or finiteness statement is inferred from (8).

A useful replacement for postcritical finiteness is the
**post-nonproper set**

\[
\mathcal P_{\mathrm{np}}(F)
=\overline{\bigcup_{j\geq0}F^j(V(D))}.
\tag{9}
\]

The immediate questions are whether (9) is a proper algebraic subset,
whether its finite truncations equal the nonproperness loci of the iterates,
and whether the degree-eight pullback in (8) belongs to a finite orbit of
boundary divisor classes on a suitable compactification.

Over \(\mathbb R\), one step has three sheets on \(D<0\) and one on \(D>0\).
For iteration, the sheet number is obtained by summing the sheet numbers of
all real ancestors.  Therefore the real bifurcation walls at level \(n\)
must be studied together with the forward/backward boundary orbit, not from
the affine Jacobian, which never vanishes.

## 5. Research queue

The next calculations are deliberately ordered so that each can end in an
exact certificate.

1. **Compactified degree proof.** Resolve the indeterminacy of the extension
   to \(\mathbb P^3\) and realize matrix (1) as pullback on a finite-rank
   divisor lattice.  Then compute the remaining dynamical degrees.
2. **Iterated monodromy.** Compute the normal closure at level two and test
   whether \(H_2=S_3\wr S_3\).  Express failure, if any, as a square-class or
   discriminant dependence among the three first-level branches.
3. **Special preimage trees.** Use the cubic marked-root equation recursively
   to classify targets with branching profiles below \(3^n\), distinguishing
   affine ancestors from length lost at the boundary.
4. **Height growth.** Search for invariant curves and surfaces, construct
   canonical-height candidates at scale \(3+\sqrt {10}\), and classify
   height-growth exceptions beginning with (6).
5. **Real sheets.** Pull back the sign condition for \(D\) through each real
   inverse branch and enumerate chambers for \(F^2\) exactly.
6. **Escaping sets.** Separate forward escaping orbits such as (6) from
   nonproper inverse branches approaching \(D=0\); these are different
   notions and should not share terminology in certificates.
7. **Invariant hypersurfaces.** Search degree by degree for irreducible
   \(H\) with \(H\mid H\circ F\), recording bounded null results as
   computations rather than proofs of nonexistence.
8. **Post-nonproper finiteness.** Compute the first images of \(V(D)\) and
   compare them on a resolved compactification.

## 6. External-work audit

The recent record is Liam Giannini,
[The Alpöge--Fable counterexample to the Jacobian Conjecture: verification,
geometry, dynamics, and an equivariant program for the plane
case](https://zenodo.org/records/21461572), Zenodo record `21461572`,
version `1.1`, created 20 July 2026.  The pinned files reported by the Zenodo
API are:

| file | Zenodo MD5 | locally computed SHA-256 |
|---|---|---|
| `jacobian-note-v1.1.pdf` | `8225a47137ed62b66ae168ae5db0245d` | `fab3ebe9e321ec8c370ae496157a92619424424ed757b3c4a86cc373659cbfb3` |
| `m-v1.1.zip` | `4c789ab70f30779c887eff9f02f2b1b7` | `e688123884caa2ab5286a8c8a3ff7167afe1ce570be1127c7ea3d7eb1598b134` |

That paper reports the six degrees displayed above and presents the
recurrence and \(3+\sqrt {10}\) as a conjecture.  The bundled
`scripts/dyn_e2.py` in version 1.1 actually loops through five iterates, so
the checked-in external script does not by itself reproduce the sixth
reported value.  Theorem 1 independently proves the recurrence and every
listed value without generic-line specialization.

The same record claims stability of the affine image under iteration.  That
claim is not used in Theorems 1 or 2 and is not registered here as proved;
its resultant case split should receive a separate denominator-safe audit
before import.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_foundational_arithmetic_dynamics.py
```

The checker verifies the top homogeneous terms, the induction matrix and its
characteristic polynomial, literal degrees through \(F^2\), the first twelve
matrix iterates, the invariant cross, the Jacobian input, and (8).  The
field-tower multiplication in (4) is a written proof.
