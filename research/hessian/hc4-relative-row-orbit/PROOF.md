# Relative-linear row orbit, gauge and the zero-fibre properness problem

13 November 2026. The exact algebraic orbit theorem is HC4MYOR1. Its
relative-Hessian and normalization consequences are HC4MYOR2. Work over
a characteristic-zero field; for geometric statements take its algebraic
closure. No theorem below uses the Tang structural reduction.

## 1. The row and an explicit polynomial normalization

The retained [relative-linear unit-pivot construction](../../foundations/HC4_MENG_YANG_QUINTIC_SCHUR_FRONTIER.md)
gives

\[
 \psi=f(x,y)+2G_1(x,y)r+2G_2(x,y)s,\qquad
 G=\gamma C,\quad C\in\mathrm{SL}_2(K[x,y]),
\]

where, with v=xy,

\[
 \gamma=(-y P(v),x Q(v)),\quad
 P=18v^5+81v^4+120v^3+60v^2-1,\quad Q=(v+1)(v+2).
\]

Put ell=(y,2x). The matrix

\[
 C_0=\begin{pmatrix}
 (3v^2+6v+2)/2 & x^2(2v+3)\\
 (3y^2/2)(18v^4+63v^3+69v^2+21v-1)&
 36v^5+108v^4+87v^3+3v^2-3v+1
 \end{pmatrix}
\]

is the already retained HC5T1 normalizer from
[the toric descent](../../foundations/HC5_NONLINEAR_TORIC_DESCENT.md),
replayed here for the orbit calculation. It has det C0=1 and
gamma*C0=ell, as literal polynomial identities. Its
maximal total entry degree, and that of C0^-1, is ten. In particular
(gamma1,gamma2)=(x,y), not just equality of their zero sets. Every row
gamma*C generates that same maximal ideal and has Jacobian -2 at zero.

## 2. Exact classification of the whole orbit

**Theorem.**

\[
 \{\gamma C:C\in\mathrm{SL}_2(K[x,y])\}
 =\{(g_1,g_2):(g_1,g_2)=(x,y),\ J(g_1,g_2)(0)=-2\}.
\]

Only the converse requires proof. Write g=ell*A by decomposing each
polynomial without constant term into a multiple of y and a multiple of
2x. The orientation assumption gives det A(0)=1. Hence
det A-1 belongs to (x,y)=(g1,g2), so choose polynomial mu,nu with

\[
 \det A-1=\mu g_1+\nu g_2.
\]

Let k_ell=(-2x,y)^T and set

\[
 A'=A+k_{\rm ell}(\nu,-\mu).
\]

Then ell*A'=g because ell*k_ell=0. Also
adj(A)*k_ell=(-g2,g1)^T, so the rank-one determinant identity gives
det A'=det A-nu*g2-mu*g1=1. Taking C=C0*A' proves the converse.

This is constructive whenever certificates expressing x,y in (g1,g2)
are supplied. One decomposes det A-1 into x and y, substitutes those
certificates, and applies the displayed rank-one correction. It does not
assume g has an inverse map or solve the Jacobian conjecture.

## 3. All redundancy over one output row

For a row g generating (x,y), put k_g=(-g2,g1)^T. Its complete
determinant-one stabilizer is

\[
 D_h=I+h k_g g,\qquad h\in K[x,y].
\]

Indeed, if gD=g, each column of D-I is a polynomial syzygy of g.
Since gcd(g1,g2)=1, those syzygies are exactly multiples of k_g.
Thus D=I+k_g(s,t); determinant one requires -g2*s+g1*t=0, which
forces (s,t)=h*g by the same coprimality. Conversely (k_g g)^2=0,
so the displayed D_h has determinant one, inverse I-h*k_g*g, and fixes g.
If C_g is one solution of gamma*C_g=g, all solutions are C_g*D_h.

For g=ell this becomes

\[
 C_h=C_0+h k_\gamma\,\mathrm{ell},\qquad
 k_\gamma=(-\gamma_2,\gamma_1)^T,
\]

because C0*k_ell=k_gamma. For any nonzero polynomial h of total degree d,
the maximal entry degree of C_h is exactly d+12, yet gamma*C_h=ell.
The degree follows from deg C0=10 and deg(k_gamma*ell)=12, so their
highest parts cannot cancel. These arbitrarily large corrections produce
the identical cotangent potential f+2yr+4xs, whose gradient inverse is
explicit: recover x from the fourth output divided by four, y from the
third divided by two, then r,s from the other two affine equations.

Matrix coefficient degree is therefore not an intrinsic measure of the
remaining plane map. Search on the output row, or quotient this complete
stabilizer, before increasing a correction-matrix bound.

## 4. A larger normalized sufficient gate

For any C write D=C0^-1*C. Then G=ell*D, so

\[
 \deg D\le99\ \Longrightarrow\ \max(\deg G_1,\deg G_2)\le100.
\]

If G is Keller, [Moh's degree-100 theorem](https://www.math.purdue.edu/~ttm/jacobian)
makes it an automorphism. This normalized gate includes every old
degree-at-most-89 matrix C, because deg C0^-1=10. It is strictly larger:
take D=[[1,0],[x^99,1]]. Then C=C0*D has degree 109, whereas
G=(y+2x^100,2x) is already a triangular automorphism.
This is an enlarged matrix-parameter region, not a new plane degree bound.
Directly bounding deg G is stronger than bounding either matrix degree.

## 5. The actual unrestricted Keller question

For every C, J(G)(0)=-2. The cotangent determinant identity gives

\[
 \det\operatorname{Hess}\psi=16J(G)^2.
\]

Thus a nonzero constant Hessian determinant is equivalent to J(G)=-2,
and then the Hessian determinant is 64. The source substitution
z'=grad f+2*DG^T*z in the two dual variables is polynomially invertible,
since DG has a unit determinant. In these source coordinates the gradient
is (z',2G), hence, up to output exchange and scaling, G times an identity
map. It follows in both directions that grad psi is an automorphism
exactly when G is. The converse can be seen by restricting an inverse
of G times identity to the zero identity slice.

A collision G(a)=G(b), a different from b, always lifts: choose any common
first-output vector w and set each dual vector to
(2*DG^T)^-1*(w-grad f) at its corresponding base point. These are
colliding gradient points. Retaining a particular previously marked pair
imposes its additional base and dual matching conditions; the orbit
classification alone does not assert that match.

Combining this with the orbit theorem shows that the unrestricted
relative-linear counterexample problem is precisely the following
restricted plane question:

> Can a nonautomorphic plane Keller map have a single reduced fibre at
> zero, with that fibre the origin and Jacobian normalized to -2?

For a Keller map, every finite geometric source fibre is reduced. Thus
the ideal condition (G1,G2)=(x,y) is exactly the singleton-zero condition.
Neither the orbit theorem nor that one fibre establishes invertibility.
No equivalence of this restricted question to unrestricted JC2 is claimed.

## 6. Properness and the exact missing length

Let G in the preceding class have generic field degree d. Normalize the
target A2 in its source function field, giving pi:N->A2. Polynomial
rings over fields have finite normalization in a finite extension.
[Zariski's Main Theorem](https://stacks.math.columbia.edu/tag/02LR)
identifies the source A2 with an open subset of N. A normal surface is
Cohen–Macaulay, so the finite map to the smooth surface is flat by the
[flatness criterion](https://stacks.math.columbia.edu/tag/00R3).
Its finite fibres all have scheme length d.

Exactly one reduced point of N_0 belongs to the affine source. The
remaining finite-normalization fibre therefore has total length d-1.
This is a scheme-length statement about the omitted part of the fibre,
not a ramification-only classification or a logarithmic point correction.

G is an automorphism if and only if it is proper over some neighborhood
of zero. For the nontrivial direction, proper and quasi-finite makes G
finite over that neighborhood. There it is finite etale of degree d;
its singleton fibre forces d=1. The birational Keller case is invertible
(or use the open-immersion form of Zariski's Main Theorem and injective
polynomial-map invertibility). The converse is immediate for an automorphism.

Consequently any relative-linear counterexample must have zero in its
nonproperness locus and missing normalization length d-1 there. This is
the concrete global boundary obligation. A one-point affine fibre does
not imply that its finite-normalization fibre is complete.

## 7. A verified source-A3 control with a singleton zero fibre

Use the known [foundational Keller map](../../foundations/verified/FOUNDATIONAL_GEOMETRY.md),
writing u=1+xy:

\[
 F=(u^3z+y^2u(4+3xy),\quad
 y+3xu^2z+3xy^2(4+3xy),\quad 2x-3x^2y-x^3z).
\]

Its determinant is -2. The finite cubic cover of the target (a,b,c) is
the root hypersurface in A3 times P1,

\[
 cT^3-2T^2S+bTS^2-2aS^3=0.                         \tag{1}
\]

This is the retained [marked-root model](../../foundations/verified/MARKED_ROOT_MODEL.md),
used here as a properness control. It is finite flat of degree three:
every fibre is a nonzero binary
cubic divisor on P1, of length three, and the coefficient of T^2S is
the fixed nonzero scalar -2. Both affine charts are smooth: on S=1
solve for a, and on T=1 solve for c. The two irreducible charts overlap,
so the cover is integral and normal.

The source maps into it by (a,b,c)=F and [T:S]=[u:x]. These homogeneous
coordinates never vanish simultaneously, since u-xy=1. The cubic's
two homogeneous partials there are 2x and -2u, so this source point is
always in the etale locus. Conversely the following formulas recover
the source from each etale chart, and their substitutions are checked.

On S=1, t=T, set D=3ct^2-4t+b and require D nonzero. Then

\[
 x=2/D,\qquad y=t-D/2,\qquad
 z=(2x-3x^2y-c)/x^3.
\]

On T=1, s=S, set E=-2+2bs-6as^2 and require E nonzero. Then

\[
 x=-2s/E,\quad y=b-3as,\quad
 z=-aE^3/8-y^2E^2/4+3y^2E/2.
\]

Thus the source is exactly the etale open of this finite normal cubic
cover, whose generic degree is three. Over zero, (1) becomes -2T^2S.
The simple point [1:0] gives the source origin; the double point [0:1]
is omitted. The source zero fibre is a single reduced point, while two
units of normalization length are missing.

There is a separate literal maximal-ideal check. Set
H=y(1+3xy)+xuz. Direct multiplication gives

\[
 y=(-3x+9x^2H/2)F_1+
       (u-x(F_2+y)/2)F_2.
\]

Also F1=z+y*K, with K=x(u^2+u+1)z+yu(4+3xy), and
2x=F3+3x^2y+x^3z. Substituting the displayed y identity gives literal
polynomial certificates for x,y,z in (F1,F2,F3). Thus the zero-fibre
ideal really is (x,y,z), independently of the cover interpretation.

There is even an explicit escaping-sheet family. For t nonzero the
complete fibre over (-t^2,0,0) is

\[
 (0,0,-t^2),\quad (-1/(2t),3t,26t^2),\quad
 (1/(2t),-3t,26t^2).
\]

Two points escape as t tends to zero and only the origin remains. This
is a control in dimension three, not a plane counterexample or a proof
that the plane singleton-zero question has either answer. It explains
why the global properness obligation cannot be replaced by local ideal
equality or source etaleness.

## Scope and reproduction

The [checker](verify.py) constructs C0, verifies the universal determinant
correction and stabilizer, exercises a constructive nontrivial orbit lift,
checks the high-degree gauge and normalized degree example, and replays
the cubic-cover charts and escaping fibres. The orbit completeness and
normal-surface flatness applications are written proofs with identified
imports. There is no finite search presented as an all-degree theorem.
The unrestricted singleton-zero plane gate remains open.
