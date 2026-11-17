# Homogeneous sweeps: the axis cover and exceptional11 decomposition

This studies every polynomial normalized sweep whose matrix
P=Dg-p df is homogeneous of a single positive ordinary degree n.
The main statements assume r>=2. They use the full differential
reconstruction of KMFI115, not a triangular ansatz.
After constant target translations,
f is homogeneous of degree n and g of degree n+1.
Write
\[
 H=g-fp,\quad T_A=H+Ap,\quad
 S(\gamma,p)=(\gamma+f,\ g+\gamma p),\quad
                        \det DS=\gamma^r.                 \tag{1}
\]
Work geometrically over algebraically closed characteristic-zero
constants. The axis is
\[
                \mathcal A=\{(A,0,\ldots,0):A\ne0\}.
\]

## 1. The whole actual axis fiber and global nonfiniteness

**KMFI118.** For every A!=0 the actual source fiber is finite,
reduced and unramified. It has cardinality
\[
                       k=1+nm,\qquad m\ge0,               \tag{2}
\]
independent of A. There is one distinguished point p0,gammaA.
Every other point has
\[
                   f=(n+1)A,\qquad\gamma=-nA.             \tag{3}
\]
For r>=2 every such nonzero axis point also has
rankP<=r-2 and rankDX=rankP+1<=r-1.
Thus positive axis orbits lie on the rank-drop locus of X,
even though S itself is etale there.
The source preimage of the axis is a finite etale cover consisting
geometrically of one section and m connected cyclic covers of
degree n. In addition, for every r>=2, this homogeneous sweep
is globally nonfinite. Rank1 is an essential comparison exception when f and g
are explicitly assumed homogeneous of degrees n and n+1;
its zero matrix P alone does not determine those degrees.

Here is the fiber proof. Euler homogeneity gives
\[
                         Pp=(n+1)H+fp.                    \tag{4}
\]
On T_A=0 this becomes Pp=(f-(n+1)A)p.
For nonzero p, nilpotency forces the eigenvalue to vanish,
giving (3). At the origin f0 and gammaA.
The exact determinant DT_A=(A-f)^r is therefore nonzero at
every root. Hence the zero scheme of T_A is zero-dimensional
and reduced; it is finite because it is affine of finite type.

To see the rank-drop assertion, (4) makes p a kernel vector
at a nonzero axis root, and df(p)=n f=n(n+1)A!=0.
If P were regular nilpotent there, KMFI116 would make X full rank.
It would then force D_pP to have nonzero eigenvalue df(p).
But homogeneity gives D_pP=nP, which is nilpotent.
Hence P is not regular and rankP<=r-2.
Since df is nonzero on kerP, adding its row raises the rank
by exactly one, proving the stated rankDX.

Scalar multiplication by n-th roots of unity preserves T_A=0
and acts freely on nonzero roots. This proves (2).
After the finite etale base change A=s^n, p=sq identifies
every fiber with T_1=0. Descent gives a finite etale cover
of the axis. Each free root-of-unity orbit gives exactly one
connected cyclic degree-n component.

For global nonfiniteness, suppose that S were finite.
Its origin fiber identifies with H^{-1}(0) by gamma=-f.
Since H is homogeneous of degree n+1, finiteness would force
H^{-1}(0) to be just the origin. Its leading forms would
then have no projective common zero.

Bezout now gives exactly (n+1)^r roots of T_1=0, all simple
by the preceding argument. The origin has determinant1;
each other root has determinant(-n)^r.
The Euler--Jacobi identity with constant numerator1 gives
\[
                     1+\frac{(n+1)^r-1}{(-n)^r}=0.        \tag{5}
\]
We use [Balanov--Kononovich--Krasnov, Theorem4.3,
printed page586](https://doi.org/10.1017/S0013091511000435);
its maximal simple zero-set hypothesis has just been proved
under the assumed absence of projective common zeros.
For even r the left side of (5) is positive.
For odd r>=3, (n+1)^r-1>n^r, so it is negative.
Both contradict zero. Thus S is nonfinite.
For r1, (5) holds identically, and the ordinary finite
rank1 examples show why that case cannot be included.

The Euler--Jacobi step is conditional on the hypothetical absence
of zeros at infinity. It is not applied to a nonproper system
by counting only its finite roots.

## 2. The complete positive-axis-fiber table for exceptional11

Assume now rankr4, full rankDX4, and either exceptional degree11
PSL2(11) action. By KMFI117, n>=2.
Let Z be the finite normalization of the target in the source field,
and take the Galois normalization as well. At the generic axis
point, write D_ax for decomposition and I_ax for inertia:
I_ax is the kernel of the action of D_ax on the residue field.
This is a codimension-four specialization, not a DVR notation
with ramification indices silently imported.

**KMFI119.** If m>=1 in (2), there are exactly four possible
decomposition rows:
\[
\begin{array}{c|c|c|c|c}
 n&D_{\rm ax}&I_{\rm ax}&
 |\operatorname{Fix}(I_{\rm ax})|&\text{actual }k\\ \hline
 2&C2&1&11&3,5,7,9\\
 2&V4&C2&3&3\\
 3&C3&1&11&4,7,10\\
 5&C5&1&11&6,11
\end{array}                                                \tag{6}
\]
These are necessary possibilities, not assertions that a
homogeneous PSL input realizing any row exists.
If n is not2,3 or5, the actual axis fiber must consist of
the distinguished point alone.

We prove the residue and subgroup restrictions explicitly.
The scaling
\[
 (\gamma,p)\mapsto(s^n\gamma,sp),\qquad
 (Y0,Y_{\rm vec})\mapsto(s^nY0,s^{n+1}Y_{\rm vec})          \tag{7}
\]
lifts to the normalizations. For the Galois closure, use the
diagonal action on the ordered fiber product and normalize:
the connected torus preserves each relevant irreducible component.

Normalize any component over the axis. Its torus orbit is
one-dimensional and covers the entire axis. The stabilizer is
a subgroup of the target isotropy mu_n, so its residue cover
is cyclic of degree q dividing n. Thus D_ax/I_ax=C_q.
The distinguished section is an actual unramified rational sheet,
so D_ax fixes it. Conjugate so that D_ax is a subgroup of its
A5 sheet stabilizer.

If m>=1, an actual degree-n axis component is cyclic and etale,
so D_ax/I_ax surjects C_n. Since q divides n, q=n.
Every actual axis point is fixed by I_ax; hence its common
fixed-sheet mark is at least1+n.

The complete A5 subgroup list is
1,C2,C3,C5,V4,S3,D10,A4,A5.
The normal subgroups with nontrivial cyclic quotient give:
\[
 (D,I,n)=(C2,1,2),(V4,C2,2),(S3,C3,2),(D10,C5,2),
                 (C3,1,3),(A4,V4,3),(C5,1,5).              \tag{8}
\]
The S3 row has inertia mark2<3, the D10 row mark1<3,
and the A4 row mark1<4; all are impossible.
For the remaining rows the action on Fix(I) has respectively
four2-cycles and three fixed sheets; one2-cycle and one fixed
sheet; three3-cycles and two fixed sheets; two5-cycles and one
fixed sheet. The actual cover retains exactly the distinguished
fixed sheet and some positive number of the n-cycles.
This proves every cardinality in (6).

The checker independently enumerates all59 subgroups of each
actual A5 stabilizer, every normal subgroup and cyclic quotient.
There are67 nontrivial cyclic quotient pairs per stabilizer.
The mark test retains exactly46 pairs in the four rows (6);
the other21 are precisely the S3,D10,A4 exclusions.
Both exceptional actions are checked in their actual11-letter
representations. This is not inferred from element orders alone.

If I_ax1, the full normalization is etale at the geometric
generic axis point. Missing sheets there are omitted unramified
normalization points. In the V4/C2 row, I_ax has three fixed
sheets and four2-orbits. All three unramified points are retained;
the other four normalization points are absent, with total
local generic rank8. That is four omitted points, not eight.

In every positive-axis-fiber row, I_ax is1 orC2.
Therefore the axis is not contained in the original visible
branch hypersurface, whose meridian has order5 by KMFI105.
The V4/C2 row can lie on a different branch component supplied
only by missing ramification. No nonempty source ramification
above the axis has been inferred.

## 3. When the axis is actually proper

**KMFI120, first assertion.** For a homogeneous sweep of generic
degree N, k<=N and equality holds exactly when S is finite etale
on a target neighborhood of the axis. Here k counts actual points,
not the length of a possibly nonflat normalization fiber.

At a geometric axis point all source points are etale and
embed as open normalization points by Zariski's Main Theorem.
After strict henselization each contributes a rank-one local
factor of the finite normalization. Any additional nonzero factor
is torsion-free over the regular target local ring and has positive
generic rank. Their ranks sum to N, so k<=N.
If k=N, no extra factor remains, so every normalization point
is present and unramified. The finite image of the omitted boundary
misses the axis, giving a neighborhood where S is the finite
etale normalization. The converse is immediate.
This argument does not assume flatness of an arbitrary finite
normalization at a codimension-four point.

Consequently a homogeneous exceptional11 input finite near the
axis must have exactly
\[
                 n=5,\quad k=11,\quad
                      D_{\rm ax}=C5,\ I_{\rm ax}=1.        \tag{9}
\]
Its homogeneous degrees are f5 and g6.
It is still globally nonfinite by KMFI118.
Its ten extra axis points must lie on the rank-drop locus of X;
a candidate with full rankDX everywhere on f!=0 cannot use this row.

Every other row in (6) forces omission along the axis.
When k1, the omitted local factors have total generic rank10;
this is not a count of divisors, primes or geometric points.
No positive-axis-fiber condition is automatic: the following
full-rank rational control has k1.

## 4. A rational full-rank homogeneous control with only one axis point

**KMFI120, second assertion.** The k1 alternative is real even
for full-rank inputs. For every r>=2 let e=r+1 and
\[
 d_r(T)=\frac{T^{r+1}-(T-1)^{r+1}}{r+1}.
\]
This is monic, degree r, separable, and has neither0 nor1 as a root.
Let M be its rational companion multiplication matrix.
In the finite etale algebra Q[T]/d_r, put
\[
 Q_p=p1I+p2M+\cdots+prM^{r-1},\qquad F(p)=\det Q_p,
\]
\[
               f=-F,\qquad g=F(M-I)p.                     \tag{10}
\]
These rational polynomials give a full-rank normalized input
with P homogeneous of degree r. Its exact geometric action is
regular cyclic C_(r+1), so it has no Keller completion even after
independent stabilization and arbitrary birational recharts.
Arithmetic over Q is C_e semidirect Gal(Q(zeta_e)/Q), with exact
constants Q(zeta_e). This is a countercontrol, not a new Keller
action.

To verify nilpotency, split the etale algebra and let a_i be the
distinct roots of d_r. In the corresponding linear coordinates z_i,
F=product z_i and H=F diag(a_i)z.
On the dense torus, P is similar to F(B-I), where
\[
                      B=\operatorname{diag}(a_i)+a\,1^t.
\]
The determinant lemma gives
\[
 \det(TI-B)=(r+1)d_r(T)-T d_r'(T)=(T-1)^r.
\]
Thus P is polynomial and nilpotent everywhere.
Since every a_i-1 is nonzero, the rank-one perturbation B-I
has nullity at most1, hence exactly1.
A kernel vector has nonzero coordinate sum; otherwise all its
entries would vanish. Consequently df is nonzero on kerP on
F!=0. KMFI116 proves full rank there.

The exact field is also simple. Set A=Y0, Z=F(p), B0=F(Y_vec).
Then
\[
 Y_{\rm vec}=(AI+ZM)p,\qquad
 Z\det(AI+ZM)=\frac{(A+Z)^{r+1}-A^{r+1}}{r+1}.
\]
Therefore
\[
             \gamma^{r+1}=A^{r+1}+(r+1)B0,\qquad
              p=(AI+(\gamma-A)M)^{-1}Y_{\rm vec}.          \tag{11}
\]
These recover both fields exactly. Over algebraically closed
constants F(Y_vec) splits into a product of independent linear
coordinates, so the radicand has a simple zero divisor and is
not a proper power. Equation(11) has degree e and the stated
cyclic group. Its splitting field adds precisely zeta_e,
giving the arithmetic constants and semidirect product.

Finally the entire axis ideal has an exact certificate.
Let V=T_1=p+FMp. Since M e1=e2 and Q_p e1=p,
\[
            \operatorname{adj}(Q_p)V=F e1+F^2 e2.          \tag{12}
\]
Its first coordinate puts F in the ideal generated by V.
Then p=V-FMp puts every p_i in that ideal. The reverse inclusion
is immediate, so (V1,...,Vr)=(p1,...,pr) exactly.
Thus k1, not merely a numerically observed fiber count.
The origin fiber contains gamma0,F0 and is positive-dimensional.

For r4, F has32 rational monomials; the cyclic coframe determinant
at p=(1,0,0,0) is1/125. The geometric group is C5, generic degree5,
and four generic sheets are absent over the axis.
This demonstrates that full rank does not force m>=1 and that
a one-point special fiber does not prove generic degree1.

## 5. Scope and replay

The homogeneous matrix hypothesis is essential. The statements
do not apply to every nonhomogeneous nonlinear matrix field
left by KMFI117. Even inside this class, the k1 PSL alternative
is not excluded. If k>1, the four rows are necessary possibilities
and not actual PSL realizations. Finite behavior near the axis
is distinct from global finiteness.

The written proofs cover all r,n in their stated ranges.
Exact checks include homogeneous Euler identities, axis ideals,
the rational full-rank norm construction, both complete A5
normal-subgroup quotient tables, their fixed-sheet actions and
primary-source bytes/hash. They do not classify all homogeneous
polynomial solutions or prove a generic action from a passed
axis test. Formal verification, independent whole-proof review
and external review are not claimed.

[Helper](homogeneous_axis.py) ·
[checker](verify_homogeneous_axis.py) ·
[model](homogeneous_axis_model.json) ·
[receipt](homogeneous_axis_receipt.json) ·
[pinned Euler--Jacobi input](sources/homogeneous-axis-manifest.json).

The nonfinite rank4 PSL2(11) field, the absolute exceptional11
action and the complete faithful-transitive Keller spectrum
remain open.
