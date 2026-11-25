# Repeated-cubic even quintic faces and the residual quartic-CM sixfold

Use the rational sphere functional
\[
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^{1}
 f(x,(1-t^2)/(2x),t)\,dt
\]
over a characteristic-zero field. The integral means its rational polynomial rule. A phase is eventually pure if all sufficiently late positive moments vanish. Consider
\[
 f=P(y,t)+xt^2(q_0+Bt^2),\quad
 P=\sum_{j+r\le5,\ r\ {\rm even}}p_{jr}y^jt^r,\quad
 A=p_{14},\ C=p_{32},\ D=p_{50},\quad BD\ne0.             \tag{1}
\]
All lower coefficients are arbitrary. Set
\[
 R(z)=B-\tfrac A2z-\tfrac C8z^2-\tfrac D{32}z^3,\qquad
 Q(z)=z^3R(1/z)=Bz^3-\tfrac A2z^2-\tfrac C8z-\tfrac D{32},
\]
\[
 \Delta=4A^2C^2+8BC^3-16A^3D-27B^2D^2-36ABCD
           =1024\operatorname{Disc}(Q).                \tag{2}
\]

**Theorem 1.** For EVERY characteristic-zero field, every phase (1) with
\(\Delta=0\) has infinitely many nonzero positive moments. This is an
all-member repeated-cubic exclusion, with no generic-End or simplicity assumption.

**Theorem 2.** Let J be a polarized abelian sixfold over a number field with
a faithful geometric \(E=\mathbb Q(\zeta_5)\)-action, Rosati conjugation,
and E-signature \((1,1,2,2)\) in character order1,2,3,4, up to relabeling.
If J is not absolutely simple with FULL geometric End^0(J)=E, it has
positive-density ordinary primes after finite extension.

The ordinary inputs are exactly the existing all-special quartic-CM
fourfold theorem, balanced-fourfold theorem, surface theorem and
potential-CM theorem in Section8. Neither result proves ordinary primes
for the remaining simple full-End sixfold.

## 1. Exact face and prime isolation

For \(x^it^h\) assign weight i+h. The top term of \(y^jt^r\) after Laurent
substitution is \((-1/2)^jx^{-j}t^{2j+r}\), weight j+r; its other terms
have smaller weight. Thus all weights in f are at most5 and its equality
face is
\[
                       xt^4R(t^2/x^2).                 \tag{3}
\]
Over a number field omit the fixed finite nonintegrality primes.
At \(p=10k+1,m=2k\), angular-constant terms in \(f^m\) have t-height at
most \(5m=p-1\). Multiplication of the rational moment by p kills all
lower even heights modulo the prime; odd heights integrate to zero.
Angular balance in (3) requires z-degree k. At sufficiently large good
absolute-degree-one primes,
\[
        p\Lambda(f^{2k})\bmod\mathfrak p
        =H_p,\qquad H_p=[z^k]R(z)^{2k}.                \tag{4}
\]
No lower coefficient contributes. Eventual purity forces every sufficiently
late H_p to vanish.

## 2. Normalized covers and the exact selected line

Take smooth projective normalizations
\[
 {\cal C}:w^5=Q(u^2),\qquad
 \overline{\cal C}:w^5=Q(z),\quad z=u^2.                 \tag{5}
\]
All nonzero cubic roots have multiplicity at most3, so these covers are
connected. The involution \(\iota(u,w)=(-u,w)\) commutes with the deck
action and has the displayed quotient.

For a fifth-root cover with positive raw branch exponents \(a_\nu\) summing
to \(5l\), the character-i eigensheaf is
\[
 {\cal O}(-il+\sum_\nu\lfloor ia_\nu/5\rfloor).          \tag{6}
\]
Locally divide \(w^i\) by the branch parameter to
\(\lfloor ia_\nu/5\rfloor\) to obtain an integral generator. These generators
give (6), and \(h^1({\cal O}(-n))=\max(n-1,0)\).
For (5), use homogeneous degree10 upstairs and degree5 downstairs:
the raw infinity exponents are4 and2. Thus, in character order1,2,3,4,
\[
\begin{array}{c|c|c|c|c}
Q\text{ roots}&{\cal C}\text{ exponents}&h^1({\cal C})_i&
 h^1(\overline{\cal C})_i&(\dim J^-,\operatorname{sig}_EJ^-)\\ \hline
1+1+1&(1,1,1,1,1,1,4)&(1,2,3,4)&(0,1,1,2)&(6,(1,1,2,2))\\
2+1&(2,2,1,1,4)&(1,2,1,2)&(0,1,0,1)&(4,(1,1,1,1))\\
3&(3,3,4)&(1,0,1,0)&(0,0,0,0)&(2,(1,0,1,0)).
\end{array}                                                   \tag{7}
\]
Here \(J^-=\operatorname{im}(1-\iota)\) up to isogeny. The plus cohomology
is the quotient pullback, so the last column is subtraction of the
preceding columns. Riemann--Hurwitz gives upstairs genera10,6,2 and
downstairs4,2,0, independently checking the totals. For example the double
row uses \(2i-1-2\lfloor2i/5\rfloor-\lfloor4i/5\rfloor\) upstairs and
\(i-1-2\lfloor2i/5\rfloor\) downstairs.

The restricted Jacobian polarization on \(J^-\) has Rosati conjugation on E:
deck transformations preserve it. The E-action is faithful because a
nontrivial order-five deck character occurs and \(\Phi_5\) is irreducible.

Every raw exponent upstairs is less than5. Its character-one eigensheaf is
therefore exactly O(-2), with Cech generator w/u. It changes sign under
\(\iota\), so the ENTIRE selected line is in the Prym cohomology, not the
quotient's. At good absolute-degree-one primes of a field containing E,
\(p=10k+1\), and Frobenius gives the scalar
\[
 [u^{p-1}]Q(u^2)^{(p-1)/5}
   =[z^{5k}]Q(z)^{2k}
   =[z^k]R(z)^{2k}=H_p.                                \tag{8}
\]
The last equality is exact polynomial reversal, valid at repeated roots
too. The Cech computation uses the normalized eigensheaf (6).
Good reduction of the chosen normalized branch data is secured outside
finitely many primes. Vanishing of H_p kills this nonzero invariant line.
An ordinary abelian variety has bijective semilinear Frobenius on H^1(O);
hence \(J^-\) is then nonordinary. The fixed isogeny splitting is harmless
outside its degree.

## 3. Repeated-root exclusion over every characteristic-zero field

In the double-plus-simple case write
\(Q=B(z-r)^2(z-s)\), with \(Brs(r-s)\ne0\).
The second row of (7) is a polarized balanced quartic-CM fourfold at EVERY
member. Apply ED0-BALANCED-FOURFOLD-POSITIVE-ORDINARY: the signature is
(1,1,1,1), the Rosati action is conjugation, and no full-End or simplicity
condition is needed. Its explicit quaternion construction produces a
rational imaginary quadratic subfield of signature(2,2), so Fite Theorem4
gives positive ordinary density after extension.

For number-field coefficients extend once to define E, the roots, factors
and relevant endomorphisms, and reapply the ordinary theorem over that
final field. An eventual zero tail and (4),(8) would make the Prym
nonordinary at every sufficiently large absolute-degree-one prime.
Prime ideals of higher absolute degree have density zero, a contradiction.

For a triple root \(Q=B(z-r)^3\), \(Br\ne0\),
\(R=B(1-rz)^3\) and
\[
          H_p=B^{2k}\binom{6k}{k}(-r)^k\ne0\bmod p.     \tag{9}
\]
Indeed \(6k<p\) and B,r are units at all sufficiently large good primes.
After adjoining E there are infinitely many absolute-degree-one primes,
whose rational characteristics satisfy p=1 mod10. Thus (9) contradicts
eventual purity directly. The genus-two surface input is an alternative.

For arbitrary characteristic-zero coefficients, adjoin the repeated roots
algebraically. A finitely generated Q-domain containing all phase
coefficients and r,s, with \(Brs(r-s)\) inverted in the double case or Br
inverted in the triple case, preserves the factorization and branch.
A maximal ideal has number-field residue field. The SAME single
specialization preserves every tail identity: each moment is a rational
polynomial in the finitely many phase coefficients. The number-field
contradictions apply, proving Theorem1. No PFP parameter theorem or
universal finite moment bound is used.

## 4. Simple sixfolds: full-End alternatives

For a geometrically simple J let \({\cal D}=\operatorname{End}^0(J)\).
This is a division algebra acting on the twelve-dimensional rational H_1.
Thus \(\dim_{\mathbb Q}{\cal D}\mid12\); the embedded E makes this dimension
a multiple of4, so it is4 or12.

Dimension4 is exactly E, the remaining case. Commutative dimension12 is a
degree12 field, hence potential CM. A noncommutative central division
algebra of dimension12 would be quaternion over a cubic center F0.
The commuting E and F0 generate a commutative DOMAIN inside the division
algebra, hence a field. Since E/Q is Galois of degree4, coprime to3, that
compositum has degree12. It fills the whole algebra and contradicts
noncommutativity. This field argument is specifically inside a DIVISION
algebra, not a claim about arbitrary commuting subfields of matrix algebras.
Every simple case other than full End=E is consequently potential CM.

## 5. Several blocks: all fourfold complements are covered

Project E into each geometric isotypic block of J. Each projection is a
unital field homomorphism and therefore injective. The block's H_1 has
dimension divisible by4, so its abelian dimension is even. Multiple blocks
have dimensions2+4 or2+2+2.

A surface block is potential CM and has a CM-type E-signature. Subtracting
its four possible CM types from (1,1,2,2) gives
\[
\begin{array}{c|c}
\text{surface signature}&\text{fourfold complement}\\ \hline
(1,1,0,0)&(0,0,2,2)\\
(1,0,1,0)&(0,1,1,2)\\
(0,1,0,1)&(1,0,2,1)\\
(0,0,1,1)&(1,1,1,1).
\end{array}                                                   \tag{10}
\]
Every conjugate pair (1,4),(2,3) of the surface sums to1. The original
polarization restricts to the blocks, and the E-adjoint relation remains
Rosati conjugation there.

The first complement is pure CM type: its rank-two E-H_1 splits into
rational E-line Hodge substructures, giving two potential-CM surfaces.
The next two are Galois relabelings/conjugates of (1,2,0,1), covered by
SLBF-FOURFOLD-ALL-SPECIAL-ORDINARY. Its verification includes the full-End,
potential-CM and nonsimple special cases; it is not a generic-End premise.
The final complement is balanced and uses
ED0-BALANCED-FOURFOLD-POSITIVE-ORDINARY.

Choose a common extension defining all factors and making every CM
surface ordinary on a density-one set. Reapply the balanced-fourfold
theorem over that common field if it occurs. Its positive-density ordinary
set intersects the CM surfaces' density-one set in positive density.
Products and isogenies preserve ordinarity away from finitely many primes.
We have not intersected arbitrary positive-density sets.

## 6. A sole nonsimple block

If \(J\sim B_0^m\) with B0 geometrically simple, \(\dim B_0\) divides6.
Dimensions1,2 are handled by the surface theorem (use an elliptic square
in dimension1); dimension6 was Section4. Only \(\dim B_0=3,m=2\) remains.

Let \(F_0=\operatorname{End}^0(B_0)\). Its division-algebra dimension divides6.
No divisor of6 has a noncommutative central-division square factor at least4,
so F0 is a field of degree \(f_0=1,2,3,6\). Degree6 is potential CM.
For \(f_0<6\), the E-action embeds E into M2(F0). If T is the image of
zeta5, its F0-minimal polynomial divides separable Phi5 and has degree at
most2 by Cayley--Hamilton. No primitive fifth root lies in F0, since
\(4\nmid f_0\). The minimal polynomial is neither linear nor split
quadratic, hence irreducible quadratic. Therefore F0[T] is an ACTUAL field
of rational degree2f0, containing E. Thus \(4\mid2f_0\), forcing f0=2 and
F0[T]=E. The scalar F0 is E's unique quadratic subfield Q(sqrt5), real.

A real quadratic field cannot act on an abelian threefold: H_1 has rank3
over it, so at a real embedding its complexified H_1 has dimension3.
Complex conjugation preserves that embedding and exchanges the (1,0)
and (0,1) parts, forcing this dimension to be even. Contradiction.
This excludes the non-CM threefold-square case. Sections4--6 prove
Theorem2.

## 7. Exact residual and structural obstructions

For a distinct cubic Q the selected Prym has signature(1,1,2,2).
Theorem2 implies: an eventual number-field phase (1) must have an
absolutely simple Prym with FULL geometric End^0=E. For arbitrary
characteristic-zero coefficients, a single specialization preserving the
distinct-root branch yields a number-field eventual phase, whose Prym
must have that remaining property. An arbitrary specialization is NOT
claimed to preserve any extra endomorphisms originally present over a
transcendental coefficient field.

The D=0 quaternion argument cannot extend unchanged. Put F=Eplus.
Each real F-embedding in Lie(J) has dimension1+2=3. A unital quaternion
F-subalgebra of End^0(J) containing the specified E would complexify to
M2(C) acting unital on that space. Every M2(C)-module has even dimension.
This is impossible. Also any invertible E-antilinear Hodge endomorphism
would identify conjugate E-Lie dimensions1 and2, again impossible.

Known symmetries do not yield a proper factor. On the anti factor iota=-1
and the nontrivial rational mu5 group-algebra block is the field E.
The algebra generated there by the deck action and iota is E, with no
nontrivial idempotents. The selected one-dimensional coherent-cohomology
character is not an abelian factor; a splitting needs extra endomorphisms
or a correspondence.

In the remaining full End=E case there is no imaginary quadratic subfield:
E's unique quadratic subfield is real. Fite's imaginary-quadratic
threefold/fourfold theorems do not cover this sixfold. The two
conjugate-pair minima in its signature sum to2, not1, so the exact
Cantoral Farfan et al. simple-signature gate also fails. This is a scope
obstruction, not a proof that every other method fails.

The smallest sufficient remaining arithmetic target is: for every actual
number-field smooth-cubic face here, after one fixed extension the selected
H_p in (8) is nonzero for infinitely many good absolute-degree-one primes.
Positive ordinary density for its simple full-End sixfold implies that
target, but the converse is not asserted. That target is NOT proved here.
Full even14/S24 and global primitive degree5 remain outside these results.

## 8. Existing input scopes and audit

The local proofs above establish normalization, selected scalar,
discriminant, division/matrix-algebra reductions and the no-quaternion
obstruction. Standard abelian-isogeny and weight-one Hodge foundations
remain imported.

The exact all-member balanced input is
ED0-BALANCED-FOURFOLD-POSITIVE-ORDINARY, in
EVEN_QUINTIC_D0_BALANCED_PRYM_ORDINARY_OBSTRUCTION.md, Sections1--3.
It assumes polarized quartic-CM multiplication, Rosati conjugation and
signature(1,1,1,1), with no simplicity or full-End hypothesis.
Its underlying exact ordinary input is Fite Theorem4; its explicit
quaternion field construction and Milne descent were already proved and
reviewed there. We reuse this theorem, not the old D=0 face theorem.

SLBF-FOURFOLD-ALL-SPECIAL-ORDINARY is in
SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md, Sections6--7.
Precomposing the E-action gives the relabelings in (10).
Its simple-signature ordinary input is used only after FULL geometric-End
equality has been checked in the appropriate fourfold case.

EXT-SLBF-SAWIN-SURFACE-ORDINARY and
EXT-SLBF-FITE-POTENTIAL-CM-ORDINARY each give density one after extension.
These are the precise inputs permitting the density intersection above.
Primary references:
[Sawin1506.04784v2](https://arxiv.org/abs/1506.04784v2),
[Fite2005.10185v3](https://arxiv.org/abs/2005.10185v3), and
[Cantoral Farfan et al.2508.11174v1](https://arxiv.org/abs/2508.11174v1).
