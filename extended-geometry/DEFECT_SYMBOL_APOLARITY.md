# Defect symbols, inverse systems, and finite-jet closure

This note extracts a common local-algebra package from the degree-forty-two
normal cone, the degree-thirty braid cells, primitive merger blocks, and the
cubic-normalization multiplication tensor.  Its main purpose is to separate
three objects which are closely related but not identical:

1. the first nonzero Kuranishi equations;
2. the full tangent-cone algebra;
3. a structural potential, such as the ternary cubic in the triple-cover
   multiplication law.

Keeping these objects separate makes the apolarity statements precise and
shows exactly when a leading symbol replaces all higher local equations.

Throughout, `k` has characteristic zero.  Divided powers give the
characteristic-free version of the inverse-system statements.

## 1. The intrinsic symbol package

Let `Z -> S` be a projected intersection problem and let `p` be a point of a
fiber.  After removing smooth variables, choose a minimal formal presentation

\[
 \widehat{\mathcal O}_{Z,p}\simeq k[[E]]/I,
 \qquad I\subset (E)^2.
\tag{1}
\]

Here `E` is the surviving cotangent space, intrinsically
`H_0(L_(Z/S) tensor k(p))`.  A minimal relation space is supplied by
`H_1(L_(Z/S) tensor k(p))`.  If `d` is the least order of a member of `I`,
write

\[
 W_d=\operatorname{in}_d(I)
 \subseteq \operatorname{Sym}^d(E).
\tag{2}
\]

Changing formal coordinates acts through `GL(E)` on (2), while changing
minimal relation generators does not change the subspace `W_d`.  Thus the
first symbol is naturally a point of

\[
 \operatorname{Gr}\bigl(\dim W_d,\operatorname{Sym}^d(E)\bigr)/GL(E),
\tag{3}
\]

not an ordered tuple of forms.

There are two associated graded algebras:

\[
 \mathcal S_p^{(d)}
   =\operatorname{Sym}(E)/(W_d),
 \qquad
 \mathcal D_p
   =\operatorname{gr}_{\mathfrak m}\widehat{\mathcal O}_{Z,p}
   =\operatorname{Sym}(E)/\operatorname{in}(I).
\tag{4}
\]

Call the first the **leading defect-symbol envelope** and the second the
**full defect-symbol algebra**.  There is a canonical surjection

\[
 \mathcal S_p^{(d)}\longrightarrow\mathcal D_p.
\tag{5}
\]

The distinction is useful even when only nilpotence is needed.  A finite
leading envelope gives an upper bound for the Hilbert function and Loewy
length of the full algebra without requiring equality in (5).

In a family, `dim E`, `d`, and `dim W_d` need not be constant.  The symbol is
therefore a Grassmannian-bundle section only after stratifying the base by
these ranks.  Rank drops and changes of the first nonzero degree are part of
the invariant, not failures of the construction.

## 2. Complete-intersection closure lemma

Suppose the initial ideal contains a homogeneous regular sequence

\[
 f_1,\ldots,f_r,\qquad
 \deg f_i=d_i,\qquad r=\dim E.
\tag{6}
\]

Then every quotient of the complete intersection

\[
 C=\operatorname{Sym}(E)/(f_1,\ldots,f_r)
\]

has maximal ideal nilpotent in degree

\[
 1+\sigma,\qquad
 \sigma=\sum_{i=1}^r(d_i-1).
\tag{7}
\]

Indeed, the Hilbert series of `C` is

\[
 \prod_{i=1}^r(1+T+\cdots+T^{d_i-1}),
\tag{8}
\]

whose top degree is `sigma`.  Consequently, if a formal defect class is
already known to lie in

\[
 I+\mathfrak m^{\sigma+1},
\tag{9}
\]

then it vanishes exactly in the completed local quotient.

This single lemma explains the relevant closure exponents:

| generator degrees | Hilbert vector | socle degree | vanishing power |
|---|---|---:|---:|
| `(2,2)` | `(1,2,1)` | 2 | `m^3` |
| `(3,3)` | `(1,2,3,2,1)` | 4 | `m^5` |
| `(5)` | `(1,1,1,1,1)` | 4 | `m^5` |
| `(4,2)` | `(1,2,2,2,1)` | 4 | `m^5` |

Thus invariant theory need not reconstruct the whole completed ideal in
order to close a finite-jet argument.

## 3. Macaulay duality and its scope

Let `Gamma(E dual)` be the divided-power algebra, with `Sym(E)` acting by
contraction.  If `A=Sym(E)/J` is Artin Gorenstein, then

\[
 J=\operatorname{Ann}(F)
\tag{10}
\]

for a divided-power polynomial `F`, unique up to a nonzero scalar.  If the
socle has dimension `t>1`, the inverse system is instead a `t`-generated
`Sym(E)`-submodule.  A single potential is therefore available precisely on
the Gorenstein locus.

This also distinguishes two uses of the word *pencil*:

- a pencil of equations is a two-plane in `Sym^d(E)` and, when it is a
  codimension-two complete intersection, has one Macaulay-dual form;
- a pencil of dual forms is a two-plane in `Gamma(E dual)` and defines a
  type-two level algebra, with a pencil of Gorenstein quotients.

The two constructions are contravariant, not identical.

## 4. Binary complete intersections

Write

\[
\begin{aligned}
 q&=a x^2+bxy+cy^2,\\
 r&=d x^2+exy+fy^2.
\end{aligned}
\tag{11}
\]

The divided-power quadratic dual to `(q,r)` is

\[
 F=(bf-ce)X^{[2]}+(cd-af)X^{[1]}Y^{[1]}
     +(ae-bd)Y^{[2]}.
\tag{12}
\]

The coefficient vector in (12) is the cross product of the two coefficient
rows in (11), so both `q` and `r` annihilate `F`.  Moreover,

\[
 \operatorname{disc}(F)
 =(cd-af)^2-(bf-ce)(ae-bd)
 =\operatorname{Res}(q,r).
\tag{13}
\]

Hence the resultant is exactly the discriminant of the dual conic, not just
an analogous invariant.

For the degree-forty-two quadrics, after removing their common factor `w_1`,
(12) gives, up to scale,

\[
\begin{aligned}
 A&=-\frac98(-e_1^2+e_2t),\\
 B&=-\frac9{16}(-e_1e_2+t),\\
 C&=\frac98e_1.
\end{aligned}
\tag{14}
\]

Therefore

\[
 B^2-AC
 =\frac{81}{256}\bigl((t+e_1e_2)^2-4e_1^3\bigr).
\tag{15}
\]

The divisor `Delta=0` is precisely the degeneracy divisor of the Macaulay
dual quadratic.  On the other branch, where the quadratic symbol vanishes,
the regular pair of binary cubics has a binary quartic inverse system.

Write that divided-power quartic in the normalized form

\[
 F_0X^{[4]}+F_1X^{[3]}Y+F_2X^{[2]}Y^{[2]}
   +F_3XY^{[3]}+F_4Y^{[4]}
\]

and use the classical invariants

\[
\begin{aligned}
 I_4&=F_0F_4-4F_1F_3+3F_2^2,\\
 J_4&=F_0F_2F_4+2F_1F_2F_3-F_0F_3^2-F_1^2F_4-F_2^3.
\end{aligned}
\tag{16}
\]

The exact terminal cubics give

\[
\begin{aligned}
 I_4&=\frac{3\cdot5^8}{2^{24}}A(A+3t^2)B,\\
 J_4&=-\frac{5^{12}}{2^{36}}A^2B^2,\\
 I_4^3-27J_4^2
   &=\frac{3^6\,5^{24}}{2^{72}}t^6A^3B^3.
\end{aligned}
\tag{17}
\]

Equivalently, the base polynomials satisfy the compact identity

\[
 (A+3t^2)^3-AB=27t^6.
\tag{18}
\]

This separates two phenomena which the cubic resultant alone does not:

- `AB=0` is where the two cubic equations fail to be a regular sequence and
  the cofactor construction of the quartic inverse system drops rank;
- on `D(AB)`, the quotient remains a `(3,3)` Artin Gorenstein complete
  intersection, but its dual quartic becomes singular exactly on `t=0`.

Thus `t=0` is an invariant-theoretic stratum inside the Gorenstein locus,
not another failure of finite-jet closure.

## 5. Known inverse systems

The primitive and braid-cell algebras have monomial inverse systems:

| algebra | inverse-system generator |
|---|---|
| `k[epsilon]/(epsilon^2)` | `X` |
| `k[epsilon_1,...,epsilon_t]/(epsilon_i^2)` | `X_1...X_t` |
| `k[u]/(u^5)` | `U^[4]` |
| `k[u,v]/(u^2,v^2)` | `U V` |
| `k[u,v]/(u^4,v^2)` | `U^[3] V` |

The Boolean merger algebra is therefore the apolar algebra of a squarefree
monomial.  Independence of Hensel clusters appears dually as the presence
of every variable in that monomial.

Two known boundary algebras show why the package must also retain
noncyclic inverse systems:

\[
 k[\ell,s]/(\ell^2,\ell s,s^3)
 \quad\longleftrightarrow\quad
 \langle L,S^{[2]}\rangle,
\tag{19}
\]

and

\[
 k[x,y]/(x^3,xy,y^2)
 \quad\longleftrightarrow\quad
 \langle X^{[2]},Y\rangle.
\tag{20}
\]

Both have socle dimension two.  They lie naturally on a boundary of the
Gorenstein classification, but neither is represented by one potential.

## 6. The ternary cubic is a structural potential

At a reduced cubic-normalization defect, the actual fiber algebra is the
square-zero length-four algebra `k plus k^3`, whose socle has dimension
three.  Its Macaulay inverse system is generated in degree one.

The ternary cubic proved in the cubic frontend instead describes the first
nonzero generalized Miranda--Tan multiplication tensor

\[
 \phi:\operatorname{Sym}^3(M)\longrightarrow\det(M).
\tag{21}
\]

It is therefore a **structural potential**, not the Macaulay dual of the
defect fiber.  One may form the auxiliary apolar algebra

\[
 \operatorname{Ap}(h)
 =\operatorname{Sym}(E)/\operatorname{Ann}(h),
\tag{22}
\]

but (22) should not be identified with the original square-zero fiber.

The complete local package should consequently be written

\[
 \boxed{(\mathcal D_p,\tau_p)}
\tag{23}
\]

where `mathcal D_p` is the full defect-symbol algebra and `tau_p` is the
first additional structural tensor.  For Gorenstein complete intersections
the two entries are linked by apolarity.  For the cubic frontend they carry
different information.

The coarse cubic division into smooth, nodal, cuspidal, and reducible loci
should be refined as follows:

- smooth cubics retain the `j`-parameter;
- singular irreducible cubics split into nodal and cuspidal types;
- reducible cubics split into conic-plus-line and three-line configurations,
  with transverse, tangent, concurrent, and nonreduced subtypes.

Aronhold invariants, the discriminant, and factorization type provide a
finite first pass before higher-order lifting and the Keller-specific gates
are imposed.

## 7. Computation programme

The next computations should be performed in this order.

1. Extract `E`, the minimal relation space, and `W_d` after unit-pivot
   elimination, and record their ranks on every base stratum.
2. Compute the leading envelope and full tangent-cone algebra separately.
3. On Gorenstein strata, solve linearly for the inverse-system generator and
   calculate catalecticants and classical invariants.
4. Refine the degree-forty-two quartic strata beyond (17), using root
   partitions and catalecticant rank on `t=0`, `A=0`, and `B=0`.
5. On nongorenstein strata, retain the complete inverse-system module and
   its generator degrees.
6. For the cubic frontend, stratify the multiplication potential by cubic
   invariants and run the higher-lift calculation separately on each
   stratum.

The exact elementary identities and the inverse-system table in this note
are checked by
[`verify_defect_symbol_apolarity.py`](../scripts/verify_defect_symbol_apolarity.py).

## 8. External coordinates

The inverse-system conventions and the distinction between cyclic
Gorenstein systems and finitely generated inverse-system modules follow
Elias--Rossi,
[The structure of the inverse system of Gorenstein k-algebras](https://arxiv.org/abs/1705.05686).
The contravariant pencil construction is the one in Iarrobino,
[Hilbert functions of Gorenstein algebras associated to a pencil of forms](https://arxiv.org/abs/math/0412361):
a two-plane of dual forms defines a type-two level algebra, and its lines
define Gorenstein quotients.  For socle degree three, the reduction of
isomorphism classes to projective cubic hypersurfaces is proved by
Elias--Rossi,
[Isomorphism classes of short Gorenstein local rings via Macaulay's inverse
system](https://arxiv.org/abs/0911.3565).

These references support the inverse-system classification layer.  The
identification of the repository's degree-forty-two resultant with the
dual-quadratic discriminant is the direct calculation (11)--(15), not an
external input.
