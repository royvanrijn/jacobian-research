# Beckmann--Black specialization with a Keller condition

## Status and main conclusion

This note separates a literal formulation, which is impossible for absolute
polynomial Keller maps, from the nonnormal point-field formulation that
matches the existing inverse-monodromy programme.

There are three conclusions.

1. If a `G-extension` means a `G`-Galois field extension \(L/\mathbb Q\)
   of degree \(|G|\), then asking for \(L\) itself to be a complete fiber of
   an absolute Keller map whose generic monodromy is the same group \(G\)
   is impossible for every nontrivial \(G\).  The generic action would be
   regular, hence the generic function-field extension would be Galois, and
   the Campbell--Razar--Wright Galois-case theorem would make the Keller map
   invertible.
2. The viable replacement fixes a core-free subgroup \(H<G\), uses the
   transitive action on \(G/H\), and asks for the degree-\([G:H]\) point
   field \(L^H\) as the complete fiber, with splitting closure \(L\).
3. In this corrected sense the **smooth-affine chart version is already
   solved** for
   \[
     A_4,\qquad D_5,\qquad F_{20},\qquad A_5.
   \]
   Each group has a two-parameter generic polynomial.  A derivative-unit
   suspension turns its root cover into one fixed determinant-one morphism
   of smooth affine threefolds, without changing its monodromy or any
   complete root fiber.  What remains open is the stronger absolute output:
   a polynomial Keller self-map of affine space with the same property.

Thus proper-monodromy Kellerization really would create a uniform
realization machine, but only after `G-extension` is interpreted through a
nonnormal transitive action.  The regular-action reading is obstructed, not
merely unfinished.

## 1. Three classical quantifiers

Let \(L/\mathbb Q\) be a Galois extension with group \(G\).

The classical Beckmann--Black problem asks whether there is a
\(\mathbb Q\)-regular Galois extension

\[
 E_L/\mathbb Q(T),\qquad \operatorname{Gal}(E_L/\mathbb Q(T))=G,
\]

whose specialization at an unramified rational point is \(L/\mathbb Q\).
The cover may depend on \(L\).

A one-parameter \(G\)-parametric extension asks for one fixed
\(\mathbb Q\)-regular \(G\)-extension \(E/\mathbb Q(T)\) which specializes
to every \(G\)-extension of \(\mathbb Q\).

An \(r\)-parameter generic polynomial is uniform in a different and stronger
base-field sense: one polynomial over
\(\mathbb Q(u_1,\ldots,u_r)\) specializes to every \(G\)-extension after
every allowed scalar extension.  It has an \(r\)-dimensional parameter base,
so it is not the same assertion as a one-parameter parametric extension.

These distinctions follow the terminology in
[Dèbes' Beckmann--Black paper](https://www.numdam.org/item/ASNSP_1999_4_28_2_273_0.pdf)
and
[Legrand's parametric-extension paper](https://arxiv.org/abs/1310.6682).
Legrand's later
[nonnormal variant](https://arxiv.org/abs/2111.07155)
is especially relevant here: it shows that removing normality changes the
problem substantially.

For the four pilot groups, the classical arithmetic input is unusually
favorable.  The Beckmann--Black property is known for alternating groups and
for the relevant dihedral case.  Black's semidirect-product theorem also
applies to

\[
 F_{20}=C_5\rtimes C_4
\]

because the two factors have coprime orders and the complement \(C_4\) has
the abelian arithmetic lifting property.  See
[Black, *On Semidirect Products and the Arithmetic Lifting Property*](https://doi.org/10.1112/S002461079900784X).
The Keller obstruction below is therefore geometric, not a lack of regular
Galois realizations.

## 2. The literal absolute question has a negative answer

Call a polynomial map

\[
 F:\mathbb A^m_{\mathbb Q}\longrightarrow\mathbb A^m_{\mathbb Q}
\]

Keller if \(\det DF\in\mathbb Q^\times\).  Its geometric degree is the
degree of

\[
 \mathbb Q(x_1,\ldots,x_m)/
 \mathbb Q(F_1,\ldots,F_m).
\]

### Proposition 2.1 -- regular-action obstruction

Let \(G\ne1\).  There is no absolute polynomial Keller map \(F\) satisfying
both:

1. the generic inverse monodromy is \(G\); and
2. some connected complete fiber is a \(G\)-Galois extension
   \(L/\mathbb Q\), viewed as a finite etale algebra.

Here “complete” means that the fiber length equals the geometric degree.

### Proof

The complete fiber has length

\[
 [L:\mathbb Q]=|G|,
\]

so \(\operatorname{gdeg}(F)=|G|\).  The transitive generic action of the
group \(G\) is therefore an action of a group of order \(|G|\) on
\(|G|\) sheets.  Its point stabilizer is trivial, so the action is regular.
Equivalently, the source function field is already its Galois closure over
the target function field.

The Galois case of the Jacobian theorem says that a characteristic-zero
Keller map with normal function-field extension has a polynomial inverse.
Consequently its geometric degree is one, contradicting \(|G|>1\).
\(\square\)

The classical sources for the last step are:

- [Campbell, *A condition for a polynomial map to be invertible*](https://doi.org/10.1007/BF01349234);
- [Razar, *Polynomial maps with constant Jacobian*](https://doi.org/10.1007/BF02764906);
- [Wright, *On the Jacobian conjecture*](https://doi.org/10.1215/ijm/1256047158).

This obstruction concerns the simultaneous requirements “complete Galois
fiber” and “same generic group.”  The repository's
[finite-etale fiber theorem](verified/FINITE_ETALE_KELLER_FIBERS.md)
can still put the field \(L\) itself into a complete fiber, but the ambient
generic group is then generally a larger nonnormal permutation group.

## 3. The corrected point-field problem

Fix a core-free subgroup \(H<G\), and put

\[
 n=[G:H].
\]

For a \(G\)-Galois extension \(L/\mathbb Q\), the point field
\(K=L^H\) has degree \(n\), and its Galois closure is \(L\).  Core-freeness
is exactly what makes the action \(G\curvearrowright G/H\) faithful.

### Definition 3.1 -- \((G,H)\)-Keller-parametric map

An absolute \((G,H)\)-Keller-parametric map is a fixed polynomial Keller map

\[
 F:\mathbb A^m_{\mathbb Q}\longrightarrow\mathbb A^m_{\mathbb Q}
\]

such that:

1. its geometric and arithmetic generic monodromy are both \(G\) in the
   action on \(G/H\);
2. for every \(G\)-Galois extension \(L/\mathbb Q\), there is a rational
   regular target \(y\) with
   \[
     F^{-1}(y)\simeq\operatorname{Spec}(L^H);
   \]
3. the displayed fiber is complete, hence has length \(n\), and its
   splitting closure is \(L\).

The first pilot actions are:

| group | point stabilizer \(H\) | degree \([G:H]\) |
|---|---:|---:|
| \(A_4\) | \(C_3\) | \(4\) |
| \(D_5\) | \(C_2\) | \(5\) |
| \(F_{20}=C_5\rtimes C_4\) | \(C_4\) | \(5\) |
| \(A_5\) | \(A_4\) | \(5\) |

In each row \(H\) is core-free and the conjugacy class of point stabilizers
is the natural one.  None of these four actions is regular, so Proposition
2.1 does not apply.

A **chart \((G,H)\)-Keller-parametric map** has the same properties but is a
determinant-one morphism between smooth affine varieties equipped with
nowhere-vanishing algebraic volume forms.  It need not be a polynomial
self-map of affine space.  This is the output proved next.

## 4. Generic polynomials give universal Keller charts

Let

\[
 P(\mathbf u,T)\in\mathbb Q(\mathbf u)[T],
 \qquad \mathbf u=(u_1,\ldots,u_r),
\]

be a monic generic polynomial for \(G\) in a fixed transitive degree-\(n\)
action.  Shrink the rational parameter space to a smooth affine open
\(B\subset\mathbb A^r\) on which:

- all coefficients of \(P\) are regular;
- its discriminant is a unit; and
- its generic splitting field is regular with geometric and arithmetic
  group \(G\).

Put

\[
 X=\operatorname{Spec}_B
   \frac{\mathcal O_B[T]}{(P(\mathbf u,T))}.
\]

Then \(X\to B\) is finite etale of degree \(n\), and
\(P_T=\partial P/\partial T\) is a unit on \(X\).

### Theorem 4.1 -- derivative-unit chart compiler

The map

\[
\begin{aligned}
 \widehat\pi:X\times\mathbb A^1_z
   &\longrightarrow B\times\mathbb A^1_Z,\\
 (\mathbf u,T,z)&\longmapsto
 \left(\mathbf u,\frac{z}{P_T(\mathbf u,T)}\right)
\end{aligned}
\tag{4.1}
\]

is one fixed finite etale degree-\(n\) morphism with the same geometric and
arithmetic generic monodromy \(G\).  With the hypersurface residue volume
form on \(X\), it has Jacobian one.  At every rational \(b\in B(\mathbb Q)\)
and \(Z_0\in\mathbb Q\), its complete fiber is

\[
 \widehat\pi^{-1}(b,Z_0)
 \simeq
 \operatorname{Spec}\frac{\mathbb Q[T]}{(P(b,T))},
 \qquad z=P_T(b,T)Z_0.
\tag{4.2}
\]

If \(P\) is generic for the action \(G/H\), every \(G\)-Galois extension
\(L/\mathbb Q\) therefore occurs through the complete point field \(L^H\).

### Proof

Because \(P_T\) is a unit in the root algebra, (4.1) is regular.  The change
of fiber coordinate

\[
 z=P_TZ
\]

identifies it, as a finite cover, with the base change of \(X\to B\) along
\(B\times\mathbb A^1_Z\to B\).  Its degree, etaleness, and generic monodromy
are therefore unchanged.

If \(\omega_B\) is a coordinate volume form on \(B\), choose the residue
form \(\omega_X\) so that

\[
 \pi^*\omega_B=P_T\omega_X.
\]

Since terms involving \(dP_T\) wedge to zero against the top form on \(X\),

\[
 \widehat\pi^*(\omega_B\wedge dZ)
 =
 P_T\omega_X\wedge
 \left(\frac{dz}{P_T}-\frac{z\,dP_T}{P_T^2}\right)
 =
 \omega_X\wedge dz.
\]

Finally, imposing \(Z=Z_0\) eliminates \(z\) by
\(z=P_TZ_0\), proving the scheme-theoretic fiber identity (4.2).
Genericity supplies, for every \(L\), a squarefree irreducible
specialization whose splitting field is \(L\) in the prescribed action.
\(\square\)

The result is a uniform cover-preserving version of the derivative-unit
suspension already used in the
[oriented quartic checkpoint](extended-geometry/A4_KELLER_INVERSE_COVER.md).
It must not be confused with the
[universal relative quadratic-gauge map](verified/UNIVERSAL_RELATIVE_KELLER_MAP.md):
that coefficient atlas realizes every presented fiber, but its moving
inverse pencil has generic group \(S_n\).  Theorem 4.1 preserves the proper
group \(G\) because it adds only a reconstructible linear coordinate.

## 5. The four first cases

Jensen--Ledet--Yui's
[*Generic Polynomials*](https://library.slmath.org/books/Book45/files/book45.pdf)
gives the required arithmetic inputs:

| group | generic-polynomial source | parameters | chart conclusion |
|---|---|---:|---|
| \(A_4\) | Theorem 2.2.9 | \(2\) | fixed determinant-one smooth affine threefold cover |
| \(D_5\) | Brumer, Theorem 2.3.5 | \(2\) | same |
| \(F_{20}\) | Lecacheux, Theorem 2.3.6 | \(2\) | same |
| \(A_5\) | Buhler construction, Theorem 2.3.7 | \(2\) | same |

The book proves that all transitive subgroups of \(S_5\) have
two-parameter generic polynomials over \(\mathbb Q\).  In particular, the
generic-polynomial input parametrizes **all** \(G\)-extensions, not merely
an infinite Hilbert family.

### 5.1 \(A_4\)

There are two chart proofs.

First, Theorem 4.1 applies to the two-parameter generic polynomial.  Second,
the repository's oriented depressed-quartic chart is directly universal.
Given an \(A_4\)-Galois extension \(L/\mathbb Q\), let

\[
 K=L^{C_3}.
\]

Choose a primitive element of the quartic field \(K\), make its minimal
polynomial monic, and translate it to

\[
 T^4+pT^2+qT+r.
\]

Its discriminant is a nonzero rational square because its splitting group
lies in \(A_4\).  Choosing \(D\in\mathbb Q^\times\) with
\(D^2=\Delta(p,q,r)\) gives a rational point of the oriented base.  Formula
(5.4) of the oriented-quartic note then gives the complete fiber
\(\operatorname{Spec}K\), whose splitting closure is \(L\).

Thus the existing determinant-one boundary-complement map is already
\((A_4,C_3)\)-parametric at chart level.  Its affine follow-up produces a
literal polynomial \(\mathbb A^3\)-map with exact \(A_4\) monodromy, but
its determinant is

\[
 4W^2K^3L,
\]

so the absolute Keller step remains open.

### 5.2 \(D_5\) and \(F_{20}\)

Brumer's and Lecacheux's generic polynomials, followed by Theorem 4.1, solve
the chart problem for every \(D_5\)- and \(F_{20}\)-extension.

The repository's
[fixed quintic Galois stratification](FIXED_QUINTIC_GALOIS_STRATIFICATION.md)
proves something complementary: one fixed absolute Keller map has
infinite rational \(D_5\) and \(F_{20}\) specialization families.  Its
ambient generic group is nevertheless \(S_5\), and the displayed Brumer
curve and De Moivre surface are not proved to parametrize all extensions.
It therefore does not solve Definition 3.1.

For an absolute construction, the natural next inputs are the actual Brumer
and Lecacheux generic covers, not the smaller non-universal subfamilies
already embedded in the \(S_5\) map.

### 5.3 \(A_5\)

The Buhler two-parameter polynomial and Theorem 4.1 solve the chart problem
for every \(A_5\)-extension.  The fixed quintic Keller map currently has
exact \(A_5\) complete fibers in both real signatures and all four
unramified cycle types.  It also has two regular Mestre source pencils.
The direct affine lift of the totally real pencil reduces to a genus-eight
curve, and a higher non-affine descent surface remains unparametrized.
Thus infinitude inside the fixed quintic map is still open.  In any case its
ambient generic group is \(S_5\), so even a successful parametrization there
would be a specialization-family theorem rather than an absolute
\((A_5,A_4)\)-Keller-parametric map.

## 6. Exact frontier

The four levels are now:

| level | \(A_4\) | \(D_5\) | \(F_{20}\) | \(A_5\) |
|---|---|---|---|---|
| classical Beckmann--Black | known | known | known | known |
| two-parameter generic polynomial | known | known | known | known |
| fixed determinant-one smooth affine chart preserving all point fields | proved | proved | proved | proved |
| absolute polynomial Keller map with generic group \(G\) and all point fields | open | open | open | open |

The remaining question is therefore:

> For a core-free \(H<G\), can the generic root cover for \(G/H\), together
> with its derivative-unit suspension, be completed to a polynomial Keller
> self-map of affine space without changing its function-field extension or
> losing any rational complete point fibers?

The phrase “all point fibers” matters.  A completion could preserve the
generic extension while deleting some rational specializations at its
boundary.  A valid construction must prove both cover preservation on a
dense open and sheet-by-sheet fullness at every parameter used by
genericity.

## 7. Recommended attack order

1. **Finish \(A_4\).**  It already has a polynomial affine monodromy core,
   an exact determinant ledger, and several increasingly sharp
   affine-modification obstructions.  It is the nearest absolute case.
2. **Kellerize Brumer's \(D_5\) cover.**  Work with the full two-parameter
   generic polynomial.  Compute \(1/P_T\bmod P\), factor its pole ledger,
   and seek an affine modification that preserves the degree-five
   nonnormal cover.
3. **Repeat for Lecacheux's \(F_{20}\) cover.**  Its nonsquare
   discriminant character distinguishes it from the alternating cases and
   tests whether the affine completion can retain a nontrivial sign
   quotient.
4. **Treat \(A_5\).**  The generic polynomial is available, but the
   non-solvable monodromy and its rational denominator boundary make it the
   strongest degree-five test.

For each case, the verifier should separately certify:

- geometric and arithmetic generic monodromy;
- the quotient-ring identity for \(1/P_T\);
- the determinant divisor before and after every modification;
- preservation of the generic function-field extension;
- scheme-theoretic completeness of specialized fibers; and
- arithmetic surjectivity onto all \(G\)-extensions, imported from the
  generic-polynomial theorem rather than inferred from Hilbert
  irreducibility.

No bounded search or infinite Hilbert family can replace the last item.
