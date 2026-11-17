# All quadratic Jacobian-power directions have explicit projective charts

The unmatched quadratic frontier is not an arbitrary coefficient search.
A two-sign midpoint argument forces a projective affine-chart isomorphism.
In rank4 every determinant-compatible primitive quadratic direction has
one of four explicit forms, with an arbitrary invertible output matrix
retained. Its source and target boundary hyperplanes need not agree.

The global multiplier-residue homotopy then bounds every nonboundary
irreducible scalar factor by degree3. These are structural reductions,
not a realization or an absolute exceptional-action obstruction.

## 1. A quadratic power Jacobian determines an affine target chart

Work over algebraically closed characteristic-zero constants.
Let V:k^r->k^r be homogeneous quadratic, r>=2, and let u be a nonzero
source linear form, with
\[
                            \det DV=\delta u^r,\quad\delta\ne0.
\]

**KMFI142.** The projective map [V] restricts to an isomorphism
\[
              \{u\ne0\}\ \simeq\ {\mathbb P}^{r-1}\setminus H_\alpha,
                                                               \tag{1}
\]
where H_alpha is a unique target hyperplane.
After scaling its equation, alpha(V)=u^2.
In the source and target affine charts, this is a polynomial
automorphism of degree at most2.
Neither alpha proportional to u nor a matched component is asserted.

There are no base points on u!=0: Euler gives DV(p)p=2V(p), so V(p)0
there would contradict invertibility of DV.
The induced projective differential is invertible, since the
invertible affine differential takes the radial direction to a nonzero
radial direction. Thus the restriction is etale.

The quadratic midpoint identity is
\[
                  V(p)-V(q)=DV((p+q)/2)(p-q).
\]
If [V(p)]=[V(q)] with both u coordinates nonzero, scale q by a square
root so V(p)=V(q). The same equality holds with -q.
At least one of u(p)+u(q),u(p)-u(q) is nonzero.
For that sign the midpoint differential is invertible, forcing p=q
or p=-q. The projective points are equal.
The proof works after every algebraically closed field extension,
so the etale morphism is universally injective and hence an open
immersion by [Stacks Tag025G](https://stacks.math.columbia.edu/tag/025G).
The retained [source manifest](sources/unmatched-quadratic-chart-manifest.json)
pins this exact imported result.

Here is why its image U, an affine space, is the complement of one
hyperplane. If the complement had no divisors, every regular function
on U would have no codimension-one pole in projective space and hence
be constant, impossible for affine space. If it had two distinct
divisors D1,D2, the ratio of defining equations raised to equal degrees
would give a nonconstant unit on U, again impossible.
Thus it has one divisor D.
A general hyperplane cut restricts to a principal divisor on U, since
affine space is factorial. Extending its defining rational function
to projective space and taking degrees gives 1+a*degD=0 for an integer a.
Consequently degD1. Additional codimension>=2 omissions outside D
are impossible: on the resulting affine chart every rational function
regular on U is polynomial, by cancellation of any denominator divisor.
Since U itself is affine, its inclusion is an isomorphism.
This proves the boundary assertion.

The homogeneous polynomial alpha(V) is nonzero and has no zero where
u!=0. Its only possible irreducible factor is u, so alpha(V)=constant*u^2.
The corresponding chart ratios are polynomials of degree<=2, with a
polynomial inverse by(1).

The affine direction map V itself is not injective: V(p)=V(-p).
In fact k(p)/k(V) has degree2. The projective inverse recovers p/u
rationally and alpha(V)=u^2 recovers the remaining radial quadratic.
This direction-map cover is not the Keller sweep field and is not
a new Keller occurrence action.

## 2. Four exhaustive primitive quadratic forms in rank4

Choose source coordinates p=(u,x,y,z).

**KMFI143.** Up to a matching constant source/target conjugacy
preserving the chosen u, every primitive rank4 quadratic direction
satisfying detDV=delta*u^4 has the form
\[
                              V=A K_\nu,\quad A\in GL_4(k),
\]
where
\[
 K_A=(u^2,ux,uy+x^2,uz+y^2),\qquad
 K_B=(u^2,ux,uy+x^2,uz+xy),
\]
\[
 K_C=(u^2,ux,uy,uz+x^2),\qquad
 K_D=(u^2,ux,uy,uz+xy).                                  \tag{2}
\]
All invertible A are allowed at this determinant-compatible stage.
Conversely every displayed A,K has primitive coordinates and
\[
                              \det DV=2\det(A)u^4.
\]
The nonprimitive affine-chart case is V=u*A*p.

Indeed KMFI142 first makes the direction into (u^2,u^2*v(p/u))
by an independent target chart choice.
Write v=M(x+H)+b with H homogeneous quadratic.
The constant Jacobian makes DH nilpotent.
The source-pinned de Bondt Theorem3.2 triangularizes it in dimension3,
as in KMFI133. Quadratic image dimension2 gives H_A or H_B there.
Image dimension1 gives (0,0,Q(x,y)); over algebraically closed constants
a nonzero binary quadratic is a square or a product of independent
linear forms, giving K_C or K_D.
The affine matrix M, constants b and the independent target chart
choice combine into the unrestricted constant matrix A.
They are not discarded by a claimed triangularization of the whole V.

The projective base loci of these forms are:
K_A has the point u=x=y=0;
K_B and K_C have the line u=x=0;
K_D has the two lines u=x=0 and u=y=0.
In particular every base locus is a finite union of linear subsets
inside u0, and each is contained in a finite union of lines in that plane.

The target boundary equation is alpha=e0^t*A^-1, while the source
boundary equation is u=e0^t. They agree exactly when A has first row
(a,0,0,0). After scalar normalization, precisely that subcase is
already closed by KMFI141.
An arbitrary A does not preserve the matched hypothesis.

For a normalized full-rank scalar F of degree s, the remaining exact
KMFI125 gate is h=lambda*u, with
\[
                              \lambda^4=(s+2)\det A.       \tag{3}
\]
The unique feedback form, nonnegative integer logarithmic residues
and full-rank coframe are still required. Formula(2) is not an assertion
that every A admits such a scalar or a desired monodromy group.

## 3. Global scalar-factor homotopy and weighted preservation

Let H=FV be any homogeneous normalized full-rank input, with V any
positive homogeneous degree and primitive as in KMFI125.
Write A0=DV, eta=dlogF and f=-Fh.

**KMFI144.** Every irreducible homogeneous factor D of F satisfies
\[
                            D(p+tV(p))\in(D)\,k[p,t].     \tag{4}
\]
No pure-power component of V is required.
This follows directly from
\[
 \eta(TI-A0)^{-1}V=1-(T-h)^r/\det(TI-A0).
\]
Its separate divisor residue makes dD*A0^j*V divisible by D for all j;
the first r powers suffice. The formal inverse proof of KMFI139,
now applied to the vector field V on the full source affine space,
proves(4). For homogeneous D it contains the projective line
joining [p] and [V(p)] whenever those points differ.

For a quadratic direction of KMFI142, each factor D not proportional
to u moreover satisfies
\[
                            D(V)=c_D u^{\deg D}D,\quad c_D\ne0. \tag{5}
\]
If alpha is not proportional to u, alpha(p) cannot be a factor of F.

To prove this, (4)'s top coefficient gives D divides D(V).
The source divisor D0 meets U={u!=0}; it is not contracted because
[V] is an isomorphism there. It cannot equal the target boundary
H_alpha, whose inverse image in U is empty.
Its image is a closed irreducible divisor of U' contained in D0,
hence the whole D0 intersect U'.
Thus inside U the pullback of D has exactly D0, with multiplicity1
by etaleness. All other homogeneous factors of D(V) are supported on u0.
Comparing degrees 2degD=degD+k gives k=degD, proving(5).
The argument does not give(5) for the boundary factor u.

## 4. Individual scalar factors now have degree at most three

**KMFI145.** For every normalized full-rank rank4 input with primitive
homogeneous quadratic direction V, choose u with detDV=delta*u^4.
Every irreducible factor D of F not proportional to u has degree<=3.
Total scalar degree and multiplicities remain unbounded.
The factor u is not excluded.

Suppose degD=m>3. Set u1 and q=D(1,x,y,z), an irreducible polynomial
of the same degree m.
The projective vector field
\[
                       W=V_{\rm high}(1,x)-V_u(1,x)x
\]
has degree at most3. It cannot vanish generically on q0:
degree comparison would make it vanish identically, which makes
V a scalar linear form times p, contrary to primitivity.

By KMFI144 the surface D0 is covered by the selected projective lines
joining p,V(p). The selected-line-family argument of KMFI140 applies
to any nonplane projective surface, without requiring the same source
and target affine charts. It makes [V] preserve each general selected
projective line. The restriction is generically injective on its
u!=0 part by KMFI142, so it is a degree1 self-map of P1.

The four quadratic sections restricted to that line therefore have
one common base point: their degree2 drops to degree1.
There are no base points off u0. Thus the point where the selected
line meets u0 lies in the base locus in(2).
Its homogeneous direction is
\[
                       R(p)=uV(p)-V_u(p)p.
\]
Each base locus is in a finite union of lines in u0.
Irreducibility forces one independent linear form ell, not proportional
to u, to satisfy ell(R)=0 on D0.
Its degree is at most3<m, so this is a global polynomial identity:
\[
                            u\,\ell(V)=V_u\,\ell(p).      \tag{6}
\]
Since u,ell are coprime, V_u=uL and ell(V)=L*ell(p) for a nonzero
linear L. Differentiating the two corresponding output coordinates
uL and ell*L shows L divides detDV: modulo L their derivative rows
are proportional. Since detDV is a nonzero power of u, L=c*u.
Consequently V_u=c*u^2, a matched pure-square component.
After scalar rescaling, KMFI141 excludes this primitive full-rank
quadratic direction. The contradiction proves the factor bound.

This does not import the old matched-class lower bound into the
unmatched class. Nonboundary linear, quadratic and cubic factors,
and powers of u, are not generally excluded here.
They must satisfy(5), the exact feedback conditions, and(3).
The eight retained permutation controls reject only those specified
matrices, not every A in(2).

## 5. Evidence and open action scope

[Helper](unmatched_quadratic.py) ·
[checker](verify_unmatched_quadratic.py) ·
[model](unmatched_quadratic_model.json) ·
[receipt](unmatched_quadratic_receipt.json) ·
[primary imports](sources/unmatched-quadratic-chart-manifest.json).

Exact controls check both signs of polarization, all four chart inverses,
every generic output matrix through determinant/adjugate identities,
the distinct source and target boundaries, global divisor homotopy,
weighted preservation, the fixed-projective-coordinate determinant factor,
and eight precisely specified unmatched feedback matrices.

The exhaustive classification here concerns determinant-compatible
directions. It is not a classification of their polynomial multipliers
or generic permutation groups. The remaining four arbitrary-A families,
their allowed low-degree factor arrangements and actual Keller field
actions remain open. Higher primitive degrees, nonhomogeneous fields,
absolute exceptional11 and the full faithful-transitive spectrum
are not decided. No new abstract occurrence group is claimed.

Formal verification, independent whole-proof replay and external review
are not claimed. The goal remains active.
