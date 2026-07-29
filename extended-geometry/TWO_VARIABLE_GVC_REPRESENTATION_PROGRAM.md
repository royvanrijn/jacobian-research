# Two-variable GVC and rank-efficient Laplacian polarization

## 1. Status

This remains a research program rather than a proof of unrestricted
\(\operatorname{GVC}(2)\), but its balanced homogeneous target is now
closed.  The
[balanced cubic theorem](TWO_VARIABLE_CUBIC_GVC_THEOREM.md) proves that
the first four pure moments force every cubic Segre point into the
one-sided nullcone, with the explicit mixed cutoff \(m>\deg Q\).
More generally, the
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) proves the GVC
conclusion in every balanced homogeneous degree: every binary operator
symbol splits into linear factors, and translated complete polarization
converts its pure values to constant terms of powers of one Laurent
polynomial.  Retaining the translation variable removes the earlier
degree restriction: homogeneous binary operators satisfy GVC for
arbitrary \(P\).

The later
[separable escape obstruction](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
closes a further nonhomogeneous class.  If \(r\) is the lowest positive
homogeneous order of a binary constant-coefficient operator \(\Lambda\),
then GVC holds for every \(P\) of degree at most \(r\).  It also proves
GVC for arbitrary \(P\) when \(\Lambda\) is a homogeneous split factor
times a differential operator with nonzero constant term.  It further proves
that separated multiplicative auxiliary specialization, rank-one
dilation, and nonlinear Segre substitution cannot convert the new
two-pair witness: the output remains rank one, whereas the witness has
rank five.  Coefficient extraction can sum rank-one channels only by
losing multiplicativity, so it supplies no automatic all-order transfer.
Moreover, fixed linear translation followed by one diagonal coefficient
represents exactly a product of powers of linear symbols.  The first
remaining operator architectures must therefore use coupled coefficients,
nonlinear translation, or an irreducible nonhomogeneous symbol.

The two continuing targets are:

1. treat genuinely nonhomogeneous pairs with \(\deg P\) greater than the
   lowest positive operator order and with no split homogeneous
factor-unit decomposition, and sharpen the split-symbol theorem to
   finite pure-moment/nullcone certificates; and
2. replace the third-order Dvorsky operator by a quadratic
   constant-coefficient operator using as few auxiliary variables as
   possible.

Neither target changes the certified dimension ledger.  In particular, the
current upper bounds remain five variables for unrestricted
constant-coefficient GVC and forty variables for ordinary-Laplacian GVC.

The first target has an additional associated-graded gate.  For arbitrary
nonhomogeneous \(\Lambda\), every fixed number of leading homogeneous
layers of \(\Lambda^m(QP^m)\) vanishes eventually under the pure premise.
Consequently a counterexample cannot live on a fixed leading face: its
mixed defect must move to unbounded depth as \(m\) grows.
If the operator has only linear and quadratic pieces and its linear part
is nonzero, the first two moments close arbitrary \(P\).  In normalized
coordinates, the unique highest \(y\)-degree term of the second moment is
the square contribution from \(C^2P_{yy}^2\), forcing \(P=ay+b\).
The same product-defect argument closes every separated drift
\(\partial_x+h(\partial_y)\), even for formal \(h\), since derivative
series act locally finitely on polynomials.  Formal Weierstrass division
factors every binary symbol with nonzero linear part as
\(U(\xi,\eta)(\xi+q(\eta))\), where \(U(\partial)\) is a locally finite
differential automorphism.  This closes the entire lowest-order-one
frontier for arbitrary \(P\).  Exact cubic and quartic jet computations
give independent finite regressions of the factorization argument.
Hence a remaining GVC(2) counterexample must have lowest positive
operator order at least two.  The first three pure equations also close
the entire \(r=2,\deg P\le3\) cell: the double-line quadratic orbit has
one apparent second-moment cancellation, killed exactly at moment three,
and every surviving branch in both quadratic orbits has a strict weighted
degree separator.  Thus the first \(r=2\) target has \(\deg P\ge4\).
The quartic analysis below closes that endpoint as well.
In that quartic target, the distinct-root quadratic orbit and the
double-line \(xy^3\) top form close through moment three and weighted
cutoffs.  Only the double-line \(y^4\) top form remains; eliminating one
coefficient from its second and third moments leaves finitely many ratios
on an explicit sextic.  Its fourth-moment octavic is coprime to that
sextic.  Moreover, weights
\(\operatorname{wt}(x)=2,\operatorname{wt}(y)=1\) show that the full pure
moments equal those of this three-term weighted face: every higher
operator jet has excess weight and every lower polynomial term has
deficient weight.  Hence the entire \(r=2,\deg P\le4\) cell is closed,
and the first quadratic-leading target has degree at least five.
For \(r=3,\deg P=4\), the triple-root cubic-symbol orbit closes through
moment four and three weighted-face reductions.  The double-root orbit
closes through moment three.  For the squarefree orbit, the first two
leading moments leave three fourth-power tips in one root-permutation
orbit; its \(x^4\) tip has the same weight-four correction face as the
double-root case and dies at moment three.  Hence every binary
constant-coefficient operator is GVC-safe for \(\deg P\le4\).  A
counterexample must have degree at least five; initially its lowest
positive operator order could only be \(2,3,\) or \(4\).  For \(r=2\),
the first two
leading equations leave only
\((\partial_x\partial_y,x^5)\),
\((\partial_x^2,xy^4)\), and
\((\partial_x^2,y^5)\).  The full correction calculation closes all
three: the distinct-root tip has a linear-factor/high-order cutoff, the
\(xy^4\) tip dies on two nested weighted faces, and the \(y^5\) tip
reduces to six ratios, four killed at moment three and two already
one-sided.  Hence the first quadratic-leading target has degree at least
six, and a degree-five counterexample must have \(r=3\) or \(4\).
For \(r=3\), the first two leading equations leave eight top-form normal
forms: four triple-root, three double-root, and one squarefree.  This is
the next finite correction calculation.  The \(r=4\) row has a
continuous squarefree-quartic cross-ratio and is the less discrete branch.
Although a squarefree tangent cone formally factors into smooth branches,
this does not transport powers: a multiplicative conjugate of a
derivation is a derivation, so a constant-coefficient umbral conjugate
must have linear symbol.  Exploiting formal branches therefore requires a
new multilinear identity, not an algebra-coordinate straightening.

The later
[two-pair SIC counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
lies in bidegree \((4,4)\), but its \(5\times5\) coefficient matrix has
determinant \(48\).  It is therefore a full-rank nonseparable point of
\(\operatorname{End}(\operatorname{Sym}^4)\), not a rank-one Segre point
\(A\otimes P\).  It falsifies unrestricted SIC(2) without changing either
GVC target above.

## 2. Balanced homogeneous GVC is a Segre slice of SIC(2)

Let \(k\) be a characteristic-zero field and let \(U\) be two-dimensional.
A homogeneous constant-coefficient operator of order \(d\) has a symbol

\[
 A\in\operatorname{Sym}^d U,
 \qquad
 \Lambda=A(\partial),
 \]

while a homogeneous polynomial of degree \(d\) is an element

\[
 P\in\operatorname{Sym}^d(U^*).
 \]

Introduce dual variables \(\zeta\) and coordinate variables \(z\), and put

\[
 f(\zeta,z)=A(\zeta)P(z).
 \tag{2.1}
\]

Then \(f\) is a rank-one point of

\[
 V_d
 =\operatorname{Sym}^d U\otimes\operatorname{Sym}^d(U^*)
 =\operatorname{End}(\operatorname{Sym}^d U).
 \tag{2.2}
\]

If \(\mathcal E_2\) is the two-pair contraction, the definitions give the
exact identities

\[
 \boxed{
 \mathcal E_2(f^m)=\Lambda^m(P^m),\qquad
 \mathcal E_2(Q(z)f^m)=\Lambda^m(QP^m).
 }
 \tag{2.3}
\]

Thus balanced homogeneous two-variable GVC is not merely analogous to the
SIC(2) problem.  It is the restriction of the SIC moments to the affine
Segre cone

\[
 \Sigma_d
 =\{A\otimes P\}\subset V_d.
 \tag{2.4}
\]

The multiplier class in GVC is smaller: \(Q\) depends only on the
coordinate variables, whereas SIC allows an arbitrary multiplier in both
sets of variables.  Consequently a moment--nullcone theorem on all of
\(V_d\) proves the balanced homogeneous GVC conclusion, but the restricted
GVC problem on \(\Sigma_d\) may be strictly easier.

After self-dualizing the binary representation, Clebsch--Gordan gives

\[
 A\otimes P
 \longmapsto
 (F_0,F_2,\ldots,F_{2d})
 \in\bigoplus_{r=0}^d\operatorname{Sym}^{2r}U.
 \tag{2.5}
\]

The components \(F_{2r}\) are the joint transvectants of \(A\) and \(P\),
and the quantities in (2.3) are their invariant contractions.  The
one-sided nullcone criterion from the
[two-pair moment--nullcone program](TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md)
becomes

\[
 F_0=0,\qquad
 F_{2r}=L^{r+1}G_{r-1}\quad(1\leq r\leq d)
 \tag{2.6}
\]

for one common linear form \(L\).  This supplies the same two-stage proof
architecture:

1. force the lowest nonzero transvectant into its binary-form nullcone; and
2. synchronize the unique high-multiplicity root through every higher
   transvectant.

The additional rank-one equations defining \(\Sigma_d\) should be imposed
before elimination.  Eliminating in the full \((d+1)^2\)-dimensional SIC
space discards the main advantage of the GVC specialization.

## 3. The cubic Segre slice is closed

The complete bidegree-\((2,2)\) theorem proves SIC(2) on all of \(V_2\), so
it proves the balanced quadratic GVC conclusion a fortiori.  It also shows
that any bidegree-preserving linear pair compression of the compact
three-pair witness to two pairs cannot work: linear pair substitution
preserves bidegree \((2,2)\), and the resulting form lies in the proven-safe
stratum.  This statement does not cover a linear substitution that mixes
dual and coordinate variables and thereby destroys the bigrading.

The next case was the cubic Segre slice.  It now satisfies the exact
equality

\[
 \boxed{
 \Sigma_3\cap
 V(\mu_1,\mu_2,\mu_3,\mu_4)
 =
 \Sigma_3\cap N_3.
 }
 \tag{3.1}
\]

Indeed, normalize the nonzero cubic symbol to its triple-root,
double-root, or squarefree orbit.  The first orbit closes after moment one,
the second after moment three, and the squarefree orbit after moment four.
In the squarefree orbit the surviving polynomial is one of three pure
cubes, each annihilated by one linear factor of the operator.  Every branch
then gives

\[
 \Lambda^m(QP^m)=0\qquad(m>\deg Q).
 \tag{3.2}
\]

The full proof and exact moment identities are in
[`TWO_VARIABLE_CUBIC_GVC_THEOREM.md`](TWO_VARIABLE_CUBIC_GVC_THEOREM.md).
Since \(\Sigma_3\) has dimension seven inside the sixteen-dimensional
\(V_3\), this also proves that the extra semistable component of the first
thirteen full \((3,3)\) moments does not meet the Segre cone.  The first
open balanced homogeneous separable case is now degree four.

The degree-four frontier has since narrowed further.  The
[low-root theorems](TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md) prove every
balanced degree when the operator symbol has at most two distinct roots,
using the one-variable Laurent constant-term theorem.  In degree four,
moments through order five also cut out the one-sided locus on the
\((2,1,1)\) three-root orbit, with five exact fourth-power certificates.
On a nonempty Zariski-open set of squarefree quartic cross-ratios, moments
through order six cut out the four annihilator lines; this follows from an
exact reduced four-point fiber at cross-ratio \(2\) and proper-family upper
semicontinuity.
The later [split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) closes the
GVC conclusion on the one-parameter squarefree quartic orbit and in every
higher balanced degree.  What remains on that quartic orbit is the stronger
finite-moment nullcone question, not eventual coordinate-multiplier
vanishing.

## 4. What the SIC(2) counterexample changes

The bidegree-\((2,2)\) theorem closes more than the natural four-parameter
compression family: it closes every two-pair form in that stratum.
The successful witness uses balanced bidegree \((4,4)\) and is genuinely
nonseparable.  Thus the earlier list of possible escape mechanisms is now
resolved as follows:

1. bidegree \((2,2)\) remains completely safe;
2. full bidegree \((3,3)\) remains undecided;
3. bidegree \((4,4)\) already fails at a full-rank nonseparable point; and
4. unequal/nonhomogeneous degree and nonlinear polarization are no longer
   needed to falsify SIC(2), though they remain possible sources of smaller
   support or degree.

A single unequal bihomogeneous component is already one-sided by its
nonzero central weight.  The viable unequal-degree option is therefore a
nonhomogeneous mixture whose different central weights can interact in
powers; unequal bidegree by itself is not an open counterexample stratum.

For ordinary homogeneous GVC, (2.1) is necessarily separable.  The new
counterexample therefore belongs only to the larger SIC problem, while the
Segre restriction remains the structural advantage available to GVC.

Moreover, the rank obstruction is stable under every separated
multiplicative conversion.  Evaluation after auxiliary polarization,
separated restriction of a rank-one dilation, and nonlinear substitutions
on the dual and coordinate factors all preserve the product form.  They
therefore cannot produce the rank-five matrix (1.6).  A sum needs at least
five rank-one channels.  Auxiliary coefficient extraction is not
multiplicative, so a construction using it must prove new identities for
every power rather than inherit the witness identities by specialization.

For a nonhomogeneous operator
\(\Lambda=\Lambda_r+\cdots+\Lambda_d\), the lowest-order theorem further
shows that \(\deg P\leq r\) is safe.  In the pure premise, all higher
operator pieces are killed by degree and one obtains
\(\Lambda_r^m(P_r^m)=0\).  In the mixed expression only a bounded number
of higher-order pieces can occur; the split-symbol Newton gap for
\(\Lambda_r\) absorbs these bounded defects.  Thus nonhomogeneity alone is
not an escape mechanism.

## 5. Rank-efficient quadraticization of the Dvorsky operator

The Dvorsky--Long witness uses

\[
 \Lambda
 =\partial_t(\partial_a\partial_d-\partial_b\partial_c),
 \qquad
 P=(t+c)(ad+bt),
 \tag{5.1}
\]

so its operator symbol is the factored cubic

\[
 \sigma_\Lambda
 =\xi_t(\xi_a\xi_d-\xi_b\xi_c).
 \tag{5.2}
\]

This does not directly give ordinary-Laplacian GVC, whose operator symbol
is a nondegenerate quadratic form.  The promising conversion problem is to
find a polynomial lift to \(N\) variables with a nondegenerate quadratic
operator \(\widetilde\Delta\), a polynomial \(\widetilde P\), and a fixed
multiplier \(\widetilde Q\), together with nonzero constants \(c_m,d_m\),
such that the all-order identities

\[
\begin{aligned}
 \widetilde\Delta^m(\widetilde P^m)
   &=c_m\,\iota\!\left(\Lambda^m(P^m)\right),\\
 \widetilde\Delta^m(\widetilde Q\,\widetilde P^m)
   &=d_m\,\iota\!\left(\Lambda^m((-c)P^m)\right)
\end{aligned}
\tag{5.3}
\]

hold for every \(m\), for an explicit injective polynomial substitution
\(\iota\).  After a linear change over an algebraic closure, a
nondegenerate quadratic constant-coefficient operator is an ordinary
Laplacian.

The optimization objective is the auxiliary rank, hence \(N\), rather than
the size of a generic polarization.  The factorization (5.2) should be
retained as a \(1+2\) tensor network: the determinant factor has rank two
and should be quadraticized as one block, not expanded and polarized
monomial by monomial.

The useful machinery from the \(HC_4\) program is methodological:

1. choose a nondegenerate quadratic pivot;
2. eliminate auxiliary blocks by a Schur complement;
3. stratify the residual coupling by rank and constant-kernel dimension;
4. use source/dual bidegrees to isolate faces that cannot cancel; and
5. synchronize any surviving moving kernel before increasing the ansatz.

A constant Schur complement of quadratic symbols is still quadratic, so it
cannot by itself reproduce the cubic operator.  Any successful
construction must act jointly on the operator and polynomial sides, or use
a nonlinear polarization.  This is the first structural gate.

The first calculation should therefore solve the universal coefficient
problem for one auxiliary quadratic block while preserving both lines of
(5.3).  That calculation is now complete for every polynomial or formal
hyperplane lift of the canonical one-pair completion.  With one new
variable \(s\), consider the
nondegenerate completion
\[
 \widetilde\Delta
 =\partial_a\partial_d-\partial_b\partial_c+\partial_t\partial_s.
\tag{5.4}
\]
For every polynomial or formal \(F\) restricting to \(P\) on \(s=0\), the
first pure equation \(\widetilde\Delta F=0\) forces
\[
 \left.
 \widetilde\Delta^2(F^2)
 \right|_{a=b=c=d=s=0}
 =12t^2-8\rho t
\tag{5.5}
\]
for one transverse first-jet coefficient \(\rho\).  Thus every such
six-variable hyperplane lift fails at moment two.  The proof and exact
certificate are in
[the one-pair obstruction](DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md).
Degree mixing cannot change the coefficient \(12\).

The next ansatz must therefore add a second auxiliary block, replace the
hyperplane restriction by a nonlinear polarization, or use a different
quadratic completion with genuinely new cross terms.
Only after all-order identities are derived should the block rank be
minimized.  A finite moment match, or preservation of the pure identity
without the mixed defect, does not lower the ordinary-Laplacian endpoint.

## 6. Laplacian promotion rule

The cubic Segre theorem and the one-pair obstruction have exact checkers.
No successful quadratic ordinary-Laplacian lift is currently attached to
the second branch.  Such a result may change the Laplacian status ledger
only after it supplies:

1. an explicit characteristic-zero construction;
2. an all-order proof of the pure vanishing;
3. an all-order or infinite-tail proof of the fixed-multiplier defect;
4. a proof that the quadratic operator is nondegenerate and hence linearly
   equivalent to the ordinary Laplacian; and
5. the exact number of variables and auxiliary-block rank.
