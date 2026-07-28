# The two-pair moment--nullcone program

## 1. Status and objective

This is a research program, not a proof of unrestricted
\(\operatorname{SIC}(2)\). It replaces further full-coefficient Gröbner
elimination in bidegree \((3,3)\) by one representation-theoretic statement
that makes sense in every balanced bidegree.

For \(d\geq 1\), let \(V_d\) be the space of two-pair forms of bidegree
\((d,d)\), and put

\[
 \mu_m(f)=\mathcal E_2(f^m),\qquad m\geq1.
 \tag{1.1}
\]

The central target is:

> **Moment--nullcone conjecture \(\mathrm{MN}_d\).** The common zero set of
> all \(\mu_m\) on \(V_d\) is the pair-linear one-sided nullcone.

The [complete bidegree-\((2,2)\) theorem](TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md)
proves \(\mathrm{MN}_2\). The
[bidegree-\((3,3)\) frontier](TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md) proves it
on every pure irreducible summand and on one mixed branch. The purpose of
this note is to organize those results into a route that can scale with
\(d\).

## 2. The balanced representation

Use contraction pairs \((W,Z),(V,Y)\). The diagonal
\(\mathrm{SL}_2\)-action preserving \(WZ+VY\) identifies

\[
 V_d
 \cong \operatorname{Sym}^d(\mathbb C^2)^*
       \otimes\operatorname{Sym}^d(\mathbb C^2)
 \cong \operatorname{End}(\operatorname{Sym}^d)
 \cong \bigoplus_{r=0}^{d}\operatorname{Sym}^{2r}.
 \tag{2.1}
\]

The last equality is the Clebsch--Gordan formula: the
\(\mathrm{SL}_2\)-module \(\operatorname{Sym}^d\) is self-dual and

\[
 \operatorname{Sym}^d\otimes\operatorname{Sym}^d
 \cong
 \operatorname{Sym}^{2d}\oplus
 \operatorname{Sym}^{2d-2}\oplus\cdots\oplus
 \operatorname{Sym}^{0}.
 \tag{2.2}
\]

Thus a form has irreducible coordinates

\[
 f=(F_0,F_2,F_4,\ldots,F_{2d}),
 \tag{2.3}
\]

where \(F_{2r}\) is a binary form of degree \(2r\). Every moment
\(\mu_m\) is an \(\mathrm{SL}_2\)-invariant homogeneous polynomial of
degree \(m\). The first moment is a nonzero scalar multiple of \(F_0\), so
moment vanishing removes the scalar summand.

The standard Hilbert--Mumford criterion gives an explicit description of
the nullcone \(N_d\). For direct sums of binary forms, the criterion is
recorded explicitly by Brouwer and Popoviciu in
[*Sylvester versus Gundelfinger*](https://sigma-journal.com/2012/075/):

\[
 f\in N_d
 \quad\Longleftrightarrow\quad
 F_0=0
 \ \text{and there is a common linear form }L
 \text{ such that }L^{r+1}\mid F_{2r}
 \text{ for every }r\geq1.
 \tag{2.4}
\]

Zero components impose no condition. In binary-form language, every
nonzero \(F_{2r}\) has a root of multiplicity strictly greater than half
its degree, and all these roots agree. Such a root is unique.

Equation (2.4) is exactly the pair-linear one-sided condition. Hence
\(\mathrm{MN}_d\) would imply \(\operatorname{SIC}(2)\) for every
bidegree-\((d,d)\) form. Unequal bihomogeneous bidegrees already have
nonzero total dual-minus-coordinate weight and are one-sided. Proving
\(\mathrm{MN}_d\) for all \(d\) would therefore settle the bihomogeneous
two-pair problem.

### Geometry of the target

The common-root description also gives the nullcone dimension without
elimination. Fixing \([L]\in\mathbb P^1\), the component \(F_{2r}\) has
the form

\[
 F_{2r}=L^{r+1}Q_{r-1},
 \qquad Q_{r-1}\in\operatorname{Sym}^{r-1}(\mathbb C^2),
 \tag{2.5}
\]

so its fiber has dimension \(r\). The resulting incidence space is a
vector bundle of rank

\[
 \sum_{r=1}^{d}r=\frac{d(d+1)}2
 \tag{2.6}
\]

over \(\mathbb P^1\). Away from the zero form, the destabilizing root of
the first nonzero component is unique, so the incidence map is generically
one-to-one. It follows that \(N_d\) is irreducible and

\[
 \boxed{\dim N_d=1+\frac{d(d+1)}2},\qquad
 \boxed{\operatorname{codim}_{V_d}N_d=\frac{d(d+3)}2}.
 \tag{2.7}
\]

For \(d=2\) this gives dimension \(4\); for \(d=3\), dimension \(7\).
These agree with the independent exact eliminations in the two frontier
notes.

### Why the nullcone gives the Mathieu conclusion

This implication is worth separating from the conjectural converse.
After a pair-linear change, a nullcone point is supported in the strict
one-sided positions

\[
 M_{ij}=W^{d-i}V^iZ^{d-j}Y^j,\qquad i>j.
 \tag{2.8}
\]

Every monomial in \(f^m\) therefore has total \(V\)-exponent minus total
\(Y\)-exponent at least \(m\). A fixed multiplier \(g\) changes this
difference by a bounded amount. For all sufficiently large \(m\), every
monomial of \(gf^m\) has more \(V\)-derivatives than \(Y\)-degree, and its
contraction is zero. Thus

\[
 f\in N_d
 \quad\Longrightarrow\quad
 \mathcal E_2(gf^m)=0
 \quad\text{for every fixed }g\text{ and all }m\gg0.
 \tag{2.9}
\]

Consequently moment--nullcone equality is sufficient for the Special
Image conclusion on the balanced stratum; it is not merely an invariant
theory reformulation of the pure moments.

## 3. Invariant-ring formulation

Let

\[
 S_d=\mathbb Q[V_d],\qquad
 R_d=S_d^{\mathrm{SL}_2},\qquad
 M_d=(\mu_1,\mu_2,\ldots)\subset R_d.
 \tag{3.1}
\]

The nullcone is the zero set of the positive-degree invariants
\((R_d)_+\). Consequently the set-theoretic content of
\(\mathrm{MN}_d\) is

\[
 \sqrt{M_d}=(R_d)_+,
 \tag{3.2}
\]

or, equivalently after extending to \(S_d\),

\[
 \sqrt{S_dM_d}=I(N_d).
 \tag{3.3}
\]

This suggests two scalable proof mechanisms.

1. Find finitely many moment orders and power certificates showing that a
   homogeneous generating set of \((R_d)_+\) lies in \(\sqrt{M_d}\).
2. Prove that \(R_d\) is integral over a finite moment subalgebra and that
   the fiber over the moment origin consists only of the invariant origin.

The exact radical certificates in bidegrees \((2,2)\) and on the pure
binary-sextic slice are finite instances of the first mechanism. No
uniform finite cutoff or integrality theorem is currently proved.

There are two useful qualifications.

First, \(R_d\) is finitely generated and therefore Noetherian. The ideal
generated by the infinite moment sequence is consequently generated by a
finite subcollection:

\[
 M_d=(\mu_{m_1},\ldots,\mu_{m_s})
 \quad\text{for some finite set of orders depending on }d.
 \tag{3.4}
\]

Thus a finite cutoff exists for each fixed \(d\) without any conjecture.
What is missing is an effective choice of the orders, a useful bound, and
proof that their radical is \((R_d)_+\).

Second, Krull dimension gives a sharp lower bound on how many moment
equations a global proof can use. For \(d\geq2\), a generic binary form in
the highest \(\operatorname{Sym}^{2d}\) summand has finite
\(\mathrm{SL}_2\)-stabilizer. Hence

\[
 \dim R_d=\dim V_d-\dim\mathrm{SL}_2=(d+1)^2-3.
 \tag{3.5}
\]

If \(s\) homogeneous moments define the nullcone, their ideal has radical
\((R_d)_+\), whose height is \(\dim R_d\). Krull's height theorem therefore
forces

\[
 \boxed{s\geq(d+1)^2-3.}
 \tag{3.6}
\]

This bound explains two computations already in the repository:

- for \(d=2\), the lower bound is six, and the first six moments attain it;
- for \(d=3\), any full-space moment--nullcone proof needs at least
  thirteen moment equations, although fewer can suffice on a proper
  irreducible slice.

If exactly \(\dim R_d\) moments define the nullcone, Hilbert's criterion
makes them a homogeneous system of parameters. This is the most economical
possible invariant-theoretic certificate. The binary-form version of this
criterion is described by Brouwer, Draisma, and Popoviciu in
[*The Degrees of a System of Parameters of the Ring of Invariants of a
Binary Form*](https://doi.org/10.1007/s00031-015-9353-8).

There is an essential degree-selection test before attempting the
zero-fiber geometry. Since \(R_d\) is Cohen--Macaulay in characteristic
zero, if homogeneous invariants of degrees \(e_1,\ldots,e_{\dim R_d}\)
form a system of parameters, then

\[
 H_{R_d}(t)\prod_i(1-t^{e_i})
 \tag{3.7}
\]

is the Hilbert series of an Artinian quotient and therefore is a
polynomial with nonnegative coefficients. Thus a single negative
coefficient rules out that degree sequence independently of the chosen
invariants.

For \(d=3\), the exact weight expansion gives

\[
 [t^{63}]H_{R_3}(t)\prod_{m=1}^{13}(1-t^m)=-2186.
 \tag{3.8}
\]

Consequently \(\mu_1,\ldots,\mu_{13}\), although algebraically
independent, cannot define the nullcone. Their zero fiber necessarily has
a semistable component. The least-total-degree replacement surviving the
same Hilbert test is

\[
 \mu_1,\ldots,\mu_{12},\mu_{14}.                          \tag{3.9}
\]

These corrected moments also have exact Jacobian rank thirteen. Their
proposed Hilbert numerator is nonnegative through degree \(100\), with
last observed nonzero term in degree \(76\) and zeros through degree
\(100\); this is a necessary-test result, not a proof that they are a
system of parameters.

This changes the all-degree architecture. One must first use the Molien
or weight series of \(R_d\) to select admissible moment degrees, and only
then prove that the selected zero fiber is the nullcone by the
first-component and synchronization lemmas. Consecutive initial moments
need not be the correct parameters.

## 4. Stratification by the first nonzero summand

The direct-sum nullcone condition separates into two problems: make the
first nonzero component unstable, then synchronize every higher component
with its unique destabilizing root.

For \(1\leq s\leq d\), consider the stratum

\[
 F_2=\cdots=F_{2s-2}=0,\qquad F_{2s}\ne0.
 \tag{4.1}
\]

A proof of \(\mathrm{MN}_d\) would follow from the following two lemmas.

> **First-component lemma.** On (4.1), vanishing of all moments forces
> \(F_{2s}\) to have a root \(L\) of multiplicity at least \(s+1\).

> **Synchronization lemma.** Once this root \(L\) exists, moment vanishing
> forces \(L^{r+1}\mid F_{2r}\) for every \(r>s\).

Because a root of multiplicity \(>s\) in a degree-\(2s\) form is unique,
the second lemma has a canonical flag to use. If \(F_{2s}=0\), one simply
moves to the next stratum. This avoids choosing a global normal form
before the moments have produced one.

### The quadratic anchor

The first stratum is the most useful starting point. Write

\[
 F_2=aX^2+2bXT+cT^2,\qquad
 \Delta_2=b^2-ac.
 \tag{4.2}
\]

The first-component lemma for \(s=1\) becomes the concrete target

\[
 \Delta_2\in\sqrt{M_d}.
 \tag{4.3}
\]

If \(F_2\ne0\), equation (4.3) gives \(F_2=L^2\), after which the
synchronization lemma asks successively for

\[
 L^3\mid F_4,\quad L^4\mid F_6,\quad\ldots,\quad
 L^{d+1}\mid F_{2d}.
 \tag{4.4}
\]

If \(F_2=0\), the same argument restarts with the quartic component.
Thus (4.3) is an anchor, not an assumption that every moment-zero point
has a nonzero quadratic component.

## 5. Current evidence

The evidence must be kept at its proved strength.

| locus | result | status |
|---|---|---|
| all of \(V_2\) | first six moments have the full one-sided nullcone radical | exact over \(\mathbb Q\) |
| all of \(V_3\) | \(\mu_1,\ldots,\mu_{13}\) are algebraically independent, but the degree-\(63\) Hilbert numerator coefficient is \(-2186\), so they cannot define the nullcone | exact over \(\mathbb Q\); an extra semistable zero component exists |
| all of \(V_3\) | \(\mu_1,\ldots,\mu_{12},\mu_{14}\) have exact Jacobian rank thirteen and pass the necessary Hilbert numerator test through degree \(100\) | exact over \(\mathbb Q\); corrected zero fiber still open |
| maximal-torus fixed diagonal slice in \(V_3\) | moments \(1,2,3,4\) have only the origin as a common zero, with seventh-power certificates for all four diagonal coefficients | exact over \(\mathbb Q\) |
| full non-null \(F_2\) branch in \(V_3\) | after \(F_2=2XT\), five residual-torus chart orbits cover the non-diagonal locus; \(\mu_2,\ldots,\mu_{12}\) have exact Jacobian rank eleven on every representative chart | exact over \(\mathbb Q\); all five affine zero fibers remain open |
| pure \(\operatorname{Sym}^2,\operatorname{Sym}^4,\operatorname{Sym}^6\) in \(V_3\) | moments cut out the corresponding binary-form nullcones | exact over \(\mathbb Q\) |
| \(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) in \(V_3\), with \(F_2=2cXT\) | moments through order six imply \(c^6=0\) | exact over \(\mathbb Q\) |
| \(\operatorname{Sym}^6\oplus\operatorname{Sym}^2\) in \(V_3\), with \(F_2=2cXT\) | even moments through order fourteen imply \(c^{25}=0\) | only over \(\mathbb F_{32003}\) |
| full mixed \(V_3\) | no moment--nullcone equality yet | open |

The normalized \(c\)-certificates say that the tested branches cannot have
\(\Delta_2\ne0\). They motivate (4.3), but the finite-field
\(c^{25}\) membership is not a characteristic-zero certificate and neither
slice proves the global quadratic-anchor statement.

## 6. What does and does not generalize

### Arbitrary direct sums of binary forms

The geometric half is already general. For

\[
 V=\operatorname{Sym}^{n_1}\oplus\cdots\oplus
   \operatorname{Sym}^{n_s},
 \tag{6.1}
\]

Hilbert--Mumford says that \((H_1,\ldots,H_s)\) is in the
\(\mathrm{SL}_2\)-nullcone exactly when all nonzero \(H_i\) have a common
root of multiplicity \(>n_i/2\). Thus the first-component and
synchronization architecture applies to any direct sum of binary forms.

What is special here is the invariant sequence \(\mu_m=\mathcal E_2(f^m)\).
An arbitrary direct sum has no canonical contraction moments with the
required properties. The generalization is therefore a reusable nullcone
geometry, not a universal moment--nullcone theorem.

### \(n\) contraction pairs

There is a natural higher-rank formulation. Let \(U=\mathbb C^n\) and

\[
 V_{n,d}
 =\operatorname{Sym}^d(U^*)\otimes\operatorname{Sym}^d(U)
 =\operatorname{End}(\operatorname{Sym}^dU).
 \tag{6.2}
\]

Pieri's rule gives the multiplicity-free \(\mathrm{SL}_n\)-decomposition

\[
 \boxed{
 V_{n,d}\cong
 \bigoplus_{j=0}^{d}
 V_{j(\omega_1+\omega_{n-1})}},
 \tag{6.3}
\]

where \(V_\lambda\) denotes the irreducible highest-weight module of
highest weight \(\lambda\). One derivation writes
\((\operatorname{Sym}^dU)^*\), up to a determinant twist, as the Schur
module of rectangular shape \((d^{\,n-1})\). Adding a horizontal
\(d\)-strip gives precisely

\[
 (d+j,d,\ldots,d,d-j),\qquad 0\leq j\leq d,
 \tag{6.4}
\]

which becomes \(j(\omega_1+\omega_{n-1})\) after removing the determinant
twist. When \(n=2\), formula (6.3) reduces to (2.1).
As a dimension check, Weyl's formula gives

\[
 \dim V_{j(\omega_1+\omega_{n-1})}
 =\frac{2j+n-1}{n-1}
   \binom{j+n-2}{n-2}^{\!2},
 \qquad
 \sum_{j=0}^{d}\dim V_{j(\omega_1+\omega_{n-1})}
 =\binom{n+d-1}{d}^{\!2}
 =\dim V_{n,d}.
\]

Define

\[
 \mu_{n,m}(f)=\mathcal E_n(f^m).
 \tag{6.5}
\]

Because \(f\) is balanced, this is a scalar
\(\mathrm{SL}_n\)-invariant. The formal higher-rank question is

\[
 \mathrm{MN}_{n,d}:\qquad
 V(\mu_{n,1},\mu_{n,2},\ldots)
 \stackrel{?}{=}\mathcal N(V_{n,d}).
 \tag{6.6}
\]

The nullcone is now described by an arbitrary destabilizing one-parameter
subgroup, or equivalently by a weighted flag in \(U\); there is generally
no reduction to one common point of \(\mathbb P^1\).

The easy direction still holds in every rank:

\[
 \mathcal N(V_{n,d})
 \subseteq V(\mu_{n,1},\mu_{n,2},\ldots),
 \tag{6.7}
\]

because every positive-degree invariant vanishes on the nullcone.
Moreover, a balanced nullcone point satisfies eventual mixed contraction
vanishing. Indeed, choose a one-parameter subgroup for which every weight
of \(f\) is positive. Weights in \(f^m\) then grow at least linearly in
\(m\). A fixed \(g\) contributes bounded weight, while every monomial that
survives \(\mathcal E_n\) has nonnegative residual coordinate exponents
whose total degree is fixed by the bidegree of \(g\), hence belongs to a
finite set of weights. Equivariance gives a contradiction for large
\(m\). Therefore

\[
 f\in\mathcal N(V_{n,d})
 \Longrightarrow
 \mathcal E_n(gf^m)=0\quad(m\gg0).
 \tag{6.8}
\]

### Exact obstruction in three or more pairs

The converse (6.6) is false as soon as \(n\geq3\), already for \(d=2\).
The repository's
[four-term three-pair counterexample](THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)

\[
 f=\tau(t-y)(wz+vt),\qquad g=y
 \tag{6.9}
\]

lies in \(V_{3,2}\) and satisfies

\[
 \mathcal E_3(f^m)=0,\qquad
 [t]\mathcal E_3(gf^m)=(-1)^{m-1}(m+1)!\,m!
 \ne0
 \quad(m\geq1).
 \tag{6.10}
\]

If \(f\) were in the \(\mathrm{SL}_3\)-nullcone, (6.8) would force the
second expression eventually to vanish. Hence \(f\) is a moment-zero
semistable point and

\[
 \boxed{\mathrm{MN}_{3,2}\text{ is false}.}
 \tag{6.11}
\]

Ignoring additional contraction pairs preserves both identities in
(6.10), so the same argument proves

\[
 \boxed{\mathrm{MN}_{n,2}\text{ is false for every }n\geq3.}
 \tag{6.12}
\]

This is the decisive boundary of the generalization. The plausible theorem
is rank two in every balanced degree, not every rank in every degree.
Rank two is special because all one-parameter subgroups are conjugate and
instability becomes a unique common-root condition for binary forms.

### Nonhomogeneous forms and positive characteristic

The balanced theorem would cover each bidegree-\((d,d)\) stratum, and a
single unequal bidegree is already one-sided. It does not automatically
cover a sum of several bidegrees: mixed products can cancel the central
dual-minus-coordinate grading. A nonhomogeneous extension therefore needs
a separate multigraded argument.

All statements in this program are over characteristic zero, with
nullcone tests made after scalar extension to an algebraic closure. In
small positive characteristic the displayed Clebsch--Gordan decomposition
need not be semisimple, so no characteristic-free version is asserted.

## 7. The next attack

The next calculations should discover small invariant or covariant
identities, rather than eliminate the full coefficient space.

1. **Construct uniform Clebsch--Gordan coordinates.** Express the
   projections \(f\mapsto F_{2r}\) and the contractions \(\mu_m\) by
   transvectants, with \(d\) retained as a parameter.
2. **Select admissible parameter degrees.** Compute the Molien/weight
   Hilbert series before attacking a zero fiber. In \(d=3\), degrees
   \(1,\ldots,13\) are now excluded and the first corrected target is
   \((1,\ldots,12,14)\). Seek a uniform rule for choosing
   \((d+1)^2-3\) moment orders whose proposed Hilbert numerator is
   nonnegative.
3. **Test the corrected minimal zero fiber.** Determine whether
   \(V(\mu_1,\ldots,\mu_{12},\mu_{14})\) equals the nullcone, preferably
   by invariant/covariant saturation. If it fails, isolate the extra
   semistable component and move only to another Hilbert-compatible
   degree set; do not merely extend a consecutive cutoff.
4. **Prove the quadratic anchor.** Search for a finite identity
   \(\Delta_2^N\in(\mu_{m_1},\ldots,\mu_{m_k})\) in invariant coordinates,
   using the explicit global projection
   \[
   r_0=(3c_{10}+2c_{21}+3c_{32})/10,\quad
   r_1=(-9c_{00}-c_{11}+c_{22}+9c_{33})/20,\quad
   r_2=-(3c_{01}+2c_{12}+3c_{23})/10.
   \]
   The torus-fixed slice is now closed exactly; the next test must allow
   generic higher-component weights.
5. **Prove one synchronization step.** On the incidence chart
   \(F_2=L^2\), show that the moments force the forbidden coefficients of
   \(F_4\) to vanish, equivalently \(L^3\mid F_4\). Formulate the
   certificate covariantly so that \(L\) can be eliminated without a large
   global Gröbner basis.
6. **Generalize the step.** With
   \(L^{j+1}\mid F_{2j}\) for \(j<r\), use the lowest remaining
   \(\mathrm{SL}_2\)-weights in the moments to force
   \(L^{r+1}\mid F_{2r}\).
7. **Handle a zero anchor.** Repeat the same construction on the strata
   \(F_2=0\), then \(F_2=F_4=0\), so the proof does not discard components
   whose lower summands vanish.

The first decisive milestone is therefore not a rational reconstruction
of the isolated \(c^{25}\) calculation. The consecutive dimension-sized
set is now excluded and the corrected set (3.9) is the first viable
target. The next milestone is a global \(d=3\) quadratic-anchor
certificate for that corrected ideal, stated in
\(\mathrm{SL}_2\)-invariant form, followed by the first common-root
synchronization certificate. Those two identities would expose the
pattern needed for arbitrary \(d\).

### Decision tree after the present tests

The best next calculation is the localized non-null-quadratic problem.
On \(\Delta_2\ne0\), use \(\mathrm{SL}_2\) to put

\[
 F_2=2cXT,\qquad c\ne0,
\]

then use overall homogeneity to set \(c=1\). The scalar component has
already been removed by \(\mu_1\), leaving twelve variables
\((F_4,F_6)\) and the twelve corrected equations
\(\mu_2,\ldots,\mu_{12},\mu_{14}\). The residual diagonal torus gives
weights

\[
 (6,4,2,0,-2,-4,-6;\;4,2,0,-2,-4).
\]

There are already \(246354\) residual-weight-zero monomials of degree at
most thirteen, and the corrected system also contains order fourteen, so
a raw expanded Gröbner calculation remains large. The weight grading is
essential: use a sparse invariant-monomial or straight-line representation
rather than first expanding in twelve ordinary variables. This slice is
nevertheless substantially smaller and more structured than saturating
the sixteen-coefficient ideal directly.

The outcomes have clear implications:

1. If the dehomogenized ideal is the unit ideal over \(\mathbb Q\), then
   \(\Delta_2\) lies in the radical of the corrected moment ideal: the
   global quadratic anchor is proved for the viable minimal set.
2. If it has a component, test whether the component survives all higher
   moments. A surviving exact point or recurrence is a candidate
   SIC(2) counterexample; a component killed by later moments only shows
   that the corrected set is still not sufficient. Its geometry should
   determine the next Hilbert-compatible replacement.
3. Once the non-null branch is excluded, move to the incidence chart
   \(F_2=L^2\) and prove \(L^3\mid F_4\). This is the first synchronization
   lemma.
4. Only after these two steps should the same identities be interpolated
   in \(d\). A full raw Gröbner basis in sixteen coefficients remains the
   least informative route.

For computation, the preferred order is: modular sparse elimination on
the \(c=1\) slice to estimate feasibility and certificate degree; exact
rational reconstruction of a unit certificate; then an independent exact
reduction. A finite-field unit ideal is evidence only until the rational
certificate is reconstructed.

## 8. Claim boundary

This program does not prove \(\mathrm{MN}_3\), \(\mathrm{MN}_d\) for
general \(d\), or unrestricted \(\operatorname{SIC}(2)\). It gives:

- an exact all-\(d\) reformulation of the desired conclusion as a common
  high-multiplicity-root condition;
- a stratified proof architecture with two explicit lemma families;
- the exact Hilbert-series obstruction to consecutive degrees
  \(1,\ldots,13\), plus the algebraically independent corrected candidate
  \(1,\ldots,12,14\);
- a precise first invariant target, \(\Delta_2\in\sqrt{M_d}\);
- a rule for using the existing exact and finite-field calculations
  without promoting experiments to theorems.

It also proves the structural statements (2.7), (3.4)--(3.8), the
corrected Jacobian independence in (3.9), and the higher-rank obstruction
(6.11)--(6.12). These are consequences of exact calculation and standard
invariant theory together with the repository's proved three-pair
counterexample; they do not improve the open SIC(2) status.

The cited literature supplies the Image/Mathieu framework, the
Hilbert--Mumford and binary-form nullcone criteria, and the invariant-ring
tools. It does not assert the moment--nullcone conjecture \(\mathrm{MN}_d\);
that conjecture and the proposed anchor/synchronization proof are specific
to this repository program.

## 9. Sources

- A. van den Essen, D. Wright, and W. Zhao,
  [*On the Image Conjecture*](https://arxiv.org/abs/1008.3962),
  J. Algebra 340 (2011), 211--224. This supplies the Image/Mathieu
  framework and the contraction setting.
- H. Derksen,
  [*Constructive Invariant Theory*](https://www.cse.iitb.ac.in/~sohoni/CS782/DerksenGIT.pdf),
  Sections 2.2 and 4. This is used for the invariant-theoretic nullcone,
  Hilbert--Mumford criterion, and finite invariant bounds.
- A. E. Brouwer and M. Popoviciu,
  [*Sylvester versus Gundelfinger*](https://sigma-journal.com/2012/075/),
  SIGMA 8 (2012), 075. This records the common high-multiplicity-root
  criterion for direct sums of binary forms and Hilbert's nullcone
  criterion for systems of parameters.
- A. E. Brouwer, J. Draisma, and M. Popoviciu,
  [*The Degrees of a System of Parameters of the Ring of Invariants of a
  Binary Form*](https://doi.org/10.1007/s00031-015-9353-8),
  Transform. Groups 20 (2015), 953--967.
- C. D. Long,
  [*Counterexamples to the \(xz\)-Conjecture and the Mathieu Conjecture for
  \(SU(2)\)*](https://arxiv.org/abs/2607.19012), arXiv:2607.19012v1
  (2026). The repository's three-pair witness is a bihomogeneous lift of
  Long's \(SU(2)\) seed.
- M. Müger and L. Tuset,
  [*The Mathieu Conjecture for \(SU(2)\) Reduced to an Abelian
  Conjecture*](https://arxiv.org/abs/2210.06582),
  Indag. Math. 35 (2024), 114--118.
