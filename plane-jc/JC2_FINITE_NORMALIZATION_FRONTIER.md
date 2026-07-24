# JC(2) finite-normalization frontier: the cubic cusp gate

## Status

The surface finite-normalization theorem is a genuinely useful unconditional
reduction.  If

\[
F=(P,Q):\mathbb A^2\longrightarrow\mathbb A^2
\]

is Keller and \(B\) is the normalization of \(k[P,Q]\) in \(k(x,y)\), then
\(B\) is finite locally free, hence free, over \(k[P,Q]\).  Every closed
fiber has the generic length.  This removes the closed-point module-flatness
defect which remains possible for a normal threefold.

This note audits the four proposed next targets.  Its main conclusion is
that the first and fourth targets need a correction:

> finite freeness, tame transverse ramification, and purity of the
> divisorial Jacobian do **not** imply that a boundary curve maps
> immersively to its (possibly singular) target image.

The obstruction is already present in the universal cubic cusp cover.  In
geometric degree three it is the minimal local model left by the generic
sheet ledger.  The next useful theorem is therefore not an abstract
Riemann--Hurwitz or unit-rank statement.  It is a global exclusion of this
cusp packet using the fact that the distinguished complement is
\(\mathbb A^2\).

## 1. What finite freeness gives

Put

\[
A=k[u,v],\qquad
B=\operatorname{Norm}_{A}k(x,y),\qquad
d=[k(x,y):k(u,v)].
\]

The verified input is:

1. \(B\) is a finite free \(A\)-module of rank \(d\);
2. every fiber of \(\operatorname{Spec}B\to\mathbb A^2\) has length \(d\);
3. the missing boundary primes freely generate \(\operatorname{Cl}(B)\);
4. \(B^\times=k^\times\);
5. every target nonproperness component has a positive affine-sheet
   contribution; and
6. every missing-boundary curve has rational projective normalization.

These statements are global and special to the distinguished open
\(\mathbb A^2\).  They must not be replaced by a merely local finite-flat
surface model.

The remaining information has three different levels:

| level | data | finite at fixed \(d\)? |
| --- | --- | --- |
| generic DVR | rows \((e_i,f_i)\) and affine residue degrees | yes |
| curve/conductor | punctures, singular branches, conductor weights and special-fiber packets | not without a curve-degree or conductor bound |
| resolved surface | boundary tree, self-intersections, discrepancies, pole vector and target grouping | not without a minimality theorem |

Thus “enumerate the degree-three ledgers” is finite only at the first level.
The word *decorated* must specify a bound or a canonical minimal model at the
second and third levels.

## 2. The cubic cusp countermodel to automatic residue immersion

Consider

\[
\pi:\mathbb A^2_{t,u}\longrightarrow\mathbb A^2_{u,v},
\qquad
v=t^3+ut.
\tag{2.1}
\]

Its coordinate algebra is

\[
k[u,v,t]/(t^3+ut-v),
\tag{2.2}
\]

which is finite free of rank three over \(k[u,v]\), with basis
\(1,t,t^2\).  The Jacobian is

\[
J_\pi=3t^2+u.
\tag{2.3}
\]

The ramification curve

\[
E=(u+3t^2=0)
\]

is smooth and occurs with the tame coefficient \(e-1=1\).  Its image is the
cusp

\[
C=(4u^3+27v^2=0),
\tag{2.4}
\]

parametrized by

\[
t\longmapsto(-3t^2,-2t^3).
\tag{2.5}
\]

The map (2.5) is not immersive at \(t=0\), although (2.3) has no residual
divisorial factor after removing \(E\).  The pullback of the discriminant is

\[
4u^3+27(t^3+ut)^2
=(u+3t^2)^2(4u+3t^2).
\tag{2.6}
\]

Hence the generic degree-three row over \(C\) is exactly

\[
\boxed{(e,f)_{\rm boundary}=(2,1),\qquad
       (e,f)_{\rm companion}=(1,1).}
\tag{2.7}
\]

At the cusp the three roots coalesce and the fiber is
\(k[t]/(t^3)\), of length three.  Finite flatness therefore does not exclude
it.

The failed local inference is easy to locate.  At a smooth point of \(C\)
one can split the logarithmic differential into a normal factor and the
tangential differential of \(E\to C\).  At the cusp there are no regular
target parameters in which \(C\) is a coordinate axis.  A unit residual
Jacobian therefore does not force the differential of the map to the
singular embedded curve to be nonzero.

This example is not a Keller counterexample: deleting \(E\) from the source
does not leave \(\mathbb A^2\).  Indeed \(J_\pi\) becomes a nonconstant unit
on that complement.  Precisely those global unit and class-group failures
show what a successful JC(2) argument still has to use.

The exact symbolic regression is
[`cas/test_cubic_cusp_local_model.py`](cas/test_cubic_cusp_local_model.py).

## 3. Correct residue statements

Let \(E^\nu\to C^\nu\) be the induced map between normalizations.  There are
two different assertions:

\[
E^\nu\longrightarrow C^\nu
\quad\text{is unramified}
\tag{3.1}
\]

and

\[
E^\nu\longrightarrow C\hookrightarrow\mathbb A^2
\quad\text{is immersive}.
\tag{3.2}
\]

Assertion (3.2) implies (3.1), but the converse fails at a singular
unibranch point of \(C\).  The cubic cusp has a degree-one map in (3.1) and
still fails (3.2).

Riemann--Hurwitz controls (3.1).  If \(E^\nu\) has \(s\) punctures and
residue degree \(f\) over \(C^\nu\simeq\mathbb A^1\), the affine
ramification of (3.1) has degree

\[
f+s-2.
\tag{3.3}
\]

Consequently (3.1), if proved, forces \(f=s=1\).  It says nothing by itself
about the conductor of \(C\) or the differential of \(C^\nu\to C\).

There are three viable corrected targets, in increasing strength:

1. prove (3.1) for every boundary normalization;
2. prove that every target nonproperness component is seminormal, or at
   least has no singular unibranch point met by a ramified boundary sheet;
3. directly exclude the completed cubic cusp packet inside a finite-free
   normalization whose distinguished complement is \(\mathbb A^2\).

The third is the sharp degree-three target.

## 4. Degree-three ledger

Let \(C\) be a target component carrying ramification.  The affine-sheet
budget and purity force the unique generic ledger

\[
3=2+1,
\tag{4.1}
\]

namely (2.7).  In particular the boundary residue degree and puncture count
are already \(f=s=1\); no general residue-immersion theorem is needed for
that numerical conclusion.

The decorations split as follows.

### 4.1 Normal target

If \(C\) is normal, then \(C\simeq\mathbb A^1\).  The known affine-line
component theorem excludes this case for a nonsingular plane polynomial
map.

### 4.2 Multibranch target singularity

If two distinct regular boundary points of transverse index two lie over
one target value, their local fiber lengths contribute at least four.  This
contradicts rank-three flatness.  This excludes the clean two-branch packet.

It does not automatically cover a singular point of the normal source,
non-Cartier boundary, or a packet whose branches meet only after contraction.
Those cases require a local normal cubic-algebra classification.

### 4.3 Unibranch target singularity

This is the essential survivor.  The cusp model (2.1) realizes:

- one boundary point;
- transverse generic index two;
- residue degree one;
- one affine companion sheet generically; and
- a curvilinear length-three special fiber.

Higher cusps and more complicated conductor weights are not bounded by the
generic ledger.  A complete degree-three enumeration should therefore be
formulated as a classification of normal finite-flat cubic algebras with
one deleted ramification prime and a distinguished étale
\(\mathbb A^2\)-open, not as a list of integer partitions.

### 4.4 Other nonproperness components

A degree-three cover can also have target components carrying only
unramified missing-boundary rows.  The six coarse signatures currently
printed by the atlas are signatures for one arbitrary target component;
only one of them is the ramified row (4.1).  A global ledger must retain all
unramified components and their intersections with the ramified component.

## 5. The combined lattice attack

The class-group, intersection, and log-canonical constraints should be
compiled into one certificate, because none is independently decisive.
On a smooth SNC completion \(S\setminus D=\mathbb A^2\), with boundary
intersection matrix \(Q\), write

\[
K_S=\sum k_iD_i,\qquad
\bar F^*L_\infty=\sum p_iD_i.
\]

The intrinsic equations include

\[
\det Q=(-1)^{r-1},\qquad
\operatorname{inertia}(Q)=(1,r-1),
\tag{5.1}
\]

\[
Qk=(-2-D_i^2)_i,\qquad
k^TQk=10-r,
\tag{5.2}
\]

\[
Qp\ge0,\qquad p^TQp=d,
\tag{5.3}
\]

\[
k+3p\ge0,\qquad k+\mathbf1+2p\ge0.
\tag{5.4}
\]

After transporting the surviving dicriticals to the finite Stein model, the
generic cubic ledger adds \(e=2,f=1\), while

\[
K+D=-2\bar F^*L_\infty+\sum e_iE_i
\tag{5.5}
\]

must agree with (5.2), including contraction differents and discrepancies.
The boundary classes must also remain a free basis; no nonzero boundary
combination may become principal.

The correct contradiction certificate is therefore:

> there is no integral weighted tree, pole vector, target grouping,
> contraction matrix, and conductor packet satisfying
> (4.1) and (5.1)--(5.5).

This has not yet been reduced to a finite enumeration.  Optional boundary
blowups create infinitely many presentations.  The missing theorem is a
canonical minimal-model or exposed-leaf reduction which preserves the
finite-normalization decorations.

## 6. Punctures and units

A plane counterexample does **not** presently force an impossible puncture
profile.  In degree three the forced boundary profile is already
\(\mathbb A^1\): one puncture and no nonconstant boundary units.  The cubic
cusp countermodel has exactly that profile.

The useful unit statement is global:

\[
B^\times=k^\times,\qquad
\operatorname{Cl}(B)=\bigoplus_i\mathbb Z[E_i].
\tag{6.1}
\]

For the local cusp cover, \(E=(u+3t^2)\) is principal and deleting it creates
the nonconstant unit \(u+3t^2\).  An actual Keller normalization is forced
to avoid this.  Thus the promising question is:

> Can the cubic cusp packet be globalized while every missing prime remains
> independent in \(\operatorname{Cl}(B)\) and the distinguished complement
> has no nonconstant units?

This is a conductor/class-group problem, not a puncture-count problem.
A useful implementation should build the localization sequence together
with the pullback factorizations

\[
\operatorname{div}_B(g_C)=2E_C+A_C
\tag{6.2}
\]

for every ramified target component, and then impose all intersections and
specializations among the affine companions \(A_C\).

## 7. External prime-degree claim

The 2024 arXiv preprint *There are no Keller maps having prime degree field
extensions* claims, in particular, to exclude geometric degree three.  It
cannot currently be used as a closure theorem for this programme.

Its first case invokes a MathOverflow answer as proving that the “rare
property” forces extension degree two.  The cited answer instead constructs
a quadratic example having the rare property; it does not prove that every
such extension is quadratic.  The implication used in the preprint is
therefore unsupported.  This does not prove the preprint's conclusion
false, but it leaves degree three open for the present boundary programme.

References:

- [the prime-degree preprint](https://arxiv.org/abs/2407.13795);
- [the cited MathOverflow question and answer](https://mathoverflow.net/questions/472877/);
- [Jelonek--Lasoń on parametric nonproperness sets](https://arxiv.org/abs/1411.5011);
- [Nguyen Van Chau on the one-place-at-infinity theorem](https://arxiv.org/abs/math/0305088);
- [Borisov's compactification/Picard frameworks](https://arxiv.org/abs/1901.04073);
- [Stacks Project miracle flatness](https://stacks.math.columbia.edu/tag/00R4).

## 8. Recommended order of attack

1. **Cubic local classification.**  Classify completed normal finite-flat
   cubic algebras over \(k[[u,v]]\) with the generic row
   \((2,1)+(1,1)\), including source singularities and non-Cartier boundary.
   Isolate the cusp packet and its possible degenerations.
2. **Globalize the divisor relation.**  Compile (6.2) for every ramified
   component in the free boundary class group, retaining the affine
   companion and conductor ideals.
3. **Minimal log model.**  Prove a contraction/minimality theorem which
   turns (5.1)--(5.5) into a finite decorated-tree enumeration.
4. **Cusp exclusion.**  Show that every surviving local cusp packet either
   creates a nonconstant unit, a forbidden principal boundary combination,
   or an overlength special fiber after the global target grouping.
5. **Only then generalize in degree.**  For \(d>3\), prove unramifiedness of
   \(E^\nu\to C^\nu\) or account for its Riemann--Hurwitz cost in the global
   residual-different budget.  Do not conflate this with immersion into the
   singular embedded curve \(C\).

The smallest high-value theorem is therefore:

\[
\boxed{\text{No cubic cusp packet can occur in the finite normalization of
a plane Keller map.}}
\]

It simultaneously attacks the degree-three ledger, the conductor escape,
and the only local countermodel to automatic residue immersion.

## 9. Clean nonimmersion consumes the cubic fiber

There is a useful local theorem beyond the countermodel.

### Proposition 9.1 -- clean tangential-defect length

Let

\[
\pi:(X,p)\longrightarrow(Y,q)
\]

be a finite morphism of smooth surface germs in characteristic zero.  Let
\(E\subset X\) be a smooth ramification curve through \(p\), of generic
transverse index two.  Assume:

1. the ramification determinant has exactly the clean factor
   \[
   \operatorname{Jac}(\pi)=z\cdot h,\qquad h(p)\ne0,
   \]
   where \(z=0\) cuts out \(E\); and
2. the restriction \(E\to\pi(E)\) is nonimmersive at \(p\).

Choose a parameter \(t\) on \(E\).  If the first nonconstant target
coordinate restricted to \(E\) has order \(m\ge2\), then

\[
\boxed{
\operatorname{length}_p\pi^{-1}(q)=m+1.
}
\tag{9.1}
\]

#### Proof

The clean Jacobian has order one at \(p\), so \(d\pi_p\) has rank one.
Choose a target coordinate \(u\) whose differential spans its image and use
\((t,u)\) as source coordinates.  Since \(d(u|_E)_p=0\), the smooth curve
\(E\) has an equation

\[
u-a(t)=0,\qquad \operatorname{ord}_t a=m\ge2.
\tag{9.2}
\]

After choosing a second target coordinate \(v\), the Jacobian is

\[
\frac{\partial v}{\partial t}
=(u-a(t))h(t,u),\qquad h(0,0)\ne0.
\tag{9.3}
\]

On the fiber \(u=0\),

\[
\frac{\partial v(t,0)}{\partial t}
=-a(t)h(t,0),
\]

so \(v(t,0)\) has order \(m+1\).  Therefore

\[
k[[t,u]]/(u,v)\simeq k[[t]]/(t^{m+1})
\]

up to multiplication by a unit, proving (9.1).

### Cubic consequence

For a finite-flat cover of rank three, (9.1) gives

\[
m+1\le3.
\]

Nonimmersion has \(m\ge2\), hence necessarily

\[
\boxed{m=2,\qquad
\operatorname{length}_p\pi^{-1}(q)=3.}
\tag{9.4}
\]

Thus a clean nonimmersive boundary point consumes the entire cubic fiber.
There is no affine point above \(q\).  Expanding (9.3) shows that the image
has leading parametrization

\[
u=c\,t^2+\cdots,\qquad v=c'\,t^3+\cdots,
\qquad cc'\ne0,
\tag{9.5}
\]

so it is analytically the ordinary cubic cusp packet.  The model in
Section 2 is not merely an example: it is the forced clean rank-three
leading form.  The length-three fiber in (9.4) is curvilinear.  Its parameter
\(t\) lifts to the finite flat local algebra and generates it by Nakayama;
comparison of ranks therefore gives a monogenic cubic presentation.  After
the usual Tschirnhaus translation, its leading equation is
\(T^3+uT-v\).  This upgrades the leading-form calculation to the standard
analytic cubic packet and determines the companion branch as well.

The affine companion has an additional forced decoration.  In the exact
model (2.6), put

\[
E=(u+3t^2),\qquad A=(4u+3t^2).
\]

Then

\[
\boxed{I_p(E,A)=2.}
\tag{9.6}
\]

Indeed, modulo \(E\), the companion equation is \(-9t^2\).  This contact
equals the normalization-conductor degree of an ordinary cusp.  Thus the
complete clean local cubic packet is

\[
\boxed{
\bigl((2,1)_E,(1,1)_A;\
k[t]/(t^3);\
I_p(E,A)=2;\
\deg\mathfrak c_C=2\bigr).
}
\tag{9.7}
\]

This explains why a naive smooth-target residual-different calculation
appears overdrawn by two: the missing two units are exactly the target
adjunction/conductor correction.  Any global lattice compiler must retain
this contact and cannot apply the smooth-target identity across the cusp.

The executable function
`clean_tangential_defect_budget` and its regression distinguish three cases:

- \(m+1<d\): the defect leaves a residual fiber budget;
- \(m+1=d\): it consumes the full fiber;
- \(m+1>d\): it is impossible by finite flatness.

### Dirty-point alternative

If the residual Jacobian factor \(h\) vanishes at \(p\), then \(p\) lies on
another ramification component or on a non-Cartier/singular-source defect.
Proposition 9.1 does not apply.  This gives a sharp dichotomy for degree
three:

\[
\boxed{
\begin{array}{ll}
\text{clean nonimmersion} &\Rightarrow
  \text{ordinary cusp and a full length-three omitted fiber},\\[1mm]
\text{non-clean nonimmersion} &\Rightarrow
  \text{a residual-ramification or singular-source vertex}.
\end{array}}
\tag{9.8}
\]

The second row is exactly the data seen by the residual-different and
intersection-matrix ledgers.  It cannot be hidden in an undecorated
puncture count.

## 10. Global consequence of the clean cubic packet

Let \(C=(g=0)\) be a ramified nonproperness component with the cubic ledger
(4.1), and suppose all points of its boundary normalization are clean in
the sense of Proposition 9.1.

At a multibranch singular value, two distinct boundary-normalization points
each contribute at least the transverse length two.  Their total is at least
four, contradicting rank-three flatness.  Hence every singularity of \(C\)
is unibranch.

At each such singularity, Proposition 9.1 consumes the full fiber.
Consequently the Keller open has no point over any singular point of \(C\).
The polynomial

\[
H=g(P,Q)\in k[x,y]
\tag{10.1}
\]

then has no critical point on its zero fiber:

\[
dH=(DF)^t(dg),
\tag{10.2}
\]

\(DF\) is invertible everywhere, and on \(C=(g=0)\) the zero locus of
\(dg\) is the singular locus of \(C\), which is omitted by \(F\).  Critical
points of the polynomial \(g\) away from \(C\) may still pull back to
critical points of \(H\), so no global-submersion claim is made.

Moreover, the generic ledger says that the divisor \(H=0\) in the Keller
open has one reduced irreducible component.  Hence \(H\) is irreducible and
\(H^{-1}(0)\) is a smooth affine rational curve.  Its normalization is
\(\mathbb A^1\), with at least one point removed for every cusp collision
with the deleted boundary.

This does not yet give a contradiction: irreducible polynomials can have a
smooth punctured-rational zero fiber, and regular fibers can still be
atypical at infinity.  It does, however, replace the vague puncture target
by a concrete global object:

> **Cubic cusp pullback problem.**  Exclude an irreducible polynomial
> \(H=g(P,Q)\), with zero a regular value, whose zero fiber is the affine
> companion of a deleted cubic ramification cusp and whose divisor closure satisfies
> \(\operatorname{div}_B(g)=2E+A\).

There are now two independent ways to close the clean cubic case:

1. prove that such a polynomial cannot admit the finite-normalization
   divisor relation \(2E+A\) with the free boundary class group (6.1); or
2. compute the topology at infinity of \(H\) from the degree-three monodromy
   and contradict the Suzuki/atypical-fiber formula.

The remaining non-clean case is routed to the lattice attack through
(9.8), rather than being incorrectly declared immersive.

## 11. Complete cubic closed-fiber algebra atlas

Finite freeness makes the scheme-theoretic closed-fiber enumeration exact.
Over an algebraically closed field, a rank-three fiber is an Artin algebra
of total length three.  Its point-length partition is one of

\[
1+1+1,\qquad 2+1,\qquad 3.
\tag{11.1}
\]

At a local point of length three there are exactly two algebra types:

\[
\boxed{
k[\epsilon]/(\epsilon^3)
\quad\text{or}\quad
k[\epsilon,\eta]/(\epsilon,\eta)^2.
}
\tag{11.2}
\]

Indeed, if \(\mathfrak m\) is the maximal ideal, then either
\(\dim\mathfrak m/\mathfrak m^2=1\), giving the curvilinear algebra, or it
equals two; length three then forces \(\mathfrak m^2=0\).

For a target component with the ramified cubic generic row, the reduced
partition \(1+1+1\) is unavailable at its generic point.  Every closed
specialization therefore belongs to:

| fiber | geometric meaning | frontier route |
| --- | --- | --- |
| \(k[\epsilon]/(\epsilon^2)\oplus k\) | simple boundary ramification plus a separate companion | track whether the companion remains affine or moves to another boundary point |
| \(k[\epsilon]/(\epsilon^3)\) | curvilinear triple collision | the clean nonimmersive case is the cusp packet (9.7) |
| \(k[\epsilon,\eta]/(\epsilon,\eta)^2\) | embedding-dimension-two triple collision | necessarily a singular-source vertex, by Proposition 11.1 below |

### Proposition 11.1 -- the square-zero row is a triple surface singularity

Let \(R=\mathcal O_{\bar X,p}\) be the normal surface local ring at a point
whose finite-flat cubic fiber is

\[
R/(u,v)\simeq k[\epsilon,\eta]/(\epsilon,\eta)^2.
\tag{11.3}
\]

Then:

1. \(R\) is not regular;
2. \(R\) is not Gorenstein; and
3. its Hilbert--Samuel multiplicity is three.

#### Proof

The target parameters \(u,v\) are an \(R\)-regular sequence because \(R\)
is Cohen--Macaulay and finite flat over the regular target.  If \(R\) were
regular, or merely Gorenstein, its quotient by this regular sequence would
be an Artin Gorenstein algebra.  The square-zero algebra in (11.3) has
two-dimensional socle, so it is not Gorenstein.  This proves the first two
claims.

The parameter ideal \((u,v)\) has colength three, so its parameter-ideal
multiplicity is three.  The maximal-ideal multiplicity of \(R\) is therefore
at most three.  It is not one because \(R\) is singular.  If it were two,
the multiplicity--embedding-dimension inequality would force embedding
dimension three.  The completed two-dimensional Cohen--Macaulay domain
would then have embedding codimension one and hence be a hypersurface, so
it would be Gorenstein.  The second claim excludes that case, leaving
maximal-ideal multiplicity three.

Thus the square-zero row is not an amorphous closed-point defect.  It is a
normal non-Gorenstein multiplicity-three surface singularity sitting on the
deleted boundary.

Miyanishi's cylinder theorem makes its rationality unconditional here:
every singularity of a normal affine surface containing a cylinderlike open
set is rational.  The finite normalization is normal affine and contains
the distinguished dense open

\[
\mathbb A^2=\mathbb A^1\times\mathbb A^1.
\]

Therefore the square-zero point is a **rational triple point**, whose
minimal resolution graphs are classically classified.  The next finite
lattice computation should:

1. import those multiplicity-three graphs;
2. mark the strict transforms of every \((2,1)\) boundary row and its affine
   companion;
3. glue the local graph into the unimodular \(\mathbb A^2\) boundary tree;
   and
4. reject graphs violating the canonical vector, pole vector, or free
   boundary class group.

References for this reduction are:

- M. Miyanishi,
  [*Singularities of normal affine surfaces containing cylinderlike open
  sets*](https://doi.org/10.1016/0021-8693(81)90264-7);
- Altıntaş--Çevik--Tosun,
  [*Nonisolated forms of rational triple point singularities and their
  resolutions*](https://arxiv.org/abs/1302.1464).

This is a genuine finite decorated enumeration at the closed-fiber-algebra
level.  The remaining infinitude is only in:

- how many target singular values occur;
- conductor/contact orders away from the clean cubic cusp;
- how the separate length-one point specializes between the affine and
  boundary loci; and
- the resolved surface tree carrying the dirty vertices.

Thus degree three has been reduced to two global frontiers:

\[
\boxed{
\begin{array}{ll}
\text{curvilinear branch} &\Rightarrow
  \text{exclude the divisor packet }(2E+A,\ I(E,A)=2),\\[1mm]
\text{square-zero/dirty branch} &\Rightarrow
  \text{exclude an order-two ramification vertex in the log lattice}.
\end{array}}
\tag{11.4}
\]

The atlas is encoded by `cubic_closed_fiber_atlas` in the boundary checker.
