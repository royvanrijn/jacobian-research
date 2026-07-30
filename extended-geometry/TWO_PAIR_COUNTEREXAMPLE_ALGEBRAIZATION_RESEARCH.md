# Algebraization tests on the local SIC2C4 deformation branch

## 1. Status and scope

This is an exact local computation at the displayed
bidegree-\((4,4)\) SIC2C4 point.  It does not assert a new component of the
all-moment-zero locus and does not prove formal isolation.  The canonical
local setup, coordinate conventions, and all-order tail identities remain in
[the local-moduli note](TWO_PAIR_COUNTEREXAMPLE_LOCAL_MODULI.md).

The main distinction needed for this experiment is:

- the two affine \(9\)-planes already known over
  \(\mathbb Q(\sqrt{41})\) are **fourth-order compatibility fibers**;
- the complete fifth-order compatibility systems on the six displayed
  components are unit ideals.

Consequently none of these components supplies a point to recurse to orders
\(6,\ldots,12\).

## 2. Directions

Write \(h=(h_0,\ldots,h_4)\) in the five-plane coordinates of the
local-moduli artifact.  The experiment uses

\[
\begin{array}{c|c}
\text{label}&h\\ \hline
\text{generic rational}&(2,-1,3,1,-2)\\
\text{previously documented}&(1,2,3,4,5)\\
\text{pure apolar-odd}&(0,1,0,0,0).
\end{array}
\]

The last vector is purely odd for the local apolar involution.  In these
coordinates the odd eigenspace is
\[
 h_0=0,
\]
and the even line is spanned by
\[
 (-210,0,18,69,64).
\]
The even line is a rank-zero fourth-order chart and was retained only as an
exceptional control; its complete fourth fiber was not decomposed here.

The known exact \(F_{a,b}\) direction needs separate interpretation.  Since
\[
 F_{a,b}=a^2b\,
 \bigl(\operatorname{diag}(t,t^{-1})\cdot F\bigr),
 \qquad t^{-2}=b/a,
\]
its tangent is orbit plus scaling.  It is therefore zero after quotienting by
both orbit and scaling, rather than a nonzero direction of the reduced
five-plane.

## 3. The exact-family control through order twelve

Set \(a=1+s,\ b=1\), and put
\[
\begin{aligned}
A_0&=R+Z,\\
B_0&=2W(R+Z)^2-2R^3-R^2Z,\\
B_1&=4W(R+Z)R-2R^3,\\
B_2&=2WR^2.
\end{aligned}
\]
Direct expansion of (6.1) in the local-moduli note gives
\[
\begin{aligned}
F_{1+s,1}
={}&\frac12A_0B_0\\
&+\frac{s}{2}(RB_0+A_0B_1)\\
&+\frac{s^2}{2}(RB_1+A_0B_2)
+s^3WR^3.
\end{aligned}
\]
Thus the coefficients at orders \(4,\ldots,12\) are exactly zero.  This is
already an exact polynomial reconstruction of degree three, hence also a
rational-function and algebraic-function reconstruction.  Its coefficient
sequence satisfies the eventually zero recurrence
\[
 C_n=0\qquad(n\geq4).
\]
The all-order identity
\(\mathcal E_2(F_{1+s,1}^m)=0\) proves that this is an algebraic family, not
merely a reconstructed formal arc.

## 4. Fourth-order fibers

For each of the three nonzero reduced directions in Section 2, restore the
complete \(11\)-parameter freedom in the second correction that remains after
the cubic equations.  Exact characteristic-zero Gröbner reduction gives the
same shape in all three cases:

\[
\begin{array}{c|c|c|c|c}
h&\operatorname{rank}M_4&\dim&\deg&
\text{discriminant square class}\\ \hline
(2,-1,3,1,-2)&2&9&2&41\\
(1,2,3,4,5)&2&9&2&41\\
(0,1,0,0,0)&2&9&2&41.
\end{array}
\]

Each Gröbner basis consists of one affine-linear equation and one rank-one
quadric.  Hence every displayed fiber splits into two conjugate affine
\(9\)-planes over \(\mathbb Q(\sqrt{41})\).  In particular, the square class
\(41\) is not special to the previously selected direction: it persists in
the generic rational and pure apolar-odd samples.

This is sample evidence, not a proof that the square class is constant on a
dense open subset of the five-plane.

## 5. Complete pointwise fifth-order test

The earlier pointwise test fixed an \(11\)-dimensional kernel in the cubic
correction.  At one exact point on one component of each paired fourth fiber,
the new experiment restores all eleven of these directions.  Their induced
changes in the fourth correction are solved from the first twelve moment
equations.  The fifth compatibility matrix then contains

- \(11\) columns from the restored lower lift, and
- \(13\) tangent columns from the new fourth correction.

The results are

\[
\begin{array}{c|c|c}
h&
\text{frozen coefficient/augmented ranks}&
\text{complete coefficient/augmented ranks}\\ \hline
(2,-1,3,1,-2)&2/3&4/5\\
(1,2,3,4,5)&2/3&4/5\\
(0,1,0,0,0)&2/3&4/5.
\end{array}
\]

All ranks are exact over \(\mathbb Q(\sqrt{41})\).  Field conjugation gives
the same ranks at the conjugate selected point.  Therefore none of these six
selected points has a fifth-order lift, even after the previously omitted
lower-kernel freedom is restored.

The pointwise calculation is strengthened component-wide below.

## 6. Component-wide fifth obstruction

Choose one of the two affine \(9\)-planes and write its coordinates as
\(q_0,\ldots,q_8\).  Restore the complete eleven-dimensional kernel in the
cubic correction.  After eliminating the rank-two image of the thirteen
new tangent corrections, the all-order fifth system consists of \(54\)
equations linear in those eleven lower-kernel parameters.

For each of the three directions, exact reduction over
\[
\mathbb Q(\sqrt{41})(q_0,\ldots,q_8)
\]
gives coefficient and augmented ranks
\[
\boxed{2\quad\hbox{and}\quad3.}
\]
Thus every \(3\)-by-\(3\) minor of the coefficient matrix vanishes
identically.  The determinantal consistency equations are the augmented
\(3\)-by-\(3\) minors.  In all three calculations, the \(52\) minors on the
selected two-column chart have parameter degree zero, and one of them is an
explicit nonzero element of \(\mathbb Q(\sqrt{41})\).  Its product with its
field conjugate is a nonzero rational number.

Therefore the augmented matrix has rank at least three at every point, while
the coefficient matrix has rank at most two at every point.  Hence
\[
\boxed{\text{the complete affine \(9\)-plane has no fifth-order lift}.}
\]
Applying \(\sqrt{41}\mapsto-\sqrt{41}\) proves the same statement for the
conjugate component.  This eliminates all six affine \(9\)-planes over the
three selected directions, not merely the six points used in Section 5.

As an independent regression check, specialization
\(\sqrt{41}\mapsto6502\) modulo \(32003\) and saturation on the selected
two-column pivot chart gives the unit Gröbner basis \((1)\) in every case.
The characteristic-zero constant minor, rather than this modular check, is
the proof.

## 7. A reconstructed factor candidate on \(h_3=h_4=0\)

The component checker accepts custom rational directions.  Nine exact
component-wide calculations on the projective plane \(h_3=h_4=0\) reconstruct
the projective ratio between the rational and \(\sqrt{41}\) parts of the
selected constant augmented minor within a quadratic projective ansatz.

Eight samples give a rank-eight interpolation system with one-dimensional
kernel.  A ninth sample at \((h_0,h_1,h_2)=(1,1,2)\) is reserved from the fit
and agrees exactly.  The reconstructed candidate, up to a nonzero rational
chart factor, is
\[
 \Theta(h_0,h_1,h_2)=h_0A_1+\sqrt{41}\,B_2,
 \tag{7.1}
\]
where
\[
\begin{aligned}
A_1={}&1155847373766150h_0
       -36297985953411000h_1
       +27223928948689200h_2,\\
B_2={}&36391011330354397h_0^2
       +6249230221583086080h_0h_1
       -6305341673135930040h_0h_2\\
&+239871821693141332800h_1^2
-485076246209760340800h_1h_2
+245240740811211768000h_2^2.
\end{aligned}
\tag{7.2}
\]
Equivalently,
\[
 \frac{\Theta_{\mathrm{rat}}}{\Theta_{\sqrt{41}}}
 =\frac{h_0A_1}{B_2}.
 \tag{7.3}
\]
The quadratic form \(B_2\) has rank three and determinant
\[
-17597199715531977537918909696315829680766290528000000.
\tag{7.4}
\]
The determinant of the full quadratic form \(\Theta\) is
\[
\begin{split}
&729245647195217455697430711533144208802994049936000000\\
&\quad
-726873922382233060747233151732018411642022039856000000\sqrt{41}.
\end{split}
\]
Its rational norm is nonzero, as checked exactly in the artifact.  Thus the
candidate exceptional locus of this selected minor is the smooth conic
\[
 h_0A_1+\sqrt{41}B_2=0
\tag{7.5}
\]
over \(\mathbb Q(\sqrt{41})\), not a union of linear directions.  The norm
is the quartic
\[
 h_0^2A_1^2-41B_2^2.
\tag{7.6}
\]

The reconstructed conic has no rational projective points.  Indeed, if
\([h_0:h_1:h_2]\in\mathbb P^2(\mathbb Q)\), then \(\Theta=0\) forces
\(h_0A_1=B_2=0\).  On the line \(h_0=0\), the binary restriction of \(B_2\)
has discriminant
\[
-6408370091478404240150045725416960000,
\]
which is negative.  On the line \(A_1=0\), after clearing denominators, its
binary discriminant is
\[
5139006856829622345071075808901549175510230116674634862952296119705720988800000,
\]
which is positive but not a rational square.

All nine evaluations and the interpolation are exact, but the finite sample
does not itself prove the quadratic degree ansatz or a universal identity on
the displayed plane.  A symbolic degree bound or direct identity check is
still required.  Even after that, vanishing of this one augmented minor on
(7.5) would not imply fifth liftability: the other augmented minors must be
restricted to the conic and checked for a common zero.

## 8. Reconstruction and stabilization conclusions

For the exact \(F_{a,b}\) control, polynomial reconstruction succeeds and the
sequence terminates at order three.

For the three reduced directions, every point of each displayed fourth-order
component is obstructed at order five.  There are therefore no coefficients
at orders \(6,\ldots,12\) to feed into rational, algebraic, periodic, or
holonomic reconstruction.  Attempting such reconstruction from a truncated
four-jet would not describe a formal arc.

The nonzero constant minors prove uniform bounded obstruction on these full
components, so no eventual-rank stabilization argument is needed there.
The remaining local question is whether other directions of the reduced
five-plane possess fourth lifts whose fifth compatibility behaves
differently.

## Reproduction

Run

```bash
.venv/bin/python scripts/research_two_pair_counterexample_algebraization.py
.venv/bin/python scripts/research_two_pair_counterexample_fifth_component.py
.venv/bin/python scripts/analyze_two_pair_counterexample_fifth_factor.py
```

The run uses the complete fourth and fifth all-order beta tails.  On the
reference machine it took about nineteen minutes.  It writes
`artifacts/generated-results/two_pair_counterexample_algebraization_research.json`.
The component-wide command runs all three directions; each subcalculation
takes about five to seven minutes and writes the corresponding
`two_pair_counterexample_fifth_component_research*.json` artifact.
The final command is a fast exact replay from the nine stored slice samples
and writes
`two_pair_counterexample_fifth_factor_plane_research.json`.
