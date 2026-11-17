# Toroidal boundary feasibility before nonlinear Kellerization

## Status and outcome

This note introduces a shared necessary-condition front end for the
inverse-Galois/Keller affine-completion problem.  It is a compiler and an
exact re-encoding of existing divisor calculations, not an affine-completion
theorem.

The front end replaces repeated family-specific ledgers by three objects:

1. a fan or cone complex carrying monomial boundary scales;
2. boundary **colors** recording distinct divisorial residue branches over
   the same ray; and
3. one integer valuation matrix whose columns are masks, target pullbacks,
   derivative divisors, conductor functions, units, and proposed affine
   modifications.

Affine completion is then screened first by exact integer identities, finite
integral feasibility, Smith factors, pole inequalities, and smooth-cone
tests.  A package that passes has status `feasible`, never `realized`.  Its
output contains a typed `nonlinear_residue`: polynomial adjugate division,
residue equations, singularities not visible in the fan, affineness outside
the certified toric category, field recovery, and finite-flat fibres.

The implementation is the optional `toroidal_boundary` block in
[`boundary_package_compiler.py`](../scripts/boundary_package_compiler.py).
The exact regression is
[`verify_toroidal_boundary_feasibility.py`](../scripts/verify_toroidal_boundary_feasibility.py).

## 1. Colored fans rather than flat divisor lists

Let \(N\) be a lattice and let \(\Sigma\) be a finite fan or toroidal cone
complex in \(N_{\mathbb R}\).  A ray \(\rho\) records a monomial scale.  A
height-one boundary valuation need not be determined by \(\rho\): several
strict branches can have the same scale and differ only by their residue on
the exceptional divisor.  The compiler therefore uses a map

\[
 \operatorname{carrier}:\{D_1,\ldots,D_r\}\longrightarrow\Sigma(1).
 \tag{1.1}
\]

The \(D_i\) are the boundary colors.  They retain the affine/boundary
coloring used by the Keller etaleness gate and may also carry a residual
label.  This prevents two common losses of information:

- collapsing tangent branches with the same Newton weight; and
- treating a monomial fan as though it already separated every residue
  component.

The code checks that every ray is primitive, every declared cone is
simplicial, and every required smooth cone has Smith diagonal all ones.  A
string certificate must still identify the cone list with the relevant
geometric fan.  The compiler does not infer a global toroidal embedding from
vectors alone.

## 2. The single valuation matrix

Choose functions \(f_1,\ldots,f_q\) which include every divisor-bearing
object used by the ansatz.  Define

\[
 V=(\nu_{D_i}(f_j))\in M_{r\times q}(\mathbb Z).       \tag{2.1}
\]

One matrix can now carry all of the following.

| former ledger | matrix form |
|---|---|
| mask orders | one mask column |
| target-boundary pullback | one pullback column |
| derivative/Jacobian divisor | one derivative column |
| conductor divisor | one conductor column |
| Cox/unit characters | a selected column block |
| affine modification | an unknown nonnegative combination of columns |
| boundary coloring | row metadata and carrier ray |

A fixed divisor identity is a vector equation

\[
 Vc=b.                                                \tag{2.2}
\]

For a family of modifications, let \(x_1,\ldots,x_s\) be bounded integral
variables and \(c_0,c_1,\ldots,c_s\in\mathbb Z^q\).  The tropical search is

\[
 V\left(c_0+\sum_jx_jc_j\right)=b,                   \tag{2.3}
\]

together with extra integral equalities or inequalities.  The latter can
encode support-function positivity, deletion choices, exponent bounds, or a
finite ansatz.  The checker enumerates the declared finite box exactly and
returns its componentwise minimal models.  An empty box is promoted to an
obstruction only when the package says that infeasibility is conclusive and
names a certificate that the box exhausts the asserted scope.  Thus a
bounded search is not silently turned into a theorem.

## 3. Affine-space screens

For a certified normal affine UFD core \(W\subset U\), let selected columns
of \(V\) be a basis of

\[
 \Gamma(W,\mathcal O_W)^*/k^*.
\]

If the selected rows are the complete codimension-one complement, the
localization sequence gives

\[
 0\longrightarrow \Gamma(U,\mathcal O_U)^*/k^*
 \longrightarrow \mathbb Z^q\mathop{\longrightarrow}^{V}
 \mathbb Z^r\longrightarrow \operatorname{Cl}(U)
 \longrightarrow0.                                  \tag{3.1}
\]

The shared Smith routine therefore computes the unit rank, free class rank,
class-group torsion, and exact unimodularity.  This is the same theorem-bearing
screen used by the normal-core boundary-lattice compiler; the toroidal front
end calls the same implementation rather than maintaining another Smith
calculation.

In the certified normal affine toric category, trivial class group and
constant units are sufficient for affine space.  Outside that category they
remain necessary only.  In particular, they do not detect every singularity,
nonaffine relatively ample failure, positive-genus horizontal divisor, or
entrywise inverse pole.

Regularity is also rowwise.  A named reconstruction function may list the
boundary colors on which poles are allowed.  A pole on an affine color is
always rejected.  A pole on any unlisted boundary color is rejected before
coefficient equations are considered.

## 4. \(A_4\): one identity replaces the three-row ledger

The pure-target lift has rows \(W,K,L\) and columns

\[
 \det D\Phi,\qquad \text{auxiliary mask},\qquad
 \Phi^*\mathcal B.
\]

Its valuation matrix is

\[
 V_{A_4}=
 \begin{pmatrix}
 2&1&3\\
 3&0&3\\
 1&1&2
 \end{pmatrix}.                                      \tag{4.1}
\]

The whole pure-target ledger is the single equality

\[
 V_{A_4}(1,1,-1)^t=0.                                \tag{4.2}
\]

The compiler verifies (4.2) coefficientwise on all three colors.  Passing
does not make the lift Keller.  The emitted nonlinear residue is exactly the
current frontier:

- realize two source-dependent masks with entrywise polynomial inverse;
- recognize source and target as affine spaces rather than boundary chart
  opens; and
- preserve the oriented quartic field and complete regular fibres.

Later \(A_4\) extractions should extend (4.1) by new exceptional rows and
selector columns.  They should not restate the \(W,K,L\) balance.  Class
content, positive-genus horizontal components, and normalized exceptional
branches then become separate matrix or residue gates attached to that one
package.

## 5. \(D_5\): the Newton fan supplies the missing exceptional rows

Recall

\[
 C=a^2-4u,\qquad
 R_\pm=a^2-\frac{3\pm\sqrt5}{2}u,\qquad
 Q=R_+R_-,\qquad
 \Delta=CQ^2.                                        \tag{5.1}
\]

The three parabolas have common weight \((1,2)\).  The smooth subdivision

\[
 (1,0),\ (1,1),\ (1,2),\ (0,1)                      \tag{5.2}
\]

has adjacent cone determinants one.  After the \((1,2)\) extraction the
three strict branches have distinct residues

\[
 \frac{u}{a^2}=\frac14,
 \quad\frac{2}{3+\sqrt5},
 \quad\frac{2}{3-\sqrt5}.                            \tag{5.3}
\]

They are therefore three colors carried by the same parabolic ray.  Include
also the exceptional colors \(E_{11}\) and \(E_{12}\).  With columns

\[
 Q,\quad \Delta,\quad C,\quad R_+,\quad R_-,
\]

the full matrix is

\[
 V_{D_5}=
 \begin{array}{c|ccccc}
   &Q&\Delta&C&R_+&R_-\\ \hline
 C   &0&1&1&0&0\\
 R_+ &1&2&0&1&0\\
 R_- &1&2&0&0&1\\
 E_{11}&2&5&1&1&1\\
 E_{12}&4&10&2&2&2
 \end{array}.                                        \tag{5.4}
\]

The last two rows are forced by the fan: each quadratic has orders \(1,2\)
at weights \((1,1),(1,2)\).  Thus \(Q\) has orders \(2,4\), and
\(\Delta=CQ^2\) has orders \(5,10\).  These exceptional checks no longer
need separate ledgers.

For a target branch order \(m\) and source orders
\((s_C,s_+,s_-)\), the compiler solves

\[
 V_{D_5}
 \begin{pmatrix}1\\-m\\s_C\\s_+\\s_-\end{pmatrix}=0.
 \tag{5.5}
\]

The exact regression for \(0\le m\le4\) has four models and the unique
componentwise minimal model

\[
 (m,s_C,s_+,s_-)=(1,1,1,1).                         \tag{5.6}
\]

This agrees with the all-\(m\) proof

\[
 (s_C,s_+,s_-)=(m,2m-1,2m-1).                       \tag{5.7}
\]

in the canonical \(D_5\) obstruction note.  The finite regression is not
offered as a second proof of (5.7).  Its purpose is to ensure that any new
candidate consumes the same colored fan and automatically checks both
exceptional valuations.

The remaining \(D_5\) work is nonlinear: moving the old graph, clearing all
adjugate entries, smoothing and recognizing the modification spaces, and
proving finite flat rank five with the correct normal closure.

## 6. Davenport: tropical feasibility isolates the new divisor

Over the generic target branch ray \(\Delta=0\), the three source colors are
\(E_3,E_6,J\).  Use the columns

\[
 \pi^*\Delta,\qquad J,\qquad E_3.
\]

The matrix is

\[
 V_{\mathrm{Dav}}=
 \begin{pmatrix}
 1&0&1\\
 1&0&0\\
 2&1&0
 \end{pmatrix},
 \qquad \det V_{\mathrm{Dav}}=1.                    \tag{6.1}
\]

Thus the integral completion problem is feasible and has no Smith
obstruction.  This conclusion is deliberately weaker than affine
completion.  Realizing the third column by the direct chart
\((T,Y)\mapsto(E_3,Y)\) creates the coprime divisor \(L(Y)=0\).  The compiler
therefore emits that divisor as nonlinear residue rather than opening new
unit-lattice subbranches.  A successful modification must absorb \(L\), fill
\(\Delta\), control all three colors, preserve the Gassmann closure, and
prove affine-space recognition.

## 7. \(F_{20}\): corrected cover and generic colors

The primary Lecacheux polynomial has

\[
\begin{aligned}
 P={}&X^5+(t^2d-2s-17/4)X^4
 +(3td+d+13s/2+1)X^3\\
 &-(td+11s/2-8)X^2+(s-6)X+1,
 \qquad d=s^2+4.                                    \tag{7.1}
\end{aligned}
\]

The square on \(t\) in the quartic coefficient is essential.  It is present
in [Lecacheux, Theorem 3.1](https://doi.org/10.4064/aa-86-3-207-216) but is
missing from the transcription in Jensen--Ledet--Yui, Theorem 2.3.6.  The
repository's previous formula inherited that typo.  The corrected exact
factorization is

\[
 \operatorname{Disc}_X(P)=\frac1{256}d^3q^2r^2,       \tag{7.2}
\]

where

\[
\begin{aligned}
 q={}&4s^2t^2+4s^2t+8st+6s-8t-5,\\
 r={}&16s^2t^3+4s^2t^2-76st-16s
       +64t^3+16t^2-164t-199.                       \tag{7.3}
\end{aligned}
\]

Also

\[
 P_t=X^2d(2X^2t+3X-1).                              \tag{7.4}
\]

These identities separate three different phenomena on the normalized root
cover.

- Modulo \(d\),
  \(P=(X-1/4)(X-1-s/2)^4\).  There is one unramified color and one tame
  index-four color, with derivative orders zero and three.
- Modulo \(q\), the generic gcd of \(P\) and \(P_X\) is the single root
  \[
   a_q=\frac{2s^2t+2s^2+3s-4}{2(s-1)}.
  \]
  Here \(P_t(a_q)=0\) as well.  The exact quadratic tangent form has nonzero
  discriminant, so the double root is a transverse crossing of two
  unramified normalization branches, not index-two inertia.  Each crossing
  color has derivative order one; the other three geometric sheets have
  derivative order zero.
- Modulo \(r\), the generic gcd has degree two and \(P_t\) is nonzero at its
  two geometric roots.  Thus there are two tame index-two colors, each of
  derivative order one, and one unramified color.

The compiler records the geometric colors in the following compressed table;
the `number` column says how many identical rows are expanded.

| colors | number | \(\nu(d)\) | \(\nu(q)\) | \(\nu(r)\) | \(\nu(P_X)\) |
|---|---:|---:|---:|---:|---:|
| \(d\)-unramified | 1 | 1 | 0 | 0 | 0 |
| \(d\)-index-four | 1 | 4 | 0 | 0 | 3 |
| \(q\)-crossing branches | 2 | 0 | 1 | 0 | 1 |
| \(q\)-residual sheets | 3 | 0 | 1 | 0 | 0 |
| \(r\)-index-two branches | 2 | 0 | 0 | 2 | 1 |
| \(r\)-unramified | 1 | 0 | 0 | 1 | 0 |

The pullback degrees over each target ray sum to five, while the derivative
columns sum to \((3,2,2)\), recovering (7.2).  This closes the generic
codimension-one handoff and prevents Cox or mask searches from confusing the
\(q\)-conductor crossing with ramification.

The first exceptional divisor can also be completed exactly.  At the node
\((s,t)=(1,-1/2)\), write

\[
 s=1+\epsilon,\qquad t=-\frac12+z\epsilon.
\]

The special root polynomial is

\[
 P_0=(X+1)(X^2-3X+1)^2.                              \tag{7.5}
\]

If \(a^2-3a+1=0\) and \(X=a+\epsilon Y\), the coefficients of
\(\epsilon^0\) and \(\epsilon^1\) vanish.  The coefficient
\(R_{a,z}(Y)\) of \(\epsilon^2\) is quadratic in \(Y\), has leading
coefficient \(5(a+1)\), and has discriminant

\[
 -\frac{25}{4}
 \bigl(796az^2-188az+11a-304z^2+72z-4\bigr),        \tag{7.6}
\]

which is nonzero in
\(\mathbf Q(a,z)\).  Moreover the coefficient of \(\epsilon\) in \(P_X\)
is \(\partial R_{a,z}/\partial Y\).  Thus the two roots \(a\), each with two
geometric slopes, give four unramified colors of derivative order one.  The
simple root \(X=-1\) gives a fifth unramified color of derivative order zero.
Since \(q=\epsilon^2(4z^2+8z-1+O(\epsilon))\), the new rows are

| exceptional colors | number | \(\nu(d)\) | \(\nu(q)\) | \(\nu(r)\) | \(\nu(P_X)\) |
|---|---:|---:|---:|---:|---:|
| node slope sheets | 4 | 0 | 2 | 0 | 1 |
| node simple sheet | 1 | 0 | 2 | 0 | 0 |

Their derivative sum is four, exactly the valuation of \(q^2\) on the
exceptional divisor.  The compiler now contains fifteen \(F_{20}\) rows:
ten generic geometric colors and these five exceptional colors.

The exact base-incidence audit makes the next blowup queue finite.  Put

\[
 h_2=16t^2+24t+13,
 \qquad h_3=8t^3+16t^2+2t-7.                         \tag{7.7}
\]

Then

\[
 \operatorname{Res}_s(d,q)=h_2^2,
 \qquad \operatorname{Res}_s(d,r)=3125h_2,           \tag{7.8}
\]

and the affine elimination polynomial for \((q,r)\) is, up to a unit,

\[
 (t-2)h_2h_3^2.                                    \tag{7.9}
\]

The apparent extra factor \(t\) in the raw resultant is a
leading-coefficient intersection at infinity.  The finite centers are:

- the ordinary node \((s,t)=(1,-1/2)\) of \(q\), whose tangent-cone
  discriminant is \(80\);
- the ramphoid \(A_4\) cusp \((11,-1/2)\) of \(r\): with
  \(v=t+1/2\) and \(w=(s-11)-50v\), treating \(r\) as a quadratic in \(w\)
  gives discriminant \(-2048v^5(2v-5)\);
- two conjugate triple points cut out by
  \(s=4t+3\) and \(h_2=0\), where \(d\) and \(q\) are tangent and \(r\) is
  transverse;
- the transverse \(q\)-\(r\) point \((7/12,2)\); and
- three \(q\)-\(r\) tangencies cut out by \(h_3=0\) and \(s=4t^2-5\).

The ramphoid cusp can likewise be completed without a family-specific
ledger.  Its standard embedded resolution has four exceptional divisors.
Exact substitutions in (7.1) give the following table; the pullback order is
the source valuation, so it already includes the displayed ramification
index.

| divisor | base \(\nu(r)\) | root profile | source \(\nu(r)\) rows | \(\nu(P_X)\) rows |
|---|---:|---|---|---|
| \(E_1\) | 2 | \((5)\) | \(10\) | \(4\) |
| \(E_2\) | 4 | \((5)\) | \(20\) | \(8\) |
| \(E_3\) | 5 | \((1,2,2)\) | \(5,10,10\) | \(2,4,4\) |
| \(E_4\) | 10 | \((1,1,1,1,1)\) | \(10,10,10,10,10\) | \(4,4,4,4,4\) |

For completeness, the Newton residuals are especially small.  With \(z\)
the generic exceptional parameter, the first two are

\[
 Y^5+\frac{25}{2}z,
 \qquad
 Y^5+\frac{25}{2}(z-180).                            \tag{7.10}
\]

On \(E_3\), one Newton segment is linear and the other has nonzero part

\[
 Y^4-50Y^2+500,                                     \tag{7.11}
\]

whose discriminant is (2000000000).  On the final exceptional divisor the
residual quintic is

\[
 Y^5-50zY^3+500z^2Y+\frac{25}{2}z^2,                \tag{7.12}
\]

with discriminant

\[
 \frac{1220703125}{16}z^8(2560z-1)^2.               \tag{7.13}
\]

Thus all residual polynomials are generically separable with precisely the
profiles in the table.  The derivative-order sums are respectively

\[
 (4, 8, 10, 20)=2(2, 4, 5, 10),              \tag{7.14}
\]

as forced by the exponent two of (r) in (7.2).  These ten cusp rows join
the ten generic and five node rows, for twenty-five exact (F_{20}) rows.

The two conjugate triple centers can be resolved by one calculation over
\(\mathbf Q(i)\).  At the center

\[
 (s_0,t_0)=\left(2i,-\frac34+\frac i2\right),
 \qquad
 P_0=(X-1/4)(X-1-i)^4,                              \tag{7.15}
\]

and complex conjugation supplies the second center.  The first exceptional
has base orders \((1,1,1)\) for \((d,q,r)\).  A two-stage Newton extraction
of the fourfold root gives one index-four color of derivative order seven;
the simple root has derivative order zero.  On the second exceptional the
base orders are \((2,2,1)\), and a separable quartic residual gives four
unramified colors of derivative order three plus the simple color.  Across
both centers the compressed rows are therefore

| exceptional colors | number | \(\nu(d)\) | \(\nu(q)\) | \(\nu(r)\) | \(\nu(P_X)\) |
|---|---:|---:|---:|---:|---:|
| triple \(E_1\), index four | 2 | 4 | 4 | 4 | 7 |
| triple \(E_1\), simple | 2 | 1 | 1 | 1 | 0 |
| triple \(E_2\), cluster sheets | 8 | 2 | 2 | 1 | 3 |
| triple \(E_2\), simple | 2 | 2 | 2 | 1 | 0 |

These fourteen rows bring the running total to thirty-nine.

Finally, the three \(q\)-\(r\) tangencies are one cubic orbit, so they also
need only one residue-field calculation.  Put

\[
 8\alpha^3+16\alpha^2+2\alpha-7=0,
 \qquad s_0=4\alpha^2-5.                             \tag{7.16}
\]

At any of the three conjugate centers,

\[
 P_0=(X+2+2\alpha)^3
     (X-8+12\alpha^2+8\alpha)^2.                    \tag{7.17}
\]

On the first exceptional, \(t=\alpha+\epsilon\) and
\(s=s_0+z\epsilon\), so the base orders are \((0,1,1)\).  The triple cluster
has one index-two color and one unramified color; the double cluster has one
index-two color.  The common tangent slope is

\[
 m=-8\alpha^2-4\alpha+10.                            \tag{7.18}
\]

Thus the next chart
\(t=\alpha+\epsilon,\ s=s_0+m\epsilon+z\epsilon^2\)
has base orders \((0,2,2)\).  Its triple and double residuals are separable
of degrees three and two, so the second exceptional has five unramified
colors.  If \(Q_j,R_j\) denote the leading pullbacks of \(q,r\) on either
chart, their residual discriminants satisfy the exact identities

\[
 \operatorname{Disc}(A_j)=C_AQ_j^2R_j,
 \qquad
 \operatorname{Disc}(B_j)=C_BR_j,                   \tag{7.19}
\]

in \(\mathbf Q(\alpha)[z]\), where

\[
 C_A=1072\alpha^2+\frac{43561}{16}\alpha
       +\frac{111205}{64},
 \qquad
 C_B=597\alpha^2+581\alpha-\frac{975}{2}             \tag{7.20}
\]

are nonzero.  The eight-row template at each center is

| exceptional colors | number per center | \(\nu(d)\) | \(\nu(q)\) | \(\nu(r)\) | \(\nu(P_X)\) |
|---|---:|---:|---:|---:|---:|
| \(E_1\), triple-cluster index two | 1 | 0 | 2 | 2 | 2 |
| \(E_1\), triple-cluster unramified | 1 | 0 | 1 | 1 | 1 |
| \(E_1\), double-cluster index two | 1 | 0 | 2 | 2 | 1 |
| \(E_2\), triple-cluster sheets | 3 | 0 | 2 | 2 | 2 |
| \(E_2\), double-cluster sheets | 2 | 0 | 2 | 2 | 1 |

The derivative sums are four on \(E_1\) and eight on \(E_2\), as forced by
\(q^2r^2\).  Repeating the single template over the cubic orbit adds
twenty-four rows.

The resulting \(F_{20}\) package has six primitive rays and sixty-three
exact color rows.  The rays added by the tangency resolutions are
\((1,1,1)\), \((2,2,1)\), and \((0,1,1)\); the compiler checks that all six
maximal cones in the corresponding star subdivisions are smooth.  The two
\(q\)-\(r\) exceptionals map to the same primitive tropical ray because
their \((q,r)\)-orders are proportional, while their distinct root-cover
rows retain the infinitely-near information.

The two geometric crossing slopes over \(q\) can now be globalized exactly.
The normalization of the affine nodal \(q\)-curve is

\[
 t=\frac{y^2-9}{8},
 \qquad
 s=\frac{4(y+2)}{(y+1)(y+3)},                        \tag{7.21}
\]

and the two conductor preimages of the node are \(y^2=5\).  On this
normalization the repeated root is \(X=2/(y+3)\), and the complete root
polynomial factors as

\[
 P=
 \frac{(X(y+3)-2)^2
       (X(y+1)^2+2(y+3))H_y(X)}
      {16(y+1)^2(y+3)^2},                            \tag{7.22}
\]

where

\[
\begin{aligned}
H_y(X)={}&16X^2+Xy^4-14Xy^2-16Xy-3X\\
         &+2y^3+10y^2+14y+6.
\end{aligned}
\]

Moving transversely in \(s\) and expanding
\(X=2/(y+3)+\epsilon Y\), the first nonzero coefficient is quadratic in
\(Y\).  Its discriminant is

\[
 \frac{y^2(y-5)(y^2+4y+5)^2}{(y+3)^5}.              \tag{7.23}
\]

Thus the two slopes do not define two rational sections over
\(\mathbf Q(y)\).  They form one connected quadratic cover

\[
 w^2=\frac{y-5}{y+3},
 \qquad
 y=\frac{5+3w^2}{1-w^2}.                             \tag{7.24}
\]

This cover is rational and is branched only at \(y=5,-3\).  The node
conductor pulls back to

\[
 c(w)=w^4+10w^2+5.                                  \tag{7.25}
\]

After deleting the affine points at infinity and the conductor, the
rational unit lattices have ranks three and four:

\[
\begin{aligned}
U_q&=
 \mathbf Q[y,(y+1)^{-1},(y+3)^{-1},(y^2-5)^{-1}]^\times/
 \mathbf Q^\times,\\
U_C&=
 \mathbf Q[w,(w-1)^{-1},(w+1)^{-1},(w^2+3)^{-1},
 c(w)^{-1}]^\times/\mathbf Q^\times .
\end{aligned}
\]

In the cover basis
\((w-1,w+1,w^2+3,c(w))\), pullback of
\((y+1,y+3,y^2-5)\) has matrix

\[
 M_q=
 \begin{pmatrix}
 -1&-1&-2\\
 -1&-1&-2\\
  1& 0& 0\\
  0& 0& 1
 \end{pmatrix}.                                     \tag{7.26}
\]

Its Smith diagonal is \((1,1,1)\), so the cokernel is free of rank one.
Adjoining the selector \(w-1\), with column \((1,0,0,0)^t\), gives
determinant \(-1\).  By contrast the anti-invariant ratio
\((w-1)/(w+1)\) gives determinant of absolute value two.  The norm identity

\[
 \operatorname{Nm}(w-1)=1-w^2=\frac{8}{y+3}          \tag{7.27}
\]

reduces selector powers modulo descended units to parity.  Thus one
selector, rather than two separately named branches, is the integral
conductor handoff.

The compiler appends four proposed mask columns to the sixty-three finite
rows:

\[
 (\mathrm{mask}_d,\mathrm{mask}_q,\mathrm{mask}_r,
   \mathrm{selector}_{w-1}).
\]

The first three repeat the \(d,q,r\) pullback columns.  The selector has zero
order on every finite row because its divisor lies over the boundary at
infinity.  This defines the first certificate-scoped Cox/mask architecture

\[
 d^a q^b r^c(w-1)^e,
 \qquad a,b,c\geq0,\quad e\in\{0,1\}.                \tag{7.28}
\]

The general
[colored divisor-span theorem](COLORED_DIVISOR_SPAN_OBSTRUCTION.md)
shows before enumeration that the derivative target is outside even the
integer lattice generated by these four columns: the generator matrix has
rank three and adjoining the target raises the rank to four.  It emits six
proportional-row witnesses, one each for the \(d\), \(q\), \(r\), triple
\(E_1\), triple \(E_2\), and \(q\)-\(r\)-tangent color classes.  Thus signs,
quotients, and unbounded exponents cannot repair this architecture.

The exhaustive box \(0\leq a,b,c\leq8\) has
\(9^3\cdot2=1458\) assignments and no models.  In fact the contradiction is
visible in two rows: the unramified \(d\)-color has derivative order zero
and forces \(a=0\), while the index-four \(d\)-color has pullback order four
and derivative order three, demanding \(4a=3\).  Hence no nonnegative
base-factor exponent works, independently of the displayed upper bound.

There are consequently no models in (7.28) on which to test entrywise
inverse-adjugate divisibility or affine-space recognition.  This is an exact
empty-survivor result for the base-factor-plus-selector architecture, not an
\(F_{20}\) affine-completion obstruction.  A continuation must add a
genuinely colored Cox divisor whose orders distinguish the unramified and
index-four \(d\)-colors; only surviving models of that broader matrix should
enter nonlinear adjugate and affine-space tests.

That continuation is now carried out at the combinatorial packet level in
[`F20_COLORED_COX_PACKET_FRONTIER.md`](F20_COLORED_COX_PACKET_FRONTIER.md).
Six primitive packet columns break the six displayed relations but leave the
augmented rank larger.  The complete positive derivative support consists of
sixteen disjoint Galois-stable packets; their indicator columns give a
rank-nineteen saturated matrix and the unique nonnegative derivative model.
They compress to three different-factor Cartier candidates satisfying
\(3D_d+D_q+D_r=\operatorname{div}(P_X)\), with unique exponent vector
\((3,1,1)\).
The \(q\)-conductor order and unit-lattice gates pass.  The residue cocycle of
the individual Cox sections and the global regular Cox algebra are not yet
certified, so this packet survivor does not reach entrywise adjugate or
affine-space testing.  Under the separate factorial-core hypothesis, the
three-column compression would retain one unit and a free class group of
rank fifty-seven; the full packet expansion lowers the latter only to
forty-four.  Both naive slices are therefore excluded from affine space.

This completes the finite base-incidence audit, the exceptional root-color
audit, and the global \(q\)-conductor handoff.  It does **not** prove a
global toroidal compactification or affine-space completion.

## 8. Compiler policy for future branches

A new \(A_4\), \(D_5\), \(F_{20}\), Davenport, Cox-fill, or
boundary-normalization ansatz should begin by extending one of these packages.

1. Add every new toroidal ray and certify its carrier geometry.
2. Split equal-scale residue branches into colors.
3. Add one row for every new boundary valuation and one column for every new
   mask, conductor, derivative, selector, or modification function.
4. State divisor cancellation as matrix identities or integral feasibility.
5. Run smooth-cone, pole, unit, class, and support inequalities.
6. Only after the tropical audit is feasible, form the nonlinear residue
   ideal and attempt coefficient elimination, adjugate division, or fibre
   reconstruction.

The existing family notes remain canonical proof sources for their exact
factorizations and obstruction theorems.  The consolidation is prospective:
new descendants should link to their matrix row or residual gate rather than
copying the same valuation argument into another ledger note.

## 9. Limits

The current compiler does **not** prove any of the following from raw matrix
input:

- that a proposed cone list comes from a global toroidal compactification;
- that a general smooth factorial constant-unit variety is affine space;
- that a relatively positive boundary divisor is ample without the required
  geometric certificate;
- that determinant cancellation clears each inverse-matrix entry;
- that conductor residue classes match after a finite jet truncation;
- that a bounded empty search excludes an unbounded family; or
- that a tropically feasible package admits a polynomial Keller map.

Conductor/contact-loss data already handled by the boundary-package compiler
remain proof-bearing inputs.  Divisorial conductor orders may be added as
columns of (V); the finite conductor quotient and its residue matching stay
in the nonlinear handoff unless separately certified by the conductor-jet
theorem.

## 10. Reproduction

Run

```bash
.venv/bin/python scripts/verify_toroidal_boundary_feasibility.py
```

The checker verifies the four fan/matrix encodings, all
smooth-cone Smith diagonals, the two exceptional \(D_5\) rows, the four
bounded \(D_5\) models and their primitive Pareto representative, the
Davenport determinant-one block, and the corrected \(F_{20}\) discriminant,
subresultant profiles, tangent crossing, ten generic colored rows, five
exact rows over the \(q\)-node exceptional divisor, and ten exact rows over
the four ramphoid-cusp exceptional divisors.  It also certifies the
finite \(F_{20}\) base-incidence atlas: one node, one
ramphoid \(A_4\) cusp, two conjugate triple tangencies, one transverse
intersection, and three further tangencies.  It checks rejection of a
perturbed \(A_4\) identity, rejection of a
synthetic class group \(\mathbb Z/4\), and the rule that an empty integral
box is obstructing only with an exhaustive-scope certificate.  Finally it
recovers the primitive positive support weight from a finite system of
integral inequalities.
