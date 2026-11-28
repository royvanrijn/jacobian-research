# The six-cell \((2,2,2)\) interface gap

**Theorem.** In the Gaussian simplex with six cells and every fixed
\(0<\rho<1\), the \((2,2,2)\) isotypic component has the gap
\[
 \mathfrak q_\rho(h)\ \ge\ \ell\sqrt{1-\rho^2}\,\|h\|^2 .
 \tag{M1}
\]
Here \(\ell=6\nu_{12}(F_{12})\). The estimate holds in any
ambient dimension \(d\ge5\).

Use the positive Young basis, the equivariant cubic field \(W\), and
the equivalent integrated residual of
[the conditional comparison](NONHOOK_CUBIC_CONDITIONAL_OBSTRUCTION.md).
The integrated comparison uses the
two-point derivative determinant bound (RC6) of
[the reflection-chamber comparison](REFLECTION_CHAMBER_DIFFUSION_COMPARISON.md).

## Truncated moments measured from the cutoff

Keep the common cutoff \(u\) fixed and write the four untied means as
\(0>m_1>\cdots>m_4\). Put \(x_i=u-m_i\), so that
\(u<x_1<\cdots<x_4\). Set
\[
 f(x)=\phi(x)/\Phi(x),\qquad
 P(x)=(a(x),b(x))=(x+f(x),\,x(x+f(x))+1).
\]
These are the first and second moments of a normal with mean \(x\)
conditioned to be nonnegative. They are also the moments of the
distance \(u-Y\) for \(Y\) with mean \(u-x\) conditioned to be at
most \(u\). In the notation of the conditional comparison, replacing
the two moment columns for \(Y\) by those for \(u-Y\) transforms
\(-X(HB_{\cdot1}\wedge HB_{\cdot2})\) into
\(X(Ha\wedge Hb)\). For any four planar moment points
\(p_1,\ldots,p_4\), let \(\mathcal U(p_1,\ldots,p_4)\) denote this
last expression with \(a,b\) their two coordinate columns.
Writing \(e_i=p_{i+1}-p_i\) and \(D_{ij}=\det(e_i,e_j)\), direct
substitution in \(X,H\) gives
\[
 \mathcal U(p_1,p_2,p_3,p_4)
  =\bigl(6D_{12}+3D_{13},\ 6D_{13}+4D_{23},\ 6D_{23}\bigr)^\top .
 \tag{M2}
\]
The identity requires no ordering for its algebra.

For \(s<t\), (RC5)–(RC6) with two starting points and the common
hard tilt to the nonnegative half-line give
\[
 0\ \le\ \det(P'(s),P'(t))\ \le\ 2(t-s).
 \tag{M3}
\]
The sign is written for increasing starting points. In fact the
matrix in (RC6) has its second derivative column \(b'/2\).
Thus for any two successive disjoint intervals on the real line,
the determinant of their increments of \(P\) is nonnegative and
at most the corresponding determinant for the curve \((x,x^2)\).
The intervals may share an endpoint. Apply this to the three
successive intervals \([x_i,x_{i+1}]\). For the untruncated moment
points the second coordinate is \(x_i^2+1\), where the constant
makes no contribution. Therefore (M2) shows
\[
                    W(m)-\mathcal W(B)\ge0
 \tag{M4}
\]
coordinatewise, with \(B\) the conditional moment array used in the
integrated residual.

## The summed adjacent contribution

Write \(f_i=f(x_i)\), \(\delta_i=f_i-f_{i+1}>0\),
\(\gamma=P(u)\), and \(d_i=P(x_{i+1})-P(x_i)\). Put
\[
 \alpha_i=\det(\gamma,d_i),\qquad \beta_{ij}=\det(d_i,d_j)\quad(i<j).
\]
In addition to nonnegativity of \(\beta_{ij}\), we use
\[
 0\le \frac{\alpha_1}{\delta_1}
       \le \frac{\alpha_2}{\delta_2}
       \le \frac{\alpha_3}{\delta_3},
 \qquad
 \frac{\beta_{12}}{\delta_2}\le\frac{\beta_{13}}{\delta_3}.
 \tag{M5}
\]
To prove this, reparametrize the curve by \(z=-f(x)\). The elementary
Mills identities proved in (S92) of
[the interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md) give
\[
 -f'(x)=f(x)a(x)\in(0,1),\qquad f''(x)>0,\qquad
 \frac{da}{dz}=\frac{1}{-f'(x)}-1>0 .
\]
So the first coordinate of \(dP/dz\) is increasing. By (M3), its
slope \(b'/a'\) is nondecreasing. The first coordinate of \(\gamma\)
is positive and its slope lies below all the derivative slopes at
\(x\ge u\). Indeed \(P(x)\to(0,0)\) as \(x\to-\infty\) and
\(P(u)=\int_{-\infty}^uP'(s)\,ds\); the slope is a weighted
average with weight \(a'(s)>0\). The endpoint limit follows
directly from the nonnegative-distance density, proportional to
\(e^{x d-d^2/2}\) on \(d\ge0\). For \(x=-L\), its normalizer is
at least \(L^{-1} e^{-1-1/(2L^2)}\), and the \(r\)-th unnormalized
moment is at most \(r!L^{-r-1}\) for \(r=1,2\).
For any vector with positive first coordinate and slope at most the
derivative slope at the start of the compared intervals, its
determinant against \(dP/dz\) is nonnegative and nondecreasing over
those intervals. Both its factors are nonnegative and nondecreasing:
the first coordinate of \(dP/dz\) and the difference of slopes.
Use the vector \(\gamma\) for the first chain in (M5) and the
vector \(d_1\) for the second. The slope of \(d_1\) is a weighted
average over the first interval, so it is at most the derivative
slopes over the second and third intervals. Averaging in \(z\)
proves (M5), since the corresponding interval lengths are
\(\delta_i\).

For adjacent term \(l\), put its replaced point first, retaining
the order of the other three points, and let \(v_2,v_3\) be the
successive increments of those three. The replacement for
\(B^{l\to u}\) in distance coordinates is \((0,0)\). That for
\(B^{(l)}\) is \(P(u)=\gamma\). The first-label matrix in the
positive Young basis is \(T_1=\operatorname{diag}(1,1,-2)\).
By (M2), the difference of these evaluations with the second
one acted on by \(T_1\) is
\[
 \left(6\det(\gamma,v_2)+3\det(\gamma,v_3),\
       6\det(\gamma,v_3),\ 18\det(v_2,v_3)\right)^\top .
 \tag{M6}
\]
For example, moving the first point from \(\gamma\) to \((0,0)\)
adds \(\gamma\) to the first increment in (M2), leaving the other
two increments unchanged. The last entry in (M6) comes from
subtracting \(-2\) times the unchanged entry.

Return the replaced point to slot \(l\) by the permutation
\(g_l=s_{l-1}\cdots s_1\), \(s_i=(i\,\,i+1)\) on the four untied
indices. This changes \(T_1\) into \(T_l\), so the \(l\)-th term
in the integrated regrouping is \(f_l R(g_l)\) times (M6).
For \(l=1,2,3,4\) the input pairs \((v_2,v_3)\) are
\[
 (d_2,d_3),\quad(d_1+d_2,d_3),\quad
 (d_1,d_2+d_3),\quad(d_1,d_2),
\]
respectively. Using the three displayed Young matrices in the
conditional comparison, the sum of the four terms has coordinates
\[
 \begin{split}
 &-3(2\delta_2+\delta_3)\alpha_1+
       6\delta_1\alpha_2+3\delta_1\alpha_3
       +12(f_3+f_4)\beta_{12}+12f_3\beta_{13},\\
 &-6\delta_3\alpha_1-4\delta_3\alpha_2+
       (6\delta_1+4\delta_2)\alpha_3
       +8\delta_3\beta_{12}+(16f_2+8f_3)\beta_{13}
       +16f_2\beta_{23},\\
 &-6\delta_3(\alpha_2+\beta_{12})
       +6\delta_2(\alpha_3+\beta_{13})
       +(18f_1+6f_2)\beta_{23}.
 \end{split}
 \tag{M7}
\]
Each row is nonnegative by (M5). In the first row use
\(\delta_1\alpha_2\ge\delta_2\alpha_1\) and
\(\delta_1\alpha_3\ge\delta_3\alpha_1\). In the second row use those
ratios for both \(\alpha_1,\alpha_2\) against \(\alpha_3\). The two
inequalities against the negative terms in the third row are
\(\delta_2\alpha_3\ge\delta_3\alpha_2\) and
\(\delta_2\beta_{13}\ge\delta_3\beta_{12}\).

The first term in the integrated regrouping is nonnegative by (M4)
and \(u+f(u)>0\). Equations (M5)–(M7) give the same conclusion for its
summed adjacent term. This proves the integrated unit-variance
coordinatewise comparison for every ordered mean.

## Passing to the full block

The shift-invariant polynomial \(W\) is homogeneous of degree three.
The Gaussian square-completion and scaling (S91c)–(S91f) of the
interface reduction apply with the transport \(T_l\) in place of the
exterior-hook transport \(\Lambda^p C_l\). The scalar row is unchanged
and equals \(a_\rho/\rho\) by (S49). Consequently the unit-variance
inequality just proved gives \(\mathsf K_\lambda W\le
\rho^2 a_\rho W\) in the ordered positive Young basis,
where \(\mathsf K_\lambda\) denotes the folded raw interaction.

By (S69) the folded kernel is entrywise nonnegative. Positive diagonal
rescaling from the rational basis to the orthonormal Young basis
preserves this assertion and all the coordinatewise inequalities.
The normalized operator \(T=a_\rho^{-1/2}\mathsf K_\lambda
a_\rho^{-1/2}\) is compact and self-adjoint, and
\(\psi=\sqrt{a_\rho}W\) is square-integrable and positive in the chamber
interior. The supersolution bound \(T\psi\le\rho^2\psi\) implies
\(\|T\|\le\rho^2\). Indeed if \(Th=\eta h\), pair
\(|\eta|\,|h|\le T|h|\) with \(\psi\); compactness then gives the norm
bound. The pointwise bound \(a_\rho/\rho\ge\ell/\sqrt{1-\rho^2}\)
gives (M1) in the active space. Each unused Hermite degree multiplies
the normalized interaction by an additional power of \(\rho\),
which preserves the lower bound in any ambient dimension.


The [exact checker](check_six_cell_222_gap.py) reproduces
the Young coordinate field, transport and coefficient identities
(M2), (M6) and (M7).

~~~sh
python3 research/plurality-stability/check_six_cell_222_gap.py
~~~

The [replay record](six_cell_222_gap_verification.json) retains
the algebraic controls. The Gaussian derivative determinant
comparison, the ratio inequalities and the spectral implication
are written proofs.
