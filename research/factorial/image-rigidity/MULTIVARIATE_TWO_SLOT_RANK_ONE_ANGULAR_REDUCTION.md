# A rank-one angular reduction for fully shared two-slot profiles

Consider two shared mirrored unit-gap pairs or a single shared
mirrored gap-two pair. At degree \(n\), their normalized array
is a circular coefficient convolution of two Gegenbauer
factors. The [negative force theorem](MULTIVARIATE_TWO_SLOT_MIRRORED_FORCE_COMPARISON.md)
gives exact matrices for the intervening symmetric
shape advances. The two residual factors then give a rank-one
or arrow matrix for the actual upper array.
This reduces the missing original interlacing inequalities
to a secular comparison at the original lower roots.

## Coordinate and force residues

For a centered reciprocal polynomial \(f(\theta)\) of degree
\(n\) in the trigonometric grading, use \(x=\cos\theta\).
If \(n=2p\), write \(f=P(x)\); if \(n=2p+1\),
write \(f=\cos(\theta/2)P(x)\), removing the central root.
In each case \(P\) has degree \(p\) and its increasing
positive-half angles correspond to decreasing roots \(x_j\).
Positive scalar factors do not matter.

Let \(\mathcal L\) denote the operator induced by
\(\partial_\theta^2\) on \(P\). Its two forms are
\[
 \mathcal L_{\rm even}P=(1-x^2)P''-xP',\qquad
 \mathcal L_{\rm odd}P=(1-x^2)P''+(1-2x)P'-P/4 .
\]
These follow by differentiating with \(x=\cos\theta\),
and by dividing out \(\cos(\theta/2)\) in the odd case.
At a root \(\theta_j\) with \(x_j=\cos\theta_j\), define
\[
                 c_j=-\sin\theta_j\,\frac{f''(\theta_j)}
                                           {f'(\theta_j)} .
\]
Then
\[
              \frac{\mathcal LP}{P}
                   =-\frac{n^2}{4}+\sum_j \frac{c_j}{x-x_j}.   \tag{1}
\]
The leading terms give the constant, including in the
odd case, and \((\mathcal LP)(x_j)/P'(x_j)=c_j\)
gives each residue. The negative force theorem gives
\(c_j>0\) for every two-slot symmetric stage in these
shared profiles.
The [Sonin barrier](MULTIVARIATE_TWO_SLOT_SONIN_POLAR_BARRIER.md#two-upper-bounds-on-the-force)
further gives \(c_j<D(1+x_j)/(B+1/2)\) in its notation.

## Native shape advances are negative rank-one changes

A matched symmetric Euler factor in either profile acts
on \(f\), up to a positive scalar, as
\(A^2+\partial_\theta^2\), where \(A=n/2+\lambda\)
for the corresponding current shape \(\lambda>0\).
For a gap-two pair the symmetric row factor raises both
Gauss shapes by \(1/2\); the \(A\) in this formula uses
the shape of the slot at its current offset for that factor.
The same formula follows directly from the matched bases
in the [finite Euler factorization](MULTIVARIATE_MIRRORED_POLAR_EULER_FACTORS.md#factorization-and-the-attached-shapes).
Normalizing the ratio (1) by its leading coefficient, the
new projected polynomial divided by the old one is
\[
                 1+\sum_j \frac{v_j^2}{x-x_j},
                  \qquad v_j^2=\frac{c_j}{A^2-n^2/4}
                               =\frac{c_j}{\lambda(\lambda+n)}.
\]
Thus the new roots are exactly the eigenvalues of
\(\operatorname{diag}(x_j)-vv^{\mathsf T}\), by the
matrix determinant lemma. The statement includes all
native steps and agrees with their strict shape
interlacing.

## The final upper lift

Write \(A_i=n/2+\lambda_{N,i}\) for the two upper shapes,
\(S_A=A_1+A_2\), and
\(Z=(n+\lambda_{N,1})(n+\lambda_{N,2})\).
Let \(P\) now denote the projected advanced degree-\(n\)
array, with roots \(x_j\) and force residues \(c_j\).
The two attached residual factors and the polar lift give
the upper centered function, up to a positive scalar, as
\[
 \cos(\theta/2)
       \left[A_1A_2 f-f''+S_A\tan(\theta/2) f'\right].
\]
This is the same centered rotation used by the
[exact two-residual sign test](MULTIVARIATE_MIRRORED_POLAR_EULER_FACTORS.md#two-decreasing-factors-can-fail-the-sign).

If \(n=2p\), remove the upper factor \(\cos(\theta/2)\).
Since \(\tan(\theta/2) f'=-(1-x)P'\), the resulting
projected upper polynomial is proportional to
\[
                   A_1A_2P-\mathcal L_{\rm even}P-S_A(1-x)P'.
\]
Its leading coefficient relative to \(P\) is \(Z\).
By (1) its monic ratio to \(P\) is
\[
                      1-\sum_j\frac{w_j^2}{x-x_j},
                   \qquad w_j^2=\frac{c_j+S_A(1-x_j)}{Z}.      \tag{2}
\]
The actual upper projected roots are therefore the
eigenvalues of \(\operatorname{diag}(x_j)+ww^{\mathsf T}\).

If \(n=2p+1\), put \(h=\cos(\theta/2)\), so \(f=hP\).
Then \(\tan(\theta/2) f'/h=-(1-x)(P'+P/(2(1+x)))\).
The projected upper polynomial is proportional to
\[
 (1+x)\left[A_1A_2P-\mathcal L_{\rm odd}P-S_A(1-x)P'\right]
                                  -\tfrac12 S_A(1-x)P .
\]
The leading coefficient relative to \(P\) is again \(Z\);
the monic ratio is
\[
                      x-d-\sum_j\frac{w_j^2}{x-x_j},
                 \qquad w_j^2=\frac{(1+x_j)[c_j+S_A(1-x_j)]}{Z}.
                                                               \tag{3}
\]
Here \(d\) is the trace of the upper roots minus the trace
of the advanced roots. Equation (3) is the characteristic
ratio for the arrow matrix with diagonal block
\(\operatorname{diag}(x_j)\), last diagonal \(d\), and
off-diagonal column \(w\).
Thus all entries in (2) and (3), apart from the projected
roots themselves, are exact shape and force expressions.


## Bidiagonal coefficient form

Put \(s=1-x\), \(n=2p+\epsilon\) with \(\epsilon\in\{0,1\}\),
and \(Q(s)=P(1-s)\). The shifted differential operator
\(\mathcal A=\mathcal L+n^2/4\) has the exact action
\[
    {\cal A}s^j=d_js^j+j(2j-1)s^{j-1},\qquad
                     d_j=p(p+\epsilon)-j(j+\epsilon).
\]
Consequently the native symmetric step, normalized by its
leading coefficient, is \(I+\mathcal A/a\) with
\(a=\lambda(\lambda+n)\). If \(Q(s)=\sum_j(-1)^j b_j s^j\)
is scaled so the coefficient magnitudes \(b_j\) are positive,
the inverse step recovers them from the advanced magnitudes
\(\widetilde b_j\) by the positive triangular recursion
\[
         b_j=\frac{a\,\widetilde b_j+(j+1)(2j+1)b_{j+1}}{a+d_j}.
\]
A common positive normalization of the two polynomials
is understood in this identity.

For even \(n\) the upper residual operator from (2) becomes
\(I-[{\cal A}+S_A(p-s\partial_s)]/Z\).
Thus from the fully advanced coefficients \(\widetilde b_j\)
the projected upper magnitudes are
\[
       b^{\rm up}_j=
       \frac{(A_1+j)(A_2+j)\,\widetilde b_j
                  +(j+1)(2j+1)\,\widetilde b_{j+1}}{Z}.
\]
Both arrays have positive coefficient transitions in this
basis when measured from the advanced array. The original
array uses the inverse of the symmetric advances and the
upper array uses one residual step. The secular test still
requires a comparison of their root positions.

## The remaining secular inequality

Write \(l_j\) for the decreasing \(x\)-roots of the
original lower profile and \(x_j\) for the advanced
ones. Native steps give \(x_j<l_j\), while the upper
\(x\)-roots \(y_j\) lie in \((x_j,x_{j-1})\),
with \(x_0=+\infty\).
The required missing ordinary inequality is \(l_j<y_j\).
It follows exactly when \(l_j\in(x_j,x_{j-1})\)
and the secular function in (2) or (3) is negative
at \(l_j\). If \(l_j\ge x_{j-1}\) for \(j>1\), the
desired ordering fails. Both secular functions are
strictly increasing between their poles. This supplies
an exact test at each lower rank.

For example, the top root decreases in one symmetric
advance by at most \(v_1^2\), while the top upper root in
the even case increases from the advanced one by at
least \(w_1^2\), directly from the secular formulas.
These are one-rank estimates. The displacements for
different stages have different rank-one directions.

There is also a trace comparison that holds for every
fully shared ordinary mirrored profile, independently of the number of
slots. If \(c_{m,1}\) is its normalized first coefficient,
the sums of the positive-half \(x\)-roots are
\(-c_{n,1}/2\) at even degree \(n\) and
\((1-c_{N,1})/2\) at odd degree \(N=n+1\).
These follow from the next-to-leading coefficients of the
Chebyshev projections, using \(U_p-U_{p-1}\) as the leading
odd-degree section after removing \(\cos(\theta/2)\).
The [exact first-factor derivative](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
gives \(c_{N,1}/N\le c_{n,1}/n\) for a fully shared profile,
and \(c_{n,1}<n\) by the strict factorial chord. Thus
the upper trace minus the lower trace is
\((1+c_{n,1}-c_{N,1})/2>0\).
The trace comparison and the one-sided rank inequalities
leave the individual secular signs above to be established.
