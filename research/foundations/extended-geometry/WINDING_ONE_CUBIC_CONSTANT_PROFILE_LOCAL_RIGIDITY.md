# Local rigidity of constant depressed cubic profiles

Fix bounds \(a,b\ge0\) on the polynomial degrees of \(p(x),q(x)\).
Write \(S_{a,b}\) for their affine coefficient space over \(\mathbb Q\)
and \(L\subset S_{a,b}\) for the linear subspace of constant pairs.
At \(u=1\), choose the branch of \(u^{1/3}\) with value one. Let
\(Y(u,w)\) be the unique formal inverse branch with leading term
\(u^{1/3}w\) in
\[
 Y^3+p\left(\frac{u}{1-u}\right)Y+q\left(\frac{u}{1-u}\right)=w^3u.
\]
Its coefficients belong to the field of formal Laurent series at
\(u=1\). Define the inverse-period equations by
\[
           \operatorname{Res}_{u=1} Y(u,w)\frac{du}{u}=0
\]
coefficientwise in \(w\).

**Theorem ConstantDepressedCubicLocalRigidity.** There is a Zariski-open
subset of \(S_{a,b}\) containing every geometric point \((p_0,q_0)\in L\)
with \(4p_0^3+27q_0^2\ne0\), on which the inverse-period equations
generate the ideal sheaf of \(L\). In particular, an integral algebraic
family of solutions to the inverse-period equations which specializes to
such a constant pair has constant profiles throughout the family. Here
constant profile means independent of \(x\). The constant values may
vary with the family parameter. If \(K=\max(a,b)\ge1\), the first
\(4K-2\) negative \(w\)-coefficients already define the same local
locus on a computable such open subset.

## 1. The coefficient equations

The branch expansion has the form
\(Y=u^{1/3}w+\sum_{j\ge1}Y_j(u)w^{-j}\). Recursion in the monic cubic
gives each \(Y_j(u)\) as a polynomial in the coefficients of \(p,q\),
with rational Laurent series coefficients at \(u=1\). This follows by
dividing at each step by the coefficient \(3u^{2/3}\) of the new
unknown in the leading derivative. Every individual residue uses
only finitely many terms of its local Laurent series. In particular,
the inverse-period equations are polynomial equations over \(\mathbb Q\).
They all vanish identically on \(L\), since for constant \(p,q\) the
branch coefficients are regular at \(u=1\).

## 2. The linearization detects every nonconstant direction

Work first at a complex point \((p_0,q_0)\in L\) satisfying the stated
nonvanishing condition. Put \(P(Y)=Y^3+p_0Y+q_0\) and \(\zeta=w^3\).
For constant \(p_0,q_0\), the inverse branch is \(Y_0(\zeta u)\), where
\(P(Y_0(z))=z\). An infinitesimal variation \((\dot p,\dot q)\) of
the coefficient pair gives, up to an overall sign, the linearized
inverse-period series
\[
 \operatorname{Res}_{u=1}
 \frac{
    \dot p(u/(1-u))Y_0(\zeta u)+\dot q(u/(1-u))
 }{uP'(Y_0(\zeta u))}\,du.                                      \tag{1}
\]
Differentiation of the inverse equation gives the formula. It is also
an analytic germ near infinity in the variable \(w\), since the inverse
branch with the chosen leading term is convergent there, locally
uniformly for \(u\) near one. In (1), the residue is a finite linear
combination of derivatives of the analytic inverse branch. A zero
linearized formal series consequently gives a zero germ. Constant
variations have zero residue.

Assume there is a nonconstant variation in the kernel. Let \(K\ge1\)
be the maximum degree of \(\dot p,\dot q\), and let \(\dot p_K,\dot q_K\)
be the two coefficients of \(x^K\), with a zero coefficient for a
polynomial of smaller degree. Since
\((u/(1-u))^K=(-1)^K(u-1)^{-K}+O((u-1)^{-K+1})\), only these terms
give a derivative of order \(K-1\) in the residue. The highest
derivative at that order is
\[
 \frac{(-1)^K\zeta^{K-1}}{(K-1)!}
 \left( \dot p_K\left(\frac{Y_0^2}{2}\right)^{(K)}(\zeta)
           +\dot q_K Y_0^{(K)}(\zeta)\right).                  \tag{2}
\]
The other terms involve derivatives of orders at most \(K-1\) of
\(Y_0\) and \(Y_0^2/2\) with coefficients regular at finite \(\zeta\).
Indeed \(1/P'(Y_0)=Y_0'\), and differentiating a function of \(\zeta u\)
in \(u\) supplies the factor \(\zeta\). When a derivative falls on
\(1/u\) it contributes only to the lower orders.

When \(p_0\ne0\), the critical points of \(P\) are the two distinct
numbers \(y_\pm\) with \(y_\pm^2=-p_0/3\). The critical values
\(c_\pm=P(y_\pm)\) are distinct and nonzero because
\(4p_0^3+27q_0^2\ne0\).
The covering of the regular values of \(P\) by their preimages is
connected, since its domain is the complex plane with finitely many
points removed. Its monodromy is therefore transitive on the three
inverse branches. For each critical value one can continue the branch
used in (1) first by such monodromy and then to a sector at that value
so that it approaches the critical point. The entire expression (1)
continues with that branch. The derivatives defining it are local in
\(u\) at one, so along any compact path of regular \(\zeta\) values they
can be continued in a common sufficiently small neighborhood of
\(u=1\). A zero germ stays zero under these continuations.

At either selected critical value the continued branch is
\[
 Y_0(\zeta)=y_\pm+\kappa_\pm(\zeta-c_\pm)^{1/2}+O(\zeta-c_\pm),
 \qquad \kappa_\pm\ne0.
\]
The coefficient of \((\zeta-c_\pm)^{1/2-K}\) in (2) is a nonzero
factor times \(\dot p_K y_\pm+\dot q_K\). The factor includes
\(c_\pm^{K-1}\), \(\kappa_\pm\) and the falling factorial of \(1/2\)
of length \(K\). None of the lower derivative terms can contribute
to this exponent. Since the germ is zero, the two distinct values
\(y_+,y_-\) force \(\dot p_K=\dot q_K=0\), contrary to the choice of
\(K\).

If \(p_0=0\), then \(q_0\ne0\) and \(Y_0(\zeta)=(\zeta-q_0)^{1/3}\)
on the chosen branch. At \(\zeta=q_0\), the two derivatives in (2)
have respective exponents \(2/3-K\) and \(1/3-K\) with nonzero
coefficients. Neither exponent occurs among the lower derivative terms,
and the factor \(\zeta^{K-1}\) is nonzero there. Their coefficients
force \(\dot p_K=\dot q_K=0\), again a contradiction.
Each geometric point over a field of characteristic zero is defined
over a finitely generated subfield that embeds in \(\mathbb C\), so the
rank statement holds in characteristic zero.

For the explicit prefix, view (1) as a rational function in
\(Y_0(\zeta)\), using \(\zeta=P(Y_0)\) and the derivation
\(D=(P')^{-1}d/dY_0\). The residue terms are linear combinations of
\(P^jD^{j+1}Y_0\) and \(P^jD^{j+1}(Y_0^2/2)\) for \(0\le j<K\).
The two derivatives of order \(j+1\) have denominator
\((P')^{2j+1}\) and numerators of degrees at most \(j\) and \(j+1\),
respectively. This follows by induction since applying \(D\) to
\(A/(P')^m\) increases the denominator exponent by two and the
numerator degree by at most one. With common denominator
\((P')^{2K-1}\), multiplication by \(P^j\) therefore gives numerators of
degrees at most \(4K-4\) and \(4K-3\). The denominator has degree
\(4K-2\). As \(Y_0(w^3)\sim w\), any nonzero such rational function has
some nonzero coefficient among \(w^{-1},\ldots,w^{-(4K-2)}\). The
injectivity thus already holds for that prefix of the equations.

## 3. From the linearization to a Zariski-open statement

Let \(I\) be the ideal generated by the inverse-period equations and
let \(J\) be the ideal of \(L\). We have \(I\subset J\). At such a
constant point, the images of the inverse-period equations span
\(J/\mathfrak m J\) in the local coefficient ring with maximal ideal
\(\mathfrak m\). This is the dual of the normal-direction injectivity
in §2. Equivalently a finite selection of the equations spans that
space, since it is finite-dimensional. Nakayama's lemma gives
\(I_{\mathfrak m}=J_{\mathfrak m}\). This argument can be made after
passing to the residue-field extension of a geometric point, since
equality of the two localized ideals descends by faithful flatness.

One can compute an explicit open subset using the finite prefix.
Write the column of its equations as \(H\boldsymbol n\), with
\(\boldsymbol n\) the column of nonconstant coefficients of \(p,q\) and
\(H\) a polynomial matrix. Such a decomposition is obtained by
assigning each term to a nonconstant coefficient dividing it. At
\(\boldsymbol n=0\), \(H\) is the normal Jacobian of the prefix. The
union of the open sets where a maximal column minor of \(H\) is
nonzero contains every constant profile in the theorem by §2. On
each such open the corresponding linear system expresses all
entries of \(\boldsymbol n\) in the ideal of the prefix. When there
are no nonconstant coefficients, take the whole space. Thus on the
indicated open the finite prefix, the full equations and \(J\)
generate the same ideal. If an integral parameter scheme maps into
the inverse-period locus and meets this open subset, the preimage
of the open is dense. The family then maps into \(L\) on a dense open
and hence everywhere. Solutions supported on the closed complement
remain to be classified.

## 4. Application to the quadruple-root septic branch

**Corollary DegreeFourSepticConstantProfileExclusion.** Let a regular
winding-one cubic phase on \(2xy+t^2=1\) of least sphere degree seven
have \(E=e(1+\lambda x)^4\), \(e\lambda\ne0\), and satisfy the
[boundary divisibilities](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root).
Rescale to \(U=1+x\) and \(e=1\). Its depressed profiles have
\(\deg p\le6\) and \(\deg q\le9\). If its pure moments vanish
eventually, its profile belongs to the closed complement in \(S_{6,9}\)
of the open subset constructed in §3 from the first 34 negative
inverse coefficients. This open contains every constant profile with
\(p_0(4p_0^3+27q_0^2)\ne0\).

To prove the corollary, write \(v=Ut\). Then
\[
 f=\frac{U}{x}\bigl(v^3+d(x)v^2+b(x)v+c(x)\bigr)
   =\frac{U}{x}(Y^3+p(x)Y+q(x)),\qquad Y=v+\frac{d(x)}3.
\]
Here the rational coordinate \(u=x/U\) satisfies
\(dx/(xU)=du/u\). On the formal inverse branch for \(f=w^3\),
\(t\,dx/x=(Y-d/3)du/u\). The term \(d/3\) is independent of \(w\),
and the leading \(w\)-coefficient on the \(Y\) side has zero residue.
Thus the [prime-height inverse-period condition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md)
for an eventually pure phase gives exactly the coefficient equations
in this theorem.

The [least sphere degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md)
gives \(\deg d\le3\), \(\deg b\le5\), \(\deg c\le7\) at degree
at most seven. Hence the depressed profiles lie in \(S_{6,9}\).
If \(p,q\) are constant, expansion of
\((v+d/3)^3+p(v+d/3)+q\) and the same degree formula give least
sphere degree \(\max(6,3\deg d)\), with degree six when \(d=0\).
This cannot equal seven. Consequently every eventually pure septic
in this branch would have its depressed profile in the closed
complement of the open set in the theorem. In particular a family
of such profiles satisfying the inverse-period equations cannot
specialize to a constant pair with \(4p_0^3+27q_0^2\ne0\).
This leaves profiles outside that open set to be decided.
