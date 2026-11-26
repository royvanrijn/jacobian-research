# Split partial ordinarity and projective p-curvature

For a fixed Hilbert–Blumenthal abelian threefold, the number of nonordinary
components at a split prime can be read from denominators of three fixed
characteristic-zero power series. These are coordinates of the Gauss–Manin
horizontal lift through the Hodge lines.

## Pointwise criterion

Let \(F\) be a totally real field of degree \(d\), and let \(S\) be a fine
Hilbert modular variety at neat level with universal abelian scheme
\(f:\mathcal A\to S\). Fix \(x\in S(K)\) over a number field \(K\) containing
the normal closure of \(F\). Work on the smooth integral model away from
the discriminant, level and any finite set needed to extend the data at \(x\).

Write \(H_i\) for the rank-two summand of \(H^1_{\mathrm{dR}}(\mathcal A/S)\)
for an embedding \(i:F\hookrightarrow K\) and \(\omega_i\subset H_i\) for
its Hodge line. Gauss–Manin preserves \(H_i\) and induces a flat projective
connection. Use the product horizontal foliation on
\(X=\prod_S\mathbf P(H_i)\), where the projective bundles parametrize lines,
and let \(P=(x,(\omega_i(x))_i)\in X(K)\).

**Proposition.** At a good prime \(\mathfrak p\) of \(K\) of characteristic
\(p\) split completely in \(F\), the \(p\)-curvature of the horizontal
foliation at \(P_{\mathfrak p}\), as a linear map from the Frobenius
pullback of the base tangent space to the vertical tangent space, has rank
equal to the \(p\)-rank of \(\mathcal A_{x,\mathfrak p}\). Its \(i\)-th
component vanishes precisely when the \(i\)-th partial Hasse invariant of
\(x_{\mathfrak p}\) vanishes.

**Proof.** In characteristic \(p\), the image \(C_i\subset H_i\) of relative
Frobenius is the conjugate line and is the kernel of Verschiebung on \(H_i\).
The split assumption makes both operators preserve the embedding labels.
The restriction of Verschiebung to \(\omega_i\) gives the \(i\)-th partial
Hasse invariant. Thus its zero condition is \(\omega_i(x_{\mathfrak p})
=C_i(x_{\mathfrak p})\), equivalent to \(p\)-rank zero in that height-two,
dimension-one factor of the \(p\)-divisible group.

In degree one the two-step conjugate filtration and Katz's formula express
the \(p\)-curvature at a Frobenius-pulled-back tangent vector as
\[
 H_i\longrightarrow H_i/C_i
       \longrightarrow C_i\longrightarrow H_i.
\]
The middle arrow is the Frobenius pullback of the \(i\)-th Kodaira–Spencer
map evaluated on this vector, up to sign. The Hilbert Kodaira–Spencer
isomorphism identifies the base tangent space with the sum of
\(\operatorname{Hom}(\omega_i,H_i/\omega_i)\). Consequently the \(p\)-curvature
of the vector bundle takes independent values in the \(d\) lines
\(\operatorname{Hom}(H_i/C_i,C_i)\).

Horizontal lift on the frame bundle descends to the projective bundle and
commutes with the passage from the restricted Lie algebra of the frame
group to its vector fields on the projective fiber. Hence, up to the sign
convention for horizontal lifts, projective \(p\)-curvature on a line \(L\)
is the infinitesimal action \(L\xrightarrow{N}H_i\to H_i/L\) of the
vector-bundle \(p\)-curvature \(N\). A nonzero \(N\) as above is rank-one
nilpotent with unique invariant line \(C_i\). The action therefore vanishes
in the \(i\)-th vertical tangent line at \(P_{\mathfrak p}\) precisely when
\(\omega_i(x_{\mathfrak p})=C_i(x_{\mathfrak p})\). Adding the ranks over the
independent directions gives the result. \(\square\)

Katz proves the comparison in [*Algebraic solutions of differential
equations (p-curvature and the Hodge
filtration)*](https://web.math.princeton.edu/~nmk/old/algsoln.pdf), §2.3 and
Theorem 3.2. See Diamond–Sasaki, [*A Serre weight conjecture for geometric
Hilbert modular forms in characteristic p*](https://arxiv.org/abs/1712.03775),
§5.1 for the partial Hasse construction and §8.2 for Kodaira–Spencer. The
split specialization of Diamond, [*Geometric weight-shifting operators on
Hilbert modular forms in characteristic p*](https://arxiv.org/abs/2011.14128),
Theorem 3.3.1 and Remark 5.2.2, also gives the required Kodaira–Spencer
isomorphism.

## A coefficient test on one fixed curve

Take a smooth algebraic curve germ through \(x\) with parameter \(t\) and
a tangent vector having all \(d\) partial Kodaira–Spencer coordinates
nonzero. An étale chart at \(x\) and a line in a direction outside the \(d\)
coordinate hyperplanes supply one. Fix a frame of the pulled-back \(H_i\)
near \(t=0\) whose first vector spans \(\omega_i(x)\) at \(0\). Write
\(\nabla_\partial=\partial+A_i(t)\), where \(\partial=d/dt\), and let
\(Y_i(t)=(u_i(t),w_i(t))^{\mathsf T}\) be the formal horizontal solution
with \(Y_i(0)=e_1\). Set \(y_i(t)=w_i(t)/u_i(t)\in K[[t]]\). These choices
are independent of the prime.

**Proposition (first nonintegral coefficients).** For all sufficiently
large good \(\mathfrak p\) as above, normalize \(v_{\mathfrak p}(p)=1\).
At an ordinary component the first nonintegral coefficients of
\(Y_i\) and \(y_i\) are both in degree \(p\). At a nonordinary
component the first degrees are \(p^2+p\) for the vector \(Y_i\)
and \(2p^2+p\) for the projective coordinate \(y_i\). Each first
nonintegral coefficient has valuation \(-1\), taking the minimum
of the two valuations for a vector.

Indeed \(A_i(t)\) has integral coefficients outside a finite set, and the
recursion gives \(n![t^n]Y_i\) integral. For the integral truncation
\(Y_i^{<p}\) through degree \(p-1\), horizontality implies
\[
 (\partial+A_i)^pY_i^{<p}(0)=-p![t^p]Y_i.
\]
Here the \(p\)-th derivative in the operator has leading coefficient \(1\),
and every lower derivative of \(Y_i-Y_i^{<p}\) vanishes at zero. Modulo
\(\mathfrak p\), the left side is \(\psi_i(\partial)e_1\) at zero because
\(\partial^p=0\) on the reduction of the formal parameter ring. Since
\(w_i(0)=0\) and \(u_i(0)=1\), the coefficients \([t^p]w_i\) and
\([t^p]y_i\) differ by an integral element. Thus \(p![t^p]y_i\) reduces to
the negative of the second coordinate of \(\psi_i(\partial)e_1\).
The chosen tangent makes the bundle curvature in each component nonzero.
Its projective action at the Hodge point detects precisely the
nonordinary components, and \(v_{\mathfrak p}(p!)=1\).

In particular, precisely one of the three \(p\)-th coefficients is integral
at a split prime with precisely one nonordinary component. The residual
density question becomes a distribution question about those denominators.

For the nonordinary case, suppress the embedding index. At a good
unramified split place, complete the parameter ring over the Witt ring
\(R\) of the residue field and use the Frobenius lift \(t\mapsto t^p\).
First crystalline cohomology gives an integral horizontal matrix
\(\Phi(t)\) from the Frobenius pullback to \(H\). Write \(\sigma\) for
Witt Frobenius and use bars for reduction modulo \(p\). The image of
\(\bar\Phi(t)\) is \(C(t)\), and its kernel is
\(\omega^\sigma(t^p)\), since Frobenius pulls back differentials to zero
and has rank one. For \(Z\) horizontal for the \(\sigma\)-twisted
connection with \(Z(0)=\Phi(0)^{-1}e_1\), horizontality and uniqueness
give
\[
                             Y(t)=\Phi(t)Z(t^p).
\]
On the fiber lattice \(M=H(0)\), write \(F=\Phi(0)\sigma\) for
semilinear Frobenius. The image of its reduction is \(C(0)=\omega(0)\).
Its square vanishes modulo \(p\) and its determinant has valuation one.
Therefore \(F^2M=pM\) and \(G=p^{-1}F^2\) is a semilinear
automorphism preserving \(C(0)\) on reduction since it commutes with
\(F\). Relative Verschiebung gives \(pM\subseteq FM\), so
\(Z(0)\) is integral. Its reduction lies outside \(C^\sigma(0)\):
otherwise applying \(\Phi(0)\) would place the primitive vector \(e_1\)
in \(F^2M=pM\).

Write \(Z_j=[t^j]Z\). The coefficients \(Z_j\) with \(j<p\)
are integral, and \(\overline{pZ_p}\) is nonzero in \(C^\sigma(0)\)
by the same \(p\)-jet calculation used in the ordinary case.
For the ratio calculation below we also need integrality of
\(pZ_j\) through \(j=2p+1\). For an integral frame covector
\(\ell\) of the twisted dual connection, the derivative of
\(\ell(Z)\) is \((D^*\ell)(Z)\), where \(D^*\) is the
dual connection in the parameter direction. On reduction
\((D^*)^p\) is the rank-one square-zero dual curvature, and hence
\((D^*)^{2p}=0\). Evaluating at zero gives
\(j!Z_j\in pM^\sigma\) for \(j\ge2p\). The factorial valuations
for \(j\le2p+1<p^2\) prove the needed bound.

Through its first two terms, \(t^{-p}\overline{pZ(t)}\) is a
horizontal section of \(C^\sigma(t)\) with nonzero initial value.
Indeed dividing by \(t^p\) preserves horizontality on reduction, and
the first derivative is determined by the initial value. The lines
\(C^\sigma(t)\) and \(\omega^\sigma(t)\) agree at zero and have
distinct tangents, since \(C\) is horizontal and the chosen partial
Kodaira–Spencer coordinate is nonzero. After the substitution
\(t\mapsto t^p\), the input section therefore leaves the kernel of
\(\bar\Phi(t)\) to exact order \(p\). Thus
\(\overline{pY(t)}\) first has a nonzero coefficient in degree
\(L=p^2+p\), in \(C(0)\). The entire residue series through degree
\(2p^2+p\) is a section of \(C(t)\).

Write \(c(t)\) for the graph coordinate of \(C(t)\) on reduction.
The integral initial part of \(Y\) projects to zero in \(H/C\)
through degree \(p^2-1\), since there it is computed from integral
\(Z_j\). In degree \(p^2\) the projection is nonzero. To see this,
lift \(\overline{pZ_p}\) as \(\Phi(0)^\sigma b\) modulo \(p\)
with \(\bar b\notin C^{\sigma^2}(0)\). Then
\[
 \Phi(0)Z_p
   =p^{-1}\Phi(0)\Phi(0)^\sigma b+\Phi(0)b_1
\]
for some integral \(b_1\). The second term projects to zero on
reduction, and the first is the matrix of \(G\) applied to \(b\)
in the twice-twisted lattice. It sends \(\bar b\) outside the line.
As the integral truncations \(u_0,w_0\) of the coordinates of \(Y\)
through degree \(L-1\) have \(u_0(0)=1\), the difference
\(c(t)-\overline{w_0/u_0}\) has exact order \(p^2\).

Up to degree \(N=2p^2+p\) write
\((u,w)=(u_0,w_0)+p^{-1}(u_1,w_1)\), with \(u_1,w_1\)
integral of order at least \(L\). As \(2L>N\),
\[
     \frac{w}{u} = \frac{w_0}{u_0}
          +\frac{w_1u_0-u_1w_0}{p u_0^2}
               \pmod{t^{N+1}} .
\]
One has \(\bar w_1=c\bar u_1\) and \(\bar u_1\) has exact
order \(L\) because the first nonzero coefficient of \(pY\) lies
in the initial line. Hence the numerator of the second term first
has a nonzero residue at order \(L+p^2=N\), as asserted.

For the integral crystalline comparison, see Lam–Litt,
[*Algebraicity and integrality of solutions to differential
equations*](https://arxiv.org/abs/2501.13175), §3.2.

The [exact Legendre check](verify_cm_newton_split_projective_jets.py)
computes the horizontal vector and its ratio from the rational
hypergeometric connection at \(\lambda=-1\), where \(j=1728\).
At \(p=7\) the first denominator degrees are \(56\) and \(105\),
and at \(p=11\) they are \(132\) and \(253\). At the ordinary
sample \(p=5\) both first degrees are five.

## Scope of available algebraicity criteria

Pointwise vanishing does not give vanishing of the reduced foliation. The
nonzero rank-one \(p\)-curvature acts nontrivially at every line other than
its kernel. This remains true on the curve near \(x\), where the chosen
tangent has nonzero components. The generic-point \(p\)-curvature conjecture
has a different hypothesis. Bost, [*Algebraic leaves of algebraic foliations
over number fields*](https://numdam.org/articles/10.1007/s10240-001-8191-3/),
Theorem 2.1(i), assumes stability by \(p\)-th powers of the foliation on the
reductions at almost all primes. Herblot's [*Algebraic points on meromorphic
curves*](https://arxiv.org/abs/1204.6336), Lemma 3.3 and Proposition 3.6, use
vanishing of the reduced foliation for the evaluation bounds at the
relevant primes. Knowing the single \(p\)-jet at \(P_{\mathfrak p}\) does
not provide those bounds.

There is a separate obstruction to Bost's Liouville hypothesis. In the
complex uniformization \(\Gamma\backslash\mathfrak H^d\), let
\(z=(z_1,\ldots,z_d)\) lift \(x\). In flat projective coordinates the
\(i\)-th Hodge line is represented by \(z_i\), using the period line
\(\mathbf C(z_i,1)\). Its horizontal leaf is the quotient of
\(\mathfrak H^d\) by the subgroup of \(\Gamma\) stabilizing that line.
On this leaf the function
\[
 w\longmapsto\left|\frac{w_i-z_i}{w_i-\overline{z_i}}\right|^2
\]
descends because the conjugated stabilizer acts by rotations of the disc.
It is bounded and plurisubharmonic, and is nonconstant on every open
neighborhood of the Hodge point. It descends also on the product leaf.

Herblot's Theorem 6.1 additionally assumes a finite-order uniformization on a
smooth affine curve at some place. The same function rules out a complex
uniformization of the lifted curve germ of this kind. Any such holomorphic
map would be tangent to the horizontal foliation by analytic continuation
wherever it maps into the projective bundle. The bounded pullback extends
across the discrete preimage of the boundary and the missing points of the
affine curve, so it would be constant on the compact Riemann surface. The
chosen Kodaira–Spencer direction gives a nonconstant restriction to the
germ. The coefficient criterion supplies no non-Archimedean uniformization.

Lam–Litt's Theorems 4.1.1 and 4.2.1 give a related argument for nonsingular
elliptic Picard–Fuchs systems. Their proof uses primes at which a chosen
initial class projects nontrivially modulo the conjugate line, supplied
by Lemmas 4.2.2 and 4.2.3. Ren's [*K3 Picard–Fuchs equations with initial
conditions in the Hodge line*](https://arxiv.org/abs/2611.31120), §1 and
Theorem 1.1, uses the corresponding argument for the Hodge line of a
family of K3 surfaces. The long integral jets above occur specifically
at primes at which the initial class belongs to the conjugate line. They
do not supply the \(\omega(p)\)-integrality hypothesis of the cited
algebraicity results, which concerns all primes outside a fixed finite set.

André's criterion in Tang, [*Algebraic solutions of differential equations
over the projective line minus three
points*](https://math.berkeley.edu/~ytang/p-curv.pdf), Theorem 3.1.5, requires
the all-coefficient bounds \(\tau(y)=0\), \(\rho(y)<\infty\) and a product
of local uniformization radii greater than one. In degree three, at least
one of the three series above has \(\tau(y_i)>0\). To see this, enlarge \(K\)
to contain the \(\ell\)-torsion of \(\mathcal A_x\) for a prime \(\ell>3\),
and use only large rational primes \(p\) splitting completely in \(K\).
At every place over such \(p\) the relative trace \(a\in\mathcal O_F\) has
\(\lvert N_{F/\mathbb Q}(a)\rvert\le 8p^{3/2}\). Two nonordinary factors
would give \(p^2\mid N_{F/\mathbb Q}(a)\) and therefore \(a=0\) for
\(p>64\). Rational \(\ell\)-torsion gives absolute trace \(6\bmod\ell\),
excluding \(a=0\). Hence at least two of the three coefficients at
degree \(p\) have valuation \(-1\) at each such place.

By Chebotarev and partial summation, the sum of \(\log p\) for rational
primes \(p\le n\) split completely in \(K\) is asymptotic to a positive
multiple of \(n\). The sum over the three series of the truncated
denominator contributions in Tang's definition of \(\tau\) is therefore
bounded below by a positive multiple of \(n\), even after deleting all
primes below any fixed bound. Dividing by \(n\) and taking the prescribed
limits gives \(\sum_i\tau(y_i)>0\), which rules out simultaneous use of
this algebraicity criterion for the three projective coordinates.
