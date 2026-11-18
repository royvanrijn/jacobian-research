# Positive convolution models and the Eisenstein angular threshold

A positive probability measure can avoid a central angular band while satisfying
the full operator bounds used separately at each Hecke index. The construction
below strengthens the scalar [EisMomentBarrier](EISENSTEIN_POWER_SHRINKING_WINDOWS.md)
to a compact-group operator with finite norm-coset tests. It identifies information
that is insufficient by itself to improve the current central exponent.
It does not give an arithmetic counterexample or a multiplicative Hecke family.

## Proposition EisPositiveOperatorBarrier

Write
\[
 S=\mathrm{SU}(2)=\left\{\begin{pmatrix}z&w\\-\bar w&\bar z\end{pmatrix}: |z|^2+|w|^2=1\right\},
 \qquad t(g)=|z|^2-|w|^2,
\]
and let \(T\) be its diagonal torus. For every \(0<h\le1/4\), there is an inverse-invariant probability measure \(\mu_h\), bi-invariant under the normalizer \(N(T)\), such that
\[
 \mu_h\{|t|\le h\}=0,\qquad
 \left\|\int_S\rho(g)\,d\mu_h(g)\right\|\le5h^{3/2}
 \tag{B1}
\]
for every nontrivial irreducible unitary representation \(\rho\) of \(S\). Convolution by \(\mu_h\) is a self-adjoint, positivity-preserving Markov operator and has norm at most \(5h^{3/2}\) on the orthogonal complement of the constants in \(L^2(S)\).

There is a finite-level version. Let \(G\) be any finite group with a surjection \(q:G\to E\) onto a finite elementary abelian 2-group. Fix \(\beta\in E\), and let \(\lambda_\beta\) be the uniform probability measure on the nonempty coset \(D_\beta=q^{-1}(\beta)\). Then
\[
                  \mathfrak m_{h,\beta}=\mu_h\otimes\lambda_\beta
 \tag{B2}
\]
is a nonnegative probability measure, is inverse-invariant and avoids the same central band. Its convolution norm on every subspace of positive archimedean weight is at most \(5h^{3/2}\), independently of \(G,q,\beta\). At archimedean degree zero, its operator has eigenvalue \(\xi(\beta)\) on the norm characters \(\xi\circ q\), \(\xi\in\widehat E\), and is zero on their orthogonal complement. Every finite test \(f:G\to\mathbb C\), including any product of local tests, has exactly its uniform conditional scalar mean
\[
                \int f\,d\lambda_\beta=\frac1{|D_\beta|}\sum_{g\in D_\beta}f(g).
 \tag{B3}
\]
For diagonal order-unit symmetry, require its finite image \(\Gamma\) to lie in \(\ker q\). The twelve units of the discriminant-3 order lie in \(N(T)\), so (B2) is bi-invariant under their diagonal image. These are conditions of the relaxation, not a claim of quaternionic realization.

Equivalently, for any proposed normalized positive-weight envelope \(0<e\le1\), all these properties coexist with an empty central band whenever
\[
                  h\le\min\{1/4,(e/5)^{2/3}\}.
 \tag{B4}
\]

## Proof of the archimedean assertion

Let \(d\sigma=dt/2\) on \([-1,1]\). Push \(\sigma\) by the map which replaces \(0<|t|<h\) by \(t+h\operatorname{sgn}(t)\) and fixes the other points. Call its law \(\nu_h\). This is exactly the measure used in EisMomentBarrier. It is a symmetric probability measure, has no atoms, and gives the closed interval \([-h,h]\) mass zero. Its density with respect to \(dt\) is zero for \(|t|<h\), one for \(h<|t|<2h\), and \(1/2\) elsewhere.

Let \(P_k\) be the normalized Legendre polynomial. For \(k\ge2\) and \(|t|\le1/2\), the elementary bound and derivative recurrence proved in [EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md) give
\[
 |P_k(t)|\le2/\sqrt{k},\qquad
 (1-t^2)P_k'(t)=k(P_{k-1}(t)-tP_k(t)),
\]
whence \(|P_k'(t)|\le6\sqrt{k}\). Here \(1/(1-t^2)\le4/3\), \(2/\sqrt{k-1}\le2\sqrt2/\sqrt k\), and \(2\sqrt2+1<4\). Orthogonality to constants gives
\[
 \begin{split}
 \left|\int P_k\,d\nu_h\right|
 &\le h\min\{6h\sqrt{k},4/\sqrt{k}\}\\
 &\le\sqrt{24}\,h^{3/2}<5h^{3/2}.
 \end{split}
 \tag{B5}
\]
The displaced points and intervening segments lie inside \([-2h,2h]\subset[-1/2,1/2]\). Odd moments, including \(k=1\), vanish by symmetry. There is therefore no low-degree exception to (B5).

For each \(t\), choose \(a_t\in S\) with first coordinate of squared modulus \((1+t)/2\), and let \(m_t\) be the probability measure which is the law of \(u a_t v\), with \(u,v\) independent Haar variables in \(T\). Define
\[
                         \mu_h=\int_{-1}^1m_t\,d\nu_h(t).
\]
Haar measure disintegrates as \(dm_S=\int m_t\,dt/2\). Thus
\(d\mu_h/dm_S\) is respectively 0, 2 and 1 on the central band,
adjacent annulus and remaining region; in particular it is absolutely
continuous with bounded density.

The coordinate \(t\) classifies the double cosets \(T\backslash S/T\), so this definition is independent of the choices. Inversion preserves each \(m_t\). Either left or right multiplication by a Weyl element replaces \(m_t\) by \(m_{-t}\); symmetry of \(\nu_h\) proves bi-normalizer invariance.

All irreducible representations of \(S\) are \(\operatorname{Sym}^m(\mathbb C^2)\). Averaging over \(T\) is the orthogonal projection \(\Pi_T\) onto its fixed vectors. If \(m\) is odd this projection is zero. If \(m=2k\), it has rank one, onto the unit vector proportional to \(X^kY^k\). The matrix-coefficient identity in [EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md), equation (4), gives the exact operator identity
\[
 \int\operatorname{Sym}^{2k}(g)\,d\mu_h(g)
       =\left(\int P_k(t)\,d\nu_h(t)\right)\Pi_T .
 \tag{B6}
\]
Equation (B5) proves (B1); Peter-Weyl decomposition gives the convolution norm assertion. This is a full operator estimate, rather than only an estimate for the one displayed zonal matrix coefficient.

## Proposition EisFiniteCosetRelaxation

Averaging over \(D_\beta\) is a translation by a representative of \(\beta\) after projection onto functions of \(q(g)\). It therefore has precisely the degree-zero description above and has norm one. Since \(\beta^{-1}=\beta\), it is self-adjoint. The tensor product with (B1) gives the level-independent positive-weight bound. The same bound holds on invariant subspaces and after orthogonal compression. Formula (B3) is the definition of the uniform coset measure. If \(\Gamma\subset\ker q\), left and right multiplication by \(\Gamma\) preserve it, giving the asserted unit symmetry.

For the prime-to-index finite groups in [EisJointCoset](EISENSTEIN_BOUNDED_PRIME_FACTORS.md), take the product of projective reduction groups, with \(q\) the joint determinant/reduced-norm squareclass map and \(\beta\) the target norm class. The ramified factor at 3 can be included with its actual norm map. The order units have norm one. Uniform conditional measure has the scalar probabilities of that finite local model; it is not being identified with an actual arithmetic shell. In particular, at an inert odd prime \(p\), a torus normalizer meets either determinant squareclass in \(p+1\) elements, whereas the squareclass has \(p(p^2-1)/2\) elements. The bad density is therefore \(2/[p(p-1)]\), exactly as required. Arbitrary exact-valuation tests at unit-target prime-square levels are covered by (B3) without a new assumption about their probabilities.

For compatible reduction levels, Haar measure on these norm cosets pushes forward to Haar measure on the smaller norm cosets. Thus the construction can be used consistently across the finitely many prime-to-index levels selected by an argument. No assertion at an index-dividing principal level is needed or made.

## Corollary EisSingleIndexInformationBarrier

Normalize the actual raw all-ideal operator by its degree \(\sigma(n)\). Its positive-weight bound is
\[
                  e_n=\tau(n)\sqrt n/\sigma(n)\ge n^{-1/2},
\]
since \(\sigma(n)\le n\tau(n)\). Consequently the model operator
\(\sigma(n)\,C_{\mathfrak m_{h,\beta}}\) satisfies the source's entire positive-weight bound whenever
\[
                        h\le5^{-2/3}n^{-1/3}.
\]
Its scalar norm-character eigenvalues are exactly \(\sigma(n)\xi(\beta)\); its remaining scalar operator is zero. Multiplying the measure by \(12\sigma(n)\) also gives the correct raw shell mass. Thus the separate source bounds, positivity, self-adjointness, unit symmetry and all these scalar coset means are compatible with an empty central band of actual difference width \(nh\). Every fixed exponent below \(2/3\) is within this relaxation for large \(n\).

For a particularly direct cofactor test, take \(n=2\ell\) with \(\ell\equiv1\pmod3\) prime. Then \(n\equiv2\pmod4\), so odd subtraction is absent; the pure condition at 3 is automatic. For any fixed tiny budget \(n^\kappa\), \(\kappa<1\), and sufficiently large \(\ell\), the only eligible target divisor is \(d=1\). All small directive and inert sieve levels are prime to \(n\) and are already covered by (B2)-(B3). The actual target has the perfectly balanced representation \(n=\ell+\ell\). Its construction uses multiplicative structure between indices, which this relaxation does not encode; the example is a test of the restricted analytic inputs, not a claim that these targets are arithmetically difficult. The positivity, norm envelope and finite-coset principal means in this relaxation, by themselves, cannot remove its angular gap. Prime-tail or cofactor arguments would need an additional property that is not captured by these inputs.

This strengthens EisMomentBarrier from scalar Legendre moments to convolution by a nonnegative measure and a finite-level relaxation. It does not strengthen the arithmetic obstruction: the model is continuous, does not obey quaternion-shell integrality, and is not a multiplicative Hecke family. A proof using relations between indices, arithmetic discreteness, or additional estimates for joint cancellation across degrees may escape it. For the moving bump the currently summed angular error remains \(\sqrt{nM/H}\); comparison with main mass \(H\) requires \(H^3>nM\) by a power. The inert annulus tail itself requires only \(H>n^{1/2}\) by a power. No improved central theorem has been proved in this pass.

## The Hecke relations are additional information

The model is not positive semidefinite in the operator-theoretic sense:
positivity here means preservation of nonnegative functions. Indeed
\(\int P_4\,d\nu_h=(105h^5-30h^3)/4<0\); the changes in the second
and fourth moments are \(2h^3\) and \(6h^5\). Its finite norm-character
eigenvalues can also be negative, just as for the Hecke operator.

The omitted prime-power recurrence is a concrete distinction. If this
absolutely continuous measure were used for a prime kernel of degree
\(p+1\), the identity \(T_{p^2}=T_p^2-pI\) would require the measure
\[
 (p+1)^2\mu_h*\mu_h-p\delta_e.
 \tag{B7}
\]
Its mass at the identity is \(-p\), since \(\mu_h*\mu_h\) is absolutely
continuous. It therefore cannot be a nonnegative kernel. This example
makes no claim that the family can satisfy those relations: they are
precisely additional arithmetic structure that the relaxation discards.

## Dependencies and scope

The complete scalar construction, the interior Legendre bound and the
zero-weight matrix-coefficient formula are proved in
[EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md). Its uniform all-ideal
Hecke bound and raw degree normalization supply the comparison with
\(e_n\). The finite groups, norm characters, unit action and coset
probabilities are described in
[EisJointCoset](EISENSTEIN_BOUNDED_PRIME_FACTORS.md).
[Moving windows](EISENSTEIN_MOVING_WINDOWS.md) gives the stated angular
error and the central-versus-endpoint distinction. Classification of the
irreducible representations of \(\mathrm{SU}(2)\), torus averaging and
Peter--Weyl decomposition are used explicitly in (B6); no new spectral
estimate is imported.

For exact classical locators, see Koornwinder,
[*Representations of SU(2) and Jacobi polynomials*, arXiv:1606.08189v1](https://arxiv.org/pdf/1606.08189v1):
Theorem 4.2 and exercise (7.9) give the zero-weight matrix coefficient
`P_k(2|z|^2-1)`, equation (6.11) gives the torus weights, and Theorem 6.4
states Peter--Weyl completeness. Averaging the torus on both sides is
orthogonal projection; its zero-weight space in degree `2k` is one
dimensional, while odd degree has none. Hence (B6) is the standard
spherical Fourier reduction, and its bounds extend to the full
`L^2` operator by orthogonal decomposition. These mechanisms are
classical; the displaced measure, its uniform `5h^(3/2)` estimate and
the failure of the additional Hecke relation still require the
specific arguments above. The source's special-function formulas and
torus calculation were read directly; its stated general Peter--Weyl
theorem remains an imported classical theorem.

The proof is analytic and exact. Internal derivation checks reconstructed
the Legendre constants, the full Fourier operator, inversion, finite
coset eigenvalues and normalization. They do not establish an arithmetic
realization, external review, formal verification or a new exponent for
Eisenstein sums. The source audit found no defect in the existing
uniform Hecke input. A proof using arithmetic discreteness, Hecke
relations between indices or joint cancellation between degrees remains
outside the obstruction proved here.
