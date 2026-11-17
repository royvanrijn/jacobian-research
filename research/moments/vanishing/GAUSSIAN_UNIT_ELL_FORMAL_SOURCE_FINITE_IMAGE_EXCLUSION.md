# A full formal original source when ell(0) is nonzero

Use one Witt Gaussian pair Z,W with E[ZW]=1 and an independent ordinary
Gaussian Y with E[Y^2]=1. Let

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

where a is constant, all other coefficients are polynomials, and
ell(0)!=0. Put Phi_P[Q]=E[Q exp(tP)].

The ORIGINAL full source extends horizontally after inverting the
Gaussian determinant, without a source-cleanliness hypothesis. If also

    deg ell>=2,
    every root of gcd(ell,c^2) has multiplicity at most two,          (F1)

then the ACTUAL scalar unit cannot have the form

    Phi_P[1]=p(t)exp(psi(t)),
    p in K[t] nonzero, psi in tK[t].                                (F2)

In particular the full actual image is not finite over K[t], and P is
neither pure nor safe. The exclusion holds over every characteristic-zero
coefficient field. The comparison part has exactly the dependencies and
internal assurance of the retained double-cover volume theorem; no
independent replay, formal verification or external review is asserted.

## 1. Completed-z Wick contraction preserves the full original source

Initially work over C and set

    R0=C[[z]][W,Y],
    Rhat=R0((t)).

A t coefficient is a formal power series in z but a POLYNOMIAL in W,Y;
its polynomial degrees may grow with the t index. Define coefficientwise

    Ehat[a(z)W^jY^k]=j![z^j]a(z) E[Y^k].                            (F3)

This is well-defined: each coefficient has finite W,Y degree, and each
summand uses only one finite Taylor jet of its formal z coefficient.
It agrees with original Wick contraction on polynomial coefficients
and is C((t))-linear.

All three Gaussian Stein identities hold on Rhat:

    Ehat[(partial_z-W)A]=0,
    Ehat[(partial_W-z)A]=0,
    Ehat[(partial_Y-Y)A]=0.                                         (F4)

For example
j![z^j]a'=(j+1)![z^(j+1)]a proves the first. For j>=1,
j(j-1)![z^(j-1)]a=j![z^(j-1)]a proves the second; for j=0
both sides vanish because a has NO negative z powers. The last identity
is ordinary polynomial Gaussian integration by parts. Parameter
differentiation commutes with Ehat coefficientwise.

Write d=1-2at and

    q=c^2+2d ell/t=c^2-4aell+2ell/t.

Because ell(0)!=0, ell is a unit of C[[z]], and

    q^-1=t/(2ell)
           *sum_(j>=0)(-t(c^2-4aell)/(2ell))^j                     (F5)

belongs to t C[[z]][[t]]. Therefore K[z,W,Y,q^-1], K=C(t), embeds
coefficientwise into Rhat. Every element has a finite lower t bound.
Also exp(tP) belongs to R0[[t]]. Thus

    L_t(A)=Ehat[A exp(tP)]                                          (F6)

is defined for EVERY element of that localized polynomial ring and
agrees with Phi_P on EVERY original polynomial marker.

Applying (F4) to A exp(tP) and differentiating (F6) proves the full
spatial and parameter identities

    L_t((partial_z-W+tP_z)A)=0,
    L_t((partial_W-z+tP_W)A)=0,
    L_t((partial_Y-Y+tP_Y)A)=0,
    partial_t L_t(A)=L_t(A_t+PA).                                  (F7)

All products and derivatives stay within the stated completion.
Consequently (F6) is an ORIGINAL horizontal functional on the full
localized spatial de Rham quotient. This is a functional extension;
no assertion that spatial localization is a quasi-isomorphism is needed.

## 2. Gaussian fibre elimination and the tracked cover seed

Over K[z,q^-1], the W,Y precision matrix, linear source and determinant are

    M=[[-2t ell,-tc],[-tc,d]],
    b=(th-z,tu)^T,
    Delta=det M=-t^2q.

The determinant is invertible. Polynomial Gaussian fibre reduction
therefore contracts all nonconstant centred W,Y monomials using M^-1;
its top quotient is free of rank one. Completing the square gives
the reduced phase

    phi=tk+(1/2)b^T M^-1 b.                                        (F8)

Differentiating the same Gaussian relations gives the connection
correction -(1/2)tr(M^-1 M_z) or -(1/2)tr(M^-1 M_t). These are exactly
the derivatives of Delta^-1/2; moving-centre terms are spatial exact
terms. This is the usual determinant-factor Gaussian elimination, with
both spatial and parameter connections retained.

Thus the localized original functional (F6) descends to the normalized
half-density system with seed Delta^-1/2 dz. This algebraic reduction is
also the one retained in the
[varying-quadratic-W source construction](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#exact-source-criterion-and-localization-of-the-original-seed).
Its former source-cleanliness requirement is replaced here by the
explicit completion argument (F3)-(F7).

On the normalized double cover

    U^2=q, U!=0,
    eta=dz/(tU),

the odd sector is this half-density module. Delta^-1/2 dz differs from
eta only by a fixed nonzero complex constant, absorbed into the functional
normalization. Extending by zero on the even sector gives a horizontal
functional on the full localized curve cohomology with

    L_curve(eta)=Phi_P[1].                                         (F9)

All ORIGINAL marker values are retained through their Gaussian fibre
reductions. This avoids substituting an untracked stationary cycle or
an alternate all-marker functional.

## 3. The scalar quotient and the complete positive-volume contradiction

Assume (F2), with p(0)=1. Gauge the curve phase to phi-psi and its
functional to exp(-psi)L_curve. The allowed seed eta/p has value one.
Horizontality on the full quotient implies, for every n>=0,

    exp(-psi)L_curve(nabla_t^n(eta/p))=partial_t^n(1),

where nabla_t is the shifted connection. Its restriction to the
K-differential span of eta/p is therefore a nonzero rational horizontal
quotient onto K, not merely a formal scalar identity.

The retained cover-volume estimates in
[finite-centre analysis](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#positive-volume-order-over-finite-centres)
and
[infinity and complete comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#infinity-and-complete-period-comparison)
give positive order for eta on every component of every sufficiently
refined marked model under (F1). These are inequalities for q and eta
ALONE. They survive subtraction of any spatially constant psi, including
all additional marked sections or components that the shifted phase
creates. For c=0 the retained fixed-cover companion gives the same
positive normalization; a geometrically disconnected normalization is
treated componentwise.

Every generic horizontal pole of eta is logarithmic or holomorphic:
fixed common factors have multiplicity at most two, moving simple branch
points are regular in their branch coordinate, and infinity is
logarithmic or holomorphic since deg ell>=2. Dividing by p changes no
spatial pole and adds deg p>=0 to every infinity order.

The shifted phase is spatially nonconstant. Its t^-1 coefficient is

    -z^2/(4ell(z)),

which has a double zero at zero and cannot be constant because ell(0)!=0.
The retained complete compact-piece, finite/finite annulus,
finite/pole descending-tail and gluing comparison therefore makes every
zero-exponential eta/p period tend to zero. Its rational cyclic unit
quotient, obtained above, would have seed value one. The retained
ordinary-parameter extension and perfect period pairing give the
contradiction exactly as in the existing cover proof.

This proves the scalar exclusion (F2). Finiteness of the FULL actual
image would supply (F2) by the
[general finite-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image),
so it too is excluded. No moving-nullsource or other-root cleanliness
condition is used.

Over any characteristic-zero field, a hypothetical (F2) uses finitely
many phase,p,psi coefficients. Embed their finitely generated field in C.
The nonzero ell(0), degree and geometric gcd-multiplicity conditions
survive extension, as do all formal original coefficients. The complex
contradiction descends. The completed-z construction itself is algebraic
over the original coefficient field.

## 4. Exact boundary and bounded independent controls

The hypothesis ell(0)!=0 is essential to this FULL-source completion
argument. If one admits a negative z power, the second identity of (F4)
already fails for A=z^-1:

    Ehat[partial_W A]=0, Ehat[zA]=1.                                (F10)

Thus allowing arbitrary Laurent-z coefficients creates an actual residue
anomaly. This explains why the simple-origin scalar bridge requires a
separate construction. No original all-marker localization is inferred
there from the present proof.

There is also an independent rationalizable scalar check. Take

    ell=(1+z)^2, h=c=a=u=k=0, t=epsilon^2.
    C=sqrt(2)epsilon(1+z),
    zeta=sqrt(2)epsilon S/(1-sqrt(2)epsilon S),
    J=1-sqrt(2)epsilon S.

Ordinary auxiliary Gaussian contraction and formal circular residue give

    Phi_P[1]=E_S[1/J]
      =sum_(n>=0)(2n)!t^n/n!.

Direct ORIGINAL contraction of ((1+Z)^2W^2)^n gives exactly the same
coefficients. In this model U=sqrt(2)/(epsilon J) has a scalar leading
coefficient, so every allowed inverse U power evaluates as a Laurent
epsilon series with polynomial S coefficients. The moving-coordinate
chain rule supplies exact checks for arbitrary displayed forms and scalar
gauges, independently of numerical moment ranks.

The [portable checker](verify_gaussian_unit_ell_formal_source.py) and
[receipt](gaussian_unit_ell_formal_source_verification.json) retain bounded
formal-germ Stein identities, q^-1 identities, localized-source controls,
the rationalizable scalar and covariant chain-rule checks, and the residue
anomaly negative control. They pin this proof, the checker and the retained
cover proof. They do not certify unrestricted comparison theory or infer
all-order conclusions from finite jets. The unrestricted original-source
extension is proved in Sections1-2, and the comparison dependencies and
pending independent-review status are unchanged.
