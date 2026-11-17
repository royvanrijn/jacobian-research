# Filtered Gaussian phases: critical action, exact generic modules and all-ring gates

16 November2026. Internal all-degree/all-ring proof; independent, formal
and external assurance outstanding; no priority claim. The
[producer/checker](verify_filtered_gaussian_critical_action.py),
[certificate](filtered_gaussian_critical_action_certificate.json) and
[source audit](filtered_gaussian_critical_action_source_audit.json) retain
finite identities. The universal statements use the arguments below.
All earlier sealed evidence is preserved.

## 1. The larger architecture

Keep a nonsingular rational Gaussian form Q0 homogeneous of weight2
for supplied NONNEGATIVE rational weights. Now allow P to have weight
AT MOST2, including lower-weight and constant terms. This strictly
enlarges the [homogeneous weighted theorem](SEMIPOSITIVE_WEIGHTED_GAUSSIAN_REDUCTION_AND_SAFE_NORMAL_FORMS.md).
It does not allow arbitrary phases without the supplied weight condition.
As before the zero-weight positions and weight2 momenta split off:

    Q0=Z dot W+(1/2)T^T B0 T,
    P=W dot H(Z)+F(Z,T), 0<w(T_i)<2.                    (1)

F now may contain T-constant and T-linear terms. Coefficients in the
zero-weight Z are unrestricted polynomials. All ring constructions hold
over EVERY commutative Q-algebra R.

## 2. Finite rational critical points and a polynomial volume form

Put v=B0^-1 grad_T F. Its i-th component has weighted degree at most w_i.
Group equal positive weights in increasing order. In a group I one has

    v_I=M_I(Z) T_I+R_I(Z,T_(lower)),                    (2)

where M_I is independent of T and R_I uses STRICTLY LOWER weights.
A nonlinear monomial cannot contain a variable of weight at least w_i.
Thus the critical point p=t v(Z,p), p(0)=0, is obtained by the FINITE
block recursion

    p_I=t(I-tM_I)^-1 R_I(Z,p_(lower)).                  (3)

Define the polynomial

    sigma(t,Z)=product_I det(I-tM_I), sigma(0,Z)=1.

It depends ONLY on the weight2 quadratic part of F. Lower-weight
quadratic terms lie strictly below the diagonal in (2); all nonlinear
terms do too, even after substitution of p. Consequently for

    S=(1/2)T^T B0 T-tF,
    B_p=Hess_T S(Z,p),

one has det(B_p)/det(B0)=sigma EXACTLY. Both p and B_p^-1 have
coefficients in R[t,Z,sigma^-1]. This determinant identity is essential;
it is not an assumption that a general moving Hessian is polynomial.

Translate T by p, then apply the ordinary-degree contractions of the
previous weighted proof with B_p in place of B. The inverse B_p is
filtered: a nonzero entry (i,j) requires w_i+w_j>=2. Hence the correction
-(1/k)B_p^-1 grad R_k has weight at most w_i in coordinate i. Since it
has ordinary degree>=2, it uses strictly lower-weight variables. The
same finite degree ceiling floor(2/min w_i) terminates the process.
We obtain a FINITE polynomial triangular automorphism Theta, including
the translation p, with polynomial inverse over R[t,Z,sigma^-1], satisfying

    det(partial Theta/partial T)=1, Theta(0,Z,T)=T,
    S(Z,Theta)= (1/2)T^T B_p T-psi(t,Z),
    psi=tF(Z,p)-(1/2)p^T B0 p.                          (4)

All identities are in the ORIGINAL potential and ORIGINAL volume.
psi is rational with denominator a power of sigma and psi(0,Z)=0.
No infinite analytic Morse map is substituted for this finite certificate.
The classical Morse/splitting mechanism and Greuel--Pfister's published
2026 theorem retain the [previous exact credit/scope audit](semipositive_gaussian_reduction_source_audit.json).
The finite filtered block construction above supplies its own ring proof.

## 3. Exact full generic image, with its exponential factor retained

Let g=tH(g), J=det(I-tH_Z(g)), L=K(t,g) over a characteristic-zero field.
Evaluate sigma and psi at Z=g and write

    a=J^-1 sigma_g^(-1/2), f=e^(psi_g) a.                (5)

All branches have constant coefficient1 or0 as appropriate. Formal
circular residues and the Jacobian1 change(4) give Phi[1]=f and
Phi[Q(Z)]=fQ(g). Every original residual amplitude is evaluated by Wick
contraction of Q(g,Theta) with covariance B_p^-1. Every W-momentum
amplitude is a finite jet in the normal source eta of
g_eta=eta+tH(g_eta). These jets require only J^-1, sigma_g^-1 and their
finite derivatives, multiplying f. Therefore the FULL actual generic
module, not just the seed or an annihilating source, is EXACTLY

    K(t)M=fL, rank_actual=[L:K(t)].                      (6)

Upper inclusion retains all amplitudes; lower inclusion is supplied by
the ORIGINAL Z monomials. f can be transcendental; multiplication by it
does not change vector-space dimension. Neither the algebraic degree
of the gauge a nor the degree of a nonminimal input polynomial is an
actual rank formula. No actual-image-kernel computation is performed.

After a finite algebraic cover the connection in(6) is a direct sum of
rank-one algebraic gauges times exponentials of the conjugates of psi.
Its geometric differential Galois identity component is therefore a
TORUS, with no additive unipotent part. This assertion does not say the
finite permutation quotient is solvable, or that the global polynomial
lattice is semisimple.

There is an exact differential-type criterion on the WHOLE projective
time line:

    actual generic connection regular everywhere
      iff psi_g=0
      iff every original marked period is algebraic.    (7)

For necessity, a nonconstant algebraic function on its smooth proper
curve has a pole. At a pole of order m>0 the rank-one exponential has
slope m in its local parameter. Algebraic gauges contribute only regular
singular factors. Restriction of scalars retains this irregular summand.
Since psi_g(0)=0, a pole-free constant is exactly0. Sufficiency is the
finite algebraic-cover argument from the previous theorem.

More precisely, at a geometric base place b, let v run over the places
of L above b, with ramification e_v and m_v=max(0,-ord_v psi_g). The
actual slope is m_v/e_v with multiplicity e_v for each v; zero values
are regular parts. Its Newton irregularity at b is sum_v m_v. In a local
coordinate t-b=s^(e_v), this is just the rank-one exponential pole
calculation and the chain rule. The d=[L:K(t)] ranks in(6) ensure no
selected-image kernel has silently removed a branch. These valuations
give an actual slope algorithm without a nonminimal scalar equation.

Thus regular TRANSCENDENTAL images cannot occur anywhere in this larger
filtered architecture. They remain possible only outside it, if at all.
Gaussian moment-preserving polynomial automorphisms with inverse preserve
the full actual image, so the same exclusion holds for its full Gaussian
symmetry orbit. Ordinary-degree cutoff data need not be unchanged there.

## 4. Field purity, all-ring purity, and the exact pure safe normal form

Over a FIELD, purity has the exact split condition

    P pure-null iff psi_g=0 AND J^2 sigma_g=1.           (8)

Indeed f=1 makes exp(psi_g)=J sqrt(sigma_g) algebraic. On a finite
proper curve dlog of an algebraic function has at most simple poles,
whereas dpsi has a pole of order>=2 at any pole of nonconstant psi.
Their equality forces psi constant, hence0. The remaining factor gives
the determinant identity. This algebraic argument works in characteristic
zero without a chosen complex embedding or a contour assertion.

Over arbitrary coefficient rings the exact identity is instead

    exp(2psi_g)=J^2 sigma_g.                             (9)

Its normalized square-root branch makes it equivalent to f=1. One MUST
NOT impose the two field equations separately on a nonreduced ring.

UNDER PURITY, all-original-marker safety over EVERY R is again equivalent
to polynomiality of the finitely many FIXED LINEAR periods g_i:

    safe iff g in R[t]^n.                               (10)

Necessity follows from Phi[Z_i]=g_i. For sufficiency g=p polynomial
makes J and sigma_g polynomials. Reduce (9) modulo each prime and pass
to its fraction field. The field result gives psi=0 and J^2sigma=1;
polynomial units and constant coefficient1 then force J=sigma=1 in
each field. Thus the finitely many nonconstant coefficients of J and
sigma_g are nilpotent in R. Their finite ideal is nilpotent, so J and
sigma_g are polynomial UNITS, without a Noetherian hypothesis.
Every critical-point/Morse/source-jet coefficient evaluated at p is now
polynomial. The common seed is1. The full-amplitude construction proves
safety, including all zero divisors.

The exact finite-polynomial pure-safe normal form is equivalently:
there is p in tR[t]^n with p=tH(p), J_p and sigma_p polynomial units,
and

    psi(t,p)=log J_p+(1/2)log sigma_p.                  (11)

These logarithms are FINITE polynomials on their nilpotent nonconstant
coefficient ideals. The action is polynomial because sigma_p is a
polynomial unit. Under (9), field reduction also makes its finitely many
coefficients nilpotent, so formal log/exp are legitimate finite inverse
operations. Conversely(11) gives purity and the full polynomial image.
Over fields(11) reduces to p polynomial,J_p=sigma_p=1 and psi(t,p)=0.

For pure inputs the exact scalar polynomial-Hom normalization ideal is

    N_P={c in R: c g_i in R[t] for every i}.             (12)

Apply(10) over R/Ann(c). This gate retains every original marked tail;
nonzero normalization alone does not certify safety. The condition is
vacuous when there are no zero-weight coordinates: every pure member
of that strictly positive filtered class is safe over every R.

### A finite ordinary residual ideal for the entire pure architecture

The coordinate criterion has an explicit FIELD degree bound and a finite
all-ring replacement. Set D_i=max(0,deg H_i) and

    B=product_i(D_i+1),
    q=g mod t^(B+1),
    I_B=(all t-coefficients of q-tH(q)) in R.           (12a)

Over a field any POLYNOMIAL g has degree at most B. Its graph is a
one-dimensional component of z_i-tH_i(z)=0: the z-Jacobian at(0,0) is
identity, so it has the expected codimension there. The classical affine
Bezout component inequality bounds its degree by product_i(D_i+1).
For completeness, the first i fixed-point equations have independent
z-Jacobian rows at(0,0), so their unique local component there has
codimension i. Their closures give a chain of components containing the
graph, with dimension dropping at every step. At each proper hypersurface
intersection degree is multiplied
by at most that equation's degree. This last assertion is the elementary
projective Hilbert-series exact sequence A(-d) --f-> A, whose quotient
has degree d*deg A; component multiplicities are nonnegative. The degree
of a polynomial graph is max(1,deg g_i), as a generic hyperplane pulls
back to a polynomial of that degree on the parameter line. This proves
the stated conservative bound; it is not a least-degree assertion.

More generally, for ANY specified truncation degree N, the ordinary
residual ideal of q_N-tH(q_N) is EXACTLY the ideal of all coefficients of
g above degree N. Modulo the tail ideal the truncated vector is already
an exact solution; conversely modulo the residual ideal formal uniqueness
makes g=q_N. These quotient arguments prove both ordinary ideal inclusions
over arbitrary R, retaining all nilpotents.

Under the FULL PURITY premise, the exact finite all-ring criterion is

    P safe iff I_B has NILPOTENT image,
    N_P=Gamma_(I_B)(R).                                (12b)

Necessity follows by prime reduction, the field degree bound and finite
generation of I_B. Conversely let I_B be nilpotent. Modulo I_B, q is an
exact POLYNOMIAL fixed point and purity still holds. The already proved
all-ring polynomial-coordinate argument makes its Jacobian a polynomial
unit there. A polynomial unit matrix lifts through a nilpotent ideal.
Finite polynomial Newton correction, starting at q, therefore yields
an exact polynomial fixed point over R: the equation and inverse-matrix
errors square at each Newton/Schulz step, so ceil(log_2 L) steps suffice
when I_B^L=0 and an inverse is supplied modulo I_B. Formal uniqueness
identifies this vector with g. Then(10) proves safety. Applying this over
R/Ann(c) proves the Gamma formula. The finite residual producer is exact;
a general inverse oracle for arbitrary abstract R is not asserted.

This uses the classical inverse-tail/nilpotent correction mechanism already
retained in [the inverse-section proof](NILPOTENT_INVERSE_SECTIONS_AND_EXACT_SCALING_IDEALS.md),
with the explicit graph-degree bound above replacing the Keller-specific
inverse-degree bound. It is not a new nilpotent Newton method.
The purity hypothesis is essential. For H=e+z over Q[e]/e^2 the residual
ideal at the field bound B=2 is (e), hence nilpotent, but g=et/(1-t) is
NOT polynomial. Its Jacobian 1-t is not a polynomial unit. Without the
pure source, residual nilpotence alone does not lift polynomiality.
Nilpotent corrections may exceed the field bound: for H=1+ez over e^4=0,
B=2 but g=t+et^2+e^2t^3+e^3t^4. Both boundaries are retained explicitly.

## 5. No zero weights: finite ordinary purity and all-ring safety WITHOUT purity

Suppose there are no Z,W pairs. Then J=1, p and psi are rational in t,
sigma in R[t], and every original period is f times a rational function
whose denominator is a power of sigma. Retain an explicit polynomial
numerator and nonnegative integer M:

    psi=nu/sigma^M, nu in tR[t].

The finite block construction supplies this data; no primitive or
antiderivative is guessed. Define

    Omega=2 sigma nu'-2M nu sigma'-sigma^M sigma',
    2 sigma^(M+1) f'=Omega f, f(0)=1.                   (13)

The FULL ORDINARY all-positive pure ideal is exactly the ideal generated
by the finitely many coefficients of Omega. If Omega has degree k,
the FIRST k+1 ORIGINAL PURE MOMENTS generate it; if Omega=0, the ideal
is0. Proof: formal inversion of f and sigma expresses every Omega
coefficient in the ideal of the first k+1 f-coefficients; modulo those
Omega coefficients, (13) gives f'=0 and hence f=1 over a Q-algebra.
Factorials are rational units. No radical replacement is made.
Changing the denominator power M merely multiplies Omega by a power
of the constant-one polynomial sigma and leaves its full coefficient
ideal unchanged. The prefix need not be optimal.

Let I be the finite ideal of all nonconstant coefficients of sigma AND
all coefficients of nu. Over EVERY R, without a purity premise,

    all-marker safety iff f polynomial iff I nilpotent,
    N_P=Gamma_I(R)=union_N Ann_R(I^N).                   (14)

For necessity reduce a polynomial f modulo each prime. Fieldwise exp(psi)
algebraic forces psi=0, and f^2sigma=1 then forces sigma=1 and f=1.
Thus every generator of I vanishes in every prime; finite generation
makes I nilpotent. Conversely I nilpotent makes sigma a polynomial unit,
all critical/Morse/covariance coefficients polynomial, psi a polynomial
with nilpotent coefficient ideal, and exp(psi) a FINITE polynomial.
Every marked image is polynomial. Applying this equivalence over
R/Ann(c) proves the normalization assertion.

The finite certificates give numerical affine bounds too. Given I^L=0,
expand every inverse using the finite nilpotent series. Let a bound the
t-degree of the coefficients of Theta, b its maximum residual ordinary
degree (use b>=1), c the t-degree of B_p^-1, and d0 the degree of f.
Then every original marker Q has the sufficient cutoff

    m>d0+ceil(a+b*c/2) deg Q.                           (15)

This follows directly from Wick contraction of Q(Theta). Under purity
d0=0 gives a linear bound. These are computed sufficient bounds, not
least-cutoff claims. The earlier homogeneous strictly positive r-moment
gate and explicit dimension/degree formula remain stronger in their
scope; neither is silently extended to lower-weight terms.

## 6. Essential retained controls and an exponential realization application

The old [delayed first-moment family](ALL_DEGREE_STRUCTURE_AND_BARRIERS.md)
X^d+Y, d>=2, with E[XY]=1 is included with weights 2/d,2-2/d. Its
seed exp(t^(d+1)) was already known; no new degree-free prefix obstruction
is claimed. Here (3)-(4) recover its exact full actual module
exp(t^(d+1))K[t], generic rank1 and infinity slope d+1. Its first d
pure moments vanish. This explicitly prevents carrying the old homogeneous
r-moment gate into the filtered class.

The ordinary nonreduced control is

    P=e*(T^2-1)/2, sigma=1-et, psi=-et/2,
    Omega=e^2*t, I=(e).                                 (16)

Over Q[e]/e^2 the seed is exactly1, though psi!=0 and sigma!=1. The
universal ordinary pure ideal is (e^2), whereas the safety/normalization
ideal topology is that of (e). This is a concrete reason to retain the
finite log identity and the distinct ordinary pure and safe ideals.

Finally EVERY algebraic psi in tK[[t]] has a seven-variable realization
of the full exponential module exp(psi)K(t,psi). Apply the previous
algebraic encoder to psi/t, obtaining a polynomial coordinate marker
C(U,Z), and add C to its pure compensated phase. All these coordinates
have weight0, so the new phase satisfies the filtered condition. Its
unit seed is exp(t C(t,g))=exp(psi), and(6) gives the stated exact
generic image. It is nonpure unless psi=0 over fields. This does NOT
realize a transcendental period inside a pure filtered phase: cancelling
the unit exponential forces the algebraic case(8).

For a ramified actual-slope example, use the five-variable compensated
lift of g=t(1+g^3) and add the original zero-position marker Z^3 to its
phase. Then psi=tg^3=g-t. On the rational curve with coordinate u=g,

    t=u/(1+u^3), psi=u^4/(1+u^3).

Above t=0 the selected point u=0 is regular, while u=infinity has
ramification2 and action pole order1. Thus the FULL rank3 actual image
has slope0 with multiplicity1 and slope1/2 with multiplicity2 at time0.
At time infinity the three roots of 1+u^3 each have ramification1 and
action pole order1, giving slope1 with multiplicity3. The corresponding
Newton irregularities are1 and3. These assertions follow from the exact
curve and full-image equality(6), not from a prefix of period jets.

The finite checker uses explicit constant-one t-adic binomial and
exponential jets. A symbolic expression such as sqrt((gamma*t-1)^2)
can select the negative constant branch; it is not a valid substitute
for the normalized formal square root in(5).

## 7. Published-source and completion boundaries

The new published primary source checked here is
[Khovanskii--Tronsgard,Izv.Math.90:1(2026),238-254](https://doi.org/10.4213/im9686e),
with its full [HTML body](https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=9686&wshow=paper).
Sections1-2 define exponentials/logarithmic derivatives and distinguish
generalized quadratures. Its Theorem2.1 on generalized special extensions
is ANNOUNCED; the paper explicitly proves a weaker result excluding
algebraic extensions. We import no unproved general theorem or Gaussian
classification from it. Our rank-one pole, Galois and slope statements
follow directly from the displayed algebraic-exponential decomposition.
Classical Euler/Newton conventions use the retained earlier source audit.

The larger filtered architecture and its full Gaussian-automorphism orbit
are now closed at the stated generic type and pure all-ring safe scopes;
the strictly positive part has finite ordinary purity,exact safe ideal,
normalization and numerical certificates. Regular transcendental Gaussian
realization must escape THIS class, not merely its old homogeneous part.
Arbitrary negative-weight/higher-filtered-degree phases,exhaustive general
safe normalforms,full actual generic connection classification and sharp
real cutoffs remain. No unrelated Mathieu work,outreach,delegation,capped
job restart,assurance promotion or goal completion. Oaku operators,torsion
witnesses and specialization thresholds stay sealed. Sabbah/Maisonobe keep
convergent rational-slope credit;nonconvergent formal/allring priority is
separately unassessed.
