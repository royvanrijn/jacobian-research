# Squarefree linear-T coefficients exclude pure Gaussian nullity in every degree

15 November 2026. Internal written all-degree theorem, using the
[logarithmic curve criterion](LOGARITHMIC_CURVE_PERIOD_CRITERION.md).
Its analytic input remains independently unreviewed. The
[checker](verify_squarefree_linear_t.py) verifies universal coordinate and
tangent identities and representative partial fractions; it does not
certify the analytic input. No external review, formal verification or
novelty claim is made.

## Theorem and structural consequence

Let E(ZW)=E(T^2)=1 and all other Wick pairings be zero. For arbitrary
complex polynomials a(Z), b(Z) and G(Z,T), set

    P=W*(a(Z)+b(Z)*T)+G(Z,T).

If b is nonconstant and squarefree, P has a nonzero positive pure Gaussian
moment. There is no degree bound on a, b or G. Thus this entire architecture
has an empty pure-null locus and satisfies GMC vacuously. The conclusion
includes every coefficient specialization preserving squarefreeness and
nonconstancy of b; no condition on a(0) or b(0) is imposed.

Equivalently, within this arbitrary-degree linear-T architecture, any
pure-null polynomial must have b constant (including zero) or b with a
repeated complex root. This is a necessary condition, not a sufficient
condition for pure nullity or for GMC failure. It leaves a sharply defined
structural boundary to classify.

Both exceptions occur for safe polynomials: b=1,a=G=0 gives P=WT,
and b=Z^2,a=G=0 gives P=WZ^2T. In the first example every monomial of
P has negative circular charge; in the second it has positive charge.
All positive pure moments vanish and mixed moments eventually vanish by
the same charge count. Hence neither constant b nor repeated-root b can
be excluded by a blanket theorem. The existing affine-W classification
retains both safe and failing examples on the b=0 boundary.

## 1. Exact formal Gaussian fixed point

Put h(z,T)=a(z)+b(z)T and solve the unique formal equation

    g=t*h(g,T),           g in t*C[T][[t]].

The retained formal Gaussian--Lagrange identity GL1 gives

    I(t)=E exp(tP)=E_T[exp(B)/D],
    B=t*G(g,T),          D=1-t*(a'(g)+b'(g)T).       (1)

This is a formal identity. It does not assume convergence of exp(tP) under
Gaussian integration. All coefficients of g,B and 1/D are polynomials in
T. We treat localization cases before constructing any quotient functional.

If b(0)!=0, b(g) is a formal unit with nonzero constant leading coefficient.
If b(0)=0 and a(0)!=0, squarefreeness gives b'(0)!=0, and

    g=a(0)*t+O(t^2),
    b(g)=b'(0)*a(0)*t+O(t^2).

Thus b(g) is still a Laurent unit with constant nonzero leading coefficient.
In both cases every q in R=C(t)[z,1/b(z)] has a well-defined substitution
q(g) in C[T]((t)) coefficientwise. Its negative t powers have a finite
lower bound. The remaining case a(0)=b(0)=0 is treated in section 4;
inverting g or b(g) there would be invalid.

## 2. The rational curve and full volume form

On the constraint z=t*h(z,T), wherever b(z)!=0, solve for T:

    T=F(t,z)=(z/t-a(z))/b(z),
    Phi(t,z)=t*G(z,F)-F^2/2.

This is a rational phase of finite degree, however large the degrees of
a,b,G are. Differentiation gives

    F_z=(1-t*(a'(z)+b'(z)*F))/(t*b(z)),
    dT/(1-t*h_z)=dz/(t*b(z)).                      (2)

The actual volume form is therefore rho*omega with

    rho=1/t,             omega=dz/b(z).

Because b is squarefree and nonconstant, omega has only simple finite
poles, at its distinct roots. At infinity its order is deg b-2 in the
coordinate y=1/z: it has a simple pole when deg b=1, is regular when
deg b=2, and vanishes when deg b>2. Its finite residues are the constants
1/b'(r), and its infinity residue is their negative sum. In particular
all residues are bounded at t=infinity.

Apply the logarithmic curve criterion on A1 minus the roots of b,
with infinity also retained in the compactification. It excludes every
nonzero constant normalized period

    (1/t)*integral exp(Phi(t,z))*dz/b(z).           (3)

No bound on phase degree, number of roots or seed cyclicity is required.
It remains to connect (1) to the actual connection of (3).

## 3. Formal horizontal functional with every tangent term retained

Use the identification of rational relative one-forms with functions q
times omega. The twisted quotient is the finite-dimensional vector space

    H=R / {b(z)*(q_z+Phi_z*q): q in R}.

The displayed set is a linear image; H is not a ring quotient. Since b is
invertible on the curve, this is exactly its twisted de Rham H1. The
normalized tangent connection on functions is

    nabla[q]=[q_t+(Phi_t-1/t)*q].                  (4)

Here omega has no t dependence. Its moving fixed-point substitution still
has full implicit derivatives. Define

    J(q)=E_T[q(g)*exp(B)/D].                       (5)

The unit checks in section 1 make this a well-defined C(t)-linear formal
functional. Implicit differentiation gives g_T=t*b(g)/D, and on the
constraint Phi(t,g)=B-T^2/2. Hence

    J(b*(q_z+Phi_z*q))
      =(1/t)*E_T[(d/dT)(q(g)*exp(B))-T*q(g)*exp(B)]
      =0,

by coefficientwise Gaussian integration by parts.

For horizontality put w=1/D and V=g_t/g_T. These expressions are valid
formal Laurent functions on the admitted branches. Since

    w=g_T/(t*b(g)),        w*V=g_t/(t*b(g)),

commutation of mixed derivatives proves the full conservation identity

    w_t+w/t=(w*V)_T.                             (6)

After differentiating (5), subtract its value on (4). The integrand of the
difference, after the common exp(B) factor is removed, is

    (w_t+w/t)*q(g)+w*V*(d/dT q(g)+(B_T-T)*q(g)).

By (6) its Gaussian integral is zero. Thus J is horizontal for the actual
connection, with no dropped h_z, implicit g_t, Jacobian or normalization
term. In particular J(1)=I(t).

If every positive pure moment vanished, I(t)=1. The cyclic-submodule
argument and local analytic extension in the curve criterion would
therefore produce the forbidden constant period (3). When Phi is constant
on the fibre, use the elementary constant-phase case of that criterion
instead of its algebraic-independence comparison. This proves the theorem
on both branches handled in section 1.

## 4. The removed formal component, a(0)=b(0)=0

Here h is divisible by z and the unique formal solution is g=0 identically.
Formula (1) becomes

    I(t)=E_T[exp(t*r(T))/(1-t*(lambda+mu*T))],
    r(T)=G(0,T),       lambda=a'(0),       mu=b'(0).

Squarefreeness at the root zero gives mu!=0. All higher coefficients of
a,b and all G terms divisible by z disappear from this exact pure series.
Do not localize this zero branch at b(g).

Instead set

    T0=(1/t-lambda)/mu,    x=T-T0,
    Phi=t*r(x+T0)-(x+T0)^2/2,
    dT/(1-t*(lambda+mu*T))=-dx/(t*mu*x).          (7)

The formal x=T+lambda/mu-1/(mu*t) is a Laurent unit with nonzero constant
leading coefficient. Thus the polynomial phase on C* in (7) has a valid
formal functional. The same integration-by-parts and conservation argument
applies with kappa=-mu and x_T=t*kappa*x/D=1, exactly the trapped-component
case already retained in the quadratic-W proof. The logarithmic form dx/x
has constant residues and the normalization -1/(mu*t) decays. Its criterion
again rules out I(t)=1, completing every case.

## Consequences and remaining boundary

This is an all-degree empty-locus theorem, not another cutoff optimization.
It extends the one- and two-puncture models to arbitrarily many roots of b.
For example b=1+Z^2 allows arbitrary-degree a and G and is already outside
the old quadratic h(Z,T) restriction. A root at zero, such as b=Z*(Z-1),
is included even when the formal branch lies entirely on that removed
component.

The residual architecture is exact: constant b or repeated-root b. The
degree-free affine-W classification handles b=0 but does not classify
every constant or repeated-root branch here. A promising structural next
question is whether the seed's higher-pole form on repeated-root branches
has a bounded logarithmic representative modulo twisted exact forms, and
which branches instead support Gaussian conjugacies or persistent mixed
failures. The examples above prohibit an indiscriminate exclusion.

The theorem depends on the internally written logarithmic-curve analytic
proof; no independent assurance is transferred from the auxiliary checks.
Least cubic prefixes, expanded radical cofactors and completed witnesses
are unchanged. [Source audit](squarefree_linear_t_source_audit.json).
