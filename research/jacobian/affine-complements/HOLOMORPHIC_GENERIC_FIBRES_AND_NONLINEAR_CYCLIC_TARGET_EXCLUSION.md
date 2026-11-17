# Holomorphic generic fibres exclude nonlinear quadratic cyclic targets

PF2RULED146. 16 November 2026. ALL degrees over C.
This excludes the entire separated coordinate family P=a*w+f(z),
a nonzero constant and f ANY nonconstant polynomial, on quadratic MX0
and classical S2. There is NO rational Jacobian mate, even in a finite
extension of the function field. No target symmetry or companion degree
bound is assumed. The new moving-label case with deg f>=2 passes the
direct degree-at-most-one test on the reached multiple-fibre line.
Arbitrary nonlinear target coordinates and full HC4/JC2 remain OPEN.

## 1. Scope audit before using the previous exclusions

PF145's complete critical-scheme formulas, mate rigidity and generic-conic
residue calculation are retained. Its regular-target exclusion for
c*z+A(u)*w+B(u) ALREADY follows from PF52: on C0={u=w=0}, the restriction
is c*z+B(0), of degree at most ONE on a reached whole source line.
If nonconstant it is injective there; if constant its actual companion
would be affine-linear there. Either forces the Keller inverse.

The same observation excludes A(u)*z+B(u)*w+C(u) for arbitrary polynomial
coefficients. A nonconstant coefficient of z by itself is not a new
surviving case. Any actual polynomial TARGET coordinate must restrict
to degree at least TWO on every reached source-coordinate line.

On classical S2, P=a*w+f(z) has an affine restriction to the OTHER source
line h=0 in PF144's chart (u,w,z)=(t^2,t(1+t^2*h),2*h+t^2*h^2).
Thus its regular-target exclusion there also has an earlier one-line
proof. The stronger NO RATIONAL MATE conclusion below is not a
consequence of that line test.

For quadratic MX0 with b of positive degree, the source x=0 instead
has w=t*sigma(t^2), z=0, so its restriction has degree 2*deg b+1.
The C0 restriction is f(linear x), of degree deg f. For deg f>=2,
neither of these old direct line arguments disposes of the candidate.
The new obstruction computes its ACTUAL required relative differential.

## 2. A general smooth hypersurface lemma

Let

    V={w^2=N(u)+M(u)*z}, M,N in C[u],
    m=deg M > ell=deg N >=1.

Assume V is smooth. Its residue volume

    Omega=du wedge dw/M

is globally regular and nowhere zero: use the equivalent residue
expressions on the other nonvanishing hypersurface-Jacobian patches.
Let

    P=a*w+f(z), a in C*, n=deg f>=1.                 (1)

Then NO rational Q in C(V) satisfies

    dP wedge dQ=Omega.                               (2)

The conclusion also holds with any nonzero constant multiple of Omega.
The proof includes EVERY infinity valuation of the geometric generic
curve; no Newton nondegeneracy or omitted-end assumption is used.

Write p for the generic value of P and work over an algebraic closure of
C(p). Eliminating w gives the entire affine generic curve

    G(u,z)=N(u)+M(u)*z-(p-f(z))^2/a^2=0.             (3)

It is smooth. For completeness, on a smooth surface in characteristic
zero the zero locus of dP has dimension at most one, and P is constant
on the normalization of each of its curve components. There are finitely
many such components and isolated critical points, hence finitely many
critical values. The geometric generic fibre is therefore smooth.
It may be treated component by component; connectedness is not needed.

Put

    G_u=N'(u)+M'(u)*z.

Differentiating the ORIGINAL surface equation gives

    2*w*dw=(N'+M'*z)*du+M*dz,
    Omega=dP wedge [dz/(a*G_u)].                    (4)

Thus (2) would force on each generic component

    dQ=omega, omega=dz/(a*G_u).                     (5)

At G_u=0 use omega=-du/(a*G_z). Smoothness of (3) makes this a regular,
nonzero differential on the ENTIRE affine generic curve, not merely
on the open where M or G_u is invertible.

## 3. Every infinity valuation makes omega holomorphic

Take a point at infinity on the smooth projective normalization of any
generic component. Its normalized discrete valuation nu is trivial on
the geometric constant field. At least one of u,z has a pole.

If z had pole order beta>0 while u were finite, the right side
(p-f(z))^2/a^2 would have pole order 2*n*beta. The left side N+M*z has
pole order at most beta, impossible. Thus u ALWAYS has a pole. Write

    alpha=-nu(u)>=1.

If z were a nonzero finite unit, M*z would have the unique dominant
pole m*alpha, since m>ell and f(z) is finite. This is impossible too.
Exactly TWO cases remain.

Case A: z vanishes. Put beta=nu(z)>0. Since p-f(0) is a nonzero generic
constant, cancellation of the N and M*z poles in (3) forces

    beta=(m-ell)*alpha, z ~ -N/M.

The leading coefficient of G_u is
(ell-m)*N_lead*u^(ell-1), which is NONZERO. Therefore

    nu(G_u)=-(ell-1)*alpha,
    ord(omega)=beta-1+(ell-1)*alpha
              =(m-1)*alpha-1 >=0.                  (6)

Case B: z has a pole of order beta>0. The M*z pole dominates N, and
balance against f(z)^2 gives m*alpha=(2*n-1)*beta.
The M'*z term strictly dominates N' in G_u, so

    nu(G_u)=-((m-1)*alpha+beta),
    ord(omega)=-beta-1+(m-1)*alpha+beta
              =(m-1)*alpha-1 >=0.                  (7)

These exact orders use characteristic zero: a function of nonzero order
r has differential of order r-1. Equations (6)--(7) exhaust ALL boundary
points, including ramified ones, with arbitrary positive alpha.
The hypothesis m>ell>=1 ensures m>=2 and the displayed nonnegativity.

Hence omega extends to a NONZERO holomorphic differential on every
smooth projective generic component. It cannot be dQ for a rational Q:
a pole of Q of order j>=1 would give a pole of dQ of order j+1.
With no poles Q is constant on the projective curve, whereas omega!=0.
This proves the lemma without a genus classification.

## 4. Finite extensions and the actual whole-plane source

The obstruction survives every finite characteristic-zero extension of
C(V). If a rational mate existed upstairs, normalized field trace would
produce a RATIONAL mate downstairs, because trace commutes with
differentiation and P belongs to the lower field. No polynomial trace
landing is used. Equivalently, pullback at any ramified curve point has

    ord(pullback omega)=e*ord(omega)+(e-1)>=0.

The nonzero holomorphic differential stays nonzero and cannot become
rationally exact. Every ramification index is retained in this formula.

For the actual cyclic plane chart q:X=A2->V, its field extension L/E is
finite and q^*Omega is a nonzero CONSTANT times the source volume.
Thus a polynomial source companion Q, even one outside E=C(V), would
give a rational mate upstairs and then downstairs. It is impossible.
This is stronger than excluding companions restricted to O(V).

Only rational exactness is traced here. We do NOT infer Tr(R) subset B,
Tr(Lambda*Lambda) subset B, or preservation of the full finite lattice.

## 5. Apply the lemma to the actual quadratic cyclic families

For MX0's g=2 family retain ANY k>=2 and positive-degree square-free b,
b(0)!=0, r=deg b. PF145's COMPLETE invariant algebra is

    V={w^2=u*sigma(u)^2+u^2*b(u)^k*z},
    sigma=-b+2*(1-k)*u*b'.

It has

    m=k*r+2, ell=2*r+1, m-ell=(k-2)*r+1>=1.         (8)

It is the full smooth affine source image, with the whole etale plane
chart retained. In that chart u=t^2, w=t*v, and

    q^*Omega=-2*dx wedge dt.

All conditions of the lemma hold. Therefore for EVERY a!=0 and
EVERY f of degree at least ONE, P=a*w+f(z) has no rational Jacobian
mate even in the original plane function field L. In particular the
entire previously surviving nonlinear class deg f>=2 is excluded,
without assuming a target-field symmetry or restricting its companion.

For example k=2,b=u-1 and P=w+z^2 pull back to

    P=t*(1-3*t^2+t^2*(t^2-1)^2*x)
        +(2*(1-3*t^2)*x+t^2*(t^2-1)^2*x^2)^2.

Its restriction at t=0 is 4*x^2, at x=0 is t-3*t^3, and on the ordinary
line t=2 has degree four in x. Those named line tests do not supply this
certificate; no claim tests every polynomial source line. More generally,
for deg f=n>=2 the restriction to EVERY line t=c has degree 2*n in x if
c^2*b(c^2)^k!=0, and degree n otherwise, since sigma never vanishes at a
bad base root. Thus the entire intrinsic vertical-line pencil passes the
old degree test. The holomorphic relative differential still excludes P.
No actual Keller mate for this P is asserted.

For classical S2, M=u^2,N=u, so m=2,ell=1 and the same result applies.
Here the geometric generic curve has the explicit hyperelliptic model

    Y=1+2*z*u,
    Y^2=1+(4/a^2)*z*(p-f(z))^2,
    omega=dz/(a*Y).                                  (9)

The polynomial on the right is square-free over C(p). A repeated root
would give z!=0, p-f=2*z*f', and
1+(16/a^2)*z^3*f'^2=0. The last nonzero polynomial forces z algebraic
over C, and then p algebraic, a contradiction. Its degree is 2*n+1,
so the complete generic curve has genus n. Formula (9) gives an
independent view of the same nonzero holomorphic differential.
The regular exclusion on S2 was already the h=0 source-line case in
Section 1; rational-mate impossibility is the stronger conclusion.

## 6. Boundaries and the full normalization

The assumption deg f>=1 is ESSENTIAL. On S2, P=w has the rational mate
Q=1/u, since dw wedge d(1/u)=Omega. That mate is polar on the reached
multiple-fibre line and is not a regular plane target companion.

The degree hypothesis is also essential: for M=1,N=u, V is A2 with
z=w^2-u and Omega=dw wedge dz. Then w+f(z),z is a regular Darboux pair.
These controls prevent extension of the lemma to arbitrary smooth
hypersurfaces or arbitrary affine-plane functions.

The exclusion survives all automorphisms of either cyclic surface and
all polynomial TARGET rechartings: unit volumes differ by scalar units,
and a transported actual target would give a mate for the original
tested coordinate. But no reduction of a GENERAL target coordinate to
a*w+f(z) is proved. Nonconstant w coefficient, simultaneous nonlinear
u,z dependence, and arbitrary higher mixed terms remain separate.

For any hypothetical original Keller pair retain its full finite
normalization S over B, full trace-dual Lambda, and EVERY ramified and
unramified boundary prime. V and E are auxiliary source-image geometry;
they do not replace S or the target field. The proof rejects the specified
first-coordinate family before needing a finite normalization of a target.
All curve infinity valuations and all ramified finite-extension orders
have been retained, without claiming a classification of the original
surface boundary.

PF142's mixed 2+3 minimal alternative remains retired. PF143's actual
generic first-velocity marked-quadratic integrality is still UNSUPPLIED.
The full HC4/JC2 objective and construction of an actual counterexample
remain OPEN.

## 7. Sources and evidence

The nonexactness of a nonzero holomorphic differential uses the elementary
pole-order argument above. It is an application of the retained canonical
boundary-differential method in HC4MYLQ1, now with the ACTUAL required
Gelfand--Leray differential computed explicitly on the cyclic family.
Exactness of an unrelated source Liouville mark or of the surface volume
is not substituted for (5).

PF52 supplies the line-criterion audit, MX0 and PF145 supply the complete
quadratic source algebra and volume, and PF142--144 remain retained.
The new checker supports the coframe identity, both complete infinity
order calculations, the nonlinear whole-plane control, the hyperelliptic
model and the boundary countercontrols. These controls support the
written universal proof; they are not a target census or an independent
proof of HC4/JC2. No predecessor checker or large search is replayed,
and no new external theorem is imported. Formal verification,
independent whole-proof replay, external review and novelty priority
are not claimed.
