# Cable genus-three curves have a nonzero action vector

Every irreducible rational fibre with affine-line normalization in the
primary genus-three CABLE ambient class has at least one nonzero
matching-endpoint action integral

    integral X*dY.

This includes repeated-contact and triple-point fibres. An individual
action can vanish; the assertion concerns the ENTIRE vector. The proof
converts the cable action to a weighted intrinsic-quartic action and
computes that weight exactly. It does not treat the intrinsic cubic
and quartic functions as a polynomial target automorphism.

The result is a pointwise gate. It does not exclude every cable family,
all conductor degree six, higher conductors, automatic target integrality,
HC4 or JC2.

## 1. The exact cable model and its intrinsic quartic

The primary form retained from Neumann, Invent. Math.98(1989),
Section8, pages482-483, is

    f1=(X^3+Y^2)^2+a*X^4+b*X^3+(a^2/4)*X^2+c*X
                              +a*X*Y^2+b*Y^2+D*Y,
    D!=0.                                               (1)

Its exact transcription and source PDFs are pinned in the
[genus-three source audit](sources/neumann-genus-three-source-audit.json).
We assume an irreducible fibre f1=value with normalization C[t].
This is the explicit cable hypothesis, not a classification of every
genus-three embedding.

Put

    H=X^3+Y^2+(a/2)*X+b/2,
    alpha=c-a*b/2, beta=b^2/4+value.

The fibre equation becomes

    F=H^2+alpha*X+D*Y-beta=0.                            (2)

In its ACTUAL coordinate ring,

    Y=(beta-H^2-alpha*X)/D,
    C[X,Y]=C[H,X].                                      (3)

Eliminating Y gives the intrinsic equation

    Psi(H,X)=(beta-H^2-alpha*X)^2+D^2*X^3
                     +(a*D^2/2)*X+b*D^2/2-D^2*H=0.     (4)

The two quotient rings in (3) are isomorphic, so Psi is irreducible.
It has H degree four and X degree three. Since the common fraction
field is C(t), its normalization functions have EXACT degrees

    deg_t H=3, deg_t X=4.                               (5)

Indeed [C(t):C(X)]=4 and [C(t):C(H)]=3 by the irreducible equation;
for a polynomial normalization function this field degree equals
its polynomial degree.

The special intrinsic coefficients in (4) are essential. Its H*X^2,
H^3 and H*X coefficients are ZERO, while its linear H coefficient
is -D^2, which is NONZERO. The change to (H,X) is not an ambient
target automorphism and does not preserve the constant Jacobian.

## 2. The cable coefficients force a restricted normalization

Translate the affine-line parameter to remove the quadratic term
of H and divide by the nonzero leading coefficients. Write

    H=h3*P, X=x4*Q, h3*x4!=0,
    P=t^3+A*t+eta,
    Q=t^4+L*t^3+B*t^2+Cphys*t+E.

For the coefficient calculation only, subtract L*P and a constant
from Q, obtaining

    p=t^3+A*t, q=t^4+B*t^2+C*t.

The exact resultant of these two polynomials is

    R(U,V)=U^4+3C*U^3+(3B-4A)*U^2*V
       +[B*(A-B)^2+3C^2]*U^2+C*(3B-5A)*U*V
       +C*[A*(A-B)^2+C^2]*U
       -V^3+2A*(A-B)*V^2-A*[A*(A-B)^2+C^2]*V.          (6)

The actual relation is R(P-eta,Q-L*P-E0), for a constant E0.
Its U*V^2 coefficient is 3L. Comparison with (4) forces L=0.
Then C=Cphys, and the U^3 coefficient is 3C-4eta, giving
eta=3C/4. The U*V coefficient becomes

    C*(2A-3B)/2.

It too must vanish. If C=0, then eta=0 and the entire resultant
is even in U, contradicting the NONZERO linear H coefficient of
(4). Thus C!=0 and B=2A/3. We have proved

    P=t^3+A*t+3C/4,
    Q=t^4+(2A/3)*t^2+C*t+E, C!=0.                      (7)

The intrinsic subtraction of L*P was only a coefficient calculation;
comparison ultimately proves L=0 in the actual normalized pair.
A constant shift E of Q has no effect on any action below.
No classification of arbitrary abstract (3,4) orders is asserted.

## 3. The ambient action is a weighted quartic action

There is an exact polynomial one-form identity

    H^2*dX-D*X*dY=F*dX+dS,
    S=-(alpha/2)*X^2-D*X*Y+beta*X.                       (8)

On the fibre F=0, the endpoint values of S agree whenever X,Y
agree. By (3), these are exactly the pairs on which H,X agree.
Therefore every matching-endpoint cycle satisfies

    integral X*dY=(1/D)*integral H^2*dX
                 =(h3^2*x4/D)*integral P^2*dQ.           (9)

This is a WEIGHTED action, not the ordinary intrinsic action
integral P*dQ. All factors in front are nonzero constants for the
individual fibre. If an ambient polynomial automorphism changes a
generator pair, all matching-endpoint actions are multiplied by
its same nonzero constant Jacobian, as in the retained nodal-action
proofs. Thus vector nonvanishing is an ambient invariant; a uniform
automorphism over a family is not required.

## 4. Exact weighted action of every identified pair

Let t_+,t_- be distinct points with equal values of P,Q in (7), and put

    sigma=t_++t_-, pi=t_+*t_-, gap=t_+-t_-.

Their equality conditions are exactly

    pi=sigma^2+A,
    g(sigma)=sigma^3+(4A/3)*sigma-C=0,
    gap^2=-3sigma^2-4A.                                 (10)

Conversely every root of g gives such a pair. Its gap cannot vanish:
if 3sigma^2+4A=0, equation g=0 would give C=0, impossible.

Direct polynomial integration and reduction using (10) give

    W(sigma)=integral_(t_-)^(t_+) P^2*dQ
             =gap*C^2*(8A*sigma+9C)/40.                 (11)

For an exact verification with no numerical roots, a primitive is

    L(t)=2t^10/5+7A*t^8/6+C*t^7+10A^2*t^6/9+2AC*t^5
       +(A^3/3+15C^2/16)*t^4+A^2C*t^3
       +9AC^2*t^2/8+9C^3*t/16.

Compute [L(t_+)-L(t_-)]/gap by the recurrence

    h0=0,h1=1, hn=sigma*h(n-1)-pi*h(n-2).

Here hn=(t_+^n-t_-^n)/gap. Substitute pi=sigma^2+A and
reduce modulo g; the result is the coefficient in (11).
This proves (11) as a universal polynomial identity, not by sampling.

If A=0, that coefficient is 9C^3/40, which is nonzero.
If A!=0, it is a nonzero linear function of sigma and can vanish
at at most ONE value. But the cubic g has at least TWO distinct
roots: a triple root would force A=C=0. Every corresponding gap
is nonzero. Hence the matching-endpoint action vector is NONZERO.

The argument includes a double root of g, which represents contact,
and A=0, where the three matching pairs share three source points
at one ordinary triple point. It does not identify the number of
pairs with the number of independent cycles. One nonzero matching
pair already gives a nonzero action on the corresponding cycle.
Equations (9),(11) prove the asserted cable gate.

## 5. Exact rational cable controls and a vanishing individual action

The nonzero-vector conclusion does not imply nonzero individual
periods. There are actual polynomial normalization controls for this
boundary. Put

    Z=t^2+A/3,
    H=t^3+A*t-3/2,
    X=Z^2-2t,
    Y=i*(Z^3-3t*Z+3/2).

They satisfy identically

    H=X^3+Y^2+A*X+3/4,
    H^2-A*X+i*Y=3/4-4A^3/27.                            (12)

Thus they are fibres of (1), with a=2A,b=3/2,c=A/2,D=i.
Their normalization is C[t]: H and X have coprime degrees three
and four, generate C[X,Y] by the first identity, and have common
fraction field C(t). The map is immersed: if H_t=0, then X_t=-2.
The monogenic trace-dual formula from
[the coordinate-order proof](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives the exact conductor ratio

    kappa=t^6+(5A/3)*t^4-4t^3+(7A^2/9)*t^2
                                  -(10A/3)*t+A^3/9+4.

It has degree six and normalization defect three. For A=0,
kappa=(t^3-2)^2, giving the three-branch triple-point control.
These are rational plane curves, not two-variable Keller maps.

In this C=-2 normalization, equation (11) is

    W/gap=(4A*sigma-9)/5.

Choose sigma^3=-5 and A*sigma=9/4. Then A^3=-729/320,
g(sigma)=0 and W=0, while gap^2=6/sigma!=0.
The cubic discriminant is -432/5, so the other pair sums are distinct.
At least one of their actions is nonzero by Section4. This is a
proper guard against an incorrect individual-action claim.

## 6. Scope and exact evidence

The checker verifies the primary cable square completion, the actual
quotient-ring elimination, all coefficient restrictions in (6),(7),
the one-form identity (8), the full primitive and symmetric reduction
(11), the nonzero-vector argument, and the rational curve controls
(12) including the zero individual action. The primary classification
remains an attributed written input, not a theorem proved by these
symbolic calculations.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_cable_genus_three_weighted_periods.py

The result excludes an ALL-ZERO action vector at a cable fibre.
It neither turns (H,X) into a target automorphism nor proves a
Jacobian-preserving reduction to the quartic branch. The later
[global sextic exclusion](GLOBAL_SEXTIC_CONDUCTOR_EXCLUSION_AND_PERIOD_POLE_REPAIR.md)
closed global conductor degree six, including its cable branch. The
[current global bound](GLOBAL_CONDUCTOR_TWENTY_EXCLUSION.md) is twenty-two.
Individual special fibres may still have degree two, four or six.
Higher conductors and automatic target/trace integrality remain open.
Independent replay, formal verification and external review are
not claimed.
