# Seven Gaussian-phase critical points and the nodal escape profiles

14 November 2026. Internal written algebraic proof and exact symbolic
checks. This algebraic note does **not** itself exclude a new nonhomogeneous GMC class.
The subsequent [center-line cut-decay theorem](NODAL_CENTER_LINE_DECAY.md)
excludes its zero-critical-line and closed-half-plane critical-value classes.
No external review, independent specialist replay or formal verification
is claimed.

The subsequent [full formal-frame theorem](NODAL_FORMAL_FRAME.md)
closes every irreducible nodal-leading cubic, including all finite
degeneracies, coincident critical values and quartic resonances. Earlier
method boundaries below describe the scope of this retained argument.

## Scope and normal form

Consider a ternary cubic whose homogeneous leading curve is irreducible
nodal, with nonisotropic node v and F2(v)=0. The retained
[height-one reduction](NODAL_CUBIC_REDUCTION.md) and binary centering give
P=T*(uv+h)+G(u,v). Absorb the mixed cubic monomials into T, then the
quadratic uv monomial into a constant translation of T. The result is

    P=T*(uv+h)+a*u^3+d*v^3+b*u^2+c*v^2+e*u+f*v+k,
    a*d!=0.                                             (1)

The T-shears change the linear and constant terms of G when h!=0; those
terms must be retained. All changes are affine and invertible. Keep the
full transformed Gaussian exponent

    Q=(A*u^2+B*v^2+m*T^2+2*C*uv+2*D*uT+2*E*vT)/2
        +l_u*u+l_v*v+l_T*T+q0,
    m!=0, det Hess Q!=0.                                 (2)

The coefficient m is nonzero because the node is nonisotropic; centering
and shearing do not change the T-axis direction. Means, cross covariances
and the constant coordinate Jacobian have not been discarded. Write

    S_t=t*P-Q.

For every t!=0, the critical scheme of S_t has length **seven**.
As t tends to infinity, its escaping branches are exactly those in the
table below, including their multiplicities. All escaping branches are
eventually Morse and have a nonzero positive-power action. If h=e=f=0,
there are no escaping branches. This completely locates the algebraic
escape possibilities; it is not a uniform period-decomposition theorem.

## The seven-point count

Homogenize the three quadratic equations grad S_t=0 in projective
three-space. At infinity their only common point is the leading node
[u:v:T:w]=[0:0:1:0]. In its chart T=1, the linear parts in (u,v,w) are

    t*v-D*w,  t*u-E*w,  -m*w.

Their Jacobian determinant is m*t^2!=0. Hence this projective point is
isolated and has intersection length one. No positive-dimensional
component can lie in the affine part: its projective closure would meet
infinity and contradict this isolated point. The three quadrics therefore
form a zero-dimensional complete intersection of degree8. Removing the
one point at infinity leaves affine length7. This is a projective Bezout
proof, not an inference from a sampled gradient basis. The checker also
retains one explicit seven-dimensional quotient as a regression.

## Finite critical points of P

For h!=0 the critical conic uv=-h is smooth. Its restricted phase is

    R(u)=a*u^3+b*u^2+e*u+k-f*h/u+c*h^2/u^2-d*h^3/u^3.

Multiplying u*R'(u) by u^3 gives

    3*a*u^6+2*b*u^5+e*u^4+f*h*u^2-2*c*h^2*u+3*d*h^3. (3)

This has degree6 and nonzero constant term. Thus P has exactly six
finite critical points with multiplicity. Locally, the regular coordinate
uv+h and a shift of T split the phase as T'*(uv+h)+R(u), so root
multiplicity in (3) is exactly the local Milnor number. Multiple roots
and coincident critical values are permitted.

For h=0, a nonzero point on v=0 is critical precisely when

    3*a*u^2+2*b*u+e=0, T=-f/u,

and the other axis has the symmetric description. The center is critical
if and only if e=f=0; then it is the entire T-axis. These facts give the
finite Milnor sums in the table. Conservation of local intersection
multiplicity under the perturbation P-Q/t retains each finite isolated
critical cluster with its original Milnor number.

## Exhaustive escaping branches

In the first row put K=h/m. In the remaining rows assume e!=0 and let
g_j*v^j be the lowest nonconstant term of G(0,v):

    j=1: g_1=f!=0;
    j=2: f=0, g_2=c!=0;
    j=3: f=c=0, g_3=d!=0.

If e=0,f!=0, exchange u,v and the corresponding coefficients. Thus these
rows exhaust h=0 with at least one nonzero linear coefficient.

| Conditions | Finite Milnor sum | Escapes n | T scale | Equation for K | Action of an escape |
|---|---:|---:|---|---|---|
| h!=0 | 6 | 1 | T~K*t | K=h/m | (h^2/(2m))*t^2+O(t) |
| h=0, e*f!=0 | 4 | 3 | T~K*t^(1/3) | K^3=e*f/m | k*t-(3m*K^2/2)*t^(2/3)+lower powers |
| h=0, e!=0, f=0, c!=0 | 3 | 4 | T~K*t^(1/4) | K^4=-2*c*e^2/m | k*t-m*K^2*t^(1/2)+lower powers |
| h=0, e!=0, f=c=0 | 2 | 5 | T~K*t^(1/5) | K^5=3*d*e^3/m | k*t-(5m*K^2/6)*t^(2/5)+lower powers |

Every K is nonzero. The common term k*t must be removed before extracting
the fractional action. It does not cancel a nonzero fractional term. For
n=4, opposite roots K and -K have the same leading action; distinct Morse
points must not be assumed to have distinct formal critical values.

For h!=0, use r=1/t and the substitutions

    u=r*U, v=r*V, T=K/r.

Multiplying each gradient equation by r and setting r=0 gives

    K*V+e-D*K=0, K*U+f-E*K=0, h-m*K=0.

Its Jacobian determinant in (U,V,K) is m*K^2!=0. The analytic implicit
function theorem gives one Laurent branch. Its action and original
Hessian determinant have leading terms

    S_t=(h^2/(2m))*t^2+O(t),
    det Hess S_t=(h^2/m)*t^4*(1+O(1/t)).              (4)

This permits cancellations in the leading u or v coordinate; the branch
and Hessian remain nondegenerate even then.

For the other rows put n=j+2, r=t^(-1/n), and

    u=r^j*U, v=r*V, T=K/r.

After subtracting k*t, the leading scaled phase is

    r^2*(S_t-k*t)=K*U*V+e*U+g_j*V^j-m*K^2/2+O(r). (5)

All terms in the remainder have positive integer r-powers near a fixed
finite (U,V,K). The critical equations of the displayed polynomial are

    K*V+e=0, K*U+j*g_j*V^(j-1)=0, U*V-m*K=0.

They imply

    K^(j+2)=(-1)^(j+1)*j*g_j*e^j/m.

There are exactly n simple nonzero roots. The critical Jacobian has
determinant n*m*K^2; the implicit function theorem gives n distinct
Puiseux branches, with

    S_t=k*t-((j+2)/(2j))*m*K^2*t^(2/(j+2))
             +lower fractional powers,
    det Hess S_t=n*m*K^2*t^(2+2/n)*(1+O(t^(-1/n))).  (6)

The finite critical multiplicities plus these constructed branches already
sum to7 in every row. The seven-point count proves that there are no
additional escaping branches, hidden axis branches or unaccounted critical
multiplicities on this finite-critical locus.

Equations(4)--(6) also identify the local Gaussian scale of a Morse
thimble: after factoring its own exponential action, its leading scalar
amplitude is a nonzero constant times t^(-1-1/n). For (5), the volume
factor r^(n-2) and three-dimensional Morse scale r^3 give r^(n+1).
For (4), the transverse variables have scale r and the longitudinal
Gaussian width stays constant, giving r^2. These are local statements.

## The center critical line: no escaping branches

Suppose h=e=f=0. An unbounded S_t-critical branch as t->infinity would
have projective direction the node. Thus T->infinity and (u,v)=o(T).
The first two gradient equations, divided by t*T, give

    |u|+|v| <= C/|t|+o(|u|+|v|).

The quadratic errors divided by T are absorbable because (u,v)=o(T).
Hence u,v=O(1/t). The last equation is then

    t*u*v=m*T+l_T+D*u+E*v,

whose left side tends to zero while the right side is unbounded, a
contradiction. Every one of the seven phase-critical points therefore
stays bounded as t->infinity, including exceptional coefficient values.

Apart from the center line, P has at most two critical points:

    (-2b/(3a),0,0),  value k+4*b^3/(27*a^2), if b!=0;
    (0,-2c/(3d),0),  value k+4*c^3/(27*d^2), if c!=0.

Each is Morse. On the center line the transverse Hessian determinant is
4bc-T^2. If this Hessian is invertible at a limiting point, u,v=O(1/t)
forces m*T+l_T=0. Thus all center-line limits lie among

    (m*T+l_T)*(T^2-4bc)=0.                           (7)

This statement does not give their separate local multiplicities or
identify a complete period basis at the colliding exceptional points.

Along a bounded critical branch put epsilon=1/t and
A(epsilon)=P(z(epsilon))-epsilon*Q(z(epsilon)). The critical equations
give A'=-Q(z(epsilon)), which is bounded. A branch tending to the center
line therefore has A=k+O(epsilon), or S_t=k*t+O(1). This controls its
critical action, not its period amplitude. In particular k=0 is the
remaining zero-action center-line case. When b=c=k=0 and Q is centered, the homogeneous case is covered by the
[Newton theorem](NODAL_NEWTON_EXCLUSION.md). A shifted Gaussian mean was
not covered by that argument. The subsequent
[center-line theorem](NODAL_CENTER_LINE_DECAY.md) supplies the stronger
result, including every k=0 member with arbitrary b,c and Gaussian means.

## Actual Gaussian width in the compact conic chart

Write Q=m*T^2/2+L(u,v)*T+Q_bin(u,v). Completing the Gaussian square is exact:

    S_t=-m/2*(T-(t*H-L)/m)^2
          +t*G-Q_bin+(t*H-L)^2/(2m), H=uv+h.          (8)

After eliminating T, put s=t*H and v=(s/t-h)/u on u!=0. The remaining
phase is

    Phi=(s-L(u,(s/t-h)/u))^2/(2m)
          +t*G(u,(s/t-h)/u)-Q_bin(u,(s/t-h)/u),
    du*dv=(1/t)*(du/u)*ds.

For u in a compact subset of C* and bounded s, Phi_ss ->1/m. Thus the
specific width-t cancellation in the earlier negative control cannot
occur on this compact conic chart. This is not uniform as u->0 or infinity,
or near s=t*h. Those are precisely the bounds an analytic proof must supply.

## What this closes, and what remains

The [checker](verify_nodal_gaussian_critical.py) checks the algebraic formulas in (3)--(8), the
projective Jacobian, each rescaled critical system, action and Hessian,
all three j-profiles, the exact Gaussian completion, the compact-chart
width and one seven-dimensional gradient quotient. The
[receipt](nodal_gaussian_critical_verification.json) retains these data.
The universal seven-point count, conservation of Milnor multiplicity,
exhaustiveness and center-line boundedness use the adjacent written proofs;
the sampled quotient is not substituted for those arguments.
The first determinant-sign assertion and the second action extractor
(which initially omitted subtraction of k*t) failed; their inputs and
terminal failure receipts are retained. Only the corrected run passes.

[Pham's 1985 paper](https://numdam.org/item/AST_1985__130__11_0.pdf),
PartII Section1.1, supplies the local isolated-critical-point power-log
framework; its Appendix constructs semialgebraic rapid-decay chains and
proves fixed-phase convergence. The [source audit](nodal_gaussian_critical_source_audit.json)
records the exact reading coverage. Neither statement automatically gives
the required uniform parameter estimate at this nodal infinity limit.

The next analytic obligation is now precise: prove that all horizontal
solutions split sectorially into the finite local Milnor blocks and the
listed escaping Morse blocks, with tails controlled uniformly in t.
The finite critical blocks have the retained local decay theorem, and
every escape has a nonzero action, but their algebraic count alone is
not a proof of that uniform splitting. Equal actions must be grouped.
A proposed semialgebraic-chain proof must establish polynomial bounds in
both t and the cut length; fixed-t convergence does not suffice.

The subsequent [center-line theorem](NODAL_CENTER_LINE_DECAY.md) supplies
Gaussian-weighted normalized cut decay throughout h=e=f=0, including the
exceptional limits in (7), and excludes every zero-critical-line member
and the larger closed-half-plane critical-value class. Bounded critical
points alone did not prove this: the new proof uses a uniform gradient
bound, constrained semialgebraic chains and integral-current compactness.
The subsequent [dominance theorem](NODAL_MORSE_DOMINANCE.md) closes
every center-line case and the finite-critical Morse/distinct-value locus
with nonzero quartic gamma. The remaining finite degeneracy, value-collision
and quartic-resonance strata remain open. This algebraic note itself supplies no new
nonhomogeneous exclusion or numerical cubic upper prefix.
