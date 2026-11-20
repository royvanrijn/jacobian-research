# Compact-complement consecutive exclusion at each integer cusp

Fix an integer s>=1 and put a=s^2, d=2s+1, D=a+d,
S_c=s(s+1), w_c=S_c^d*(s+1)/s. Define

    psi_w(S)=S^a*(w+S^d)*exp(-S),
    I_m(w)=integral_0^infinity psi_w(S)^m dS,
    A_m(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k).

The order m is always a positive integer in the period identities. Thus
psi_w(S)^m is entire in S. Direct expansion and the Gamma integral give

    I_m(w)=(a*m)!*A_m(m^d*w)/m^((a+d)*m+1).       (G1)

All constants below may depend on s and on the specified compact set.
No effective cutoff, uniform-in-s estimate, or all-order coprimality is
claimed.

**Theorem CompactComplementCuspConsecutiveExclusion.**
For every compact K subset of C minus {w_c}, there is an unevaluated
N_0(s,K) such that, for all integers N>=N_0, the actual polynomials
A_N and A_(N-1) have no common root x with x/N^d in K.

The essential new step is a complete one-Gaussian contour near a
potential Stokes wall where the critical magnitudes are different. It
is built by a finite low-level detour. Neither a Picard--Lefschetz jump
formula nor an assumption that both saddles are active throughout the
parameter plane is needed.

## 1. Explicit imported statements and their scope

The proof imports the following already proved statements from
CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md, equations(T1)--(T15), and the
distinguished sector theorem it cites.

For Im(w)>0, the physical wedge

    C={S: -pi/d<arg(S)<0}

has exactly two simple critical points S_L,S_R of a local phase
Phi=log(psi_w), and exactly one interior zero
zeta=(-w)^(1/d). The real-part separation is

    Re(S_R-S_L)>0.                                  (G2)

The closed wedge, compactified at0 and infinity, is a disk X with
three marked decaying ends E={0,zeta,infinity}. All finite closed
f-bands for f=log|psi_w| are compact away from the ends. The normalized
downward flow is

    V_w=-1/Phi_S=-psi_w/(partial_S psi_w),
    d f/dt=-1.                                      (G3)

It points strictly inward on both angular sides. Negative end caps
are contractible. Passing each of the two simple harmonic critical
points attaches one index-one handle; the two attachments form a
tree connecting the three end caps. This finite-level statement
holds also when the two critical factors have the same argument.

At a parameter where lambda_R/lambda_L is not positive real, where
lambda_i=psi_w(S_i), the actual descending thimbles are complete tree
edges. They are an integral relative basis and the original0-to-
infinity path is their unique tree path, with coefficients0 or +/-1,
not both zero. On compact subsets of one such chamber their complete
periods have uniform Gaussian expansions

    P_i(m,w)=lambda_i(w)^m*sqrt(2*pi/m)*c_i(w)
                 *(1+O(1/m)),
    c_i(w)=1/sqrt(-Phi_SS(S_i(w),w)) !=0.            (G4)

The square-root orientations agree with the corresponding paths.

Every upper-half-plane point with equal critical magnitudes lies on
the distinguished co-dominant arc. On that arc an explicit global
logarithm D_crit=Phi_R-Phi_L satisfies

    -pi<Im(D_crit)<0.

Consequently no possible Stokes wall
lambda_R/lambda_L in R_(>0) meets the equal-magnitude locus.
The envelope identity is

    d*w*partial_w D_crit=-(S_R-S_L).                (G5)

Finally, those imported statements already prove compact regular
non-Stokes consecutive exclusion using(G2), (G4), and(G5). We repeat
the comparison in Section5 below to make the real-boundary extension
explicit. The proof of(G4) does not give a uniform estimate on a
Stokes wall; Section2 supplies the missing complete contours.

All exact relative-path integral identities used here follow by
Cauchy's theorem after truncation. Small caps contribute
O(r^(am+1)) at0 and O(r^(m+1)) atzeta; large closed-sector caps
contribute O(R^(Dm+1)*exp(-mR/2)). Thus these caps vanish. A finite
sublevel cap is not silently declared to have zero integral.

## 2. A Stokes-wall neighborhood with a complete dominant packet

**Lemma StokesFiniteSublevelDetour.**
Let Im(w_0)>0 and lambda_R(w_0)/lambda_L(w_0)>0. There is a
neighborhood U of w_0 and m_0 such that I_m(w) is nonzero for
w in U and every integer m>=m_0.

Proof. By Section1 the two critical real parts are unequal at w_0.
Label them p_l,p_h so that

    h_l=log|psi_w(p_l)| < h_h=log|psi_w(p_h)|.

Shrink a relatively compact simply connected U so that the critical
points, critical factors, and local Morse coordinates are analytic,
the labels stay ordered, and their real-height gap has a fixed
positive lower bound. Choose a fixed level c and a positive eta
so that throughout U,

    h_l+4*eta < c < h_h-4*eta.                    (G6)

Every finite band between the two critical neighborhoods is compact,
uniformly after shrinking U. All endpoint caps can likewise be
chosen uniformly; around the moving zeta the cap moves analytically.

The lower critical point has a complete genuine downward thimble
throughout U, even at a possible Stokes parameter. After leaving its
local Morse chart, its f-value is strictly below h_l, so it cannot
approach the higher critical point by(G6). It cannot return to
itself after a strict phase loss. Compact finite-band continuation
and the three stable end caps give two complete arms, with fixed
end labels. Its local Gaussian arc and complete tails therefore
give

    P_l(m,w)=lambda_l^m*sqrt(2*pi/m)*c_l(w)
                        *(1+O_U(1/m)).           (G7)

This is the same elementary local-Morse and complete-tail proof as
(G4), with height separation replacing argument separation.

For the higher point, use its analytic Morse coordinate u with

    Phi(S,w)=Phi(p_h,w)-u^2/2.

Take a symmetric real local arc -r<=u<=r, oriented once and for all.
Continue its two descending arms only until the level c. No critical
point lies in this upper band. The stopping times are the exact
differences of f-values by(G3); ODE continuation is uniform on the
compact band. Both cutpoints therefore depend continuously, indeed
real-analytically locally, on w.

Below the lower handle but above all negative caps there are three
components, and just above that handle there are exactly two.
Thus X_c={f<=c} has two components. One contains two ends and the
other contains the remaining end. The two cutpoints of the upper
handle belong to different components: otherwise the upper handle
would not connect the final space X, whereas there is no remaining
handle and X is connected. This is an incidence assertion from the
two-handle tree, not an inference from rank alone.

At w_0, continue each cutpoint a little farther to a regular level
below c but still above h_l. In each of the two components, join it
to a chosen end by a piecewise smooth path in that component,
truncated at a stable end cap. Such paths exist because the
component is path-connected, and can be taken to satisfy f<c on
their finite portions. Choose one end of the two-ended component
and the sole end of the other component. Shrink U once more.
Keep the finite connector paths fixed, attach short paths from
the moving cutpoints, and move the terminal cap at zeta when it
occurs. Uniform continuity of |psi_w| on their compact union gives

    f <= c+eta < h_h-3*eta                       (G8)

on these connector pieces. One may first leave a strict margin
below c at the original cutpoints, so that the short attachments
obey the same estimate. Near finite ends use short radial pieces
inside the stable caps; at infinity use a fixed decaying sector ray.
Their weighted integrals are uniformly convergent.

This produces an actual admissible contour Q_h through p_h whose
two endpoints are the chosen ends. Its only part without a strict
loss from h_h is the symmetric local Morse arc. The finite flow
pieces outside that arc have the fixed loss r^2/2, and the
detoured pieces have(G8). Their finite lengths are uniformly
bounded. At a finite end, |psi_w| is bounded by C r^a or C r;
at infinity it is at most C |S|^D exp(-Re S). For any contour
piece with sup |psi_w|<=M<|lambda_h|, its full integral is bounded by

    M^(m-1)*integral |psi_w(S)| |dS|,             (G9)

and the last integral is uniformly bounded. Therefore the complete
detoured period has the uniform expansion

    Q_h(m,w)=lambda_h^m*sqrt(2*pi/m)*c_h(w)
                        *(1+O_U(1/m)).          (G10)

No claim is made that Q_h is the complete higher downward thimble
at the wall; it generally is a finite low-level detour.

The lower thimble joins the two ends of the two-ended component.
The detoured higher contour joins one of those ends to the third.
These two endpoint edges form a tree. Hence their end-relative
classes are an integral basis. The original0-to-infinity class is
the unique tree path and, exactly as in the imported cap argument,

    I_m=n_l*P_l+n_h*Q_h,  n_l,n_h in {0,+1,-1}.   (G11)

They are not both zero and are constant on U. Their constancy can
also be seen directly: all selected endpoint labels are fixed, so
the unique tree path is fixed. It is independent of m.

If n_h is nonzero, the uniform height gap(G6) makes the second
term dominate the first, and(G10) excludes zeros for all large m.
If n_h=0, then n_l=+/-1 and(G7) excludes zeros. This proves the
lemma across the wall, without assuming a jump formula or
nonvanishing of a particular coefficient on both sides.

## 3. The entire nonnegative real boundary, including zero

**Lemma NonnegativeRealOnePacketNeighborhood.**
For every real w_0>=0 with w_0!=w_c, the actual I_m is nonzero
throughout some complex neighborhood of w_0 for every sufficiently
large integer m.

For w_0>0 the stationary map

    Omega(S)=S^d*(D-S)/(S-a)

is a decreasing bijection from(a,D) to positive w, with its sole
vanishing derivative at S=S_c, corresponding to w=w_c. The real
phase has exactly one positive stationary point. It is a
nondegenerate strict global maximum when w_0!=w_c, since the
positive integrand tends to zero at both real ends and the
stationary point is simple. For w_0=0,

    psi_0(S)=S^D*exp(-S)

has its unique nondegenerate global maximum at S=D. Thus this
case has the same one-packet property; no lower root-modulus
bound is required to exclude escape to w=0.

Here is a uniform complex-neighborhood proof, including the
noncompact parts of the actual cycle. On the whole positive ray,

    |psi_w(S)-psi_(w_0)(S)|
      <=|w-w_0|*S^a*exp(-S).                    (G12)

The right-hand factor has finite supremum and integral. Also
integral_0^infinity |psi_w(S)|dS is uniformly bounded for w near
w_0. Select a small real interval J about the unique maximum,
so that its complement has a strict supremum gap from the
maximum. This gap persists for complex w near w_0 by(G12).
The critical point and its local analytic Morse arc continue
holomorphically; near this interval psi_w is nonzero and a
logarithm exists even at w_0=0. Deform only J to that Gaussian
arc, adding short endpoint connectors which retain a strict
gap. The leading coefficient is analytic and nonzero. The
complement of J is bounded exponentially by(G9), now using the
uniform L1 bound on the full positive ray. It follows that

    I_m(w)=lambda(w)^m*sqrt(2*pi/m)*c(w)
                        *(1+O(1/m))

uniformly on a smaller complex neighborhood. This is a complete
period estimate and proves the lemma. It uses only local
logarithms; possible zeros of psi_w near the original endpoint
when w is close to0 do not obstruct the entire-integer-integrand
deformation or the direct absolute tail bound.


## 4. Negative real boundary: two complete packets, not zero-freeness

**Lemma NegativeRealPairNeighborhood.**
For every real w_0<0 there is a complex neighborhood U and an
unevaluated N_0 such that A_N,A_(N-1) have no common root x with
x/N^d in U for N>=N_0.

Write zeta_0=(-w_0)^(1/d)>0. On each of the real intervals
(0,zeta_0) and(zeta_0,infinity), |psi_(w_0)| has exactly one
nondegenerate maximum. Indeed Omega is a strictly decreasing
bijection(0,a) to(-infinity,0) and(D,infinity) to(-infinity,0).
These give the two positive stationary points p_L in(0,a) and
p_R>D. On the left interval psi is negative, on the right it
is positive; both maximum absolute values are nonzero.

Choose the analytic branch zeta(w) near w_0 and the exact endpoint
periods

    E_L(m,w)=integral_0^zeta(w) psi_w(S)^m dS,
    E_R(m,w)=integral_zeta(w)^infinity psi_w(S)^m dS.
    I_m=E_L+E_R.                                  (G13)

Use the original positive ray except for a finite moving path
to zeta. Because the integrand is entire for integer m, these
periods are independent of their finite path choices and
holomorphic in w. The infinite part stays on a fixed rightward
ray after a compact connector.

Both real maxima have local analytic Morse coordinates. On
the remaining parts of their respective real intervals there
is a strict maximum-modulus gap, including at all finite ends.
The same argument as(G12), on fixed paths plus short moving
endpoint connectors, preserves those gaps for w near w_0.
For the infinite tail the uniform L1 bound is available.
Thus deform near each maximum only, obtaining two complete
uniform expansions

    E_i(m,w)=lambda_i(w)^m*sqrt(2*pi/m)*c_i(w)
                          *(1+O_U(1/m)), i=L,R.    (G14)

The factors lambda_i and c_i are nonzero and analytic on a
small simply connected U. No logarithm crossing the zero
zeta is needed: the Gaussian logarithms are local near the
two separate maxima, and the remaining pieces are bounded
directly in absolute value.

The two stationary branches are analytic on U and
Re(p_R-p_L)>0 there after shrinking U. The envelope identity
(G5) holds on U as an identity of local logarithmic derivatives.
Choose local analytic logarithms and set

    D_crit=log(lambda_R/lambda_L),
    B=log|c_R/c_L|.

Even though the factors have opposite real signs at w_0,
a logarithm exists on a small complex disc because neither
vanishes. At an actual zero of I_m, (G14) implies

    m Re(D_crit(w))+B(w)=O_U(1/m).                (G15)

The exact adjacent-order comparison in Section5 therefore
applies on a smaller disc. It proves pair noncollision.

It would be false to replace this conclusion by uniform
zero-freeness of each I_m near the negative equimodular
crossing. Two packets may cancel there; this lemma concerns
their exclusion at the two different consecutive orders
and own-order scaled parameters.

## 5. The exact adjacent-order comparison and finite covering

We give the argument used both on a regular non-Stokes cell
and in Section4. Let K_0 be a compact small neighborhood
where a two-packet complete expansion holds and

    Re(S_R-S_L)>=c>0.

In a regular flow chamber the two oriented coefficients may
be0 or +/-1. If only one is nonzero, its Gaussian expansion
already excludes zeros. Otherwise their absolute values
are1 and (G15) holds; a minus sign merely changes an imaginary
logarithmic lift and does not affect(G15).

Let n=N-1 and suppose x is a common root. Set w=x/N^d.
By(G1),

    I_N(w)=I_n((N/n)^d*w)=0.

Work in a slightly larger compact collar and put

    m=N-theta,
    w_theta=(N/(N-theta))^d*w, 0<=theta<=1.         (G16)

For N large the full path stays in this collar. Equation(G5)
gives the exact identities

    d D_crit(w_theta)/dtheta=-(S_R-S_L)/m,
    d[m Re(D_crit(w_theta))]/dtheta
              =-Re(D_crit(w_theta))-Re(S_R-S_L). (G17)

All leading amplitudes are bounded and bounded away from0,
with bounded derivatives. The upper zero and(G15) give
Re(D_crit(w))=O(1/N). The first identity in(G17) keeps this
estimate uniform on the short path. Also

    d B(w_theta)/dtheta=O(1/N).

Therefore, writing F_m=m Re(D_crit)+B and integrating(G17),

    F_n(w_1)-F_N(w)<=-c+O(1/N),

whereas the two actual zeros give F_n(w_1)=O(1/N) and
F_N(w)=O(1/N). This is impossible for large N. Only values
of uniform Gaussian errors were used; no unspecified error
term was differentiated.

We now prove the main compact-complement theorem. At each
upper-half-plane point there are two possibilities:
if it is non-Stokes, the imported regular-cell theorem and
the preceding comparison give a neighborhood with eventual
pair noncollision; if it is a possible Stokes point,
Section2 gives a neighborhood where I_N itself is nonzero
for every sufficiently large N. Equal-magnitude points
cannot be Stokes by the explicit logarithmic argument in
Section1. Lower-half-plane neighborhoods follow by real
coefficients and conjugation. Sections3 and4 give all real
boundary neighborhoods except the positive cusp itself;
zero is included in Section3.

For each point of K choose a smaller neighborhood whose
closure lies in the corresponding controlled neighborhood.
The eventual pair statement on these smaller sets allows
the radial path(G16) when needed. Finitely many smaller
neighborhoods cover K. The maximum of their finitely many
cutoffs proves CompactComplementCuspConsecutiveExclusion.

## 6. Root capture using only the retained global upper bound

Theorem thm:product-annulus in GAMMA_PRODUCT_ROOT_MODULUS.tex
proves, for T_(m,gamma,d)(t)=sum_k binom(m,k)*(gamma)_(dk)*t^k,
that every reciprocal root satisfies

    |1/t| <= (gamma+d(m-1))_d.

Apply it to gamma=a*m+1 and to A_m(x)=x^m*T_(m,am+1,d)(1/x).
Every root x of A_m is nonzero and

    |x| <= ((a+d)*m-d+1)_d <= ((a+d)*m)^d,
    |x/m^d| <= (a+d)^d.                          (G18)

No lower annulus or new all-gamma improvement is needed.
For any epsilon>0, apply the compact-complement theorem to

    K={w: |w|<=(a+d)^d, |w-w_c|>=epsilon}.

It follows that for all sufficiently large N every possible
common root x of A_N,A_(N-1) satisfies

    |x/N^d-w_c|<epsilon.                         (G19)

This is common-root localization in every fixed physical
neighborhood, not all-root localization and not localization
in a fixed or growing quartic scaled window.

As a separate conditional deduction, if a proved neighborhood
theorem supplies r_s>0 and N_s such that no common root has
|x/N^d-w_c|<r_s for N>=N_s, then choose epsilon=r_s in(G19).
Combining the two exclusions proves eventual consecutive
coprimality for this fixed integer cusp. The present note
does not itself prove that inner-neighborhood theorem.
No finite N_0 is evaluated, and no assertion about the
remaining smaller orders follows.

## 7. Exact scope and proof-source record

The Stokes construction above is a finite sublevel detour
proved directly from the imported actual inward flow and
two-handle topology. The cycle method, Gaussian asymptotics
and one-dimensional Morse framework are classical; no
priority claim is made. An external Stokes jump theorem is
not used as an unverified premise.

Necessary sources and locators:

* CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md: T1--T4 for finite-band
  compactness, inward flow, end caps and handles; T6 for
  uniform end velocities; T8--T10 for exact path periods
  and regular Gaussian packets; T12--T15 for the regular
  consecutive comparison and absence of a Stokes point
  on the equal-magnitude arc.
* CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md:
  the actual two distinguished sector branches and their
  strict real-part separation, imported by the first
  source; no all-stationary dominance inference.
* GAMMA_PRODUCT_ROOT_MODULUS.tex:
  thm:product-annulus, used only for its old global upper
  bound in(G18).

All parameters here are fixed integer cusps. The result
does not assert a uniform bound in s, an effective cutoff,
all-order fixed-cusp SFC, complete binomial SFC, or an
independent replay or formal verification of the analytic
proof. Canonical files and status records were not changed.

