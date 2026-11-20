# Strong Factorial at gap two for all sufficiently large lower exponents

The conclusion has an EXISTENTIAL threshold a_0; no numerical value, effective algorithm for a_0, or claim for every lower exponent is made. All coefficient-plane and actual-cycle hypotheses are stated below.

**Theorem LargeLowerExponentGapTwo.** There is an integer a_0 such that for every integer a>=a_0, every N>=2, and every characteristic-zero field, the ordinary factorial moments of X^a(A+B*X^2) at orders N-1 and N cannot both vanish unless A=B=0.

This theorem concerns the actual coupled shapes aN+1 and a(N-1)+1. It is not an assertion for a common fixed Gamma shape.

## 1. Actual normalization and global localization

For D=a+2 put g(s,w)=s^a(w+s^2) and

I_N(w)=integral_0^infinity g(s,w)^N exp(-Ns) ds,
J_N(w)=integral_0^infinity g(s,w)^(N-1) exp(-Ns) ds.

These are the actual consecutive factorial-moment polynomials up to nonzero factors and the common coefficient rescaling w=(A/B)/N^2. Monomial cases are immediate. The [retained linear lower-exponent criterion](GAMMA_PRODUCT_ROOT_MODULUS.tex) already proves the required noncollision when a>=2(2N-3). It remains to treat

N>a/4+3/2.                                                   (E1)

Let delta=2/sqrt(a), v=(s-a)/sqrt(a), w=a^2*(-1+delta*z). The exact phase, after removing a nonzero common factor, is

F_delta(v,z)=4*delta^(-2)*[log(1+delta*v/2)-delta*v/2]
             +log(z+v+delta*v^2/4).                          (E2)

It extends analytically to delta0 on every fixed bounded v region avoiding the displayed second logarithm's zero. Its limit is F_0=-v^2/2+log(z+v).

The independently proved [global Grace--Laguerre rectangle](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md) gives, for every root of the actual upper-order polynomial and sufficiently large a,

|Re z|<=8/sqrt(a)+18/a^(3/2), |Im z|<4.                     (E3)

The bound is uniform in N. Its proof uses the exact duplication into two binomial Grace convolutions and the two enclosing Laguerre matrix disks; it is not an asymptotic density statement.

The exact actual-order rescaling is

J_N(w)=(n/N)^(D*n+1)*I_n((N/n)^2*w), n=N-1.                (E4)

Thus the two endpoint parameters are z and

z_n=z+rho*(z-1/delta), rho=(N/n)^2-1>2/N.                  (E5)

On the hard region(E1), |z_n-z|=O(delta) uniformly for |z|<=5, while N,n tend uniformly to infinity with a.

## 2. Complete contours away from the turning points

Fix any eta in(0,1/10). There is a threshold a_eta such that all common roots in the hard region(E1) lie in

||Im z|-2|<eta.                                           (E6)

Here are the full contour facts behind this statement. In the limiting bulk 0<=y<=2-eta, the selected saddles at z=iy are

v_+/-=+/-sqrt(1-y^2/4)-iy/2.

Their real-part separation, positive horizontal Gaussian curvature, and the centre barrier are uniformly bounded away from zero on this compact y interval. On the shared horizontal line Im v=-y/2, the exact real profile is

H(t)=-t^2/2+y^2/8+(1/2)*log(t^2+y^2/4).

On each half-line this rises to its selected peak and then falls. The centre gap is -u/2-log(1-u)/2>0, u=1-y^2/4; a logarithmic zero at y0 only improves the gap. Uniform analytic perturbation for small delta and |Re z|=O(delta) preserves two local Gaussian packets, their nonlocal horizontal losses, and the short-connector gap. Deform the WHOLE original contour, cancel the common central vertical, and assign only the short remaining connector to the stronger packet. The real tails and the far verticals at Re v=+/-3 have uniform strictly smaller phase; the estimates are the same elementary derivatives used in [the retained full-bulk proof](COALESCING_HERMITE_BULK_CONTOUR_EXCLUSION.md). The resulting complete relative packet errors are O_eta(N^(-1/2))+O_eta(sqrt(N)*exp(-c_eta*N)), uniformly in small delta. No N*delta smallness is required: only the geometric estimates are perturbed, and the actual phase is retained in each exponent.

For 2+eta<=y<=4 choose the Hermite stable saddle -ib, b+1/b=y and b<1. Along its horizontal line the entire profile has

H''(t)<=-1+1/(y-b)^2=-1+b^2<0,

uniformly on this compact interval. The same complete contour with one horizontal and two far verticals gives a single nonzero Gaussian packet, so the individual upper moment cannot vanish there. Conjugation treats negative y.

Finally, in the bulk the exact action derivative along the actual radial interpolation is

d/dtheta[Delta_theta/(N-theta)]
 =-N*Re(S_theta,+-S_theta,-)/(N-theta)^2.

The uniform saddle separation in v gives an action drop at least c_eta*sqrt(a) between the actual orders. At an upper-order zero the within-order Gaussian prefactor ratio and small relative errors bound |Delta_0| by a fixed constant; the drop makes lower-order cancellation impossible for large a. This proves(E6).

All assertions in this section are uniform compact statements with positive margins. They can also be quantified by the explicit bulk/outside-endpoint guard packets. In the present existential theorem no uncalculated numerical threshold is substituted for a proof of positivity.

## 3. Jointly analytic endpoint normal form and its actual cycle

It suffices to treat the upper endpoint. At(delta,v,z)=(0,-i,2i),

F_v=F_vv=0, F_vvv=2i, F_vz=1.

The implicit-function determinant for(F_v,F_vv) in(v,z) is nonzero. Hence the critical point v_c(delta) and critical parameter z_c(delta) are jointly analytic near delta0. The [classical regular cubic theorem](gap_two_large_lower_endpoint_primary_audit.md) supplies jointly analytic coordinates

F_delta(v,z)=Phi(z,delta)+t^3/3-zeta(z,delta)*t,

with zeta(z_c(delta),delta)=0, v_t=i at delta0 and the critical point, and zeta_z=-i there. Shrink a fixed outer polydisc so that the transformation is one-to-one, its Jacobian stays nonzero, and

|zeta_z+i|<1/10.                                         (E7)

The extra delta parameter causes no new uniformity assumption: the phase and all relevant derivatives are bounded holomorphic functions on one fixed closed polydisc, and the cubic theorem's local inverse construction uses those same bounds. Equivalently, apply its proof with z-z_c(delta) as unfolding parameter; all nonzero derivatives have uniform positive lower bounds after shrinking the delta disk.

The actual cycle is the ordinary Ai cycle, up to one nonzero orientation factor. This is determined explicitly at the limiting critical point. Put v=u-i with u real. Then

F_0(u-i,2i)-Phi_0=i*u^3/3-u^4/4+O(u^5).

With v_t=i, the canonical coordinate satisfies t=-i*u+u^2/4+O(u^3). Thus the endpoints u=-r and u=+r of a sufficiently small fixed horizontal segment map into the upper and lower Ai DECAY sectors respectively, with positive real parts. Extending those endpoints within their decay sectors gives the standard Ai contour, with reversed orientation. The sign is a nonzero common factor and has no effect on zeros. No Bi or arbitrary Airy combination is introduced.

The full original contour has an actual nonlocal gap. At the limiting critical parameter its horizontal profile is

Re(F_0(u-i,2i)-Phi_0)=-u^2/2+(1/2)*log(1+u^2)
 <=-u^4/[4*(1+u^2)].                                    (E8)

For |u|>=r on the compact central part this is at most -c_r<0. The far verticals and original real tails are still strictly lower. These inequalities persist uniformly for delta and z-z_c(delta) in a sufficiently small fixed polydisc. The canonical endpoints remain in the same two decay sectors and may be extended radially after shrinking |zeta|; the tails of these extensions are exponentially small. This proves both the cycle identification and a uniform absolute nonlocal bound C*exp(m*(Re Phi-c)) for fixed c>0.

## 4. Uniform Airy vector remainder

The regular cubic transformation is the classical input. The first-order remainder needed here follows from one exact integration by parts.

First record an elementary Airy-contour envelope fact with contours compatible with a fixed local analytic disk. For every complex xi, there is a contour representing the standard Ai cycle whose absolute exponential integral is at most

C*(|Ai(xi)|+|Ai-prime(xi)|).                               (E9)

Here phi(tau)=tau^3/3-xi*tau. We give explicit constructions, so no uniform compact-thimble assertion is needed.

In the sector |arg xi|<=5*pi/6, let r=|xi|>=1, s=sqrt(xi)=p+iq, so p>sqrt(r)/4 and |q|<=sqrt(r). Set Y=6*sqrt(r), and parameterize the contour by

tau=s+x(y)+iy, y real,
x(y)=0 for |y|<=Y, x(y)=(|y|-Y)/sqrt(3) for |y|>=Y.

The exact phase difference is

Re[phi(s+x+iy)-phi(s)]
 =p*(x^2-y^2)-2*q*x*y+x^3/3-x*y^2.

On the two bends, 0<=x<=|y|/sqrt(3) and |y|>=6*sqrt(r). Hence

p*(x^2-y^2)<=-(2/3)*p*y^2,
x^3/3-x*y^2<=-(8/9)*x*y^2,
|-2*q*x*y|<=(1/3)*x*y^2.

The loss is at least (2/3)*p*y^2+(5/9)*x*y^2. On the central vertical it is exactly p*y^2. Since the arc-length derivative is at most2/sqrt(3), the whole absolute integral is at most

sqrt(2*pi/p)*exp[-Re(2*xi^(3/2)/3)]
 <6*r^(-1/4)*exp[-Re(2*xi^(3/2)/3)].

The arms asymptote to the Ai rays +/-pi/3 and, crucially, have cubic negative phase at infinity. For bounded xi use a fixed version of this contour through tau=1 with bends beginning at |y|=6; its absolute integral is uniformly bounded on a compact xi disk.

Near the negative axis write xi=-r*exp(i*theta), r>=1, |theta|<=pi/6, and scale tau=sqrt(r)*exp(i*theta/2)*u, lambda=r^(3/2). For psi(u)=u^3/3+u use the explicit upper curve

u=x(y)+iy, y>0,
x(y)=(y-1)*sqrt((y+2)/(3*y)).

A direct calculation gives

Im psi=2/3,
-Re psi=D(y)=2*(y-1)^2*(2*y+1)^2/(9*y)*sqrt((y+2)/(3*y)).

It goes from the left infinity sector as y tends to0, through u=i at y1, to the upper Ai infinity sector as y tends to infinity. The lower conjugate curve has Im psi=-2/3 and the same D. Multiplication by exp(3i*theta/2) preserves a loss at least D/sqrt(2), apart from the constant upper/lower saddle factors, since cos(3theta/2)>=1/sqrt(2).

The resulting one-dimensional absolute integral is bounded explicitly by50*lambda^(-1/2) times its saddle exponential factor. For completeness, on y in[1/2,3/2], D>=4*(y-1)^2/9 and |du/dy|<4. On y>=3/2, D>=2*y^2/27 and |du/dy|<3. On y<=1/2, substitute t=y^(-1/2): D>=t^2/18 and |du|<=4*dt. Gaussian integration bounds the three constants respectively by15, 27/2 and21, whose sum is less than50. After returning to tau, the two curves give total envelope

50*r^(-1/4)*(exp(Re Xi)+exp(-Re Xi)), Xi=2*xi^(3/2)/3.

Their common left-infinite pieces cancel homologically; their other ends have angles theta/2+pi/3 and theta/2-pi/3, uniformly inside the correct Ai decay sectors.

The uniform single-exponential expansions [DLMF9.7.5--6](https://dlmf.nist.gov/9.7#ii) compare the first contour's bound with the right side of(E9). The cosine/sine two-vector [expansions9.7.9--10](https://dlmf.nist.gov/9.7#ii) give the corresponding comparison near the negative axis: their leading matrix is invertible and |cos w|+|sin w| dominates half the sum of the two exponential moduli. On bounded xi the joint norm has a positive minimum, because Ai and Ai-prime have no common zero. These facts prove(E9) on the whole plane.

The constructions also settle a necessary truncation issue. In the physical cubic coordinate t=m^(-1/3)*tau, the positive-sector bends begin at distance O(sqrt(|zeta|)), or at O(m^(-1/3)) in the bounded-xi case. The negative-sector curves similarly enter fixed cubic decay subsectors once |t| is large compared with sqrt(|zeta|). Choose sqrt(|zeta|) much smaller than the fixed analytic-disk radius, and then m large. At its fixed outer truncation radius all arms therefore have uniformly negative CUBIC phase, not merely the vanishing quadratic loss of an unbent vertical. Connecting arcs can be drawn inside those fixed decay subsectors. This is what permits the contour envelope to be used for a locally analytic divided amplitude.

The same argument gives the weaker global lower bound

|Ai(xi)|+|Ai-prime(xi)|
 >=c_1*(1+|xi|)^(-1/4)*exp[-(2/3)*|xi|^(3/2)].             (E10)

Now let h(t,z,delta)=v_t be the analytic transformed amplitude on a fixed closed t disk. Uniform analyticity gives a constant M. Write exactly

h=A(z,delta)+B(z,delta)*t+(t^2-zeta)*G(t,z,delta).

The even/odd interpolation formulas at +/-sqrt(zeta) make A,B analytic in zeta, including zeta0. At the central critical point A=i; shrink the outer polydisc so |A|>=c_0>0 and |B|<=B_0, choosing B_0 positive by enlargement if needed. Indeed G(t)=(1/(2*pi*i))*integral_|u|=R h(u)/[(u-t)*(u^2-zeta)] du; differentiate under that fixed contour. On a smaller fixed disk this gives

|G_t|<=R*M/[(R-r)^2*(R^2-r^2)],                           (E11)

when the enclosing amplitude disk has radius R, the smaller one radius r, and |zeta|<=r^2. Thus G_t has a uniform bound independent of m, a and delta.

On the actual local cycle exact integration by parts yields

integral exp(m*phi)*h
 =A*integral exp(m*phi)+B*integral t*exp(m*phi)
  +(1/m)*[G*exp(m*phi)]_boundary
  -(1/m)*integral G_t*exp(m*phi), phi=t^3/3-zeta*t.          (E12)

Choose the local portion of the canonical descent cycle inside the smaller analytic disk. Here is why that replacement has only exponentially small connecting pieces. The original two endpoints and the corresponding canonical outer points lie in the same upper/lower cubic decay sectors; within either sector join them on a fixed outer annulus. After shrinking |zeta|, the cubic term dominates on those arcs and their real phase has a fixed negative bound. Near the negative xi axis the two descent pieces may also reach the third, left decay sector. Connect their two truncated left endpoints on that outer annulus: their common infinite tails cancel exactly, and the connector again has a negative cubic bound. All contour changes remain inside the analytic amplitude disk. Thus only the local portions carry the bounded divided amplitude; it is never assumed to extend to infinity. The added arcs belong to the already exponentially small boundary contributions. The bounded G_t and(E9), with tau=m^(1/3)*t, bound the last integral by

C*m^(-1/3)*(|Ai(xi)|+|Ai-prime(xi)|), xi=m^(2/3)*zeta.

The two complete canonical basis integrals are 2*pi*i*m^(-1/3)*Ai(xi) and -2*pi*i*m^(-2/3)*Ai-prime(xi). Their finite-endpoint differences and the explicit boundary term in(E12) belong to the same nonlocal exponential estimate from Section3.

Shrink the zeta disk once more so (2/3)*|zeta|^(3/2)<c/2, with c the nonlocal phase loss. The lower bound(E10) shows that every nonlocal exponential term, even after multiplying by the polynomial m factors, is at most C/m times the joint Airy norm for all sufficiently large m. This step explicitly handles sectors where Ai itself is exponentially small.

After dividing by the nonzero orientation and common factors, the actual period therefore has the exact representation

P_m=A*Ai(xi)+m^(-1/3)*B_*Ai-prime(xi)+R_m,
|R_m|<=C_R/m*(|Ai(xi)|+|Ai-prime(xi)|),                    (E13)

uniformly on one fixed endpoint polydisc for all m>=m_0; choose C_R positive by enlargement if needed. Here B_* differs from B by an irrelevant fixed sign. This derivation supplies the required remainder power directly. The explicit compatible cubic contours above supply its uniform contour bound.

## 5. Uniform endpoint zero band

For m>=2*C_R/c_0, a zero of(E13) cannot have Ai-prime(xi)=0. Dividing by that nonzero quantity gives

|Ai(xi)/Ai-prime(xi)|<=K*m^(-1/3),
K=2*(B_0+C_R)/c_0.                                      (E14)

The Airy--Robin lemma is particularly useful here. If Ai(xi)+lambda*Ai-prime(xi)=0, |lambda^2*xi|<=1/4, and |lambda|^3<=1/16, continue the zero backward along t*lambda. It obeys

dxi/dt=-lambda/(1-t^2*lambda^2*xi).

The bootstrap |xi(t)-xi(1)|<=2*|lambda| keeps the denominator at least5/8, so continuation reaches t0 and ends at a [negative real Ai zero](https://dlmf.nist.gov/9.9#i). Ai-prime never vanishes on this zero curve, by Airy-equation uniqueness. Thus |Im xi|<=2*|lambda|, uniformly in the Airy zero index. The escaping Robin root of scale lambda^(-2) is excluded by the explicit product condition.

Compute K on the outer polydisc before shrinking the actual endpoint window. Choose an inner zeta radius at most1/(8*K^2); then, for both actual orders lying in that window and m>=16*K^3, the ratio in(E14) satisfies both Robin conditions. We obtain the uniform actual zero band

|Im zeta(z,delta)|<=2*K/m.                                (E15)

There is no assertion that an individual Ai factor dominates throughout the window. The joint norm and the root continuation are essential.

## 6. The actual radial shift contradicts both zero bands

Choose the actual endpoint window with its closure strictly inside the INNER zero-band window from Section5. By(E1),(E5), the entire segment from z to z_n stays in that inner window once a is large enough, so both endpoints satisfy the Robin product guard as well as the outer uniform CFU estimates. Use the weaker bound |zeta_z+i|<=1/4, which follows from(E7), and |z|<=3. For delta<=1/15, the exact integral mean value formula gives

Im[zeta(z_n,delta)-zeta(z,delta)]
 >=rho*[(delta^(-1)-Re z)-(1/4)*|z-delta^(-1)|]
 >=rho*[(3/4)*delta^(-1)-15/4]
 >=rho/(2*delta)>1/(delta*N)=sqrt(a)/(2*N).                (E16)

If both actual moments vanished, their two bands(E15) would instead give

|Im[zeta(z_n,delta)-zeta(z,delta)]|
 <=2*K/N+2*K/n<=6*K/N.

This contradicts(E16) as soon as sqrt(a)>12*K. The argument holds for every hard N, including N much larger than a^(3/2); no amplitude Taylor approximation on an Airy scale is used.

## 7. Threshold order and exact conclusion

All choices are made in a noncircular order. First choose a fixed local cubic disk and a fixed outer parameter polydisc at the Hermite endpoint, with its nonlocal gap. Obtain the uniform analytic amplitude and first-order remainder constants there. Then shrink the inner endpoint window to satisfy the Robin product condition. For the two compact complementary y intervals, use Section2 with that smaller fixed window. Finally choose a_0 large enough for the global root rectangle, the required delta smallness, all Gaussian/CFU lower-order thresholds m_0 and16*K^3, the actual radial-segment containment, and sqrt(a)>12*K. The hard inequality N>a/4+3/2 makes every lower-order requirement on N,n automatic after enlarging a_0.

Every possible common root is now excluded: the large-lower-exponent regime by the retained exact theorem; the hard regime outside the turning windows by Section2 and the all-root rectangle; the two turning windows by Sections3--6 and conjugation. Therefore the actual rational moment polynomials are coprime over C for each integer a>=a_0 and every N>=2. Their nonzero rational resultants, or Bezout identities over Q, transfer the conclusion to every characteristic-zero field.

This establishes existence of a uniform lower-exponent threshold at fixed gap two. No explicit numerical a_0 is supplied here. The finite set of lower exponents below that threshold is not solved by this theorem, and no general-gap, fixed-shape-Gamma, or full binomial SFC claim follows.


## 8. Genuine threshold and existential endpoint equivalence

**Corollary EventualGapTwoEndpointEquivalence.**

For EACH integer A>=4*10^16 define the following assertions:

* T(A): for every integer a>=A and every N>=2, the actual ordinary
  gap-two consecutive moment polynomials are coprime.
* E(A): for every integer a>=A and HARD order N>a/4+3/2, the actual
  consecutive pair has no common complex zero in the SAME common-N
  endpoint windows |Re z|<=10^(-7),19/10<|Im z|<21/10.

Then T(A) holds if and only if E(A) holds. The forward implication is
restriction. For the reverse, use the pointwise hard-order endpoint
equivalence in
[the retained endpoint reduction](GAP_TWO_EXTENDED_BULK_ENDPOINT_REDUCTION.md).
For every a>=A the complementary soft orders are covered by the
separate linear lower-exponent criterion a>=2(2N-3), including equality.
This is the same common-N argument for BOTH actual polynomials; no
fixed-shape or independent lower-order scaling is substituted.

Consequently
\[
 \exists\,A\ge4\cdot10^{16}\ T(A)
 \quad\Longleftrightarrow\quad
 \exists\,A\ge4\cdot10^{16}\ E(A).
\]
Both existential assertions are proved by the main theorem: choose
A=max(a_0,4*10^16). Conversely any such E(A) implies T(A), hence the
existence of an ordinary gap-two all-terminal threshold.

The FIXED assertions T(4*10^16) and E(4*10^16) are not proved by this
existential argument; they are established separately in
[the effective gap-two theorem](GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md).
This argument itself gives no numerical bound, and the a>=40000
assertion remains open. No finite moment-order box for lower exponents
follows from the existential theorem.

## Exact reproduction and source scope

Run from the repository root:

    python3 research/factorial/image-rigidity/verify_gap_two_endpoint_contour_identities.py \
      --repo-root . --output /tmp/gap_two_endpoint_contour_review.json

The standard-library checker verifies universal polynomial identities
for the explicit Airy contours, stationary phase, divided amplitudes,
Robin flow and exact radial drift, together with the stated rational
guards and corruption controls. The analytic compactness and contour
arguments are the written proof, not a conclusion from finite fixtures.

The mandatory source ledger pins ALL fourteen listed proof/checker/
primary-context inputs, including the four actual cached NIST pages.
The [primary applicability audit](gap_two_large_lower_endpoint_primary_audit.md)
separates the imported regular cubic normal form from the first-order
remainder proved here. There is no claimed downloaded-byte hash for
the currently unavailable CFU full-PDF URL.

This existential theorem does not itself bound a_0 numerically. The
separate effective theorem supplies the a>=4*10^16 result; the broader
a>=40000 assertion remains open. Neither the existential theorem nor
that effective bound supplies a finite moment-order box for the
remaining lower exponents.

Updated 20 November 2026.
