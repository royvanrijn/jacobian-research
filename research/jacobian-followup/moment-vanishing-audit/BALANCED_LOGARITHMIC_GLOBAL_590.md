# Balanced logarithmic seeds: conic336, W-linear425 and global590

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

15 November 2026. Internal written improvement of the retained
[logarithmic conic proof](QUADRATIC_W_LOG_PERIODS.md) and
[rescaled singular cover](RESCALED_SINGULAR_GLOBAL_888.md).

| Class | Sufficient initial zero-pure-moment gate | Current least-prefix interval |
|---|---:|---:|
| Every ternary polynomial of degree at most three | **590** | **[20,590]** |
| All homogeneous ternary cubics | **336** | **[14,336]** |
| All sixteen W-linear cubic coefficients | **425** | Optimal value not computed |
| Nondegenerate quadratic W coefficient, cubic G | **336**, a pure exclusion | Optimal value not computed |

The first590 zero pure moments force all-order pure nullity and the
retained Gaussian frame, with sharp mixed cutoff m>3*deg Q. No separate
frame condition is imposed on the global590 test. The odd-only interval
is[20,590]. Smooth-leading and homogeneous-smooth intervals remain[20,94]
and[14,56]. Both lower-bound coefficient witnesses remain complete with
first nonzero orders20 and14. No numerical leastness or20/14 sufficiency
is asserted.

The [checker](verify_balanced_logarithmic_seed_bounds.py),
[certificate](balanced_logarithmic_seed_bounds_certificate.json),
[receipt](balanced_logarithmic_seed_bounds_verification.json) and
[source audit](balanced_logarithmic_seed_bounds_source_audit.json) distinguish
exact algebra from the retained analytic inputs. The new all-coefficient
degree and class-cover arguments are written mathematics. Independent
replay, external review, formal verification and expanded global radical
cofactors remain outstanding.

## 1. Balance the conic coordinate without changing the logarithmic form

Write

    P=W*h(Z,T)+G(Z,T),
    h=a+bZ+cT+dZ^2+eZT+fT^2,
    e^2-4df!=0, deg G<=3.

The retained logarithmic proof gives normalized periods

    I(t)=constant*(1/t)*integral exp(Phi(t,x))*dx/x. (1)

Use y=t*x. Unlike an ordinary volume form, dx/x becomes exactly dy/y;
the prefactor1/t and its connection term -1/t stay unchanged.

In the f!=0 chart set eta=sqrt(e^2-4df), rho=(e+eta)/(2f),
sigma=(e-eta)/(2f), and

    ubar=((rho*c-b)*t+1)/eta,
    vbar=((b-sigma*c)*t-1)/eta,
    kbar=ubar*vbar-a*t^2/f.

Then the rescaled conic coordinates are

    Zhat=tZ=(f/eta)*(y-kbar/y-ubar+vbar),
    That=tT=y-ubar-rho*Zhat.                     (2)

They have y-orders in[-1,1], and their coefficient of y^k is a polynomial
in t of degree at most1-k. This includes the quadratic coefficient kbar
at the negative end; replacing all coordinate coefficients by a single
worst pole bound loses the useful asymmetry.

For f=0,e!=0 put z0=-c/e and A(z)=a+bz+dz^2. On the smooth component,

    Zhat=y+t*z0,
    That=(y-d*y^2-t*(b+2d*z0)*y+t*z0-t^2*A(z0))/(e*y). (3)

The same coefficient bounds hold. When a=c=f=0, the old formal branch
has g=0 and must not be discarded by inverting g. Its separate logarithmic
coordinate instead has

    That=y+(1-bt)/e,
    Phi=t*G(0,T)-T^2/2,

again satisfying the bounds. The old Laurent-unit checks remain valid
after multiplying the coordinate by t. Thus the trapped formal component
is included, and no singular t=0 Gaussian integral is introduced.

Let Psi(t,y)=Phi(t,y/t) and define

    W=t^2*Psi=t^3*G(Zhat/t,That/t)-That^2/2.

The coordinate bounds prove, for every coefficient point,

    W=sum_(k=-3)^3 A_k(t)*y^k,
    deg_t A_k<=3-k.                              (4)

Indeed a degree-m monomial of G contributes t^(3-m) times m coordinate
factors, with coefficient degree at most3-k. The Gaussian term satisfies
the stronger bound2-k. Cancellations can lower the actual Laurent orders;
all such endpoints are handled below.

The normalized source is

    Psi_t-1/t=S/t^3,
    S=t*W_t-2W-t^2,
    deg_t S_k<=3-k.                              (5)

The new source differs from the pulled-back old one by the exact
logarithmic twisted image -Dlog(1)/t, where

    Dlog(q)=y*q_y+(y*W_y/t^2)*q.                 (6)

There is no extra volume source. This contrasts with the ordinary binary
rescaling in the previous global888 proof, whose two-dimensional volume
factor was t^-2.

## 2. Keep the two endpoint denominators separate

Suppose the actual positive and negative orders p,q are both nonzero.
The full logarithmic quotient basis is y^-q,...,y^(p-1), of rank p+q.
This includes the regular-end and noncyclic-seed conventions from the
retained proof, not an ordinary Laurent ideal quotient.

Clear(6) by t^2. A high monomial is reduced using p*A_p; a low one using
-q*A_-q. Put a_+=3-p, a_-=3+q, the maximum endpoint coefficient degrees.
For an input of weighted degree M with wt(t)=wt(y)=1, one high reduction
cleared by A_p has degree at most M+a_+; one low reduction cleared by
A_-q has degree at most M+a_-. The Euler term t^2*y*q_y is one degree
lower than the phase-gradient contribution. These bounds remain valid
when M or a monomial exponent is negative.

A high reduction never produces an exponent below -q, and a low reduction
never exceeds p-1. Consequently their denominators do not alternate on
one reduction chain. For an input with support[-qj,pj], the necessary
endpoint depths are

    h_j=p(j-1)+1, ell_j=q(j-1).                  (7)

Using (A_p*A_-q)^max(h_j,ell_j) would be a valid but looser denominator.
Keeping A_p^h_j*A_-q^ell_j separately improves the result.

The actual seed derivatives have unreduced amplitudes R_j/t^(3j), where

    R_0=1,
    R_(j+1)=t^3*R_j'-3j*t^2*R_j+S*R_j.         (8)

Their support lies in[-qj,pj], and every coefficient of y^k has t-degree
at most3j-k. Reducing gives N_j/H_j with

    H_j=t^(3j)*A_p^h_j*A_-q^ell_j,
    deg_t N_(k,j)<=b_j-k,
    deg_t H_j<=b_j,
    b_j=3j+a_+*h_j+a_-*ell_j.                   (9)

The lower scalar coefficients are estimated at the actual seed rank,
not by imposing generic rank at special parameter points. A chosen
nonzero row minor includes the row k=0, since the initial seed is1.
If its row exponents sum to K0 and its rank is s, the zeroth polynomial
cofactor has degree at most

    sum_(j=1)^s b_j-K0.                         (10)

The higher scalar coefficients contain H_j and hence vanish to order
at least3j at t=0. If I=1+O(t^(L+1)), their derivative errors have order
at least L+1-j+3j>=L+1. Thus taking L equal to(10) contradicts the
nonzero zeroth coefficient. Nonvanishing uses the retained logarithmic
growth/period-comparison proof: a nonzero constant normalized period is
impossible. The algebraic bounds alone do not establish that input.

## 3. Regular ends and the exhaustive cubic bound

For q=0,p>0 the basis is1,...,y^(p-1), and h_j=p(j-1)+1,ell_j=0.
For p=0,q>0 it is1,y^-1,...,y^(-(q-1)), with h_j=0 and
ell_j=q(j-1)+1. The extra1 at the negative irregular end is essential:
the two-irregular-end basis had included y^-q, whereas this basis includes1.
For a constant phase the basis is1, H_j=t^(3j), and the gate is at most3.
Here the source cannot be zero: t*W_t-2W has zero coefficient at t^2,
while S includes -t^2.

All actual endpoints satisfy0<=p,q<=3. The worst bound over every seed
rank and every row subset containing0 is at full quotient rank. Adding
the next seed degree b_j>=3j exceeds any added row exponent, which is
at most2. Applying(7)--(10) and their one-end versions gives the complete
finite table:

| p \ q | 0 | 1 | 2 | 3 |
|---|---:|---:|---:|---:|
| 0 | 3 | 7 | 30 | 93 |
| 1 | 5 | 20 | 63 | 164 |
| 2 | 12 | 39 | 108 | 255 |
| 3 | 15 | 52 | 145 | **336** |

For p=q=3, b_j=21j-18 and the row exponents -3,-2,-1,0,1,2 sum to-3.
Hence21*21-18*6+3=336. This proves the nondegenerate conic exclusion
for cubic G, including all endpoint cancellations and the trapped branch.

The f=0 face therefore has sufficient gate max(336,320)=336, using the
retained fourteen-coefficient theorem when e=0. The other boundary
f!=0,e^2-4df=0 has the retained rank-one gate425. Thus the full W-linear
cubic gate is max(336,320,425)=425. Its sharp mixed cutoff3*deg Q is unchanged.
The older all-degree conic exclusions remain valid;336 is asserted here
only in the cubic coefficient scope, not for arbitrary degree G.

## 4. A homogeneous nodal complement bound129

For a homogeneous irreducible nodal cubic with nonisotropic node, use
P=Tuv+a*u^3+d*v^3 with a*d!=0 and the full homogeneous Gaussian quadratic.
The coordinate changes are linear, so its Gaussian means stay zero.
In the prior binary elimination/rescaling, D=m!=0 is now constant and
the rescaled quartic W is independent of t. Explicitly, before rescaling
its numerator is

    m*(t*(a*u^3+d*v^3)-q2(u,v))
         +(tuv-ju*u-jv*v)^2/2.

Substitute u=(x+jv)/t,v=(y+ju)/t and multiply by t^2. Every t cancels.
The transformed prefactor is still t^-2*m^-1/2, giving source
-2*m*t*W-2*m^2*t^3 of total degree at most5.

The four retained ordinary-node divergence projections, with constant
phase coefficients and delta=m*t^2, preserve total degree in t,x,y after
clearing their constant principal divisors. In particular the highest
polynomial terms in the fourth operator cancel as in its original identity;
one cannot estimate its two summands independently as ordinary products.
The actual R_j therefore reduce with coefficient degrees at most5j-d_i.
Their denominators have order4j. At seed rank s<=7, the zeroth scalar
degree is at most5s(s+1)/2 minus the selected basis grade sum. Its maximum
is5*28-11=129. The same retained nonconstant-period argument gives the
finite exclusion.

If the node is isotropic, homogeneity makes F2=0 and puts the polynomial
in the nondegenerate quadratic W-coefficient branch. That branch now has
gate336. Thus every homogeneous nodal cubic is excluded by336. Combining
smooth homogeneous56, irreducible cusp10 and homogeneous reducible10
gives a pure-only gate336 for *all* homogeneous ternary cubics. The lower
witness14 remains unchanged. An even336 bound is not reduced to335 by parity.

## 5. Updated global cover and assurance

| Leading type, arbitrary lower terms | Full-class sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | max(353,336)=353 |
| Irreducible cuspidal cubic | 558 |
| Transverse line-smooth-conic | max(334,336)=336 |
| Tangent line-smooth-conic | 590 |
| Three independent lines | max(304,336)=336 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The non-frame complements353/334/304 and weighted558/590 bounds are retained
from the preceding rescaled proof. Their frame branches now use conic336
or rank-one425 as appropriate. The maximum is590, proving the global
pure-only gate. The complete cubic theorem supplies its frame and mixed
consequences. In the rational twenty-coefficient ring,
radical(I590)=radical(I_infinity); no ordinary ideal equality or expanded
radical replay is claimed. At most ceil(log2(590-19))=10 exact monotone
search rounds would determine the least K in[20,590]. None has been run.

The checker covers three formal-coordinate charts, all16 endpoint cases,
nine exact Laurent models with37 connection and27 direct seed certificates,
45 constant-phase nodal reductions and three homogeneous normalization
identities. It preserves the initial Laurent-parser failure and the earlier
successful conic checkpoint. No high-order moment expansion is used.

The comparison source remains the retained Hien--Roucairol theorem; the
current version metadata was rechecked, without claiming a new body review
or importing a new theorem. The logarithmic growth proof, all-coefficient
degree arguments, full frame/class cover and analytic zeroth-coefficient
nonvanishing remain written and independently unreviewed. Leastness,
20/14 sufficiency, expanded global cofactors, exceptional schemes,
term-count minima and lower-power GVC remain open.
