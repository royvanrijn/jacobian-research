# Uniform actual-phase, matched-data and complete-contour estimates

This account proves uniform estimates for every integer gap d>=2 at fixed small delta. The effective gap-two proof supplies numerical contour templates whose hypotheses are explicitly reverified for this actual family. Its d=2-scoped theorem is not simply applied to another gap.

Let d>=2 be an integer, a>0, delta=d/sqrt(a), epsilon=delta/d, and

h(v)=((1+delta*v/d)^d-1)/delta,
U(v)=d^2/delta^2*[log(1+delta*v/d)-delta*v/d],
F(v,z)=U(v)+log(z+h(v)).                                  (U1)

The removable value at delta0 is h=v,U=-v^2/2. Parameter bounds below allow complex delta unless actual original contours are specified.

## 1. Universal bounded-v majorants

**Lemma UniformActualPhaseMajorants.** For |v|<=R, t=|delta|R and q=t/d<1,

|h-v|<=(|delta|R^2/2)*exp(t),
|h'-1|<=|delta|R*exp(t),
|h''|<=|delta|exp(t), |h'''|<=|delta|^2 exp(t),
|partial_delta h|<=(R^2/2)*exp(t),                        (U2)

and

|U+v^2/2|<=|delta|R^3/[3d(1-q)],
|partial_delta U|<=R^3/[3d(1-q)^2],
U'=-v/(1+delta*v/d), U''=-(1+delta*v/d)^(-2),
|U'''|<=(2|delta|/d)/(1-q)^3.                             (U3)

Every bound is independent of d after replacing d by2 in the denominator. Proof: write

h=sum_(k=1)^d c_(d,k)*delta^(k-1)*v^k/k!,
c_(d,k)=product_(j=0)^(k-1)(1-j/d) in[0,1].

The exponential series gives(U2). In particular(k-1)/k!<=1/[2(k-2)!] gives the parameter derivative. More generally |partial_delta^ell h|<=R^(ell+1)*exp(t)/(ell+1). The logarithm series for U has coefficient1/(k*d^(k-2)); summing its absolute tail gives(U3). For its parameter derivative use(k-2)/k<=(k-2)/3 and the differentiated geometric series. No limit d fixed is taken.

The EXACT stationary numerator is

P(v)=z*v+(v-delta)*h(v)-1,
F_v=-P(v)/[(1+delta*v/d)*(z+h(v))].                       (U4)

Indeed h'=(1+delta*v/d)^(d-1) and(1+delta*v/d)^d=1+delta*h. These identities will also supply the exact coupled-order action derivative.

## 2. A uniform selected quadratic factor without a d=2 cubic formula

Assume |z-2i|<=11/100 and |delta|<=1/100000. Put u=v+i and P0=v^2+z*v-1=u^2+(z-2i)*v. On |u|=1, |P0|>=39/50 and(U2) gives

|P-P0|<=(2+|delta|)*2|delta|exp(2|delta|)+2|delta|
        <(601/100)*|delta|.                              (U5)

The same estimate at |u|=1/2, |v|<=3/2, is below(16/5)|delta|, whereas |P0|>=17/200. Rouche gives exactly TWO stationary roots inside |u|<1/2 and exactly two inside |u|<1. All other stationary roots lie outside the larger disk.

Let Q(v)=(v-v1)(v-v2) be their monic quadratic, with symmetric coefficients defined by the Cauchy power sums. Thus Q is jointly holomorphic in parameters, including at a double root, and P=Q*L with L analytic and nonzero on |u|<1. At delta0, Q=P0 and L=1.

Here is a quantitative factor bound which does not assume an explicit formula for the other d-1 stationary roots. Set

nu=(601/100)*(1/100000)/(39/50)<.000078.

On |u|=1, log(P/P0) is the convergent logarithm near0, bounded by nu/(1-nu). The boundary decomposition is

log(P/P0)=log L+log(Q/P0).

The second term is analytic outside the selected root disk and vanishes at infinity: Q and P0 are monic quadratics with both root pairs inside1/2. Its Laurent expansion contains only negative powers, so its interior Cauchy projection is0. Consequently

log L(u)=(1/(2*pi*i))*integral_(|s|=1)
                         log(P(s)/P0(s))/(s-u) ds.

This argument also constructs the jointly holomorphic branch of log L, normalized by0 at delta0. It is valid through the selected double root because Q is defined symmetrically. On |u|<=1/2,

|log L|<=2nu/(1-nu), |L'/L|<=4nu/(1-nu),
|L-1|<161/1000000, |L'/L|<32/100000.                     (U6)

Define K=-L/[(1+delta*v/d)*(z+h(v))], so F_v=Q*K. On this same disk let R0=61/100+(9/8)*10^(-5)*exp(3/(2*10^5)). Then |z+h(v)-i|<=R0<.610012. The multiplicative error in L/(1+delta*v/d) is below

(161/1000000+3/(4*10^5))/(1-3/(4*10^5)).

The disk for1/[1-i(z+h-i)] gives

Re(K/i)>3/5, |K|<13/5, |K'/K|<3.                        (U7)

For the last bound use(U6), |delta/d|/(1-|delta|*3/(2d)), and exp(3|delta|/2)/(1-R0). These are strict numerical inequalities, not an asymptotic claim.

## 3. Uniform matched data and the true cubic inverse

The selected-root construction from the quantitative endpoint-data proof now applies with exactly the shared constants(U7), on the OUTER complex disk |delta|<=10^(-5), |z-2i|<=.11. To spell out its symmetric construction, put c=(v1+v2)/2,q=(v1-v2)/2 and

J=(3/4)*integral_-1^1(1-s^2)K(c+q*s)ds, mu=J/i,
b=-i*q*mu^(1/3), zeta=b^2=-q^2*mu^(2/3).

The right-half-plane bound fixes the roots of mu without a branch ambiguity. The matched mean Phi is the average of the two critical values. The auxiliary H(v)=i*mu^(1/6)*(K(v)/i)^(-1/2) agrees with the TRUE inverse Jacobian at the two saddles only. Its interpolation supplies A,B; it is not used as the full Jacobian. From(U7), |q|<1/2 and |q/b|<5/4,

|H|<2, |H'|<3, |A|<2, |B|<4, |zeta|<.49.                 (U8)

These are exactly the outer bounds used by the retained quantitative endpoint-data proof. At delta0 the Hermite data are unchanged: on |z-2i|<=.101,

|A|>.989, |B|<500/989,
|zeta|<((118/117)^2)/9,
|zeta_z+i|<.03.

Taking ACTUAL |delta|<=10^(-10) reproduces the SAME Cauchy transfer ratio

|delta|/(10^(-5)-|delta|)<=1/99999.

Cauchy in z uses the same .009 outer margin. Hence all integer d>=2 satisfy

|A|>.98, |B|<.52, |zeta|<1/8, |zeta_z+i|<1/4             (U9)

on the actual endpoint rectangle |Re z|<=10^(-5),1.8999<=Im z<=2.1001. In particular the Hermite-to-actual zeta transfer is at most .49/99999, independently of the chosen physical endpoint.

The full-window regular inverse also persists with the SAME numerical radii and local amplitude bounds. Here are its quantitative inputs. On |u|=.94, v=u-z/2 has |v|<2 and |z+v|>=.00995. Equations(U2)--(U3) bound |partial_delta F|<203; at a selected critical point |v_j|<3/2, the identity

z+h(v_j)=(1+delta*v_j/d)^d/v_j

and Bernoulli's inequality give |partial_delta F(v_j)|<2.3, hence |partial_delta Phi|<2.3. Therefore the actual normalized phase differs from its Hermite value by less than206*10^(-10), which is smaller than the old certified110*10^(-8) allowance. The same Hermite outward-ball boundary .09, complex-z allowance.00052 and cubic radius.455 inequality therefore apply. The exact critical count does not use a degree-three stationary polynomial: the actual physical |u|<.94 disk is contained in |v+i|<.9901<1, and contains the selected pair because |v_j+i|<.5 implies |v_j+z/2|<.5501. The polynomial projection count in Section2 therefore gives exactly two physical critical points and no others. Both actual logarithm denominators are nonzero there, since |delta*v/d|<10^(-10) and |z+v|>=.00995 while |h-v|<3*10^(-10). The phase boundary is strictly above.0894, has three zeros by homotopy from the Hermite triple zero, and dominates the cubic on radius.455. This explicitly rechecks every numerical/count premise of the family-independent cover proof. The polynomial projection on |v+i|=1 never required the logarithm phase to be analytic on that larger circle; only P and P0 are used there. The written three-sheet/node/monodromy proof uses only these inputs and is unchanged. It proves the TRUE regular inverse on |t|<.455 with |u|<.94, h<38 on.425 and G<20000,G_t<1300000 on.41.

At the physical endpoints |u|=.4, the point derivative is below2.7, so the matched phase derivative is below5. Thus the actual W endpoint perturbation is at most2*10^(-5)+5*10^(-10), smaller than the exact ball certificate's2*10^(-5)+4*10^(-8). The zeta perturbation is no larger than its original allowance because the Cauchy transfer ratio and zeta_z bound in(U9) are identical. The complete branch-labelled endpoint cubic certificate and its decreasing radial connectors therefore apply to every d, with no sampling in d and no need to assume an auxiliary Jacobian equals the true one.

## 4. A uniform one-packet theorem at EVERY outside height

**Lemma UniformAllHeightOnePacket.** Let d>=2,a=d^2/delta^2 be positive integers, 0<delta<=10^(-10), m>10^12 an integer, |Re z|<=10^(-4), and |Im z|>=21/10. Then the ACTUAL normalized moment integral for(U1) is nonzero. There is no upper bound on |Im z| in this lemma.

Conjugation reduces to y=Im z>=2.1. Let b be the smaller positive solution of b+1/b=y. Then0<b<3/4 and1/b-b>7/12. On |v+i*b|=1/1000, the Hermite numerator P0=v^2+i*y*v-1 has modulus above(1/1000)*(7/12-1/1000). With R=751/1000, the actual perturbation is bounded by

10^(-4)*R+delta*[(R+delta)*R^2*exp(delta*R)/2+R].

It is strictly smaller than that boundary lower bound. Hence a unique simple stationary root V satisfies |V+i*b|<.001. This is uniform as y tends to infinity; no lower bound on b is used.

On the entire horizontal v=t+i*Im V, |t|<=3, |v|<16/5 and

Im(z+h(v))>2.1-.751-(delta/2)*(16/5)^2*exp(16delta/5)>4/3.

By(U2)--(U3),

Re F_vv<=-1+(2q+q^2)/(1-q)^2
              +delta*exp(16delta/5)/(4/3)
              +exp(32delta/5)/(4/3)^2<-2/5,
q=8delta/5.

At V this also gives |A|<8/5 for A=-F_vv(V). On the complex radius.05 around V the logarithmic denominator exceeds5/4, |v|<.801, and

|F_vvv|<=delta/(1-delta*.801/2)^3
 +delta^2*exp(delta*.801)/(5/4)
 +3*delta*exp(2delta*.801)/(5/4)^2
 +2*exp(3delta*.801)/(5/4)^3<6/5.

Thus the same Gaussian leading term, local |u|^3/5 remainder, and global horizontal decay u^2/5 are available uniformly in d and y.

Here is the necessary far-phase comparison, which replaces the old bounded-height assumption. Put S=Re F(V). The general estimate Re U(V)>-.283 and

|z+h(V)|>y-.752>(3/5)*y

give S>log y-1. At a real endpoint v=+/-3,

Re U(v)<-4.49,
|z+h(v)|/y<7/4,

so Re F(v)<log y-3.4. On the far verticals Re v=+/-3, |Im v|<.751,

Re U(v)<-4.2,
|z+h(v)|/y<1+(16/5)/y<8/3,

so Re F(v)<log y-3. Every far vertical lies at least2 below the stationary peak, independent of y. This comparison subtracts log y from BOTH phases; an absolute bound like S>0 is neither asserted nor needed at unbounded height.

The original entire integrand is

exp(-m*d*v/delta)*(1+delta*v/d)^(a*m)*(z+h(v))^m.

Only the finite original interval [-3,3] is deformed; the remaining original contour goes from -d/delta to infinity on the real axis. For v>=3, h>=v and

h'/h<=(1/v+delta)/(1+delta*v/d).

Since h/(h-|Re z|)<11/10,

Re F_v<=[-v+(11/10)*(1/v+delta)]/(1+delta*v/d)<-2.

For v=-r<=-3, down to the original endpoint, 0<=h'<=1 and -h is increasing. Also

-h(-3)>=(1-exp(-3delta))/delta>=3/(1+3delta).

Therefore Re F_v>=r-1/[3/(1+3delta)-10^(-4)]>2. Both unbounded/long tails are controlled by integration of these derivative bounds. The left original endpoint vanishes because a*m is a positive integer. No other stationary contribution is left unassigned by this whole-contour deformation.

The Gaussian leading modulus exceeds exp(m*S)/sqrt(m). The same exact Gaussian integrals as in the retained outside proof yield total relative error at most

[5+100*exp(-m/2000)]/sqrt(m)
 +6*sqrt(m)*exp(-m/2000)+3*sqrt(m)*exp(-m),

which is less than1/100 for m>10^12. In the far terms the new relative gap2 is stronger than the gap1 used in this display. A single leading packet cannot cancel with this error. This proves the lemma for all y>=2.1, not only for a bounded interval depending on d.

## 5. The actual central two-packet contour is uniform as well

Assume |Re z_theta|<=2*10^(-7), |Im z_theta|<=1.91 for the coupled radial family, and actual delta<=10^(-10). The Hermite pair has modulus1, real separation>74/125, and root disks of radius r=10^(-6). The general stationary perturbation on those disks is below

2*10^(-7)*(1+r)+(8/5)*delta,

which is smaller than r*(74/125-r). Both actual roots are within r of their respective Hermite roots. Their derivatives are nonzero, and their real separation exceeds.59.

The following explicit numeric bounds reverify the retained extended-bulk proof for the new phase. Set R=1+r, t=2*10^(-7)+r+(delta*R^2/2)*exp(delta*R). Then the curvature difference is bounded by

(2q+q^2)/(1-q)^2
 +delta*exp(delta*R)/(1-t)
 +(2e+e^2)/(1-t)^2+(2t+t^2)/(1-t)^2<3*10^(-6),
q=delta*R/2, e=delta*R*exp(delta*R).

Hence Re A>.17 and |A|<3. On the complex radius.05 about each actual saddle, the logarithmic denominator exceeds.94. The general third-derivative bound from(U2)--(U3) is below2.5, including the nonzero h''' term absent at d=2. The actual local horizontal decay is therefore at least3*u^2/50.

On the compact whole horizontal/connector region |v|<=16/5, the first-phase remainder is below6delta and |h-v|<6delta. These are literally the same coarse perturbation allowances used in the old complete horizontal comparison. The two actual heights differ by at most r from the common Hermite height. If the Hermite logarithmic denominator is at least.1, the real-phase difference is below1/50000; otherwise the actual phase is below-1.5 while each actual peak is above-1. The exact Hermite centre barrier is above1/600 and its profile outside local distance.049 has the old strict loss. Consequently the remaining horizontal loss is above1/10000 and the short connector loss above1/1000. The two common central verticals are cancelled BEFORE taking absolute values.

The far real tails have derivative bounds <-2 and>2 from Section4, now using |Re z|<=2*10^(-7) and -h(-3)>2.9 rather than a lower bound on Im z. At their endpoints the phase is below-2. On the far verticals |Im v|<.955001, |z+h(v)|<6 and the Gaussian comparison gives phase below-2, while each selected peak exceeds-1. All entire-contour nonlocal estimates are uniform in d.

The resulting complete packet errors are the old numerical bounds

[120+240*exp(-m/5000)]/sqrt(m),
3*sqrt(m/2)*exp(-m/10000),
4*sqrt(m/2)*exp(-m),
2*sqrt(m/2)*exp(-m/1000).

For m>10^12 each is below.001. Thus each complete packet has relative error below.01 and curvature prefactor ratio below sqrt(300/17)<5. These statements bound the ACTUAL whole contour, so extra stationary roots of the degree-(d+1) equation are not omitted.

The coupled action derivative is also EXACT and uniform in d. At a stationary point, v*(z+h)=1+delta*h implies

(z-delta^(-1))*F_z=1-v/delta.

Let D(z)=F(V_+(z),z)-F(V_-(z),z). Along the actual radial path

z_theta-delta^(-1)=(1-theta/N)^(-d)*(z_0-delta^(-1)),

stationarity cancels the root derivatives and gives

dD/dtheta=-d/[delta*(N-theta)]*(V_+-V_-).                (U10)

Define Delta_0=N*Re D(z_0), Delta_1=n*Re D(z_1). Since d/delta=sqrt(a), Re(V_+-V_-)>.59 and n*log(N/n)>5/6 for N>=6, integration gives Delta_1<(n/N)*Delta_0-.4*sqrt(a). At an upper-order zero the packet prefactor/error bounds give |Delta_0|<log6. The lower packet ratio is then below40*exp(-.4*sqrt(a))<1 under actual delta<=10^(-10), hence sqrt(a)>=2*10^10. Thus the actual adjacent pair has no common bulk zero.

## 6. Reverified original endpoint nonlocal bounds and the complete remainder

For the actual endpoint window, the original contour is shifted by exactly the same Hermite horizontal heights as before: D=Im z/2, s=0 for D<=1 and s=sqrt(D^2-1) otherwise. Its physical local endpoints have |u|=.4. The limiting normalized horizontal loss is at least.0007 below -kappa; the matched critical values have limiting real parts0,0 for D<=1 and +/-kappa for D>=1. These are facts about the common Hermite limit, independent of d.

Here are all phase-dependent perturbation steps needed to transfer the full numerical nonlocal proof, rather than applying its d=2-scoped statement. On its horizontal |v|<16/5, |z+v|>=18999/20000, equations(U2)--(U3) give |partial_delta F|<11. The matched mean has derivative below2.3. Thus the normalized phase perturbation is below

(11/5)*10^(-5)+14*10^(-10)<3*10^(-5).

On the physical |u|=.4 critical-count circle, |P0| has the old margin.16-.103 and |P-P0|<3.2delta. Exactly two selected roots lie there, so |q|<.4. At a stationary point

F_z=v/(1+delta*v/d)^d.

The v derivative of this last rational expression is

[1-(d-1)*(delta/d)*v]/(1+delta*v/d)^(d+1).

On the joining segment |v|<=1.5 its modulus is below1.001, uniformly in d. Hence each individual matched critical action has z derivative below.401. Its delta derivative is below4.6, from the two2.3 stationary/mean bounds. The actual critical-action perturbation is therefore below

.401*10^(-5)+4.6*10^(-10)<10^(-5).

This also describes the correct selected Puiseux continuation through the merger; the unordered critical pair and its mean are analytic, and either action can be continued piecewise with these bounded derivatives. No derivative of a singularly labelled root is inserted.

On the far vertical rectangle the point delta derivative is below12, rather than the old9. Since the ACTUAL delta is at most10^(-10), the strict bounds -2.5+12delta<-1.9 and -3+12delta<-2 still hold. The mean lower bound Re Phi>.4 follows with its2.3delta allowance. The uniform real-tail derivative estimates in Section4 handle both original tails. Thus every final numeric margin in the original nonlocal proof is unchanged: its horizontal loses.00067 relative to-kappa, the actual critical comparison costs at most.00001, and the far pieces have the same stronger losses.

Consequently, for EVERY d>=2 and every actual order m>=10^16, its normalized original nonlocal contribution still obeys

|R_nonlocal|<60000*sqrt(m)*exp(-33*m/50000)*M,
M=|Ai(m^(2/3)zeta)|+|Ai'(m^(2/3)zeta)|.                  (U11)

The proof of(U11) has now been supplied for the new actual family: it follows from the same explicit Hermite inequalities with the new uniform perturbation/entire-tail estimates just listed. It is not inferred merely because its formula matches a previous theorem.

Every other complete-period input depends only on the now-reverified matched data, true inverse, endpoint EW/EZ boxes and canonical cubic, not on d. In detail: |Im zeta|<2*10^(-5) follows from the real Hermite value and the identical transfer; the old compact contour and optimized positive/negative Airy contours have unchanged joint envelopes and weighted tails; the right/left circle connectors lie in the same certified sublevel ranges; the physical disk is simply connected and the local original chain is inside it; and h38,G20k,G_t1.3e6 hold on the same circles. The exact Cauchy replacement and one integration by parts therefore yield, for the ACTUAL normalized F-period,

P_m=A*Ai(xi)-m^(-1/3)B*Ai'(xi)+R_m,
|R_m|<=10^14*M/m if |xi|<=4,
|R_m|<=10^9*M/m if |xi|>=4, xi=m^(2/3)zeta.               (U12)

For clarity the numerical aggregation is the same five proved contributions: normalized original60000sqrt(m)exp(-.00066m), weighted canonical tails5*10^14 exp(-m/40000), leading connectors10^9 sqrt(m)exp(-m/40000), G_t connectors5*10^13 m^(-1/2)exp(-m/40000), and both physical G endpoints2*10^11 m^(-1/2)exp(-m/40000), all multiplied by M. Their sum is below10^16*m*exp(-m/100000)*M<M/m for m>=10^16. The local core term has coefficient1300000*C_A/(2pi), with C_A=2*10^7 or300; this proves exactly the two constants in(U12).

The exact Robin-flow product/cubic guards consequently give, at EVERY actual period zero in the endpoint window,

|Im zeta|<5200/m.                                       (U13)

All constants in(U11)--(U13) are independent of d. This is a new uniform family bridge for those complete-period estimates, proved by verifying their hypotheses; it is not a claim that the old uniform C_R=10^10 question has been resolved.

