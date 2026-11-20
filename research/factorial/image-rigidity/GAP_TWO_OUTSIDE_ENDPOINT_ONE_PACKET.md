# One-packet exclusion beyond the Hermite endpoint for the actual gap-two moments

This is a separate actual one-packet lemma. It leaves the compact turning-point window open.

Let a>=4*10^16 be an integer, epsilon=1/sqrt(a), delta=2/sqrt(a)<=10^(-8), and m>10^12 an integer. Set

w=a^2*(-1+delta*z), s=a+sqrt(a)*v,
h(v)=v+delta*v^2/4,
F(v,z)=epsilon^(-2)*(log(1+epsilon*v)-epsilon*v)+log(z+h(v)).

**Lemma OutsideEndpointOnePacket.** If |Re z|<=10^(-4) and 21/10<=|Im z|<=4, then the actual integral integral_0^infinity [s^a(w+s^2)]^m exp(-ms) ds is nonzero.

Proof. Conjugation reduces to z=x+iy with y>=21/10. Let b be the smaller positive solution of b+1/b=y. Then 1/4<b<3/4, because 1/4+4>4 and 3/4+4/3<21/10. The Hermite stationary polynomial at z=iy has roots -ib and -i/b, separated by more than7/12.

The exact stationary polynomial for the actual phase is

P(v,z)=v^2+z*v-1+delta*(v^3/4-v)-delta^2*v^2/4.

On |v+ib|=1/1000, the Hermite polynomial has modulus greater than (1/1000)*(7/12-1/1000). The remaining terms have modulus at most

(1/10000)*(751/1000)
 +10^(-8)*[(751/1000)^3/4+751/1000]
 +10^(-16)*(751/1000)^2/4,

which is strictly smaller. Rouche isolates a simple stable stationary branch V within1/1000 of -ib. Direct differentiation also bounds P_v away from zero on that disk, so the selected branch varies analytically over the indicated rectangle. In particular |V|<751/1000 and |Re V|<1/1000.

Write eta=Im V. On the ENTIRE horizontal segment v=t+i*eta, |t|<=3,

Im[z+h(v)]=y+eta*(1+delta*t/2)>4/3.

Put q=epsilon*(16/5). The identities h'=1+delta*v/2, h''=delta/2 and U''=-(1+epsilon*v)^(-2), for U=epsilon^(-2)(log(1+epsilon*v)-epsilon*v), give

Re F_vv<=-1+(2q+q^2)/(1-q)^2
          +(delta/2)/(4/3)+(1+(8/5)*delta)^2/(4/3)^2<-2/5.

Thus the actual phase, not merely its Hermite approximation, is strictly concave on the whole horizontal segment. At V its curvature A=-F_vv has Re A>2/5 and |A|<8/5. On the complex radius1/20 around V the logarithmic denominator exceeds5/4; the third derivative has modulus below6/5. The whole segment consequently satisfies

Re[F(V+u)-F(V)]<=-u^2/5.

The stationary peak S=Re F(V) is positive: Re(-V^2/2)>31/1000, the Taylor remainder of U+v^2/2 is below10^(-8), and |z+h(V)|>4/3.

Deform the entire original contour to real [-1/epsilon,-3], a vertical to -3+i*eta, the single horizontal to3+i*eta, a vertical back to3, and the real right tail. The integrand is entire, since a,m are integers. No uncancelled centre vertical occurs.

For the outer right real tail, Re z is bounded by10^(-4), so

Re F_v<=[-v+(10001/10000)*(1/v+delta)]/(1+epsilon*v)<-2,

using h>=v and h'/h<=(1/v+delta)/(1+epsilon*v). The endpoint phase is below -900/203+log6<-2. The left tail has h'<1, |h|>=591/200 and hence Re F_v>=3-1/(591/200-1/10000)>2; its endpoint phase is below -9/2+log6<-2. Each tail therefore has absolute size below exp(-2m)/(2m).

On the far verticals |Im v|<751/1000, |z+h(v)|<8, and

Re F<=-(9-(751/1000)^2)/2+6*10^(-8)+log8<-21/10.

Their lengths total less than2. They too are exponentially smaller than the selected peak S>0.

The Gaussian leading term has modulus greater than exp(mS)/sqrt(m). On the local radius1/20, comparison with the Gaussian along the real displacement gives absolute error at most 5*exp(mS)/m: the Taylor remainder is at most |u|^3/5, both interpolated real phases have decay at least u^2/5, and integral_R |u|^3 exp(-m*u^2/5) du=25/m^2. Truncating the Gaussian contributes at most100*exp(mS-m/2000)/m. The remaining horizontal length is below6 and has loss at least1/2000 from the peak. Therefore the complete relative error is bounded by

[5+100*exp(-m/2000)]/sqrt(m)
 +6*sqrt(m)*exp(-m/2000)+3*sqrt(m)*exp(-m),

which is less than1/100 for m>10^12. The one leading packet cannot cancel with this error. This proves the lemma.

Combining this statement with the independent all-root Grace rectangle bounds shrinks every actual gap-two moment zero, in the present huge-a range and hard terminal regime, to |Im z|<21/10. It still leaves the endpoint neighborhoods and gives no all-gap-two Strong Factorial theorem. The separate central whole-contour packet only excludes simultaneous zeros for |Im z|<=29/20.


The independently proved all-root result referred to above is
[the exact global gap-two Grace rectangle](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md).
The earlier central comparison is
[the complete bulk contour theorem](COALESCING_HERMITE_BULK_CONTOUR_EXCLUSION.md).
These references do not turn the remaining endpoint noncollision
assertion into a conclusion.

## Exact reproduction

    python3 research/factorial/image-rigidity/verify_gap_two_endpoint_reduction_guards.py \
      --repo-root . --output /tmp/gap_two_endpoint_reduction_review.json

The common supporting ledger also pins the distinct thin-strip proof
and every imported local/global/contour proof.

Updated 20 November 2026.
