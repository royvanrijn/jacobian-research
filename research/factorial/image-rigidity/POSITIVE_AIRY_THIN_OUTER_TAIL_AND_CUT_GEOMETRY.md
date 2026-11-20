# Thin positive-sector Airy tails and the radius-0.41 cuts

This scratch result concerns canonical Airy contours and their cut geometry.
It does not itself bound an original period or supply a numerical Strong
Factorial threshold.

## 1. A compatible shallow bend and its loss

Let m>=1, |zeta|<=1/8, |Im zeta|<=2*10^(-5), and
|arg zeta|<=5*pi/6. Write s=sqrt(zeta)=p+i*q with p>0. The case zeta0
is included by continuity or by taking p=q=0 explicitly. We have

    p>sqrt(|zeta|)/4, |q|/p<4,
    q^2<=4*10^(-5), |q|<7/1000, |s|<177/500.            (P1)

The q bound follows from 2p|q|=|Im zeta|. Set Y=6|q| and choose

    t=s+x(y)+i*y, x(y)=max(0,|y|-Y)/sqrt(3), y real.

These are standard outgoing cubic decay directions. Their arc element
is at most (2/sqrt(3))*dy<(6/5)*dy. For phi=t^3/3-zeta*t,

    phi(t)-phi(s)=s*(x+i*y)^2+(x+i*y)^3/3.

The exact real loss is

    L=Re(phi(s)-phi(t))
     =p*(y^2-x^2)+2qxy+x*y^2-x^3/3.

On a bent part, |y|>=6|q| and x<=|y|/sqrt(3). The identity

    L-(2/3)*p*y^2-(5/9)*x*y^2
      =(p+x/3)*(y^2/3-x^2)+x*(y^2/3+2q*y)

therefore proves

    L>=(2/3)*p*y^2+(5/9)*x*y^2.                         (P2)

On the vertical part x0, (P2) holds directly. This bend depends on the
imaginary part of s, rather than its full modulus.

## 2. An absolute outer loss

The vertical part cannot reach rho=3/8, because there

    |t|^2<=p^2+49q^2<=1/8+48*(4/100000)
                         =3173/25000<9/64.

On each bent tail write a=|y|>=Y and y=epsilon*a. Its squared radius
has derivative

    d|t|^2/da=2*(p+x)/sqrt(3)+2*(a+epsilon*q)>0

away from the harmless origin when p=q=0. Thus there is exactly one
cut on each tail at every radius at least rho; no additional outer
components are omitted.

For every point with |t|>=rho,

    |y|>(5/6)*(rho-|s|).

If p>=1/4, this implies |y|>7/400. The first term of(P2) gives

    L>49/960000>1/20000.

If p<1/4, then |s|<251/1000, so |y|>31/300. Using Y<21/500 and
sqrt(3)<9/5 gives x>23/675. The second term of(P2) gives

    L>(5/9)*(23/675)*(31/300)^2>1/5000.

Thus, on the complete outer set,

    L>1/20000.                                          (P3)

## 3. Uniform unweighted and linear-weight tail integrals

Let T_(m,j) be the absolute integral over |t|>=3/8 with the extra
factor |t|^j, j=0 or1. Extract half of(P3), retaining half of(P2):

    T_(m,j)<=exp(m Re phi(s))*exp(-m/40000)
         * integral_contour |t|^j exp(-(m/2)*L) |dt|.      (P4)

Write xi=m^(2/3)*zeta and M(xi)=|Ai(xi)|+|Ai'(xi)|.

If |xi|>=4, the remaining Gaussian exp(-m*p*y^2/3) and the arc
bound imply

    I_0<9*m^(-1/2)*|zeta|^(-1/4),
    I_1<9*m^(-1/2)*|zeta|^(1/4)
                       +18*m^(-1)*|zeta|^(-1/2).

Here I_j denotes the last integral in(P4). These follow directly from
the whole-line Gaussian zeroth and first absolute moments; for I_1 use
|t|<=sqrt(|zeta|)+(6/5)|y|. The joint Airy lower bound from
EXPLICIT_AIRY_CONTOUR_ENVELOPE.md is

    M(xi)>=(3/(8sqrt(pi)))*|xi|^(-1/4)*exp(m Re phi(s)).

Since |zeta|>=4*m^(-2/3), both resulting normalized tails are less
than200*m^(-1/3)*exp(-m/40000).

If |xi|<=4, scale tau=m^(1/3)*t, S=sqrt(xi)=P+i*Q, and
w=m^(1/3)*y. The bend starts at6|Q|<=12. On |w|>=24, its real
increment X is at least |w|/(2sqrt(3))>|w|/4, so the remaining loss
is at least |w|^3/16. This exceeds |w|. On the middle interval discard
the nonnegative loss. The elementary estimates

    integral exp(-L_scaled/2)|d tau|<60,
    integral |tau|exp(-L_scaled/2)|d tau|<1900

follow as follows. The middle interval has width48; |tau|<=2+(6/5)|w|
there. On each outer interval majorize by exp(-|w|), whose zeroth and
first moments starting at24 are exp(-24) and25exp(-24). In particular
the weighted bound is below

    (6/5)*(48*(2+(6/5)*24)+2*(2+(6/5)*25))<1900.

Now |m Re phi(s)|<=16/3, exp(16/3)<225, and the compact joint lower
bound is M(xi)>1/36000. Thus each j=0,1 satisfies the single bound

    T_(m,j)<2*10^10*m^(-1/3)*exp(-m/40000)*M(xi).         (P5)

This is stronger than the requested bound with constant10^15. It
handles the canonical linear amplitude by multiplication by |A|+|B|.
It does not extend any original amplitude outside its inverse domain.

## 4. Cut positions and a circular-connector consequence

Set r_c=41/100. On the positive contour |t|=r_c forces
|y|>7/150>|q|, and it is beyond the vertical part. Therefore the upper
cut has Re t>0, Im t>0; the lower cut has Re t>0, Im t<0.

For the negative sector, write zeta=-q0*exp(i*theta), 0<q0<=1/8,
|theta|<=pi/6, and |Im zeta|<=2*10^(-5). Use exactly the explicit
negative curves of NEGATIVE_AIRY_CURVE_OUTER_TAIL.md and rotate by
theta/2. At radius r_c their incoming branch has y<3/8 and their
outgoing branch has y>1. Before rotation write the physical coordinates
as X=sqrt(q0)*x(y), Y=sqrt(q0)*y.

On the outgoing branch, (y-1)^2<=3x(y)^2, and the reverse triangle
inequality gives

    sqrt(R)-1<=|u-i|=sqrt(x^2+(y-1)^2)<=2x.

Consequently the following slightly weaker bounds hold:

    X>=sqrt(3)/4*(r_c-sqrt(q0)),
    Y>=sqrt(3)*r_c/2.

For the Y bound use R<=1+(4/3)*(y^2-1), which gives
y^2>=(3R+1)/4. The two negative-curve radial cuts are unique as well:
R is strictly decreasing for y<3/8 and strictly increasing for y>1,
and the entire set R>=9/8 lies in these two ranges.

Use cos(theta/2)>24/25 and |sin(theta/2)|<1/3. If q0>=1/1000,
the thin bound improves the latter to1/50. Then sqrt(q0)<177/500
gives Re t>71/5000>.01. If q0<=1/1000, use sqrt(q0)<4/125;
the same estimate gives an even larger positive margin. The upper
outgoing cut has Im t>(7/15)*r_c>0, and the lower cut the opposite
sign. Thus both outgoing cuts are strictly in their right quadrants.

On the incoming branch,

    |X|^2>=r_c^2-9/512>(19/50)^2,
    |Y|<(177/500)*(3/8)=531/4000.

After rotation both incoming cuts consequently satisfy Re t<-8/25.
There is no unproved assumption that these incoming branches lie in
the right half-plane.

Finally on the right circular quadrant |t|=r_c let p0=Re t>=0.
For REAL zeta the real phase is

    g(p0)=p0*((4/3)*p0^2-r_c^2-zeta).

Since g''(p0)=8p0>=0, any subarc between two points stays below the
larger endpoint phase. For complex zeta in the stated thin strip, the
total two-endpoint comparison costs at most

    2*r_c*|Im zeta|<=41/2500000=0.0000164.

Therefore endpoints lying below the SAME real reference phase by
1/20000 yield a full right-quadrant circular connector lying below
that reference by more than1/30000. This is a conditional geometric
statement: the endpoints must be certified relative to the same
reference phase. It does not certify any unspecified original
endpoint or close the original-period contour argument by itself.

The negative incoming cuts admit an analogous LEFT-circle connector.
On p0<=-8/25,

    g'(p0)=4p0^2-r_c^2-Re zeta
       >=4*(8/25)^2-(41/100)^2-1/8
        =233/2000>0.

Moving from either incoming endpoint towards -r_c decreases p0 and
therefore decreases g. The two arcs through -r_c stay in this range.
Restoring the imaginary part costs at most the same41/2500000, so
the certified negative-tail endpoint loss1/7000 leaves a whole-arc
loss greater than1/10000. If the two saddle reference phases differ,
use their maximum as the common reference; this only improves both
endpoint inequalities.

## 5. The physical-disk Cauchy replacement rule

Here is the precise rule for using these compact canonical chains.
Suppose an integrand f(v) is holomorphic on a physical disk D, a true
regular inverse v(t) takes the canonical disk |t|<.455 into D, and two
canonical endpoints t_-,t_+ are CERTIFIED to satisfy v(t_+)=v_+ and
v(t_-)=v_-. Let Gamma be any canonical chain joining t_- to t_+ inside
|t|<=.41, including the circular connectors above. Let gamma be any
physical chain in D with endpoints v_-,v_+. Then

    integral_gamma f(v) dv
      =integral_Gamma f(v(t))*v_t(t) dt.                  (P6)

Indeed f has a holomorphic primitive on the physical disk, and the
two integrals are the same difference of endpoint values. For the
actual local phase f(v)=exp(m*(F_delta(v,z)-Phi)), its pulled-back
integrand is exp(m*(t^3/3-zeta*t))*h(t). An additional holomorphic
physical amplitude is handled in the same way.

The rule requires the PHYSICAL chain to remain in D and the selected
canonical chain to remain in the inverse disk. It does not require
the canonical image of the entire physical chain to lie in the
inverse disk. In particular it supplies no permission to assume that
the image of an actual horizontal segment stays within that disk.
For the proposed application D is |v+z/2|<.94. The actual endpoint
matching must come from its separate certificate.

Reproduce the literal algebra/rational guards with
`python3 verify_positive_airy_thin_outer_tail.py`.
Scratch20November2026. No canonical files or status declarations changed.
