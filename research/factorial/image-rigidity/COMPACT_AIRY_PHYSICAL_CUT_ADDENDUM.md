# The compact Airy representative also has certified physical cuts

This additive scratch lemma closes the geometric obligation for the OLD
compact contour used in EXPLICIT_AIRY_CONTOUR_ENVELOPE.md. It does not
replace that contour by the optimized contour through sqrt(xi).

Let m>=10^16, epsilon=m^(-1/3)<1/200000, |xi|<=4 and
zeta=epsilon^2*xi. Thus |zeta|<=4epsilon^2<10^(-10).
The old scaled bent vertical, in physical coordinates, is

    t=epsilon+x(y)+iy,
    x(y)=max(0,|y|-6epsilon)/sqrt(3).

Its full absolute integral has the proved C_A=2*10^7 joint-norm bound.
For |t|>=rho=3/8, the elementary estimate
|t|<=epsilon+(6/5)|y| gives |y|>3/10. Hence it is on a bent part.
Writing p=Re t gives

    p=|y|/sqrt(3)-(2sqrt(3)-1)*epsilon,
    |y|/2<p<|y|/sqrt(3).

For the lower bound use1/sqrt(3)>5/9, 2sqrt(3)-1<3, and
|y|/18>3epsilon. Therefore, with r=|t|,

    |y|>5r/6,
    -Re(t^3/3-zeta*t)
      >=(4/9)*|y|^3-|zeta|*r
      >=(125/486)*r^3-4epsilon^2*r
       >r^3/5.                                          (C1)

The last inequality follows from
(125/486-1/5)*rho^2>4/(200000)^2.
The real parts of the two canonical saddle phases have modulus at most
16/(3m), so C1 gives a loss greater than1/100 from either saddle
reference on this outer set. The exact guard is

    rho^3/5-16/(3*10^16)>1/100.

The two radius-.41 cuts are therefore in the right upper/lower
quadrants, with more than the1/20000 loss required for the certified
right-circle connectors. Radius is strictly increasing along each
outer bent tail: d(r^2)/d|y|=2p/sqrt(3)+2|y|>0. Thus there is exactly
one cut per tail, and the whole retained compact contour chain stays
inside the circle between those cuts.

For completeness the compact tail integrals also satisfy an explicit
bound. The derivative of radius is greater than5/6, while the arc
element is less than(6/5)d|y|. Thus |dt|<(3/2)dr. For j=0,1 the two
outer tails have absolute integral at most

    3*integral_rho^infinity r^j exp(-m*r^3/5)dr
       <=5*rho^(j-2)/m*exp(-m*rho^3/5).

Since the compact joint norm M(xi)>1/36000, each tail is less than

    2*10^6*m^(-1/3)*exp(-m/40000)*M(xi).

This includes all complex xi in the compact disk, so no sector switch
or Ai-alone lower bound is needed there. These are canonical tails;
an original period still requires the physical-disk replacement rule
and separately bounded original nonlocal pieces.

Research scratch20November2026; no canonical files changed.
