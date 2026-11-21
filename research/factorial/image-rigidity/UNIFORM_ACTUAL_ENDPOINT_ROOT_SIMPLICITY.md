# Uniform simplicity inside the actual Airy endpoint window

This note proves the missing open-collar and derivative estimate for the
actual endpoint expansion in
research/factorial/image-rigidity/UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md,
(U9), (U12) and (U13). It concerns zeros in that endpoint window only.

**Theorem UniformActualEndpointRootSimplicity.** Let d>=2 and a>=1 be
integers, delta=d/sqrt(a)<=10^(-10), and m>=10^16 an integer. Define

    A_m(x)=sum_(j=0)^m binom(m,j)(am+1)_(dj) x^(m-j),
    x=m^d*a^d*(-1+delta*z).

Every zero of A_m in the closed upper endpoint rectangle

    |Re z|<=10^(-7),       1.9<=Im z<=2.1                 (E1)

is simple. The conjugate lower endpoint rectangle has the same conclusion.
The constants are independent of a,d,m and of the root index. The
statement does not assert that every root lies in these rectangles or
that roots outside them are simple.

## 1. The actual holomorphic expansion being used

For fixed a,d,m, the complete normalized period in (U12) is

    P_m(z)=A(z)*Ai(xi)-epsilon*B(z)*Ai'(xi)+R_m(z),
    xi=mu*zeta(z),   mu=m^(2/3),   epsilon=m^(-1/3),       (E2)

on the interior of the expanded physical rectangle

    Omega={|Re z|<10^(-5), 1.8999<Im z<2.1001}.           (E3)

This is one actual normalized period, not independent choices of a
period at different parameters. Its normalization is a nonzero constant
times exp(-m*Phi(z)); Phi is the analytic matched mean. The actual
original integral is a polynomial in z. The regular cubic inverse and
the symmetrically constructed matched A,B,zeta in (U6)--(U9) are
holomorphic on this physical domain. Consequently P_m and the remainder
defined by subtraction in (E2) are holomorphic there. Piecewise choices
of estimating contours do not alter this identity.

Write y=Ai and use the sign convention B_tilde=-B below. The source gives
throughout Omega

    |A|>49/50,  |B_tilde|<13/25,  |A|<2,
    |zeta|<1/8,  |zeta_z+i|<1/4,                        (E4)

and, with M=|y(xi)|+|y'(xi)|,

    |R_m|<=10^14*M/m     everywhere,
    |R_m|<=10^9*M/m      where |xi|>=4.                 (E5)

The first inequality follows also in the second region by enlarging its
constant. Thus no artificial discontinuity in the error is introduced
at |xi|=4. The thin |Im zeta| restriction used to prove (U12) holds
throughout the stated physical Omega, not just at actual roots.

For completeness the exact original-polynomial normalization is as
follows. With D=a+d and w=a^d*(-1+delta*z),

    I_m^orig(w)=int_0^infinity [S^a(w+S^d)]^m exp(-mS)dS
      =Gamma(am+1)/m^(Dm+1)*A_m(m^d*w).

The change S=a+sqrt(a)*v gives

    I_m^orig(w)
      =sqrt(a)*[a^(a+d)*delta*exp(-a)]^m I_m^v(z).

All displayed constant factors are nonzero. Multiplication by the
Airy normalization m^(1/3)exp(-m*Phi)/(2*pi*i), with either fixed
orientation, is also locally nonzero. This records the precise
multiplicity-preserving link from P_m to A_m.

## 2. A common open zeta collar, with amplitude derivative bounds

Fix a zero z0 in (E1). Its closed physical disk of radius

    r_z=1/125000=8*10^(-6)

is strictly inside Omega: its real boundary has modulus at most
8.1*10^(-6)<10^(-5), and its imaginary boundary is more than
9*10^(-5) from either expanded horizontal edge.

The derivative estimate in (E4), integrated along a segment in this
disk, gives

    |zeta(z)-zeta(z')+i(z-z')|<|z-z'|/4,
    |zeta(z)-zeta(z')|>3|z-z'|/4.                       (E6)

In particular zeta is one-to-one. For any
|tau-zeta(z0)|<3r_z/4, Rouche on the physical boundary, comparing with
-i(z-z0), gives exactly one solution of zeta(z)=tau. The inverse is
holomorphic and |z_tau|<4/3. Thus its image includes a common disk of
radius 6*10^(-6), and we reserve the smaller OPEN collar

    |tau-zeta(z0)|<3*10^(-6).                            (E7)

For tau in (E7), (E6) places its inverse within r_z/2 of z0. A physical
Cauchy circle of radius r_z/2 is available there. Using |A|<2 and
|B_tilde|<13/25 in (E4), and then |z_tau|<4/3, gives on all of (E7)

    |A_tau|<2*(2/r_z)*(4/3)<10^6,
    |B_tilde_tau|<(13/25)*(2/r_z)*(4/3)<10^6.             (E8)

The complete period and remainder in (E2) are holomorphic on (E7), and
(E4)--(E5) hold there. This is a common admissible complex collar, not
an estimate confined to a zero curve.

If one also uses (U13), the root itself satisfies
|Im zeta(z0)|<5200/m. The Airy-spacing disks used next therefore lie
very deeply inside the thin strip: their imaginary modulus is below
5200/m+1/(4mu)<10^(-10). The physical collar already proves their
admissibility, so (U13) is not needed to establish that collar.

## 3. The Airy vector on a whole root-spacing disk

Put xi0=mu*zeta(z0) and D0=sqrt(1+|xi0|). Since m>=10^16,

    epsilon<1/200000,        mu>4*10^10,
    D0*epsilon<3/8,          D0/m<10^(-11).             (E9)

Indeed D0^2*epsilon^2=m^(-2/3)+|zeta(z0)| and
D0/m<=1/m+sqrt(1/8)*m^(-2/3). The elementary strict bounds
sqrt(1/8)<3/8 and mu>4*10^10 give (E9).

At P_m(z0)=0, (E4)--(E5) show that y'(xi0) is nonzero. Otherwise the
zero equation would force y(xi0)=0, contradicting uniqueness for the
nontrivial Airy ODE solution. The same equation gives, with
t0=|y(xi0)/y'(xi0)|,

    t0<[(13/25)*epsilon+C/m]/[49/50-C/m].                (E10)

For all xi0 one may take C=10^14, whence t0<1/90. If |xi0|<=5,
D0<5/2 and hence D0*t0<1/4. If |xi0|>5, use C=10^9 at the root
instead. By (E9),

    D0*t0<
       [(13/25)*(3/8)+10^9*10^(-11)]
       /[49/50-10^(-7)] <1/4.                           (E11)

On the entire complex s disk |s|<=1/4 define

    Y(s)=D0*y(xi0+s/D0)/y'(xi0),
    V(s)=y'(xi0+s/D0)/y'(xi0).

They satisfy

    Y'=V,       V'=[(xi0+s/D0)/D0^2]*Y.

The coefficient has modulus at most5/4. Their initial l1 norm is less
than5/4, by (E11) and its compact case. Gronwall on each radial segment
therefore gives

    |Y(s)|+|V(s)|
       <(5/4)*exp(5/16)<20/11<2.

Since D0>=1, throughout this ENTIRE disk,

    |y(xi)|+|y'(xi)|<2*|y'(xi0)|.                       (E12)

No all-index Airy zero asymptotics or separation theorem is used.

The corresponding tau disk has radius

    rho=1/(4*mu*D0)<=1/(4mu)<1/(16*10^10).

It lies in the common physical collar (E7). If |xi0|<=5 use C=10^14
on its boundary. If |xi0|>5, the boundary has |xi|>5-1/4>4, so C=10^9
is valid throughout it. Cauchy's derivative estimate and (E12) imply

    |(R_m)_tau(zeta(z0))|/[mu*|y'(xi0)|]
       <=8*C*D0/m <1/5.                                (E13)

For the compact case use D0<5/2 and m>=10^16. In the other case
8*C*D0/m<8/100<1/5. This supplies the derivative of the actual
remainder from an open-domain estimate, rather than differentiating
a pointwise big-O statement.

## 4. The differentiated period cannot vanish

For notational simplicity write B for B_tilde in this section. From
the exact zero equation and y''(xi)=xi*y(xi), differentiation of (E2)
gives at tau0=zeta(z0)

    (P_m)_tau
      =mu*[A-(B^2/A)*tau0]*y'(xi0)
        -mu*epsilon*(B*xi0/A)*R_m
        +A_tau*y(xi0)+epsilon*B_tau*y'(xi0)+(R_m)_tau.    (E14)

The relation epsilon^2*xi0=tau0 is exact. The leading coefficient has
the uniform lower bound

    |A-(B^2/A)*tau0|
      >49/50-(13/25)^2/[8*(49/50)]
       =4633/4900>47/50.                                (E15)

The normalized amplitude-derivative error is, by (E8)--(E10),

    |A_tau*y+epsilon*B_tau*y'|/[mu*|y'|]
       <10^6/(90*4*10^10)+10^6/10^16<10^(-6).            (E16)

For the mixed remainder term in (E14), (E10) gives
|R_m|/|y'(xi0)|< (C/m)*(91/90). If |xi0|<=5, its normalized modulus
is less than

    (1/200000)*(26/49)*5*(1/100)*(91/90)<1/400.

If |xi0|>5, use C=10^9, |tau0|<1/8 and
epsilon*|xi0|/m=|tau0|/mu to get the same bound from

    (26/49)*(10^9/[8*(4*10^10)])*(91/90)<1/400.           (E17)

Combining (E13)--(E17) proves the explicit nonvanishing estimate

    |(P_m)_tau(tau0)|
      >[47/50-1/5-1/400-10^(-6)]*mu*|Ai'(xi0)|
      >(7/10)*mu*|Ai'(xi0)|>0.                           (E18)

Finally zeta_z and x_z=m^d*a^d*delta are both nonzero, and the
normalization in Section1 is nonzero. Multiplicity therefore transfers
from P_m back to A_m: the actual polynomial root is simple. Real
coefficients supply the conjugate endpoint statement.

## 5. Exact scope of the conclusion

This closes the common holomorphic-collar and endpoint derivative
obligation using the actual (U12) vector constants on their full stated
physical reserves. It treats the entire old endpoint window and
therefore does not discard an annular part of that window by shrinking
a fold disk. The threshold here is the already available m>=10^16;
no numerical claim about all roots or about smaller orders follows.

A proof of global within-order squarefreeness still has to handle every
other root region with its own derivative/count argument. Neither
support-uniform squarefreeness nor all-order SFC is inferred in this note.

