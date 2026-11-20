# Certified cubic images of the actual horizontal endpoints

The computation below proves only the locations of the selected cubic images of two actual physical endpoints. It does not prove a complete-period contour connection, a period remainder constant, or a numerical Strong Factorial threshold.

Write z=x+2iD with

|x|<=10^(-5), 18999/20000<=D<=21001/20000,
|delta|<=10^(-8).                                         (H1)

The physical coordinate is u=v+z/2. Let W_delta=F_delta-Phi and zeta be the matched data of QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md, and let T(u,z,delta) be the regular cubic coordinate of FullActualEndpointInverse, so

T^3/3-zeta*T=W_delta.

Define the Hermite horizontal endpoints at the same real D by

u_+=2/5, u_-=-2/5,                      D<=1;
s=sqrt(D^2-1), X=sqrt(4/25-s^2),
u_+=X+i*s, u_-=-X+i*s,                 D>=1.              (H2)

They have |u_+|=|u_-|=2/5. For nonzero x or delta retain exactly these u endpoints, setting v=u-z/2.

**Theorem HorizontalEndpointCubicLocation.** Throughout(H1), the selected coordinates T_+=T(u_+,z,delta), T_-=T(u_-,z,delta) satisfy

379/1000<|T_+|,|T_-|<403/1000,
Re T_+,Re T_->1/100,
Im T_+<-1/5, Im T_->1/5.                                (H3)

Writing T=X_t+iY_t, both coordinates also satisfy

Y_t^2-X_t^2/3>1/1000,
3*(Y_t^2-X_t^2/3)+Re zeta>1/10.                          (H3a)

Consequently both endpoints lie inside |t|<17/40 and outside |t|<=3/8, with their imaginary parts on opposite sides. Their stationary branch labels are part of the certificate, rather than inferred from a numerical choice of a cubic root.

## 1. A real analytic formula for zeta through D=1

Put e=D-1 and

A(e)=sum_(k>=0) binom(1/2,k)*3/(2k+3)*(e/2)^k.

The selected real Hermite value is

zeta_0(D)=2e*A(e)^(2/3),                                (H4)

where A(e)>0 near the entire interval in(H1). To verify this, for D>1 the matched critical action gives

zeta_0(D)^(3/2)=(3/2)*(D*sqrt(D^2-1)-arcosh D)
 =3*integral_1^D sqrt(t^2-1)dt
 =2*sqrt(2)*e^(3/2)*A(e).

The signed form(H4) analytically continues across e=0 to D<1, and selects the branch zeta_0=2e+O(e^2). It agrees there with

zeta_0(D)=-[(3/4)*(2*theta-sin(2*theta))]^(2/3),
theta=arccos D.

The verifier evaluates the first33 terms k=0..32 using exact rational coefficients and outward-rounded Arb arithmetic. On its slightly enlarged domain, |e|/2<51/2000. Since |binom(1/2,k)|<=1 and3/(2k+3)<=1, the absolute omitted tail is bounded by

(51/2000)^33/(1-51/2000).

It checks99/100<A<101/100 on every parameter cell, so the positive real cube root and its square are unambiguous. In particular no acosh/arccos interval singularity or fractional-power sign ambiguity occurs at D=1.

At a Hermite endpoint the exact normalized phase is

W_0=-(u-iD)^2/2+log(u+iD)-D^2+1/2-i*pi/2,                (H5)

with the logarithm in the upper half-plane. For the endpoint u_-=-conjugate(u_+), equation(H5) gives W_-=conjugate(W_+) and zeta_0 is real.

## 2. Exact perturbations covering the actual family

Along the x segment at fixed u, |u|=2/5 and |D|<1051/1000. The exact moving-centre derivative is

W_(0,z)=(u+z/2)/2+1/[2*(u+z/2)].

Its modulus is below2, using |u+z/2|>=18999/20000-2/5. At these endpoints |v|<3/2. Direct differentiation of the actual phase gives

|F_delta derivative in delta|
 <=(3/2)^3/[6*(1-|delta|*(3/2)/2)^2]
   +(3/2)^2/[4*(18999/20000-2/5-|delta|)]
 <2.

The selected matched mean has |Phi_delta|<2, by the retained data theorem. Thus

|W_actual-W_Hermite|<2*10^(-5)+4*10^(-8)
 =501/25000000=:E_W.                                    (H6)

The same theorem supplies |zeta_z|<5/4 and the delta transfer bound(49/100)/99999. Therefore

|zeta_actual-zeta_Hermite|
 <(5/4)*10^(-5)+(49/100)/99999
 =139199/7999920000=:E_zeta.                             (H7)

Both bounds apply to both endpoints, even though their actual complex perturbations need not remain conjugate.

## 3. The certified closed cover and cubic-root test

The verifier covers the D<=1 interval directly. For D>=1 it uses the nonsingular parameter s in(H2), on the slightly enlarged rational interval0<=s<=801/2500. Exact guards show that this covers D<=21001/20000, keeps X positive, and keeps D<2101/2000. The resulting matched-data domain remains |z-2i|<101/1000. Only the subset(H1) is used in the theorem.

Each cell has rational endpoints, an exact rational complex centre t_0, and radius r=1/1000. Floating-point Newton iteration proposes t_0 only; it has no certifying role. At the proposed centre define

P(t_0)=t_0^3/3-zeta_Hermite*t_0-W_Hermite.

Using outward-ball values over the ENTIRE closed parameter cell, the verifier accepts only if

|P(t_0)|+E_W+|t_0|*E_zeta+|t_0|*r^2+r^3/3
 < (|t_0^2-zeta_Hermite|-E_zeta)*r.                       (H8)

On |t-t_0|=r the linear term dominates all other terms. Rouche therefore gives exactly one actual cubic root in that disc, for every parameter in the cell and every actual perturbation in(H6)--(H7). This root is simple. The program additionally verifies

|t_0|-r>379/1000, |t_0|+r<403/1000,
Re t_0-r>1/100, Im t_0+r<-1/5,

for the plus endpoint. It additionally verifies, with L=-Im(t_0)-r and U=Re(t_0)+r,

L^2-U^2/3>1/1000,
3*(L^2-U^2/3)+zeta_Hermite-E_zeta>1/10.

These inequalities prove(H3a) over the whole actual root disc. Conjugate discs give the minus bounds using its own actual perturbation estimate.

The completed run uses7262 closed cells,14522 attempted cells and maximum bisection depth13. The parameter cells on each side are sorted and checked to join exactly without gaps. The full stream SHA256 is

3daab939d21f7735b25c06f48b19bda42beb79754951d4cde49f0af1a9f5bf00.

## 4. Certifying the branch, not merely an arbitrary root

At each shared cell endpoint the verifier independently certifies the Hermite cubic root in a point disc of radius10^(-9), again by(H8), now with E_W=E_zeta=0. It proves that this smaller disc lies strictly inside BOTH neighboring certified discs. This supplies7264 shared-point certificates, including the two side endpoints. Hence the cell roots continue as one branch, without inferring branch identity merely because two large discs overlap.

At D=1, u=2/5, put

S(u)=log(1-iu)+iu-u^2/2.

The known critical coordinate is

T_base=-i*u*[3*S(u)/(i*u^3)]^(1/3),

using the cube root near1. Both parameter chains meet a point disc containing this exact base branch. The critical inverse normalization therefore fixes the branch along the full Hermite path. For every fixed D, continuous interpolation of x and delta from0 stays in the unique root disc(H8). The full regular-sheet theorem makes the selected physical coordinate continuous and holomorphic along that interpolation. It therefore remains the same root, proving(H3) for the actual family.

The negative controls reject a proposed single coarse parameter cell, a false global upper radius.38, and the wrong conjugate branch in the actual plus-endpoint cubic. The source and receipt prove locations and labels only. They do not establish that a prescribed original contour connects to any canonical Airy contour inside a quantitative decay sublevel set.

## 5. A strictly decreasing radial connector to a common cubic radius

**Corollary EndpointRadialCubicDecay.** For either selected endpoint let t(lambda)=lambda*T, with1<=lambda<= (41/100)/|T|. This path stays within the proved inverse disk, from radius<.403 to the common radius.41. Along it

d/dlambda Re(t(lambda)^3/3-zeta*t(lambda))<-99/100000.       (H9)

Indeed writing T=X+iY with X>1/100 and a_t=Y^2-X^2/3>1/1000, the exact derivative is

-X*(3*lambda^2*a_t+Re zeta)+(Im zeta)*Y.

By(H3a) its first term is below-1/1000. Hermite zeta is real, so |Im zeta|<E_zeta; and |Y|<403/1000. The exact rational inequality

1/1000-(139199/7999920000)*(403/1000)>99/100000

proves(H9). In particular the radial connector cannot raise the normalized phase. This is only a local canonical connector; the rest of the original-period deformation and its full remainder bound are not inferred from this statement.
