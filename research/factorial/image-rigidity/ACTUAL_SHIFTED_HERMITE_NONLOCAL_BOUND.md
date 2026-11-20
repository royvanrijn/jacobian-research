# An explicit nonlocal bound for the actual shifted-Hermite period

Let 0<delta<=10^(-8), a=4/delta^2 be a positive integer, m>=1 an integer,
and z=x+2iD, with |x|<=10^(-5) and
18999/20000<=D<=21001/20000. Use the actual phase

 F(v,z)=4/delta^2*[log(1+delta*v/2)-delta*v/2]
           +log(z+v+delta*v^2/4),

and its matched mean Phi and unfolding zeta from the retained matched-data
proof. The analytic phase is used on the bounded deformed contour; the
original integrand is the entire function

 E_m(v)=exp(-2m*v/delta)*(1+delta*v/2)^(am)
                            *(z+v+delta*v^2/4)^m.

Put s=0 if D<=1, s=sqrt(D^2-1) otherwise, H=D-s,
and X=2/5 if D<=1, X=sqrt(4/25-s^2) otherwise.
Then 0<H<=1 and X^2>57/1000.

Define I_m as the integral of E_m on [-2/delta,infinity), and L_m as its
left-to-right integral on v=t-z/2+i*s, -X<=t<=X. Set N_m=I_m-L_m.
The nonlocal normalized contribution is

 R_nonlocal=(m^(1/3)/(2*pi*i))*exp(-m*Phi)*N_m.

**Theorem ActualShiftedHermiteNonlocalBound.** With
M(xi)=|Ai(xi)|+|Ai-prime(xi)| and xi=m^(2/3)*zeta,

 |R_nonlocal| < 60000*sqrt(m)*exp(-33m/50000)*M(xi).       (N1)

In particular it is below
10^15*m^2*exp(-m/100000)*M(xi) for m>=10^16.
This is a bound for the specified original-period nonlocal pieces only.
It does not identify the local canonical chain, bound its complete
one-integration-by-parts remainder, or prove a numerical Strong Factorial
threshold.

## 1. The exact contour and limiting horizontal losses

Deform only the finite original interval [-3,3] to its horizontal translate
at height -H, using verticals at real parts -3 and3. The horizontal parameter
is t=Re(v)+x/2 and therefore runs from -3+x/2 to3+x/2.
This exact endpoint convention avoids silently replacing those limits by
-3,3. In particular |t|<=3.000005. The selected local subinterval lies
strictly between these endpoints.

All logarithmic denominators stay away from zero in this finite
rectangle: |v|<16/5, Re(1+delta*v/2)>0 and
Im(z+v+delta*v^2/4)>.8998. Alternatively the displayed entire E_m
makes Cauchy deformation immediate. The remaining nonlocal chain consists
of the horizontal pieces outside [-X,X], the two verticals, and the two
original real tails. Its total horizontal length is at most6 and its
total vertical length at most2, hence below3.

Write W0=F0-Phi0, Phi0=-z^2/4-1/2+i*pi/2. On the limiting x=delta=0
horizontal, if D<=1,

 Re W0=-1/2*[q-log(1+q)], q=t^2+D^2-1.

Outside the cut, q>39/625. For q>0,

 q-log(1+q)=integral_0^q r/(1+r)dr >=q^2/[2(1+q)].

Hence the loss is greater than
(39/625)^2/[4(1+39/625)]>9/10000.

If D>=1, put V=D-s, T=D+s=1/V. The stable saddle is v=-iV.
Its real normalized action is -kappa(D), where

 kappa(D)=D^2-1/2-V^2/2-log(T)>=0.

On its horizontal,

 Re W0+kappa=-t^2/2+(1/2)*log(1+t^2/T^2)
             <=-t^4/[4(1+t^2)].

Outside the cut, t^2>=X^2>57/1000, so this loss is greater than
(57/1000)^2/[4(1+57/1000)]>7/10000.
Set kappa=0 for D<=1. For D>1 write D=cosh(eta).
Then kappa=sinh(2eta)/2-eta=(2/3)*zeta0^(3/2), using the selected
positive Hermite unfolding. The matched-data bound zeta0<1/8 gives
0<=kappa<3/100.

Thus the entire limiting nonlocal horizontal satisfies

 Re W0 <=-kappa-7/10000.                                 (N2)

## 2. Actual horizontal and critical-action perturbations

Throughout the horizontal, |z+v|>=D+s>=18999/20000 and |v|<16/5.
At fixed physical u=t+i*s,

 partial_z W0=(u+z/2)/2+1/[2(u+z/2)].

Here |u+z/2|<33/10: use D+s<1.372 and |t+x/2|<=3.00001.
Thus |partial_z W0|<11/5. The delta derivative of the first log
term is bounded by

 |v|^3/[6*(1-delta*|v|/2)^2],

and that of the second by

 |v|^2/[4*(|z+v|-delta*|v|^2/4)].

Their sum is below9 on the horizontal. The selected matched mean has
|partial_delta Phi|<2. Consequently

 |W_delta(t-z/2+i*s,z)-W0(t-iH,2iD)|
       <(11/5)*10^(-5)+11*10^(-8)<3/100000.               (N3)

For the relative Airy norm we also need an actual critical-action bound.
The two selected critical roots lie in |v+z/2|<2/5. Indeed their polynomial
is P_delta=v^2+z*v-1+delta*(v^3/4-v)-delta^2*v^2/4.
On that circle |1+z^2/4|<.103, while |v+z/2|^2=.16 and
the delta perturbation is below3*10^(-8).
Rouche gives precisely the two selected roots there. Hence |v_j|<3/2
and |q|<2/5 for q=(v1-v2)/2, c=(v1+v2)/2.

Let h(v)=v+delta*v^2/4 and p(v)=1+delta*v/2.
At a critical point, z+h(v_j)=p(v_j)^2/v_j.
If w_j=F(v_j)-Phi, then

 partial_z w1 =-q*(1+delta*c/2)*v1*v2/[p(v1)^2*p(v2)^2],
 partial_z w2 =-partial_z w1.

It follows that |partial_z w_j|<901/1000.
The stationary identity also gives |partial_delta F(v_j)|<1.41;
together with |partial_delta Phi|<2 it gives
|partial_delta w_j|<7/2.
The individual critical values may be Puiseux functions at a merger;
no holomorphy of each labelled action in the base parameter is assumed.
Along either straight parameter segment they can be continued
piecewise across isolated mergers with the displayed bounded
derivatives. Their action derivatives have finite limits there.
The merger zero set on a segment is discrete unless it is identically
zero; in the latter case both normalized critical values are zero.
Integrating on the remaining pieces and passing to the endpoints
therefore gives a matching with the two limiting critical actions such that

 |w_j(z,delta)-w_j(2iD,0)|
       <(901/1000)*10^(-5)+(7/2)*10^(-8)<10^(-5).         (N4)

The limiting real critical actions are 0,0 if D<=1 and
-kappa,+kappa if D>1. Since the actual actions are
+/-2*zeta^(3/2)/3, (N4) gives

 |Re(2*zeta^(3/2)/3)| <=kappa+10^(-5).                    (N5)

The sign of the fractional power is immaterial.

## 3. A uniform lower bound for the joint Airy norm

The retained explicit Airy envelope proof, Sections2--4, gives the
following lower estimates, not only an upper contour envelope.
For |xi|<=4, M(xi)>1/36000.
For |xi|>=4, its finite n=1 errors and rotated connection argument give

 M(xi)>=(3/16)*|xi|^(-1/4)
                          *exp(-|Re(2*xi^(3/2)/3)|).

The latter weaker lower bound holds in both sectors: on the negative
axis sector the proof actually gives the sum of both exponentials.
Using |zeta|<1/8 and (N5), both cases imply the single bound

 M(m^(2/3)*zeta) >
  m^(-1/6)/36000 *exp[-m*(kappa+10^(-5))].                (N6)

This compares to the actual critical action, so no unjustified condition
m*delta small is used.

## 4. Original tails and far verticals

On v=+/-3+i*r, -H<=r<=0, |v|<16/5 and
|z+v|^2<(3+10^(-5))^2+(2.1001)^2<14.
At delta0, Re F0<=-4+(1/2)*log14<-5/2, since exp3>14.
The actual delta derivative is below9 throughout these verticals:
the same derivative formula applies with |z+v|>=1.8999-1>.8998.
Thus

 Re F_delta<-19/10                                      (N7)

there. At the two real endpoints v=+/-3, the limiting bound is
Re F0<-3, so the actual endpoint phases are below-2.

For real v>=3, let p=1+delta*v/2 and h=v+delta*v^2/4.
Here h-|x|-2p=v-2-|x|+delta*v*(v-4)/4 is increasing for v>=3 and at3 equals1-|x|-3*delta/4>0. Thus h-|x|>2p, and the second-log real derivative is at most1/2.
The first derivative is -v/p<=-3/(1+3delta/2), so the full
real derivative is below-2. For -2/delta<v<=-3, 0<p<=1,
-v/p>=3 and the second-log real derivative has absolute value
at most1/Im(z)<1. The full derivative is above2.
The two original tail integrals are therefore bounded in absolute
value by exp(-2m)/(2m) each. The endpoint at -2/delta has E_m=0;
the infinite endpoint is exponentially dominated, so no omitted
endpoint contribution is present.

Finally Re Phi>D_min^2-x^2/4-1/2-2delta>2/5. Thus all
vertical and original-tail contributions together satisfy

 |N_far| <4*exp(-19m/10),
 exp(-m*Re Phi)|N_far|<4*exp(-23m/10).                    (N8)

## 5. Combined original-period bound

By (N2)--(N3) the nonlocal horizontal contributes at most

 6*exp[m*(Re Phi-kappa-67/100000)].

Combining with (N6), its ratio to exp(m*Re Phi)*M is below
216000*m^(1/6)*exp(-33m/50000).
By kappa<.03 and (N8), the far contribution is below
144000*m^(1/6)*exp(-2.26999m) in the same ratio.
Their sum is below360000*m^(1/6)*exp(-33m/50000).
Multiplication by m^(1/3)/(2*pi), using pi>3, proves (N1).

The indicated weaker target follows from
60000<10^15, sqrt(m)<=m^2 and33/50000>1/100000 for m>=1.
No local-canonical-chain or complete-period remainder conclusion
is inferred from this nonlocal estimate.


## Exact reproduction and source scope

    python3 verify_actual_shifted_hermite_nonlocal.py \
      --repo-root /path/to/checkout --output /path/to/review.json

The verifier uses exact rational margins and checks all eight mandatory
source pins. It imports the selected matched data and the finite joint-Airy
lower estimates from their stated retained proofs. The full-window cubic
inverse is not a premise of this nonlocal theorem: the local canonical
chain and its complete remainder are separate. The analytic integral and
contour inequalities remain written proof obligations discharged above,
not consequences of numerical sampling.
