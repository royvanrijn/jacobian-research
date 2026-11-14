# Complete exclusion of the quartic pole pair three and four

14 November 2026. PF2RULED4 excludes the full reduced pole pair (3,4)
under PF2RULED1's compatible-completion hypothesis. Together with
PF2RULED2, no compatible quartic finite map has any target-coordinate
presentation of pole height at most four.

The argument covers all leading cubic data, their positions relative
to the marked point of the conic, all constant lifts, and both roots
of the remaining quadratic equation. The full finite-map problem at
higher pole height and the actual existence of a compatible completion
remain open.

## A normalized cubic relation

Use the [four-section normal form](POLE_FOUR_FRONTIER.md). Replacing its
functions by A=3X, H=3Y and
C=2*kappa*X+ell*T+(4/3)*Z gives

\[
 \begin{aligned}
 P&=T^3+3TX+3Y,\\
 Q&=T^4+4T^2X+4TY+2X^2+\kappa T^2+2\kappa X+\ell T+\tfrac43Z.
 \end{aligned}\tag{1}
\]

Here T,X,Y,Z lie in V1=H0(S,O(B)). The scalar ell absorbs the integration
constant of the leading differential equation. Their symbols
t,u,v,r are cubic sections on the conic B and satisfy

\[
 \boxed{u^3-3tuv+t^2r=0.}\tag{2}
\]

The required differential has divisor E-3B. Thus every pole above three
must cancel, and its pole-three coefficient must cut out a fourth-power
contact quartic on B.

## The proportional branch is globally impossible

If u=mu*t, (2) gives r=3mu*v-mu^3*t. Since V1/C maps isomorphically
to the cubic sections, the functions satisfy

\[
 X=\mu T+\chi,\qquad Z=3\mu Y-\mu^3T+\zeta
\]

for constants chi,zeta. The map (1) factors polynomially through (T,Y).
Its outer Jacobian is

\[
 (12\mu^2-12\chi-6\kappa)T-12Y
 -6\kappa\mu-3\ell+4\mu^3.
\]

The coefficient of Y is -12, so this is a nonconstant polynomial.
The [prime-ramification factorization theorem](PRIME_RAMIFICATION_FACTORIZATION.md)
rules it out. The two intermediate functions are algebraically independent
because P,Q are independent. No finiteness of the intermediate map is
assumed in this divisor argument.

If t has three simple zeros, (2) forces u to vanish at all of them,
hence u=mu*t. This excludes the whole squarefree-t branch as well.

## All nonproportional cubic data have two root types

Suppose u is not proportional to t. Over C there are two remaining
possibilities, up to a nonzero scalar in t.

If t=L^2 M for distinct linear forms L,M, the orders in (2) force
u=L M N. Nonproportionality means N is not proportional to L. Solving
for r gives

\[
 \boxed{
 t=L^2M,\quad u=LMN,\quad
 v=\tfrac13MN^2+LD,\quad r=3ND,}
 \tag{3}
\]

where D is an arbitrary binary quadratic. Indeed
r=N(3v-MN^2)/L and gcd(N,L)=1, which proves both necessity and sufficiency
of the displayed polynomial divisibility.

If t=L^3, the order of u at L cannot be one: its cube would have order
three while the other terms of (2) have larger order. Thus u=L^2N,
with N not proportional to L. Solving (2) then gives

\[
 \boxed{
 t=L^3,\quad u=L^2N,\quad v=LD,\quad r=3ND-N^3.}
 \tag{4}
\]

These forms include every nonproportional solution of (2).
The zero section u is in the proportional branch already excluded.

## The marked point is preserved in the atlas

Normalize S as the blowup of p0=[0:0:1] in P2, with B the strict
transform of XZ-Y^2=0. Write [s:t] for the conic parameter, so p0 is
s=0. Projective automorphisms of the conic fixing p0 lift to
automorphisms of S preserving B; they give the full PGL2 stabilizer
of that marked point.

For (3), there are exactly three positions relative to p0. We may use

\[
 (L,M)=(s,t),\quad(t,s),\quad(s-t,t).
\]

They correspond to the double root at p0, the simple root at p0, and
neither root at p0. In each chart write N=n0*L+n1*K, with K independent
of L; nonproportionality is exactly n1!=0.

For (4), the triple root is either p0 or away from p0. We may use
L=s or L=t. These are the two triple-root charts. A weighted target
rescaling (P,Q)->(c^3 P,c^4 Q) normalizes the remaining nonzero scalar
of t. All of these transformations preserve the compatible-completion
class, the B valuation, generic degree and contact-four requirement.

## Mixed double and simple roots fail the next coefficient

Use the four global functions

\[
 g_1=y/b,\quad g_2=y^2/b,\quad
 g_3=xy/b,\quad g_4=x^2/b,\qquad b=x-y^2.
\]

Together with1 they form a basis of V1. Their cubic symbols are
1,w,w^2,w^3 in the usual trivialization. For a cubic
f=f0+f1*w+f2*w^2+f3*w^3, every lift is

\[
 \widetilde f=
 \frac{wf(w)}z+f_2w+2f_3w^2+f_3z+c,\qquad z=b,\ w=y,
 \tag{5}
\]

with an arbitrary constant c.

Substituting the full families (3), all constant lifts and kappa,ell
into (1), the exact pole-four coefficient F4(w) has these consequences:

| Chart | Coefficient functional | Result |
| --- | --- | --- |
| Double root at p0 | coefficient of w^4 | -(20/3)*n1^4 |
| Simple root at p0 | coefficient of w^9 | (20/3)*n1^4 |
| Neither root at p0 | derivative F4'(1) | -(20/3)*n1^4 |

Every result contradicts n1!=0. These are polynomial identities in every
remaining parameter, not sampled specializations. The retained
[full coefficient atlas](pole4_symbol_atlas.json) and
[localized unit certificates](pole4_mixed_certificates.json) verify them.
For example, if a coefficient functional gives c*n1^4, adjoining
h with h*n1=1 yields the explicit identity

\[
 1=\frac{h^4}{c}(c\,n_1^4)
 +(1-hn_1)(1+hn_1+(hn_1)^2+(hn_1)^3).
\]

No mixed-root chart survives.

## Triple-root normalization loses no candidate

In (4), first write N=n0*L+n1*K. The choice of lift T of its leading
symbol is free. Replacing T by T+c changes the symbol u to u-c*t
when the normal form is reconstructed. Taking c=n0 therefore makes
n0=0. One may see the first-component identity directly from

\[
 T'=T+c,\quad X'=X-cT,\quad Y'=Y-cX-c^3/3,
\qquad T'^3+3T'X'+3Y'=P.
\]

Reconstructing Q's normal form may apply Q->Q-dP; that does not change
the differential. This is a representation choice and an allowed
target change, not an equivalence of different root partitions.

Next normalize n1 to1. The conic-preserving automorphism
(x,y)->(a^2x,ay), a!=0, acts on (g1,g2,g3,g4) by
(a^-1,1,a,a^2). In the chart L=s, a weighted target rescaling with
c=a^-2 keeps t=s^3 normalized and sends n1 to n1*a^-3.
In the chart L=t, use c=a; it sends n1 to n1*a^2.
Cube or square roots exist over C, so either nonzero n1 can be
normalized to1. The zero case was excluded separately.

Thus all triple data are covered by

\[
 \begin{array}{c|c|c}
 \text{chart} & (L,N) &D\\ \hline
 \text{at }p_0 &(s,t)&a s^2+b st+\delta t^2\\
 \text{away from }p_0 &(t,s)&\delta s^2+b st+a t^2 .
 \end{array}\tag{6}
\]

## Solve every pole-four equation in the triple charts

In either chart, the pole-four coefficient forces

\[
 f(\delta)=6\delta^2-12\delta+5=0.
 \tag{7}
\]

The other three equations are triangular in the lift constants t0,x0,ell.
Their relevant coefficients are nonzero multiples of
3delta-2, delta-1 and1. No root is lost: f(1)=-1 and f(2/3)=-1/3.

Let xbar,ellbar denote the solutions with a=y0=z0=0. Modulo f,
the two charts have

\[
 t_0=(9-6\delta)b
\]

and

\[
\begin{array}{c|l|l}
 &\bar x&\bar\ell\\ \hline
 \text{at }p_0&
 (-21+15\delta)b^2+(1-\tfrac32\delta)\kappa&
 (84-60\delta)b^3+(-6+6\delta)b\kappa+4\delta-\tfrac43\\
 \text{away}&
 (-21+15\delta)b^2+(1-\tfrac32\delta)\kappa+2\delta-1&
 (84-60\delta)b^3+(-6+6\delta)b\kappa+(8-8\delta)b .
\end{array}
\tag{8}
\]

The general solutions are x0=xbar-a and
ell=ellbar+4a*t0-4y0. The constants a,y0,z0 only translate P,Q.
Indeed their differences from the normalized pair are

\[
 P_{\rm gen}-P_{\rm norm}=3(y_0-a t_0),\qquad
 Q_{\rm gen}-Q_{\rm norm}=2a^2-4a\bar x-2\kappa a+\tfrac43z_0.
\]

Thus their differentials, and the contact quartic, are unchanged.
This removes three parameters by exact identities rather than an
unproved symmetry. The remaining parameters b,kappa are unrestricted,
and (7) retains both algebraic values of delta simultaneously.

## Contact four is a fourth-power binary quartic

A pole-three two-form on S corresponds to a section of
K_S+3B=3L_plane-2C0. In the plane chart it is

\[
 \Omega=F(x,y)\,dx\wedge dy/(x-y^2)^3,
\]

where F is a cubic with multiplicity at least two at p0. Restriction
to the strict conic is the binary quartic obtained by homogenizing

\[
 G(w)=F(w^2,w)/w^2.
 \tag{9}
\]

The division by w^2 accounts for the double point at the blown-up
marked point; it must not be omitted. The contact condition B.E=4
supported at one point says G is a nonzero fourth power of a linear form,
up to scalar. Homogenization includes the point at infinity.

Write G=sum q_i*w^i and c_i=q_i/binomial(4,i). Every fourth power satisfies

\[
 F_2=c_0c_2-c_1^2=0,\quad
 F_3=c_0^2c_3-c_1^3=0,\quad
 F_4=c_0^3c_4-c_1^4=0.
 \tag{10}
\]

For the away chart, reverse the five coefficients first; fourth-power
contact is invariant under that reversal. These are necessary equations;
their insufficiency for a general quartic cannot weaken an exclusion.

## Both remaining contact systems have exact contradictions

The [retained quartics](pole4_triple_candidates.json) give, modulo f,
the following simple consequences.

In the marked chart,

\[
 F_2=(70-120\delta)\kappa,\qquad
 F_3=(2100-3600\delta)b\kappa+(1900-3200\delta)/3.
\]

Since f(7/12)=1/24, the first equation forces kappa=0. The second
then forces delta=19/32, but f(19/32)=-5/512. This is impossible.
An explicit polynomial unit certificate is

\[
 1=-\tfrac{512}{5}f+
 b\left(\tfrac{432}{25}\delta-\tfrac{243}{10}\right)F_2+
 \left(\tfrac{81}{100}-\tfrac{72}{125}\delta\right)F_3 .
\]

In the away chart,

\[
 F_2=(70-120\delta)(\kappa-5/3).
\]

Thus kappa=5/3. After this substitution the next equations are

\[
 F_3=-\tfrac{100}{3}(18\delta-11)b,\qquad
 F_4=\tfrac{500}{9}\big((432\delta-264)b^2+71-120\delta\big).
\]

Now f(11/18)=-5/54, so F3 forces b=0. Then F4 forces
delta=71/120, but f(71/120)=1/2400. Again this is impossible.
The [certificate file](pole4_contact_certificates.json) gives a full
rational polynomial identity expressing1 in (f,F2,F3,F4).

Both are characteristic-zero exact contradictions. No approximation
to the quadratic roots is used.

## Coverage and consequence

Every leading cubic t is squarefree, double-plus-simple, or triple.
The squarefree case and every proportional u are excluded by the
global factorization obstruction. The three marked mixed-root positions
are excluded by exact coefficient identities. Both normalized triple
positions are excluded by contact, with every lift accounted for above.

Therefore the full pair (3,4) is impossible. Combined with the earlier
pole reduction,

\[
 \boxed{\text{Every compatible quartic finite map has target pole
 height at least five in every polynomial target presentation.}}
\]

The subsequent [all-height coordinate bound](LOW_COORDINATE_POLE.md)
excludes every (2,n) pair. At height five the first reduced pairs are
(6,8),(7,8), after the
[source-coordinate transfer](POLE_SEVEN_EXCLUSION.md) excludes all heights
through seven. The intrinsic [complete (3,5) proof](POLE_THREE_FIVE_EXCLUSION.md)
is retained separately.
Those pairs, higher poles and actual ruled compatibility remain open.
The factorization theorem also applies directly to both actual quartic
packets, independent of that compatibility assumption.

The [checker](verify_pole4_exclusion.py) reconstructs the coefficient
atlas by direct quotient differentiation, checks normalization identities
and all unit certificates, and independently reconstructs the two contact
quartics. Root coverage, geometric normalization, canonical divisors and
the contact interpretation are written proofs, not formal certification.
No independent whole-proof implementation, external review or novelty
priority is claimed.


The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
