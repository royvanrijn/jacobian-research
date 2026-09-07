# F2 `k=1` carrier-jet factorization

> **Status.**  Exact in both the normalized and fixed target charts.  The
> `k=1` implicit quintic and its fixed-coordinate transforms are irreducible,
> so factoring the target polynomial cannot reveal a missing component.  On
> the carrier, however, its normalized seven-jet lies on a prime
> codimension-three complete intersection: the first four jet coefficients
> recover `a,b,c,d`, and the last three are forced.  This gives three finite
> compatibility equations for maximal carrier contact once the target
> normalization parameters are fixed.  The inverse fixed
> target normalization is an explicit weighted triangular carrier
> automorphism, so the equations now apply directly to the raw fixed carrier
> centers and five target-normalization parameters.  Its raw seven-jet
> Jacobian is three times `Res(p',q')`; hence on the immersed locus the free
> normalization parameters absorb the three residuals and carrier jets alone
> give no generic obstruction.  On the `E_6+A_1` subfamily, Lagrange
> inversion instead gives the complete binomial carrier sequence; its
> normalized seven-jet is a one-parameter monomial curve and remains
> codimension three after all three target transports.  Its `E_8` endpoint
> is the hyperplane section `J_1=0`, a prime codimension-four complete
> intersection with four scale-free equations in the raw centers.  Apparent full
> factors `C^3` and
> `s^2` in the zero-translation chart are coordinate artifacts.  The
> terminal leading fiber
> `s^9(135s^3+405s^2+396s+125)` is invariant and reproduces the already
> known terminal passport; it is not a new affine-row factor.

The identities in this note are replayed by
[`verify_f2_affine_k1_carrier_jet_factorization.py`](../scripts/verify_f2_affine_k1_carrier_jet_factorization.py).

## 1. What can and cannot factor

Write the normalized target parametrization as

\[
 p(t)=t^3+at,\qquad
 q(t)=t^5+bt^4+ct^2+dt,                          \tag{1.1}
\]

and let `F(P,Q)` be the twelve-support implicit quintic from
[`F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md`](F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md).
For every specialization of `a,b,c,d` over a characteristic-zero field,

\[
 [k(t):k(p,q)]\mid 3,
 \qquad
 [k(t):k(p,q)]\mid 5.                            \tag{1.2}
\]

Thus `k(p,q)=k(t)`.  The parametrization is birational onto its image, and
its resultant has exponent one.  Consequently `F` is irreducible.  The
fixed-coordinate polynomial `G` of `PF2K1PB1` is obtained from `F` by an
invertible affine target change and multiplication by a unit, so it is also
irreducible.

This does **not** prove that `G(P_s,Q_s)` is irreducible for a Keller map.
An étale inverse image of an irreducible curve may be disconnected.  It does
show that any useful factorization must occur in the source pullback or in a
boundary initial form, not in the target quintic itself.

As a regression experiment, exact factorization over `QQ` found one
irreducible quintic for all

\[
 (a,b,c,d)\in\{-2,-1,0,1,2\}^4                  \tag{1.3}
\]

(625 specializations).  The experiment is not used in the proof; (1.2) is
the general argument.

## 2. The common-power edge restriction

Normalize the fixed F2 common-power edge by replacing `-Q` with
`(5/9)(-Q)`.  Its leading pair is then

\[
 P=C^3,\qquad Q=C^5.                              \tag{2.1}
\]

Put `z=C^{-1}` and

\[
 S=a^4+a^3b^2-2a^2bc+2a^2d+ac^2+d^2.           \tag{2.2}
\]

Direct substitution gives the exact finite weighted restriction

\[
\begin{aligned}
z^{15}F(z^{-3},z^{-5})={}&
3bz-5az^2+(ab+b^3+3c)z^3\\
&+(-5a^2-4ab^2+3bc+3d)z^4
 +2a(ab-c)z^5\\
&+(-abc+4ad+3b^2d+3c^2)z^6\\
&+(a^3b-3a^2c-5abd+3cd)z^7\\
&+(a^3c+a^2b^2c-a^2bd-2abc^2+5acd+3bd^2+c^3)z^9\\
&-aS z^{10}+dS z^{12}.                           \tag{2.3}
\end{aligned}
\]

There are structural gaps at weights `8` and `11`.  In the translated
normal form, (2.3) is divisible by `z^3` after returning to `C`; equivalently
`F(C^3,C^5)` is divisible by `C^3`.  This is not invariant.  It only records
that the target curve and the pure cusp both pass through `(0,0)` in this
normalization.  For example,

\[
 F(C^3,C^5)|_{(a,b,c,d)=(1,0,0,0)}
 =-C^5(5C^8+5C^6+1).                             \tag{2.4}
\]

After restoring target translations, `G(0,0)` is not identically zero, so
there is no universal `C` factor at all.  The robust information in (2.3)
is its weighted coefficient pattern, not the origin factor.

## 3. The normalized carrier graph

At infinity put `u=t^{-1}` and use the normalized carrier coordinates

\[
 \pi=\frac{p^3}{q^2},\qquad h=\frac{p^5}{q^3}.   \tag{3.1}
\]

Then

\[
\pi=u\frac{(1+au^2)^3}{(1+bu+cu^3+du^4)^2},
\qquad
h=\frac{(1+au^2)^5}{(1+bu+cu^3+du^4)^3}.         \tag{3.2}
\]

Because `pi=u+O(u^2)`, there is a unique graph expansion

\[
 h=1+H_1\pi+H_2\pi^2+\cdots.                    \tag{3.3}
\]

Through the seven carrier centers, the coefficients are

\[
\begin{aligned}
H_1={}&-3b,\\
H_2={}&5a,\\
H_3={}&14ab-b^3-3c,\\
H_4={}&-20a^2+46ab^2-3b^4-12bc-3d,\\
H_5={}&-187a^2b+160ab^3+32ac-9b^5-45b^2c-18bd,\\
H_6={}&175a^3-1204a^2b^2+574ab^4+328abc+41ad\\
     &\quad-28b^6-168b^3c-84b^2d-12c^2,\\
H_7={}&2754a^3b-6630a^2b^3-442a^2c+2100ab^5\\
     &\quad+2250ab^2c+500abd-90b^7-630b^4c\\
     &\quad-360b^3d-135bc^2-30cd.               \tag{3.4}
\end{aligned}
\]

The first four equations are triangular.  Their inverse is

\[
\boxed{\begin{aligned}
a={}&H_2/5,\\
b={}&-H_1/3,\\
c={}&(5H_1^3-126H_1H_2-135H_3)/405,\\
d={}&(5H_1^4-90H_1^2H_2-540H_1H_3
      -324H_2^2-405H_4)/1215.
\end{aligned}}                                    \tag{3.5}
\]

Thus four jet coefficients contain the entire target-normal-form packet.
The remaining coefficients are compatibility conditions, not new free
Laurent data.

## 4. Three exact seven-jet equations

Write `x_j=H_j`.  Substitution of (3.5) into the last three rows of (3.4)
gives

\[
\boxed{\begin{aligned}
x_5={}&\frac{1}{2025}\bigl(
10x_1^3x_2-2025x_1^2x_3-2223x_1x_2^2
-4050x_1x_4-4320x_2x_3\bigr),                   \tag{4.1}\\
x_6={}&\frac{1}{54675}\bigl(
-100x_1^6+1575x_1^4x_2+118800x_1^3x_3
+118746x_1^2x_2^2\\
&\qquad+170100x_1^2x_4+63180x_1x_2x_3
-43011x_2^3-149445x_2x_4-72900x_3^2\bigr),       \tag{4.2}\\
x_7={}&\frac{1}{820125}\bigl(
4375x_1^7-78750x_1^5x_2-2868750x_1^4x_3
-2685960x_1^3x_2^2\\
&\qquad-3543750x_1^3x_4+2126250x_1^2x_2x_3
+3736854x_1x_2^3+6561000x_1x_2x_4\\
&\qquad+455625x_1x_3^2+2646270x_2^2x_3
-2733750x_3x_4\bigr).                            \tag{4.3}
\end{aligned}} 
\]

Give `x_j` weight `j`.  Equations (4.1)--(4.3) are weighted homogeneous of
weights `5,6,7`.  Each equation is monic linear in the new variable
`x_5,x_6,x_7`, respectively.  Therefore the normalized seven-jet locus is

\[
 \mathbb A^4\hookrightarrow\mathbb A^7          \tag{4.4}
\]

as a prime codimension-three complete intersection, with coordinate ring
`QQ[x_1,x_2,x_3,x_4]`.

This pattern is not special to order seven.  Formal reversion of (3.2)
gives every `H_j` as a polynomial in `a,b,c,d`.  Since (3.5) is polynomial
over `QQ`, for every `N>=4` the normalized `N`-jet locus is a prime
codimension-`N-4` complete intersection

\[
 x_j-R_j(x_1,x_2,x_3,x_4)=0,
 \qquad 5\le j\le N,                              \tag{4.5}
\]

where `R_j` is weighted homogeneous of weight `j`.  The carrier has an
infinite triangular recurrence with only four genuine moduli.

## 5. Exact fixed-coordinate transport

The compatibility equations above are useful only if the fixed F2 target
translations are retained.  This transport has a compact closed form.
Use the inverse target coordinates from `PF2K1PB1`,

\[
 U=\frac{P-P_0}{A},\qquad
 V=\frac{A(Q-Q_0)-\Gamma(P-P_0)}{AB},             \tag{5.1}
\]

and put

\[
 x=\frac{U^3}{V^2},\quad y=\frac{U^5}{V^3},
 \qquad
 \kappa=\frac{A^3}{B^2},\quad
 \lambda=\frac{A^5}{(-B)^3},                    \tag{5.2}
\]

\[
 X=\frac{\pi}{\kappa},\quad Y=\frac{h}{\lambda},
 \qquad
 \mu=\frac{P_0}{A},\quad
 \eta=\frac{\Gamma}{B},\quad
 \nu=\frac{Q_0}{B}.                             \tag{5.3}
\]

The elementary identities

\[
 U=\frac{y^2}{x^3},\qquad V=\frac{y^3}{x^5}      \tag{5.4}
\]

turn the entire five-parameter affine target change into the exact
three-parameter carrier automorphism

\[
\boxed{\begin{aligned}
X={}&x\frac{(1+\mu x^3/y^2)^3}
 {(1+\eta x^2/y+\nu x^5/y^3)^2},\\
Y={}&y\frac{(1+\mu x^3/y^2)^5}
 {(1+\eta x^2/y+\nu x^5/y^3)^3}.
\end{aligned}}                                    \tag{5.5}
\]

Thus two of the five fixed-target parameters contribute only the leading
scales `kappa,lambda`; the jet transport itself depends only on
`mu,eta,nu`.

Write the normalized and dimensionless fixed graphs as

\[
 y=1+\sum_{j\ge1}H_jx^j,
 \qquad
 Y=1+\sum_{j\ge1}J_jX^j.                         \tag{5.6}
\]

If the raw fixed carrier centers are `zeta_1,...,zeta_7` in

\[
 h-\lambda-\sum_{j=1}^7\zeta_j\pi^j,             \tag{5.7}
\]

then

\[
 \boxed{J_j=\frac{\kappa^j}{\lambda}\zeta_j.}   \tag{5.8}
\]

Expanding the inverse of (5.5) gives the following triangular
**untransport**:

\[
\begin{aligned}
H_1={}&J_1,\\
H_2={}&J_2+3\eta,\\
H_3={}&J_3-2J_1\eta-5\mu,\\
H_4={}&J_4+2J_1^2\eta+8J_1\mu-4J_2\eta-6\eta^2,\\
H_5={}&J_5-2J_1^3\eta-11J_1^2\mu+6J_1J_2\eta
 +15J_1\eta^2+11J_2\mu\\
&\quad-6J_3\eta+30\eta\mu+3\nu,                \tag{5.9}\\
H_6={}&J_6+2J_1^4\eta+14J_1^3\mu-8J_1^2J_2\eta
 -28J_1^2\eta^2-28J_1J_2\mu\\
&\quad+8J_1J_3\eta-104J_1\eta\mu-8J_1\nu
 +4J_2^2\eta+28J_2\eta^2\\
&\quad+14J_3\mu-8J_4\eta+28\eta^3-35\mu^2,     \tag{5.10}\\
H_7={}&J_7-2J_1^5\eta-17J_1^4\mu+10J_1^3J_2\eta
 +45J_1^3\eta^2+51J_1^2J_2\mu\\
&\quad-10J_1^2J_3\eta+240J_1^2\eta\mu+15J_1^2\nu
 -10J_1J_2^2\eta-90J_1J_2\eta^2\\
&\quad-34J_1J_3\mu+10J_1J_4\eta-120J_1\eta^3
 +153J_1\mu^2-17J_2^2\mu\\
&\quad+10J_2J_3\eta-160J_2\eta\mu-10J_2\nu
 +45J_3\eta^2+17J_4\mu\\
&\quad-10J_5\eta-225\eta^2\mu-30\eta\nu.       \tag{5.11}
\end{aligned}
\]

Equations (5.8)--(5.11) are the fixed-coordinate compiler.  Apply (3.5) to
the first four untransported coefficients, then form

\[
 \boxed{\mathcal R_j
 =H_j-R_j(H_1,H_2,H_3,H_4),\qquad j=5,6,7,}       \tag{5.12}
\]

where `R_5,R_6,R_7` are the right sides of (4.1)--(4.3).  The fixed carrier
jet belongs to the `k=1` target locus exactly when all three `mathcal R_j`
vanish.

This remains a weighted triangular complete intersection.  Assign

\[
 \operatorname{wt}(J_j)=j,\qquad
 \operatorname{wt}(\eta,\mu,\nu)=(2,3,5).       \tag{5.13}
\]

Then `H_j` and `mathcal R_j` have weights `j`.  Each residual is linear in
the new variable `J_j`.  Consequently restoring the fixed target does not
increase the intrinsic number of carrier equations or destroy their
triangular order.

### 5.1 The raw fixed seven-jet is dominant

There is a decisive qualification.  The three residuals in (5.12) are
obstructions only after `mu,eta,nu` have been fixed or related by independent
global equations.  Let

\[
 \Phi:(a,b,c,d,\mu,\eta,\nu)
 \longmapsto(J_1,J_2,\ldots,J_7)                 \tag{5.14}
\]

be the fixed seven-jet map obtained from (3.4) and (5.5).  Exact
differentiation and factorization give

\[
\boxed{
 \det(d\Phi)=3\operatorname{Res}_t(p'(t),q'(t)).} \tag{5.15}
\]

For the normal form (1.1), the resultant is

\[
\begin{aligned}
\Delta_{\mathrm{imm}}={}&25a^4+48a^3b^2-144a^2bc
 +90a^2d+108ac^2+81d^2.                         \tag{5.16}
\end{aligned}
\]

Thus on the immersed locus `Delta_imm!=0`, the map (5.14) is étale and in
particular dominant.  The union of the normalized codimension-three jet
loci over the three free target-normalization parameters fills a dense open
subset of the raw fixed seven-jet space.

This changes the obstruction interpretation:

- with `mu,eta,nu` fixed, (5.12) gives three genuine carrier equations;
- with them free, the same equations generically determine normalization
  parameters and do **not** exclude a `k=1` target curve; and
- rank can drop only on `Res(p',q')=0`, precisely the nonimmersion locus
  where the generic conormal/Chern packet already requires separate point
  corrections.

Therefore the fixed carrier seven-jet by itself cannot supply the missing
generic obstruction.  It must be coupled to global equations fixing the
target normalization, to the four affine singular fibers, or to the
boundary normalization/`Fitt_1` module.

## 6. Contact interpretation

Suppose the seven carrier centers, after transport through the **exact
inverse fixed-target normalization** of `PF2K1PB1`, are
`c_1,...,c_7`.  For

\[
 w=h-1-\sum_{j=1}^7c_j\pi^j,                    \tag{6.1}
\]

one has

\[
 \operatorname{ord}_u(w)\ge r
 \quad\Longleftrightarrow\quad
 H_j=c_j\ \text{for }1\le j<r,                  \tag{6.2}
\]

through the order-eight truncation.  Hence:

- matching `c_1,...,c_4` uniquely reconstructs `a,b,c,d` by (3.5);
- the next match is exactly (4.1);
- the next two matches are exactly (4.2) and (4.3); and
- maximal truncated contact `8` requires all three compatibility equations.

These are not three equations on the seven raw fixed-coordinate carrier
coefficients alone.  Equations (5.8)--(5.12) supply the required scalings,
translations, and triangular shear explicitly.  Ignoring that transport is
exactly what produces the false `C^3` factor in Section 2.

For the `(75,125)` compiler this replaces a blind seven-step coefficient
descent by a four-parameter reconstruction followed by three explicit
tests.  It does not supply the missing fixed-coordinate source Laurent pair.

In raw fixed coordinates, use (5.8)--(5.11) first.  Since (5.5) has
`X=x+O(x^2)` and is formally invertible, it preserves the first mismatch
order and hence the contact number.

### 6.1 Exact carrier pattern on the `E_6+A_1` escape

The complement-monodromy stratification isolates the first one-component
topological escape as

\[
 p=t^3,\qquad q=t^5+\beta t^4,qquad\beta\ne0.   \tag{6.3}
\]

On this stratum the carrier functions collapse to

\[
 \pi=\frac{u}{(1+\beta u)^2},\qquad
 h=\frac1{(1+\beta u)^3}.                       \tag{6.4}
\]

Put `z=beta*u` and `X=beta*pi`.  Then

\[
 X=\frac{z}{(1+z)^2},\qquad h=(1+z)^{-3}.        \tag{6.5}
\]

Lagrange inversion gives the complete coefficient pattern

\[
 \boxed{
 H_j=-\frac3j\binom{2j-4}{j-1}\beta^j\quad(j\ge1).}          \tag{6.6}
\]

In particular,

\[
 (H_1,\ldots,H_7)=
 (-3\beta,0,-\beta^3,-3\beta^4,-9\beta^5,
  -28\beta^6,-90\beta^7).                    \tag{6.7}
\]

Thus the normalized seven-jet is the explicit monomial curve cut out by

\[
\boxed{
 H_2=0,\quad27H_3=H_1^3,\quad27H_4=-H_1^4,
 \quad27H_5=H_1^5,
 \quad729H_6=-28H_1^6,
 \quad243H_7=10H_1^7.}                         \tag{6.8}
\]

After adjoining the three carrier transports `(mu,eta,nu)`, elimination is
still triangular.  Four fixed jets recover all four parameters:

\[
\boxed{\begin{aligned}
\beta={}&-J_1/3,\qquad \eta=-J_2/3,\\
\mu={}&-(J_1^3-18J_1J_2-27J_3)/135,\\
\nu={}&-(6J_1^5-109J_1^3J_2-297J_1^2J_3-27J_1J_2^2
       +297J_2J_3+135J_5)/405.
\end{aligned}}                                                   \tag{6.9}
\]

The remaining jets satisfy exactly

\[
\boxed{
 E_4=J_1^4-18J_1^2J_2-72J_1J_3-30J_2^2-45J_4=0,}              \tag{6.10}
\]

\[
\boxed{\begin{aligned}
E_6={}&-187J_1^6+3186J_1^4J_2+11178J_1^3J_3
 +6480J_1^2J_2^2\\
&-972J_1J_2J_3-9720J_1J_5+3780J_2^3
 -5103J_3^2-3645J_6=0,
\end{aligned}}                                                   \tag{6.11}
\]

\[
\boxed{\begin{aligned}
E_7={}&89J_1^7-1604J_1^5J_2-5181J_1^4J_3
 -1314J_1^3J_2^2\\
&+6408J_1^2J_2J_3+3375J_1^2J_5+4131J_1J_3^2
 -1125J_2^2J_3\\
&-2250J_2J_5-675J_7=0.
\end{aligned}}                                                   \tag{6.12}
\]

These weighted equations have degrees `4,6,7` and are monic linear up to a
nonzero scalar in `J_4,J_6,J_7`.  Hence they define a prime
codimension-three complete intersection with free coordinates
`J_1,J_2,J_3,J_5`; the exact `7 by 4` parametrization Jacobian has rank four
at `(beta,mu,eta,nu)=(1,0,0,0)`.  Unlike the immersed chart, free target
transport does **not** make the `E_6` locus dominant.

The leading scale `kappa` can also be eliminated completely.  Since
`J_j=kappa^j*zeta_j/lambda` and `E_w` has carrier weight `w`, define

\[
 \boxed{
 \widehat E_w(\zeta;\lambda)
 =\lambda^w E_w
 \left(\frac{\zeta_1}{\lambda},\ldots,
       \frac{\zeta_7}{\lambda}\right),
 \qquad w=4,6,7.}                               \tag{6.13}
\]

Then `E_w(J)=kappa^w*widehat E_w/lambda^w`, so the exact raw-center gate is

\[
 \boxed{\widehat E_4=\widehat E_6=\widehat E_7=0.}             \tag{6.14}
\]

For example,

\[
\widehat E_4=
 \zeta_1^4-18\lambda\zeta_1^2\zeta_2
 -72\lambda^2\zeta_1\zeta_3-30\lambda^2\zeta_2^2
 -45\lambda^3\zeta_4.                         \tag{6.15}
\]

The checker expands and verifies all three equations.  On the special F2
carrier point one substitutes the already fixed residue
`lambda=125/729`; no unknown leading scale remains.

This is the promised finite carrier gate for the first cusp escape.  To turn
it into an exclusion, the actual scaled fixed F2 center vector
`J_1,...,J_7` must be supplied by the missing global source Laurent pair;
with those data, membership is the three direct tests
`E_4=E_6=E_7=0`, not a 927-variable descent.

### 6.2 The `E_8` endpoint

At `beta=0`, the separated node in the `E_6+A_1` packet coalesces with the
cusp and

\[
 p=t^3,\qquad q=t^5.                              \tag{6.16}
\]

This is the affine `E_8` endpoint.  Its normalized carrier graph is simply
`h=1`, so `H_j=0` for every `j`.  After the three target transports, the
parameters are recovered by

\[
 \boxed{
 \eta=-J_2/3,\qquad \mu=J_3/5,\qquad
 \nu=-\frac{11J_2J_3+5J_5}{15}.}                 \tag{6.17}
\]

The remaining four fixed centers obey

\[
\boxed{\begin{aligned}
 J_1={}&0,\\
 2J_2^2+3J_4={}&0,\\
 140J_2^3-189J_3^2-135J_6={}&0,\\
 5J_2^2J_3+10J_2J_5+3J_7={}&0.
\end{aligned}}                                                   \tag{6.18}
\]

These equations are just `J_1=0` and the specialization of
(6.10)--(6.12), with content removed.  They are linear with nonzero
coefficient in `J_1,J_4,J_6,J_7`, respectively.  Consequently the `E_8`
fixed-jet locus is a prime codimension-four complete intersection with free
coordinates `J_2,J_3,J_5`.  The exact parametrization Jacobian has rank
three at the zero transport point.

The leading scale again cancels.  In the raw fixed centers the gate is

\[
\boxed{\begin{aligned}
 \zeta_1={}&0,\\
 2\zeta_2^2+3\lambda\zeta_4={}&0,\\
 140\zeta_2^3-189\lambda\zeta_3^2
                 -135\lambda^2\zeta_6={}&0,\\
 5\zeta_2^2\zeta_3+10\lambda\zeta_2\zeta_5
                 +3\lambda^2\zeta_7={}&0.
\end{aligned}}                                                   \tag{6.19}
\]

For the special F2 carrier one substitutes `lambda=125/729`.  Thus both
nonimmersed one-component topological escape strata now have finite raw
**maximal-contact** tests: three equations for `E_6+A_1`, and four for
`E_8`.  A packet leaving the carrier earlier only satisfies the initial
center matches preceding its contact order.  What is still missing is the
actual center vector, to the order reached, produced by the global
fixed-coordinate Laurent pair.

## 7. Terminal factor audit

Use a terminal transverse parameter `rho` and residue coordinate `s`.  The
certified leading terminal pair is

\[
 P=\rho^{-3}s^2(1+s),\qquad
 -Q=\rho^{-5}s^3(1+3s+\tfrac95s^2).              \tag{7.1}
\]

In the zero-translation normalized target chart, clearing `rho^{-15}`
produces a universal factor `s^2`.  It is again an origin artifact.  In the
fixed target chart with `A=1,B=-9/5`, one instead has

\[
 \left.\rho^{15}G(P,Q)\right|_{s=0}
 =\rho^{15}G(0,0),                               \tag{7.2}
\]

and `G(0,0)` is not identically zero.  Thus there is no universal `s`
factor in the full fixed-coordinate pullback.

On the special-residue locus `A^5/(-B)^3=125/729`, the leading terminal
fiber is invariant up to its nonzero leading-scale factor:

\[
\boxed{
 \left.\rho^{15}G(P,Q)\right|_{\rho=0}
 =-\frac{B^3s^9}{729}
 (135s^3+405s^2+396s+125).}                      \tag{7.3}
\]

After the parameter normalization `A=1,B=-9/5`, the scalar in (7.3) is
`1/125`.

The cubic is the finite simple part of the known terminal `h=125/729`
fiber: it is irreducible over `QQ` and has discriminant `-98415`.  Together
with the order-three point at infinity it recovers the
existing terminal `A_6` passport.  Equation (7.3) is a valuable regression
for the fixed-coordinate compiler, but it does not create an affine-purity
divisor at the resolved terminal node.

## 8. Consequence for the logarithmic `ch_2` route

The factor audit sharpens the boundary-budget strategy:

1. do not spend more algebra trying to factor the target quintic—it is
   irreducible;
2. do not count the normalized `C^3` or `s^2` factors as source components;
3. use (4.1)--(4.3) as the finite carrier-contact gate after restoring all
   fixed-target parameters;
4. regard (7.3) as the already-accounted terminal packet; and
5. reserve any remaining localized `ch_2` point budget for the unresolved
   upstream attachment, outgoing tail, affine-purity row, or uncompiled
   centers.

The three carrier equations detect contact with the normalized target
curve.  They do not by themselves compute a normalization/conductor or
`Fitt_1` length of the source boundary pullback.  That local module still
requires the missing source factor or an equivalent completed logarithmic
matrix.

The subsequent
[`coprime carrier-jet discriminant pattern`](COPRIME_CARRIER_JET_DISCRIMINANT_PATTERN.md)
explains why this packet misses an invariant by exactly one order.  A
primitive `(m,n)` slice has `m+n-4` normalized coefficients and the carrier
transport has three parameters, so the raw saturation order is
`N_*=m+n-1`; here `N_*=7`, exactly the number of prescribed centers.  The
first raw invariant occurs at jet eight, which the present contact test does
not constrain.  On `Res(p',q')=0`, the generic target packet is one cusp
plus three nodes, the conductor splits `2+6`, and the raw jet map has corank
one.  The unibranch attachment theorem shows that this rank loss is not the
point length: a minimal transverse boundary attachment over the ordinary
cusp contributes `2q_p`, with `q_p` its local residue index; a complete
residue-degree-`f` fiber contributes `2f` when all attachments are minimal.
Locating those attachments and verifying the hypotheses remain open.

<!-- status-consumer: PCJDP1 d4c16bb71dfc6b80 -->
<!-- status-consumer: LUAF1 b0279670ffbd3fa5 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_carrier_jet_factorization.py
```

The command checks (2.3), the seven coefficients (3.4), the inverse (3.5),
the three equations (4.1)--(4.3), their weights and complete-intersection
form, the fixed-coordinate transport (5.5)--(5.12), the terminal factor
(7.3), including its arbitrary
leading scale, and the disappearance of the translation-dependent origin
factors.  On the `E_6+A_1` stratum it checks the closed binomial jet formula,
the triangular recovery (6.9), and the prime fixed-coordinate complete
intersection (6.10)--(6.12).  It also checks the scale-free `E_8` endpoint
complete intersection (6.18)--(6.19) and runs the explicitly bounded
625-target irreducibility regression from (1.3).

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->
