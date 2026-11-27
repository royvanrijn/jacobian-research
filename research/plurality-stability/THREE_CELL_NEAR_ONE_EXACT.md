# Exact three-cell simplex optimality near correlation one

There is a $\rho_*<1$ such that for every $\rho_*<\rho<1$ and
every $d\ge2$, the regular simplex maximizes noise stability among
measurable Gaussian partitions into three cells of measure $1/3$.
Equality forces the simplex, up to orthogonal transformations, labels
and null sets. The proof compares a clipped moving-junction graph and
the remaining exterior changes through a bounded annulus where the
actual interfaces are small normal graphs.

The argument uses the [near-one bounds](THREE_CELL_NEAR_ONE.md),
the [moving-junction estimate (MJ0)](THREE_CELL_MOVING_JUNCTION.md),
the [exterior tube estimates](THREE_CELL_EXTERIOR_TUBES.md), and
the [three-cell interface gap (S37)](SIMPLEX_INTERFACE_SYMMETRY.md).

## Relative Taylor remainder for a clipped moving tripod

Let $T$ be a central planar simplex, with any fixed orientation, $\rho=1-e^2/2$,
$\sigma^2=1-\rho^2$, and $a^2=(1+\rho)/2$. Write $q_T(x,y)=1_{\{T(x)\ne T(y)\}}$,
$w=(x+y)/2$, and $z=(y-x)/e$. Use outward ray coordinates $r_k=x\cdot\tau_k,\
t_k=x\cdot n_k$ and cyclically oriented normals. Fix two small absolute angular
widths $0<\eta_0<\eta_1$, and a smooth scalar cutoff $\zeta$ equal to one on
$[-\eta_0,\eta_0]$, zero off $[-\eta_1,\eta_1]$. The three positive-ray angular
supports are disjoint, and their only incident reference labels are those of the
corresponding ray.

Write $\phi_a$ for the centered planar Gaussian density with covariance
$a^2 I$, and $\phi,\Phi$ for the scalar standard density and distribution
function. The norm $\|\cdot\|_\nu$ of a ray field sums over all three rays.

Let $R_*>0$ be fixed. For each ray let $h_k(0)=0$,
$\operatorname{Lip}(h_k)\le\delta$, and $h_k$ be constant on $r\ge R_*$ (the
constants may differ). Let $|p|\le\delta$. Define
\[
 V(x)=\sum_{k:r_k>0}n_k h_k(r_k)\zeta(t_k/r_k),\qquad
 F_\lambda(x)=x+\lambda(p+V(x)),\qquad 0\le\lambda\le1.
\]
Set $C=F_1(T)$ and $f_k=h_k+p\cdot n_k$. Thus the three new interfaces are
$p+r\tau_k+h_k(r)n_k$. Write $h^0_k$ for extension by zero to the negative
tangential half-line, $\mathcal D_\rho(g)=\langle
g,(I-T_\rho^{(1)})g\rangle_{\gamma_1}$, and
\[
 Z_e(h)=e^{-1}\sum_k\int_0^\infty e^{-c_0 r^2/e^2}|h_k(r)|^2 dr
\]
for a sufficiently small absolute $c_0>0$. All norms of half-line speeds below use
$d\nu=\phi(0)\phi(r)\,dr$. For small $\delta$ and all small $e$, the precise
relative remainder is
\[
 \left|N_\rho(C)-N_\rho(T)-\rho\,\mathfrak q_\rho(f)\right|
 \le C_{R_*}\delta\left(e|p|^2+e\|h\|_\nu^2+
                  e^{-1}\sum_k\mathcal D_\rho(h^0_k)+Z_e(h)\right).
 \tag{R1}
\]
The ambient translation $p$ is allowed to be much larger than $e$.
The quadratic form uses the full kernel in (B6), extended to
unbalanced speeds as in (S39). The volume constraint is used
when applying its positive gap.

## Differentiating the pulled-back pair density

Since $|h|\le\delta r$ and $|h'|\le\delta$ almost everywhere, each angular term of
$V$, extended by zero outside its tube, is Lipschitz with constant $C\delta$. In
particular $|V(x)|\le C\delta|x|$, and $F_\lambda$ is bi-Lipschitz for small
$\delta$. At a point in a tube the derivative of $V$ has a single image direction
$n_k$, so the exact area factor is
\[
 J_\lambda(x)=1+\lambda j(x),\qquad
 j(x)= h_k(r_k)\zeta'(t_k/r_k)/r_k.
\]
It equals one outside the tubes, and $|j|\le C\delta$. The support of $j$ lies in
the angular gluing regions $\eta_0r_k\le |t_k|\le\eta_1r_k$. There are no terms
$j(x)^2$ from a single endpoint.

By the area formula the new pair kernel evaluated on reference pairs is
\[
 e^{-2}J_\lambda(x)J_\lambda(y)
      \phi_a(w+\lambda m)\phi_1(z+\lambda \xi),
 \quad
 m=p+\tfrac12(V(x)+V(y)),\quad \xi=(V(y)-V(x))/e .
\]
In particular $p$ cancels from $\xi$ identically. Put $G_\lambda=\phi_a(w+\lambda
m)\phi_1(z+\lambda\xi)$. We have $|m|\le C\delta(1+|w|+e|z|)$ and $|\xi|\le C\delta
|z|$. The two deformed arguments have combined squared norm at least
$c(|w|^2+|z|^2)-C$, uniformly in $0\le\lambda\le1$. For $1\le j_0\le3$, Gaussian
differentiation consequently gives
\[
 |G_\lambda^{(j_0)}|
 \le C (|m|+|\xi|)^{j_0}(1+|w|+|z|)^{j_0}
                    \exp[-c(|w|^2+|z|^2)].
\]
Indeed each directional derivative inserts a polynomial of degree at most one more
than the preceding derivative and a factor in $m,\xi$, with coefficients uniformly
bounded since $a$ is bounded away from zero. The deformed arguments have the
indicated polynomial bound in $w,z$.

The third derivative of the full product, with $e^{-2}$ omitted, is exactly
\[
 J_\lambda(x)J_\lambda(y)G_\lambda'''
 +3(j(x)J_\lambda(y)+J_\lambda(x)j(y))G_\lambda''
 +6j(x)j(y)G_\lambda'.
\]
Use the small Lipschitz bounds in one factor of each term and absorb the polynomial
into a wider Gaussian. The resulting absolute bound is
\[
 C\delta\,\mathcal K_e(x,y)
                 (|m|^2+|\xi|^2+|j(x)j(y)|),\qquad
 \mathcal K_e=e^{-2}\exp[-c(|w|^2+|z|^2)]
 \tag{R2}
\]
for a fixed absolute $c>0$ after decreasing it. The following three estimates
integrate (R2) only on differing-label pairs.

**Midpoints.** For fixed $z$, a differing-label pair must cross a supporting line,
so its midpoint lies in one of three strips of width $e|z|$. Their Gaussian weight
is at most $C\min(e|z|,1)$, giving $\iint q_T\mathcal K_e\le Ce$. This controls the
$p^2$ term. At a point $x$ in an angular support, the distance to a different
reference label is at least $c|t_k|$. Since $w=x+ez/2$,
\[
 \int q_T(x,y)\mathcal K_e(x,y)\,dy
              \le C\exp(-c'|x|^2-c't_k^2/e^2).
\]
This uses $|x-y|\ge c|t_k|$ for differing labels, and the Gaussian tail in the $z$
integral. Multiply by $h_k(r_k)^2$, integrate in the normal coordinate to gain a
factor $e$, and sum the rays. For $r\le R_*$ the widened marginal $e^{-c'r^2}$ is
bounded by $C_{R_*}\phi(r)$. For $r>R_*$ the height is constant, so the widened
integral is still bounded by $C_{R_*}$ times the usual Gaussian norm there. We
obtain $\iint q_T\mathcal K_e |m|^2\le C_{R_*}e(|p|^2+\|h\|_\nu^2)$.

**Differences, pairs in the same pure tube.** In the pure tube $|t|\le\eta_0 r$,
$\zeta=1$. If both endpoints lie in the same pure tube and have different labels,
their normal coordinates have opposite signs. The displacement difference is
$n_k(h_k(s)-h_k(r))$. The integral in the two normal coordinates of
$\exp[-c(t-u)^2/e^2]$ over $tu\le0$ is at most $Ce^2$: on each of the two quadrants
set $v=|t|+|u|$ and integrate $v e^{-cv^2/e^2}$. Including the $e^{-2}$ kernel
prefactor and the $e^{-2}$ from $|\xi|^2$, the integral on these pairs is at most
\[
 C e^{-2}\iint_{r,s>0}
 e^{-c((r+s)^2+(r-s)^2/e^2)}(h_k(r)-h_k(s))^2\,dr\,ds .
 \tag{R3}
\]
Here constants in the two exponents may be decreased.

For completeness, the nonlocal weighted comparison used in (R3) is as follows. Let
$L>1$ be a sufficiently large fixed integer and $\rho'=1-L(1-\rho)$, so the noise
length is $\sqrt L e$. The Gaussian tangential joint density $p_{\rho'}(r,s)$ has
midpoint exponent $-(r+s)^2/[4(1+\rho')]$, difference exponent $-(r-s)^2/(2Le^2)$,
and prefactor comparable to $e^{-1}$. Wherever $h_k(r)\ne h_k(s)$ with $r,s>0$,
either $r$ or $s$ is at most $R_*$, so $(r+s)^2\le C R_*^2+C(r-s)^2$. Choose
$1/(2L)<c/2$, then choose $e$ small. The difference exponent in (R3) absorbs the
cost of changing the midpoint exponent, and its kernel is bounded by
$C_{R_*}e^{-1}p_{\rho'}(r,s)$. Extend by zero and enlarge the domain to the full
tangential line. On Hermite degree $n$,
\[
 1-(\rho')^n=(1-\rho')\sum_{j=0}^{n-1}(\rho')^j
             \le L(1-\rho)\sum_{j=0}^{n-1}\rho^j=L(1-\rho^n).
\]
Thus (R3) is at most $C_{R_*}e^{-1}\mathcal D_\rho(h_k^0)$.

**Differences, remaining pairs.** Use $|V(y)-V(x)|^2\le2|V(x)|^2+2|V(y)|^2$. If
$V(x)\ne0$ and the endpoints are not both in that same pure tube, then every
differing-label pair of this type satisfies $|x-y|\ge c r_k(x)$. If
$|t_k(x)|\ge\eta_0 r_k(x)/2$, the phase depth gives that distance. Otherwise the
distance from $x$ to the complement of the pure tube is at least $c r_k(x)$. The
latter alternative includes the distance to its vertex. Conditional integration of
$\mathcal K_e$ in $y$ gives $Ce^{-c r_k(x)^2/e^2}$. The normal support at $x$ has
length $O(r_k(x))$. Consequently the contribution is at most $C e^{-2}\sum_k\int r
e^{-c r^2/e^2}|h_k(r)|^2dr\le C Z_e(h)$.

**Jacobian product.** Both endpoints are now in angular gluing regions. Their phase
depths and differing labels imply $|x-y|\ge c\max(r_k(x),r_\ell(y))$. Hence
$\mathcal K_e\le C e^{-2}\exp[-c(r_k(x)^2+r_\ell(y)^2)/e^2]$. Integrating each
normal coordinate cancels $1/r$ in its corresponding $j$. This bounds the integral
by
\[
 C e^{-2}\left(\sum_k\int_0^\infty e^{-c r^2/e^2}|h_k(r)|\,dr\right)^2
 \le C Z_e(h).
\]
The last inequality is Cauchy–Schwarz with a Gaussian integral of order $e$. No
bound on a point trace of $h/r$ in terms of energy is used.

Taylor's theorem with (R2) and these estimates gives (R1) with the first and second
variations in place of $0$ and $\rho\mathfrak q_\rho(f)$. For smooth vector fields
those are the standard interface variations: every pair score ties on the
corresponding reference facet, so the first variation and the acceleration term in
the second variation vanish. The half second derivative is $\rho\mathfrak
q_\rho(f)$, with no mass constraint in that identity. For the present Lipschitz
field, convolve it with compact smooth mollifiers. The mollifications converge
uniformly, preserve Lipschitz and linear-growth bounds, and their derivatives
converge almost everywhere. One may also multiply by a smooth spatial cutoff equal
to one on a ball of increasing radius, with support in the doubled ball and
inverse-radius gradient bound. This gives smooth compactly supported fields with
uniform Lipschitz and linear-growth bounds. The pullback derivatives at zero up to
second order converge by the same Gaussian domination. Facet speeds converge
pointwise with dominated Gaussian bounds, so the interface identity at fixed $\rho$
passes to the limit. This proves (R1).

## Relating the remainder to sector energies

Use the standard and sign squared-difference energies $E_+(P_{\rm st}f),E_-(P_{\rm
sg}f)$ of (S18), with the $1/2$ normalization of those formulas. Write $E$ for their
sum. Translation adds a standard constant $q_k=p\cdot n_k$, and leaves these
differences unchanged. With the reflected kernels $N=M_1+M_{-1}$, $D=M_1-M_{-1}$,
and $c_\rho=\phi(0)/\sigma$, denote their even and odd ground-state energies by
$E_N,E_D$ as in (S31). The full-line zero extension satisfies the exact identity
\[
 \frac{\phi(0)^2}{\sigma}\mathcal D_\rho(h^0)
   =c_\rho\|h\|_\nu^2-\langle h,M_1 h\rangle_\nu
   =E_N(h)+\langle h,M_{-1}h\rangle_\nu .
\]
In addition,
\[
 E_N(h)-E_D(h)=(1-\rho)c_\rho\|h\|_\nu^2
                                  -2\langle h,M_{-1}h\rangle_\nu .
\]
These follow from the row identities for $N$ and $D$ in (S30); the first uses the
independent normal joint density at $(0,0)$. Since $M_{-1}\ge0$ pointwise and its
row is $c_\rho\Phi(-\rho r/\sigma)$, the absolute value of its pairing is at most $C
Z_e(h)$ by $2|h(r)h(s)|\le h(r)^2+h(s)^2$. Then (S32) controls $E_N$ in the standard
sector and $E_D$ in the sign sector. Thus
\[
 e^{-1}\sum_k\mathcal D_\rho(h^0_k)
       \le C(E+e\|h\|_\nu^2+Z_e(h)).
 \tag{R4}
\]
The constants here are absolute.

## Cutoff rotation gauge and balance

One may take the reference $T_\theta$ itself rotated before forming the two cutoffs,
and use the same $T_\theta$ in the exterior construction. Fix $R$. In the chart
centered at the actual vertex $p$ use the actual heights $h_k^\theta(r)$, with
central tangential coordinate $s=r+p\cdot\tau_k^\theta$. Choose $\chi=1$ for $s\le
R$, $\chi=0$ for $s\ge2R$, and a complementary mask $\psi$ with $\chi^2+\psi^2=1$,
both nonnegative. The inner total normal speed in these coordinates is
\[
 f_{C,k}^\theta(r)=\chi(r+p\cdot\tau_k^\theta)
                      (h_k^\theta(r)+p\cdot n_k^\theta).
\]
The inner shifted height equals $f_C^\theta-p\cdot n^\theta$, so it has constant
tails. At fixed $r$, differentiation of the change of coordinates from one ray to
its rotation gives
$\partial_\theta h_k^\theta=-r-h_k^\theta\partial_r h_k^\theta$
with the chosen sense of rotation. On the bounded chart this is $-r+O_R(\delta^2
r)$. Differentiating the sign rotation moment $\sum_k\int r
f_{C,k}^\theta(r)d\nu(r)$ gives
\[
 -\sum_k\int r^2\chi(r+p\cdot\tau_k^\theta)d\nu(r)
                         +O_R(|p|+\delta^2).
\]
The main term is bounded above by a negative absolute constant for $R\ge2$ because
$\chi=1$ on a fixed inner interval. The moment at the initial orientation is
$O_R(\delta+|p|)$. It follows that an $O_R(\delta+|p|)$ change of common orientation
imposes the exact orthogonality to the sign rotation mode. The exterior cutoff is
constructed with that same orientation; actual tube containment persists under this
small rotation when the barrier is chosen with an angular margin.

For an unbalanced inner $C$, let $b=P_{\rm st}\int f_Cd\nu/\int d\nu$ be its
standard constant projection. Linearized volume change on standard constants is an
invertible fixed map onto centered label masses. Expanding the masses of the
translated simplex and integrating its signed strips gives
\[
 |b|\le C|\boldsymbol\gamma(C)-(1/3)\mathbf 1|
                    +C\delta(|p|+\|h_C\|_\nu).
 \tag{R5}
\]
For the strip expansion the tangential weight is $\phi(r+p\cdot\tau_k)$, whose ratio
to $\phi(r)$ differs from one by $O(|p|)$ in $L^2(\gamma_1)$. The normal integral
differs from $\phi(0)h_k(r)$ by at most $C (|p|+|h_k(r)|)|h_k(r)|$, using the
Lipschitz property of $\phi$. Since $|h_k(r)|\le\delta r$, integration gives the
asserted error by Cauchy–Schwarz. The translation of the reference masses has
$O(|p|^2)$ Taylor remainder and linear normal speed $p\cdot n$.

The centered standard mean projection $f_C-b$ is balanced and has the same
squared-difference energy $E$. With the above rotation gauge, (S37) supplies an
absolute $ce$ gap on its $L^2(\nu)$ norm for $\rho$ close to one. In (S18) the
standard row potential is bounded by an absolute $Ce$ times the squared norm, and
the sign sector is exactly a nonnegative squared-difference form. Combining these
two inequalities gives $\mathfrak q_\rho(f_C-b)\ge c(E+e\|f_C-b\|_\nu^2)$ with
absolute $c>0$. Replacing $f_C-b$ by $f_C$ changes only the standard row potential,
because a standard constant has zero edge difference. Bound its mixed term by
$\varepsilon e\|f_C-b\|_\nu^2+C_\varepsilon e|b|^2$. This yields
\[
 \mathfrak q_\rho(f_C)\ge c(E+e\|f_C\|_\nu^2)-C e|b|^2
 \tag{R6}
\]
with absolute $c,C$. Equations (R1), (R4), (R5), (R6) give the needed inner estimate
once the anchored $Z_e(h_C)+e(|p|^2+\|h_C\|_\nu^2)$ is bounded by
$C_R(E+e\|f_C\|_\nu^2+\text{relative exterior error})$. The moving-junction estimate
supplies this bound, using the localized score residual from replacing $A$ by $C$.
In particular the positive constants that must compete with the annular localization
error come from (S37) and (R6); the fixed-$R$ constants enter only the coefficients
multiplied by the small deformation parameter $\delta$.

For the quadratic annular IMS cutoff error itself, take $\chi=\cos\vartheta,\
\psi=\sin\vartheta$ with $|\vartheta'|\le C/R$ supported on $[R,2R]$. The pair
identity for scaled heights has error
$g(r)g(s)[(\chi(r)-\chi(s))^2+(\psi(r)-\psi(s))^2]$. For pairs with both tangential
coordinates in $[R/2,3R]$, its absolute integrated value with the binary quadratic
edge coefficient of order $e^{-1}$ is at most
\[
 C e^{-1}R^{-2}\int_{R/2}^{3R}|g(r)|^2
         {\mathbb E}[(S-S')^2\mid S=r]\,d\gamma_1(r)
 \le C e R^{-2}\int_{R/2}^{3R}|g(r)|^2\,d\gamma_1(r).
\]
Here the exact conditional moment is $\sigma^2+(1-\rho)^2 r^2$, bounded by an
absolute $Ce^2$ for small $e$ at every fixed $R$. Thus the IMS coefficient is
absolute and can be made small by choosing $R$ large before sending $\rho\uparrow1$.
Pairs outside this local class need the separate Gaussian-separated mixed
comparison.

## The clipped competitors and their energies

Take a sequence of actual balanced planar maximizers $A$ with
$\rho\uparrow1$. After subsequence and alignment, the simplex $T$
is the limit from (P3). Fix a sufficiently narrow tube aperture for
the exterior comparison. By (P17)–(P18) with a smaller aperture,
all large-radius changes are in those tubes, also after an $o(1)$
rotation of $T$. Choose a fixed large $R$ beyond the confinement
radius. By bounded-region graph alignment, on the three ray
rectangles with $s\in I=[R/2,6R]$ and $|t|$ small but fixed
the only changes are the signed binary strips of heights $g_k(s)$.
Their heights and slopes tend uniformly to zero. The chart about
the actual vertex $p$ extends throughout this bounded region,
with $p\to0$. There are no additional changes in the bounded
interior away from the strips and the moving vertex, by bounded
winner-set convergence.

In the central tangential coordinate $s$ choose a smooth cutoff
angle equal to zero for $s\le R$, to $\pi/2$ for $s\ge2R$, and
with derivative bounded by $C/R$ for an absolute $C$. Let
$\chi,\psi$ be its cosine and sine. In the three tubes the inner
partition $C$ makes the $\chi g_k$ graph changes on $[R,2R]$
and agrees there with $T$ beyond $2R$. The exterior partition
$O$ makes the $\psi g_k$ graph changes on $[R,2R]$ and takes
all the actual changes in those tubes beyond $2R$. Off the
exterior tube portions $s\ge R$, take $C=A$ and $O=T$.
The definitions give $C=T$ outside a fixed ball. On $I$,
$C,O$ have respective heights $\chi g_k,\psi g_k$.

In vertex coordinates the inner interfaces have heights $h_C$
obtained from $f_C=h_C+p\cdot n$ in the cutoff rotation formula.
In particular their endpoints at $r=0$ are at $p$, and the
heights have constant tails. The ambient shear in (R1) maps the
reference tripod rays to these interfaces. Since each is a
Lipschitz graph and the shear is one-to-one, it maps the reference
cells to the corresponding inner cells. Choose the common
orientation by the cutoff rotation gauge above, and use this
reference for both $C$ and $O$. All graph and containment
conclusions hold at that orientation.

Let $E_C=E(f_C)+e\|f_C\|_\nu^2$ and let $L_A,L_C,L_O$ be the
three regrets relative to $T$. Define $E_O$ as the sum of the
following nonnegative quantities for the three tubes:

- $e$ times the sum of the two post-merge squared profile norms $V$.
- The sum of the two profile energies $\mathcal B$ of (CT2).
- The sum of the two sorting slacks integrated with $q_R$, and
  the third-label mass integrated with $q_R$.
- $e^2R^2 L_O$, counted once over all the tubes.

There is an absolute $c_o>0$ such that for sufficiently large fixed
$R$ and then sufficiently small $e$,
\[
 N_\rho(O)-N_\rho(T)\ge c_o E_O.
 \tag{R7}
\]
Indeed apply (CT2) and (CT3) to each of the two global binary
merges, retaining a fixed fraction of the $\mathcal B$ term.
Average the resulting sorting and diagonal estimates over the two
possible merges. The average third-label cut decrease gives the
$q_R$ mass term in (CT6). The compression estimate (CT10) also
applies with either of the two merges and retains the $e^2 R^2 L_O$
reserve. The regret error and the cross-ray interactions are
absorbed by (CT7)–(CT9) exactly as in the tube proof, using
$L_s^*\le Cp(s)^2/\sigma$ for each merge. The constants can be
chosen independently of $R$ once $R$ exceeds a sufficiently
large absolute threshold, using (CT3) with that threshold to
absorb the negative diagonal cost.

The two merged binary profiles are, up to signs, the signed
changes in the two incident fiber masses before merging. The
third mass change is their negative sum. Hence the squared norm
of the full signed volume change of $O$ is at most
$C\gamma_1([R,\infty))\sum V$, by Cauchy–Schwarz in each
tube. On the transition annulus the mass error in adding the
clipped changes has absolute fiber magnitude at most
$C|\psi g_k|$. The original-strip minus inner-strip length
is $(1-\chi)|g_k|\le\psi|g_k|$, and the outer strip has
length $\psi|g_k|$. The outer
binary profile on the transition annulus has magnitude at least
$c|\psi g_k|$ since the graphs are uniformly small. Balance of
$A$ therefore gives
\[
 e|\boldsymbol\gamma(C)-(1/3)\mathbf 1|^2
    \le C\gamma_1([R,\infty)) E_O .
 \tag{R8}
\]

## The anchored core estimate

The transition cuts preserve a fixed fraction of the regret. At
each fiber one of the cuts retains at least $|g_k|/\sqrt2$. The
pointwise margin on a graph strip is comparable to
$\min(|t|/\sigma,1)$, by (P14) and the upper bound $C/\sigma$
on the score derivatives. The normal Gaussian density is
comparable to a constant on these strips. Integration in the
normal coordinate gives $L_A\le C(L_C+L_O)$, since outside the
transition one clipped partition makes the actual transfer.

The inner regret obeys
\[
 L_C\le C_R e^{-1}(|p|^2+\|h_C\|_\nu^2).
 \tag{R9}
\]
The ball of radius a sufficiently large fixed multiple of $|p|$
about the origin contributes at most $C|p|^2$. Outside it the
defects from translation of the reference are in strips of width
$C|p|$ along the incident lines. Their area is at most
$C|p|$, and their regret per point at most $C|p|/\sigma$.
For the inner graph strips in vertex coordinates, the distance
from the central pair line is at most $|p|+|h_{C,k}(r)|$.
For small slopes the three strips stay separated outside the
small ball and all changes there are between the corresponding
incident labels. Integrate the upper pointwise margin against
the Gaussian strip area bound $C_R\phi(r)|h_{C,k}(r)|$ on the
bounded support of the actual defect. This proves (R9).

Apply (MJ0) on a fixed inner chart, for example $0\le r\le R/4$,
where $C=A$. Its local edge forms are bounded by $C_R E(f_C)$
since the translation adds a standard constant. Extending $Z_e$
beyond $R/4$ costs at most $C_R e^{-c_R/e^2}\|h_C\|_\nu^2$,
using the bounded region and constant tails. Also
$e\|h_C\|_\nu^2\le C(e\|f_C\|_\nu^2+e|p|^2)$.
The area $d$ in (MJ0) obeys $d^2\le C(L_A+|p|^2)$, by (P14)
and the translation area bound. Use $L_A\le C(L_C+L_O)$,
(R9), and $e^2R^2L_O\le E_O$. Absorbing the exponentially
small inner terms yields
\[
 Z_e(h_C)+e(|p|^2+\|h_C\|_\nu^2)
     \le C_R E_C+C_R e^{-c_R/e^2}E_O.
 \tag{R10}
\]
In particular every exponentially small multiple of
$L_A+L_C+L_O$ is $o_R(1)(E_C+E_O)$, where $o_R(1)$ tends to
zero along the maximizing sequence at fixed $R$.

By (R1), (R4)–(R6), (R8) and (R10),
\[
 N_\rho(C)-N_\rho(T)
   \ge c_c E_C-\bigl(C\gamma_1([R,\infty))+o_R(1)\bigr)E_O,
 \tag{R11}
\]
with absolute $c_c>0$ and $C$. In fact the relative Taylor
error and the squared linearized volume error are bounded by
$C_R(\delta+\delta^2) E_C+o_R(1)E_O$, and $\delta$ tends to
zero at fixed $R$ by bounded graph alignment.

## Localizing the comparison in the graph annulus

Write $F(U)=N_\rho(U)-N_\rho(T)=2L_U-Q(1_U-1_T)$.
The cutoff error satisfies
\[
 |F(C)+F(O)-F(A)|
 \le \bigl(C/R^2+o_R(1)\bigr)(E_C+E_O)
 \tag{R12}
\]
with absolute $C$. Use a binary comparison with the full flat
halfspace for each ray. Its original, inner and outer heights
are $g_k,\chi g_k,\psi g_k$ on $I=[R/2,6R]$, and zero off $I$.

For the linear terms the cutoff commutator vanishes outside
the transition $[R,2R]$. There the tripod and halfspace regret
differ by a relative $C\exp(-cR^2/e^2)$, by the conditional
third-wedge bound preceding (CT7). The sum of those errors is
at most $C\exp(-cR^2/e^2)(L_A+L_C+L_O)$.

For the quadratic terms compare on three disjoint rectangles
with tangential window $I$ and normal width $\kappa R$ for a
fixed small $\kappa>0$. The reference and the actual changes
in each rectangle are the binary model, and the auxiliary
binary defects have no support outside it. For both endpoints
in one rectangle the quadratic commutators agree. Every
remaining noncanceling pair either has a transition endpoint
or joins the two distinct pure cutoff parts. A transition
defect is at distance at least $cR$ from the complement of
its rectangle. A pure inner defect has radius at most
$\sqrt{1+\eta^2}R$ by confinement and $s\le R$ in the tubes,
whereas a pure outer defect has radius at least $2R$.
The unaccounted pairs thus have separation at least $cR$.
All three defect vectors are either zero or the original
transfer at each point, since clipping only shrinks a binary
strip. Extract a Gaussian factor at this separation and use
the thin-strip and deep-transfer bounds of (P14) as in (CT11).
This bounds the residual by $C\exp(-cR^2/e^2)L_A$. By (R10)
the exponential errors are $o_R(1)(E_C+E_O)$.

Apply (P12) to the auxiliary binary commutators. If one
tangential coordinate is outside $I$, a nonzero commutator
requires the other to be in the transition. The normal change
against zero height is bounded by $C|g_k|^2/e$ from (P12),
and the conditional tangential probability is exponentially
small at the fixed separation. Also
\[
 e\sum_k\int_I g_k^2\,d\gamma_1\le C(E_C+E_O).
 \tag{R13}
\]
Indeed $g_k^2=(\chi g_k)^2+(\psi g_k)^2$. The first term
is controlled by the inner $L^2$ norm after translation
of the tangential coordinate by $p\cdot\tau_k$. On $I$
the resulting Gaussian weights differ by at most a factor
two for all sufficiently small $p$ at fixed $R$. The
second term is controlled by the binary profile norm in
$E_O$, since $|\Phi(u)-1/2|\ge c|u|$ on the bounded
small heights after clipping.

On $I^2$, the normal diagonal term in (P12) has expansion
$D_\rho(u)=d_\rho u^2+O(eu^4)$ with $|d_\rho|\le Ce$.
For example differentiate the two equal-threshold probabilities
to get $D_\rho'(u)=2\phi(u)[1-2\Phi((1-\rho)u/\sigma)]$,
and expand at zero. The quadratic diagonal contributions cancel
by $\chi^2+\psi^2=1$, leaving $o_R(1)(E_C+E_O)$ by (R13).

In the double normal integral in (P12), replace the normal
pair density by its value at $(0,0)$. The corresponding
quadratic term in the commutator is that value times
\[
 g_k(s)g_k(t)\bigl((\chi(s)-\chi(t))^2+(\psi(s)-\psi(t))^2\bigr).
\]
The bound following (R6) applies on $I$ using the same
conditional moment, with an absolute constant for small $e$
at fixed $R$. By (R13), the integrated error is at most
$C R^{-2}(E_C+E_O)$.

It remains to estimate the error in the density replacement.
Let $G$ be any of the three heights on $I$, with height and
Lipschitz bounds $\varepsilon\to0$. At normal positions $a,b$
in the interval between $G(s)$ and $G(t)$, the density
relative to its value at the origin is
$\exp[-(a+b)^2/(4(1+\rho))-(a-b)^2/(4(1-\rho))]$.
Its difference from one is at most
$C\varepsilon^2(1+|s-t|^2/e^2)$. Integrating gives an error
bounded by $C\varepsilon^2/e$ times a tangential difference
integral with the pair density multiplied by
$1+|s-t|^2/e^2$.

At fixed $R$ this last integral is at most $C_R$ times
the same integral unmultiplied and restricted to $|s-t|\le e$
in $I^2$. To prove this, on the bounded interval the midpoint
density has positive upper and lower bounds depending only
on $R$. In a displacement band
$2^j e<|s-t|\le2^{j+1}e$, split the segment into
$m=2^{j+1}$ equal steps. Cauchy–Schwarz bounds the squared
difference by $m$ times the sum of the squared step
differences. Changing variables from $(s,t)$ to the first
endpoint of each step and the step length has Jacobian $m$.
All step endpoints remain in $I$ and are separated by at
most $e$. The Gaussian density and polynomial supply a
factor at most $C_R(1+4^j)e^{-c4^j}$ relative to the
near-diagonal density. The product with $m^3$ is summable.

The restricted difference integral for the inner profile,
divided by $e$, is bounded by $C_R(E_C+o_R(1)E_O)$ by
translation of the tangential coordinates and (R4), (R10).
For the outer profile it is at most $C\mathcal B$ on that
ray since the binary mass profile difference is comparable
to the graph difference and, by the small slope, has magnitude
at most $e$ on the restricted pairs. For the original profile
the quadratic cutoff identity bounds the integral by the sum
of the two clipped integrals plus the already bounded quadratic
cutoff term. Hence the full normal density replacement costs
at most $o_R(1)(E_C+E_O)$. This proves (R12).

## Exact comparison and equality

Combine (R7), (R11), and (R12) to obtain
\[
 F(A)\ge (c_c-C/R^2-o_R(1))E_C
     +(c_o-C/R^2-C\gamma_1([R,\infty))-o_R(1))E_O.
 \tag{R14}
\]
Choose $R$ fixed sufficiently large to leave a strictly
positive coefficient on both energies, and then take $\rho$
sufficiently close to one along the sequence. Since $A$
maximizes and the simplex is balanced, $F(A)\le0$. Both
energies therefore vanish. By (R10), $p=0$ and $h_C=0$,
so the inner partition is the reference. The outer partition
agrees with the reference by the regret term in its energy.
The two cutoffs together cover all changes, since at each
transition fiber at least one multiplier is positive. Thus
$A=T$ up to null sets. A sequence of nonsimplex planar
maximizers at correlations tending to one is impossible.

Heilman–Tarter, Lemma 2.3 and Theorem 1.9, give existence and
dimension reduction for these fixed masses at every positive
correlation. For three cells and ambient dimension $d\ge2$,
every maximizer is cylindrical over a planar maximizer up to
rotation and null sets. Thus the same interval and equality
conclusion hold in all such dimensions.
[Heilman–Tarter](https://stevenheilman.org/~heilman/papers/plurthreeB.pdf)

## A uniform quadratic penalty for unequal masses

There are absolute \(C>0\) and \(\rho_0<1\) such that for any ordinary
three-cell partition \(A\) in \(d\ge2\), \(\rho_0<\rho<1\), and
\(e=\sqrt{2(1-\rho)}\),
\[
 S_\rho(T)-S_\rho(A)\ge
 -C e\big|\boldsymbol\gamma(A)-(1/3)\mathbf1\big|^2 .
 \tag{R15}
\]
Here \(T\) is the balanced simplex. For a contradiction take
\(\rho_j\to1\) and partitions with the deficit less than
\(-j e_j |m_j|^2\), where \(m_j=\boldsymbol\gamma(A_j)-(1/3)\mathbf1\).
Write \(N_\rho=1-S_\rho\). The bounds \(N_\rho(A_j)\ge0\) and
\(N_\rho(T)\le C e_j\) imply \(|m_j|^2\le C/j\). The masses are
eventually positive. Replace each competitor by an actual maximizer
at its own masses and use Heilman–Tarter to work in the plane.

In (P2)–(P3) of the [near-one reduction](THREE_CELL_NEAR_ONE.md)
the noise bound is still \(N_{\rho_j}(A_j)\le N_{\rho_j}(T)\),
and limiting masses are \(1/3\) by \(L^1\) convergence and \(m_j\to0\).
The equal-mass perimeter equality forces an aligned tripod and the
same energy-measure convergence. For (P4) the two compact volume
variations still have independent limiting mass columns. They
repair to the masses of \(A_j\). Global maximality at those masses
gives the winning-score system by the linearized assignment argument.
Fixed-volume stationarity with the bumps gives the \(o(e_j)\)
multiplier bound.
The repairs for replacements and blow-ups in (P6)–(P8) use the
same target masses. Exterior confinement (P17)–(P18) uses this
stationarity and the limiting tripod. In (P18), identically equal
winning scores for two labels would make the OU image of their
indicator difference constant. Injectivity on the nonconstant
Hermite modes would make their difference constant, impossible for
disjoint cells both of positive mass. The
[moving-junction construction (MJ0)](THREE_CELL_MOVING_JUNCTION.md)
then applies with the same rank-two repair and the adjusted score tie.

For fixed large \(R\), (R1)–(R7) and (R9)–(R12) apply with one
assembly change. The signed-strip estimate preceding (R8) gives
\[
 e_j|\boldsymbol\gamma(C)-(1/3)\mathbf1|^2
 \le 2e_j |m_j|^2+C\gamma_1([R,\infty))E_O .
\]
It bounds \(\boldsymbol\gamma(C)-\boldsymbol\gamma(A_j)\) by the
same outer changes, and \(\boldsymbol\gamma(A_j)-(1/3)\mathbf1=m_j\).
The anchored core estimate already permits unequal masses, so
(R11) and (R14) acquire only the error \(-C e_j|m_j|^2\).
Fix \(R\) as in (R14) and let \(j\) grow. The two energy
coefficients are positive, giving the claimed contradiction.

## Discrete plurality on the near-one interval

Fix $\rho_*<\rho<1$ and $\epsilon>0$. There are an integer
$D$ and a $\tau>0$ such that, for every $m\ge1$, every
balanced rule $f:[3]^m\to\Delta_3$ under the uniform law with
$\operatorname{Inf}^{\le D}_i(f_j)\le\tau$ for every coordinate
$i$ and choice $j$ satisfies
\[
 \mathbb E\langle f(X),f(Y)\rangle
 \le \lim_{n\to\infty}\mathsf S_\rho(\operatorname{Plur}_{3,n})+\epsilon.
\]
Here $Y$ independently retains each vote from $X$ with
probability $\rho$ and resamples it uniformly otherwise.
Balance means $\mathbb E f_j=1/3$ for every $j$.
The degree-$D$ influence sums the squared product-Fourier
coefficients involving the coordinate and of degree at most
$D$. The ordinary influence condition with threshold
$\tau$ implies this hypothesis. Isaksson–Mossel, Theorem 7.1,
provides an ordinary Gaussian realization of the same masses
in dimension $2m$ with stability within $\epsilon$.
Their Lemma 7.3 identifies the plurality limit with the
planar simplex stability, including arbitrary fixed tie
breaking since count-tie probabilities tend to zero.
[Isaksson–Mossel](https://arxiv.org/pdf/0903.3362)

## Two endpoints of any strict failure interval

If the simplex is strictly suboptimal at some positive
correlation, each connected component $(a,b)$ of the strict
failure set has $21/50<a<b<1$. Both $a$ and $b$ have
nonsimplex balanced global maximizers tied with the simplex.
Indeed a strict improvement by a fixed partition persists by
continuity, so the strict failure set is open. The known
small-correlation interval and the present near-one interval
confine its components away from zero and one. At each
endpoint the simplex is optimal. Take improving planar
maximizers along a sequence in the failure component tending
to that endpoint. By the positive-correlation compactness
argument they have a subsequential ordinary maximizer limit.
If the limit were a simplex, the fixed-correlation local growth
(B11), with (S37), would rule out the sequence. The argument
is two-sided in correlation. Finally the equality result at
$21/50$ excludes that value as an endpoint. This implication
addresses strict failure intervals and leaves open nonsimplex
ties in the optimality set.
[Continuation and local growth](POSITIVE_CORRELATION_CONTINUATION.md)
[Equality at the certified lower interval](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md)
