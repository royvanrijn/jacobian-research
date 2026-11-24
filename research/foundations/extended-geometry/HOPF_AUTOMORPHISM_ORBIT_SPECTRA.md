# Exact degree spectra of Hopf-quadric automorphism orbits

Let \(K\) be a field of characteristic zero and
\[
A_K=K[x,y,t]/(t^2+2xy-1).
\]
Write \(\deg_A f\) for the least total degree of a polynomial representative
of nonzero \(f\in A_K\). All automorphisms below are polynomial
\(K\)-algebra automorphisms. Every assertion remains valid after every
coefficient-field extension, including for automorphisms defined only over
that extension.

The retained seed families are
\[
\begin{aligned}
f_4&=(1+x)(y-t^2-xt^2/2),\\
xg_j&=(1+x)^4(t+1-x^j)^3-4(1+x)^2(t+1-x^j),\qquad j\ge1,\\
p_{n,\lambda}&=\frac{tC_{n,\lambda}(x)}x
                       (1-t^2U_{n,\lambda}(x)^2),
                       \qquad n\ge1,\quad\lambda\in K^*,\\
C_{n,\lambda}(x)&=(1+\lambda x)^{n+1},\\
U_{n,\lambda}(x)&=\frac{1+\lambda x}{2n}
                  \sum_{r=0}^{2n-1}\binom{2n}{r+1}(\lambda x)^r .
\end{aligned}                                                     \tag{1}
\]
Here \(g_1\) is the sextic seed; \(j\) is a parameter, not its degree.
The \(p_{n,\lambda}\) in this statement are the cubic-height hyperbolic
family; higher odd starting heights are not included in that notation.

**Theorem (full orbit spectra and minimum degrees).** For every
\(\sigma\in\operatorname{Aut}_K(A_K)\), put
\[
a=\deg_A\sigma(x),\qquad b=\deg_A\sigma(t),\qquad
h_j=\deg_A\sigma(t+1-x^j).
\]
These are positive integers, and
\[
\deg_A\sigma(f_4)=2(a+b),\qquad
\deg_A\sigma(g_j)=3(a+h_j),\qquad
\deg_A\sigma(p_{n,\lambda})=5na+3b.                    \tag{2}
\]
Consequently:

1. The full orbit of \(f_4\) has degree spectrum \(\{4,6,8,\ldots\}\).
2. All \(g_j\) lie in the orbit of \(g_1\), whose degree spectrum is
   \(\{6,9,12,\ldots\}\).
3. The minimum degree of the orbit of \(p_{n,\lambda}\) is \(5n+3\).
   Two such seeds \(p_{n,\lambda}\) and \(p_{m,\mu}\) are equivalent
   under an automorphism and multiplication by a nonzero scalar if and
   only if \(n=m\), over every common coefficient-field extension.

These families give infinitely many inequivalent pure unsafe phases with
unbounded minimum orbit degree. None of their automorphism images has least
degree five or seven. This does not exclude a witness in another orbit or
classify every pure unsafe phase.

## 1. The degree function and regularity

The associated graded ring is
\[
\operatorname{gr}A_K=K[x,y,t]/(t^2+2xy).                \tag{3}
\]
Indeed, the leading homogeneous part of a nonzero multiple
\(H(t^2+2xy-1)\) is the leading part of \(H\) times \(t^2+2xy\).
The latter quadratic has rank three, whereas a product of two linear
forms has rank at most two. It is irreducible even over an algebraic
closure. Thus (3) is a domain after every field extension, giving
\[
\deg_A(uv)=\deg_Au+\deg_Av\quad(u,v\ne0).               \tag{4}
\]
A sum of terms of different degrees has the maximum degree. If \(u\) is
nonconstant and \(B\in K[X]\) has degree \(r\), its leading term strictly
dominates the others, so
\[
\deg_A B(u)=r\deg_Au.                                 \tag{5}
\]
The degree-zero part is \(K\). An automorphism cannot send a nonconstant
element to a constant, since its inverse fixes \(K\). This proves
\(a,b,h_j\ge1\), even when leading terms in
\(\sigma(t)-\sigma(x)^j\) cancel.

The definitions in (1) are regular in \(A_K\). For \(g_j\), its numerator
at \(x=0\) is \((t^2-1)(t+3)\). Subtract this before dividing by \(x\);
replace the remaining \((t^2-1)/x\) by \(-2y\).
For the hyperbolic seed, \(U_{n,\lambda}(0)=1\), so
\[
p_{n,\lambda}=tC_{n,\lambda}(x)
 \left(2y-t^2\frac{U_{n,\lambda}(x)^2-1}{x}\right)      \tag{6}
\]
is a polynomial representative. Here \(\deg C_{n,\lambda}=n+1\) and
\(\deg U_{n,\lambda}=2n\).

The quotient embeds in \(K[x,x^{-1},t]\): the relation is primitive and
linear in \(y\), with coprime coefficients \(2x,t^2-1\), so the quotient
is a domain and localization at \(x\) is injective. All power conclusions
below concern \(A_K\), not merely the localization.

## 2. The quartic orbit

In the quotient one has
\[
2xf_4=(1+x)-(1+x)^3t^2.                               \tag{7}
\]
After applying \(\sigma\), the two right-hand terms have degrees \(a\)
and \(3a+2b\). They differ strictly. By (4),
\(a+\deg_A\sigma(f_4)=3a+2b\), proving the first formula in (2).

For every \(q\in K[x]\), the substitutions
\[
\tau_q(x)=x,\qquad \tau_q(t)=t+xq(x),\qquad
\tau_q(y)=y-tq(x)-xq(x)^2/2                            \tag{8}
\]
preserve \(t^2+2xy\), and direct substitution gives
\(\tau_q^{-1}=\tau_{-q}\). For \(q=x^{r-1}\), \(r\ge1\), one has
\(a=1,b=r\), including \(r=1\). Formula (2) gives degree \(2r+2\).
This realizes every even degree at least four without a classification
of the automorphism group.

## 3. The sextic orbit, including cancellations

Write \(V=\sigma(t+1-x^j)\). The transformed numerator identity for \(g_j\)
has terms of degrees \(4a+3h_j\) and \(2a+h_j\), whose difference is
\(2a+2h_j>0\). Hence
\[
a+\deg_A\sigma(g_j)=4a+3h_j.
\]
This proves its formula in (2). Retaining the actual degree \(h_j\) is
essential: \(\max(b,ja)\) can overestimate it when leading terms cancel.

Choose \(q=1-x^{j-1}\) in (8). Then
\[
\tau_q(t+1-x)=t+1-x^j,\qquad \tau_q(g_1)=g_j.           \tag{9}
\]
For the identity automorphism, \(a=1,h_j=j\), including \(j=1\).
Thus \(\deg_Ag_j=3j+3\), realizing every multiple of three at least six.
Formula (2) excludes every other degree in that orbit.

## 4. Hyperbolic minimum degrees and inequivalent orbits

For \(X=\sigma(x)\), \(T=\sigma(t)\), equation (5) gives
\[
\deg_A C_{n,\lambda}(X)=(n+1)a,\qquad
\deg_A U_{n,\lambda}(X)=2na.
\]
Therefore
\[
\deg_A(1-T^2U_{n,\lambda}(X)^2)=2b+4na.
\]
The numerator identity in (1) and (4) yield
\[
a+\deg_A\sigma(p_{n,\lambda})=b+(n+1)a+2b+4na.
\]
This is the last formula in (2). Since \(a,b\ge1\), the minimum is
\(5n+3\), attained by the identity.

The minimum degree of a full automorphism orbit is invariant under
automorphisms and nonzero scalar multiplication. Distinct \(n\) therefore
give inequivalent phases. Conversely, the diagonal automorphism
\[
\delta_c(x)=cx,\qquad \delta_c(y)=c^{-1}y,\qquad \delta_c(t)=t
\]
gives \(\delta_c(p_{n,\lambda})=c^{-1}p_{n,c\lambda}\).
Taking \(c=\mu/\lambda\) proves the converse. The same proof in \(A_L\)
for any extension \(L/K\) shows that inequivalence persists after extending
constants. The quartic and sextic orbit minima are four and six, so those
orbits are distinct from each other and from the hyperbolic orbits.

## 5. Moment and original Image interpretation

The normalized formal Hopf functional is
\[
\Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)\,dt.                \tag{10}
\]
Integration is coefficientwise polynomial integration. A phase is pure
when \(\Lambda(f^m)=0\) for all \(m\ge1\); it is unsafe when one fixed
polynomial marker \(g\) has \(\Lambda(gf^m)\ne0\) for infinitely many \(m\).

The retained all-order inputs are:

- [The quartic proof](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md), Theorem 1.1
  and sphere normalization, gives
  \(\Lambda(xf_4^m)=m!/(2m+1)!!\ne0\) for every \(m\ge1\).
- [The off-center proof](OFFCENTER_CUBIC_ALL_MULTIPLES_OF_THREE_IMAGE_FAMILY.md),
  N6--N8, gives \(\Lambda(xg_j^m)=(-8)^m\kappa_m\ne0\), with
  \(\kappa_m=2^m m!/\prod_{r=0}^m(m+2r+1)\).
- [The hyperbolic proof](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md),
  H6--H9 and its endpoint argument, gives zero odd marked moments and
  nonzero even marked moments for every \(p_{n,\lambda}\).

All positive pure moments of these seeds vanish. These all-order results
are premises for the present interpretation, not consequences of finite
tests. The [rank-one trace theorem](RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md),
Sections 3--5, proves \(\Lambda\circ\sigma=\Lambda\) for every polynomial
automorphism over every characteristic-zero extension. Hence
\[
\Lambda(\sigma(x)\sigma(f)^m)=\Lambda(xf^m).            \tag{11}
\]
The marker \(\sigma(x)\) is fixed independently of \(m\). No claim is made
that the unchanged \(x\) works for every automorphism; the survival orders
are precisely those provided by the selected seed.

All these seeds are absolutely scalar-nonproper-powers in \(A_K\).
Specialize on the polynomial curve \(x=1,y=(1-t^2)/2\). The quartic becomes
\(1-4t^2\), with a simple root at \(1/2\); every \(g_j\) becomes
\(16t^3-16t\), with a simple root at zero. For \(p_{n,1}\), the result is
\(2^{n+1}t(1-U_{n,1}(1)^2t^2)\), with derivative \(2^{n+1}\ne0\) at zero.
A scalar proper power cannot specialize to a nonzero polynomial with a
simple root. Diagonal scaling handles every \(\lambda\ne0\).
The argument survives extensions and is preserved by automorphisms.

In the original pairs, put
\[
R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad
W=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2.
\]
The balanced subring is \(K[R,Z,W,T]/(T^2+2ZW-R^2)\). Restricting \(R=1\)
identifies its degree-\(D\) piece with \(A_{\le D}\); denote the unique
corresponding balanced homogenization by \(H_D(f)\).
[The degree and power-descent note](HOPF_LND_ORBIT_DEGREE_PROPAGATION_AND_POWER_DESCENT.md),
Sections 4--5, proves
\[
R\nmid H_D(f)\quad\Longleftrightarrow\quad D=\deg_A f
\]
and proves that the least lift is a scalar proper power in the **full**
original polynomial ring if and only if \(f\) is one in \(A_K\).
Both grading-extreme arguments are needed for that statement.
It holds after every coefficient extension. Thus the least lifts of
the three orbit families are radial-primitive and absolutely
scalar-nonproper-powers.

For \(D=\deg_A\sigma(f)\), \(e=\deg_A\sigma(x)\), set
\(F=H_D(\sigma(f))\), \(Q=H_e(\sigma(x))\). Exact contraction gives
\[
\mathcal E_2(F^m)=0,\qquad
\mathcal E_2(QF^m)=(mD+e+1)!\Lambda(xf^m).              \tag{12}
\]
The factorial uses the degrees after transport. This realizes the degree
laws as balanced degrees of actual primitive Image failures. Inequivalence
here means scalar multiplication and automorphisms of the Hopf quadric,
not arbitrary transformations of the original four-variable polynomial ring.

## 6. A separate consequence: the coarse balanced Image boundary

Combining retained results gives the following exact **existence** boundary.

| Number \(n\) of contraction pairs | Balanced degrees admitting a counterexample |
|---|---|
| \(n=1\) | None |
| \(n=2\) | Exactly \(D\ge4\) |
| \(n\ge3\) | Exactly \(D\ge2\) |

Radial primitivity is not imposed in the middle row. For \(n=1\), a nonzero
balanced form is \(c(\xi z)^D\), with first contraction \(cD!\ne0\).
In every number of pairs, pure balanced-degree-one forms are safe for all
original polynomial markers by
[LinearSampledVanishingNilpotenceAndSafety](PRIME_BALANCED_IMAGE_FAILURES_ARE_ABSOLUTELY_NONPOWER.md).

The [complete two-pair cubic theorem](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md),
its opening statement and Section "Complete classification and the original
Image conclusion," supplies safety through balanced degree three. That
section explicitly treats unbalanced original polynomial markers,
characteristic-zero descent and the lower degrees. Its geometric
classification is an imported input for this corollary only.
The quartic gives failure in degree four; the
[direct higher-pair family](ALL_BALANCED_DEGREES_IMAGE_WITNESSES_AT_LEAST_THREE_PAIRS.md)
gives primitive nonpower failures for every \(n\ge3,D\ge2\).

For the remaining two-pair entries, use the all-order radial identity
\[
\mathcal E_n(R_n^kB)=(n+d)_k\mathcal E_n(B),\qquad
R_n=\sum_i\xi_i z_i,                                  \tag{13}
\]
where \(B\) is balanced of degree \(d\) and \((u)_k\) is the rising
factorial. For a matching monomial, the one-step sum is
\(\sum_i(\alpha_i+1)\alpha!=(n+d)\alpha!\); for a nonmatching balanced
monomial both sides vanish. Iteration proves (13).
Replace the quartic witness \(F_4\) by \(R^{D-4}F_4\).
Apply (13) to its powers and their products with the fixed balanced marker
\(Z\). Pure contractions remain zero and marked contractions remain nonzero
at every positive order. Padding is not radial-primitive when \(D>4\).

Primitive degrees five and seven are therefore a finer orbit/existence
problem, not missing entries in the coarse boundary.
This corollary inherits the cubic theorem's existing written-proof and
review scope; it is not a new verification of that theorem.
The degree-zero case is harmless: a pure constant is zero and safe.

## 7. Antecedents, exact controls, and scope

Sections 1--4 are self-contained degree proofs. The earlier
[single-shear note](HOPF_LND_ORBIT_DEGREE_PROPAGATION_AND_POWER_DESCENT.md),
Section 9, gives particular progressions; formula (2) controls every
automorphism without classifying coordinate-degree pairs.

For primary context, Blanc and Dubouloz, *Automorphisms of
\(\mathbb A^1\)-fibered affine surfaces*, Trans. Amer. Math. Soc. **363**
(2011), 5887--5924, [Theorem 5.4.5(7), pp. 5915--5918](https://algebra.dmi.unibas.ch/blanc/articles/autoa1fibrations.pdf),
give an amalgamated-product description for automorphisms of \(uv=P(w)\),
including the split smooth quadric. Here \(u=x,v=-2y,w=t\) gives
\(uv=w^2-1\). The actual theorem and factorization argument were inspected.
It is context rather than a premise of the degree laws and does not state
these seed spectra. No worldwide-novelty claim is made.

The exact checker reconstructs regular representatives and numerator
identities, checks formal shear inverses, quotient degrees, low-order
transported moments, original balanced contractions, simple-root controls
and the radial identity. It actively checks a degree-nine member descending
to degree six and rejects changed shear signs, changed seed terms,
wrong degree formulas and a false simple-root certificate for a square.
The deterministic cases are in the checker; there is no separate input file.

Run from the repository root:

~~~sh
python3 research/foundations/extended-geometry/verify_hopf_automorphism_orbit_spectra.py --repo-root .
~~~

[The source ledger](hopf_automorphism_orbit_spectra_source_ledger.json)
pins the proof, checker, arithmetic helper and precisely scoped mathematical
inputs. [The receipt](hopf_automorphism_orbit_spectra_verification.json)
records the bounded replay. Computation supports the displayed identities;
the written arguments prove the unrestricted automorphism and field
assertions. The all-order Image interpretation has the explicit premises
in Section 5. No external review, formal verification, publication,
independent replay or solution of primitive degrees five and seven is
asserted by this evidence packet.

## 8. All odd starting exponents: degree laws and family classification

The family in [the retained hyperbolic construction](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md),
Section 6, allows every odd starting exponent \(h\ge1\).
Using \(C=C_{n,\lambda}\) and \(U=U_{n,\lambda}\) from (1), write
\[
p_{h,n,\lambda}
 =\frac{t^h C(x)U(x)^{h-1}}x(1-t^2U(x)^2)
 =(tU(x))^{h-1}p_{n,\lambda}.                          \tag{14}
\]
Here \(n\ge1\), \(\lambda\ne0\), and \(p_{1,n,\lambda}=p_{n,\lambda}\).
The exponent \(h\) is the starting power of \(t\); the actual polynomial
degree in \(t\) on the chart \(x\ne0\), or height for \(x\partial_t\), is
\(h+2\). It is not \(h\).

**Theorem HigherOddStartHyperbolicOrbits.** Over every characteristic-zero
field and every coefficient extension, for every polynomial automorphism
\(\sigma\), set \(A=\deg_A\sigma(x)\), \(B=\deg_A\sigma(t)\). Then
\[
\deg_A\sigma(p_{h,n,\lambda})
  =(2h+3)nA+(h+2)B.                                   \tag{15}
\]
Its full automorphism orbit has minimum degree
\[
d_{h,n}=(2h+3)n+h+2.                                  \tag{16}
\]
For the displayed family, two seeds \(p_{h,n,\lambda}\) and
\(p_{h',n',\mu}\) are equivalent under a polynomial automorphism and
nonzero scalar multiplication if and only if
\[
(h,n)=(h',n').                                       \tag{17}
\]
Inequivalence persists after every coefficient-field extension.
All these phases are absolutely scalar-nonproper-powers; their least
balanced lifts are primitive Image failures using the retained all-order
moments and trace invariance.

### 8.1. The degree law

Regularity follows from (6) and the second expression in (14).
By (4)--(5), the transformed numerator \(xp_{h,n,\lambda}\) has degree
\[
hB+(n+1)A+2n(h-1)A+2B+4nA
  =\bigl((2h+3)n+1\bigr)A+(h+2)B.
\]
Subtracting the degree \(A\) of \(\sigma(x)\) proves (15).
Since \(A,B\ge1\), the identity automorphism attains the minimum (16).
The same diagonal scaling as in Section 4 gives
\[
\delta_c(p_{h,n,\lambda})=c^{-1}p_{h,n,c\lambda}.        \tag{18}
\]
Thus equal \((h,n)\) imply equivalence for every nonzero pair of parameters.
Minimum degree alone does not prove the converse in (17).

### 8.2. Geometric zero-component multiplicities

Work over an algebraic closure of the coefficient field. The affine
quadric is smooth: its gradient \((2y,2x,2t)\) has no common zero on
\(t^2+2xy=1\). It is therefore normal. At an irreducible divisor its
local ring is a discrete valuation ring, and the multiplicity of a
function along that divisor is its order in this ring. An automorphism
identifies these local rings and preserves those orders. Multiplying a
function by a nonzero constant changes none of them. Consequently the
multiset of multiplicities of the geometric irreducible components of
its zero divisor is a scalar-and-automorphism invariant.

On the chart \(x\ne0\), the coordinate ring is the unique factorization
domain \(\overline K[x,x^{-1},t]\). Put \(u=1+\lambda x\); then
\[
U(x)=\frac{u}{2n}\sum_{j=0}^{2n-1}u^j
     =\frac{u(u^{2n}-1)}{2n(u-1)}.                     \tag{19}
\]
Its \(2n\) roots are simple: \(u=0\), and the \(2n-1\) nonidentity
\(2n\)-th roots of unity. None has \(x=0\), since \(U(0)=1\).
The numerator in (14), with \(x\) a unit on this chart, factors as
\[
t^h\,u^{n+1}\,U(x)^{h-1}\,(1-tU(x))(1+tU(x)).
\]
The factors \(1-tU\) and \(1+tU\) are irreducible there: each is primitive
and linear in \(t\) over \(\overline K[x,x^{-1}]\), with constant
coefficient one, so Gauss's lemma applies. They are distinct, mutually
coprime, and have no common factor with \(t\) or \(U\). Their closures
are distinct irreducible divisors of the quadric, each with multiplicity
one. The component \(t=0\) has multiplicity \(h\), the component \(u=0\)
has multiplicity \(n+h\), and each other root of \(U\) has multiplicity
\(h-1\), provided \(h>1\).

No component has been missed on the omitted chart boundary. The two
irreducible divisors of \(x=0\) have \(t=1\) and \(t=-1\), respectively.
In the regular representative
\[
p_{h,n,\lambda}
 =t^h C U^{h-1}\left(2y-t^2\frac{U^2-1}{x}\right),
\]
restriction to these divisors is, respectively,
\[
2\bigl(y-U'(0)\bigr),\qquad
-2\bigl(y-U'(0)\bigr).
\]
Both are nonzero polynomials in \(y\), since \(h\) is odd and \(C(0)=U(0)=1\).
Thus neither boundary divisor belongs to the zero divisor.

The complete geometric multiplicity multiset is therefore
\[
\mathcal M_{h,n}
 =\{\,n+h,\ h,
       \underbrace{h-1,\ldots,h-1}_{2n-1\text{ entries}},\ 1,\ 1\,\},
 \quad\text{with zero entries discarded}.             \tag{20}
\]
For \(h=1\), this is \(\{n+1,1,1,1\}\).
For every allowed \(h,n\), the largest multiplicity is \(n+h\), and
the second-largest **distinct** multiplicity is \(h\). These recover
\(h\), then \(n\). This proves the necessary direction in (17).
Any automorphism over an extension further extends to its algebraic
closure, so the same invariant proves the assertion over all extensions.

The multiplicity-one components also prove absolute nonpower directly:
every component of the zero divisor of a scalar \(r\)-th power has
multiplicity divisible by \(r\). This argument applies for \(h>1\),
where the simple root at \(t=0\) used for the cubic-height seeds no
longer exists.

### 8.3. What this adds, and its boundary

The degree-\(23\) seeds \(p_{1,4,1}\) and \(p_{3,2,1}\) have the same
minimum orbit degree but different multiplicities:
\[
\mathcal M_{1,4}=\{5,1,1,1\},\qquad
\mathcal M_{3,2}=\{5,3,2,2,2,1,1\}.
\]
They are therefore inequivalent even after coefficient extension.
This is information not supplied by orbit minimum degrees.

**Corollary HyperbolicFamilyOrbitCount.** For every integer \(D\ge2\),
let \(c_{\rm hyp}(D)\) count the scalar-and-automorphism orbits represented
by the displayed family whose minimum degree is \(D\). Put
\(\tau_1(m)=\#\{v>0:v\mid m,\ v\equiv1\pmod4\}\). Then, over every
characteristic-zero field and after every coefficient extension,
\[
c_{\rm hyp}(D)=\tau_1(2D-1)-1-\mathbf1_{D\ {\rm odd}}.   \tag{21}
\]

**Proof.** Equations (16)--(17) identify these orbits with parameter pairs
\((h,n)\) satisfying
\[
2D-1=(2n+1)(2h+3),\qquad h\ge1\text{ odd},\quad n\ge1.
\]
They are in bijection with divisors \(v\mid 2D-1\) such that
\(v\ge5\), \(v\equiv1\pmod4\), and \(v<2D-1\), through
\[
h=(v-3)/2,\qquad n=\bigl((2D-1)/v-1\bigr)/2.
\]
The quotient \((2D-1)/v\) is odd and greater than one, hence at least
three. These are exactly the required parameter conditions.
Among the divisors counted by \(\tau_1\), discard \(v=1\), and also
\(v=2D-1\) when \(2D-1\equiv1\pmod4\), equivalently when \(D\) is odd.
Since \(D\ge2\), these excluded divisors are distinct whenever both occur.
This proves (21). \(\square\)

In particular there are minimum degrees supporting arbitrarily large
finite collections of distinct orbits from this family.
For any integer \(r\ge2\), set
\[
d_r=\frac{5^r+1}{2},\qquad
h_j=\frac{5^j-3}{2},\qquad
n_j=\frac{5^{r-j}-1}{2}\quad(1\le j<r).
\]
Each \(h_j\) is a positive odd integer and each \(n_j\ge1\).
The pairs are distinct and
\((2n_j+1)(2h_j+3)=5^r=2d_r-1\).
All \(r+1\) divisors of \(5^r\) are congruent to one modulo four.
Formula (21) therefore gives **exactly** \(r-1\) displayed-family orbits
with minimum degree \(d_r\), represented by these pairs.
This exact count places no upper bound on other orbits at that degree.

For the Image interpretation, the precise additional all-order input is
HyperbolicOddHeightImage, H18--H22 and Section 6 of the retained
hyperbolic note. Its literal lift is
\((2n)^{h+1}H_{d_{h,n}}(p_{h,n,\lambda})\).
Every positive pure moment is zero; the marker \(x\) has zero odd moments
and nonzero even moments. Equations (11)--(12) transport those properties
and give primitive, absolutely nonpower least lifts in every orbit here.
All the minima (16) are at least eight, so no member of any of these
higher-start orbits has primitive degree five or seven.

Statement (17) classifies equivalence **among this displayed family**.
It does not classify all functions, pure phases, unsafe phases,
automorphisms, or attainable coordinate-degree pairs \((A,B)\).
In particular (15)--(16) do not assert a complete set of degrees in
each higher-start orbit. The bounded checker tests representative degree
laws, factorization and boundary identities, multiplicity recovery and
the equal-minimum inequivalent pair, and bounded independent divisor and
parameter enumerations for (21). The all-parameter classification
is the written argument above, not an extrapolation from those tests.
