# A universal p-power obstruction and the separable polynomial classification

15 November 2026. **QPC1.** In the quadratic coefficient chart
[a,v]=h/v, every finite scalar operator G(a)+mu*v with nonzero mu
and nonlinear polynomial G has no finite Weyl mate, of any order,
over any finite algebraic coefficient cover in characteristic zero.

The proof isolates one coefficient in the p-th power after good reduction.
It is independent of the lower coefficients of G and of the proposed
mate's order and poles. With [QVL1](ALGEBRAIC_INFINITY_LIMITS.md) and
[QZM1](QUANTUM_ZERO_MODE_TRANSPORT.md), it completes the finite-mate
classification of the entire pure polynomial form v*H(a)+G(a).

This closes a family, not all Long representatives. Genuine mixed
v,a coefficients, additional negative-v terms, other infinity faces,
nonpositional first images and general projective carriers remain open.

## Statement and linear trace descent

Let K be any characteristic-zero field, h in K nonzero, and
delta=v^-1*partial_v on F=K(v). For a finite algebraic extension E/F,
extend delta uniquely and form the finite Ore operator rings
\[
 D(F)=F[a;h\delta]\subset D(E)=E[a;h\delta],\qquad
 S=G(a)+\mu v,\quad
 G=\sum_{j=0}^n g_j a^j,\quad n\ge2,\quad g_n\mu\ne0.          \tag{1}
\]
There is no finite T in D(E) with [S,T]=kappa for nonzero kappa in K.
Scalars may belong to a Laurent-series field such as K0((h)).
The leading coefficient g_n need not be an h-valuation unit.

Normalized field trace tau=Tr_(E/F)/[E:F] commutes with delta.
Coefficientwise extension
\[
 \tau\left(\sum_i t_i a^i\right)=\sum_i\tau(t_i)a^i,\qquad
 \tau(P T Q)=P\tau(T)Q\quad(P,Q\in D(F))                     \tag{2}
\]
is a D(F)-bimodule retraction. Left multiplication by a uses
delta*tau=tau*delta; right multiplication uses the same identity
in the Ore Leibniz formula. Multiplication by base coefficients
commutes with tau. These generator checks prove the assertion.
Thus [S,T]=kappa implies [S,tau(T)]=kappa, reducing existence to F.
There is no bound on the extension degree or the mate's order.

This linear descent also works for matrix mates of the fixed scalar
S*I_m: take the normalized matrix trace first, then the coefficient
field trace. It does not descend a simultaneously unknown nonlinear
tuple or an idempotent. For example the quadratic-cover projector
e=(1/2)*[[1,w],[w^-1,1]], w^2=v, has e^2=e but coefficient trace
tau(e)=I_2/2, which is not idempotent.

## The finite-characteristic necessity

Let A,B be N-by-N matrices over a field of characteristic p and
[A,B]=kappa*I_N, with nonzero kappa. Over dual numbers epsilon^2=0,
\[
 (I-\epsilon B)A(I+\epsilon B)=A+\epsilon\kappa I,\qquad
                       \chi_A'(z)=0.                       \tag{3}
\]
Indeed conjugacy preserves the characteristic polynomial, whereas
translation replaces chi_A(z) by chi_A(z-epsilon*kappa).
Consequently chi_A is a polynomial in z^p. When N=p it is
z^p+c, so Cayley--Hamilton makes A^p scalar.

Over a perfect field k of odd characteristic p, D=k(v)[a;h*delta]
has center and central degree
\[
 Z(D)=k(v^p)[a^p],\qquad
 \dim_{\operatorname{Frac}Z(D)}
   \left(D\otimes_Z\operatorname{Frac}Z(D)\right)=p^2.        \tag{4}
\]
For completeness put b=v*a. Then [b,v]=h and a=v^-1*b.
The ring D is the localization of the rank-one Weyl algebra
k[v,b] at all nonzero polynomials in v. Such polynomials are
invertible after central localization since f(v)^p is central.
The usual Weyl center is k[v^p,b^p], and the generic algebra is
central simple of degree p. Also delta^p=0 on k(v): each Laurent
monomial acquires p successive factors j,j-2,...,j-2(p-1),
one of which is zero. The p-th power of a derivation is a
derivation, so vanishing on v and k gives vanishing on k(v).
Thus a^p is central. The p^2 monomials v^i*a^j, 0<=i,j<p,
give the displayed center and dimension by PBW and commutation
with v and a. Generic splitting is faithful into p-by-p matrices.

Hence a scalar mate in D would make S^p central by (3).
This is a necessary condition, with no prescribed center map.
The standard Weyl-center and splitting framework is described in
[Lauritzen--Thomsen, Section1.2](https://arxiv.org/html/1508.02499)
and their [2026 paper, Section3.1](https://arxiv.org/html/2601.23110).
Neither a lift modulo p^2 nor their Poisson-center lifting theorem
is used in the argument here.

## One coefficient in every degree

Suppose p>2n+1. In normal order, the exact coefficient is
\[
 [v^{-2p}a^{J}]\,S^p
  =\frac{n^2-1}{24n}\,g_n^{\,p-2}\mu^2h^{\,p+1},
 \qquad J=(n-1)p-(2n+1).                                   \tag{5}
\]
It is nonzero: p>2n+1 excludes every prime divisor of
24*n*(n^2-1). Also J>0 and J is not divisible by p.
Therefore S^p is not central.

Here is a direct computation valid for every n and prime in that range.
A word with ell copies of v has p-ell polynomial factors.
If normal ordering uses c derivative contractions, then
\[
 \deg_v=\ell-2c,\qquad
 \operatorname{ord}_a\le n(p-\ell)-c.                       \tag{6}
\]
To reach v^-2p one needs ell even and c=p+ell/2.
With ell=0 no coefficients are differentiated and no such term occurs.
For ell>=4 the upper a order is strictly below J.
For ell=2 equality forces all p-2 polynomial factors to contribute
g_n*a^n. Thus lower G terms and words with any other v count
cannot cancel the coefficient in (5).

Put m=(p+1)/2 and c_r=product_(t=0)^(r-1)(1-2t).
The two-v words are a^(ni)*v*a^(nj)*v*a^(nk), i+j+k=p-2.
Their normal-order contribution with r+s=p+1 is
\[
 \binom{ni}{r}\binom{n(i+j)-r}{s}
       \delta^r(v)\delta^s(v),\qquad
 \delta^r(v)=c_r v^{1-2r}.                                  \tag{7}
\]
Since c_r=0 for r>=m+1, only r=s=m survives.
After removing the scalar factor g_n^(p-2)*mu^2*h^(p+1),
the desired coefficient is
\[
 c_m^2\sum_{u=0}^{p-2}
       \left(\sum_{i=0}^{u}\binom{ni}{m}\right)
                          \binom{nu-m}{m}.                 \tag{8}
\]

Work now in F_p[u]. Set
A(u)=binom(nu,m), B(u)=binom(nu-m,m), and let F(u) be the
polynomial summatory function sum_(i=0)^u A(i). It has degree m+1<p;
it exists by polynomial finite summation since m+1 is invertible.
Furthermore F(p-1)=sum_(i in F_p)A(i)=0, because deg A=m<p-1.
Thus the outer sum in (8) extends to all of F_p.
Since deg(F*B)=2m+1=p+2<2p-2, finite-field power sums give
\[
                    \sum_{u\in\mathbb F_p}F(u)B(u)
                         =-[u^{p-1}]F(u)B(u).              \tag{9}
\]

Only the four highest coefficients of each factor matter:
the target degree p-1=2m-2 is three below 2m+1.
Normalize A and B by the leading coefficient L=n^m/m!.
The top four descending coefficients, on reducing m=1/2 in F_p, are
\[
 \begin{split}
 A/L &: \left(1,\frac1{8n},\frac1{128n^2},
                                      -\frac5{1024n^3}\right),\\
 B/L &: \left(1,-\frac1{8n},\frac1{128n^2},
                                       \frac5{1024n^3}\right).
 \end{split}                                               \tag{10}
\]
These formulas follow by elementary symmetric sums of the roots
0,...,m-1 for A and m,...,2m-1 for B, divided by n.
For example the first three symmetric sums for A are
m(m-1)/2, m(m-1)(m-2)(3m-1)/24, and
m^2(m-1)^2(m-2)(m-3)/48.
For B, use the power sums
m(3m-1)/2, m(2m-1)(7m-1)/6, and
m^2(15m^2-14m+3)/4, then Newton's identities.
All displayed denominators are invertible for p>2n+1.

The top four descending coefficients of F/L are
\[
 \left(
 \frac23,\quad
 \frac12+\frac1{4n},\quad
 \frac1{24}+\frac1{16n}-\frac1{64n^2},\quad
 -\frac1{192n}+\frac1{256n^2}+\frac5{1536n^3}
 \right).                                                  \tag{11}
\]
Indeed if A/L has top coefficients A0,A1,A2,A3, finite summation
gives A0/(m+1), A0/2+A1/m,
A0*m/12+A1/2+A2/(m-1), and
A1*(m-1)/12+A2/2+A3/(m-2).
These expressions give (11) on putting m=1/2.

Multiplying the appropriate four coefficient pairs yields
\[
 \frac{[u^{p-1}]F B}{L^2}
                       =\frac{1-n^2}{96n^3}.               \tag{12}
\]
Finally c_m=2^m*m! modulo p, since
c_m=2^m*m!*binom(1/2,m) and 1/2=m in F_p.
Fermat's identity gives 2^(2m)=4 and n^(2m)=n^2. Combining
(8), (9), and (12) gives
\[
 -c_m^2 L^2\,\frac{1-n^2}{96n^3}
        =-4n^2\,\frac{1-n^2}{96n^3}
        =\frac{n^2-1}{24n},                                \tag{13}
\]
which proves (5).

## Return to characteristic zero

Suppose a finite rational mate exists in characteristic zero.
Adjoin to Z the finitely many scalar coefficients in S,T,h,kappa,
and invert g_n,mu,h,kappa and one nonzero coefficient from each
denominator polynomial in v. This gives a nonzero finitely
generated characteristic-zero coefficient ring R.
Clear rational denominators in the finitely many coefficients of
[S,T]=kappa to obtain polynomial identities over R.

There are good reductions to finite fields in arbitrarily large
characteristic. One elementary proof specializes the finite-type
Q-algebra R tensor Q at a maximal ideal. Zariski's lemma makes
the residue a number field. The image of R lies in O_L[1/N]
for some integer N, and the chosen units remain units after
excluding finitely many prime ideals. Reducing at primes of
arbitrarily large rational characteristic preserves the degree
of G, every required nonzero scalar, every denominator polynomial
and the commutator identity. This is the usual finite-type
reduction framework; see
[Lauritzen--Thomsen, Section1.4](https://arxiv.org/html/1508.02499).

Choose a good p>2n+1. The reduced pair would make S^p central
by (3)--(4), contradicting (5). Together with (2), this proves
the theorem for all finite coefficient covers.

This does not specialize h-adic series by evaluating h at 1.
Only finitely many field elements occur as coefficients of the
hypothetical finite operators, even if those elements themselves
are Laurent series. Finite generation and good reduction apply
to that scalar ring.

## Complete pure separable classification

For S=v*H(a)+G(a), H,G in K[a], a finite scalar mate with
commutator h exists over a finite coefficient cover exactly
in the following cases, where it already exists over K(v).

- H=mu is constant, G=lambda*a+c is affine, and S is nonconstant.
  If mu is nonzero, an explicit mate is
\[
 T=-\frac{\lambda a^2}{2\mu^2}
       -\frac{va}{\mu}-\frac{h}{2\mu v},
                       \qquad [S,T]=h.                    \tag{14}
\]
  If mu=0 and lambda is nonzero, use
\[
                         T=\frac{v^2}{2\lambda}.           \tag{15}
\]
- H=lambda*a+mu with lambda nonzero and G constant. Then
\[
                              T=\frac v\lambda.            \tag{16}
\]

Every other case is excluded. Constant S has zero commutator.
For H constant nonzero and G nonlinear use the present theorem.
For H=0 and G nonlinear use QVL1's ordinary infinity lemma,
with coordinate u=v^2/2. For deg H>=2 use QVL1's complete
v-linear bounded-lower theorem. For deg H=1 with G nonconstant
use QZM1. These cases exhaust all polynomial H,G.

In particular no bounded full exact-R scalar Long lift can have
this pure separable polynomial form: its classical G is
4a-3r*a^2, so its generic G remains nonlinear. Every possible
finite momentum degree in H and G is covered, including
h-small highest coefficients. Finite matrix padding with S
unchanged as S*I does not help, by normalized matrix trace.
This is a structural closure of the ansatz, not a requirement
merely to increase a rank or order threshold.

## Matrix boundary, controls and remaining target

For N=m*p, equation (3) gives chi_A(z)=f(z^p).
It does not force A^p scalar when m>1.
The actual [QPH1 phase matrices](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md),
reduced at p=7 on q=1+t, a=partial_t, t^7=0, give
[S,T]=I_14 and chi_S(z)=z^14-1, while S^7 is nonscalar.
Thus applying the scalar p-central criterion indiscriminately
to matrix carriers would reject an existing valid representation.
The characteristic-polynomial derivative test remains necessary
for arbitrary finite matrix fibers. It is not sufficient.

[Exact controls](check_p_power_obstruction.py) check full sparse
Laurent Ore powers in declared degrees and good primes, the
independent two-v word sum, the rational coefficient identity
(12)--(13), scalar p-by-p split fibers, all positive classification
families, linear trace identities, the projector countercontrol,
and the genuine phase-matrix boundary. A small-prime boundary
case has central S^p and shows why the stated prime gate matters.
Finite controls are not a replay of universal all-degree coverage.

The proof gives an obstruction for a fixed operator and, with the
previous results, a complete scalar classification in one entire
form. It does not prove every bounded Long carrier is equivalent
to that form. Additional negative-v or genuinely mixed coefficient
terms, alternative leading faces, deformed first images and
projective/matrix actions with a genuinely nonscalar S remain
possible targets. Any positive model still needs all full relations,
the correct classical graph, global admission and finite-left
generation. A scalar branch computation alone does not supply those
finite-carrier conditions.

No prescribed arithmetic support, p^2 lift, all-representative
nonlifting, finite-left Long carrier, bounded full Long lift,
scalar or matrix counterexample, independent theorem replay,
formal verification, external review, novelty priority or DC2
solution is claimed.
