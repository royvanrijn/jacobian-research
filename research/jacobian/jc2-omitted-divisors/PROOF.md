# Reconstruction divisors and pointwise missing lengths

14 November 2026. Work over C. The finite-normalization statements below
apply to dominant quasi-finite polynomial maps; flatness of the plane
normalization is used only in the plane assertions. The plane boundary
structure and its purity consequence are the explicit imported results of
[PF2BATT1](../jc2-boundary-attachment/PROOF.md), with its source pins and
quantifiers unchanged. The two-class r7 application retains every hypothesis
of PF2R7BND1. The parameter n in the marked-root control is unrelated to r=7.

## 1. A complete omission test from reconstruction

Let F:U=A^d -> Y=A^d be dominant and quasi-finite. Write A=C[Y],
K=C(U), B=the integral closure of A in K, N=Spec B and pi:N -> Y.
Finiteness of normalization and Zariski's Main Theorem identify U with an
open of N, as in the [singleton-zero proof](../../hessian/hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length).
Regard the original source coordinates x_i as rational functions on N.
Then

\[
 U=\{q\in N:x_i\in\mathcal O_{N,q}\text{ for every }i\}.       \tag{1}
\]

Indeed on this open the x_i define a morphism h to A^d. The identities
F(x)=pi hold as rational identities and hence as regular identities.
On a neighborhood W of q, the maps jh and the inclusion W -> N agree on
their dense intersection with U. Since N is separated they agree on W.
Thus q=j(h(q)) belongs to U. The other implication is immediate.

Normality now makes the test divisorial, using [Stacks, Lemma 031T](https://stacks.math.columbia.edu/tag/031T).
For every prime divisor D of N,

\[
 D\subset N\setminus U\quad\Longleftrightarrow\quad
           \min_i v_D(x_i)<0.                               \tag{2}
\]

Here containment of a prime means that its generic point is omitted;
closedness of the complement then omits the whole divisor. A normal local
domain is the intersection of its height-one valuation rings. Consequently
any point at which a coordinate fails to be regular lies on a negative
valuation divisor. Thus (2) lists the entire complement, including all its
special points; there is no separate isolated omitted set.

There is also an ideal calculation that does not depend on a chosen
denominator expression. Form I_i=(B:x_i)={b in B:b*x_i in B}. Each x_i is
regular on the open D(I_i), so

\[
 N\setminus U=V(\prod_i I_i).                              \tag{3}
\]

Compute the height-one minimal primes of this ideal, with valuations in B.
Cancellation of a displayed fraction must be performed in B; raw resolvent
discriminants or denominator factors in a nonnormal order cannot certify
this calculation. To certify completeness one needs B itself, the inverse
identities and equality (3), not just a list of ramification primes.

For each target divisor H the full finite-cover DVR identity is

\[
 \sum_{D\mid H} e_D f_D=[K:C(Y)].                         \tag{4}
\]

This sum includes retained primes. It is a useful generic check but has
different content from (3). The ramification coefficient e_D−1, the omission
flag, the image H, and f_D must occupy separate ledger fields.

## 2. Exceptional fibres require local Artin lengths

For a finite flat pi of degree d and a closed target point y, define

\[
 \ell_q(y)=\dim_C\bigl(\mathcal O_{N,q}/\mathfrak m_y
                                  \mathcal O_{N,q}\bigr).
\]

For a Keller source every point retained in U is etale and has length one.
The finite fibre splits into its local Artin algebras, giving

\[
 d=\#F^{-1}(y)+\sum_{q\in(N\setminus U)\cap N_y}\ell_q(y).
                                                               \tag{5}
\]

Count each omitted point once, even if several divisors meet there. A
boundary curve's own fibre length is not in general the ambient local
length ell_q. In particular, replacing ell_q by generic e_D*f_D at an
exceptional value has no justification. Without flatness replace d in (5)
by length(N_y). The normal finite surface over smooth A2 is flat: it is
Cohen–Macaulay and miracle flatness applies. Hence (5) holds for plane
Keller maps at every finite y.

## 3. An explicit marked-root control

Fonollosa's [Theorem 1.3 and Sections 5.2–5.4](https://arxiv.org/html/2609.33491v1#S5)
give the following model. Put m=n−2 and
P(T)=qT^n+rT^(n−1)+sum(lambda_j*T^j)+T−b. The finite normalization has charts

\[
 P(t)=0,\quad a^m=P'(t);\qquad
 Q(u)=u^nP(1/u)=0,\quad v^m=-Q'(u),
 \quad u=t^{-1},\ v=ua.                                   \tag{6}
\]

In the infinity chart the reconstruction is x=1/v, y=u−v. At t=0 the
branches have a in mu_m; reconstruction y=(1−a)/t retains a=1 only.
The omitted primes are Ecrit=(v=0) and E_zeta=(t=0,a=zeta), zeta!=1.
Their respective (e,f) are (2m,1) over the binary discriminant and (1,1)
over b=0. They are pairwise disjoint. The simple root at infinity is
retained. These source facts are used as the control, not as a plane result.

Here is a direct fibre-algebra derivation of all exceptional lengths from
(6). If r is a finite root of multiplicity k, localize C[t]/(P) at t=r.
It has length k. Adjoining a through a monic equation of degree m makes
the marking algebra free of rank m over it. Thus its total length is mk.
For k>=2 the only marking in its reduced support is a=0, so all mk units
are concentrated at one omitted point. For k=1 the derivative is nonzero
and the m markings are distinct, each of length one. If r=0, P'(0)=1
forces k=1: one marking remains affine and m−1 are omitted. The same
argument on Q handles infinity, including degree drops of P. Therefore,
writing s for the number of simple roots of the degree-n binary form,

\[
 \begin{split}
 \operatorname{length}(N_P)&=nm,\\
 L_{\rm omitted}(P)&=m\sum_{r:k_r\ge2}k_r
                    +(m-1)\mathbf1_{P(0)=0},\\
 \#U_P&=ms-(m-1)\mathbf1_{P(0)=0}.                         \tag{7}
 \end{split}
\]

The sums include infinity. This also distinguishes distinct points from
scheme lengths: a k-fold root gives one omitted point of length mk, not
mk reduced points. The algebra checker uses Groebner standard monomials,
not a root counter, for explicit finite and infinity charts.

At the zero target P(T)=T, the finite zero root is simple and infinity
has multiplicity n−1. Equation (7) therefore gives the singleton-zero control

\[
 nm=1+\underbrace{m(n-1)}_{\text{one ramified omitted point}}
          +\underbrace{(m-1)}_{\text{unramified omitted points}}. \tag{7a}
\]

For n=4 this is 8=1+6+1; for n=5 it is 15=1+12+2. The ramification-only
account misses a positive amount even at the singleton value. This applies
to the full source A^n and, using the cited coefficient descent, to the
three-dimensional zero-lower-coefficient slice. It supplies no plane
counterexample.

The n=4 unramified loss is also visible in the defining polynomial, without
counting its roots numerically. With p=1+xy and source coordinates
(x,y,z3,z4), put

\[
 \Phi_4(T)=(1-yT)^3(pT-x)
 +(pT-x)^2[4y(yT+1/2)+z3(pT+x/2)+z4T^2].
\]

For b!=0, both of the following source points have exactly Phi4(T)=T−b:

\[
 (b,0,0,0),\qquad (-b,2/b,12/b^2,-8/b^3).                \tag{7b}
\]

The second has t=b, a=−1, x=t/a and y=(1−a)/t. As b tends to zero it
escapes along the unramified omitted branch, whereas the first tends to
the source origin. Direct substitution checks (7b); the local algebras
give the persistent length-six infinity point and the extra length-one
point at b=0. Thus the missing unramified contribution is independently
exhibited by a rational source path.

## 4. A pointwise completion of the conditional r7 ledger

Assume precisely the PF2R7BND1 setting: geometric degree six and the entire
reduced ramification divisor one prime R of generic (e,f)=(2,1). Its image
C is one of the two classified nonimmersed curves. PF2BATT1 says that all
other omitted primes D_j are disjoint from R and from each other, etale
at every finite point, and map by the normalization of an immersed
polynomial image C_j. Each has e=f=1 and C_j!=C; different j may share an
image. PF2R7BND1 gives t=3−epsilon such primes, with epsilon in {0,1}.

Let b_j(y) be the number of points of C_j^nu=A1 over y; equivalently it is
the degree of gcd(p_j(T)−y_1,q_j(T)−y_2) for a birational polynomial
normalization of C_j. The gcd is squarefree because the parametrization
is immersive. Put B(y)=sum_j b_j(y), counting each PRIME separately even
when images coincide. All these points have ambient length one.

The ramified contribution L_R and all remaining fibre lengths are:

| Target stratum | Points on R | Their local ambient lengths | L_R | Source length |
|---|---:|---|---:|---:|
| y outside C | 0 | none | 0 | 6−B(y) |
| smooth point of C | 1 | 2 | 2 | 4−B(y) |
| A1 node or the A5 value | 2 | 2, 2 | 4 | 2−B(y) |
| A2 cusp, epsilon=0 | 1 | 2 | 2 | 4−B(y) |
| A2 cusp, epsilon=1 | 1 | 3 | 3 | 3−B(y) |

The local orbit assertions are exactly those proved with their global
separation hypotheses in PF2R7BND1. The table then follows from (5), not
from a generic e*f substitution. It includes all intersections with C_j
and their self-intersections through b_j(y).

In particular every two-branch singular value of C permits B(y)<=2.
No collection of extra omitted primes may contribute three branches there.
The cusp permits B(y)<=4−epsilon. At a point outside C the bound is six.
These are necessary conditions; zero source points are allowed and do not
by themselves imply a contradiction. If a point additionally has a
singleton affine fibre, then

\[
 B(y)=5-L_R(y):\quad 5\text{ off }C,\quad 3\text{ on }C_{sm},
 \quad 1\text{ at a two-branch value},\quad
 3-\epsilon\text{ at the cusp}.                            \tag{8}
\]

For example, if all extra images pass through the cusp with one branch
each, their contribution t=3−epsilon is numerically compatible with a
singleton source fibre there. The count is not a realization theorem.
The images, actual intersection parameters, pole valuations and infinity
attachments are still unknown. The table is an exact conditional ledger,
not a certified enumeration of the global boundary. It does not cover
the remaining immersed normal-r7 coefficient family.

## 5. The singleton-zero ledger

For a plane Keller G with (G1,G2)=(x,y), let d be its generic degree.
The retained zero fibre has length one, hence

\[
 \sum_{q\in(N\setminus U)\cap N_0}\ell_q(0)=d-1.           \tag{9}
\]

The plane boundary theorem makes its primes disjoint. We may split (9)
without double counting into

\[
 \sum_{D:e_D>1}\ \sum_{q\in D\cap N_0}\ell_q(0)
       +\sum_{D:e_D=1} b_D(0)=d-1.                        \tag{10}
\]

In the second sum, e=f=1, the image is immersed and every local length is
one by PF2BATT1. The first sum requires actual local fibre algebras; e*f
alone is insufficient. Empty ramified contribution over zero does not
make the second sum vanish. Unramified image curves can self-intersect at
zero, producing more than one length-one point on a single omitted prime.

Since pi is finite, pi(N\U) is closed. Thus G is proper over a
neighborhood of zero exactly when this omitted fibre is empty. In the
singleton setting that forces d=1 and invertibility, by HC4MYOR2. To close
properness, one must rule out every omitted point of BOTH types in (10).
No such exclusion is supplied by a ramification-only calculation.

## 6. Completeness audit and remaining inputs

The [machine-readable ledger](ledger.json) has three scopes: the explicit
marked-root control with complete omission and fibre formulas; the
conditional two-class normal-r7 table; and the unrestricted singleton-zero
length obligation. A complete application must supply the finite normal
algebra, the reconstruction ideal (3), all resulting primes and images,
their e and f, and local Artin algebras at every exceptional image value.
Checking the DVR degree sums and (5) is necessary but does not prove that
the candidate normal algebra or omission list was complete.

Unknown images and exceptional local algebras remain null, never empty
lists interpreted as a proof of absence. The checker deliberately rejects
completeness for those unresolved applications. It also tests an immersed
nodal image, where one omitted prime contributes two etale fibre points.
This prevents confusing a prime count with an exceptional missing length.

No new large elimination is required or authorized by these ledgers.
The next mathematical input is global reconstruction and attachment data
for the extra r7 images, or a theorem eliminating the entire omitted fibre
at zero. The source pins and written imports are not an independent
reproof or formal verification of the cited boundary theorems.
