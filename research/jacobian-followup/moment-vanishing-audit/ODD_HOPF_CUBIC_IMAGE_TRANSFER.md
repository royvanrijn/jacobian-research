# The entire odd-harmonic two-pair (3,3) Image slice

16 November 2026. Internal corollary of the retained ternary cubic
theorems. The polar and multiplier transfers below are exact; their
Gaussian input retains its outstanding independent specialist review.

## Statement

**SIC2B33ODD20261116.** Let P be a complex biform of dual degree three
in (W,V) and coordinate degree three in (Z,Y). Suppose its Hopf-sphere
restriction is

\[
 F=H_3+H_1,
\]

where H_j is a degree-j spherical harmonic. If every positive pure
contraction vanishes, then, for every polynomial multiplier Q,

\[
 \boxed{\mathcal E_2(P^mQ)=0\quad(m>\deg_{Z,Y}Q).}       \tag{1}
\]

Here coordinate degree is the largest such degree of any bihomogeneous
component of Q. In particular, the proposed bound m>3 deg_(Z,Y) Q holds.
The coefficient one in (1) is sharp, including for phases with both
harmonic summands nonzero. This closes the entire ten-dimensional
Sym^6 + Sym^2 slice, including every mixture and degeneration.

No assertion about nonreduced coefficient algebras is made. The full
mixed Sym^6 + Sym^4 + Sym^2 problem remains open.

## Hopf normalization and homogeneous cubic lift

Write

\[
 \begin{split}
 s&=WZ+VY,\\
 h_1&=WY+VZ,\qquad h_2=(WY-VZ)/i,\qquad h_3=WZ-VY.
 \end{split}
\]

Then h_1^2+h_2^2+h_3^2=s^2. On the unit complex three-sphere,
W=conjugate(Z), V=conjugate(Y), and s=1. Its normalized measure pushes
forward by h to normalized surface measure sigma on S^2. Thus for any
balanced biform A of bidegree (D,D),

\[
 \mathcal E_2(A)=(D+1)!\int_{S^3} A\,d\sigma_3.       \tag{2}
\]

For example, W^r V^(D-r) Z^r Y^(D-r) has contraction r!(D-r)!
and sphere integral r!(D-r)!/(D+1)!; unmatched monomials vanish by
the circle phases. This fixes the normalization without a convergence
assumption. Balanced biforms descend through the Hopf circle; a form
of bidegree (b,b) gives spherical harmonics of degrees at most b.

Use the homogeneous harmonic polynomial representatives of H_3 and H_1
and set

\[
 G(x)=H_3(x)+(x\cdot x)H_1(x).                         \tag{3}
\]

It is a homogeneous ternary cubic and G|_(S^2)=F. In these coordinates
P=G(h)=H_3(h)+s^2 H_1(h). The restriction determines P: real radial
homogeneity recovers a balanced form from its values on S^3, and the
real locus (W,V)=(conjugate(Z),conjugate(Y)) determines its polynomial
coefficients. The harmonic dimensions are seven and three.

For an independent standard real Gaussian X=RU in R^3, R and U are
independent and U has law sigma. Consequently

\[
 \begin{split}
 \mathbb E[G(X)^m]
   &=c_{3m}\int_{S^2} F^m\,d\sigma,\\
 \mathbb E[G(X)^m K_\ell(X)]
   &=c_{3m+\ell}\int_{S^2} F^m K_\ell\,d\sigma,\\
 c_k&=\mathbb E[R^k]
     =2^{k/2}\frac{\Gamma((k+3)/2)}{\Gamma(3/2)}>0,
 \end{split}                                                   \tag{4}
\]

where K_ell is any homogeneous harmonic polynomial of degree ell.
All integrands are polynomials and all displayed moments are finite.
Equations (2) and (4) transfer spherical pure nullity to Gaussian pure
nullity, with no loss of indices.

The [full ternary cubic theorem](TERNARY_CUBIC_GMC_CLOSURE.md)
G3CUBICALL20261115 gives marked Gaussian vanishing for m>3 ell.
This alone proves the proposed spherical cutoff m>3b for every marker
of spherical degree at most b. There is a stronger retained input:
[G3HOMOGENEOUS20261114](NODAL_NEWTON_EXCLUSION.md) gives

\[
 \mathbb E[G(X)^m K(X)]=0\quad(m>\deg K)               \tag{5}
\]

for every pure-null homogeneous ternary cubic. It also follows from
the [odd cubic classification](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md).
Apply (5) to each K_ell and divide by the nonzero radial factor in (4).
Harmonic decomposition now gives

\[
 \int_{S^2} F^m q\,d\sigma=0\quad(m>b)                \tag{6}
\]

for every spherical polynomial q of degree at most b. Decomposition
is essential: an inhomogeneous marker need not have one radial factor.

## From balanced markers to the full polynomial output

Let Q have dual degree a and coordinate degree b. If a>b, contraction
is zero by total degree for every m. If a<=b, the output
E_2(P^m Q) is homogeneous of coordinate degree b-a. Write
z=(Z,Y), w=(W,V). For every multiindex gamma of size b-a,

\[
 [z^\gamma]\mathcal E_2(P^mQ)
   =\frac{\mathcal E_2(w^\gamma P^mQ)}{\gamma!}.        \tag{7}
\]

Indeed E_2(w^gamma A)=partial_z^gamma E_2(A); on an output homogeneous
of degree |gamma| this derivative is the constant gamma! times its
gamma coefficient. The balanced multiplier w^gamma Q has bidegree
(b,b), independent of m. Its Hopf restriction q has degree at most b.
Therefore (2) gives

\[
 \mathcal E_2(w^\gamma P^mQ)
    =(3m+b+1)!\int_{S^2}F^m q\,d\sigma=0\quad(m>b).
\]

This proves every output coefficient zero. Bihomogeneous decomposition
of Q proves (1). This is the same balancing argument already recorded in
Section 3 of the [Rodrigues-survivor note](../../archive/non-elliptic/extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RODRIGUES_SURVIVOR.md).

## Sharpness within mixed H_3 + H_1

Let A=x+iy and take G=A z^2. Its harmonic decomposition is

\[
 H_3=A z^2-\tfrac15(x^2+y^2+z^2)A,\qquad H_1=\tfrac15 A.
\]

Both components are nonzero. The biform P=2WY(WZ-VY)^2 has this Hopf
restriction. Its pure moments vanish by rotation about the third axis.
For b>=1 take Q_b=(2VZ)^b, a balanced marker of coordinate degree b.
At m=b, (2) gives the strictly positive value

\[
 \mathcal E_2(P^b Q_b)
  =(4b+1)!\,\frac12\int_{-1}^1(1-t^2)^b t^{2b}\,dt
  =\frac{(4b+1)!}{2}B(b+\tfrac12,b+1)>0.              \tag{8}
\]

Thus the strict inequality and coefficient one in (1) cannot be improved
uniformly, even inside the mixed slice.

## Remaining structural boundary

The first pure moment removes H_0. After this closure, an unresolved
phase must have H_2 nonzero. Pure H_2 is already safe; the new structural
question is its coupling to H_3+H_1, after retaining the other known
closed subfamilies.

For the full phase H_3+H_2+H_1, a homogeneous lift of degree three
would be H_3(x)+sqrt(x.x) H_2(x)+(x.x)H_1(x). The middle term is not
polynomial. More generally, ordinary radial homogenization uses
(x.x)^((d-j)/2) on a degree-j harmonic, so one homogeneous polynomial
degree can accommodate only one parity of j. This identifies a failure
of this transfer method, not an obstruction to some other proof.
Another arbitrary support-size threshold does not address that parity
coupling. No dense mixed-parity classification is claimed here.

## Exact checks and assurance

Run from the repository root:

```sh
.venv/bin/python research/jacobian-followup/moment-vanishing-audit/verify_odd_hopf_cubic_transfer.py
```

The [checker](verify_odd_hopf_cubic_transfer.py) verifies the Hopf identity,
the ten-dimensional cubic substitution, balanced sphere normalization on
a polynomial basis, Gaussian polar moments, coefficient extraction and
the sharp mixed-summand example. Its [receipt](odd_hopf_cubic_transfer_verification.json)
records these finite exact algebra checks. The all-order transfer is the
written argument above; neither finite checks nor this corollary supply
independent verification of the retained Gaussian classification.
