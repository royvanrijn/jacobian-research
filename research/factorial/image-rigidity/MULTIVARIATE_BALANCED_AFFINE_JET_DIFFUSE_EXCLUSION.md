# Exact exclusion near isolated diffuse contacts of balanced affine templates

Fix an ordinary two-monomial core \(C\) and adjacent orders
\(n,N=n+1\), \(n\ge2\). Let the appended row types be
\[
 \left(\frac{\beta_i}{u}+a_i,\frac{\beta_i}{u}+b_i\right),
 \qquad \beta_i>0,\quad \beta_i,a_i,b_i\in\mathbb Q,\qquad
 d_i=b_i-a_i\in\mathbb Z\setminus\{0\},
\]
with multiplicities \(r c_i+s_i\), where the \(c_i\) are fixed
nonnegative integers, at least one is positive, and the \(s_i\)
are fixed integers. The admissible \(u>0\) and integer \(r\)
must give nonnegative integer endpoints and multiplicities.
Assume \(\sum_i c_i d_i=0\) and put
\[
 A=\sum_i c_i d_i^2/\beta_i>0,\quad
 B=\sum_i c_i d_i/\beta_i,\quad
 h=\sum_i s_i d_i,\quad
 \tau=Aru,\quad \Lambda=N/n .
\]
The slopes \(\beta_i\), offsets, and orientations are arbitrary
under these hypotheses. At bounded \(\tau\) and \(u\to0\),
the ordinary diffuse parameters satisfy
\[
 \sum_i\frac{(r c_i+s_i)d_i^2}{\beta_i/u+\max(a_i,b_i)}
           =\tau+O(u),\qquad
 \sum_i\frac{(r c_i+s_i)d_i}{\beta_i/u+\max(a_i,b_i)}
           =\frac BA\tau+O(u).
\]
The signed gap is the fixed integer \(h\).

Write the raw core coefficients as \(C_{j,k}=[z^k]Q_j^C\).
For \(j=n,N\) and \(0\le k\le j\), put
\[
 f_{j,k}^0=\frac{k(k-j)}{2j}
               +\frac{B}{A}\frac{k}{2}\left(\frac1j-\frac1n\right).
\]
Let \(\tau_*>0\) be a physical common-root strength of the
reference pair with coefficients
\[
 C_{j,k}^*=C_{j,k}\Lambda^{h k\,\boldsymbol 1_{\{j=N\}}}
                         \exp(\tau_*f_{j,k}^0).
\]
The second rational vector \(f_{j,k}^1\) is calculated below
from the exact row factors. Take an integer \(D>0\) clearing
the denominators of both rational vectors. Test whether
\((X,Y)=(1,1)\) is an isolated complex point of the locus
where the two polynomials with coefficients
\[
                 C_{j,k}^* X^{D f_{j,k}^0}Y^{D f_{j,k}^1}
                                                              \tag{1}
\]
have a common factor of degree at least two. This is an
algebraic test on the two-dimensional torus, with algebraic
coefficient constants.

**Theorem.** If the test passes, there are computable positive
\(\epsilon,\delta\) such that the original ordinary moments of
\(C\) with the indicated tail are coprime at \(n,N\) whenever
\(0<u<\epsilon\) and \(|\tau-\tau_*|<\delta\). The test permits
multiple contacting roots and singular contacts. It includes
unmatched orientations, different slopes, and a nonzero limiting
signed drift. The bounded additive multiplicities \(s_i\) affect
the reference test through \(h\) alone. Once their actual values
are specified they enter the computable neighborhood bounds.

## Exact physical normalization and its first jet

In sums with a signed integer upper index \(q\), use the convention
\(\sum_{\ell=1}^q=-\sum_{\ell=q+1}^0\) when \(q<0\).
Define analytic germs with rational Taylor coefficients by
\[
 R_{j,k}^{(i)}(u)=
    \sum_{\ell=1}^{k d_i}
              \log\left(1+\frac{a_i+\ell/j}{\beta_i}u\right),
 \qquad
 g_{j,k}^{(i)}(u)=R_{j,k}^{(i)}(u)-\frac{k}{n}R_{n,n}^{(i)}(u).
\]
Here each logarithm has value zero at \(u=0\). Cancelling the
row's factorial at index zero gives exactly the factor
\((j\beta_i/u)^{k d_i}\exp R_{j,k}^{(i)}(u)\), since the
signed finite product multiplies or divides by the appropriate
consecutive factors.

Divide each ordinary moment by its scalar tail coefficient at
index zero. Then rescale the original physical argument once by
dividing by the lower tail endpoint scale, namely
\[
 \prod_i (n\beta_i/u)^{(r c_i+s_i)d_i}
             \exp\left(\frac1n\sum_i (r c_i+s_i)R_{n,n}^{(i)}(u)\right).
\]
The ratio of the two orders' explicit power factors is
\(\Lambda^{\sum_i(r c_i+s_i)d_i}=\Lambda^h\). It follows that
the actual rescaled coefficients are
\[
 C_{j,k}\Lambda^{h k\,\boldsymbol 1_{\{j=N\}}}
                \exp\bigl(\tau f_{j,k}(u)+G_{j,k}(u)\bigr),
 \quad
 f_{j,k}(u)=\frac{\sum_i c_i g_{j,k}^{(i)}(u)}{Au},\qquad
 G_{j,k}(u)=\sum_i s_i g_{j,k}^{(i)}(u).                      \tag{2}
\]
These analytic functions have rational Taylor coefficients.
In particular, all the \(G_{j,k}\) vanish at zero.

The identity \(\sum_{\ell=1}^q\ell=q(q+1)/2\) holds with the
signed convention. Applying it to the linear log term gives
\[
 \frac{d}{du}R_{j,k}^{(i)}(0)=
       \frac{k d_i a_i+(k^2d_i^2+k d_i)/(2j)}{\beta_i},
\]
which proves \(f_{j,k}(0)=f_{j,k}^0\). The quadratic log term
gives the explicit first jet
\[
 f_{j,k}'(0)=f_{j,k}^1
  =-\frac{1}{2A}\sum_i\frac{c_i}{\beta_i^2}
     \left[\sum_{\ell=1}^{k d_i}(a_i+\ell/j)^2
                      -\frac{k}{n}\sum_{\ell=1}^{n d_i}(a_i+\ell/n)^2\right].
                                                              \tag{3}
\]
For example, \(|B|\le A\) follows from \(|d_i|\le d_i^2\);
the reference signed drift at the upper order is
\(\exp[-(B/A)\tau k/(2nN)]\).

The two rational vectors in (1) are linearly independent.
Indeed (3) is cubic in \(k\) with cubic coefficient
\(-\sum_i c_i d_i^3/\beta_i^2/(6Aj^2)\). If this coefficient
is nonzero, degree \(N\ge3\) detects it on at least four index
values, whereas \(f^0\) is quadratic. If the coefficient
vanishes, compare the quadratic coefficients at \(j=n,N\).
The ratio of the quadratic coefficient in (3) to that in
\(f^0\) equals
\[
                 -\frac1A\sum_i\frac{c_i d_i^2}{\beta_i^2}
                                   \left(a_i+\frac1{2j}\right).
\]
These two ratios differ since the weights are nonnegative
with a positive total. The isolation test additionally
examines the common-factor locus under these directions.

## Algebraic isolation and exact exclusion

Put \(q=\exp[-\tau/(2nN)]\). As all leading exponents in
the reference are rational, the reference resultant is a
Laurent polynomial in \(q^{1/D_0}\) for some fixed positive integer \(D_0\),
with rational coefficients. It is nonzero. If \(B=0\), its
limit at \(q=0\) is the resultant of two endpoint binomials.
At a common root their respective evaluations of \(z^{nN}\)
would have opposite signs because one of \(n,N\) is odd and
the other even. If \(B\ne0\), remove the upper argument
factor \(\Lambda^h q^{B/A}\) when comparing root moduli.
Before that argument factor, the Gaussian deformations of
both core polynomials have fixed nonzero constant and leading
coefficients and uniformly bounded remaining coefficients for
real \(0<q\le1\). Both sets of root moduli are bounded above
and away from zero by the monic root bound and reversal.
Restoring the upper factor therefore separates the moduli
for all small positive \(q\).

At a positive reference contact, \(q_*\) is algebraic by this
nonzero Laurent polynomial. Thus every \(C_{j,k}^*\) is algebraic,
and \(\tau_*\) is transcendental by Hermite–Lindemann.
The latter claim uses that \(q_*=\exp(-\tau_* /(2nN))\)
and \(\tau_*\ne0\).
[Hermite–Lindemann input](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)

Introduce analytic coefficient arrays
\[
  C_{j,k}^*\exp\left(
       G_{j,k}(u)+v f_{j,k}(u)+t\frac{f_{j,k}(u)-f^0_{j,k}}{u}
              \right).                                      \tag{4}
\]
The quotient is extended at zero. The actual normalized arrays
(2) use \(v=\tau-\tau_*\) and \(t=\tau_*u\).
At \(u=0\), formula (4) agrees with (1) for
\(X=\exp(v/D)\) and \(Y=\exp(t/D)\). This is a biholomorphic
coordinate change locally at the claimed isolated point.

Let \(P\) be the Sylvester determinant of (4) and \(M_i\) all
the minors one size smaller. The arrays keep their degrees
near the point. At a physical common root the even ordinary
moment excludes real zeros, so the conjugate root gives a
common factor of degree at least two. Thus \(P\) and all
\(M_i\) vanish there. The determinant \(P(0,0,v)\) has finite
vanishing order at \(v=0\), since it is the same nonzero
reference resultant. Choose a sufficiently small complex
\(v\)-circle containing just that zero at \(u=t=0\).

Use the contour elimination construction of the
[isolated reflected-template proof](MULTIVARIATE_ISOLATED_TWO_CLOCK_DIFFUSE_EXCLUSION.md#exact-elimination-of-an-accumulating-contact).
For small \(u,t\), multiply the formal linear forms
\(\sum_i w_i M_i(u,t,v_\ell)\) in independent weights \(w_i\)
over all the enclosed zeros \(v_\ell\) of \(P\), counting
multiplicities. The product is polynomial in the weights
and analytic in \(u,t\). Its Taylor coefficients are algebraic:
power sums of the factors are contour integrals of
\((\sum_i w_i M_i)^a P_v/P\), whose successive parameter
coefficients are computable residues of Laurent series in
\(v\) with algebraic coefficients. Newton identities then
recover the product.

An accumulating sequence of actual contacts with
\((u,\tau)\to(0,\tau_*)\) makes that product identically zero
in the weights at every associated \((u,\tau_*u)\).
The analytic identity theorem gives vanishing on that germ.
Its successive coefficients there are polynomials in the
transcendental number \(\tau_*\) with algebraic coefficients.
Each such polynomial must be zero, forcing the entire
two-variable analytic product to vanish. At \(u=0\) this
gives arbitrarily close contacts with \(t\ne0\), since the
enclosed resultant zeros all approach \(v=0\) as \(t\to0\).
By (1) this contradicts isolation.

The test for isolation is finite real quantifier elimination
on the real and imaginary parts of the two complex coordinates
and the algebraic minors. A sufficient test is that the Laurent
resultant of (1) is squarefree after discarding factors that
do not vanish at \((1,1)\). Corank at least two annihilates
both partial derivatives of the determinant by the cofactor
formula. A reduced complex plane curve has isolated singular
points, so this implies the required isolation.

The effectivity proof with a rational
contour radius, a polydisc, and residue computation of the
first nonzero Taylor coefficient proceeds as for the reflected
template. The nonzero coefficient on \(t=\tau_*u\) is found
by a finite search guaranteed to terminate by isolation.
Its nonvanishing is tested before evaluating at \(\tau_*\)
by the algebraic polynomial criterion. Bounds on (4) and
Cauchy bounds on the remaining terms certify a rational
punctured disk, giving computable \(\epsilon,\delta\).

As an all-degree consequence, suppose every positive reference
contact at every paired order below the
[absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
passes the test, and the zero-strength reference pairs
\(Q_n^C(z),Q_N^C(\Lambda^h z)\) are coprime there. Then for
all sufficiently small admissible \(u\) every allowable \(r\)
gives ordinary consecutive coprimality in every degree.
Indeed a hypothetical sequence of collisions has bounded
\(\tau\) by the
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md):
for small \(u\) the exact tail curvature is at least
\(\tau/2-O(u)\) with the constants depending only on the
template. Pass to a convergent subsequence of \(\tau\) at one
of the finitely many low orders. Formula (2) and the
resultant force a reference contact at its limit. The assumed
zero-strength tests and the proved positive exclusion rule
out all such limits. The reference clock roots can be isolated algebraically, and
their logarithms can be enclosed effectively. The uniform
analyticity bounds on compact strength sets then give a
computable threshold.
