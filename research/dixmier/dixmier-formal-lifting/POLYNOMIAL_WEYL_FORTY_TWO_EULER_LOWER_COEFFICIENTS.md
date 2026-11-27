# Lower Euler components at the simple degree-forty-two roots

**Theorem PolynomialWeylFortyTwoEulerLowerCoefficients.** Use the
[actual degree-forty-two root constraints](POLYNOMIAL_WEYL_FORTY_TWO_SECONDARY_ROOT_CONSTRAINTS.md).
At any simple root \(\delta^3=v\), use their exact shear with parameter
\(\gamma\), their Euler variable \(C=3t^2(E-\delta)\), and the resulting
Laurent coefficients \(\widehat P_k(E)\). Put \(A=\alpha H'(\delta)^2\).
The leading Euler component is \(t^2 f(C)\), where
\[
 f(C)=\frac A9 C(C-2)+\frac{\widehat P_4'(\delta)}3 C+
              \widehat P_2(\delta).
\]
Let \(c\) be the mean of the two roots of \(f\), and \(m\in\{1,3,5,7\}\)
their positive integer gap. The polynomial of the next Euler component
is
\[
 p(C)=\widehat P_1(\delta)+\frac{\widehat P_3'(\delta)}3 C+
              \frac{\widehat P_5''(\delta)}{18} C(C-2).
\]
Then
\[
             p(c+1/2)+\frac{m^2-1}{16}\,p''=0.              \tag{1}
\]
For \(m\in\{3,5,7\}\), the third lower Euler component must also
satisfy the linear test (7) after the explicit correction there.
Equation (1) determines \(P_1(\delta)\) from the higher Laurent
coefficients and the previously allowed \(P_2(\delta)\). For any
of the three gaps covered by (7), that test determines
\(P_{-1}(\delta)\) from the higher coefficients.

## A linear obstruction for the Euler components

We use \(x=C-r\) with \(r=c-m/2\), so that the roots of the
monic leading polynomial are \(0,m\). In this section take more
generally any odd \(b\ge3\) with \(m\) odd and \(1\le m\le b\).
The commutation in the Euler algebra is
\[
 [x,t]=t,\qquad a(x)t^j=t^j a(x+j) .
\]
Suppose leading components \(t^2 f_0(x), t^b g_0(x)\) commute
and \(f_0(x)=x(x-m)\), with \(g_0\) monic of degree \(b\).
The root sets of \(g_0\) are the odd and even progressions
\[
                  2-b,4-b,\ldots,m,\qquad
                  m+2-b,m+4-b,\ldots,0.                    \tag{2}
\]
The even progression is empty if \(m=b\). Indeed, the
multiplicity polynomial of its roots in the integer orbit is
\[
 z^{2-b}\frac{(1+z^m)(1-z^b)}{1-z^2}.
\]
Its roots as a multiplicity record are given by the two
progressions in (2), each with multiplicity one. The resulting
monic polynomial satisfies
\(f_0(x+b)g_0(x)=f_0(x)g_0(x+2)\).
Every other monic polynomial solution agrees with it because
the quotient is a rational function of period two. A rational
function with nonzero period has no poles and is constant.

For odd \(d\) with \(1\le d\le m\), consider polynomial solutions
\(a,q\) of the linearized equation
\[
 f_0(x+b-d)q(x)-f_0(x)q(x+2)+
       a(x+b)g_0(x)-a(x)g_0(x+2-d)=0.                     \tag{3}
\]
Then \(a\) is in the image of
\[
       L_d(h)=h(x+2)f_0(x)-f_0(x-d)h(x).                  \tag{4}
\]
For every polynomial \(h\), the pair of commutators of
\(t^{-d}h(x)\) with the two leading components gives a solution
of (3) with \(a=L_d(h)\), by the Jacobi identity. If \(h\)
has degree \(e\), \(L_d(h)\) has degree \(e+1\), with
leading coefficient \(2(e+d)\) times that of \(h\). Thus
every polynomial \(a\) can be reduced modulo this image to a
constant. It suffices to rule out \(a=1\).

Set \(x_0=m-b+d\). For \(x=x_0,x_0+2,\ldots,m\), (2) gives
\(g_0(x)=0\), so at \(a=1\) the recurrence is
\[
         f_0(x+b-d)q(x)-f_0(x)q(x+2)=g_0(x+2-d).
\]
The initial term \(f_0(x_0+b-d)q(x_0)\) vanishes.
Up to \(x<d\) the right side vanishes by the even progression
in (2), and \(f_0(x)\ne0\) if \(x<m\) because \(x\) is odd.
Thus the successive values of \(q\) obtained in this part are
zero. For the remaining \(x=d,d+2,\ldots,m-2\),
\(f_0(x)<0\), and \(f_0(x+b-d)>0\) except when \(x=x_0\),
where it is zero. The numbers \(g_0(x+2-d)\) at these
positive even arguments are nonzero and alternate signs,
also alternating in sign at the next argument for \(x=m\).
At each recurrence step the contributions
\[
   \frac{f_0(x+b-d)}{f_0(x)}q(x)
             -\frac{g_0(x+2-d)}{f_0(x)}
\]
have the same sign when both are nonzero. Hence, if \(m>d\),
\(q(m)\) has the sign opposite to \(g_0(m+2-d)\).
If \(m=d\) and \(b>d\), the first part gives \(q(m)=0\).
Both contradict the equation at \(x=m\), since
\(f_0(m+b-d)>0\) when \(b>d\).
If \(m=d=b\), the equation at \(x=m\) already requires
\(0=g_0(2)\), impossible by (2).
All recurrence values whose signs were used are determined
by (3) and the rational monic \(f_0,g_0\), so the sign argument
applies over \(\mathbb C\).

## Application to the first and third lower components

At the simple root of the degree-forty-two pair the Euler
weights of \(t,E-\delta\) are \(1,-2\). They give an exact
grading, and the leading weights are \(2,7\), so \(b=7\).
Divide \(P\) and \(Q\) by their respective leading Euler
polynomial coefficients and translate \(C=x+r\).
Write their components at weights \(2-i\) and \(7-i\) as
\(t^{2-i}p_i(x)\) and \(t^{7-i}q_i(x)\). In particular
\(p_0=f_0\), \(q_0=g_0\), and \(p_1(x)=9p(x+r)/A\).
The weight-eight part of the commutator is (3) with
\(d=1\), \(a=p_1\), and \(q=q_1\).
Since \(\deg p_1\le2\), (4) gives \(p_1\) in the span of
\[
 L_1(1)=2x-(m+1),\qquad L_1(x)=4x^2-(3m+1)x .
\]
Set \(n=(m+1)/2\). On polynomials \(u_0+u_1x+u_2x^2\) this
plane is the kernel of
\(u_0+nu_1+n(3n-1)u_2/2\).
Since \(r+n=c+1/2\), this is exactly (1).

To apply (3) at \(d=3\), first normalize the lower components.
All conjugations in this argument can be made in the Euler
completion, with weights bounded above and polynomial
coefficients at each weight. For each weight the computations
involve finitely many terms. Choose the unique polynomial \(h_1\)
of degree at most one with \(L_1(h_1)=-p_1\), and apply
\(\exp(\operatorname{ad}(t^{-1}h_1(x)))\) to both elements.
This kills the first lower \(P\) component. For a polynomial
\(u\) in a component of weight \(a\), set
\[
           D_a(u)=h_1(x+a)u(x)-u(x-1)h_1(x).              \tag{5}
\]
In the conjugated \(P\) the second and third lower
polynomials are
\[
 p_2+\tfrac12 D_1(p_1),\qquad
 N_3=p_3+D_0(p_2)+\tfrac13D_0(D_1(p_1)).
\]
These follow by expanding the exponential and using
\(D_2(f_0)=-p_1\). Remove the second lower component
as follows. Any polynomial is a scalar modulo the image of
\(L_2\), by the leading-degree calculation for (4).
A conjugation by \(\exp(\operatorname{ad}(t^{-2}h_2(x)))\)
leaves the first and third lower \(P\) components unchanged
and places the second in the scalars. Subtract that scalar
from \(P\). These operations leave \(N_3\) unchanged.

With the first lower \(P\) component zero, equation (3)
at \(d=1\) says that the first lower \(Q\) component
commutes with \(t^2 f_0(x)\).
The weight-six homogeneous polynomial solutions of this
commutation are multiples of the cube of \(t^2 f_0(x)\).
Indeed that cube is a solution, and the quotient of any
other polynomial solution by its coefficient is a rational
function of period two.
Subtract a scalar multiple of the cube of the normalized
\(P\) from \(Q\) to remove this component, preserving the
scalar commutator. At \(d=3\), the lower terms in its
equation now vanish, so (3) applies to \(a=N_3\).
If \(m\ge3\) this gives \(N_3\in\operatorname{im} L_3\).

For \(m\in\{3,5,7\}\), a linear function annihilating this
image is obtained by weighting the evaluations at
\(x=3,5,\ldots,m\) with \(\omega_x\) satisfying
\[
       \omega_{x+2} f_0(x-1)=\omega_x f_0(x).             \tag{6}
\]
The sum telescopes in (4), with zero coefficients at the
two endpoints. Use the unnormalized weights \(1\) for
\(m=3\), \(1,1\) for \(m=5\), and \(5,6,5\) for \(m=7\).
The third-component test is therefore
\[
               \sum_{x=3,5,\ldots,m}\omega_x N_3(x)=0.   \tag{7}
\]
For these first four \(P\) components, the original coefficients
enter as
\[
 p_i(x)=\frac9A\sum_{j\ge0}
     \frac{\widehat P_{2-i+2j}^{(j)}(\delta)}{3^j j!}
            \prod_{\ell=0}^{j-1}(x+r-2\ell),\qquad 0\le i\le3.
\]
This follows from \((E-\delta)^j=t^{-2j}\prod_{\ell=0}^{j-1}
(C-2\ell)/3^j\). Only \(p_3\) uses \(\widehat P_{-1}(\delta)\).
Its coefficient in \(N_3(x)\) is the constant \(9/A\), so
(7) determines it because the weights sum to a nonzero number.
The shear lowers Laurent degree away from the diagonal Taylor
term, so \(\widehat P_{-1}(\delta)-P_{-1}(\delta)\) depends only on
higher Laurent coefficients. The same reasoning at degree one
gives the assertion about \(P_1(\delta)\) from (1).
With the [fourth-root scalar restrictions](POLYNOMIAL_WEYL_FORTY_TWO_SCALAR_FACE_COEFFICIENT_DESCENT.md),
equation (1) fixes the value at \(v\) of the polynomial
\(C_1\) that also satisfies \(C_1(u)+uC_1'(u)=0\).

## All odd components up to the Euler gap

**Corollary PolynomialWeylFortyTwoEulerOddHierarchy.** In the same
Euler completion, after a scalar shift of the normalized \(P\),
conjugations of strictly negative weights put it in a form
whose components at \(2-d\) vanish for all odd \(1\le d\le m\).
For \(m=5,7\), the resulting test at \(d=5\) fixes
\(P_{-3}(\delta)\) from higher Laurent coefficients.
For \(m=7\), the test at \(d=7\) also fixes \(P_{-5}(\delta)\).
At each spectral value the original sheared \(P\) has two
formal Laurent eigenvectors at the simple root on the lattice
of exponents \(r+\mathbb Z\), with largest exponents \(r\)
and \(r+m\), respectively.

The normalization above kills the first two lower components
of \(P\), and the first lower component of \(Q\). Iterate
at odd \(d\le m\), starting from \(d=3\). When all lower
odd \(P\) and \(Q\) components vanish, the part of the
commutator at \(2+b-d\) is (3). Each middle term has
index pair \(i+j=d\), so one index is a lower odd index.
Consequently \(p_d\in\operatorname{im}L_d\). Conjugate by
\(\exp(\operatorname{ad}(t^{-d}h_d(x)))\) with
\(L_d(h_d)=-p_d\). This kills \(p_d\) and preserves all
components at earlier indices. The component \(q_d\) then
commutes with \(t^2 f_0(x)\). Its weight \(b-d\) is even
and nonnegative, so it is a scalar multiple of the leading
component of \(P^{(b-d)/2}\). The quotient argument using
a rational function of period two proves this uniqueness,
including weight zero. Subtract that multiple of the actual
power from \(Q\). The commutator and the earlier components
are unchanged. This completes the induction since
\(2+b-d\ge2\) for all \(d\le m\).

For a test at a particular odd \(d\), let \(N_d\) be
the \(P\) polynomial at \(2-d\) just before its conjugation.
The functional killing \(\operatorname{im}L_d\) is
\[
  \sum_{x=d,d+2,\ldots,m}\omega_x N_d(x)=0,\qquad
  \omega_d=1,\qquad
  \omega_{x+2} f_0(x+2-d)=\omega_x f_0(x).                 \tag{8}
\]
The weights are positive since both entering \(f_0\) values
are negative when \(d\le x<m\).
In particular the new tests are \(N_5(5)=0\) at \(m=5\),
\(N_5(5)+N_5(7)=0\) at \(m=7\), and \(N_7(7)=0\) at
\(m=7\). Since the prior conjugations strictly lower weight,
they add to \(N_d\) terms formed from components before
index \(d\). An original component earlier than \(d\)
involves only \(\widehat P_k\) at \(k>2-d\). The dependence
on \(\widehat P_{2-d}(\delta)\) in \(N_d\) is thus exactly
the constant coefficient \(9/A\) times that value. This
proves the assertions for \(d=5,7\), using again that the
shear changes a Laurent coefficient only by terms from
higher coefficients.

For the eigenvectors, use the formal module with basis
\(t^{r+j}\) for integer \(j\), allowing series with
exponents bounded above. Here \(x t^{r+j}=j t^{r+j}\).
Write a normalized eigenvector with largest exponent
\(r+m\) as \(\sum_{n\ge0}u_n t^{r+m-n}\), \(u_0=1\).
In the normalized equation \((P-\mu)u=0\), the coefficient
at exponent \(r+m+2-n\) has diagonal term \(f_0(m-n)u_n\).
The spectral term \(-\mu\) has index two. Once the odd
components through index \(m\) are zero, all odd \(u_n\)
with \(0<n<m\) vanish in the recursion. At \(n=m\) the
diagonal coefficient vanishes, but the equation only uses
the already zero odd coefficients or the zero odd components.
The remaining diagonal coefficients are nonzero. The series
with largest exponent \(r\) has no obstruction in its
descending recursion. These two series are independent,
and every series satisfying the equation is their linear
combination, by its largest exponent and the same recursion.
The negative-weight conjugations act on this module by
triangular series transformations that preserve the largest
exponent. Undoing them, the scalar shift and leading
coefficient rescaling proves the statement for the
original sheared \(P\).
