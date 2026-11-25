# Shared terminal-nine separation with outer mean at least .95

Let the unequal rows of the ordinary exponent columns be fully shared.
Let \(u\) be the order-eight outer geometric mean after endpoint normalization
and removal of the binomial factors. If \(u\ge19/20\), the order-eight and
order-nine factorial moment polynomials have disjoint root arguments under
independent positive rescalings. They are coprime over every
characteristic-zero field.

## Common scale and the two outer means

Put \(p=1-u=\lambda^2\) and \(q=1-w=rp\), where \(w\) is the order-nine
outer mean. The checker uses
\[
 \lambda\le\frac{22361}{100000},\qquad
 \frac{507}{500}\le r\le\frac{10177}{10000},\qquad q\le q_0=\frac{13}{250}.
\]
For the matched reciprocals \(x=d/(8a+h)\), \(x'=d/(9a+h)\), one has
\(x'=8x/(9-\gamma x)\) with \(\gamma=h/d\in(0,1]\).
Also \(x\le(u^{-2}-1)/7\).
The ratios \(8x'/(7x)\) are at least \(64/63\). Expansion of the
outer product gives \(w\le u(1+(1-u^2)/63)^{-1/2}\), hence
\[
     r\ge1+\frac{u(1+u)}{126(1+(1-u^2)/63)^{3/2}}.
\]
The ratios \(8x'/(7x)\) are at most \(64/[7(9-(u^{-2}-1)/7)]\).
Apply \((1+b)^c\ge1+cb\) to the outer product at that upper ratio \(c\).
This gives \(w\ge u^c\ge1-c(1-u)\) and \(r\le c\).
The rational checks for both endpoints use the smallest allowed \(u\).
These bounds also hold when all reciprocals in the two matched lists are
scaled down by the same positive factor. Their ratios stay fixed, and the
lower bound on the new \(u\) only improves.

Use scale \(\lambda\) at order eight. At order nine use its own scale
\(\mu=\sqrt q=\sqrt r\,\lambda\) with \(\mu\le113/500\).
In the relevant scale \(z_0=\lambda\) or \(\mu\), use
\(t=(z_0Z-1)/(z_0Z+1)\) and multiply the normalized polynomial
by \((1+z_0Z)^j/(2^jz_0^j)\).
The order-eight coefficients and their error enclosures come from the
[within-order octic bounds](MULTIVARIATE_TERMINAL_EIGHT_NEAR_FLAT.md#octic-linked-means-and-skews)
at deficit \(p\).

## Nonic means, signed skews and their coupling

Write the nonic multipliers at \(i,9-i\) as
\(w^{e_i}\exp(D_i\pm\alpha_i)\) for \(0\le i\le4\), with
\(e=(0,1,7/4,9/4,5/2)\). Let \(s=\tfrac12\log(1+8x')\) on a single
Gauss factor. The mean-loss formula from the
[shared coefficient envelope](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#joint-coefficient-bounds)
gives
\[
 D=H U+(0,0,0,J_3,J_4),\qquad U=(0,0,1,15/7,20/7),
\]
where \(G_a(s)=as-\operatorname{artanh}(a\tanh s)\),
\(H=\sum G_{3/4}\), \(J_3=\sum(G_{1/2}-\tfrac87G_{3/4})\), and
\(J_4=J_3+\sum(G_{1/4}-\tfrac57G_{3/4})\).
Writing \(C_a=a(1-a^2)\), the derivative
\(G'_a=C_a\tanh^2 s/(1-a^2\tanh^2 s)\) gives
\(0\le H\le C_Hq^3\), \(J_i\le0\), and
\[
 |J_3|\le C_{1/2}'q^5,\qquad
 |J_4|\le(C_{1/2}'+C_{1/4}')q^5,
\]
with \(s_*=q_0/(1-q_0)\) and
\[
 C_H=\frac{C_{3/4}}{3(1-(3/4)^2s_*^2)(1-q_0)^3},\qquad
 C_a'=\frac{C_a((3/4)^2-a^2)}
   {5(1-(3/4)^2s_*^2)(1-a^2s_*^2)(1-q_0)^5}.
\]
Indeed \(\sum s=-\log(1-q)\le q/(1-q_0)\), and the differences
of derivatives have their quadratic terms canceled. Each \(D_i\)
is a sum of nonnegative mean defects.

For skew, center the Gauss logarithms at \(m=4\) and write
\(\tau=x'/(1+4x')\). The logarithm at \(m\) is the logarithm at
\(4\) plus \(\log(1+(m-4)\tau)\).
The skew weights are symmetric in \(m,8-m\), so only even powers survive.
For either orientation choose the sign \(\sigma\) so that the skew of
this factor at index \(i\) has the expansion \(\sigma\sum_{h\ge1}K_{h,i}\tau^{2h}\),
where
\[
 K_{h,i}=\frac1{4h}\left(
 \frac{9-2i}{9}\sum_{m=0}^8(m-4)^{2h}
                 -\sum_{m=i}^{8-i}(m-4)^{2h}\right).
\]
Put \(B=(2-q_0)/(8(1-q_0)^2)\), \(A_i=K_{1,i}\), \(L_i=K_{2,i}\),
and \(d_h=\sum\sigma\eta^h\) with \(\eta=\tau/(Bq)\).
Since \(\sum x'\le((1-q)^{-2}-1)/8\), one has \(\sum\eta\le1\).
The first two skew terms are
\(\alpha_{1,i}=A_i B^2q^2d_2\) and
\(\alpha_{2,i}=L_i B^4q^4d_4\).
The remaining absolute error is at most \(C_\alpha q^6\), where
\[
                 C_\alpha=\frac{9(4B)^6}{12(1-(4Bq_0)^2)} .
\]
This follows from \(|K_{h,i}|\le9\cdot4^{2h}/(4h)\) and summing
the geometric tail. The index zero has identically zero skew.

The signed moments obey \(|d_2|,|d_4|\le1\) and
\[
 |d_4-d_2|\le\sum(\eta^2-\eta^4)\le2/5.
\]
For the last step, \(0\le\eta\le1\) and
\(\eta-\eta^3\le2/(3\sqrt3)<2/5\).
There is also a lower bound on the mean defect at fixed skew:
\[
   \frac{H}{q^3}\ge C_{\rm link}d_2^2,\qquad
   C_{\rm link}=\frac{C_{3/4}}3(1-s_*^2/3)^2(4B)^3 .
\]
In fact \(\tanh t\ge t-t^3/3\) for \(0\le t\le s_*\), and
\(s=\operatorname{artanh}(4\tau)\ge4\tau\). The derivative formula
bounds \(G_{3/4}(s)\) below by
\(C_{3/4}(1-s_*^2/3)^2(4\tau)^3/3\).
Cauchy–Schwarz gives \(\sum\eta^3\ge(\sum\eta^2)^2\ge d_2^2\).

## Cancellations before Cayley scaling

For a mean contribution to the Cayley polynomial, expand
\[
 e^{D_i}\cosh\alpha_i
                  =1+D_i+\tfrac12\alpha_{1,i}^2+q^6\epsilon_i .
\]
For a skew contribution, expand
\[
 e^{D_i}\sinh\alpha_i
         =\alpha_{1,i}+\alpha_{2,i}+H U_i\alpha_{1,i}
              +q^6\delta_i .
\]
The mixed term \(H U_i\alpha_{1,i}\) is kept as one shared direction
across the multipliers. With \(v=i(9-i)\), one has
\(U_i=v(v-8)/84\) and \(A_i=v(9-2i)/12\). The mixed vector
\(U_i A_i\) has degree seven in the coefficient index.

Here are the error bounds used by the checker. For each \(i\) let
\(c=U_i C_H\), \(a=A_i B^2\),
\(b=L_i B^4+C_\alpha q_0^2\), and \(a_*=a+bq_0^2\).
Thus \(0\le D_i\le c q^3\),
\(|\alpha_{1,i}|\le a q^2\), \(|\alpha_i-\alpha_{1,i}|\le bq^4\),
and \(|\alpha_i|\le a_*q^2\).
Write \(j_i\) for the bound of \(|D_i-HU_i|/q^5\).
The checker uses
\[
 \begin{aligned}
 |\epsilon_i|&\le c^2+(a_*+a)b/2+a_*^4q_0^2/20
                              +(501/500)c a_*^2q_0,\\
 |\delta_i|&\le C_\alpha+a_*^3/5+q_0(j_i a+cb)
       +q_0^2c^2(a_*+a_*^3q_0^4/5)+q_0^3c a_*^3/5 .
 \end{aligned}
\]
The index-zero errors vanish. All arguments of the exponentials and
hyperbolic functions have modulus below \(1/100\).
Use \(|e^{D}-1-D|\le D^2\), \(|e^D-1|\le2D\),
\(\cosh\alpha-1\le(501/1000)\alpha^2\),
\(|\cosh\alpha-1-\alpha^2/2|\le\alpha^4/20\), and
\(|\sinh\alpha-\alpha|\le|\alpha|^3/5\).
In the mean error, separately bound \((\alpha^2-\alpha_1^2)/2\)
and the product of the exponential and hyperbolic corrections.
In the skew error use
\(D_i\alpha_i-HU_i\alpha_{1,i}=(D_i-HU_i)\alpha_{1,i}
+D_i(\alpha_i-\alpha_{1,i})\).
The remaining terms are the skew-log tail, the hyperbolic cubic error
with its product by \(D_i\), and the exponential quadratic error times
\(\sinh\alpha_i\).

Before reflection pairing, the Cayley weight of coefficient index
\(i\) in degree \(h\) is
\[
 2^{-9}\binom9i [Z^h](Z-1)^i(Z+1)^{9-i}.
\]
Pairing doubles this weight in odd degrees for the means and in even
degrees for the skews. Each direction is thus a rational combination
of \((1-q)^{e_i}\) times its displayed \(q\)-power.
The checker forms the first nine Taylor coefficients of the combination
exactly, removes the initial zeros, and bounds the tail by the sum of
the absolute ninth coefficients times a geometric series. This bound
uses the decreasing absolute binomial coefficients after degree nine
for \(0\le e_i\le5/2\). Substituting \(q=\mu^2\) before division shows
that every resulting power of the nonic scale is nonnegative. The base
direction gives the only terms at zero. The independent error arrays
have raw power \(q^6\), so they do not require cancellation.

## Exact ray exclusions

Leading coefficient dominance confines the roots to \(|Z|<2\) in their
respective scales. Use seven rays with arguments
\(\pi-2\arctan(\lambda b)\), where the rational levels
\(b=b_0+b_1p\) depend affinely on \(p\).
Set \(b_*=b\) in the octic coordinate. In the nonic coordinate
use \(b_*=(b_0+b_1q/r)/\sqrt r\), which gives the same angle since
\(\mu b_*=\lambda b\). For deficit \(h=p\) or \(q\) and the
corresponding ray level \(b_*\), the ray has equation
\[
                    (Y-b_*)(1+h b_*Y)=-h b_*x^2 .
\]
The checker encloses \(Y\) by three interval updates. Keeping the nonic
coefficients in their own scale makes the ratio interval enter only this
ray equation. Reversing both columns replaces \(x\) by \(-x\) at fixed
\(Y\), so the nonnegative half suffices.

For a nonic projection, use \(d_4=d_2+(d_4-d_2)\) and
\(|d_4-d_2|\le2/5\). This leaves a real quadratic in \(d_2\) together
with the projected term \(H/q^3\). Where that term's sign calls for
a lower bound, include \(C_{\rm link}d_2^2\) in the quadratic.
The quadratic range on \([-1,1]\) is computed by its two half-intervals,
checking endpoints and any interior vertex. All intervals and tests
of strict projected sign are rational. A proposed projection may also
be taken perpendicular to the midpoint of the leading skew direction
on the cell. This only selects the real linear functional, and the
strict sign on the whole cell is checked separately. For an interval
rectangle centered at \(Z_c\), the checker can evaluate each coefficient
family through \(P(Z_c)+P'(Z_c)(Z-Z_c)\) with absolute remainder at most
\(\sup|P''|\,|Z-Z_c|^2/2\). The supremum uses the termwise
absolute coefficient envelope and the maximum radius of the rectangle.
The same bound applies to both real and imaginary parts.

At zero deficit the continuous scaled anchors are the imaginary
Hermite polynomials of degrees eight and nine with variance
\(1/14\) and \(1/16\) in their respective scales. The latter
variance is \(r/16\) in the octic scale. The checker obtains the
exact nonic anchor coefficients from the leading coefficients of
the nominal Taylor families. Every correction family vanishes at zero.
The checker isolates the positive squared imaginary heights and compares them
with the squared ray heights at zero deficit. The upper-root order
from the smallest argument is \(9,8,9,8,9,8,9,8\).
Scaling both reciprocal lists together to zero preserves the ratio
bounds and yields the scaled Hermite anchors, so the ray exclusions
preserve the counts. In the sector joined across the negative real
axis, the sole order-nine root is real by conjugacy. The order-eight
pair is nonreal because the even factorial moment is strictly
positive on real arguments. This proves disjointness of arguments.

The [checker](verify_multivariate_terminal_nine_near_flat.py) visits
53,298 cells through all fourteen ray checks. Run it from the repository
root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_nine_near_flat.py
~~~
