# First linked corrections at simple exceptional cubic limits

Use the profile parameters \(\alpha=\ell A_0'(0)\) and
\(\beta=\ell^4 E(0)\) from the
[actual cubic expansion](MULTIVARIATE_PRIVATE_CUBIC_LINKED_DOUBLE_ZERO.md#the-actual-first-two-corrections).
In the one-private endpoint model, write \(\epsilon=n^{-1/3}\),
\(T=\ell^2(1-\sigma)/\epsilon\), \(H=\ell\xi/\epsilon^2\).
The profile is fixed within an adjacent pair. The exact drift
of the physical field at a common root is
\(\xi_{n+1}-\xi_n=\Delta\log(1+1/n)\).

The following conditions apply to bounded scaled fields
and temperatures, with all limiting fields specified below
assumed to be simple zeros of the cubic model. They refine
the [model obstructions](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md)
using the actual first-order row coefficients.

## Actual zero-branch coefficient

On a complex neighborhood of a simple model branch
\(H=h(T)\), let
\[
 p=F_H(T,h(T)),\qquad v=h'(T)+T/2=1/(2p).
\]
The first-order actual interpolating root is
\[
       h_{\rm act}(T,\epsilon)
         =h(T)-\epsilon\bigl(\alpha+\beta Q(T)\bigr)
                                    +O(\epsilon^2),
 \qquad Q=F_{HHHH}/F_H .
 \tag{1}
\]
This is the implicit expansion of the first two terms of
\(J=F+\epsilon(\alpha F_H+\beta F_{HHHH})+O(\epsilon^2)\).
The error and its \(T\)-derivative are uniform on smaller
disks, and its \(\epsilon\)-derivative is \(O(\epsilon)\),
by the differentiated expansion there.
Integer-order roots are exponentially close to these
interpolating roots, since the lattice error is exponentially
small and the local derivative has a fixed nonzero limit.

The cubic equation and backwards heat equation give,
on a zero branch,
\[
 F_{HH}=1-Tp,\qquad
 Q=2(h+T^2)v+2-2Th-T^3,qquad
 v'=-2v^3+Tv^2+h v/2 .                                  \tag{2}
\]
For example, differentiating \(p\) along the branch gives
\(p'=-hp/2+1/(2p)-T/2\), from which the equation for \(v\)
follows.

## Distinct limits at the same imaginary height

Suppose the lower and upper scaled fields of a common root
approach distinct simple zeros \(h_1(T_0),h_2(T_0)\).
Write
\[
 d(T)=h_2(T)-h_1(T),\qquad Q_d(T)=Q_2(T)-Q_1(T).
\]
The fields have the same limiting imaginary height,
so \(d(T_0)\) is real and nonzero. The exact shift
and (1) give the complex necessary condition
\[
               d(T)-\epsilon\beta Q_d(T)
                           =\ell\Delta\epsilon+O(\epsilon^2).
 \tag{3}
\]
Indeed \(T_{n+1}/T_n=(1+1/n)^{1/3}\) when \(T_n\ne0\),
and the change of a bounded temperature is \(O(n^{-1})\)
also at zero. Similarly \(\epsilon_{n+1}/\epsilon
=(1+1/n)^{-1/3}\). These effects on the bounded smooth
branches are \(O(n^{-1})\) in the scaled difference.
The exact drift divided by \(\epsilon^2/\ell\) is
\(\ell\Delta\epsilon+O(|\Delta|\epsilon/n)\).
Here \(|\Delta|\epsilon=O(1)\) because both scaled fields
are bounded. The term \(-\alpha\epsilon\) cancels between
the two actual branches, leaving (3).

In particular \(\ell\Delta\epsilon\to d(T_0)\), and the
imaginary part of (3) requires
\[
                 \Im d(T)
                      =\epsilon\beta\,\Im Q_d(T)+O(\epsilon^2).
 \tag{4}
\]
The real-analytic function \(\Im d(T)\) has a finite order
of vanishing at \(T_0\), since
[distinct model branches cannot persist at a common height](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#isolated-exceptional-temperatures).
If the order is \(m\), (4) gives
\(|T-T_0|=O(\epsilon^{1/m})\). At a transverse crossing,
let \(g'=\Im d'(T_0)\ne0\). Then the necessary conditions
sharpen to
\[
 \begin{split}
 T-T_0&=\epsilon\frac{\beta\Im Q_d(T_0)}{g'}+O(\epsilon^2),\\
 \ell\Delta&=\frac{d(T_0)}\epsilon
       +\beta\left(\frac{\Re d'(T_0)\Im Q_d(T_0)}{g'}
                      -\Re Q_d(T_0)\right)+O(\epsilon).
 \end{split}                                             \tag{5}
\]
All terms here use the same actual profile within the pair.
In particular the drift must have order \(n^{1/3}\) and its
integer value is subject to the second estimate in (5).

## A coincident simple limit of zero order speed

Suppose instead both fields approach the same simple zero
\(H_0=h(T_0)\), and let \(L_0(T)=2h-T h'\).
For nonzero \(\Delta\), the actual same-branch displacement
is \(O(n^{-5/3})\) by (1), which is smaller than the exact
integer drift. For \(\Delta=0\), only a limiting zero of
\(L_0\) needs the following comparison.

At every such zero of \(L_0\) with real \(T_0\) one has
\(T_0\ne0\). There are no real model zeros, and at \(T_0=0\)
the speed is \(2h\ne0\).
Using \(h'=2h/T\) and (2) gives
\[
 L_0'=v-Tv'
   =\frac{(4h+T^2)(16h^2+3hT^2+2T)}{4T^2}
       \quad\text{at }L_0=0 .                            \tag{6}
\]
This derivative is nonzero. The first factor cannot vanish
at a nonreal \(h\). If \(T<0\), the quadratic in the second
factor has real roots by its discriminant \(9T^4-128T>0\).
For positive \(T\), any nonreal roots of that quadratic
have real part \(-3T^2/32<0\), whereas all model zeros
have positive real part by the
[monotone sine-transform bound](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#all-zeros-at-large-positive-temperature).

The true physical root has coordinate
\(\xi=\epsilon^2h_{\rm act}(T,\epsilon)/\ell\) with
\(T\epsilon\) fixed across orders. Its derivative, after
division by \(\epsilon/\ell\), is
\[
     L_0(T)+\epsilon L_1(T)+O(\epsilon^2),\qquad
     L_1=-3\alpha-\beta(3Q-TQ').
]
At a model speed zero, (2) implies
\(3Q-TQ'=6+2(h+T^2)L_0'\), by substituting \(Th'=2h\).
Thus
\[
             L_1(T_0)
              =-3\alpha-6\beta-2\beta(H_0+T_0^2)L_0'(T_0).
 \tag{7}
\]
When \(\Delta=0\) and both original values vanish,
integrate the physical derivative over the adjacent
order interval. The temperature moves by only
\(O(n^{-1})\), and the coefficients in this expansion
are bounded with derivatives near the limiting zero.
After dividing by the nonzero integral of
\(\epsilon/\ell\), the exponentially small zero errors
therefore require
\[
                 L_0(T)+\epsilon L_1(T)+O(\epsilon^2)=0 .
]
By (6) this forces \(T-T_0=O(\epsilon)\).
If \(\alpha_*,\beta_*\) denote the coefficient limits along
a subsequence, the imaginary part of the last equation
divided by \(\epsilon L_0'(T_0)\) gives the necessary
condition
\[
 -2\beta_*\Im H_0+
    3(\alpha_*+2\beta_*)\frac{\Im L_0'(T_0)}
                                  {|L_0'(T_0)|^2}=0 .
 \tag{8}
\]
In particular, at an upper zero, if \(\beta_*<0\) and
\((\alpha_*+2\beta_*)\Im L_0'(T_0)\ge0\), the actual
first-order correction excludes the speed collision.
When \(\beta_*=0\), the sole retained-row case of the
double-zero proof gives \(\alpha_*=1\), so (8) requires
\(\Im L_0'(T_0)=0\).

Equations (3)–(8) are conditions on the actual linked
profile at a simple limiting zero. They leave the
possibility of vanishing transverse corrections at
the remaining exceptional temperatures.
