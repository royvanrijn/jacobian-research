# A Sylvester normal test at multiple Gaussian contact roots

The second-order signed-tail test for several simultaneous simple Gaussian roots can be formulated directly in the quotient by the reference gcd. The resulting criterion allows repeated roots of that gcd. It also allows a candidate actual common factor to converge to a proper divisor that uses only part of a repeated reference block. It uses actual ordinary shared tails of arbitrary heights, gaps, and multiplicities.

## Linearized rank of the full common-factor condition

Let \(F_*,G_*\in\mathbb R[z]\) be fixed polynomials of degrees \(n,N\), with monic gcd \(H_*\) of degree \(d_*\). Write
\[
                         F_*=H_* A,\qquad G_*=H_* B .
\]
In this section \(F=F_*+f\) and \(G=G_*+g\) have the same degrees as the references and are close in coefficient norm. For a polynomial \(U\), write \(M_U\) for multiplication by \(U\) on \(\mathbb R[z]/(H_*)\). There is a choice of fixed domain and range bases for the Sylvester map in which the exact rank in excess of \(n+N-d_*\) is the rank of a \(d_*\)-dimensional matrix
\[
                        M_{Bf-Ag}+O((\|f\|+\|g\|)^2).         \tag{1}
\]
The remainder is effective in a neighborhood of the fixed pair \(F_*,G_*\). For the general normal-block rank condition, see the [real determinantal tangent cone](https://doi.org/10.1007/s11228-026-00796-4). The multiplication operator in (1) comes from the Sylvester structure.

To prove this, use the Sylvester map \((u,v)\mapsto Fu+Gv\) with multiplier degrees below \(N,n\), respectively. At the reference its kernel consists of \((B\phi,-A\phi)\) with \(\deg\phi<d_*\), and its image is the space of multiples of \(H_*\) of degree below \(n+N\). Use the reference kernel for the last domain block and identify the reference cokernel with remainders modulo \(H_*\). The upper left block remains invertible near the reference. The linear part of its Schur complement is the induced change from kernel to cokernel, namely \(\phi\mapsto (Bf-Ag)\phi\bmod H_*\). The correction is quadratic because both off-diagonal blocks vanish at the reference. This proves (1).

Multiplication by \(U\) on the quotient has rank \(d_*-\deg\gcd(U,H_*)\), including when \(H_*\) has repeated roots. Indeed its kernel consists of the multiples of \(H_*/\gcd(U,H_*)\) in degree below \(d_*\). For actual \(\deg\gcd(F,G)\ge d\), the matrix in (1) must have rank at most \(d_*-d\). This yields a condition on the polynomial \(Bf-Ag\) without requiring any individual perturbed root branch. When \(d>d_*\), continuity of the invertible reference block already excludes the proposed degree near the reference.

## The factorial normal directions

Fix a finite ordinary core, adjacent orders \(n,N=n+1\), and the signed integer gap \(h\) of the appended shared rows. Consider a positive physical Gaussian reference \((W_*,\lambda_*)\). Write its lower and upper polynomials in the common lower-endpoint-normalized variable as \(F_*,G_*\), using the relative upper scale \((N/n)^h e^{-\lambda_*/(2nN)}\). The reference polynomials, and the real monic gcd \(H_*\), have algebraic coefficients at the [finite positive reference contacts](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md). This normalization is a common positive argument scaling of the actual parents, so it preserves common-factor degrees. The even reference has no real roots by the positive Gaussian mixture, and \(H_*\) has no real roots.

Use the [midpoint signed tail statistics](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md) \(W,\lambda,K,Z,B_4,L_2,\eta\). In particular \(Z>0\) for a nonempty unequal tail, and near the reference
\[
                  |K|,Z=O(\eta),\quad
                  B_4\le WZ,\quad |L_2|\le Z .
\]
Put \(x=W-W_*,\ y=\lambda-\lambda_*\). The signed expansion for the two coefficient arrays, with the relative upper scaling included, has terms linear in \(x,y,K,Z,B_4,L_2\), with coefficient remainder
\[
                   O\bigl(\eta Z+(|x|+|y|+|K|+Z)^2\bigr).     \tag{2}
\]
All constants are uniform over the ordinary tail in a fixed neighborhood of this reference. For the lower and upper array derivatives, the log weights on coefficient \(k\) at order \(j\) are
\[
 \begin{array}{c|cc}
   &j=n&j=N\\ \hline
 x& F_j(k)&F_j(k)\\
 y&0&-k/(2nN)\\
 K& H_j(k)&H_j(k)\\
 Z& G_j(k)&G_j(k)\\
 B_4& I_j(k)&I_j(k)\\
 L_2&0&c k
 \end{array}
\]
where \(F_j(k)=(k^2-jk)/(2j)\), \(G_j=-F_j/(2j)\),
\(H_j(k)=(k/j-1/2)k(j-k)/(6j)\),
\(I_j(k)=j((k/j-1/2)^4-1/16)/12\), and
\(c=(n^{-2}-N^{-2})/12\).
These are the coefficient formulas of the cited signed expansion. Formula (2) follows by expanding the exponentials around the fixed reference weights. The error \(\eta Z\) is their coefficient-log remainder, and the squared term bounds the remaining Taylor terms.

For \(t\in\{x,y,K,Z,B_4,L_2\}\), let \((F_t,G_t)\) be its two derivative arrays obtained by multiplying each reference coefficient by the given weight. Define
\[
                    \Omega_t=B F_t-A G_t\bmod H_* .
\]
Here \(A,B\) are the reference cofactors. In particular
\(\Omega_{L_2}=-2nN c\,\Omega_y\).

## A finite rank and interval exclusion test

Fix an even target common-factor degree \(d\le d_*\). For every monic real divisor \(J\mid H_*\) of degree \(d\), write brackets for reduction modulo \(J\). These divisors include all choices of allowable partial multiplicities of repeated conjugate pairs. The list is finite and computable with real algebraic coefficients. Make the following two checks for each \(J\).

1. The three vectors \([\Omega_x],[\Omega_y],[\Omega_K]\) are real-linearly independent in \(\mathbb R[z]/J\).
2. Let \(P_J\) be the orthogonal projection off their real span, using the standard coefficient basis. With \(M_*=W_*-|\lambda_*|\), the segment
   \[
              \left\{P_J[\Omega_Z]+uP_J[\Omega_{B_4}]:
                                     0\le u\le M_*/2\right\} \tag{3}
   \]
   has positive distance from zero.

**Theorem.** When both checks pass at the reference, no ordinary shared tail sufficiently near it and sufficiently diffuse can have common-factor degree at least \(d\). The neighborhood and diffuseness threshold are effective. This applies to multiple reference roots and to a candidate that uses only a proper part of a repeated reference root cluster.

First suppose an exact common factor of degree at least \(d\) occurs and set \(s=|x|+|y|+|K|\). For \(s>0\) consider
\[
      M_{\Omega_x(x/s)+\Omega_y(y/s)+\Omega_K(K/s)} .
\]
Its rank exceeds \(d_*-d\) uniformly on the compact unit coordinate sphere \(|x/s|+|y/s|+|K/s|=1\). Otherwise the multiplier polynomial would be divisible by some real degree-\(d\) divisor \(J\) of \(H_*\). Such a divisor exists whenever the common divisor has degree at least \(d\), because it has no real roots and \(d\) is even. This contradicts check 1. The first forbidden singular value therefore has a positive minimum. By (1)–(2) and the exact rank bound one obtains \(s\ll Z+(s+Z)^2\), hence \(s=O(Z)\) in a small neighborhood.

Now \(K^2/Z=O(Z)\). The pure moment constraints in the [second-order signed-tail polytope](MULTIVARIATE_TRANSVERSE_TAIL_RESPONSE_POLYTOPE.md#the-joint-skew-quartic-and-center-polytope) give
\(0\le B_4/Z\le M_*/2+O(\eta)\), also when \(M_*=0\).
This use of its inequalities does not require simple reference roots. After dividing the Schur complement by \(Z\), (1)–(2) give, with error \(O(\eta)\), multiplication by
\[
        \Omega_x x/Z+\Omega_y y/Z+\Omega_K K/Z
                    +\Omega_Z+\Omega_{B_4}B_4/Z+\Omega_{L_2}L_2/Z .
                                                               \tag{4}
\]
The first three normalized coordinates are bounded. The last term vanishes after applying every \(P_J\). Check 2 therefore says that all the limiting multipliers in the resulting compact set have rank exceeding \(d_*-d\). The positive forbidden singular-value margin contradicts the exact rank bound for sufficiently small \(\eta\).

The two strict margins and the coefficient bounds are effective. For the margins one can use real algebraic rank tests and compact semialgebraic bounds, enclosing \(W_*,\lambda_*\) at a precision certifying the given strict distances. The Schur remainder constants follow from an invertible reference minor and its bounded inverse near the reference.

For comparison, at simple reference roots reduction modulo \(J\) is evaluation at distinct roots and recovers the simultaneous simple-branch condition. Here the quotient retains derivatives up to the selected multiplicities. In particular, if all positive references of a fixed ordinary core on a compact strength range pass these two checks for \(d=4\) whenever \(d_*\ge4\), sufficiently diffuse ordinary tails on that range can have at most one irreducible quadratic common factor. This follows by the finite reference theorem and the absolute eventual cutoff. No simple-contact hypothesis is needed for this conditional quadratic reduction. The criterion does not itself exclude the remaining degree-two factors.

## Coefficient response for a whole common block

The normal rank test also yields a shape interval when the selected
limiting divisor \(J\mid H_*\) takes whole primary blocks, meaning
\(\gcd(J,H_*/J)=1\). This does not require the roots of \(J\) to be
simple. Suppose a monic actual common factor \(H\) tends to \(J\).
The two reference cofactors \(A_0=F_*/J,\ B_0=G_*/J\) generate the unit
ideal modulo \(J\): at each root of \(J\), at least one reference
multiplicity was exhausted. Choose \(U,V\) with
\(UA_0+VB_0=1\bmod J\). Then the monic factor coefficients satisfy
\[
                  H-J=\operatorname{rem}_J(Uf+Vg)
                            +O((\|f\|+\|g\|)^2).              \tag{5}
\]
Indeed write \(F=H A_1,\ G=H B_1\). In the quotient modulo \(J\),
the right-hand remainder before the error term is
\((H-J)(U A_1+V B_1)\). The second factor tends to one and is
invertible there, giving \(H-J=O(\|f\|+\|g\|)\). Monic division
now bounds \(A_1-A_0,\ B_1-B_0\) at that same rate, proving (5).

The nonzero constant of \(J\) makes the projective factor shape a
rational coefficient function regular near \(J\). Thus (5)
gives its linear response and an explicit quadratic remainder, even
at a multiple reference root. If check 1 holds for all possible
degree-\(d\) divisors, the rank argument already makes \(x,y,K=O(Z)\),
whether or not (3) excludes the reference. In the whole-block case
this makes the real shape interval \(O(Z)\) wide for comparison with
the [grouped residue reconstruction](MULTIVARIATE_GROUPED_SHAPE_CONGRUENCES.md).

## Exact repeated-reference algebra

For a synthetic reference at orders \(6,7\), take
\[
          h_1=z^2+z+1,\qquad h_2=z^2+z+2,\qquad
          F_*=h_1^2h_2,\qquad G_*=(z+1)F_* .
\]
There are two real monic degree-four divisors of its gcd, namely \(h_1^2\) and \(h_1h_2\). The latter takes a proper part of the repeated \(h_1\) block. In both quotients the first three normal directions have rank three. A row annihilating them, in increasing coefficient order, is
\[
 (493/130,\ 726/65,\ -131/10,\ 1)
        \quad\text{or}\quad
 (698/1061,\ 363/1061,\ -1303/1061,\ 1),
\]
respectively. Their evaluations on \(([\Omega_Z],[\Omega_{B_4}])\) are
\[
     (-121/1820,\ -498713/4815720)
       \quad\text{and}\quad
     (-121/12732,\ -2048203/78607368).
\]
Both rows give a strict sign in (3) for every nonnegative \(u\). The [exact checker](verify_multivariate_sylvester_contact_normal_polytope.py) verifies these values, the two ranks, and the multiplication-rank identity at repeated roots. This is a test of the repeated-reference algebra, not a claim of an ordinary factorial core with this Gaussian reference.
