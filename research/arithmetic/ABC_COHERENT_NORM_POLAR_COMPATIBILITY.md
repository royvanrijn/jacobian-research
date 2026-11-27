# Pair compatibility of marked norm vectors

For a projective binary cubic with imaginary quadratic Hessian, two marked norm vectors determine a third by polarization. An exact determinant gcd decides whether the primitive third vector has the required norm and completes a coherent prime path. A single norm vector can also be tested by one equality of ideal cube classes. A separate count bounds the sector in which the opposite primes divide two or three normalized coefficients.

## Completion from two norm vectors

Let \(F(S,T)=aS^3+3bS^2T+3cST^2+dT^3\) come from an invertible ideal \(I\) of an imaginary quadratic order \(\mathcal O_D\) and a generator \(\Theta\) of \(I^3\). Use an oriented integral basis to identify \(w\in\mathbb Z^2\) with \(\gamma_w\in I\), and set \(n_0=NI\). Extend to rational coordinates when needed. The form and its normalized Hessian are
\[
 F(w)=\frac{\gamma_w^3/\Theta-\overline{\gamma_w}^{\,3}/\overline\Theta}{\sqrt D},
 \qquad H(w)=\frac{N\gamma_w}{n_0}
       =(b^2-ac)S^2+(bc-ad)ST+(c^2-bd)T^2.
\]
This includes the [cofactor cubic](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#exact-cubic-and-quadratic-forms). The ideal description for projective integer-matrix cubics is given by [Bhargava–Varma, Theorem 9 and §2.2](https://arxiv.org/pdf/1401.5875).

Let \(T_F\) be the symmetric trilinear form with \(T_F(w,w,w)=F(w)\). In particular its four entries on basis triples are \(a,b,c,d\). Write \(K_F(x,y)\) for the integer vector determined by
\[
             \det(K_F(x,y),w)=T_F(x,y,w)\quad(w\in\mathbb Z^2).
\]
Thus if the row of the polar on the right is \((t_1,t_2)\), the vector on the left is \((t_2,-t_1)\). For a nonzero vector \(z\), let \(\operatorname{cont}(z)\) be the positive gcd of its two coordinates.

**Proposition (Polar completion and determinant test).** Let \(r_1,r_2,r_3\) be distinct odd primes split and prime to \(D\), and \(N=r_1r_2r_3\). Suppose \(x,y\) are primitive with \(H(x)=Nr_1\) and \(H(y)=Nr_2\). Set \(z=K_F(x,y)\), \(g=\operatorname{cont}(z)\). Then
\[
 g=\frac{N}{\gcd(N,\det(x,y))},\qquad H(z/g)=\frac{N^2r_1r_2}{g^2}.
 \tag{1}
\]
The marks \(x,y,z/g\) have the three required norm values and the coherent product of their representing ideals exactly when \(\gcd(N,\det(x,y))=r_3\). In that case they reconstruct the path. Its cofactor order is \(\mathcal O_D\) when the three edge coefficients relative to \(\sqrt D\) are primitive as a triple. The only primitive choices for the completing mark with norm \(Nr_3\) and the coherent product are \(\pm z/g\) when the field radicand is greater than three.

To prove this, polarize the ideal formula for \(F\). In the oriented basis, the determinant is
\[
 \det(v,w)=\frac{\overline{\gamma_v}\gamma_w-\gamma_v\overline{\gamma_w}}{n_0\sqrt D}.
\]
Therefore \(\gamma_z=n_0\overline{\gamma_x\gamma_y/\Theta}\). For \(J_w=(\gamma_w)I^{-1}\), this implies
\[
                         J_z=\overline{J_x}\,\overline{J_y}.          \tag{2}
\]
Indeed \((\overline\Theta)=\overline I^{\,3}\) and \((n_0)=I\overline I\). Taking norms gives \(H(z)=H(x)H(y)\).

For each \(p\mid N\), primitivity implies \(J_x,J_y\) are supported on one of the two primes above \(p\) each. The two exponents are \(2,1\), \(1,2\), or \(1,1\), so their minimum is one. Opposite orientations give exactly one factor of the rational integer \(p\) in \(J_z\); agreeing orientations give none. A rational divisor of \(J_z\) is exactly a divisor of the coordinates of \(z\) in the basis of \(I\). Also the reductions of \(x,y\) in \(I/pI\) are both nonzero, each lying on the line specified by its vanishing component at one prime over \(p\). Thus \(p\) divides \(\det(x,y)\) exactly when they agree in orientation. This proves the first formula of (1). Homogeneity of \(H\) gives the second.

The condition \(\gcd(N,\det(x,y))=r_3\) gives \(g=r_1r_2\) and \(H(z/g)=Nr_3\), and (2) then yields \(J_xJ_yJ_{z/g}=(N^2)\). The [three-mark inverse](ABC_COHERENT_NORM_MARKED_REPRESENTATIONS.md#exact-reconstruction-and-automatic-orientation) recovers the witnesses, with the stated condition on the normalized edge coefficients for the cofactor order. Conversely any completion by the coherent product must have \(J_{w}=(N^2)(J_xJ_y)^{-1}\). In the successful case this is \(J_{z/g}\), so their generators differ only by a unit. In an imaginary field with radicand greater than three the units are \(\{\pm1\}\). If the gcd test fails there can be no primitive completion with the coherent product, since at any differing prime the required orientations of \(J_x,J_y\) at that prime already fail.

For the [cofactor example in discriminant \(-4027\)](ABC_COHERENT_NORM_COFACTOR_CUBIC.md), the divided coefficients are \((-1,4,1,-15)\). Take \(x=(-103,4)\), \(y=(-71,-5)\) and \((r_1,r_2,r_3)=(29,13,17)\). The determinant is \(799=17\cdot47\), and the rotated polar is \((29783,6409)=29\cdot13(79,17)\), so the third primitive mark is \((79,17)\). For a control from the [identity-coset example](ABC_COHERENT_NORM_MARKED_REPRESENTATIONS.md#failure-of-automatic-orientation-in-the-identity-coset), use divided coefficients \((1,4,-3,-11)\), with \(H=19S^2-TS+53T^2\) and \((r_1,r_2,r_3)=(1009,13,17)\). The primitive vectors \((-3168,-835)\), \((-363,83)\) have norms \(Nr_1,Nr_2\), respectively. Their determinant has gcd \(1009\cdot17\) with \(N\), and their polar has content \(13\), ruling out their joint completion.

This test writes the compatibility condition as an incidence on two points of a single Hessian, with the third determined by them. In the rank-one count, vary the projective cofactor cubic and retain the recovered coefficient box, field window and lattice-rank test.

## Completion from one mark and an ideal class comparison

**Corollary (One marked norm and a ring-class condition).** With \(F,H,I\) as above, take a primitive vector \(w\) with \(H(w)=r_1^2r_2r_3\), with distinct odd primes split and prime to \(D\). Write \(J_w=\mathfrak a_1^2\mathfrak a_2\mathfrak a_3\) for the prime ideals selected by its local orientations. The two other primitive marks with norms \(Nr_2,Nr_3\) and coherent product exist exactly when
\[
                 [(\mathfrak a_1\mathfrak a_2)^3]=1
                    \quad\text{in }\operatorname{Pic}(\mathcal O_D). \tag{3}
\]
The three orientations and the edge data are then determined up to the choices of signs of the two generators if the field radicand is greater than three. The test works with \(\mathfrak a_3\) in place of \(\mathfrak a_2\). The recovered coefficient gcd still determines whether the cofactor order is precisely \(\mathcal O_D\).

In any completion the representing ideal of the squared mark at \(r_2\) must be \(J_2=\overline{\mathfrak a}_1\,\overline{\mathfrak a}_2^{\,2}\mathfrak a_3\). It comes from a vector of \(I\) exactly when \(I J_2\) is principal. But \(I J_w\) is principal and the ratio \(J_w J_2^{-1}\) has ideal class \([(\mathfrak a_1\mathfrak a_2)^3]\), since \([\overline{\mathfrak a}_i]=[\mathfrak a_i]^{-1}\). This gives (3). A generator of \(I J_2\) lies in \(I\) and is primitive by the single orientation at each prime of its ideal norm. Apply the determinant test to \(w\) and this generator. Their orientations differ at \(r_1,r_2\) and agree at \(r_3\), so the third mark is given by (1). The two alternatives for the class test are equivalent because \(J_w^3\) is principal, and \(J_w^3=(\mathfrak a_1\mathfrak a_2)^3(\mathfrak a_1\mathfrak a_3)^3\). Thus in the varying-order incidence count one can start with one vector, factor its norm to recover the three primes, and test one cube class in the actual order.

## Counting opposite-prime divisibility

Use the reordered opposite primes \(r=(p_3,p_1,p_2)\) and the absolute normalized coefficients \(q=(|x|,|y|,|z|)\) of the [complete-modulus congruences](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#exact-full-moduli-and-the-lattice-count). In this section the three \(q_i\) are pairwise coprime, \(q_i\) is prime to the incident primes \(r_j,r_k\), and \(r_j^\ell q_j^2\equiv r_k^\ell q_k^2\bmod q_i\). These properties hold also at \(\ell=3\) by [signed reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients).

**Proposition (Two or three opposite primes in the moduli).** Fix an odd prime \(\ell\). For \(r_i\in[Z,2Z)\), \(q_i\ll V\), the number of coherent witnesses, counted by the prime triple and signed normalized coefficients, with \(r_i\mid q_i\) for at least two indices is
\[
                         O_\varepsilon\big((ZV)^\varepsilon ZV^2(1+V/Z^2)\big).       \tag{4}
\]
With \(r_i\mid q_i\) for all three indices it is
\[
                         O(Z^3+V^2+V^3/Z^3).                         \tag{5}
\]
These bounds apply to either field sign and to arbitrary further restrictions on coherent paths.

For (4), fix the two indicated moduli \(q_i,q_j\). Each has \(O(\log(2V)/\log(2Z))\) possible marked primes. Fix \(r_k\) at cost \(O(Z)\). The congruences modulo \(q_i,q_j\) impose a unit square residue on \(q_k\) modulo their product, giving at most \(O_\varepsilon((ZV)^\varepsilon)\) root classes by the local unit square-root bound and CRT. Each class contains \(O(1+V/(q_iq_j))\) coefficients. The product of the marked moduli is at least \(Z^2\). Summing gives (4).

For (5), fix the three primes and set \(q_i=r_i s_i\). Modulo \(r_i\), the two incident terms give \(r_j^{\ell+2}s_j^2\equiv r_k^{\ell+2}s_k^2\). The incident coefficients are units. Thus all \(s\) lie in at most eight lattices with equations \(s_j\equiv c_i s_k\bmod r_i\) for fixed unit roots \(c_i\). Each lattice has determinant \(N\) by CRT. Every candidate in the box of side \(O(A)\), \(A=V/Z\), is primitive by pairwise coprimality of the \(q_i\). If the vectors in the box span rank at most one, there are \(O(1)\) primitive candidates. In rank two, the cross product of any two lattice vectors has its \(i\)-th coordinate divisible by \(r_i\), so the lattice in their spanning plane has covolume at least \(Z\). Its second minimum is \(O(A)\); the two-dimensional successive-minima count gives \(O(A^2/Z)\) points. In rank three the same count with all three minima \(O(A)\) gives \(O(A^3/N)\) points. Summing over the prime triples gives (5).

At \(V=Z^{5/4}\), (4) and (5) are \(O_\varepsilon(V^{14/5+\varepsilon})\) and \(O(V^{12/5})\). At \(V=Z^{7/4}\), they are \(O_\varepsilon(V^{18/7+\varepsilon})\) and \(O(V^2)\). For an unrestricted coherent example with all three opposite primes dividing their moduli, take \(r=(5,7,17)\) and \(q=(305,147,17)\) at \(\ell=3\). The signed reconstruction with all positive coefficient signs gives traces \((-1119,1087,-393)\) in edge order \(12,23,13\), with \(E=59\). These bounds leave open the occupied rank-one records with at most one opposite prime dividing its modulus.
