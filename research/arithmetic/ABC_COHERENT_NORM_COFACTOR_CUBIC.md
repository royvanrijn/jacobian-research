# The binary cubic on a coherent prime-path lattice

A coherent triple determines an integral binary cubic of discriminant \(27d\delta^2\) on its cofactor lattice. Its primitive Hessian is the norm form of the ideal with class \(2c_2-c_1-c_3\), where \(c_i\) are the oriented prime classes. For two independent edge differences this is the fourth reflected line, alongside the three edge lines.

## Exact cubic and quadratic forms

Use the [signed reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients) and write \((v_{12},v_{23},v_{13})=\delta(x,y,z)\), \(E=d\delta^2\), and \(D=-E\). Assume \(d>3\) is squarefree and the three primes are distinct, odd and unramified. The order \(\mathcal O_D\) of discriminant \(D\) contains the witnesses \(\alpha_{ij}=(u_{ij}+v_{ij}\sqrt{-d})/2\). Indeed \(-E\) is a discriminant by the integral parity conditions, and the conductor is prime to the three primes since \((E,p_1p_2p_3)=1\). Write \(\mathfrak p_i\) also for the corresponding invertible primes of this order.

Let \(p=(p_1,p_2,p_3)\), \(Q=|xyz|\), and use the exact coefficient moduli:
\[
 \Lambda=\{w\in\mathbb Z^3:
     p_2 w_1=p_1 w_2\bmod |x|,\quad
     p_2 w_3=p_3 w_2\bmod |y|,\quad
     p_3 w_1=p_1 w_3\bmod |z|\}.
\]
The incident units and pairwise coprimality give \([\mathbb Z^3:\Lambda]=Q\). Compared with the [original full-modulus lattice](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#full-moduli-including-the-factor-2), this uses the same odd parts and refines the conditions at \(2\) by a bounded index. Write \(L=p\times\Lambda\) and use any integral basis of \(L\) as the two coordinates of the forms below.

**Proposition (Cofactor cubic).** Let \(I=\overline{\mathfrak p}_1\mathfrak p_2^{\,2}\overline{\mathfrak p}_3\), \(n_0=NI=p_1p_2^2p_3\), and \(\Theta=\alpha_{12}\overline{\alpha}_{23}\), so \(I^3=(\Theta)\) in \(\mathcal O_D\). Then \(n=(n_1,n_2,n_3)\mapsto\gamma\) below is an isomorphism of \(\mathbb Z\)-lattices \(L\to I\):
\[
 \gamma=\frac{p_1p_2^2(n_1/y)+\alpha_{12}(n_2/z)}{x}.
                                                                       \tag{1}
\]
The cofactor cubic and its Hessian are
\[
 \begin{split}
 \mathcal B(n)&=-\frac1{xyz}
     \left(\frac{n_1^3}{y^2}+\frac{n_2^3}{z^2}+\frac{n_3^3}{x^2}\right)
       =\frac{\gamma^3/\Theta-\overline{\gamma}^{\,3}/\overline{\Theta}}{\sqrt D},\
 H_{\mathcal B}(n)&=9\,\frac{N\gamma}{n_0}.
 \end{split}                                                          \tag{2}
\]
The cubic has integer coefficients, with the two mixed coefficients divisible by \(3\). The quadratic form \(H_{\mathcal B}/9\) is primitive of discriminant \(D\). In particular \(\operatorname{disc}(\mathcal B)=27E\).

To prove (1), put \(a_1=n_3/x,\ a_2=n_1/y,\ a_3=n_2/z\). These are integers by the congruences defining \(\Lambda\). They satisfy \(p_3x a_1+p_1y a_2+p_2z a_3=0\). Since \(p\) is primitive, extend \(p\) to a basis \(p,w,w'\) of \(\Lambda\). The two basis vectors \(p\times w,p\times w'\) of \(L\) have cross product \(\pm Qp\). Thus the determinant on coordinates \((a_2,a_3)\) is \(\pm xp_3\). The identity \(y\alpha_{12}+x\overline{\alpha}_{23}=p_2^3z\) gives
\[
                 y\gamma=-p_2^2p_3 a_1-\overline{\alpha}_{23} a_3.
\]
Both \(x\gamma\) and \(y\gamma\) are in \(\mathcal O_D\) and \((x,y)=1\), so \(\gamma\in\mathcal O_D\). The first expression puts it in \(\overline{\mathfrak p}_1\mathfrak p_2^2\), since \(x\) is a unit at both primes. The second puts it in \(\overline{\mathfrak p}_3\), since \(y\) is a unit there. The unscaled map to \(x\gamma\) has determinant of absolute value \(p_1p_2^2 |x|\cdot|xp_3|\) in the integral basis of \(\mathcal O_D\). Division by \(x\) gives index \(n_0\), the norm of \(I\), proving equality. The generators \(\alpha_{12},\alpha_{23}\) have their expected prime ideal factorizations in \(\mathcal O_D\), because their norms are prime to its conductor.

For the first identity of (2), use \(A=p_1^3 y^2,\ B=p_2^3 z^2,\ C=p_3^3 x^2\), \(U=p_1n_1,\ V=p_2n_2,\ W=p_3n_3\), and \(h=xyz\sqrt D\). Then \(U+V+W=0\). Set \(\rho=(A+B-C+h)/(2B)\); the trace and norm identities give
\[
 \rho\overline\rho=A/B,\quad
 (1-\rho)(1-\overline\rho)=C/B,\quad
 \rho-\overline\rho=h/B,\quad
 \frac{\gamma^3}{\Theta}=\frac{(U+\rho V)^3}{Q^2\rho(1-\rho)}.
\]
In the last fraction, the partial fractions in \(\rho\) of \((U+\rho V)^3/(\rho(1-\rho))\) have terms \(U^3/\rho\), \((U+V)^3/(1-\rho)\) and polynomial part \(-\rho V^3-3UV^2-V^3\). Subtracting the conjugate proves (2). In the basis of \(I\), coefficients before division by \(\sqrt D\) come from products of three basis elements divided by \(\Theta\), hence from \(\mathcal O_D\). Extracting the signed imaginary coefficient gives integers, with the binomial factor \(3\) on the mixed coefficients. The coefficient extraction is the ideal-to-form construction of [Bhargava, Theorem 13 and equations (21)–(24)](https://annals.math.princeton.edu/wp-content/uploads/annals-v159-n1-p03.pdf). Formula (1) identifies the ideal input of that construction directly on the prime-path lattice.

For the Hessian, use \(H_f=-\tfrac14\det(\partial_i\partial_j f)\). In the linear coordinates \((\gamma,\overline\gamma)\) the cubic is diagonal, and the basis determinant in those coordinates has square \(D n_0^2\). The diagonal Hessian and \(N\Theta=n_0^3\) give \(H_{\mathcal B}=9N\gamma/n_0\). Its normalized discriminant is \(D\) by the same basis determinant. To see primitivity at a rational prime \(\ell\), choose the omitted coefficient and opposite prime to be units at \(\ell\), relabeling the triple if necessary. At least one choice is possible, since the coefficients are pairwise coprime and the primes are distinct. For the choice \((x,p_3)\), the determinant \(\pm xp_3\) is a unit and in \((a_2,a_3)\) the norm form is
\[
 \frac{p_1p_2^2 a_2^2+u_{12}a_2 a_3+p_1^2p_2 a_3^2}{p_3x^2}.
\]
An outer coefficient is a unit if \(\ell\nmid p_1p_2\); the middle one is a unit otherwise by the edge norm and its incident units. Permuting the primes leaves the sum of the three terms in \(\mathcal B\) unchanged up to signs and coordinate permutation, so the same reasoning covers all \(\ell\). Now \(\operatorname{disc} H_f=-3\operatorname{disc} f\) gives the cubic discriminant. The Hessian describes the corresponding cubic order's shape by [Bhargava–Shnidman, Theorem 9 and Propositions 11–12](https://arxiv.org/pdf/1206.4746#page=7).

There are also three distinguished primitive vectors of \(L\), one on each coordinate hyperplane \(n_i=0\). The image in \(I\) of these three lines, given by (1), has primitive generators up to sign
\[
   p_3\alpha_{12}\ (n_1=0),\qquad n_0\ (n_2=0),\qquad
   p_1\overline{\alpha}_{23}\ (n_3=0).
\]
Indeed the scaled lattice of the \(a_i\) is the entire integer plane \(p_3x a_1+p_1y a_2+p_2z a_3=0\): membership follows from \(p\times\mathbb Z^3=p^\perp\cap\mathbb Z^3\), and the incident units make the coordinate divisibilities equivalent to the defining congruences on a preimage. The three coefficients of the plane are pairwise coprime by the incident units. With \(n=p_1p_2p_3\), the three primitive values of \(H_{\mathcal B}/9\) on these lines are \(np_1,np_2,np_3\), respectively. The norm form therefore retains three marked representations with the same product of primes.

## An exact example in discriminant \(-4027\)

For \((p_1,p_2,p_3)=(13,17,29)\), take prime cube generators with doubled coefficient pairs \((69,1),(125,1),(182,4)\) in \(\mathbb Q(\sqrt{-4027})\). Their three coherent edges have pairs \((6326,-28),(19429,159),(14333,47)\) in order \(12,23,13\), so \(\delta=1\). In the scaled \(a\)-coordinates of (1), use the basis \((-4,5,-17),(-103,-71,79)\); its cross product is the primitive normal of the plane. The cubic and its primitive Hessian are
\[
    \mathcal B(S,T)=-S^3+12S^2T+3ST^2-15T^3,\qquad
    H_{\mathcal B}(S,T)/9=17S^2-11ST+61T^2.
\]
Reduction of \(\mathcal B(S,1)\) modulo \(2\) is \(S^3+S+1\), which has no root. The reflected trinomials on the three edge classes are irreducible modulo \(5,2,2\), respectively. Thus this example has all four lines nonzero, with the fourth supplied by the cofactor cubic.

## Irreducible path selection

**Corollary (Irreducibility and the cofactor height).** In the original tail construction, the same lower bound (2) for selected paths holds with \(\mathcal B\) irreducible. On each retained path every nonzero cofactor \(k\in L\) satisfies
\[
       \|k\|\gg (|xyz|\min(|x|,|y|,|z|)^2)^{1/3}.
                                                                       \tag{3}
\]
For a balanced occupied rank-one record at \(V=Z^{5/4}\), this gives \(\lambda_2\gg ZV^{1/15}\), where \(\lambda_2\) is the second minimum of the original full-modulus lattice.

Indeed \(\mathcal B\) has a rational root exactly when \(\Theta\) is a cube in \(K\). If \(\gamma^3/\Theta\) is rational and nonzero, all valuations of this quotient in the maximal order are multiples of \(3\). Since the ramification indices over \(\mathbb Q\) are \(1\) or \(2\), it is a rational cube. Conversely a cube root of \(\Theta\) gives a rational line in \(I\otimes\mathbb Q\). Unique ideal factorization in the maximal order and the units \(\{\pm1\}\) then show that reducibility is equivalent to \(2c_2-c_1-c_3=0\) in \(\mathrm{Cl}(K)\). This test takes place in the maximal order; the Hessian in (2) retains the conductor of \(\mathcal O_D\).

For an irreducible form, its splitting field is the reflected field for \(\Theta\), from [Kishi's trinomial](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#cubic-fields-detect-the-first-color) \(t^3-3n_0t-\operatorname{Tr}(\Theta)\). Over \(K(\zeta_3)\), the three ratios \(\gamma/\overline\gamma\) at the roots of \(\mathcal B\) are the cube roots of \(\Theta/\overline\Theta=\Theta^2/n_0^3\). These generate the same Kummer extension as a cube root of \(\Theta\), and both splitting fields are real. When the edge differences \(c_2-c_1,c_3-c_2\) are independent they give two lines, their sum gives the third, and their difference gives the line of this cubic.

For the path count in \([X,2X)\), take \(X\) large and use \(Z\le X^{1/3}\). No two distinct selected primes can have the same full ideal class. Otherwise the product \(\overline{\mathfrak p}_i\mathfrak p_j\) has a nonrational generator of norm \(p_i p_j<4Z^2\), whereas its nonzero imaginary coefficient forces norm at least \(d/4\). Thus the classes in each quotient coset are distinct points of an affine \(\mathbb F_3\)-space. An unordered pair of classes determines at most one third class on its line, so the number of collinear triples in a coset of size \(s\) is at most \(\binom{s}{2}/3\). The original four-point argument gives at least \(s(s-1)(s-2)/24\) unequal-color adjacent pairs. Discarding all adjacent pairs from collinear triples costs at most \(\binom{s}{2}\). The surplus is \(\gg s^3\) for \(s\ge28\), and raising the minimum coset size to \(28\) has the same negligible loss in the witness selection. The remaining paths have independent class differences, hence irreducible cubics.

For such a path, \(|\mathcal B(k)|\ge1\) for \(0\ne k\in L\); the first expression in (2) gives (3). If an occupied rank-one lattice has shortest direction \(p\), then a second independent vector of length \(O(\lambda_2)\) in the refined lattice gives a nonzero \(p\)-cofactor of norm \(O(Z\lambda_2)\). The refinement at \(2\) changes successive minima by at most a fixed factor. When \(|x|,|y|,|z|\asymp V\), (3) gives \(Z\lambda_2\gg V^{5/3}\), which yields the claimed power at the original crossover. The divisor bound for rank-two counts uses individual values of the cofactor cubic. The present identification gives a norm form and its class on the whole cofactor lattice, including the long second directions, where a uniform occupied-record bound is still needed.
