# Joint native residue and denominator bounds for common factors

Local quadratic residue compatibility can hold simultaneously at every prime of the linked odd-corner moduli. A global common quadratic has one rational angular invariant across all those reductions. An exact divisor bound on its denominator turns the joint residue constraints into one arithmetic progression in a bounded real interval. The interval-length criterion is linear in the divisor bound. It applies to original factorial rows with unrestricted counts, gaps, and heights. It also extends to higher common-factor degrees when modular data determine the factor at the tested primes.

## Obtaining the common quadratic at native primes

Let \(P_n,P_N\in\mathbb Z[z]\) be the primitive parents of a given original binomial, with \(N=n+1\). Their coefficient scales can be chosen independently. The even parent is strictly positive on the real line. Every nonconstant rational common factor therefore has even degree, including after a common reversal of the variables or any common nonzero rational argument scaling.

Suppose that at a prime \(p\) one of the parents has its leading coefficient a unit, the two reductions have gcd of degree two, and that gcd has unit constant coefficient. Any primitive rational common factor preserves its degree modulo \(p\), by Gauss's lemma and the unit leading coefficient of that parent. Thus if a common factor exists, its primitive representative \(H=A z^2+Bz+C\) has degree two and both \(A,C\) are units at \(p\). Its monic reduction is the unique monic gcd there. This applies even if that modular gcd is inseparable.

The linked corner anchors supply the same information at full depth. At one selected \(p^e\), assume the marked parent reduces modulo \(p^e\) to a quadratic with leading and constant coefficients units and the unmarked parent has leading coefficient a unit. The scalar remainder must vanish for any common factor. At \(p\), the marked quadratic contains the reduction of the entire common factor, so the degree is two and \(A,C\) are units. Over \(\mathbb Z/p^e\mathbb Z\), divide the marked parent by the monic reduction of \(H\), which has degree two because \(A\) is a unit. This forces the normalized marked quadratic itself to equal the monic reduction of \(H\) modulo \(p^e\). The argument uses monic division and is valid at an inseparable residual quadratic as well. [Linked corner anchors](MULTIVARIATE_ODD_CORNER_LINKED_QUADRATIC_ANCHORS.md)

At a general prime power one may also certify a unique monic local factor of degree two by supplying a monic generator of the ideal of both reduced polynomial arrays. Given a degree-two modular gcd with unit leading and constant coefficients, it is enough to show that the monic local degree-two factor divides both parents and is an integral polynomial combination of them modulo \(p^e\). The reduction of a global degree-two common factor then divides the generator and equals it after monic normalization.

## One rational shape across the good moduli

Write
\[
            \xi=\frac{B^2}{AC}=\frac{U}{D},\qquad
            \gcd(U,D)=1,\qquad 0\le U<4D .
\]
The positivity of the even parent gives \(B^2-4AC<0\) and \(AC>0\). Equivalently \(\xi=4\cos^2(\arg\alpha)\) for either root \(\alpha\) of \(H\). At a good modulus \(p^e\) with recovered monic factor \(z^2+b_pz+c_p\), one has
\[
                  \xi\equiv b_p^2/c_p\pmod {p^e}.                \tag{1}
\]
All fractions in this congruence have denominators prime to \(p\). Combine the congruences from any list of good prime powers into \(\xi\equiv r\bmod M\), where \(M\) is their least common multiple. Inconsistent congruences are already an exclusion. Use at least one good prime to force the global common degree to be two.

Let \(D_0\) be a positive integer divisible by every allowable \(D\). For example, for each prime take the maximum spread of two allowed common Newton root valuations and multiply the corresponding prime powers, or use the [oriented factorial block gcd bound](MULTIVARIATE_ORIENTED_BLOCK_DENOMINATOR_GCD.md). Two distinct root valuations of a rational quadratic are integers since the unequal-slope quadratic splits over \(\mathbb Q_p\). In the equal-slope case the spread is zero. The exact denominator identity \(v_p(D)=|s-t|\) for the two root valuations follows from the monic trace and norm, or from the [projective denominator formula](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#exact-projective-denominators). The full Newton branch compatibility tests can be used to reduce the divisor budget.

Remove from \(D_0\) every prime power supported on primes dividing \(M\), and call the result \(D_1\). Since \(A,C\) are units at every good prime, \(D\mid D_1\). Consequently the integer \(T=D_1\xi\) satisfies
\[
                          T\equiv D_1 r\pmod M .                 \tag{2}
\]
If every possible \(\xi\) lies in a certified real interval \(I\), failure of this arithmetic progression to intersect \(D_1 I\) excludes a common root. This test is exact using a rational interval enclosing the required angle ranges. It always applies with \(I=[0,4)\). An interval of length less than \(M/D_1\) contains at most one permitted \(\xi\); for a union, apply the test in each interval. Even if multiple integers remain, (2) leaves only the listed rational shapes. Only removal of the primes dividing \(M\) from a specified divisor budget is required, which can also be computed by repeated gcds.

For comparison, rational reconstruction with a denominator upper bound \(D_{\max}\) alone gives uniqueness when \(D_{\max}^2 |I|<M\): the cross determinant of two distinct candidates is a nonzero multiple of \(M\) and has absolute value at most \(D_{\max}^2 |I|\). Formula (2) uses the divisibility information from the actual factorial arrays and gives the linear condition \(D_1|I|<M\). [Wang, Guy and Davenport, rational reconstruction](https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA13/GuyDavenWang.pdf)

An isolated candidate \(\xi\in(0,4)\) leaves only a rational scale. Write \(\rho=B/A\ne0\). The monic candidate is \(z^2+\rho z+\rho^2/\xi\). Reduce \(P_j(\rho x)\), for \(j=n,N\), modulo the fixed \(x^2+x+1/\xi\). Its two coefficients are rational polynomials in \(\rho\). The shape occurs as a common factor exactly when all four scalar polynomials have a common nonzero rational root. This is direct polynomial arithmetic and lets one use bounds on the physical root scale before division. If \(\xi=0\), write the monic candidate as \(z^2+c\) with \(c>0\) rational. Evaluate the even and odd coefficient polynomials of each parent at \(-c\), obtaining the corresponding four scalar conditions.

## Degree greater than two

Suppose other modular or algebraic constraints force any common factor to have degree \(d\), and the monic degree-\(d\) factor is uniquely specified modulo each tested \(p^e\), with unit constant. For example a good modular gcd irreducible of degree \(d\) with unit constant, together with a parent's unit leading coefficient there, forces the degree. A higher-depth local ideal certificate is then used in the same way as for the quadratic.

For the monic global candidate \(H\) with roots \(\alpha_i\), put \(c=\prod_i\alpha_i\) and
\[
             B_H(Y)=\prod_{i=1}^d(Y-\alpha_i^d/c)
                     =c^{-d}\operatorname{Res}_z(H(z),cY-z^d).
\]
This is the projective shape polynomial, with rational coefficients. Its reduction is explicitly determined by the local factor because the formula only uses integer polynomial operations and inversion of the unit \(c\). Let \(L_0\) be one of the available exact divisor budgets for the common denominator of the coefficients of \(B_H\), for instance the bound from the [oriented factorial blocks](MULTIVARIATE_ORIENTED_BLOCK_DENOMINATOR_GCD.md). Strip the prime powers of the joint modulus \(M\) from \(L_0\), obtaining \(L_1\). For every coefficient \(b_k\) of \(B_H\), its integer lift \(L_1b_k\) lies in the progression specified by the joint reductions of \(L_1 B_H\). Any certified real interval for a coefficient applies independently, with the same linear uniqueness criterion \(L_1|I_k|<M\). Rough intervals follow already from root-radius bounds, since \(|\alpha_i^d/c|\le(r_{\max}/r_{\min})^d\).

## An original profile whose two modular shape tests pass separately

Take three rows \((0,1)\) and three \((1,0)\). At orders \(5,6\), dividing by the endpoint coefficient makes the coefficient at index \(k\) equal to \(\binom jk^{-2}\). The primitive integral arrays in increasing coefficient order are
\[
                 (100,4,1,1,4,100),\qquad
                 (3600,100,16,9,16,100,3600).
\]
Shared Newton root-valuation widths permit degree two at every prime. The only possible nonzero spread is two at \(p=2\): the shared slopes are \(1,-1\) of width two each there, and only zero is shared at \(3,5\). At all other primes the endpoints are units. Thus the width-only denominator budget is \(D_0=4\).

Modulo 7 and 31 the monic gcds are \(z^2+5z+1\) and \(z^2+3z+1\). At 7 this gives \(\xi\equiv4\), and at 31 it gives \(\xi\equiv9\). The modulus 7 alone permits \(T=2,9\) in \(0\le T<16\), and the modulus 31 alone permits \(T=5\). Both local gcds have degree two and unit endpoints. Together they give \(r=102\bmod217\) and force \(T\equiv191\bmod217\), excluding a common quadratic by (2). The example isolates the simultaneous shape obstruction relative to these two modular reductions and the Newton width budget. Other modular factor certificates may also exclude this profile. The [exact checker](verify_multivariate_joint_native_shape_reconstruction.py) recomputes the arrays from the ordinary factorial formula, checks the full Newton widths at the endpoint primes, reconstructs the two modular gcds, and verifies the separate and joint residue progressions.

The criterion leaves a genuine scope condition for unrestricted closure. The original valuation budget \(D_1\) can grow with the row multiplicities, and a linked anchor product and a given angular interval need not yet isolate the rational shape.
