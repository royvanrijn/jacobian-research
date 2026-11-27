# Fixed-level unweighted Frey monodromy and abc

For a primitive positive triple \(a+b=c\) put \(R=\operatorname{rad}(abc)\),
\(A=\log R\), and \(h=\log c\). Fix any finite set of primes and any
required multiplicities in that set. Full abc is equivalent to abc on
primitive positive triples
\[
                  a=3x^2,\qquad b=y^2,\qquad c=z^3
\tag{FL1}
\]
with those prime powers all in \(a\), and with \(32\mid a\). The Frey
curve \(E:y^2=x(x-a)(x+b)\) is semistable and has split multiplicative
reduction at every prime dividing \(a\). In particular one can force
\(3,5,7,11\mid a\). On that family the sharp bound is equivalent to an
unweighted primitive Hecke norm bound at the single pivot \(5\), with
Shimura levels \(110\) and \(770\).

## Moving prescribed prime powers into one letter

Take a primitive positive input \(a_0+b_0=c_0\), with \(a_0\ne b_0\).
Let \(K\) be divisible by \(16\), \(3\), and all prescribed prime powers.
Choose a fixed even \(n\ge4\), divisible by \(\varphi(K)\), with
\(n\ge\max_{p\mid K}v_p(K)\). Set \(\kappa=1\) if \(c_0\) is even and
\(\kappa=0\) otherwise, and use Ellenberg's operation in positive
coordinates:
\[
 U=(|a_0-b_0|/2^\kappa)^n,\qquad
 W=(c_0/2^\kappa)^n,\qquad V=W-U .
\]
The two bases are coprime. Their gcd before dividing by \(2^\kappa\)
divides \(2\), and when \(c_0\) is even the halved bases have opposite
parity because their sum is \(a_0\) or \(b_0\), both odd. Thus
\(U+V=W\) is positive and primitive. If a prime \(p\mid K\) divides
either base, then \(p^{v_p(K)}\) divides its \(n\)th power. If both
bases are units, their \(n\)th powers agree modulo \(p^{v_p(K)}\).
Consequently \(K\mid P=UVW\). Factoring the even-power difference by
\(c_0^2-(a_0-b_0)^2=4a_0b_0\) gives
\[
 \log W=n\log c_0+O_n(1),\qquad
 \operatorname{rad}(P)\ll_n \operatorname{rad}(a_0b_0c_0)c_0^{n-1}.
\tag{FL2}
\]
Indeed the remaining homogeneous factor has degree \(n-2\) and
\(\operatorname{rad}(U)\le |a_0-b_0|\).
[Ellenberg, Theorem 3 and Lemmas 4–5](sources/ellenberg-congruence-abc.pdf).

Put \(L=U^2+UV+V^2\) and \(J=(U-V)(2U+V)(U+2V)\). Expansion gives
\[
                         27P^2+J^2=4L^3 .
\]
This is the homogeneous identity for the inverse Legendre \(j\)-map
\(1728/j=27t^2(1-t)^2/[4(1-t+t^2)^3]\), with \(t=U/W\).
Its numerator vanishes at each of the three boundary points.
Since a prime dividing \(P\) cannot divide \(L\), and \(3\mid P\),
no odd prime divides both \(27P^2\) and \(4L^3\). At \(2\), \(L\) is odd,
and \(J\) is even by the displayed identity. The common divisor of
the three terms is exactly \(4\). We obtain
\[
                 (a,b,c)=\big(3(3P/2)^2,(J/2)^2,L^3\big).
\tag{FL3}
\]
Here \(J\ne0\) because \(U/W\) is an \(n\)th power of a rational
number and cannot be \(1/2\). The sole positive primitive input
with \(a_0=b_0\) is \((1,1,2)\).
For odd \(p\mid K\) one has \(v_p(a)=2v_p(P)+3\,\mathbf 1_{p=3}\).
At \(2\) one has \(v_2(a)=2v_2(P)-2\ge6\).
Thus all the required prime powers divide the designated letter.

Write \(h_0=\log c_0\) and \(A_0=\log\operatorname{rad}(a_0b_0c_0)\).
Since \(3W^2/4\le L\le W^2\), and
\(\operatorname{rad}(abc)\le \operatorname{rad}(P)\,|J|L\) with
\(|J|\le4W^3\), (FL2) yields
\[
              h=6nh_0+O_n(1),\qquad
              A\le A_0+(6n-1)h_0+O_n(1).
\tag{FL4}
\]
More precisely, for \(d=6n\), \(\varepsilon>0\), and
\(\varepsilon'=\varepsilon/[d+(d-1)\varepsilon]\),
\[
 h-(1+\varepsilon')A
 \ \ge\
 \frac{d}{d+(d-1)\varepsilon}\bigl(h_0-(1+\varepsilon)A_0\bigr)
       -O_{n,\varepsilon}(1).
\tag{FL5}
\]
A uniform bound for the left side therefore gives abc on all
primitive positive inputs. The reverse implication is restriction.
The use of \(O(1)\) constants is equivalent to the finiteness
formulation for all positive exponent slacks.

At any odd prime dividing \(abc\), reduction of the Frey cubic has
a double root and a distinct simple root. At a prime dividing \(a\)
the tangent cone at the node is \(y^2=bx^2\), which splits because
\(b\) is an odd-prime unit square there. At \(2\), substitute
\(x=4X,\ y=8Y+4X\) to get the integral model
\[
 Y^2+XY=X^3+\frac{b-a-1}{4}X^2-\frac{ab}{16}X .
\]
Here \(b\) is an odd square and \(32\mid a\), so its reduction is
\(Y^2+XY=X^3\). The node has rational tangents \(Y=0\) and
\(Y+X=0\). At primes outside \(abc\) the original model has good
reduction. Hence \(N_E=R\), and all primes on \(a\) are split
multiplicative. On the whole family (FL1) with \(32\mid a\), a
prime \(q\mid c\) is odd and different from \(3\), and
\(y^2\equiv-3x^2\bmod q\) implies \(q\equiv1\bmod3\).
At an odd \(q\mid b\) the tangent square class is \(-a\), so the
reduction is split precisely when \((-3/q)=1\). At \(q\mid c\) the
tangent square class is \(a\), which equals the class of \(-1\)
because \(a\equiv-b\bmod q\). Thus the split primes there are
precisely \(q\equiv1\bmod12\).

## Two fixed levels with unit monodromy weights

Use the family (FL1) with \(32\cdot3\cdot5\cdot7\cdot11\mid a\).
For the map above one can take \(K=16\cdot3\cdot5\cdot7\cdot11\)
and \(n=60\): \(60\) is divisible by the unit group exponents
modulo \(16,3,5,7,11\) separately, which is all the preceding
prime-by-prime argument requires. Let \(k\) be the number of odd
prime divisors of \(R\), and set
\[
       M=\begin{cases}110& k\ \text{even},\\770&k\ \text{odd},\end{cases}
       \qquad D=R/M,\qquad p=5 .
\]
Then \(D\) is odd and is a product of evenly many distinct
multiplicative primes. The factorization is admissible and \(M\)
contains at least two odd primes. The remaining variable
quaternion discriminant in the definite description is \(5D\).
The Eichler level \(M/5\) is \(22\) or \(154\).

The integral monodromy at \(5\) is isometric to the degree-zero
ideal class module of an Eichler order of this level and discriminant.
Its weights are \(|\mathcal O_i^\times|/2\), with \(\mathcal O_i\) the
left orders of locally principal right ideals. Each local order at
\(11\) is conjugate to the integral matrices upper triangular modulo
\(11\). A noncentral unit in \(\mathcal O_i\) has norm \(1\), and its
trace is an integer strictly between \(-2\) and \(2\) by definiteness.
Its reduced characteristic polynomial would then be
\(T^2+1\) or \(T^2\pm T+1\). These polynomials are irreducible
modulo \(11\): their discriminants \(-4\) and \(-3\) are nonsquares.
But an upper triangular matrix has split characteristic polynomial
over \(\mathbb F_{11}\). Thus every unit is \(1\) or \(-1\), and every
weight is \(1\).
[Kohel, §3.2 and Theorem 4.3](sources/kohel-hecke-module-structure.pdf).

Let \(s\) be the number of classes and \(v\in\mathbb Z^s\) the primitive
generator up to sign of the rational Hecke line of \(E\). In these
coordinates
\[
                 \sum_i v_i=0,\qquad
                 \gcd_i v_i=1,\qquad
                 H_{D,M,5}=\sum_i v_i^2 .
\tag{FL6}
\]
The lattice has Gram matrix \(I_{s-1}+\mathbf 1\mathbf 1^{\mathsf T}\)
in the basis \(e_i-e_s\). Its dual modulo the lattice is cyclic of
order \(s\): adjoining
\((1/s)\mathbf 1-e_s\) gives a dual vector of order \(s\), and the
determinant is \(s\). Hence the geometric component group at \(5\)
is cyclic of that order. The pairing gcd on the line is
\(\gcd_{i,j}(v_i-v_j)\) and divides \(s\), since summing the coordinate
differences with one index fixed gives \(s v_i\) and the vector is
primitive. The Eichler mass formula gives
\[
           s=\frac{\varphi(5D)\psi(M/5)}{12}
             =\begin{cases}
                12\,\varphi(D) & M=110,\\
                96\,\varphi(D) & M=770,
               \end{cases}
\tag{FL7}
\]
where \(\psi(m)=m\prod_{q\mid m}(1+q^{-1})\).
[Voight, Theorem 25.3.18](https://doi.org/10.1007/978-3-030-56694-4_25).

## The remaining uniform estimate

By the [monodromy height comparison](ABC_FREY_MONODROMY_NORMS.md#two-vertices-suffice-for-the-norm-bound),
equation (MN4), uniformly over the family just defined,
\[
   \log H_{D,M,5}=h+A+O\!\left(\frac{h}{\log(h+2)}\right).
\]
Consequently full abc is equivalent to
\[
   \boxed{\quad \sum_i v_i^2\ll_\varepsilon R^{\,2+\varepsilon}
                \quad\text{for every }\varepsilon>0.\quad}
\tag{FL8}
\]
In the reverse direction absorb the height error with exponent
slack in the family, then apply (FL5). It suffices to bound only
the explicitly constructed images (FL3) with the single choice
\(n=60\). Conversely a sequence in that image with
\(H_{D,M,5}/R^{2+\varepsilon}\) unbounded for a fixed
\(\varepsilon>0\) disproves abc.

The same setup implies \(s\ll R\) by (FL7). Thus a uniform bound on
the largest primitive coordinate of size \(O_\eta(R^{1/2+\eta})\)
for each \(\eta>0\) is sufficient. A bound for the eigenvector in
unit norm alone leaves the primitive integral normalization
uncontrolled. The open step in (FL8) is to use the arithmetic of the
Frey Hecke line to control that normalization as \(D\) varies.
