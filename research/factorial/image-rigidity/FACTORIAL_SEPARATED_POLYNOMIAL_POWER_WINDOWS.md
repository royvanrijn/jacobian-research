# Factorial windows for separated polynomial sums and their powers

Let
\[
             F=c+\sum_{i=1}^m q_i(X_i)
\]
over a characteristic-zero field, with arbitrary one-variable
polynomials \(q_i\). Absorb constant summands into \(c\), and assume
at least one remaining summand is nonconstant. Put
\(d_i=\deg q_i\), \(R=\prod_i d_i\), and let \(l\) be the number
of linear summands. For \(a,b\ge1\) write
\[
                       W(a,b)=a+b\,a(a+3)/2 .
\]

**Theorem.** For every \(k\ge1\), every positive starting order,
and every window of length
\[
                       \left\lceil W(kR,\,2+kl)/k\right\rceil ,
 \tag{1}
\]
at least one of the factorial moments of \(F^k\) in the window is
nonzero. The bound is independent of all phase coefficients. Also, for
any \(k\) and \(0\le h<k\), if the sequence of moments of \(F\) at
orders \(kn+h\) is eventually zero, it is identically zero on that
progression. The same finite-support exclusion applies to finite
complex linear combinations of exponential generating series of
such phases. The window (1) can exceed the number of monomials.

The same theorem holds for a product Gamma functional with fixed positive
rational shapes \(\alpha_i\), namely
\(\Lambda(x_1^{e_1}\cdots x_m^{e_m})=\prod_i(\alpha_i)_{e_i}\).
Ordinary factorial is the case \(\alpha_i=1\).

## A finite system for one summand

Work over \(\mathbb C\), using formal moment series at the origin.
For the one-variable functional with shape \(\alpha>0\), put
\[
 I_\phi(s)=\sum_{n\ge0}\Lambda(\phi q^n)\frac{s^n}{n!},
 \qquad u=1/s.
\]
For any polynomial \(a\), the identity
\(\Lambda(x a)=\Lambda(\alpha a+x a')\) follows on monomials
from \((\alpha)_{n+1}=(\alpha+n)(\alpha)_n\).
Applied coefficientwise to \(a e^{s q}\) it gives
\[
             s I_{x a q'}=I_{(x-\alpha)a-xa'}.                  \tag{2}
\]
For \(q(x)=c_0+a x\) with \(a\ne0\), the series is
\(I_1=e^{c_0s}(1-a s)^{-\alpha}\). Its scalar differential system
has leading action \(c_0\) and, after removal of that exponential,
positive exponent \(\alpha\) in \(u\).

For \(d=\deg q\ge2\), use the \(d\) coordinates \(I_\phi\) for a basis
of polynomials of degree at most \(d-1\). Divide
\[
     q\phi=xq'a_\phi+r_\phi,\qquad
     (x-\alpha)a_\phi-xa'_\phi=b_\phi+\gamma_\phi xq',
\]
where \(\deg r_\phi,\deg b_\phi\le d-1\) and
\(\gamma_\phi\) is constant since \(\deg a_\phi\le d-1\).
Then (2), used twice, gives the exact system
\[
            I_\phi'=I_{r_\phi}+u I_{b_\phi}
                                  +u^2\gamma_\phi I_{x-\alpha}. \tag{3}
\]
The last amplitude has degree one, so it belongs to the basis span.
The matrix has the form \(A_0+u A_1+u^2 A_2\).

Identify the amplitudes with \(\mathbb C[x]/(xq')\).
The leading action is multiplication by \(q\).
At every nonzero critical point \(\xi\) with \(q-q(\xi)\)
of order \(h\), \(xq'\) has order \(h-1\). At the endpoint \(0\),
if \(q-q(0)\) has order \(h\), \(xq'\) has order \(h\).
On each local quotient the leading action is the scalar \(q(\xi)\).
Thus \(A_0\) is diagonalizable, including at critical points of
arbitrary multiplicity and when critical values coincide.

For an amplitude in the eigenspace of value \(v\),
\(r_\phi=v\phi\) and \(a_\phi=(q-v)\phi/(xq')\) exactly.
At a nonzero root \(\xi\) in that block, the quotient
\((q-v)/(xq')\) is \((x-\xi)/(h\xi)+O((x-\xi)^2)\).
If the leading jet of \(\phi\) has degree \(j\), the leading
multiplier of \((x-\alpha)a_\phi-xa_\phi'\) on that jet
is \(-(j+1)/h\).
At the endpoint the quotient is \(1/h+O(x)\),
and the corresponding multiplier is \(-(\alpha+j)/h\).
Reduction modulo \(xq'\) preserves these jets. Projecting onto the
same leading action and using its separate local jet filtrations
therefore gives triangular residue blocks with strictly
negative rational diagonal entries in the \(s\)-equation.
Distinct roots with one value cause no mixing of these jets:
the quotient expression is a local multiple of \(\phi\),
and the indicated derivative preserves the orders just computed.

## Formal solutions at infinity

Here and below derivatives denoted by primes are in \(s\).
After a constant change of basis diagonalizing \(A_0\), there
is a formal gauge \(U(u)=1+O(u)\) which separates distinct
eigenvalues of \(A_0\) into blocks. At each power \(u^j\) the
off-block coefficients of the gauge are solved using the
nonzero difference of the two eigenvalues. Differentiating
\(u^j\) in \(s\) raises its order, so the derivative term
does not enter this step. In particular the \(u\)-coefficient
inside a block is unchanged. On removing the exponential
of the leading action \(v\), the block is then regular
singular in \(u\), with residue the negative of the
diagonal block of \(A_1\).

For clarity, such a block has a formal basis of the form
\(u^\lambda\sum_{n\ge0}u^n v_n(\log u)\) with each \(v_n\) a
polynomial vector and \(\lambda\) ranging over residue
eigenvalues. If the residue matrix is \(B\), the leading
vectors are obtained by exponentiating the nilpotent part
on its generalized eigenspace at \(\lambda\).
Subsequent coefficients satisfy equations with left-hand
operator \((\lambda+n)I-B+\partial_{\log u}\).
On the nonresonant eigenspaces it is invertible on
polynomial vectors by a finite derivative expansion.
On the resonant eigenspace it is a nilpotent constant
matrix plus ordinary differentiation, solved by multiplying
by its polynomial exponential integrating factor.
The leading generalized eigenvector bases give an invertible
matrix of formal columns. Thus the solutions have the
required Laurent–Puiseux and logarithmic form with no
analytic convergence assumption.

For (3) all the residue eigenvalues in the \(u\)-equation
are strictly positive. It follows that **every** formal
fundamental output column for \(I_1\) has positive \(u\)-order
after removing its exponential action. This proof includes
arbitrary critical multiplicities and root collisions with the
endpoint. The conclusion holds for linear \(q\) by the explicit
one-coordinate system.

## Independent sums, sampled orders, and scalar elimination

For independent variables the exponential generating series multiply:
\[
                G(s)=\sum_{n\ge0}\Lambda(F^n)\frac{s^n}{n!}
                     =e^{cs}\prod_i I^{(i)}_1(s).
\]
This is a rational-system output using tensor products of the
systems (3), the one-coordinate systems for linear terms, and
the constant action \(c\). Tensor products of formal fundamental
columns add their exponential actions and multiply their scalar
observable amplitudes. The product still has positive \(u\)-order,
even when some exponential actions cancel. Finite direct sums of
such systems preserve this property for the sum of their outputs.
If \(s\) is replaced by \(\zeta s\) for \(\zeta\ne0\), the action
is multiplied by \(\zeta\) and the same output order holds.

We give the scalar consequence algebraically. Write a system as
\(Y'=A(s)Y\) with rational output row \(b(s)\).
Suppose it has an invertible formal matrix of columns \(v_j\)
with exponential factors removed, so \(v_j'=A v_j-a_jv_j\)
for constant actions \(a_j\). They lie in
\(\mathbb C(\log u)((u^{1/L}))\) for some common \(L\); differentiation
is specified by \(d u/ds=-u^2\) and \(d\log u/ds=-u\).
Suppose that each \(b v_j\) has positive \(u\)-order or is zero.
Define \(b_0=b,\ b_{i+1}=b_i'+b_i A\), and let \(r\) be
the first dependent index of these rows over \(\mathbb C(s)\).
For a nonzero \(b\) this gives a scalar equation of order \(r\).

The first \(r\) rows applied to the formal fundamental matrix
have rank \(r\), so select \(r\) columns to get an invertible
jet matrix \(H\). If \(C(s)\) is the scalar companion
matrix, differentiation yields \(H'=CH-H\operatorname{diag}(a_j)\)
on the selected columns. If a polynomial \(P(s)\) solves the
scalar equation, its jet vector \(J\) satisfies \(J'=CJ\).
It follows that \(\gamma=H^{-1}J\) has
\(\gamma'_j=a_j\gamma_j\). Differentiation raises \(u\)-order
by at least one, so \(\gamma_j=0\) if \(a_j\ne0\).
When \(a_j=0\), \(\gamma_j\) is constant. Indeed the equation
on a coefficient \(R(\ell)\) of \(u^v\), with \(\ell=\log u\),
is \(R'(\ell)+vR(\ell)=0\). A nonzero rational function cannot
have a nonzero constant logarithmic derivative, leaving only
\(v=0\) and constant \(R\). Thus the first component of \(J\)
is a constant linear combination of the positive-order outputs
\(b v_j\) of action zero. This precludes a nonzero polynomial
in \(s=1/u\). If the output is identically zero, the assertion
that it cannot be a nonzero polynomial needs no equation.

For the progression \(kn+h\) use the output
\[
       G_{k,h}(s)=\frac1k\sum_{\zeta^k=1}\zeta^{-h}G(\zeta s).
\]
It is the direct sum output that retains precisely those
degrees. It follows that an eventually zero progression
has all moments zero on the progression. This proves also
the claim about finite linear combinations. For \(h=0\)
the output has constant term one and is nonzero.

For an explicit bound on the scalar width, take a nonzero
common polynomial denominator \(Q\) of degree \(q\) with
\(Q A\) and \(Q b\) of degree at most \(q\). In the tensor
system above it suffices to use \(s^2\) times the product of
the \(l\) linear denominators \(1-a_i s\). For \(G_{k,0}\)
one can use \(s^2\prod_{i,\zeta}(1-a_i\zeta s)\), of degree
\(2+kl\), and the dimension is \(kR\).
The numerator rows \(Q^{j+1}b_j\) have degree at most
\((j+1)q\), by differentiating and using \(QA\).
Signed \(r\)-row minors in \(r\) columns on which the
first \(r\) rows are independent give the scalar equation
with polynomial coefficients of degree at most
\(K=q((r+1)(r+2)/2-1)\).
In fact the minor for row \(j\) omits its degree bound
\((j+1)q\), and translating the relation on numerator rows
to the relation on the \(b_j\)'s restores at most \(jq\)
through the multiplier \(Q^j\).
The relation on these columns holds on the whole row space:
restriction is injective on the span of the first \(r\)
rows and the last row lies in that span.

With \(\theta=s\partial_s\), multiply the scalar equation
by \(s^r\) to write it as
\(\mathcal L=\sum_{i=0}^{K+r}s^i P_i(\theta)\).
Its width is at most \(K+r\le W(kR,2+kl)\) for \(G_{k,0}\).
If consecutive coefficients of \(G_{k,0}\) vanish for at least
that width at a positive start, the Taylor truncation \(B\)
before the block satisfies \(\mathcal L B=0\): the terms from
\(G_{k,0}-B\) start above the possible degree of \(\mathcal L B\).
But \(B(0)=1\), contrary to absence of nonzero polynomial
solutions. A block of \(N\) consecutive zero moments of \(F^k\)
gives \(kN\) consecutive zero coefficients of \(G_{k,0}\),
starting at the corresponding multiple of \(k\).
Taking \(N\) as in (1) proves the window.
All field extensions used in the proof preserve the formal
factorial moments; the field generated by the finitely many
phase coefficients embeds into \(\mathbb C\).

## Adding quadratic blocks in disjoint variables

The same argument permits arbitrary nonconstant inhomogeneous
quadratics \(P_j\) on disjoint blocks of \(m_j\) ordinary factorial
variables, in addition to the separate \(q_i\).
The [quadratic face system](FACTORIAL_QUADRATIC_EVENTUAL_NONVANISHING.md)
has dimension at most \(2^{m_j}\), with a common denominator
of degree at most \(2^{m_j}+1\) for its matrix and output.
It also has the needed positive order for each output after
its exponential action is removed. The only action block that
can have residue exponent zero is the block of the constant
face. Its constant basis vector contributes zero at \(u=0\)
to a nonconstant output, by the face elimination and the
vanishing column at equal action in that proof.

Thus for \(F=c+\sum_i q_i(X_i)+\sum_j P_j(X_{I_j})\)
with all blocks disjoint, the conclusions hold with
\(R'=(\prod_i d_i)\prod_j 2^{m_j}\),
\(q'=2+l+\sum_j(2^{m_j}+1)\), and window
\(\lceil W(k R',k q')/k\rceil\).
Indeed products of the individual denominators clear both the
tensor matrix and its output. After scaling by \(k\)-th roots of
unity one can use the product of the \(k\) resulting denominators,
of degree at most \(kq'\). All output orders remain positive
since the symbols are on disjoint variables and the exponential
moment series multiply.
