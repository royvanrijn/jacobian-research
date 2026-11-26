# Arbitrarily deep fixed-prime obstruction with one shared row

Fix a prime \(p\), a positive precision \(L\), and a positive integer
order bound \(J\ge2\). For every nonempty open interval
\(I\subset(0,\infty)\) and every positive relative-gap cap, there is
an ordinary fully shared profile with exactly one unequal row
\((a,a+d)\) of curvature \(d^2/(a+d)\in I\) and relative gap below
the cap, with the following property. One common positive rational
rescaling of the variable, and division of each moment by its constant
coefficient, gives polynomials congruent to \((1+z)^m\) modulo \(p^L\)
for every \(1\le m\le J\). The congruences are in \(\mathbb Z_{(p)}[z]\).

There is also a balanced version with precisely the two mirrored rows
\((a,a+d),(a+d,a)\) and total curvature in \(I\).
In that version **no rescaling of the argument is needed**.
The relative gap can again be arbitrarily small.
In either version, for fixed \(p,J,I\), the maximum exponent may
be taken \(O_{p,J,I}(p^{6L})\) if a separate relative-gap cap is
not prescribed.

The examples concern the original factorial moment polynomials.
The one-row profiles are coprime at every consecutive pair by the
[ordinary one-variable theorem](ORDINARY_BINOMIAL_RESIDUAL_REGION.md).
In the mirrored version one may also arrange all-degree angular
separation, by taking curvature approaching any positive rational
strength in \(I\) and using the
[diffuse theorem for the empty core](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md).
Thus arbitrarily deep congruence to a common-root pair at one prime
can occur even in minimum dimension and in a prescribed positive
curvature band. This is a fixed-prime obstruction. It makes no
assertion about simultaneous congruences at different primes with
one or two rows.

## Digit complementation and modular units

Choose \(h=p^M\) with \(M\ge L-1\), and \(T\) with
\(p^{T-L}>J(1+h)\). Set
\[
                D=p^T-1,\qquad d=hD,\qquad a=D+p^{2T}A ,
\]
where \(A\ge1\) is an integer to be chosen. Then for \(0\le k\le m\le J\)
write
\[
       ma+dk=X+D s,\qquad X=p^{2T}mA,\qquad s=m+hk .
\]
The formula \(Ds=(s-1)p^T+(p^T-1-(s-1))\) shows that its low
\(T\) digits are the complements of the digits of \(s-1\).
Since \(1\le s<p^{T-L}\), the digit sum of \(Ds\) is
\((p-1)T\), and \(Ds<p^{2T}\). Addition of \(X\) creates
no carry into the high digits. Legendre's factorial valuation
formula therefore gives exactly
\[
 v_p((ma+dk)!)-v_p((ma)!)=\frac{dk}{p-1}.                    \tag{1}
\]
Here \(d/(p-1)\) is an integer. Legendre's formula follows by
counting multiples of each \(p^j\).

For a positive integer \(y\) let
\(\operatorname{unit}_p(y)=y/p^{v_p(y)}\), and denote the same
operation on \(y!\) by \(U(y!)\).
The key identity is that there is a unit \(\xi\bmod p^L\),
independent of \(m,k,A\), such that
\[
                    U((ma+dk)!)/U((ma)!)\equiv \xi^k
                    \pmod{p^L}.                            \tag{2}
\]
To prove it, put \(x_0=X/p^T=p^T mA\) and define the block
product
\[
       B(j)=\prod_{r=1}^{p^T-1}\operatorname{unit}_p(jp^T+r),
                 \qquad j\in\mathbb Z_{\ge0}.
\]
Complete the factorial \((X+Ds)!\) to
\((p^T(x_0+s)-1)!\). The additional factors are
\(p^T(x_0+s)-r\) for \(1\le r\le s-1\). Each has unit
congruent to \(-\operatorname{unit}_p(r)\pmod{p^L}\).
Indeed \(v_p(r)<T-L\), while \(p^T(x_0+s)\) is divisible
by \(p^T\). Next separate the multiples of \(p^T\)
from the complete factorial. Their unit product is
\(U((x_0+s-1)!)\). Since \(x_0\) is also divisible by \(p^T\),
the quotient \(U((x_0+s-1)!)/U((s-1)!)\) is congruent to
\(U(x_0!)\) by the same margin \(T-v_p(r)>L\).
Altogether,
\[
 U((X+Ds)!)\equiv
       (-1)^{s-1} U(x_0!)\prod_{j=0}^{x_0+s-1}B(j)
                        \pmod{p^L}.                         \tag{3}
\]

For \(1\le r<p^T\) one has
\(v_p(jp^T+r)=v_p(r)\le T-1\).
Increasing \(j\) by \(p^{L-1}\) changes the unit by
a multiple of \(p^L\). Thus \(B(j)\) is periodic with
period \(p^{L-1}\) modulo \(p^L\).
In (3) each step \(k\mapsto k+1\) adds \(h\) consecutive
blocks, and \(h\) is a multiple of that period.
Their unit product is independent of the starting index.
Together with the sign \((-1)^h\), it gives the same
geometric multiplier \(\xi\) at every step. This proves (2).

## The two residue statements and the size choice

For one increasing row, the scalar-normalized moment at order
\(m\) has \(k\)-th coefficient
\(\binom mk(ma+dk)!/(ma)!\).
Choose any positive integer representative \(\widehat\xi\) of
the unit \(\xi\). Substitution of
\(z/(p^{d/(p-1)}\widehat\xi)\) for its argument puts the
coefficients in \(\mathbb Z_{(p)}\) and gives the specified
congruence by (1) and (2). The scale is independent of \(m\).

For the two-row profile the factorial coefficient ratio
to the constant contribution is
\[
      \frac{(ma+dk)!(ma+d(m-k))!}{(ma)!(ma+dm)!}.
\]
Its valuation is zero by (1). Its unit reduces to
\(\xi^k\xi^{m-k}/\xi^m=1\pmod {p^L}\) by (2).
The binomial coefficient is unchanged, giving the assertion
in the original argument.

For a target one-row curvature \(S_*>0\), take
\(A=\lfloor h^2/S_*\rfloor\) for sufficiently large \(h\).
Then
\[
 \frac{d^2}{a+d}
  =\frac{h^2(1-p^{-T})^2}
         {A+(h+1)(p^{-T}-p^{-2T})}\longrightarrow S_*,
 \qquad \frac d{a+d}\longrightarrow0 .
\]
For the paired construction use half the desired total
curvature as \(S_*\). For the size bound one can take
\(M=L+c\) and \(T=M+L+c_J\), with fixed \(c,c_J\)
depending on \(p,J,I\). Choose \(p^{c_J}>2J\) and then
choose \(c\) large enough to place all resulting curvatures
in \(I\). In this choice \(a+d=O_{p,J,I}(p^{6L})\).

In particular, for any consecutive pair \(n,N\le J\),
the resultant after the indicated normalizations has
valuation at least \(nL\). Its Sylvester matrix is the exact matrix
for the pair \((1+z)^n,(1+z)^N\), plus \(p^L\) times an
integral matrix over \(\mathbb Z_{(p)}\). The binomial matrix
has rank \(N\) over \(\mathbb Q_p\), since its kernel consists
of pairs \((-B(1+z),B)\) with \(\deg B<n\).
In the multilinear determinant expansion a nonzero term
therefore uses at least \(n\) columns from the error.
This explains why fixed-prime coefficient screens can
require arbitrarily deep precision even in these dimensions.
The [simultaneous prime-screen obstruction](MULTIVARIATE_QUANTITATIVE_PRIME_SCREEN_OBSTRUCTION.md)
instead treats many primes at once with a larger row count.


## Exact unit replay

The code evaluates factorial units by first taking the product of
the factors not divisible by \(p\) and then iterating with the
quotient by \(p\). Complete blocks modulo \(p^L\) are counted by
division, so no large factorial is formed. The chosen parameters
exercise both signs, small primes, and several precisions.

~~~python
from math import comb

def check(p, L, J, M, T, A):
    h = p**M
    D = p**T-1
    d = h*D
    a = D+p**(2*T)*A
    modulus = p**L
    assert M >= L-1 and p**(T-L) > J*(1+h)
    prefix = [1]
    for i in range(1, modulus):
        prefix.append(prefix[-1]*(i if i % p else 1) % modulus)

    def unit_factorial(n):
        value = 1
        while n:
            q, r = divmod(n, modulus)
            value = value*pow(prefix[-1], q, modulus)*prefix[r] % modulus
            n //= p
        return value

    def factorial_valuation(n):
        value = 0
        while n:
            n //= p
            value += n
        return value

    xi = unit_factorial(a+d)*pow(unit_factorial(a), -1, modulus) % modulus
    for m in range(1, J+1):
        arguments = [m*a+d*k for k in range(m+1)]
        vals = [factorial_valuation(x) for x in arguments]
        units = [unit_factorial(x) for x in arguments]
        for k in range(m+1):
            assert vals[k]-vals[0] == d*k//(p-1)
            assert units[k]*pow(units[0], -1, modulus) % modulus == pow(xi,k,modulus)
            assert vals[k]+vals[m-k] == vals[0]+vals[m]
            assert units[k]*units[m-k]*pow(units[0]*units[m], -1, modulus) % modulus == 1
        normalized = [comb(m,k)*units[k]*pow(units[0]*pow(xi,k,modulus), -1, modulus) % modulus
                      for k in range(m+1)]
        assert normalized == [comb(m,k) % modulus for k in range(m+1)]
    return xi

for inputs in [(2,1,10,0,6,1), (2,7,10,7,19,239),
               (3,4,10,4,12,851), (5,3,10,4,11,972),
               (11,2,10,3,8,25)]:
    check(*inputs)
print("Exact digit valuations and common modular unit transport passed.")
~~~
