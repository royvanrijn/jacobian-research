# A common root under monotone block comparison for unequal Jacobi pencils

The leading and trailing block comparison used for symmetric
Gegenbauer matrices does not by itself extend to the complex
[quadratic pencil for unequal linked beta moments](LINKED_BETA_VISHNYAKOVA_PENCIL.md).
There are exact degree-three and degree-four pencils with the
ordinary \((2,1)\) beta diagonals, ordered block weights in the
physical directions, and a common nonreal root. Both polynomials
have positive coefficients, strictly log-convex binomial-normalized
coefficients, and simple roots. The degree-four polynomial is
strictly positive on the real axis. The off-diagonal magnitudes are
chosen independently of the factorial law.

## Which block comparisons are preserved

For the opposed unit rows \((a+1,a),(b,b+1)\), the degree-\(m\)
moment has the determinantal pencil with diagonal
\(ma+j+z(mb+j)\), \(j=1,\ldots,m\), and off-diagonal entries
\(\sqrt z\sqrt{k(k+m(a+b)+1)}\), \(k=1,\ldots,m-1\).
Its normalized factorial coefficients, with the constant
normalization arbitrary, satisfy
\[
 \frac{c_{m,k+1}}{c_{m,k}}
       =\frac{m-k}{k+1}\frac{mb+k+1}{ma+m-k}.
 \tag{1}
\]
Fix \(z_0=r e^{2i\theta}\) with \(r>1\) and \(0<\theta<\pi/2\),
and set \(\zeta=\sqrt r e^{i\theta}\) and
\(p=\Re(\zeta+\zeta^{-1})>0\). Divide the pencil by \(m\zeta p\)
at \(z_0\). Its diagonal at \(x=j/m\) is the affine function
\[
 d(x)=x+\frac{a+rb}{r+1}
                +i q(x-\delta),\qquad
 q=\frac{r-1}{r+1}\tan\theta,\qquad
 \delta=\frac{a-rb}{r-1}.
 \tag{2}
\]
The exact positive squared weights of the normalized pencil are
\[
                  w_{m,k}=\frac{k(k+m(a+b)+1)}{m^2p^2}.
 \tag{3}
\]

Let \(N=n+1\) and \(a+b>0\). The consecutive physical weights
satisfy \(w_{N,k}<w_{n,k}<w_{N,k+1}\) for \(1\le k<n\).
The leading comparison ratio \(w_{n,k}/w_{N,k}\) increases
strictly with \(k\), and the trailing ratio \(w_{N,k+1}/w_{n,k}\)
decreases strictly. In the leading ratio, the only varying factor
is \((a+b+(k+1)/n)/(a+b+(k+1)/N)\), which is increasing.
In the trailing ratio, both \((k+1)/k\) and
\((a+b+(k+2)/N)/(a+b+(k+1)/n)\) decrease. The strict ordering
of the weights follows directly from (3), since
\((k+2)/N\ge(k+1)/n\) and \((k+1)/N>k/n\) on this range.

The diagonal samples \(x=j/n\) lie between the corresponding
leading and trailing samples \(j/N\) and \((j+1)/N\). At a null
vector \(v\) of a normalized block, the exact current identity
\[
  \sqrt{w_{m,k}}\,\Im(\overline v_k v_{k+1})
       =-q\sum_{j=1}^{k}(j/m-\delta)|v_j|^2
 \tag{4}
\]
holds for every internal edge, by taking the imaginary part of
the sum of the first \(k\) row equations against \(\overline v_j\).
The complete sum on the right is zero at \(k=m\).
When \(q>0\), the right-hand side in (4) is strictly positive.
Indeed the factors \(j/m-\delta\) change sign at most once, and
the first and last coordinates of \(v\) are nonzero by the end
row equations. The total vanishing then forces the first and
last nonzero masses to have opposite signs in the sum. Every
proper prefix sum is strictly negative. This proof of the
current sign uses positivity of the weights and (2).

## Exact obstruction

Take \(a=2,\ b=1,\ r=66/41,\ \tan\theta=107/5\). Then (2) gives
\[
                  d(x)=x+\frac{148}{107}+5i(x-16/25).
\]
Keep every diagonal entry. Use squared weights \(P,Q\) in size
three and \(A,B,C\) in size four, with the exact rational values
computed below. They obey
\[
 0<A<P<B<Q<C,\qquad 1<P/A<Q/B,\qquad B/P>C/Q>1.
 \tag{5}
\]
Thus both physical block ratio directions and the diagonal
sampling survive, as does the current sign (4). The size-three
determinant and size-four determinant nevertheless both vanish
at \(z_0\). In the native variable this yields the common factor
\(z^2-2\Re(z_0)z+r^2\).

For the construction, let \(v_j=d(j/3)\), \(u_j=d(j/4)\).
Choose \(P,Q\) to solve the two real equations
\(P/v_1+Q/v_3=v_2\), and set \(C=(26/25)Q\). Choose real \(A\)
so that
\((u_2-A/u_1)(u_3-C/u_4)\) is real, then use this product as
\(B\). The two determinant identities are
\[
 v_1v_2v_3-Pv_3-Qv_1=0,\qquad
 (u_1u_2-A)(u_3u_4-C)-B u_1u_4=0 .
\]
The exact replay also verifies all polynomial positivity and
simplicity claims. The counterexample violates the physical
within-grid identity (3) already at size three, since \(Q<2P\)
whereas \(w_{3,2}>2w_{3,1}\). The coefficient identities (1)
also fail. Any proof by this block comparison must use additional
information about the linked factorial weights.

~~~python
import sympy as sp

z = sp.Symbol("z")
a, b = 2, 1
r, tangent = sp.Rational(66, 41), sp.Rational(107, 5)
d0, delta = (a+r*b)/(r+1), (a-r*b)/(r-1)
slope = tangent*(r-1)/(r+1)
assert (d0, delta, slope) == (sp.Rational(148, 107), sp.Rational(16, 25), 5)
p2 = (r+1)**2/(r*(1+tangent*tangent))
z0 = r*(1-tangent*tangent+2*sp.I*tangent)/(1+tangent*tangent)
common = z*z-2*sp.re(z0)*z+r*r

def diagonal(x):
    return x+d0+sp.I*slope*(x-delta)

v1, v2, v3 = [diagonal(sp.Rational(j, 3)) for j in range(1, 4)]
u1, u2, u3, u4 = [diagonal(sp.Rational(j, 4)) for j in range(1, 5)]
Pvar, Qvar = sp.symbols("P Q", real=True)
equation = sp.expand_complex(Pvar/v1+Qvar/v3-v2)
solution = sp.solve([sp.re(equation), sp.im(equation)], [Pvar, Qvar])
P, Q = solution[Pvar], solution[Qvar]
C = sp.Rational(26, 25)*Q
tail = u3-C/u4
A = sp.cancel(sp.im(u2*tail)/sp.im(tail/u1))
B = sp.cancel((u2-A/u1)*tail)
assert 0 < A < P < B < Q < C
assert 1 < P/A < Q/B and B/P > C/Q > 1
assert Q < 2*P
assert sp.expand(v1*v2*v3-P*v3-Q*v1) == 0
assert sp.expand((u1*u2-A)*(u3*u4-C)-B*u1*u4) == 0

def native_polynomial(m, weights):
    def diag(j):
        return m*a+j+z*(m*b+j)
    previous, current = sp.S.One, diag(1)
    for j, squared_weight in enumerate(weights, 2):
        previous, current = current, sp.expand(
            diag(j)*current-z*m*m*p2*squared_weight*previous)
    return current

for m, weights in ((3, (P, Q)), (4, (A, B, C))):
    f = native_polynomial(m, weights)
    coefficients = [sp.Poly(f, z).nth(k) for k in range(m+1)]
    binomial_normalized = [coefficients[k]/sp.binomial(m, k)
                           for k in range(m+1)]
    assert all(c > 0 for c in coefficients)
    assert all(binomial_normalized[k]**2
               < binomial_normalized[k-1]*binomial_normalized[k+1]
               for k in range(1, m))
    assert any(coefficients[k+1]*(k+1)*(m*a+m-k)
               != coefficients[k]*(m-k)*(m*b+k+1) for k in range(m))
    other_factor, residue = sp.div(f, common, z)
    assert residue == 0 and sp.gcd(f, sp.diff(f, z)) == 1
    if m == 4:
        assert sp.discriminant(common, z) < 0
        assert sp.discriminant(other_factor, z) < 0
        assert sp.Poly(other_factor, z).LC() > 0
print("Exact common factor, two strict block comparisons, and positivity passed")
~~~
