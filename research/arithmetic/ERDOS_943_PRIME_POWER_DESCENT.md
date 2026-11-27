# Ideal-class descent for square-full sums at prime powers

Let \(N\) be a positive integer. A coprime representation of \(N\) is
\[
                  a^2b^3+c^2d^3=N,                         \tag{1}
\]
where \(a,b,c,d\) are positive integers, \(b,d\) are squarefree and the
two summands are coprime. A square-full summand has a unique such
square-times-squarefree-cube decomposition.

A split prime ideal can have order much smaller than the exponent in
a norm equation. In that case (1) descends to a smaller prime-power
target. The descent preserves the squarefree kernels after removing
factors supported on the exponent of the power substitution. This
makes the nonfundamental part negligible in the balanced range of
Erdős #943. The count of fundamental representations remains open.

The class-order mechanism is related to Heuberger–Le, Section 6,
where the least exponent for a quadratic-form class and composition
control representations of powers. Their setup there uses an odd
base coprime to the discriminant. The proof below works directly in
the maximal order, includes the dyadic correction for \(N=2^m\),
and keeps the divisibilities needed to recover square-full summands.
[Heuberger–Le, Theorems 6.1–6.2 and Lemma 6.3](https://www.math.aau.at/heuberger/publications/pdf/nagell.pdf)

## The oriented ideals

Suppose \(N>4\) is odd or divisible by \(8\), and \(\max(b,d)\ge5\).
Write \(K=\mathbb Q(\sqrt{-bd})\), \(\pi=\sqrt{-bd}\).
The coprimality in (1) makes \(bd\) squarefree and all four variables
units modulo \(N\). The field has unit group \(\{\pm1\}\).
Let \(t=1\) when \(N\) is odd and \(t=2\) when \(8\mid N\), and put
\[
          n=N/t^2,\qquad \alpha=(b^2 a+dc\pi)/t,
          \qquad \mathrm{N}_{K/\mathbb Q}(\alpha)=bn.        \tag{2}
\]
If \(t=2\), both summands are odd and \(b+d\equiv0\bmod8\)
from (1), since \(a,c\) are odd. Thus \(bd\equiv7\bmod8\),
\(2\) splits in \(K\), and \(\alpha\) is integral in its maximal order.
Its two coefficients in the basis \(1,\pi\) are half odd integers.
Every element of \(2\mathcal O_K\) has integral coefficients in that
basis, so \(\alpha\notin2\mathcal O_K\).

Let \(\mathfrak b\) be the product of the ramified prime ideals at
rational primes dividing \(b\), with \(\mathfrak b=\mathcal O_K\) if
\(b=1\). Its norm is \(b\) and \(\mathfrak b^2=(b)\).
At an odd prime \(p\mid n\), \(\alpha\notin p\mathcal O_K\) since
\(b,d,t\) are units at \(p\) and \(a,c\) cannot both be divisible
by \(p\). No prime dividing \(n\) ramifies because \((bd,n)=1\).
Thus each prime dividing \(n\) splits, and only one of its two prime
ideals divides \(\alpha\): in the inert case divisibility of the
norm would give \(\alpha\in p\mathcal O_K\), and in the split case
both primes would also give that containment. The dyadic argument
above already handles \(p=2\). Taking ideal norms now gives
\[
       (\alpha)=\mathfrak b\prod_{p\mid n}
                             \mathfrak p_p^{v_p(n)},         \tag{3}
\]
where \(\mathfrak p_p\) is the prime above \(p\) selected by \(\alpha\).

## Descent with the kernel divisibilities

**Theorem SquareFullIdealClassRepair.** In (1)–(3), suppose
\(n=M^k\) with \(M>1,k\ge2\), and \(k\) is odd or \(b=1\). Let
\[
 \mathfrak J=\prod_{p\mid M}\mathfrak p_p^{v_p(M)},\qquad
 f=\gcd(b,\operatorname{rad}k),\qquad g=\gcd(d,\operatorname{rad}k).
\]
If \(\mathfrak b\mathfrak J\) is principal, there are nonzero
integers \(u,v\) with
\[
            f(b/f)^3u^2+g(d/g)^3v^2=4M.                     \tag{4}
\]
The two summands in (4) have gcd dividing \(4\). For fixed \(k,t,f,g\)
the ordered pair in (4) recovers the original representation in
(1) up to an absolute bounded multiplicity. For fixed \(N,k\), the
number of representations to which the descent applies is
\(O_\varepsilon(N^{3/(10k)+\varepsilon})\).

Choose a generator \(\zeta\) of \(\mathfrak b\mathfrak J\). It is
integral, of norm \(bM\), and (3) gives
\[
              \alpha=\pm \zeta^k/b^{(k-1)/2}                \tag{5}
\]
if \(k\) is odd. For \(b=1\) the formula is \(\alpha=\pm\zeta^k\)
for any \(k\). Indeed \(\mathfrak b^k=\mathfrak b (b)^{(k-1)/2}\)
for odd \(k\), and the two elements with the same ideal differ by a unit.

Write \(2\zeta=X_0+y\pi\) with integral coefficients. The norm
makes \(b\mid X_0^2\), hence \(b\mid X_0\). Put \(X_0=bx\), so
\[
                              bx^2+dy^2=4M.                 \tag{6}
\]
Neither \(x\) nor \(y\) is zero, because then \(\zeta\) would be real
or purely imaginary and its ideal would be invariant under conjugation.
Its ideal includes just one prime above each prime of \(M>1\).

Let \(q\) be an odd prime dividing \(b\). Its exponent in \(b\) is one,
and \(q\nmid 2Md\) implies \(q\nmid y\) by (6). In this case \(k=2e+1\)
is odd. In the real coefficient of (5), the term with \(2e\) factors
of \(\pi\), after multiplication by \(2^k\) and up to sign, is
\[
                        k b x(-d)^e y^{2e}.
\]
Every other real term is a multiple of \(b^2\): the term with \(2j\)
factors of \(\pi\) contains \(b^{e+1-j}\) after division by \(b^e\).
Since the real coefficient of \(\alpha\) is \(b^2 a/t\), comparing
modulo \(q^2\) gives \(q\mid kx\). Thus \(q\mid x\) if \(q\nmid k\).

For an odd prime \(q\mid d\), (6) gives \(q\nmid x\). In the coefficient
of \(\pi\) in (5), all terms except the one with a single factor of
\(\pi\) vanish modulo \(q\) after multiplying by \(2^k\), because
\(\pi^2=-bd\). The remaining term is a unit times \(ky\) modulo \(q\),
as \(b,x\) are units there. The original coefficient \(dc/t\) vanishes
modulo \(q\), so \(q\mid ky\). This argument also applies for even \(k\)
when \(b=1\). Finally, when \(bd\) is even, \(1,\pi\) is an integral
basis, so \(X_0,y\) are even. If \(b\) is even, the norm of the integer
basis expression shows \(b\mid X_0/2\), and \(x\) is even as well.
Consequently \(b/f\mid x\) and \(d/g\mid y\), including the prime two.
Substitution into (6) gives (4).

An odd prime dividing both summands in (4) would divide \(M\) and
would divide \(x,y\), since \(b,d\) are coprime to \(M\). That puts
\(\zeta\) in both prime ideals above that prime, contrary to its ideal.
If \(2\mid M\), \(b,d\) are odd, so divisibility of both summands by
\(8\) would give \(4\mid x,y\). This would put \(\zeta\) in
\(2\mathcal O_K\), again impossible. If \(M\) is odd the gcd divides
\(4\) by (4). These cases prove the stated gcd bound.

To reconstruct, divide the first summand in (4) by \(f\) and the second
by \(g\). Unique squarefree cube decomposition gives \(b/f,|u|\) and
\(d/g,|v|\). These give \(b,d,|x|,|y|\), and the bounded number of
sign choices gives \(\zeta\) and then \(\alpha\) by (5). The positive
\(a,c\) follow from (2).

For the count, divide both summands of (4) by their gcd \(w\mid4\).
In each resulting summand, the product of the primes occurring to
exponent one divides \(2f\) or \(2g\), respectively. Removing that
product leaves a square-full integer. The resulting primitive
weighted equation at target \(4M/w\) is subject to the
[coefficient-uniform three-tenths bound](ERDOS_943_DUAL_KERNEL_GRAPHS.md#exact-optimization).
The choices of \(f,g,w\) and these two coefficients are bounded by a
fixed power of \(k\), and \(k\ll\log N\).
Summing the bounds \(O_\varepsilon(M^{3/10+\varepsilon})\)
with a smaller intermediate epsilon proves the claimed count.

## Prime powers and the fundamental case

**Theorem PrimePowerNonfundamentalRepresentations.** Let \(N=p^m\) for a
prime \(p\) and \(m\ge1\). In any square-full representation of \(N\),
first divide both summands by their gcd \(p^j\). The two quotients are
still square-full and coprime, with target \(N_0=p^{m-j}\).
Write their decomposition as (1), now with target \(N_0\).
If both normalized kernels are at most \(4\) there are
\(O_\varepsilon(N^\varepsilon)\) original pairs in this sector.

Otherwise orient the normalized pair so that \(d\ge5\).
Use (2)–(3) with \(n=p^e=N_0/t^2\) and write \(r\) for the order
of \([\mathfrak p_p]\) in the ideal class group of \(K\).
If \(b>1\), \(r\) is even and \(e=(r/2)k\) with \(k\) odd.
If \(b=1\), then \(e=rk\).
Call the representation fundamental if \(k=1\).
Uniformly in \(p,m\), the number of original nonfundamental pairs
with both normalized kernels greater than \(1\) is
\(O_\varepsilon(N^{1/10+\varepsilon})\). With one normalized
kernel equal to \(1\) it is \(O_\varepsilon(N^{2/15+\varepsilon})\).
The same bounds hold for ordered pairs.

The gcd of the original summands divides \(p^m\), and the two
quotients are \(p\)-units. Their valuations at every other prime
remain unchanged, so they are square-full. There are at most
\(m\) choices of the common exponent, each of which lifts a
normalized pair uniquely. The bounded-kernel sector follows by
the [fixed-kernel ideal norm bound](ERDOS_943_DUAL_KERNEL_GRAPHS.md#fixed-pairs-and-mixed-coordinates).
The cases \(N_0=2,4\) can also be included in this sector.

When \(b>1,d\ge5\), the ideal \(\mathfrak b\) is nonprincipal.
Indeed a generator of norm \(b\) would give
\(X_0^2+bd y^2=4b\), forcing \(y=0\) and then \(b\) a square.
Its class therefore has order two. By (3),
\([\mathfrak p_p]^e=[\mathfrak b]\), and the unique order-two element
in the cyclic group generated by \([\mathfrak p_p]\) is its
\(r/2\)-th power. This proves the assertion about \(e,k\) and shows
that \(\mathfrak b\mathfrak p_p^{e/k}\) is principal. If \(b=1\),
\([\mathfrak p_p]^e=1\) instead, so \(r\mid e\) and
\(\mathfrak p_p^{e/k}\) is principal. Thus the preceding theorem
applies whenever \(k>1\), with \(M=p^{e/k}\). The number of possible
\(k\) is at most \(m\). In the first case \(k\ge3\), giving the bound
\(N^{1/10+\varepsilon}\) after summing.

In the case \(b=1\), the first summand in (4) is a square.
After division by \(w\mid4\), it is a square or twice a square.
The primitive repaired equation therefore has the shape
\[
                     \kappa_1 x_1^2+\kappa_2 c_1^2d_1^3=m_1,
                      \qquad m_1=4M/w .
\]
Consider a dyadic box \(c_1\asymp C,d_1\asymp D\) with
\(C^2 D^3\ll m_1\). Fixing \(d_1\) gives the fixed-kernel bound,
hence \(O_\varepsilon(m_1^\varepsilon D)\) points.
Fixing \(c_1\) leaves an absolutely irreducible plane cubic
\(A x_1^2+B d_1^3=m_1\) with \(A,B,m_1>0\).
The cubic in \(d_1\) has simple complex roots, so is no square
in \(\mathbb C(d_1)\). Both variables are at most \(m_1^{1/2}\).
The coefficient-uniform Bombieri–Pila bound used in the same
[pointwise argument](ERDOS_943_DUAL_KERNEL_GRAPHS.md#fixed-pairs-and-mixed-coordinates)
therefore gives \(O_\varepsilon(m_1^{1/6+\varepsilon} C)\)
points. Taking the minimum gives
\[
             m_1^\varepsilon\min(D,m_1^{2/3}D^{-3/2})
                          \ll m_1^{4/15+\varepsilon}.
\]
Summing boxes and coefficient choices and using \(k\ge2\)
gives \(N^{2/15+\varepsilon}\), as claimed.

For a numerical example, \(157^3=2^3\,686^2+5^3\,29^2\).
Cubing \(8+5\sqrt{-10}\) and dividing by \(2\) gives the
generator in (2) up to sign. The smaller equation (4) is
\(128+500=4\cdot157\).
A dyadic example is \(2^{17}=3^3\,31^2+5^3\,29^2\).
Here \((3+\sqrt{-15})/2\) has norm \(6\).
Its fifteenth power divided by \(3^7\) gives (2) up to sign
and conjugation. The selected prime above \(2\) has order two:
its class is that of the ramified ideal above \(3\) because
the displayed element has norm \(6\), and that ramified ideal
is nonprincipal by the norm argument above.
The smaller repaired equation is \(3+5=4\cdot2\).

## Large kernels and class-group torsion

**Corollary LargeKernelsRequireTorsion.** Under (1)–(3), if at least
one of \(b,d\) exceeds
\[
                      (4M\,\operatorname{rad}(k)^2)^{1/3}, \tag{7}
\]
then \(\mathfrak b\mathfrak J\) represents a nonidentity class of
order dividing \(k\) whenever \(n=M^k\), \(M>1\), and \(k\) is odd
or \(b=1\). For odd prime \(k\) the order is exactly \(k\).

In fact (3) makes \((\mathfrak b\mathfrak J)^k\) principal in
these cases. If the class itself were principal, the two terms
in (4) would imply \(b^3/f^2,d^3/g^2\le4M\), which gives (7)
as an upper bound for both kernels. This proves the claim.

For an odd perfect power \(N=M^k\) of odd exponent \(k\ge3\),
a primitive representation in a balanced box
\(a,b,c,d\asymp N^{1/5}\) therefore requires nonzero \(k\)-torsion
for large \(N\): the threshold in (7) is at most
\(O(N^{1/9}(\log N)^{2/3})\).
For arbitrary prime-power targets and primitive balanced pairs,
both kernels are large, so every sufficiently large balanced
pair is fundamental by applying the same threshold to the
prime-power descent. Neither the torsion condition nor the
nonfundamental count supplies the remaining fundamental or
all-target representation bound in Erdős #943.
