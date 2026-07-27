# An integral Jacobian-one Hasse failure

This note gives a genuinely integral form of the degree-five Hasse failure.
The map is étale over all of `Spec Z`, the target is integral, and the local
source points lie in `Z_p^3`, including at the primes where the inverse
polynomial has bad reduction.

The paper manuscripts are not used as verification sources here.  The
construction is certified by the independent SymPy, Singular,
dependency-free Python, and PARI/GP checks listed in Section 7.

## 1. The map and the main theorem

Put

\[
 t=1+2xy,\qquad q=t^2z-y^2(1+3t)
\]

and define \(F=(F_1,F_2,F_3)\) by

\[
\begin{aligned}
F_1={}&x(1-3xy)+2x^3z-3x^4q^4+3x^5q^5,\\
F_2={}&y-6xq+6t^2x^2q^4-5t^2x^3q^5,\\
F_3={}&tq.
\end{aligned}
\]

Every displayed coefficient is integral, and exact differentiation gives

\[
\boxed{\det JF=1.}
\]

Take

\[
\mathbf y=(3,-1,1)
\]

and let

\[
\begin{aligned}
C(S)&=S^3-2S^2+8,\\
Q(S)&=S^2-S+6,\\
P(S)&=C(S)Q(S)\\
&=S^5-3S^4+8S^3-4S^2-8S+48.
\end{aligned}
\]

### Theorem 1.1

The complete fiber over \(\mathbf y\) is

\[
\boxed{
F^{-1}(\mathbf y)_{\mathbb Q}
\simeq \operatorname{Spec}\mathbb Q[S]/(P).
}
\]

Moreover,

\[
\mathbf y\in F(\mathbb Z_p^3)
\quad\text{for every prime }p,
\qquad
\mathbf y\in F(\mathbb R^3),
\]

but

\[
\boxed{\mathbf y\notin F(\mathbb Q^3).}
\]

Thus \(F:\mathbb A^3_{\mathbb Z}\to\mathbb A^3_{\mathbb Z}\) is an
everywhere-étale polynomial map with a genuinely integral Hasse-failing
fiber.  Its geometric fiber degree is five, which is the minimum permitted
by the degree-four arithmetic obstruction in
[`MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md`](MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md).

## 2. Integral unimodular normalization

Start from the root-engineered quadratic gauge for

\[
G(S)=P(S)-48
=S^5-3S^4+8S^3-4S^2-8S.
\]

Its first coefficients are

\[
(g_1,g_2,g_3,g_4,g_5)=(-8,-4,8,-3,1).
\]

If its coordinates are denoted by \((\Pi,B,C_0)\), it has determinant
\(-2\).  Apply the source dilation

\[
x_{\rm old}=2x
\]

and the target lattice change

\[
(\Pi,B,C_0)\longmapsto
\left(\frac{C_0}{4},\,B-\Pi,\,\Pi\right).
\]

The quarter is integral after the source dilation: the low term becomes

\[
\frac{2x(5-3(1+2xy))}{4}=x(1-3xy),
\]

and every higher term has the required power of two.  The determinant is

\[
(-2)\cdot 2\cdot\left(-\frac14\right)=1.
\]

The target \((3,-1,1)\) corresponds to

\[
(\Pi,B,C_0)=(1,0,12).
\]

The inverse equation is therefore

\[
G(S)-\frac{g_1}{2}(BS^2+C_0)
=G(S)+48=P(S).
\]

The independent Singular certificate constructs both quotient-ring maps
and verifies that they compose to the identity.

## 3. Arithmetic Hasse certificate

The two factors are irreducible.  The cubic has no integer root, and the
quadratic has negative discriminant.  Their discriminants are

\[
\operatorname{Disc}(C)=-1472=-23\cdot 8^2,
\qquad
\operatorname{Disc}(Q)=-23.
\]

Thus the cubic has Galois group \(S_3\), and \(Q\) cuts out its discriminant
field.  At every prime unramified in the common splitting field, each
Frobenius element either fixes a cubic root or acts trivially on the
quadratic field.  Hence one of the two factors has a local root.

The factor resultant is

\[
\operatorname{Res}(C,Q)=392=2^3\,7^2.
\]

Consequently the product-discriminant primes requiring separate integral
checks are \(2,7,23\):

\[
\begin{array}{c|c|c}
p&\text{root source}&\text{check}\\ \hline
2&Q,\ S\equiv0&Q'(0)\equiv1\pmod2\\
7&P,\ S\equiv1&P'(1)\equiv1\pmod7\\
23&C,\ S\equiv7&C'(7)\equiv4\pmod{23}.
\end{array}
\]

At every odd prime, a simple root \(s\) gives an integral source through

\[
\begin{aligned}
d&=-\frac{P'(s)}8,&t&=d^{-1},&
x&=\frac{s}{2d},\\
\beta(s)&=1-4s+\frac32s^2-\frac58s^3,&
y&=-\beta(s)-s,\\
z&=d^2\bigl(d+y^2(1+3t)\bigr).
\end{aligned}
\]

At \(2\), lift the root of \(Q\) with \(s=2u\).  Then

\[
d=-(4u-1)(u^3-u^2+1)\in\mathbb Z_2^\times,
\]

and

\[
x=\frac{u}{d},\qquad
y=-1+6u-6u^2+5u^3,\qquad
z=d^2\bigl(d+y^2(1+3d^{-1})\bigr)
\]

are all \(2\)-adic integers.  This proves
\(\mathbf y\in F(\mathbb Z_p^3)\) at every prime.

The odd-degree polynomial \(P\) has a real root, and the reconstruction is
regular because \(P\) is squarefree.  Neither irreducible factor has a
rational root.  The quotient-ring fiber description therefore proves the
real and rational assertions in Theorem 1.1.

## 4. The integral nonproper boundary

Let

\[
K_3=\mathbb Q(\alpha),\quad C(\alpha)=0,
\qquad
K_2=\mathbb Q(\beta),\quad Q(\beta)=0.
\]

Both maximal orders have field discriminant \(-23\).  The Zariski-main
finite normalization compactification of the integral fiber is

\[
\overline X
=\operatorname{Spec}\mathcal O_{K_3}
\sqcup
\operatorname{Spec}\mathcal O_{K_2}.
\]

The actual étale fiber \(X=F^{-1}(\mathbf y)\) is the open subscheme obtained
by deleting exactly the prime ideals in the rows marked `deleted` below.
The valuation columns record

\[
(v(d),v(t),v(x),v(y),v(z)).
\]

| component | prime branch | \(e,f\) | valuations | status |
|---|---:|---:|---:|---|
| \(K_3\) | \(2\), inert | \(1,3\) | \((0,0,0,0,0)\) | retained |
| \(K_3\) | \(7\), \(\alpha\equiv1\) | \(1,1\) | \((0,0,0,0,0)\) | retained |
| \(K_3\) | \(7\), degree two | \(1,2\) | \((1,-1,-1,0,1)\) | deleted |
| \(K_3\) | \(23\), \(\alpha\equiv7\) | \(1,1\) | \((0,0,0,0,0)\) | retained |
| \(K_3\) | \(23\), ramified \(\alpha\equiv9\) | \(2,1\) | \((1,-1,-1,0,1)\) | deleted |
| \(K_2\) | \(2\), \(\beta\equiv0\) | \(1,1\) | \((0,0,0,3,0)\) | retained |
| \(K_2\) | \(2\), \(\beta\equiv1\) | \(1,1\) | \((-3,3,2,-3,-12)\) | deleted |
| \(K_2\) | \(7\), inert | \(1,2\) | \((1,-1,-1,0,1)\) | deleted |
| \(K_2\) | \(23\), ramified | \(2,1\) | \((1,-1,-1,0,1)\) | deleted |

Thus exactly five prime ideals are missing.  The three retained degree-one
branches supplying the exceptional integral local points are

\[
\beta\equiv0\pmod2,\qquad
\alpha\equiv1\pmod7,\qquad
\alpha\equiv7\pmod{23}.
\]

Finally,

\[
N_{K_3/\mathbb Q}(d)=-1127,\qquad
N_{K_2/\mathbb Q}(d)=\frac{1127}{8},
\qquad
1127=7^2\cdot23.
\]

Together with the displayed powers of two in the reconstruction, this proves
that no prime outside \(2,7,23\) supports a deleted boundary point.

## 5. More failures in the same integral map

For an arbitrary integral target

\[
(u,v,1)\in\mathbb Z^3,
\]

the inverse polynomial of the same fixed map is

\[
\boxed{
E_{u,v}(S)
=S^5-3S^4+8S^3+4vS^2-8S+16u.
}
\]

Besides \((u,v)=(3,-1)\), two further targets have exact Hasse certificates:

\[
\begin{array}{c|l|c|l}
(u,v,1)&E_{u,v}(S)&
\text{discriminant squareclass}&
\text{bad primes}\\ \hline
(-66,-1,1)&
(S^3-8S^2+36S-88)(S^2+5S+12)&
-23&
2,5,23,8677\\
(-12865,-637,1)&
(S^3-18S^2+112S-1240)(S^2+15S+166)&
-439&
2,5,31,439,555253.
\end{array}
\]

In each row both factors are irreducible, the cubic is \(S_3\), and the
quadratic is its discriminant field.  At every listed odd bad prime the
inverse polynomial has an explicitly recorded simple root; at \(2\) the
artifact records an actual source point modulo \(2\).  Since \(\det JF=1\),
these residue source points lift to \(\mathbb Z_p^3\).  The same Frobenius
covering handles all remaining primes.  Hence both rows are additional
integral Hasse failures for the one fixed map \(F\).

## 6. Scope of the bounded search

Write a monic \(3+2\) factorization as

\[
(S^3+aS^2+bS+c)(S^2+\delta S+\epsilon).
\]

Coefficient comparison with \(E_{u,v}\) gives

\[
\begin{aligned}
\delta&=-3-a,\\
\epsilon&=a^2+3a+8-b,\\
c\delta+b\epsilon&=-8,\\
16u&=c\epsilon,\\
4v&=c+b\delta+a\epsilon.
\end{aligned}
\]

The exact search exhausts

\[
|a|\le2000,\qquad |b|\le2000
\]

subject to these identities, irreducibility, a common cubic/quadratic
discriminant field, and all bad-prime local conditions.  It returns exactly
the original target and the two new targets in Section 5.

This is an exhaustive computation only in the stated coefficient box.  It
is not a proof that the fixed map has only three such targets, and it is not
an infinitude theorem.  The complete machine-readable result, including
simple roots and source points at every bad prime, is
[`integral_hasse_pencil_search.json`](../artifacts/generated-results/integral_hasse_pencil_search.json).

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_integral_hasse_jacobian_one.py
/opt/homebrew/bin/Singular -q scripts/verify_integral_hasse_jacobian_one_independent.sing
python3 scripts/verify_integral_hasse_arithmetic_independent.py
/opt/homebrew/bin/gp -q scripts/verify_integral_hasse_boundary.gp
.venv/bin/python scripts/search_integral_hasse_pencil.py
```

The first checker verifies the integral map, determinant, quotient-ring
reconstruction, and local integrality.  Singular independently computes the
Jacobian and verifies mutually inverse quotient-ring maps.  The
dependency-free checker proves the finite \(S_3\) covering and exceptional
Hensel data without importing SymPy.  PARI/GP computes the complete
normalization-boundary valuation table.  The final command reproduces the
bounded fixed-pencil search and compares it exactly with the archived JSON
result.
