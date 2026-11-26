# Aggregate support depth in polynomial models

A test at each original support prime can be replaced by a bound on
the total prime-power mass that survives the tests. In a binary
polynomial family a Wronskian bounds the formal part of this mass.
A least common multiple of resultants measures the remaining part
at primitive integer parameters. The bound allows the residual
linear combination to change from prime to prime. A counterexample
shows why the resultant term cannot simply be dropped when the list
of residuals grows.

## AS1. A power saving from the aggregate test

Fix a positive integer \(d\). At scale \(H\), let \(A,B,C\)
be pairwise coprime nonconstant homogeneous binary forms in
\(\mathbb Z[s,t]\), each of degree \(d\), with coefficient heights
\(H^{o(1)}\). Put \(F=ABC\). Let \(F_1,\ldots,F_\ell\) be its
distinct primitive irreducible factors over \(\mathbb Q\), with
\[
                     S=\sum_{i=1}^{\ell}\deg F_i .
\]
Suppose each factor has degree at most four. Let
\(E_1,\ldots,E_J\) be a nonempty finite list of nonzero forms, each
an integral linear combination of \(A,B,C\) with all three
coefficients nonzero. The coefficient heights of the \(E_j\) are
uniformly \(H^{o(1)}\). The list may depend on \(H\).

Let \(G_j\) be the primitive formal gcd of \(F\) and \(E_j\), and let
\(R_j\) be the positive absolute homogeneous resultant of
\(F/G_j\) and \(E_j/G_j\). If the second quotient is a nonzero
constant \(k\), use \(R_j=|k|^{\deg(F/G_j)}\). Set
\(\mathcal R=\operatorname{lcm}_{1\le j\le J}R_j\). Assume that
\[
             \mathcal R\le H^{\rho+o(1)}
                    \qquad\text{for fixed }\rho\ge0.
\]
For primitive integer pairs \((s,t)\) with
\(H\le\max(|s|,|t|)<2H\) and \(F(s,t)\ne0\), define
\[
 \mathcal D(s,t)=\operatorname{lcm}_{1\le j\le J}
                  \gcd\bigl(F(s,t),E_j(s,t)\bigr).
\]
Here gcds are nonnegative and \(\gcd(n,0)=|n|\).
If \(\eta\ge0\) is fixed and \(\rho+\eta<3\), the count satisfying
\[
                   \mathcal D(s,t)\ge
                     H^{-\eta}\operatorname{rad}(F(s,t))^3
                                                               \tag{AS1a}
\]
is at most
\[
           H^{\,2-(1-(\rho+\eta)/3)/(2\ell+1)+o(1)} .
\]
At fixed \(d,\rho,\eta\) the allowances are uniform under the
stated subpower height and resultant bounds. In particular a
fixed number of residuals automatically allows \(\rho=0\).
The original third-support-layer condition uses \(J=1\),
\(E_1=A+B-C\ne0\) and implies (AS1a) with \(\eta=0\).
Condition (AS1a) also allows some support primes to contribute
less than three layers when other primes compensate.

### The formal and evaluated gcds

Write \(G=\operatorname{lcm}_jG_j\), taken primitive. Suppose first
that \(A,B,C\) have rank three over the constants. Choose a
rational projective chart for the degree argument whose point at
infinity is disjoint from the roots of \(F\). In this chart use
the nonzero Wronskian with derivative orders \(0,1,2\). At
\(P=F_i\) belonging to \(A\), choose an \(E_j\) maximizing
\(v_P(E_j)\) and use the basis \(A,B,E_j\). Its Wronskian is a
nonzero constant multiple of the original Wronskian. In each
determinant term the derivative orders on \(A,E_j\) sum to at
most three. Thus, writing \(e_i=v_{F_i}(F)\) and
\(m_i=\max_j v_{F_i}(G_j)\), the Wronskian has valuation at
least \(e_i+m_i-3\). For a factor belonging to either other
letter, use the basis given by that letter, one other letter
and the chosen residual. The hypotheses on the coefficients
ensure each chosen list is a basis. Summing and using the degree
bound \(3d-3\) gives
\[
          3d+\deg G-3S\le 3d-3,\qquad
                         \deg G\le 3S-3.                    \tag{AS1b}
\]

If \(A,B,C\) have rank two, their unique constant relation up
to scaling uses all three letters, since nonconstant coprime
forms cannot be proportional. Its coefficients are nonzero
even when some \(E_j\) is proportional to a single original
letter. In the chart used for (AS1b), the two-letter Wronskian
\(A'B-AB'\) is nonzero. The relation allows the pair \(A,B\)
to be replaced by any pair of letters, changing the Wronskian
only by a nonzero scalar. Thus it is divisible by
\(\prod_i F_i^{e_i-1}\) in that chart. Its degree is at most
\(2d-1\), so \(3d-S\le2d-1\) and \(d\le S-1\).
It follows that \(\deg G\le\deg F=3d\le3S-3\).
Rank one is excluded by the same nonproportionality.

The homogeneous resultant of a coprime pair of integral
binary forms gives Bézout identities with the resultant
times a power of \(s\) and a power of \(t\) on their respective
right sides. These follow by applying the adjugate to the
Sylvester map in monomial coordinates with each pure power
as the desired output. Evaluate at a primitive pair.
The evaluated gcd divides the resultant, since for any
prime at least one coordinate is a unit. The constant quotient
case satisfies the same conclusion directly. Gauss's lemma
makes the quotients by the primitive \(G_j\) integral.
Therefore \(\gcd(F(s,t),E_j(s,t))\) divides
\(|G_j(s,t)| R_j\). Since every \(G_j\) divides \(G\),
\[
                        \mathcal D(s,t)\mid |G(s,t)|\mathcal R.
                                                               \tag{AS1c}
\]
All formal divisor heights are subpower by the bounded degrees
and the subpower coefficient heights of the input. For example,
dehomogenize in an affine chart and bound the roots by the
usual polynomial root bound. The leading coefficient of an
integral factor is bounded by that of the original polynomial,
and the elementary symmetric functions of its roots give a
subpower bound for its other coefficients. A homogeneous
coordinate factor lost in the chart has height one.
For fixed list length each \(R_j\) is subpower by its determinant
formula. Equations (AS1a)–(AS1c) give, with
\(b=1-(\rho+\eta)/3>0\),
\[
                     \operatorname{rad}(F(s,t))
                                 \le H^{S-b+o(1)}.          \tag{AS1d}
\]

### Counting parameters with (AS1d)

Put \(\vartheta=b/(2\ell+1)\). If
\(\prod_i|F_i(s,t)|<H^{S-\vartheta}\), choose the real chart
with denominator of absolute value at least \(H\). The
dehomogenized product is separable and has subpower height.
The affine roots are mutually separated by \(H^{-o(1)}\),
by the nonzero integral discriminant, the coefficient bounds
and the bound on their absolute values. The value in the
chart has magnitude \(O(H^{-\vartheta})\). Root factorization
first puts the parameter within \(H^{-\vartheta/S+o(1)}\) of
one root. All the other factors are bounded below by
\(H^{-o(1)}\), so the distance to that root is
\(H^{-\vartheta+o(1)}\). The close root is real by the same
separation bound applied to conjugates. Counting the
integer numerator and denominator at this distance from
one of the bounded number of roots costs
\(O(H^{2-\vartheta+o(1)}+H)\).

Otherwise use the homogeneous resultants of the distinct
\(F_i,F_j\). Common prime factors of their values at a
primitive pair divide those resultants. Therefore
\[
 \prod_i \operatorname{rad}(F_i(s,t))
     \le \operatorname{rad}(F(s,t))
                \prod_{i<j}|\operatorname{Res}(F_i,F_j)|
     \le H^{S-b+o(1)} .
\]
At least one \(F_i(s,t)\) has
\[
    |F_i(s,t)|/\operatorname{rad}(F_i(s,t))
                             \ge H^{(b-\vartheta)/\ell-o(1)} .
\]
The largest square divisor of an integer is at least its
absolute value divided by its radical, prime by prime.
Thus this \(F_i(s,t)\) is divisible by \(r^2\) with
\(r\ge H^{\vartheta-o(1)}\).

For \(r\le H^{7/4}\), the congruence \(F_i(s,t)\equiv0\bmod r^2\)
in primitive pairs gives \(H^{o(1)}\) projective congruence
classes. In either unit-coordinate chart, the derivative
resultant of the fixed-degree separable dehomogenization is
nonzero of subpower size. At a prime power above twice the
valuation of that resultant, each root lifts to a simple
root over \(\mathbb Q_p\) and at most the resultant valuation
of its last digits are free. At smaller depths use the
trivial prime-power bound. CRT, the two chart choices per
prime and the divisor bound give the stated subpower count.
A constant chart can be handled by its nonzero constant
value. In each projective class the cross determinant of
two parameter vectors is divisible by \(r^2\). Separate
opposite signs and use angular spacing \(\gg r^2/H^2\)
for distinct primitive vectors. Summing
\(H^{o(1)}(1+H^2/r^2)\) over this range costs
\(H^{2-\vartheta+o(1)}+H^{7/4+o(1)}\).

For \(r>H^{7/4}\) only degree four can occur. Set
\(k=F_i(s,t)/r^2\ne0\), so \(|k|\le H^{1/2+o(1)}\).
Writing \(F_i(s,t)=Q_i(s^2,st,t^2)\) with \(Q_i\) quadratic,
map each pair to \([s^2:st:t^2:r]\) on
\[
                       UW=V^2,\qquad k Z^2=Q_i(U,V,W).
\]
In the chart \(W=1\) the equation is \(kZ^2=F_i(V,1)\)
and is geometrically irreducible since \(F_i(V,1)\) is
separable of degree four. The hyperplane \(W=0\) contributes
only finitely many points. The complete intersection has
a geometrically irreducible curve of degree four. The
projective height is \(H^{2+o(1)}\), so Heath-Brown's
coefficient-uniform curve estimate gives \(H^{1+o(1)}\)
points for each \(k\). Primitive pairs are recovered up to
sign from their projective ratio. This branch costs
\(H^{3/2+o(1)}\) and proves the asserted count, since
\(\vartheta\le1/7\).
[Heath-Brown, *The density of rational points in curves and surfaces*, Theorem 5](https://arxiv.org/abs/math/0405392)
· [Earlier depth-three parameter count](ABC_POLYNOMIAL_SUPPORT_DEPTH.md#extension-through-degree-four)

## AS2. A growing list can give full depth on false candidates

Fix integers \(M\ge e\ge3\). Use the fixed binary forms
\[
                       A=s^e,\qquad B=t^e,\qquad C=(s+t)^e .
\]
There are \(H^{2-o(1)}\) positive primitive parameter pairs
with \(s+t\le H\) and a single list, of cardinality and coefficient
height \(H^{o(1)}\) at this scale, for which
\[
     \mathcal D(s,t)=A(s,t)B(s,t)C(s,t)
                   \quad\text{and}\quad A(s,t)+B(s,t)<C(s,t).
                                                               \tag{AS2a}
\]
No residual in the list has value zero at a retained pair.
The actual radical budget is \(\lambda=3/e\),
\(\operatorname{rad}(ABC)<C^\lambda\), and
\(\log_C\operatorname{rad}(W)=\lambda/3+o(1)\) for each
\(W=A,B,C\). The count of distinct minimally repaired target
values at cutoff \(Y=H^M\) is at least \(Y^{1/M-o(1)}\).

Put \(L=\log H\) and \(y=\exp(\sqrt L)\).
Harper's estimate gives \(\gg \Psi(H,y)^3/H\) positive solutions
\(s+t=w\le H\) with all prime factors of \(stw\) at most \(y\).
Here \(\Psi(H,y)\) counts \(y\)-smooth integers at most \(H\).
Already the products of \(\lfloor\sqrt L\rfloor\) distinct primes
in \([y/2,y]\) give \(\Psi(H,y)\ge H\exp(-O(\sqrt L\log L))\)
by the prime number theorem and the binomial count. Thus
the count is \(H^2\exp(-O(\sqrt L\log L))\).
[Harper, *Minor arcs, mean values, and restriction theory for exponential sums over smooth numbers*, Corollary 1](https://arxiv.org/abs/1408.1662)

Discard rows with \(\gcd(s,t)>\exp(L^{2/3})\). Summing
\((H/d)^2\) over integers \(d\) above that threshold proves
a loss \(O(H^2\exp(-L^{2/3}))\). Divide each remaining row by
its gcd, preserving smoothness and the additive identity.
At most \(\exp(L^{2/3})\) rows normalize to one primitive
row, so at least \(H^2\exp(-O(L^{2/3}))\) distinct
primitive rows remain.

Set \(U=\exp(L^{3/4})\). Discard rows where a base value
\(s,t,w\) is smaller than \(H/U\), has radical smaller
than \(H/U^2\), or has a full prime power larger than \(U\).
There are \(O(H/U)\) small values. For a value at least
\(H/U\) with radical at most \(H/U^2\), a square divisor
has square root at least \(U^{1/2}\). The number of those
values is \(O(H/U^{1/2})\) by summing multiples of the
squares. For each \(p\le y\), base values with a full
\(p\)-power larger than \(U\) are multiples of the least
such power. They number \(O(H/U)\) for that prime.
There are \(O(yH/U)\) exceptional values across those
primes. Each exceptional value in a prescribed role
accounts for at most \(H\) rows. All three deletions are
lower order.

Let \(T_H=\lceil U^e\rceil\) and use the residuals
\[
            A+B-jC,\qquad jA+B-C,\qquad 1\le j\le T_H.
\]
Each list entry has all three coefficients nonzero, and
the three original forms are linearly independent for
\(e\ge3\). For a retained pair and \(p\mid A\), choose
\(j\equiv B C^{-1}\bmod p^{v_p(A)}\) in the first list.
For \(p\mid B\) choose \(j\equiv A C^{-1}\bmod p^{v_p(B)}\)
in the same list. For \(p\mid C\) choose
\(j\equiv -B A^{-1}\bmod p^{v_p(C)}\) in the second list.
All three chosen classes are units by primitivity.
The full prime powers of \(A,B,C\) are at most \(U^e\),
so the list contains all the indicated representatives.
This proves (AS2a)'s gcd identity. To exclude zero
evaluations entirely, discard the zeros of all residuals.
Each nonzero homogeneous binary form provides \(O_e(H)\)
integer pairs in the positive height box, by fixing one
positive parameter and solving a nonzero univariate equation.
The cost \(O_e(T_H H)\) is lower order.

Continuing to write \(w=s+t\), the strict inequality in
(AS2a) follows from \(s^e+t^e<(s+t)^e\).
The radical caps follow from the retained base sizes and
radicals, including \(\operatorname{rad}(ABC)\le stw<w^3\).
The actual least repair satisfies
\[
     \operatorname{rad}(w)^M\le H_M(w^e)
                  =\prod_{p\mid w}p^{\max(M,e v_p(w))}
                  \le w^M\le Y.
\]
There are \(H^{1-o(1)}\) distinct \(w\), since there are
at most \(H\) representations of any \(w\) as \(s+t\).
The fixed-support bound gives at most \(H^{o(1)}\) choices
of \(w\le H\) per radical. Hence \(H^{1-o(1)}\) different
radicals occur among the false target representatives,
and different radicals give different repaired values.
[Fixed-support bound and exact support occurrence](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-reduction-to-occurring-support-unions)

The counterexample uses per-prime choice of the linear residual.
It concerns relaxed acceptance with no true additive identity,
despite passing at the complete original depths in those tests.
The fixed-list saving applies to bounded-format polynomial
carriers. The occurrence problem for unrestricted primitive
integer triples still requires a separate argument.
