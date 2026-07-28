# Normal coverings and component complexity

## Status and purpose

This is a standalone research memo.  It is not part of the
`common-arithmetic-fibers` manuscript, it is not a canonical proof source,
and it does not change any entry of `MATH_STATUS.json`.

The purpose is to record a possible later connection between:

1. everywhere locally soluble finite-etale Keller fibers;
2. normal coverings of finite groups;
3. the normal covering number \(\gamma(G)\);
4. Nicolas Banks' degree-\(5\)-through-\(10\) classification of possible
   Galois groups of strongly intersective polynomials.

The clean theorem below is proved here, but the tables assembled from Banks
must retain the distinction between a candidate group and a group for which
an explicit strongly intersective polynomial is known.

## 1. Source audit

The main external sources are:

- Sean Eberhard and Connor Mellon,
  [*Normal covering numbers for \(S_n\) and \(A_n\) and additive
  combinatorics*](https://arxiv.org/abs/2410.06999), especially the
  introduction.  They define \(\gamma(G)\), explain its relation to
  intersective polynomials, state the factor-count bound, and use
  \[
    (T^3-19)(T^2+T+1)
  \]
  as their \(S_3\) example.
- Nicolas Banks,
  [*Classification Results for Intersective Polynomials With No Integral
  Roots*](https://uwspace.uwaterloo.ca/items/baef0bb0-3712-4117-81b5-c068944ae100),
  PhD thesis, University of Waterloo, 2025.  The classification theorem is
  Theorem 1.1.3; the degree-\(5\)-through-\(10\) candidates and displayed
  examples are in Table C.1, printed pages 74--82.
- Banks'
  [companion algorithm repository](https://github.com/N2Banks/Intersective-Polynomials-Algorithms),
  which contains the GAP conjugate-covering test, the Sage subdirect-product
  search, and the ramification-degree checker.
- D. Berend and Y. Bilu,
  [*Polynomials with roots modulo every
  integer*](https://doi.org/10.1090/S0002-9939-96-03210-8),
  Proc. Amer. Math. Soc. 124 (1996), 1663--1671.

The relevant internal inputs are:

- [`verified/FINITE_ETALE_KELLER_FIBERS.md`](verified/FINITE_ETALE_KELLER_FIBERS.md),
  for the determinant-one realization of a squarefree polynomial quotient;
- [`verified/MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md`](verified/MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md),
  for the explicit quintic Hasse failure;
- [`papers/common-arithmetic-fibers/sections/01-keller-fibers.tex`](papers/common-arithmetic-fibers/sections/01-keller-fibers.tex),
  for the current manuscript statement of the realization theorem.

## 2. Normal coverings

Let \(G\) be a finite noncyclic group.  A **normal covering** of \(G\) is a
family of proper subgroups \(H_1,\ldots,H_m<G\) such that
\[
  G=\bigcup_{i=1}^m\ \bigcup_{\sigma\in G}\sigma H_i\sigma^{-1}.
  \tag{2.1}
\]
The **normal covering number** \(\gamma(G)\) is the least possible \(m\).

For a finite group, Jordan's derangement theorem implies
\(\gamma(G)>1\).  A cyclic finite group has no normal covering by proper
subgroups: a generator belongs to none of them.

## 3. Component-complexity theorem

### Theorem 3.1

Let \(X\) be a finite etale \(\mathbb Q\)-scheme.  Assume
\[
  X(\mathbb Q)=\varnothing
  \tag{3.1}
\]
and
\[
  X(\mathbb Q_p)\ne\varnothing
  \quad\text{for all but finitely many rational primes }p.
  \tag{3.2}
\]
Write
\[
  X=\coprod_{i=1}^m\operatorname{Spec}K_i
  \tag{3.3}
\]
with the \(K_i/\mathbb Q\) finite separable fields.  Let \(L/\mathbb Q\)
be the minimal splitting field of \(X\), and put
\[
  G=\operatorname{Gal}(L/\mathbb Q).
  \tag{3.4}
\]
Then \(G\) is noncyclic and
\[
  \boxed{m\geq\gamma(G).}
  \tag{3.5}
\]

In particular, the conclusion holds for an everywhere locally soluble
finite-etale scheme without a rational point, and hence for every full
Keller fiber with that local-global behavior.

### Proof

Choose a geometric point \(\omega_i\) on the \(i\)-th connected component
and define
\[
  H_i=\operatorname{Stab}_G(\omega_i).
  \tag{3.6}
\]
The orbit of \(\omega_i\) is \(G/H_i\), and
\[
  [G:H_i]=[K_i:\mathbb Q].
  \tag{3.7}
\]
Condition (3.1) says that no component is \(\operatorname{Spec}\mathbb Q\),
so every \(H_i\) is proper.

Fix \(g\in G\).  Enlarge the finite exceptional set in (3.2) to include
all primes ramified in \(L\).  Chebotarev supplies a prime outside this set
whose Frobenius conjugacy class is the conjugacy class of \(g\).
By (3.2), a \(\mathbb Q_p\)-point of \(X\) exists.  Such a point fixes a
geometric point under the corresponding decomposition group, and hence
under its Frobenius generator.  Therefore a conjugate of \(g\) lies in one
of the \(H_i\).  Since \(g\) was arbitrary, (2.1) holds.

Thus the \(H_i\) form a normal covering.  This rules out cyclic \(G\) and
gives \(m\geq\gamma(G)\).

### Remarks

1. Fullness and the Keller condition are not used in the proof.  They enter
   only because every Keller fiber under discussion is finite etale and the
   realization theorem transfers finite-etale examples to full fibers.
2. Solubility at the real place and at the finitely many ramified primes is
   not needed for (3.5).  It is needed for an actual Hasse failure.
3. The theorem is an unramified necessary condition.  At a ramified prime,
   the stronger Sonn--Banks condition requires the entire decomposition
   group to lie in a conjugate of some \(H_i\).

## 4. The faithful \(G\)-set dictionary

The geometric points of \(X\) form the faithful intransitive \(G\)-set
\[
  \Omega=X(\overline{\mathbb Q})
  \simeq\coprod_{i=1}^mG/H_i.
  \tag{4.1}
\]
The three relevant arithmetic properties become:

\[
\begin{array}{c|c}
\text{arithmetic property}&\text{\(G\)-set property}\\
\hline
X(\mathbb Q)=\varnothing&\Omega^G=\varnothing\\
\text{local points at almost all primes}&
  \Omega^g\ne\varnothing\text{ for every }g\in G\\
L\text{ is the minimal splitting field}&
  \displaystyle\bigcap_i\operatorname{core}_G(H_i)=1.
\end{array}
\tag{4.2}
\]

There are consequently two different complexity measures:

\[
  m=\text{number of orbits}=\text{number of connected components},
  \tag{4.3}
\]
and
\[
  \deg X=|\Omega|=\sum_{i=1}^m[G:H_i].
  \tag{4.4}
\]

The invariant \(\gamma(G)\) controls (4.3), but it does not control the
index sum (4.4).  For later use one could introduce the provisional
weighted invariant
\[
 \delta(G)=
 \min\left\{
   \sum_i[G:H_i]:
   \begin{array}{l}
   H_i<G,\ \{\!H_i\!\}\text{ normally covers }G,\\
   \bigcap_i\operatorname{core}_G(H_i)=1
   \end{array}
 \right\}.
 \tag{4.5}
\]
Then every finite-etale Hasse failure with splitting group \(G\) satisfies
\[
  \#\pi_0(X)\geq\gamma(G),
  \qquad
  \deg X\geq\delta(G).
  \tag{4.6}
\]
Banks' factorization patterns are precisely the small-index data that the
unweighted invariant \(\gamma(G)\) forgets.

The notation \(\delta(G)\) in (4.5) is only a proposal for this project; no
claim is made that it is standard terminology.

## 5. The quintic \(S_3\) certificate

Consider
\[
  P_5(T)=(T^3-19)(T^2+T+1).
  \tag{5.1}
\]
The cubic is irreducible and has discriminant
\[
  \operatorname{disc}(T^3-19)
  =-27\cdot19^2
  =-3(3\cdot19)^2.
  \tag{5.2}
\]
It therefore has splitting group \(S_3\).  Its unique quadratic subfield is
\(\mathbb Q(\sqrt{-3})\), which is the splitting field of \(T^2+T+1\).
Consequently the product (5.1) also has splitting group \(S_3\), rather than
\(S_3\times C_2\).

The quadratic orbit has stabilizer
\[
  H_2=A_3,\qquad [S_3:H_2]=2,
  \tag{5.3}
\]
and the cubic orbit has a root stabilizer
\[
  H_3\simeq C_2,\qquad [S_3:H_3]=3.
  \tag{5.4}
\]
The subgroup \(A_3\) contains the identity and both 3-cycles.  The three
conjugates of \(C_2\) contain the transpositions.  Hence
\[
  S_3=A_3\cup\bigcup_{\sigma\in S_3}\sigma C_2\sigma^{-1}.
  \tag{5.5}
\]
Jordan's theorem excludes a one-subgroup normal covering, so
\[
  \boxed{\gamma(S_3)=2.}
  \tag{5.6}
\]
The quintic is therefore optimal both in rank and in number of components:
\[
  \deg P_5=5,\qquad
  \#\{\text{irreducible factors}\}=2=\gamma(S_3).
  \tag{5.7}
\]

Eberhard--Mellon use this exact example in their introduction and state the
general factor-count inequality in the same discussion.

## 6. What Banks' theorem does and does not classify

Banks defines \(P\in\mathbb Z[T]\) to be strongly intersective when it has
a root modulo every positive integer but has no integer root.  For a monic
polynomial this implies that it has no rational root.

Theorem 1.1.3 of the thesis is a necessary classification:

> If \(P\) is strongly intersective, \(5\leq\deg P\leq10\), and \(G\) is
> its Galois group, then \(G\) occurs in Table C.1.

It is not an if-and-only-if realization theorem.  Banks gives examples for
many entries, but Table C.1 contains dashes, especially in degrees \(9\)
and \(10\).  The thesis explicitly describes filling those gaps as future
work.  Therefore a later Keller table must use at least the statuses

- `candidate`: survives the group-theoretic filters;
- `realized`: an explicit strongly intersective polynomial is supplied;
- `Keller-realized`: the displayed polynomial has been passed through the
  finite-etale realization theorem.

Normal covering is only the first filter.  Banks additionally uses:

1. the possible transitive groups of the irreducible factors;
2. subdirect-product constraints between their splitting fields;
3. faithfulness of the combined root action;
4. ramified-prime decomposition or ramification checks.

## 7. Minimum-component slice in degrees \(5\)--\(10\)

The following table retains only factorization patterns with the smallest
number of irreducible factors at each total degree.  Thus it is a
component-complexity slice of Banks' Table C.1, not a replacement for the
full table.

| degree | minimum \(m\) | minimum factor-orbit patterns | possible \(G\) for those patterns |
|---:|---:|---|---|
| \(5\) | \(2\) | \((2,3)\) | \(S_3\) |
| \(6\) | \(3\) | \((2,2,2)\) | \(C_2^2\) |
| \(7\) | \(2\) | \((2,5)\) | \(D_5,\ F_5\) |
|  |  | \((3,4)\) | \(A_4,\ S_4,\ C_3{:}S_4\) |
| \(8\) | \(3\) | \((2,2,4)\) | \(D_4,\ C_2\times D_4\) |
|  |  | \((2,3,3)\) | \(D_6,\ C_3\times S_3,\ C_3{:}S_3,\ S_3^2,\ C_6{:}S_3\) |
| \(9\) | \(2\) | \((2,7)\) | \(D_7,\ F_7\) |
|  |  | \((4,5)\) | \(F_5,\ C_2\times F_5,\ C_{20}{:}C_4,\ D_{10}{:}C_4,\ A_4{:}F_5\) |
| \(10\) | \(2\) | \((3,7)\) | \(C_7{:}C_3,\ F_7,\ C_{21}{:}C_6\) |
|  |  | \((4,6)\) | 37 candidate group IDs, listed below |

Banks uses \(D_n\) for the dihedral group of order \(2n\), and \(F_n\) for
the Frobenius group of order \(n(n-1)\).

The 37 degree-\(10\), pattern-\((4,6)\) candidate entries have the following
LMFDB abstract-group IDs:

```text
12.3
24.12  24.13
48.30  48.31  48.48  48.49  48.50
72.42  72.43  72.44
96.186 96.187 96.194 96.195 96.197 96.226 96.227 96.229
144.183 144.189
192.1470 192.1472 192.1488 192.1538
216.164
288.1024 288.1026
432.743 432.747 432.748
576.8653 576.8657
864.4669 864.4670
1440.5847
8640.o
```

These IDs are preferable to compressed semidirect-product names in a
machine-facing ledger because Table C.1 contains distinct groups with the
same abbreviated structural notation.

## 8. Low-height displayed certificates

There is no presentation-independent notion of a smallest polynomial:
affine changes of the primitive element alter the coefficients.  For a
finite comparison, define the expanded coefficient height of a monic
integral polynomial by
\[
  H_{\mathrm{coeff}}(P)=\max_j|[T^j]P|.
  \tag{8.1}
\]

Among the actual certificates displayed by Banks for a minimum-component
pattern, the following are low-height representatives.  In each degree,
the displayed choice has the smallest value of (8.1) among the relevant
examples printed in Table C.1.  This is not a global minimality claim, and
an as-yet-unrealized candidate group could eventually yield a smaller
certificate.

| degree | orbit pattern | \(G\) | displayed strongly intersective polynomial \(P_d(T)\) | \(H_{\mathrm{coeff}}(P_d)\) |
|---:|---|---|---|---:|
| \(5\) | \((2,3)\) | \(S_3\) | \((T^3-19)(T^2+T+1)\) | \(19\) |
| \(6\) | \((2,2,2)\) | \(C_2^2\) | \((T^2-2)(T^2-17)(T^2-34)\) | \(1156\) |
| \(7\) | \((3,4)\) | \(S_4\) | \((T^4-T+1)(T^3-4T-1)\) | \(4\) |
| \(8\) | \((2,2,4)\) | \(D_4\) | \((T^4+T^2-2T+1)(T^2-T-4)(T^2+1)\) | \(7\) |
| \(9\) | \((4,5)\) | \(C_2\times F_5\) | \((T^5-T^4+T^2+3T+1)(T^4-T^3-6T^2+T+1)\) | \(18\) |
| \(10\) | \((4,6)\) | \(S_4\) | \((T^6-T^5-T^4+T^3-T^2-T+1)(T^4+T^2-T+1)\) | \(2\) |

Their expanded forms are:

\[
\begin{aligned}
P_5={}&T^5+T^4+T^3-19T^2-19T-19,\\
P_6={}&T^6-53T^4+680T^2-1156,\\
P_7={}&T^7-4T^5-2T^4+T^3+4T^2-3T-1,\\
P_8={}&T^8-T^7-2T^6-4T^5-4T^4+4T^3-5T^2+7T-4,\\
P_9={}&T^9-2T^8-5T^7+8T^6+2T^5-9T^4-18T^3-2T^2+4T+1,\\
P_{10}={}&T^{10}-T^9-T^7-2T^4+T^3+T^2-2T+1.
\end{aligned}
\tag{8.2}
\]

The selected examples are normal-covering optimal:

| degree | selected \(G\) | component count \(m\) | \(\gamma(G)\) |
|---:|---|---:|---:|
| \(5\) | \(S_3\) | \(2\) | \(2\) |
| \(6\) | \(C_2^2\) | \(3\) | \(3\) |
| \(7\) | \(S_4\) | \(2\) | \(2\) |
| \(8\) | \(D_4\) | \(3\) | \(3\) |
| \(9\) | \(C_2\times F_5\) | \(2\) | \(2\) |
| \(10\) | \(S_4\) | \(2\) | \(2\) |

For \(C_2^2\), all three order-two subgroups are required because the group
is abelian.  For \(D_4\) of order eight, the rotation subgroup and the two
reflection classes force three covering types.  The other selected groups
have an exhibited two-subgroup normal covering, and Jordan excludes one.

## 9. Keller targets are coordinate-dependent

The phrase "smallest Keller target" is not intrinsic.

First, if \(F^{-1}(y)\simeq X\), postcomposing \(F\) with translation by
\(-y\) gives a Keller map with the same fiber over the origin.  Thus, across
unrestricted map presentations, every certificate has target
\[
  (0,0,0).
  \tag{9.1}
\]

Second, even if the normalized realization-theorem shape
\[
  y_{P,a}=\left(1,0,-\frac{2P(a)}{P'(a)}\right)
  \tag{9.2}
\]
is retained, an affine change of primitive element makes the target uniform.
Choose \(b\in\mathbb Q\) such that
\[
  P'(b)P'''(b)\ne0.
  \tag{9.3}
\]
Because \(P\) has no rational root, \(P(b)\ne0\).  Put
\[
  r=-\frac{2P(b)}{P'(b)},
  \qquad
  Q(U)=P(b+rU).
  \tag{9.4}
\]
Then
\[
  Q'(0)=rP'(b)\ne0,\qquad
  Q'''(0)=r^3P'''(b)\ne0,
  \tag{9.5}
\]
and
\[
  -\frac{2Q(0)}{Q'(0)}
  =-\frac{2P(b)}{rP'(b)}
  =1.
  \tag{9.6}
\]
The realization theorem applied to \(Q\) at \(0\) therefore gives
\[
  y_{Q,0}=(1,0,1),
  \tag{9.7}
\]
while the affine substitution \(T=b+rU\) gives
\[
  \mathbb Q[U]/(Q)\simeq\mathbb Q[T]/(P).
  \tag{9.8}
\]

For the six representatives in Section 8, exact convenient pairs \((b,r)\)
are:

| degree | \(b\) | \(r=-2P_d(b)/P_d'(b)\) |
|---:|---:|---:|
| \(5\) | \(0\) | \(-2\) |
| \(6\) | \(-2\) | \(195/152\) |
| \(7\) | \(0\) | \(-2/3\) |
| \(8\) | \(0\) | \(8/7\) |
| \(9\) | \(0\) | \(-1/2\) |
| \(10\) | \(0\) | \(1\) |

Thus a later manuscript table should use the column heading
`normalized Keller target`, not `smallest Keller target`.  With the present
realization theorem every row may use \((1,0,1)\), and an unrestricted
target translation moves every row to the origin.

## 10. Exact checks performed for this memo

The following finite algebraic checks were performed over
\(\mathbb Q\):

1. expansion of all six displayed products in (8.2);
2. calculation of their coefficient heights;
3. squarefreeness via \(\gcd(P_d,P_d')=1\);
4. verification that the pairs in Section 9 satisfy
   \(P_d'(b)P_d'''(b)\ne0\);
5. verification of (9.6) for every row.

These are small exact SymPy calculations.  No claim of global
coefficient-height minimality was tested, and no missing candidate row in
Banks' table was realized.

No repository checker, `make check`, or manuscript build is required for
this memo because it changes no theorem source, checker, generated status,
or manuscript input.

## 11. Possible later integration

If this material is eventually folded into the manuscript, the safest
sequence is:

1. add Theorem 3.1, preferably immediately after the finite-etale
   local-global setup;
2. identify the quintic stabilizers \(A_3\) and \(C_2\), proving
   \(2=\gamma(S_3)\);
3. cite Eberhard--Mellon for the motivating normal-covering viewpoint;
4. add only the compact minimum-component table from Section 7;
5. maintain separate `candidate` and `realized` statuses for Banks' rows;
6. call the final column `normalized target`;
7. avoid claiming global polynomial or target minimality without a fixed
   normalization and a separate exhaustive search.

The strongest conceptual summary is
\[
\boxed{
\begin{array}{c}
\text{finite-etale Hasse failure}\\
\Downarrow\\
\text{faithful fixed-point-everywhere intransitive \(G\)-set}\\
\Downarrow\\
\#\pi_0(X)\geq\gamma(G),\quad
\deg X=\sum_i[G:H_i]\\
\Downarrow\\
\text{Banks' low-degree orbit-pattern constraints}\\
\Downarrow\\
\text{explicit full Keller realization.}
\end{array}}
\tag{11.1}
\]

