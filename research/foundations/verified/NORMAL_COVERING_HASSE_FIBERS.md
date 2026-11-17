# Normal coverings and the arithmetic components of Hasse-failing fibers

This note puts the repository's Hasse-failing complete Keller fibers into
the normal-covering framework for intersective polynomials.  The main bound
is intrinsic to the finite etale fiber; Keller realization then transports
the same arithmetic object without changing its component fields.

## 1. Arithmetic components

Let \(K\) be a number field and let \(X\) be a finite etale \(K\)-scheme.
Write its canonical decomposition into \(K\)-connected components as

\[
 X=\coprod_{i=1}^{r}\operatorname{Spec}K_i.             \tag{1.1}
\]

Here \(r\) is the **arithmetic component count**.  It is not the number of
geometric components: after base change to an algebraic closure, a
rank-\(N\) finite etale scheme has \(N\) geometric points.

Choose a finite Galois extension \(L/K\) containing the normal closures of
all \(K_i/K\), minimal with this property, and put

\[
 G=\operatorname{Gal}(L/K).
\]

After choosing one embedding of each \(K_i\) into \(L\), put

\[
 H_i=\operatorname{Gal}(L/K_i).
\]

The \(G\)-set of geometric points of \(X\) is

\[
 \Omega=\coprod_{i=1}^{r}G/H_i,                         \tag{1.2}
\]

so

\[
 [K_i:K]=[G:H_i],\qquad
 \operatorname{rank}_K X=\sum_{i=1}^{r}[G:H_i].         \tag{1.3}
\]

Minimality of \(L\) says that this action is faithful:

\[
 \bigcap_{i=1}^{r}\operatorname{core}_G(H_i)=1.         \tag{1.4}
\]

## 2. The normal-covering bound

Recall that a normal covering of a finite noncyclic group \(G\) is a family
of proper subgroups whose conjugates cover \(G\).  Its minimum cardinality
is the normal covering number \(\gamma(G)\).

### Theorem

Assume

\[
 X(K)=\varnothing,\qquad X(K_v)\ne\varnothing
\quad\text{for every place }v\text{ of }K.              \tag{2.1}
\]

Then the component stabilizers form a faithful normal covering:

\[
 G=\bigcup_{i=1}^{r}\ \bigcup_{g\in G}gH_i g^{-1},
 \qquad
 \bigcap_{i=1}^{r}\operatorname{core}_G(H_i)=1.         \tag{2.2}
\]

Consequently \(G\) is noncyclic and

\[
 \boxed{r\ge\gamma(G).}                                 \tag{2.3}
\]

Since (2.1) also makes every component degree at least two,

\[
 \boxed{\gamma(G)\le r\le
 \frac{\operatorname{rank}_K X}{2}.}                    \tag{2.4}
\]

### Proof

There is no \(K\)-point precisely when no \(K_i\) equals \(K\).  Hence each
\(H_i\) is proper.

For every \(g\in G\), Chebotarev supplies an unramified finite place \(v\)
whose Frobenius conjugacy class is the conjugacy class of \(g\).  A
\(K_v\)-point of \(X\) is a geometric point fixed by a decomposition group
at \(v\), hence in particular fixed by a Frobenius generator.  In the
coset description (1.2), \(g\) therefore belongs to a conjugate of some
\(H_i\).  This proves the covering identity.  Equation (1.4) proves
faithfulness, and (2.3) follows from the definition of \(\gamma(G)\).
Finally, every \([G:H_i]\ge2\), which proves (2.4).

This proof is the finite-etale formulation of the factor-stabilizer
argument for intersective polynomials.  It applies to every complete
regular Keller fiber, independently of the construction used to realize
that fiber.

## 3. The ramified-prime condition

The normal covering controls unramified primes but is not sufficient for
intersectivity.  Let \(w\) be a place of \(L\) over a finite place \(v\) of
\(K\), and let \(D_w\le G\) be its decomposition group.  The exact local
condition is

\[
 X(K_v)\ne\varnothing
 \quad\Longleftrightarrow\quad
 D_w\le gH_i g^{-1}
 \text{ for some }i\text{ and }g\in G.                 \tag{3.1}
\]

At an unramified place \(D_w\) is cyclic and the normal covering handles
its Frobenius generator.  At ramified places (3.1) is an additional finite
audit.  For an irreducible polynomial factor, the equivalent field-level
test is that some prime in the corresponding root field has ramification
degree and inertia degree both equal to one.

This is the division between the group-theoretic and arithmetic filters in
Banks' extension of Sonn's criterion:

- Nicolas Banks,
  [*Classification Results for Intersective Polynomials With No Integral Roots*](https://uwspace.uwaterloo.ca/items/baef0bb0-3712-4117-81b5-c068944ae100),
  PhD thesis, University of Waterloo, 2025;
- the accompanying
  [SageMath and GAP implementations](https://github.com/N2Banks/Intersective-Polynomials-Algorithms).

## 4. A degree-sensitive covering cost

The normal covering number remembers the number of arithmetic components,
but not their degrees.  For the search programme it is useful to introduce
the following repository terminology:

\[
\delta_{\mathrm{fc}}(G)=
\min_{\mathcal H}
\sum_{H\in\mathcal H}[G:H],                             \tag{4.1}
\]

where \(\mathcal H\) ranges over normal coverings by proper subgroups such
that

\[
\bigcap_{H\in\mathcal H}\operatorname{core}_G(H)=1.
\]

Call \(\delta_{\mathrm{fc}}(G)\) the **faithful normal-covering degree**.
This name is not asserted to be standard.  Equations (1.3) and (2.2) give

\[
 \operatorname{rank}_K X\ge\delta_{\mathrm{fc}}(G).     \tag{4.2}
\]

More precisely, the component-degree multiset of \(X\) is the subgroup
index multiset

\[
 \{[G:H_1],\ldots,[G:H_r]\}.                            \tag{4.3}
\]

Thus a low-degree search should enumerate faithful normal covers of the
required index sum before searching polynomial coefficients.  A
ramification-decorated refinement of (4.1) may further restrict the
covering families by requiring (3.1) for prescribed decomposition groups.

For \(S_n\) and \(A_n\), Eberhard and Mellon prove linear asymptotic bounds
for \(\gamma(G)\), exact formulas on several infinite subsequences, and a
translation of the remaining cases into additive-combinatorial avoidance
problems:

- Sean Eberhard and Connor Mellon,
  [*Normal covering numbers for \(S_n\) and \(A_n\) and additive combinatorics*](https://doi.org/10.1112/blms.70154),
  *Bulletin of the London Mathematical Society* 57 (2025), 3307--3325.

The parameter \(n\) in \(S_n\) or \(A_n\) is not automatically the total
fiber rank in (1.3); the actual Keller-fiber rank is the sum of subgroup
indices.

## 5. Low-degree candidate data

Banks' Table C.1 supplies necessary factorization-shape and abstract-group
candidates in total degrees five through ten.  It is not an
if-and-only-if coefficient classification: several rows have no displayed
example, and every candidate still requires the ramified-prime audit.

The versioned transcription

[`arithmetic/banks_degree_5_10_candidates.json`](../arithmetic/banks_degree_5_10_candidates.json)

records 163 group rows across the following shapes:

| degree | possible factorization shapes |
|---|---|
| \(5\) | \((2,3)\) |
| \(6\) | \((2,2,2)\) |
| \(7\) | \((2,2,3),(2,5),(3,4)\) |
| \(8\) | \((2,2,2,2),(2,2,4),(2,3,3)\) |
| \(9\) | \((2,2,2,3),(2,2,5),(2,3,4),(2,7),(4,5)\) |
| \(10\) | \((2,2,2,2,2),(2,2,2,4),(2,2,3,3),(2,2,6),(2,3,5),(2,4,4),(3,3,4),(3,7),(4,6)\) |

The data record the thesis PDF hash, the pinned commit of Banks' code, and
whether Table C.1 displays an example.  The verifier checks the schema,
shapes, group labels, row counts, and necessary-only status; it does not
promote unverified thesis examples to repository theorems.

[`scripts/normal_covering_certificate.g`](../scripts/normal_covering_certificate.g)
is the GAP front end for general finite permutation groups.  It computes the
component stabilizers, their conjugate union, their common core, and the
exact normal covering number from conjugacy classes of maximal subgroups.
The dependency-free Python implementation in
[`jcsearch/normal_covering.py`](../jcsearch/normal_covering.py) independently
replays the checked-in small certificates without requiring a GAP
installation.

## 6. Exact covering-minimal examples

### The quintic \(S_3\) fiber

For

\[
 f_5(T)=(T^3-19)(T^2+T+1),                             \tag{6.1}
\]

the splitting group is \(S_3\).  The quadratic component is the sign orbit
and the cubic component is the natural three-point orbit.  Their
stabilizers have indices two and three and cover \(S_3\); their common core
is trivial.  Exhaustive subgroup enumeration gives

\[
 r=2=\gamma(S_3),\qquad \operatorname{rank}X=5.         \tag{6.2}
\]

The ramified primes are \(3\) and \(19\).  Strong Hensel at \(-2\) in the
cubic covers \(3\), and the simple quadratic root \(7\) covers \(19\).
This replays the arithmetic fiber used in
[the minimal degree-five Keller construction](MINIMAL_HASSE_PRINCIPLE_KELLER_FIBER.md).

### A degree-six \(C_2^2\) Keller fiber

Put

\[
 f_6(T)=(T^2-2)(T^2-17)(T^2-34).                       \tag{6.3}
\]

Its splitting field is \(\mathbb Q(\sqrt2,\sqrt{17})\), with group
\(C_2^2\).  The three quadratic component stabilizers are the three
order-two subgroups.  They cover \(C_2^2\), have trivial common
intersection, and exhaustive subgroup enumeration gives

\[
 r=3=\gamma(C_2^2),\qquad \operatorname{rank}X=6.       \tag{6.4}
\]

At \(2\), the factor \(T^2-17\) satisfies the strong Hensel inequality at
\(T=1\):

\[
 v_2(1-17)=4>2v_2(2)=2.
\]

At \(17\), \(T=6\) is a simple root of \(T^2-2\) modulo \(17\).  The normal
cover handles every unramified prime, and \(f_6\) has real roots and no
rational root.

The public quadratic-gauge compiler uses the translation \(T=1+S\).  It
produces the seed

\[
 G(S)=S^6+6S^5-38S^4-192S^3+377S^2+1154S             \tag{6.5}
\]

and the determinant-one target

\[
 \boxed{\left(1,0,\frac{528}{577}\right)}.              \tag{6.6}
\]

The complete inverse polynomial is

\[
\begin{aligned}
f_6(1+S)
 &=(S^2+2S-1)(S^2+2S-16)(S^2+2S-33),                 \tag{6.7}
\end{aligned}
\]

and the compiled coordinate degrees are

\[
 (7,38,36).
\]

The exact checker expands the map, proves its Jacobian determinant is one,
and verifies (6.7).  Therefore (6.6) is a degree-six complete regular
Keller fiber with points over \(\mathbb R\) and every \(\mathbb Q_p\), no
rational point, splitting group \(C_2^2\), and covering-minimal arithmetic
component count three.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_normal_covering_certificates.py
python3 scripts/verify_banks_degree_5_10_candidates.py
.venv/bin/python scripts/verify_degree_six_normal_cover_keller.py
```

The first command independently enumerates the finite groups and all their
subgroups for the \(S_3\) and \(C_2^2\) certificates.  It verifies the
normal covers, common cores, exact covering numbers, polynomial
factorizations, discriminants, real witnesses, and ramified-prime Hensel
witnesses.  The second validates the pinned necessary-candidate
transcription.  The third compiles and checks the degree-six
determinant-one Keller map.
