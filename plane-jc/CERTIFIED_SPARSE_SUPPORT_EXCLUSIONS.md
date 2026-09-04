# Certified sparse support exclusions for normalized JC(2)

## Status and scope

Let \(k\) be a field of characteristic zero and write a tangent-to-identity
plane map as

\[
F=(x+P,\ y+Q),\qquad P,Q\in (x,y)^2.
\]

Count support **with coordinate multiplicity**:

\[
\sigma(F)=|\operatorname{Supp}(P)|+|\operatorname{Supp}(Q)|.
\]

The exact Keller equation is

\[
E(P,Q):=P_x+Q_y+P_xQ_y-P_yQ_x=0. \tag{1}
\]

This note proves and certifies the following fixed-coordinate statement.

> **Sparse normalized exclusion.**
>
> 1. If one of \(P,Q\) is a monomial and the other has at most five
>    monomials, then (1) is inconsistent on the exact-support chart, except
>    for
>    \[
>    P=a y^m,\qquad
>    Q=b x^2+cxy^m+d y^{2m},\qquad m\geq 2,
>    \]
>    its transpose.  Coefficient-zero loci are delegated to smaller exact
>    supports; only the triangular boundary charts survive.  On the
>    exceptional exact chart the coefficient ideal forces
>    \(c=2ab\) and \(d=a^2b\); the map is a composition of two monomial
>    shears and is polynomially invertible.
>    With four monomials on the other side, the unique survivor is
>    \[
>    P=ay^m,\qquad Q=b(x+ay^m)^3,
>    \]
>    the quartic chain \(Q=b(x+ay^m)^4\), and their transposes.
> 2. If \(P\) and \(Q\) each have exactly two monomials, in arbitrary
>    degrees, then the exact-support Keller ideal is the unit ideal.
> 3. If the support split is \(2+3\) or \(3+2\), in arbitrary degrees, then
>    the exact-support Keller ideal is the unit ideal.
> 4. The \(2+4\) and \(4+2\) exact-support Keller ideals are also unit
>    ideals.  For \(3+3\), the only non-singleton exponent chart has both
>    supports equal to \(\{x^2,xy,y^2\}\); its coefficient ideal forces a
>    directional quadratic shear.
> 5. Consequently every normalized Keller map with
>    \(\sigma(F)\leq 6\) is an automorphism, without a degree bound.

The statement is not invariant under affine changes of source coordinates:
translations and nonmonomial linear changes can enlarge or shrink support.
It does not improve the universal larger-coordinate frontier \(125\), and
it is not a low-degree counterexample search.  It gives a certified
fixed-coordinate support-cardinality frontier at seven nonlinear monomial
occurrences and a reusable certificate format.

## 1. Exact-support coefficient ideals

For fixed finite supports \(S_P,S_Q\), attach one coefficient variable to
each monomial and let \(A\) be their product.  The exact-support chart is
encoded by the Rabinowitsch equation

\[
zA-1=0.
\]

Let \(I_{S_P,S_Q}\) be generated over \(\mathbb Q\) by all coefficients of
\(E(P,Q)\) and by \(zA-1\).  If one monomial of \(E\) receives the single
contribution \(nM\), where \(n\in\mathbb Z\setminus\{0\}\) and \(M\) is a
coefficient monomial dividing \(A\), then

\[
\frac{z(A/M)}{n}(nM)-(zA-1)=1. \tag{2}
\]

Thus (2) is an explicit characteristic-zero unit-ideal certificate.  No
claim rests merely on a Gröbner basis containing \(1\).

Boundary charts are not lost: setting one of the displayed coefficients to
zero gives a smaller support and is handled by its own row.  If one of
\(P,Q\) is empty, (1) forces the other to be a polynomial in the transverse
variable, giving a triangular automorphism.

## 2. Arbitrary-degree singleton-versus-five classification

Put

\[
P=a x^r y^s,\qquad
Q=\sum_{j=1}^q b_jx^{u_j}y^{v_j},\qquad q\leq5,
\]

where all displayed monomials have total degree at least two and all
coefficients are nonzero.  Set \(A=(r,s)\), \(B_j=(u_j,v_j)\), and
\(\Delta_j=rv_j-su_j\).  Equation (1) has the contributions

\[
\begin{aligned}
&ar\,x^{r-1}y^s,\\
&b_jv_j\,x^{u_j}y^{v_j-1},\\
&ab_j\Delta_j\,
  x^{r+u_j-1}y^{s+v_j-1}.
\end{aligned} \tag{3}
\]

The bracket exponents in the last row are distinct.  A nonzero bracket term
cannot meet the \(P_x\) term: equality would force \(B_j=(0,1)\), which is
linear.  Therefore, unless (2) already applies, every
\(\Delta_j\neq0\) forces another \(Q\)-exponent

\[
B_h=B_j+(r-1,s). \tag{4}
\]

This turns the support into directed chains.

### The case \(r>0\)

The nonzero \(P_x\) term must meet a \(Q_y\) term.  The required exponent is

\[
B_*=A+(-1,1).
\]

It has

\[
\det(A,B_*)=r+s\neq0.
\]

Relation (4) therefore forces \(B_1=B_*+(r-1,s)\).  Its determinant with
\(A\) is \(r+s+s\neq0\), so it forces \(B_2\).  Inductively,
\(\det(A,B_i)=r+s+is\neq0\), so \(B_4\) forces a sixth distinct
\(Q\)-monomial.  This contradicts \(q\leq5\).

### The case \(r=0\)

Write \(A=(0,m)\), with \(m\geq2\).  Relation (4) is

\[
(u,v)\longmapsto(u-1,v+m).
\]

A root of a chain has no incoming bracket term.  Its \(Q_y\) contribution
would then be a singleton unless \(v=0\).  Every node with \(u>0\) has
\(\Delta=-mu\neq0\) and must have a successor; a terminal node has \(u=0\).
Hence a complete chain is

\[
(u,0),(u-1,m),\ldots,(0,um).
\]

The cases \(u=0\) and \(u=1\) begin with a constant or linear monomial.
For \(q=3\), the unique admissible chain has \(u=2\):

\[
(2,0),\quad(1,m),\quad(0,2m). \tag{5}
\]

Substitution into (1) gives

\[
E=m(c-2ab)xy^{m-1}+m(2d-ac)y^{2m-1}.
\]

The exact coefficient ideal is therefore

\[
(c-2ab,\ 2d-ac),
\]

and its reduced Gröbner consequences include

\[
c=2ab,\qquad d=a^2b.
\]

With \(X=x+ay^m\), the forced map is

\[
X=x+ay^m,\qquad Y=y+bX^2.
\]

It has the two-sided polynomial inverse

\[
y=Y-bX^2,\qquad x=X-a(Y-bX^2)^m. \tag{6}
\]

Swapping both source and target coordinates gives the transposed class.
The Singular replay computes the saturated coefficient Gröbner basis and
checks both compositions in the representative \(m=5\) member; the argument
above proves the exponent-uniform statement.

Every admissible component has at least three nodes, and there is only one
chain with a given root \(u\).  Thus four nodes cannot split into multiple
admissible components.  For \(q=4\), the unique chain instead has \(u=3\):

\[
(3,0),\quad(2,m),\quad(1,2m),\quad(0,3m).
\]

Write

\[
Q=bx^3+cx^2y^m+dxy^{2m}+ey^{3m}.
\]

Then

\[
\begin{aligned}
E={}&m(c-3ab)x^2y^{m-1}\\
   &+2m(d-ac)xy^{2m-1}\\
   &+m(3e-ad)y^{3m-1}.
\end{aligned}
\]

The coefficient ideal forces

\[
c=3ab,\qquad d=3a^2b,\qquad e=a^3b,
\]

so, with \(X=x+ay^m\),

\[
X=x+ay^m,\qquad Y=y+bX^3. \tag{7}
\]

Its two-sided inverse is

\[
y=Y-bX^3,\qquad x=X-a(Y-bX^3)^m. \tag{8}
\]

The same Singular replay checks the saturated Gröbner consequences and both
inverse compositions.

For \(q=5\), the same argument leaves only the chain beginning at \(u=4\):

\[
(4,0),\quad(3,m),\quad(2,2m),\quad(1,3m),\quad(0,4m).
\]

Writing

\[
Q=bx^4+cx^3y^m+dx^2y^{2m}+exy^{3m}+fy^{4m},
\]

the coefficient ideal of (1) is

\[
(c-4ab,\ 2d-3ac,\ 3e-2ad,\ 4f-ae).
\]

After saturation by \(abcdef\), its Gröbner consequences are

\[
c=4ab,\qquad d=6a^2b,\qquad e=4a^3b,\qquad f=a^4b.
\]

Thus the chart is exactly the quartic shear chain

\[
X=x+ay^m,\qquad Y=y+bX^4,
\]

with inverse

\[
y=Y-bX^4,\qquad x=X-a(Y-bX^4)^m. \tag{9}
\]

The Singular replay checks these saturated relations and both inverse
compositions at \(m=5\).  Hence the singleton-versus-five classification is
arbitrary-degree; it is not inferred from a bounded exponent search.

## 3. Two broad ray classes

These classes are useful controls because their supports can be arbitrarily
large.

### Separated axes

For \(P=p(y)\) and \(Q=q(x)\),

\[
E=-p'(y)q'(x).
\]

If both exact supports are nonempty, each product of one coefficient of
\(p\) and one of \(q\) occurs at its own exponent and gives (2).  If one
side is zero, the map is triangular.

### One common monomial ray

Let \(h=x^ry^s\), \(P=p(h)\), and \(Q=q(h)\).  Their bracket vanishes.
If \(r,s>0\), the \(P_x\) and \(Q_y\) supports are disjoint and every
nonzero contribution is a singleton.  If \(r=0\), equation (1) forces
\(Q=0\) and leaves the triangular shear \(x\mapsto x+p(y^s)\).  If \(s=0\),
it forces \(P=0\) and leaves the transposed triangular shear.

## 4. Arbitrary-degree balanced \(2+2\) exclusion

Write

\[
\begin{aligned}
P&=a x^{A_x}y^{A_y}+b x^{B_x}y^{B_y},\\
Q&=c x^{C_x}y^{C_y}+d x^{D_x}y^{D_y},
\end{aligned}
\]

where \(A,B,C,D\in\mathbb N^2\), every vector has coordinate sum at least
two, \(A\neq B\), and \(C\neq D\).  There are at most eight contributions
to (1), in the fixed order

\[
\begin{array}{c|c|c}
i&\text{source}&\text{exponent}\\ \hline
0&P_{A,x}&A-(1,0)\\
1&P_{B,x}&B-(1,0)\\
2&Q_{C,y}&C-(0,1)\\
3&Q_{D,y}&D-(0,1)\\
4&[A,C]&A+C-(1,1)\\
5&[A,D]&A+D-(1,1)\\
6&[B,C]&B+C-(1,1)\\
7&[B,D]&B+D-(1,1).
\end{array} \tag{10}
\]

The first four are present exactly when
\(A_x,B_x,C_y,D_y\) are positive.  The last four are present exactly when
the corresponding determinant

\[
\Delta_{AC},\Delta_{AD},\Delta_{BC},\Delta_{BD}
\]

is nonzero.

Assume no coefficient of (1) is a singleton.  Every present term in (10)
must then share its exponent with another present term.  Encode the eight
presence decisions by one byte, with bit \(i\) corresponding to row \(i\).
There are only \(256\) presence masks.

An exact integer-linear collision sieve gives:

1. \(85\) masks are compatible with nonnegative, distinct nonlinear
   supports and their declared determinant-zero pattern;
2. only the following \(15\) masks can satisfy the linear exponent
   collisions after temporarily forgetting the determinant equations:
   \[
   \mathtt{0f,3f,5f,66,69,6f,7f,95,9a,9f,af,bf,cf,df,ef};
   \]
3. their canonical equality partitions reduce to the following twenty
   rows.

Here `02|13` means \(T_0=T_2\), \(T_1=T_3\), with different blocks having
different exponents.

| mask | collision partitions | forced contradiction |
|---|---|---|
| `0f` | `02|13`, `03|12` | respectively \(\Delta_{AC}=C_x+C_y>0\) or \(\Delta_{BC}=C_x+C_y>0\), but all four brackets are declared zero |
| `3f` | `02|15|34`, `03|14|25` | the declared-zero \(\Delta_{BC}\) is respectively \((3D_x+D_y+1)/2\) or \(3D_x+2D_y-1\), both positive |
| `5f` | `02|14|36`, `06|12|34` | the declared-zero \(\Delta_{AD}\) is respectively \((D_x+3D_y)/3\) or \((2D_x+3D_y)/3\), both positive |
| `66` | `12|56` | the two zero determinants force \(D_x=C_y=1\), making \(A=(0,1)\) linear |
| `69` | `03|56` | the two zero determinants force \(D_x=0,D_y=1\), making \(D\) linear |
| `6f` | `02|13|56` | the declared-zero \(\Delta_{AC}=C_x+C_y>0\) |
| `7f` | `02|134|56` | the declared-zero \(\Delta_{BD}=D_x+D_y>0\) |
| `95` | `02|47` | the two zero determinants force \(D_x=C_y=1\), making \(C=(0,1)\) linear |
| `9a` | `13|47` | the two zero determinants force \(D_x=0,D_y=1\), making \(A\) linear |
| `9f` | `03|12|47` | the declared-zero \(\Delta_{BC}=C_x+C_y>0\) |
| `af` | `03|15|27`, `07|13|25` | the declared-zero \(\Delta_{AC}\) is respectively \(D_x+3D_y-2\) or \(2D_x+3D_y-2\), both positive |
| `bf` | `03|125|47` | the declared-zero \(\Delta_{BC}=2D_x+2D_y-1>0\) |
| `cf` | `06|13|27`, `07|12|36` | the declared-zero \(\Delta_{AC}\) is respectively \(3D_x+2D_y-1\) or \((3D_x+D_y+1)/2\), both positive |
| `df` | `036|12|47` | the declared-zero \(\Delta_{AD}=D_x+D_y>0\) |
| `ef` | `027|13|56` | the declared-zero \(\Delta_{AC}=2D_x+2D_y-1>0\) |

The positivity assertions use only nonnegativity, the declared active
derivatives, and \(D_x+D_y\geq2\).  Thus all twenty rows are impossible.
Every balanced \(2+2\) support has a singleton contribution \(nM\), and
(2) supplies its explicit unit-ideal certificate.

The checker performs this finite proof twice.  First it enumerates the
presence masks and the twenty canonical linear partitions, then checks each
residual row with exact nonlinear integer arithmetic.  Second it submits the
unpartitioned no-singleton formula directly and obtains `unsat`.  Both use
the pinned Z3 `4.15.3` arithmetic engine.  This is an arbitrary-degree
proof: the eight exponent coordinates are unbounded nonnegative integers.

## 5. Arbitrary-degree \(2+3\) and \(3+2\) exclusion

For the split \(2+3\), write the exponents of \(P\) as \(A,B\) and those of
\(Q\) as \(C,D,E\).  Equation (1) has at most eleven contributions:

1. the two \(x\)-derivative terms from \(P\);
2. the three \(y\)-derivative terms from \(Q\); and
3. the six bracket terms indexed by
   \((A,C),(A,D),(A,E),(B,C),(B,D),(B,E)\).

As before, a term is present exactly when its derivative coordinate or
determinant is nonzero.  If no coefficient is a singleton, every present
one of these eleven exponent vectors must equal another present vector.
The verifier constructs this assertion over ten unbounded nonnegative
integer exponent coordinates, together with the five degree-at-least-two
conditions and the pairwise support inequalities.

The exact collision certificate is:

\[
\begin{array}{c|r}
\text{all presence masks}&2048\\
\text{compatible derivative/determinant masks}&321\\
\text{masks surviving the linear collision sieve}&98\\
\text{global nonlinear integer no-singleton system}&\mathrm{unsat}.
\end{array}
\]

The pinned artifact records SHA-256 digests of the ordered lists of all
\(321\) and \(98\) masks.  Since the no-singleton formula is unsatisfiable,
every exact \(2+3\) support has a coefficient \(nM\) occurring alone.
Formula (2), now with \(A=abcde\), makes its coefficient ideal the unit
ideal.  Swapping both coordinates proves the \(3+2\) case.

There are therefore no coefficient-scheme survivors to classify in the
balanced support-five rows.  The only nontriangular exact support-five
families are the cubic shear chain (7) and its transpose.

## 6. Arbitrary-degree support-six classification

The support-six splits, up to transposition, are \(1+5\), \(2+4\), and
\(3+3\).  The first is classified in Section 2: its only no-singleton chart
is the quartic shear (9).

For the other two splits, the verifier constructs the same exact
no-singleton formula used above, now without enumerating Boolean presence
masks.  A \(2+4\) support has at most fourteen contributions: six divergence
terms and eight brackets.  Over the twelve unbounded nonnegative exponent
coordinates, Z3 proves

\[
\boxed{\text{the \(2+4\) global no-singleton formula is unsatisfiable}.}
\]

Thus every \(2+4\) chart contains a singleton coefficient \(nM\), and (2)
is its explicit unit-ideal certificate with coefficient product
\(A=abcdef\).  Transposition gives \(4+2\).

A \(3+3\) support has at most fifteen contributions.  Its global
no-singleton formula is satisfiable, but adding the negation of

\[
\operatorname{Supp}(P)=\operatorname{Supp}(Q)
  =\{(2,0),(1,1),(0,2)\} \tag{12}
\]

makes the exact integer formula unsatisfiable.  Hence (12) is the unique
exponent survivor, up to reordering the three monomials in each coordinate.
This is an arbitrary-degree classification: the solver variables have no
upper bounds.  A separate census through degree six checks all
\(5\,290\,000\) ordered \(3+3\) support pairs and finds the same single
collision support, but that census is only a regression and is not needed
for the theorem.

Write the survivor as

\[
\begin{aligned}
P&=ay^2+bxy+cx^2,\\
Q&=dy^2+exy+fx^2.
\end{aligned}
\]

The exact coefficient ideal of (1), saturated by \(abcdef\), is generated
by

\[
b+2d,\quad 2c+e,\quad bf+2c^2,\quad
2af+bc,\quad 4ac-b^2. \tag{13}
\]

Its Gröbner basis gives

\[
c=\frac{b^2}{4a},\quad d=-\frac b2,\quad
e=-\frac{b^2}{2a},\quad f=-\frac{b^3}{8a^2},
\qquad ab\ne0. \tag{14}
\]

Set \(v=1\), \(u=b/(2a)\), and \(\lambda=a\).  Then

\[
(P,Q)=\lambda(v,-u)(ux+vy)^2.
\]

The linear form \(ux+vy\) is invariant under \(F\), so the nonlinear part
\(H=(P,Q)\) satisfies \(H\circ(\mathrm{id}+H)=H\).  Therefore

\[
F^{-1}=\mathrm{id}-H.
\]

The Singular replay verifies (13)--(14), the exact-support saturation, and
both inverse compositions.  Consequently no support-six counterexample
survives:

\[
\boxed{\sigma(F)\leq6\ \text{and}\ \det JF=1
       \quad\Longrightarrow\quad F\ \text{is an automorphism}.}
\]

## 7. Independent bounded \(2+2\) regression

There are \(88\) nonlinear monomials of total degree \(2,\ldots,12\), hence

\[
\binom{88}{2}=3828
\]

two-term supports for each coordinate and

\[
3828^2=14\,653\,584
\]

ordered \(2+2\) support pairs.  For every pair, the checker expands the at
most eight contributions in (1), groups them by the exponent of \(x,y\),
and selects the lexicographically first coefficient with exactly one
contribution.  Formula (2) is then its exact unit certificate.

The generated artifact records the census parameters, count, certificate
identity, Singular version, and SHA-256 digest of every selected certificate.
The digest input is the concatenation of big-endian signed-short records

```text
(Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,ex,ey,n,mask,source).
```

The pinned result has no survivor.  This is an exhaustive rational
calculation, not a finite-field screen.  It is retained as a large
independent regression of the arbitrary-degree collision proof.

## 8. Reproduction

Replay the pinned artifact and the exceptional Gröbner/inverse certificate:

```bash
.venv/bin/python plane-jc/cas/verify_sparse_support_exclusions.py
```

For cleanup and ledger verification, pin-check the committed JSON and Singular
certificate and validate the arbitrary-degree/bounded-regression distinction
without enumerating either bounded support census or invoking a solver:

```bash
.venv/bin/python plane-jc/cas/verify_sparse_support_exclusions.py \
  --audit-existing-only
```

The balanced census takes about one minute on a laptop.  Run the Singular
piece alone with:

```bash
Singular -q plane-jc/cas/sparse_support_exceptional.sing
```

Intentional regeneration is separate:

```bash
.venv/bin/python plane-jc/cas/verify_sparse_support_exclusions.py --refresh
```

The pinned artifact is
[`artifacts/generated-results/jc2_sparse_support_exclusions.json`](../artifacts/generated-results/jc2_sparse_support_exclusions.json).

## 9. Affine-normalized support complexity

The fixed-coordinate theorem has an immediate affine-invariant
reformulation.  For a plane Keller map \(F\), a point \(p\in k^2\), and
\(L\in\operatorname{GL}_2(k)\), define its tangent-to-identity
normalization

\[
N_{p,L}(z)
  =(JF(p)L)^{-1}\bigl(F(p+Lz)-F(p)\bigr). \tag{11}
\]

Let

\[
\sigma_{\mathrm{aff}}(F)
  =\min_{p,L}\sigma(N_{p,L}).
\]

Every \(N_{p,L}\) is affinely left--right equivalent to \(F\), so it is
polynomially invertible exactly when \(F\) is.  The support-six theorem
therefore gives the coordinate-invariant consequence

\[
\boxed{\text{a noninvertible plane Keller map has }
       \sigma_{\mathrm{aff}}(F)\geq7.}
\]

This does not yet translate into a larger coordinate-degree lower bound.
It supplies the affine-normalization gate that a sparse standard/minimal
pair must pass before any coefficient elimination.

## 10. What this changes, and the next approach

This closes every normalized exact-support stratum with at most six
nonlinear monomial occurrences in arbitrary degree.  The falsification-first
search found two support-six survivor components, but both are certified
automorphisms: the quartic shear chain and the directional quadratic shear.
There is no unresolved support-six coefficient component.

This is the stopping point for sequential support-cardinality escalation.
Repeating the same program for support seven, eight, and nine would increase
the affine lower bound but would not create a bridge to the degree or
geometric frontiers.  The next approach should instead seek a structural
link between affine-normalized support and one of:

1. Newton-polygon edge count and mixed area;
2. geometric degree or the nonproperness curve;
3. the standard/minimal-pair normal form at infinity.

A useful theorem would force a hypothetical minimal counterexample into
affine-normalized support at most six.  Combined with the present exclusion,
that would rule it out.  Without such a bridge, another finite sparse layer
has diminishing mathematical value.  The direct bridge audit in
[`AFFINE_SUPPORT_NEWTON_BRIDGE.md`](AFFINE_SUPPORT_NEWTON_BRIDGE.md) proves
that coarse Newton/boundary data cannot supply such an upper bound and
replaces it with a Kummer-character-resolved lower-band gate.

For context, this fixed-support approach is complementary to geometric
conditions such as Chau's
[plane JC for simple polynomials](https://arxiv.org/abs/0711.3894) and to
coefficient-ideal methods for
[two-variable \(d\)-linear maps](https://arxiv.org/abs/2111.10739).  Neither
external result is used in the certificates above.
