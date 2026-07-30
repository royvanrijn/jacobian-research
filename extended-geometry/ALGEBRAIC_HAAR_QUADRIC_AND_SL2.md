# Algebraic Haar functionals on the Hopf quadric and `SL(2)`

This note removes compact integration from the `SO(3)` and `SU(2)` moment
reproductions.  The main result is an algebraic theorem over a
characteristic-zero field.  For the statements involving the standard
coordinates \(r_{ij}\) of \(SO_3\), assume additionally that
\(\sqrt{-1}\in k\); alternatively use the split orthogonal group of
\(UV+T^2\).

The exact checkers are

- [`verify_algebraic_quadric_haar.py`](../scripts/verify_algebraic_quadric_haar.py);
- [`verify_algebraic_sl2_haar.py`](../scripts/verify_algebraic_sl2_haar.py).

Their bounded loops are regressions.  The all-exponent proofs are the
recurrences below.

## 1. The quadric functional

Put

\[
A=k[U,V,T]/(UV+T^2-1).
\]

For nonnegative \(a,b,c\), define

\[
\lambda(U^aV^bT^c)=0\quad(a\ne b),
\tag{1.1}
\]

\[
\lambda((UV)^nT^{2j+1})=0,
\tag{1.2}
\]

and

\[
\lambda((UV)^nT^{2j})
=
\frac{2\,4^n n!(2j)!(n+j+1)!}
     {j!(2n+2j+2)!}.
\tag{1.3}
\]

In particular,

\[
\lambda(T^{2j})=\frac1{2j+1},\qquad
\lambda((UV)^n)=\frac{4^n(n!)^2}{(2n+1)!}.
\tag{1.4}
\]

Every ambient monomial is reduced exactly as requested.  If \(a=b\), retain
\((UV)^aT^c\).  If \(a>b\), expand

\[
U^aV^bT^c
=U^{a-b}(1-T^2)^bT^c
=\sum_{q=0}^b(-1)^q\binom bq U^{a-b}T^{c+2q};
\tag{1.5}
\]

the case \(b>a\) is symmetric with \(V\).  Thus the only normal-form families
used by the checker are

\[
U^rT^c,\qquad V^rT^c,\qquad (UV)^nT^c.
\]

### 1.1 Descent to the quotient

Let \(L_{n,c}=\lambda((UV)^nT^c)\).  Formula (1.3) satisfies

\[
L_{n+1,c}+L_{n,c+2}-L_{n,c}=0.
\tag{1.6}
\]

If \(a\ne b\), all three terms obtained by multiplying
\((UV+T^2-1)U^aV^bT^c\) have the same nonzero weight and vanish.  If \(a=b\),
their values are exactly (1.6).  Hence

\[
\lambda((UV+T^2-1)f)=0
\]

for every \(f\in k[U,V,T]\), so (1.1)--(1.3) define a functional on \(A\).

### 1.2 Infinitesimal invariance

Use the three derivations

\[
D_0=U\partial_U-V\partial_V,\qquad
D_+=U\partial_T-2T\partial_V,\qquad
D_-=V\partial_T-2T\partial_U.
\tag{1.7}
\]

They annihilate \(UV+T^2-1\).  Their brackets are

\[
[D_0,D_+]=D_+,\qquad [D_0,D_-]=-D_-,
\qquad [D_+,D_-]=-2D_0.
\]

Thus they span the split form of \(\mathfrak {so}_3\simeq\mathfrak {sl}_2\).
The \(D_0\) identity is immediate from weight.  Applied to a monomial,

\[
D_+(U^aV^bT^c)
=cU^{a+1}V^bT^{c-1}-2bU^aV^{b-1}T^{c+1}.
\tag{1.8}
\]

Both terms have zero \(\lambda\)-value unless \(b=a+1\).  In that remaining
case, \(\lambda(D_+f)=0\) is precisely

\[
cL_{a+1,c-1}=2(a+1)L_{a,c+1},
\tag{1.9}
\]

which follows directly from (1.3); the odd cases are zero.  The \(D_-\)
calculation is symmetric.  Therefore

\[
\lambda(Df)=0\qquad
(D\in\langle D_0,D_+,D_-\rangle,\ f\in A)
\tag{1.10}
\]

is proved without an integral.

### 1.3 Lie invariance uniquely forces the formula

There is also a derivation of (1.1)--(1.3), rather than merely a verification.
The quotient has the vector-space basis

\[
\{T^c:c\ge0\}
\cup\{U^rT^c:r\ge1,c\ge0\}
\cup\{V^rT^c:r\ge1,c\ge0\}.
\tag{1.11}
\]

The equation \(\lambda(D_0f)=0\) kills the last two families.  Write
\(\mu_c=\lambda(T^c)\).  In \(A\),

\[
D_+(VT^c)
=cUVT^{c-1}-2T^{c+1}
=cT^{c-1}-(c+2)T^{c+1}.
\]

Consequently

\[
\mu_1=0,\qquad
\mu_{c+1}=\frac{c}{c+2}\mu_{c-1}.
\tag{1.12}
\]

Normalization \(\mu_0=1\) now gives

\[
\mu_{2j+1}=0,\qquad \mu_{2j}=\frac1{2j+1}.
\tag{1.13}
\]

Finally, \(UV=1-T^2\) gives

\[
L_{n,c}=\lambda((1-T^2)^nT^c),
\]

and the binomial expansion with (1.13) is (1.2)--(1.3).  Hence:

> **Quadric uniqueness theorem.**  
> The functional (1.1)--(1.3) is the unique normalized functional on \(A\)
> annihilated by the infinitesimal \(\mathfrak {so}_3\)-action.

For a connected characteristic-zero group, infinitesimal invariance is
equivalent to invariance under the algebraic group action.  Thus this is
exactly the normalized algebraic invariant, or Reynolds, functional on the
quadric.

## 2. The all-order Mathieu witness

Set

\[
P=(1+U)\bigl(V-(2+U)T^2\bigr),\qquad Q=U.
\tag{2.1}
\]

Using \(V=(1-T^2)/U\) in the Laurent localization gives

\[
P=\frac{1+U}{U}\bigl(1-T^2(1+U)^2\bigr).
\tag{2.2}
\]

On Laurent polynomials in \(U,T\), (1.1) and (1.13) say: take the constant
term in \(U\), kill odd powers of \(T\), and send \(T^{2j}\) to
\((2j+1)^{-1}\).  Define the polynomial

\[
J_m(X)=\sum_{j=0}^m
 \frac{(-1)^j\binom mj}{2j+1}X^{2j+1}.
\tag{2.3}
\]

This definition is algebraic in characteristic zero, and

\[
J_m'(X)=(1-X^2)^m.
\tag{2.4}
\]

The constant-term expansion of (2.2) identifies

\[
\lambda(P^m)=[s^m]\,(1+s)^{m-1}J_m(1+s),
\tag{2.5}
\]

\[
\lambda(UP^m)=[s^{m-1}]\,(1+s)^{m-1}J_m(1+s).
\tag{2.6}
\]

Equation (2.4) has a zero of order \(m\) at \(X=1\).  Therefore every
derivative \(J_m^{(q)}(1)\) with \(1\le q\le m\) is zero.  The right side of
(2.5) is zero, while (2.6) is \(J_m(1)\).

The value of \(J_m(1)\) also needs no analytic integral.  Let
\(\mathcal I(f)\) be the endpoint difference \(F(1)-F(0)\), where \(F'=f\)
and \(F(0)=0\).  This is a purely algebraic functional on \(k[X]\).  Applying
it to

\[
\frac{d}{dX}\bigl(X(1-X^2)^m\bigr)
=(1-X^2)^m-2mX^2(1-X^2)^{m-1}
\]

gives

\[
(2m+1)J_m(1)=2mJ_{m-1}(1),\qquad J_0(1)=1.
\]

Hence

\[
J_m(1)=\frac{4^m(m!)^2}{(2m+1)!}.
\tag{2.7}
\]

Combining (2.5)--(2.7) proves:

> **Algebraic quadric counterexample.**  
> For every \(m\ge1\),
> \[
> \lambda(P^m)=0,\qquad
> \lambda(UP^m)=\frac{4^m(m!)^2}{(2m+1)!}\ne0.
> \]
> Consequently \(\ker\lambda\) is not a Mathieu subspace of \(A\).

This is the requested proposition.  The proof uses only quotient algebra,
derivations, coefficient extraction, and formal polynomial antiderivatives.

## 3. The algebraic `SL(2)` functional

Write

\[
g=\begin{pmatrix}a&c\\b&d\end{pmatrix},\qquad
B=k[a,b,c,d]/(ad-bc-1).
\]

Define

\[
\mathcal H(a^rb^sc^td^u)
=(-1)^s\delta_{r,u}\delta_{s,t}
\frac{r!s!}{(r+s+1)!}.
\tag{3.1}
\]

### 3.1 Descent and the six derivations

For every monomial \(f=a^rb^sc^td^u\), the three terms of
\((ad-bc-1)f\) can contribute only under the same balance conditions
\(r=u,s=t\).  In that case (3.1) gives

\[
\mathcal H(adf)-\mathcal H(bcf)-\mathcal H(f)=0.
\tag{3.2}
\]

The left derivations are

\[
\begin{aligned}
L_H&=a\partial_a-b\partial_b+c\partial_c-d\partial_d,\\
L_E&=b\partial_a+d\partial_c,\\
L_F&=a\partial_b+c\partial_d,
\end{aligned}
\tag{3.3}
\]

and the right derivations are

\[
\begin{aligned}
R_H&=a\partial_a+b\partial_b-c\partial_c-d\partial_d,\\
R_E&=a\partial_c+b\partial_d,\\
R_F&=c\partial_a+d\partial_b.
\end{aligned}
\tag{3.4}
\]

All six preserve \(ad-bc-1\).  Direct substitution in (3.1) proves
\(\mathcal H(Df)=0\) for every \(D\) in (3.3)--(3.4).  The checker performs
this on the sparse derivative outputs; it never invokes conjugation,
\(SU(2)\), or compact measure.

### 3.2 The formula is forced

The two Cartan equations have weights

\[
r+t-s-u,\qquad r+s-t-u.
\]

Unless \(r=u\) and \(s=t\), at least one weight is nonzero and invariance
forces the monomial value to vanish.  Put

\[
h_{r,s}=\mathcal H(a^rb^sc^sd^r).
\]

Apply \(L_E\) to \(a^{r+1}b^sc^{s+1}d^r\).  Invariance gives

\[
(r+1)h_{r,s+1}+(s+1)h_{r+1,s}=0.
\tag{3.5}
\]

The determinant relation gives

\[
h_{r+1,s}-h_{r,s+1}=h_{r,s}.
\tag{3.6}
\]

Starting from \(h_{0,0}=1\), equations (3.5)--(3.6) uniquely yield

\[
h_{r+1,s}=\frac{r+1}{r+s+2}h_{r,s},\qquad
h_{r,s+1}=-\frac{s+1}{r+s+2}h_{r,s},
\]

which is exactly (3.1).  Thus (3.1) is the unique normalized left-and-right
infinitesimally invariant functional on \(k[SL_2]\).

### 3.3 Long's witness without compact integration

Let

\[
F=(1+c)(ad+b),\qquad G=-c.
\tag{3.7}
\]

In the expansion of \(F^m\), formula (3.1) retains only equal choices of
\(b\) and \(c\).  Therefore

\[
\mathcal H(F^m)
=\frac1{m+1}\sum_{j=0}^m(-1)^j\binom mj=0.
\tag{3.8}
\]

For \(GF^m\), the surviving indices differ by one, and the same calculation
gives

\[
\mathcal H(GF^m)
=\frac1{m+1}\sum_{j=1}^m(-1)^{j+1}\binom m{j-1}
=\frac{(-1)^{m-1}}{m+1}.
\tag{3.9}
\]

This is a proof in the coordinate Hopf algebra \(k[SL_2]\).  The older
\(SU(2)\) calculation remains an independent analytic reproduction, not a
dependency of (3.8)--(3.9).

## 4. The quotient map and its variance

There is a variance correction to the proposed ring map.  The third-column
morphism

\[
\pi:SO_3\longrightarrow X
\]

induces the pullback

\[
\pi^*:k[X]\longrightarrow k[SO_3],
\tag{4.1}
\]

not a canonical algebra map in the reverse direction.  In standard
orthogonal coordinates,

\[
\pi^*(U)=r_{13}+i r_{23},\qquad
\pi^*(V)=r_{13}-i r_{23},\qquad
\pi^*(T)=r_{33}.
\tag{4.2}
\]

The third-column norm relation gives \(UV+T^2=1\).  Right multiplication by
the stabilizer \(SO_2\) fixes the third column, and

\[
k[X]\simeq k[SO_3]^{SO_2}\hookrightarrow k[SO_3].
\tag{4.3}
\]

Equivalently,

\[
SO_3/SO_2\simeq X.
\]

A reverse algebra homomorphism compatible with (4.1) would be extra section
or algebra-retract data.  It is not produced by “taking the third column.”
The right-\(SO_2\) Reynolds operator does give a linear projection
\(k[SO_3]\to k[SO_3]^{SO_2}\), but a Reynolds projection is not generally an
algebra homomorphism.

The double-cover model makes (4.1) completely explicit without \(3\times3\)
matrix elimination.  For
\(g=\left(\begin{smallmatrix}a&c\\b&d\end{smallmatrix}\right)\), conjugation
of \(\operatorname{diag}(1,-1)\) gives

\[
U=-2ac,\qquad V=2bd,\qquad T=ad+bc.
\tag{4.4}
\]

Then

\[
UV+T^2=(ad-bc)^2,
\tag{4.5}
\]

and the three expressions are invariant under the right diagonal torus.
Thus (4.4) realizes

\[
SL_2/T\simeq PGL_2/T\simeq SO_3/SO_2\simeq X.
\]

The `SL(2)` checker verifies from (3.1), monomial by monomial, that

\[
\mathcal H_{SL_2}(\pi^*f)=\lambda(f).
\tag{4.6}
\]

Abstractly, (4.6) follows even faster: the left side is a normalized
\(\mathfrak {so}_3\)-invariant functional on \(A\), so the uniqueness theorem
of Section 1.3 identifies it with \(\lambda\).

## 5. The transfer theorem

Let \(G\) be a connected reductive group over an algebraically closed
characteristic-zero field, let \(H\subseteq G\) be reductive, and put
\(Y=G/H\).  The affineness of \(Y\) is the relevant case of
[Matsushima's criterion](https://encyclopediaofmath.org/wiki/Matsushima_criterion);
an algebraic proof is given by Arzhantsev in
[*Invariant Ideals and Matsushima's Criterion*](https://arxiv.org/abs/math/0506430).
Let \(\Lambda_Y:k[Y]\to k\) be the normalized \(G\)-invariant functional.

> **Equivariant-quotient transfer theorem.**  
> Suppose there is a surjective equivariant morphism
> \[
> \pi:Y\longrightarrow X=SO_3/SO_2,
> \]
> where equivariance is taken with respect to a homomorphism
> \(\alpha:G\to SO_3\).  Then
> \[
> \Lambda_Y\circ\pi^*=\lambda
> \]
> and \(\ker\Lambda_Y\) is not Mathieu.

Indeed, \(\pi^*\) is injective.  Its composite with \(\Lambda_Y\) is
normalized and invariant under the transitive image action on \(X\), hence
equals the unique functional \(\lambda\).  The elements
\(\pi^*P,\pi^*U\) then satisfy

\[
\Lambda_Y((\pi^*P)^m)=0,\qquad
\Lambda_Y(\pi^*U(\pi^*P)^m)
=\frac{4^m(m!)^2}{(2m+1)!}\ne0
\]

for every \(m\ge1\).

More generally, surjectivity and homogeneous-space language are stronger
than necessary.  Any injective algebra map

\[
\phi:A\hookrightarrow B
\]

with normalized functionals satisfying
\(\Lambda_B\phi=\lambda\) transfers the same non-Mathieu witness.  A retract
of \(\phi\) is a convenient sufficient certificate of injectivity, but is
not needed by the transfer argument.

## 6. Classification of equivariant quotients

The first two proposed sufficient conditions admit a sharp structural
classification.

Fix a homomorphism \(\alpha:G\to SO_3\).  An \(\alpha\)-equivariant morphism

\[
G/H\longrightarrow SO_3/SO_2
\]

is determined by the image \(x\) of \(eH\), and it exists exactly when

\[
\alpha(H)\subseteq\operatorname{Stab}_{SO_3}(x),
\tag{6.1}
\]

equivalently when \(\alpha(H)\) is contained in a conjugate of \(SO_2\).
It is surjective exactly when \(\alpha(G)\) acts transitively on the
quadric.

Since \(SO_3\simeq PGL_2\), every proper connected reductive subgroup has
dimension at most one and cannot act transitively on the two-dimensional
quadric.  Therefore, for connected reductive \(G\), surjectivity forces

\[
\alpha:G\twoheadrightarrow SO_3.
\tag{6.2}
\]

On Lie algebras this is a quotient
\(\mathfrak g\twoheadrightarrow\mathfrak {sl}_2\).  A reductive Lie algebra
has such a quotient exactly when its semisimple part has a simple ideal of
type \(A_1\).  Conversely, the adjoint action on such an ideal supplies a
surjection to \(PGL_2\).  Thus:

> **Quotient classification.**  
> For connected reductive \(G\), a surjective equivariant map
> \(G/H\to SO_3/SO_2\) exists exactly when \(G\) has an \(A_1\) quotient
> \(\alpha:G\twoheadrightarrow PGL_2\) for which \(\alpha(H)\) lies in a
> conjugate maximal torus.

This is readily automatable from the root datum of \(G\) and the image of
\(H\).  It also shows why “every noncommutative reductive group” is too
optimistic at this level: a root \(SL_2\) is a subgroup of \(G\), whereas
(6.2) requires an \(A_1\) quotient of \(G\).  Higher-rank simple groups have
many root \(SL_2\)-subgroups but no quotient onto \(PGL_2\).

## 7. What a root `SL(2)` orbit does and does not give

Let \(S\simeq SL_2\) be a root subgroup of \(G\).  Its orbit through \(gH\)
is

\[
S/(S\cap gHg^{-1}).
\]

If the stabilizer is a maximal torus and the orbit has the required closed
embedding, that orbit is a copy of the quadric \(S/T\).  But a subvariety
produces a restriction map in the direction

\[
k[G/H]\longrightarrow k[S/T].
\tag{7.1}
\]

The witness transfer needs an injection in the opposite direction.  Thus an
\(SL_2/T\) orbit, even a closed one, does not by itself provide the claimed
quadric subalgebra or a Haar-compatible retract.  One still needs an
extension map, quotient morphism, or explicit splitting of (7.1).

This variance obstruction is the precise gap between the abundant-root-
\(SL_2\) observation and a theorem for all non-torus reductive groups.

## 8. Revised computational search

For an affine orbit \(Y=Gv\), the symbolic search should target:

1. elements \(U,V,T\in k[Y]\) satisfying \(UV+T^2=1\);
2. the three infinitesimal transformation rules (1.7), under a chosen
   \(\mathfrak {sl}_2\)-quotient or action;
3. injectivity of \(A\to k[Y]\), for example by proving the image has
   transcendence degree two or that \(Y\to X\) is dominant;
4. normalized functional compatibility, proved from Lie invariance and the
   uniqueness theorem rather than inferred from finitely many moments.

Bounded Reynolds moments remain useful reconnaissance, but they are not a
certificate.  A discovered triple is promoted only after the quadric
relation, action identities, injectivity, and all-order functional
compatibility are proved.

For condition (6.2), no orbit-equation search is needed: inspect the root
datum for an \(A_1\) quotient and test the stabilizer containment (6.1).
For a root-subgroup strategy, the missing object to search for is not merely
an \(SL_2/T\) orbit but a compatible algebra section or quotient map.

## 9. Reproduction

Run

```bash
python3 scripts/verify_algebraic_quadric_haar.py
python3 scripts/verify_algebraic_sl2_haar.py
```

The first command checks quotient descent, all three infinitesimal
derivations, the forcing recurrences, Long's first twenty exact moment pairs,
and agreement with the existing spherical checker.  The second checks the
determinant ideal, all six left/right derivations, the factorial recurrences,
Long's first twenty `SL(2)` moment pairs, the explicit \(SL_2/T\) pullback,
and agreement with the older compact beta checker.

Long's external `SU(2)` source remains
[*Counterexamples to the xz-Conjecture and the Mathieu Conjecture for
SU(2)*](https://arxiv.org/abs/2607.19012).  The present note changes the
local proof mechanism, not the provenance of his witness.
