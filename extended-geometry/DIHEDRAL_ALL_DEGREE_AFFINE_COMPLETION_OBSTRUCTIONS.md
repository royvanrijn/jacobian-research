# All-degree dihedral affine-completion obstructions

## 1. Outcome

Let \(k\) have characteristic zero and contain the required dihedral
reflection constants.  For \(n\ge3\), define

\[
P_0=2,\qquad P_1=a,\qquad
P_n=aP_{n-1}-uP_{n-2},
\]

\[
J_n=\partial_aP_n,\qquad
C=a^2-4u,\qquad
\Delta_n(u,v)=v^2-4u^n.
\]

The uniform Dickson identity is

\[
\boxed{
\Delta_n(u,P_n)
=C\left(\frac{J_n}{n}\right)^2.
}                                                   \tag{1.1}
\]

This note generalizes the
[\(D_5\) two-mask audit](D5_TWO_MASK_BLOWDOWN_OBSTRUCTIONS.md)
to every degree \(n\ge3\).  It proves:

1. the complete odd and even branch-supported source-Jacobian ledgers;
2. a canonical two-mask blowdown of determinant \(\Delta_n\);
3. a polynomial Dickson-base obstruction in every degree;
4. all-degree rigidity of automorphic mask mixing, using the positive-genus
   generic curve \(y^2=4x^n+\delta\);
5. an all-degree contraction-divisor obstruction to the first
   nonautomorphic normalized-cusp chart; and
6. an all-degree no-go for every tangential chart affine in one normal
   coordinate; together with
7. the exact normal-degree resonance required of every nonlinear
   one-normal-coordinate continuation, and the resulting complete
   one-normal no-go in both parities.

The remaining stratum is uniform:

\[
\boxed{
\begin{gathered}
\text{two independent nonlinear tangent/normal}\\
\text{directions or a different determinant-\(\Delta_n\) blowdown.}
\end{gathered}}
\]

No Keller map and no complete affine-completion obstruction is proved.

## 2. Uniform source boundary

In splitting variables \(a=x+y,\ u=xy\),

\[
P_n(a,u)=x^n+y^n,\qquad
\frac{J_n}{n}=\frac{x^n-y^n}{x-y}.                 \tag{2.1}
\]

Therefore

\[
P_n^2-4u^n
=(x^n-y^n)^2
=(x-y)^2\left(\frac{x^n-y^n}{x-y}\right)^2,
\]

which is (1.1).

Over a fully split reflection field, write the distinct irreducible
ramification colors as

\[
\frac{J_n}{n}=R_1\cdots R_r,\qquad
r=\left\lfloor\frac n2\right\rfloor.
\]

The source boundary basis is

\[
([C],[R_1],\ldots,[R_r]).
\]

For odd \(n\), the target branch is irreducible and

\[
d=(1,2,\ldots,2),\qquad
j=(0,1,\ldots,1),                                  \tag{2.2}
\]

where \(d\) is the pullback of \(\Delta_n\) and \(j\) is the divisor of
\(J_n/n\).

## 3. Complete branch-supported valuation ledgers

### 3.1 Odd degree

Suppose \(n\) is odd and a target exceptional determinant has order
\(\ell\) along \(\Delta_n=0\).  If no other divisor meets the generic source
boundary points, constant Jacobian forces the source rechart orders \(s\)
to satisfy

\[
j+s=\ell d.
\]

Hence

\[
\boxed{
s=(\ell,2\ell-1,\ldots,2\ell-1),\qquad \ell\ge1.
}                                                   \tag{3.1}
\]

The primitive case is

\[
\ell=1,\qquad
\det D\alpha\sim C\frac{J_n}{n}.                   \tag{3.2}
\]

### 3.2 Even degree

Let \(n=2m\).  Put

\[
A_m=\frac{x^m-y^m}{x-y},\qquad
B_m=x^m+y^m=P_m(a,u).
\]

Both are symmetric polynomials in \(x,y\), and

\[
\frac{J_n}{n}=A_mB_m.                              \tag{3.3}
\]

The target branch has two components

\[
D_-:\ v-2u^m=0,\qquad
D_+:\ v+2u^m=0.
\]

Their pullbacks are

\[
\boxed{
P_n-2u^m=CA_m^2,\qquad
P_n+2u^m=B_m^2.
}                                                   \tag{3.4}
\]

In the composite boundary basis \(([C],[A_m],[B_m])\),

\[
d_-=(1,2,0),\qquad
d_+=(0,0,2),\qquad
j=(0,1,1).                                         \tag{3.5}
\]

If the target exceptional determinant has orders
\((\ell_-,\ell_+)\) along \((D_-,D_+)\), the unique source ledger is

\[
\boxed{
s=
(\ell_-,\,2\ell_--1,\,2\ell_+-1),
\qquad \ell_-,\ell_+\ge1.
}                                                   \tag{3.6}
\]

When \(A_m\) or \(B_m\) splits further, the displayed order applies to each
of its irreducible factors.  The primitive case
\(\ell_-=\ell_+=1\) again gives

\[
\det D\alpha\sim CA_mB_m=C\frac{J_n}{n}.           \tag{3.7}
\]

Thus (3.2) and (3.7) are one uniform primitive source ledger.  Even degree
also permits unequal higher orders on the two target components; an
all-degree obstruction must not silently collapse this two-parameter
possibility.

## 4. Canonical two-mask blowdown

For every \(n\ge3\), define

\[
\boxed{
\begin{aligned}
\beta_n(U,V,A,B)
=\bigl(
U,V,\ VA+2UB,\ 2U^{n-1}A+VB
\bigr).
\end{aligned}}                                      \tag{4.1}
\]

Its mask matrix is

\[
M_n=
\begin{pmatrix}
V&2U\\
2U^{n-1}&V
\end{pmatrix},
\]

so

\[
\det D\beta_n=\det M_n=V^2-4U^n=\Delta_n.          \tag{4.2}
\]

It is generically birational, with

\[
\boxed{
A=\frac{VS-2UT}{\Delta_n},\qquad
B=\frac{-2U^{n-1}S+VT}{\Delta_n}.
}                                                   \tag{4.3}
\]

The determinant is exactly the primitive total target branch factor.
Polynomiality still requires both adjugate numerators in (4.3) to be
divisible.

## 5. Polynomial Dickson-base obstruction

Let \(a_0,u_0,A,B\) be arbitrary polynomials on an affine source of
dimension four, and put

\[
F=(u_0,P_n(a_0,u_0),A,B),\qquad
\widetilde F=(a_0,u_0,A,B).
\]

The chain rule gives

\[
\boxed{
\det DF=-J_n(a_0,u_0)\det D\widetilde F.
}                                                   \tag{5.1}
\]

Since \(J_n\) is nonconstant for \(n\ge2\), \(F\) cannot have nonzero
constant Jacobian.

This excludes every factorization through (4.1) which makes the old
Dickson source coordinates \(a_0,u_0\) polynomially visible.  Nonlinear
dependence of those coordinates on the masks does not help.  A successful
stable field certificate must recover \(a_0\) only rationally.

## 6. Automorphic mask-mixing obstruction

Constant-linear coefficient comparison already forces the standard cusp
scalings.  If linear forms \(p,q\) satisfy

\[
q^2-4p^n=\lambda(V^2-4U^n),\qquad\lambda\ne0,
                                                               \tag{6.1}
\]

then mask \(n\)-th powers first force \(p\) to be mask-independent, mask
squares force the same for \(q\), and comparison of \(V^n,U^2,V^2,U^n\)
gives

\[
p=cU,\qquad q=dV,\qquad d^2=c^n.                   \tag{6.2}
\]

The first inverse numerator in (4.3) then has degree at most two, below
\(\deg\Delta_n=n\), so divisibility forces

\[
VS-2UT=0.
\]

For linear mask rows,

\[
S=2eU,\qquad T=eV,
\]

making the rechart singular.

The all-degree nonlinear version is geometric.

### Theorem 6.1

Let

\[
p,q\in k[U,V,Z_1,\ldots,Z_s]
\]

satisfy (6.1).  Then \(p,q\) are independent of every mask \(Z_i\).

### Proof

Hold all variables but \(Z_i\) in the coefficient field \(K\).  Equation
(6.1) gives a polynomial map

\[
\mathbb A^1_{Z_i}\longrightarrow
\mathcal C_{n,\delta}:\ y^2=4x^n+\delta,
\qquad
\delta=\lambda(V^2-4U^n)\in K^\times.              \tag{6.3}
\]

The smooth projective completion has genus

\[
g=\left\lfloor\frac{n-1}{2}\right\rfloor\ge1.
\]

The map extends to \(\mathbb P^1\), and Riemann--Hurwitz excludes a
nonconstant morphism from \(\mathbb P^1\) to a positive-genus curve.
Repeat for every \(Z_i\). \(\square\)

Therefore no polynomial automorphic rechart whose Keller condition reduces
to (6.1) can mix masks into the branch coordinates, in any degree
\(n\ge3\).

## 7. First nonautomorphic cusp obstruction

Put

\[
D=V^2-4U^n,\qquad
p=V^2+DS,\qquad
q=2V^n+DT.                                         \tag{7.1}
\]

Then

\[
E=\frac{q^2-4p^n}{D}
\]

is polynomial.  Completing the rechart while retaining \(U\) would require
a polynomial \(s(U,V,S,T)\) satisfying

\[
\det\frac{\partial(p,q,s)}{\partial(V,S,T)}=E.      \tag{7.2}
\]

But

\[
\nabla p\times\nabla q
=
\left(
D^2,\,
-2V(S+1)D,\,
-2V(T+nV^{n-2})D
\right),                                           \tag{7.3}
\]

up to simultaneous sign.  Hence every left side in (7.2) is divisible by
\(D\).  On the other hand,

\[
\boxed{
E\equiv
4V^n(T-nSV^{n-2})
\pmod D,
}                                                   \tag{7.4}
\]

which is nonzero.  No polynomial completion exists.

This failure is uniform and independent of the degree assigned to \(s\).

## 8. Affine-normal tangential obstruction

Let

\[
\begin{aligned}
h&=V+S,\\
p&=h^2+A(U,V,S)DT,\\
q&=2h^n+B(U,V,S)DT,\\
\rho&=(p,q,U,S).
\end{aligned}                                      \tag{8.1}
\]

The Jacobian \(\det D\rho\) has \(T\)-degree at most one.  If

\[
q^2-4p^n=\lambda D\det D\rho,                      \tag{8.2}
\]

the \(T^n\) coefficient on the left is

\[
-4A^nD^n.
\]

Since \(n\ge3\), the right side has lower \(T\)-degree, so \(A=0\).
The \(T^2\) coefficient is then

\[
B^2D^2,
\]

while the right side is independent of \(T\).  Hence \(B=0\), after which
\(\det D\rho=0\).

Therefore:

\[
\boxed{
\text{for every \(n\ge3\), no tangential cusp chart affine in one
normal coordinate is nondegenerate and log-crepant.}
}                                                   \tag{8.3}
\]

## 9. Arbitrary normal-degree gate

The remaining one-normal-coordinate search has an exact degree constraint
before any coefficients are solved.

### Proposition 9.1

Let \(K\) be a characteristic-zero field with a derivation
\(\partial\), extended to \(K[T]\) by \(\partial T=0\).  Suppose

\[
p,q\in K[T],\qquad
q^2-4p^n
=c\bigl((\partial p)q'-p'(\partial q)\bigr),
\qquad c\in K^\times.                              \tag{9.1}
\]

Then either \(p,q\) are both independent of \(T\), or, writing

\[
r=\deg_Tp,\qquad s=\deg_Tq,
\]

one has

\[
\boxed{2s=nr,\qquad q_s^2=4p_r^n,}                 \tag{9.2}
\]

where \(p_r,q_s\) are the leading coefficients.

### Proof

The right side of (9.1) has \(T\)-degree at most \(r+s-1\).  If just one
of \(r,s\) is zero, the left side has degree \(nr\) or \(2s\), strictly
larger than that bound.

Now let \(r,s>0\).  If \(2s>nr\), then

\[
2s>r+s-1.
\]

If \(nr>2s\), then

\[
nr>r+s-1
\]

because \(n\ge3\).  In either case the unique leading term on the left
cannot be matched by the right.  Thus \(2s=nr\), and cancellation of the
common top degree gives \(q_s^2=4p_r^n\). \(\square\)

For odd \(n\), (9.2) forces

\[
r=2d,\qquad s=nd,\qquad
p_r=c_0^2,\qquad q_s=\pm2c_0^n.                   \tag{9.3}
\]

The smallest unresolved normal-degree pair is therefore
\((r,s)=(2,n)\).  For even \(n=2m\),

\[
s=mr,\qquad q_s=\pm2p_r^m,                         \tag{9.4}
\]

so the smallest pair is \((1,m)\).  This parity split is essential:
a search beginning uniformly at bidegree \((2,n)\) would miss the first
even-degree stratum.

In the geometric chart, take

\[
K=k(U,V,S),\qquad \partial=\partial_V,\qquad
c=\lambda(V^2-4U^n).
\]

Then the bracket in (9.1) is precisely the two-variable Jacobian of
\((p,q)\) with respect to \((V,T)\), up to sign.  Proposition 9.1 therefore
applies to every one-normal-coordinate log-crepant rechart retaining
\((U,S)\), not only to the affine-normal form of Section 8.

### Theorem 9.2

Let \(n=2m\ge4\).  Under the hypotheses of Proposition 9.1, every solution
with nonconstant \(p,q\) has zero bracket.  In particular, there is no
nondegenerate log-crepant one-normal-coordinate rechart retaining
\((U,S)\).

### Proof

By (9.2), \(s=mr\) and the leading coefficients satisfy
\(q_s=\varepsilon2p_r^m\), where \(\varepsilon\in\{1,-1\}\).  Put

\[
b=q-\varepsilon2p^m.
\]

The second factor \(q+\varepsilon2p^m\) has degree \(mr\).  Since their
product is the left side of (9.1), whose degree is at most

\[
r+s-1=(m+1)r-1,
\]

one has

\[
\deg_Tb\le r-1.                                    \tag{9.5}
\]

The bracket with \(p^m\) vanishes, so

\[
(\partial p)q'-p'(\partial q)
=(\partial p)b'-p'(\partial b),                    \tag{9.6}
\]

of degree at most \(r+\deg_Tb-1\).  If \(b\ne0\), then

\[
q^2-4p^{2m}
=\varepsilon4bp^m+b^2
\]

has degree \(mr+\deg_Tb\), which is strictly greater than the bound in
(9.6) because \(m\ge2\).  This contradicts (9.1).  Hence \(b=0\),
so \(q=\varepsilon2p^m\) and both sides of (9.1), including the bracket,
vanish. \(\square\)

For even degree, the canonical one-normal-coordinate programme is
therefore closed at every polynomial degree.  For odd degree, Proposition
9.1 leaves precisely the resonant pairs

\[
(\deg_Tp,\deg_Tq)=(2d,nd),\qquad d\ge1.             \tag{9.7}
\]

### Theorem 9.3

Let \(n\ge3\) be odd.  Under the hypotheses of Proposition 9.1, every
solution with nonconstant \(p,q\) has zero bracket.  Hence there is no
nondegenerate log-crepant one-normal-coordinate rechart retaining
\((U,S)\).

### Proof

Write \(\deg_Tp=2d\) and \(\deg_Tq=nd\).  By (9.3), the leading
coefficient of \(p\) is a square.  In

\[
L=K(T)(w),\qquad w^2=p,
\]

choose the place at infinity for which

\[
w\sim c_0T^d.
\]

Put

\[
y=\frac{q}{2w^n}.
\]

Its leading residue is some \(\varepsilon\in\{1,-1\}\).  Direct
substitution into (9.1) gives

\[
\boxed{
w^{n-1}(y^2-1)
=c\bigl((\partial w)y'-w'(\partial y)\bigr).
}                                                   \tag{9.8}
\]

Suppose \(y\ne\varepsilon\), and let

\[
k=\operatorname{ord}_\infty(y-\varepsilon)>0.
\]

The left side of (9.8) has exact order

\[
k-(n-1)d.                                          \tag{9.9}
\]

The derivation extends by

\[
\partial w=\frac{\partial p}{2w},
\qquad
w'=\frac{p'}{2w}.
\]

Consequently

\[
\begin{aligned}
\operatorname{ord}_\infty(\partial w)&\ge-d,\\
\operatorname{ord}_\infty(w')&\ge-d+1,\\
\operatorname{ord}_\infty(y')&\ge k+1,\\
\operatorname{ord}_\infty(\partial y)&\ge k.
\end{aligned}
\]

The right side of (9.8) therefore has order at least

\[
k-d+1.                                             \tag{9.10}
\]

But

\[
k-(n-1)d<k-d+1
\]

for \(n\ge3\) and \(d\ge1\), a contradiction.  Thus
\(y=\varepsilon\), so \(q^2=4p^n\).  Since \(n\) is odd and \(K[T]\) is
a unique factorization domain, there is \(z\in K[T]\) with

\[
p=z^2,\qquad q=\varepsilon2z^n.
\]

The bracket vanishes. \(\square\)

Combining Theorems 9.2 and 9.3:

\[
\boxed{
\text{the canonical one-normal-coordinate log-crepant programme is
degenerate for every \(n\ge3\).}
}                                                   \tag{9.11}
\]

## 10. Uniform obstruction table

| construction class | all-\(n\) verdict |
|---|---|
| one retained polynomial Dickson base | excluded by (5.1) |
| constant-linear base/mask rechart | excluded by coefficient rigidity and adjugate dependence |
| polynomial automorphic mask mixing | excluded by Theorem 6.1 |
| symmetric normalized-cusp corrections \(DS,DT\) | excluded by (7.3)--(7.4) |
| tangential chart affine in one normal coordinate | excluded by (8.3) |
| nonlinear dependence on one normal coordinate, even \(n\) | excluded by Theorem 9.2 |
| nonlinear dependence on one normal coordinate, odd \(n\) | excluded by Theorem 9.3 |
| two independent nonlinear tangent/normal directions | open |
| different determinant-\(\Delta_n\) blowdown | open |

For odd \(n\), the valuation ledger (3.1) covers the whole
single-branch-supported class.  For even \(n\), the independent orders in
(3.6) must also be considered if the two branch components are treated
asymmetrically.

## 11. Search compiler and false-counterexample triage

No one-normal coefficient system should now be sent to Gröbner,
interpolation, or numerical solving.  For a claimed candidate over
\(K[T]\), compute

\[
r=\deg_Tp,\qquad s=\deg_Tq.
\]

The exact routing chain is:

| condition | certificate | verdict |
|---|---|---|
| exactly one of \(r,s\) is zero | degree of (9.1) | reject |
| \(r,s>0\) and \(2s\ne nr\) | Proposition 9.1 | reject |
| \(2s=nr\), \(n\) even | cusp factorization, Theorem 9.2 | zero bracket |
| \(2s=nr\), \(n\) odd | infinity valuation, Theorem 9.3 | zero bracket |

Thus resonance is a necessary intermediate condition, never a surviving
search branch.

Two simple false positives are useful regressions.  Over
\(K=k(V)\), with \(\partial=\partial_V\), the even resonant-looking pair

\[
p=T+V,\qquad q=2p^{n/2}+V
\]

has bracket \(-1\), but

\[
q^2-4p^n=4Vp^{n/2}+V^2
\]

still depends on \(T\).  For odd \(n\), the pair

\[
p=T^2+V,\qquad q=2T^n
\]

has the resonant degrees \((2,n)\), but the cusp difference has degree
\(2n-2\), while its bracket has degree \(n-1\).  These examples prevent a
search from treating the leading resonance as sufficient evidence for a
counterexample.

A candidate advances beyond this compiler only if it supplies one of:

1. two genuinely active normal variables, so no two coordinates can be
   retained to reduce the equation to \(K[T]\); or
2. a different birational blowdown whose exceptional determinant is
   \(\Delta_n\).

Only then should the expensive adjugate-divisibility, stable-field,
factoriality, and complete-fibre checks run.

## 12. Next theorem

The most useful continuation is no longer a degree-by-degree coefficient
search.  It is the following uniform statement.

> **Two-normal dihedral obstruction.**  Every polynomial log-crepant
> rechart with two coupled normal coordinates either has degenerate
> Jacobian, retains an adjugate pole, or changes the stable Dickson
> extension.

Proving it would close the canonical two-mask programme for all \(n\ge3\).
Failure would produce a sharply structured candidate rather than an
uninterpreted coefficient component.  A different determinant-\(\Delta_n\)
blowdown remains logically separate.

## 13. Reproduction

Run

```bash
.venv/bin/python scripts/verify_dihedral_all_degree_affine_completion_obstructions.py
```

The checker replays the Dickson branch identity, the even-component
factorizations, the canonical blowdown, the positive-genus condition, the
first nonautomorphic cusp remainder, and the affine-normal coefficient
locks for \(3\le n\le12\).  It also checks the normal-degree inequality
behind Proposition 9.1, the even factorization gap, and the odd
valuation-at-infinity gap in bounded ranges.  It also checks the two
false-positive resonant pairs and the search-routing compiler.  The bounded
replay is a regression certificate; the proofs above are uniform.
