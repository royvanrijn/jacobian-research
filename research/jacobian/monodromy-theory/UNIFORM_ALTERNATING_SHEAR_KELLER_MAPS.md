# Uniform polynomial Keller maps with alternating generic monodromy

Work over `K=Q(i)` and its algebraic closure, in characteristic zero. For every odd integer `m>=3`, the formulas below give a polynomial Keller map of affine three-space with full geometric monodromy `A_(m(m+1))`.

The member `m=3` is the alternating-shear construction of Kistner–Shaska, [arXiv:2608.02863v1](https://arxiv.org/pdf/2608.02863v1), Section 5.1, PDF pp.15–17. Its twelve displayed coefficient terms, lift conditions and actual three-variable Jacobian have been checked exactly here. The uniform formulas and proofs below extend the same odd-shear mechanism. This is not a claim that generic alternating Keller realizations were previously unknown: the cited paper already supplies them through several constructions. Priority of the uniform formulation has not been assessed.

## 1. Explicit uniform construction

Let

\[
p=\frac{m-1}{2},\qquad
c=\frac{m(m-i)}{m^2-1},\qquad d=(1-c)^{-m}.
\]

For independent variables `Lambda,w`, put

\[
C_1=1+c\Lambda^p w^m,\qquad
C_2=w+d\Lambda^p C_1^m,\qquad A=C_1C_2.
\]

Define the polynomial

\[
\begin{split}
B={}&C_1^2
-2c\sum_{j=0}^{m}
 \frac{\binom mj(-d)^j}{(p+1)j-p}
 \Lambda^{p(j+1)}C_1^{1+mj}C_2^{m-j}\\
&+c^2\sum_{j=0}^{2m}
 \frac{\binom{2m}j(-d)^j}{(p+1)j-m}
 \Lambda^{p(j+2)}C_1^{mj}C_2^{2m-j}.
\end{split}
\tag{1}
\]

None of the denominators vanishes: `p+1` divides neither `p` nor `m=2p+1`. All exponents in this formula are nonnegative.

On source coordinates `x,y,z`, set

\[
\Lambda=1+xy,\qquad w=xz-1,
\]

and define

\[
\boxed{F_m(x,y,z)=
\left(\frac{A(1+xy,xz-1)}{x},\
      \frac{B(1+xy,xz-1)}{x},\
      x(1+xy)\right).}
\tag{2}
\]

**Theorem.** Formula (2) is a polynomial map `A^3_K -> A^3_K` with determinant `-1`. Its generic degree is

\[
N=m(m+1),
\]

and its geometric generic monodromy, in its natural degree-`N` action, is `A_N`. The arithmetic group over `K` is `A_N` when `m` is a perfect square and `S_N` otherwise; in the latter case the normal closure has constant field `K(sqrt(m))`.

For `m=3`, the coefficients are

\[
c=\frac{9-3i}{8},\qquad d=\frac{1664+1152i}{125},
\]

exactly the indicated primary model. Its component degrees are `(33,159,3)`. The subsequent uniform geometric degrees are `30,56,90,132,...`.

No rational-coefficient model or descent to `Q` is asserted by this theorem. It is a construction in dimension three, not a polynomial plane Keller realization.

## 2. Polynomiality and the lift point

Set `a=1-c`. At `(Lambda,w)=(1,-1)` one has `C_1=a` and

\[
C_2=-1+da^m=0.
\]

Thus `A(1,-1)=0`. In (1), only the terms with zero `C_2` exponent survive at that point. Writing `D_1=(m^2+1)/2`, they give

\[
B(1,-1)=a^2+\frac{2ca}{D_1}+\frac{c^2}{m^2}.
\tag{3}
\]

The quadratic on the right has roots `c=m(m±i)/(m^2-1)` when expressed in `c`, so our choice makes it zero. More explicitly its coefficients are

\[
\frac{(m^2-1)^2}{m^2(m^2+1)}c^2
-\frac{2(m^2-1)}{m^2+1}c+1.
\]

Substitution proves the identity for every `m`. Both numerators in (2) therefore vanish on `x=0`, where `(Lambda,w)=(1,-1)`, and are divisible by `x` in `K[x,y,z]`. This proves polynomiality everywhere.

One can avoid even a displayed quotient by using polynomial divided differences. For either `G=A` or `B`, with `G(1,-1)=0`, write

\[
G(\Lambda,w)=(\Lambda-1)D_\Lambda G(\Lambda,w)
 +(w+1)D_wG(1,w).
\]

The corresponding coordinate in (2) is the polynomial
`y D_Lambda G(1+xy,xz-1)+z D_wG(1,xz-1)`.

## 3. A birational shear calculation gives the Jacobian

Adjoin `alpha,beta` with

\[
\Lambda=\alpha^2,\qquad w=\beta/\alpha.
\]

The involution changes both signs and has invariant field `K(Lambda,w)`. The two odd polynomial shears

\[
\alpha_1=\alpha+c\beta^m,\qquad
\beta_1=\beta+d\alpha_1^m
\]

are an automorphism of determinant one and commute with that involution. Since

\[
\alpha_1=\alpha C_1,\quad\beta_1=\alpha C_2,
\quad X=\alpha_1^2=\Lambda C_1^2,
\quad Y=\alpha_1\beta_1=\Lambda A,
\]

one has `K(X,Y)=K(Lambda,w)`. Moreover

\[
\operatorname{Jac}_{\Lambda,w}(X,Y)=X.
\tag{4}
\]

Indeed the Jacobians of `(alpha,beta)->(Lambda,w)` and `(alpha_1,beta_1)->(X,Y)` are `2` and `2X`, respectively.

For the rational function

\[
\begin{split}
P(X,Y)={}&X
-2c\sum_{j=0}^{m}
 \frac{\binom mj(-d)^jY^{m-j}X^{(p+1)j-p}}{(p+1)j-p}\\
&+c^2\sum_{j=0}^{2m}
 \frac{\binom{2m}j(-d)^jY^{2m-j}X^{(p+1)j-m}}{(p+1)j-m},
\end{split}
\tag{5}
\]

termwise differentiation gives

\[
P_X=
\left(1-cX^{-(p+1)}(Y-dX^{p+1})^m\right)^2
=\frac{\Phi(X,Y)^2}{X^{m+1}},
\qquad
\Phi=X^{p+1}-c(Y-dX^{p+1})^m.
\tag{6}
\]

There is no logarithmic term because the two excluded denominators are precisely the possible residues. Substituting `X=Lambda C_1^2,Y=Lambda C_1C_2` in (5) gives `P=Lambda B`. For example a term from its first sum becomes

\[
\Lambda^{p(j+1)+1}C_1^{1+mj}C_2^{m-j},
\]

and a term from its second becomes
`Lambda^(p(j+2)+1) C_1^(mj) C_2^(2m-j)`, proving (1) directly.

The inverse shear gives

\[
\frac{\alpha}{\alpha_1}
=1-cX^{-(p+1)}(Y-dX^{p+1})^m.
\]

Consequently `P_X=Lambda/X`. Combining with (4) proves

\[
\operatorname{Jac}_{\Lambda,w}(P,Y)=\Lambda.
\tag{7}
\]

Since `P=Lambda B` and `Y=Lambda A`, this is exactly

\[
\Lambda(B_\Lambda A_w-B_wA_\Lambda)+BA_w-AB_w=1.
\tag{8}
\]

For the full three-dimensional map, put `S=F_3=x Lambda`. The source change to `(Lambda,w,x)` has Jacobian `x^2`, the target change to `(Y=F_1S,P=F_2S,S)` has Jacobian `S^2`, and (7) gives the Jacobian of `(Y,P,S)` with respect to `(Lambda,w,x)` as `-Lambda^2`. Hence `det DF_m=-1` on a dense open set. All components are polynomial, so the identity holds everywhere.

## 4. The generic extension has degree N

The element `S=x Lambda` is independent of `Lambda,w`, and reconstruction gives

\[
x=S/\Lambda,\qquad y=(\Lambda-1)/x,\qquad z=(w+1)/x.
\]

Thus the source and target fields are respectively

\[
K(X,Y,S),\qquad K(P,Y,S).
\tag{9}
\]

They are obtained by adjoining the same independent transcendental element to the rational cover (5), preserving its normal closure group.

For generic nonzero `Y`, the poles of `P` have orders `m` at zero and `m^2` at infinity, with leading coefficients

\[
-\frac{c^2}{m}Y^{2m}X^{-m},\qquad
\frac{c^2d^{2m}}{m^2}X^{m^2},
\]

respectively. Both are nonzero. Hence the cover degree and the degree of (9) are `m+m^2=N`.

## 5. Individual three-cycle inertia and distinct branch values

For generic `Y`, the `N/2=m(p+1)` zeros of `Phi` are nonzero and simple. To see simplicity, put `T_0=X^(p+1)`. The equation is

\[
T_0=c(Y-dT_0)^m.
\]

A repeated root requires
`1+cmd(Y-dT_0)^(m-1)=0`, which allows only finitely many `Y`; `Y=0` is also excluded. Equation (6) then makes each finite critical point have local degree three.

We must check that these critical points have different branch values for generic `Y`. The critical curve has the rational parametrization

\[
X=c^2t^m,\qquad
T_0'=Y-dX^{p+1}=ct^{p+1},\qquad
Y_c(t)=ct^{p+1}(1+c^mdt^h),\qquad h=2p(p+1).
\tag{10}
\]

Its rational inverse is `t=(T_0')^2/X`, since `2(p+1)-m=1`. The binomial critical curve is irreducible because `gcd(m,p+1)=1`.

Substitution into (5) gives a polynomial

\[
P_c(t)=c^2t^m\Xi(c^mdt^h),
\]

where

\[
\begin{split}
\Xi(Z)={}&1-2\sum_{j=0}^{m}
 \frac{\binom mj(-1)^jZ^j(1+Z)^{m-j}}{(p+1)j-p}\\
&+\sum_{j=0}^{2m}
 \frac{\binom{2m}j(-1)^jZ^j(1+Z)^{2m-j}}{(p+1)j-m}.
\end{split}
\]

The two required endpoint values are

\[
\Xi(0)=1+\frac2p-\frac1m
=\frac{2(p+1)^2}{pm}\ne0,
\qquad
\Xi(-1)=1-\frac2{D_1}+\frac1{m^2}
=\frac{(m^2-1)^2}{m^2(m^2+1)}\ne0.
\tag{11}
\]

Therefore `P_c,Y_c` have exactly one common zero, namely `t=0`, and their orders there are `m,p+1`, which are coprime.

This proves that `t -> (P_c,Y_c)` is birational onto its image. Indeed the map is finite onto that affine image, since `t` is integral over the ring generated by the nonconstant polynomial `Y_c`. The induced finite map to the normalization has only the point `t=0` above the image of the origin. Its ramification index divides both displayed orders, so it is one; with algebraically closed residue fields and a unique point over that place, the extension degree is one. Infinity is not another point above the origin, because `Y_c` has a pole there. Thus generic critical values are distinct.

The finite inertia permutations are consequently individual three-cycles, one for each of `N/2` critical points. The infinity permutation has cycle type `(m^2,m)`. Riemann–Hurwitz checks the whole branch count:

\[
\frac N2(3-1)+(m^2-1)+(m-1)=2N-2.
\]

## 6. Full geometric monodromy and the arithmetic constants

A connected cover of the punctured projective line has transitive monodromy generated by its finite branch permutations; the infinity permutation is the inverse of their product. A transitive group generated by individual three-cycles is the full alternating group. Here is an elementary proof of that group fact.

The three-element supports form a connected hypergraph, since their connected components are exactly the group orbits. Start with one support, on which the group contains `A_3`, and add supports meeting the preceding union. A three-cycle with two old symbols and one new enlarges `A_r` to `A_(r+1)` by conjugation; the case `r=3` is checked by its three unordered pairs, and for `r>=4` the old group is two-transitive. For a support with one old and two new symbols, conjugate `(a u v)` to `(b u v)` using an old element. Their quotient is a three-cycle on the old symbols `a,b` and the new `u`. Apply the preceding case to add `u`, then use the original cycle to add `v`. Induction gives the full alternating group on the connected support union.

This gives `A_N` after algebraically closing the parameter field `K(Y)`. To establish the ambient geometric group over `Kbar(P,Y)` as well, it remains to rule out an odd quadratic twist in `Y`. The following exact discriminant calculation does that and also determines arithmetic monodromy.

Use an independent target variable `T` and put

\[
f(X)=X^m(P(X,Y)-T),\quad
 a_N=\frac{c^2d^{2m}}{m^2},\quad
 a_0=-\frac{c^2}{m}Y^{2m}.
\]

It is a degree-`N` polynomial and satisfies

\[
Xf'(X)-mf(X)=\Phi(X,Y)^2.
\tag{12}
\]

At its roots `alpha`, one has `alpha f'(alpha)=Phi(alpha)^2`. Using `deg Phi=N/2`, the discriminant-resultant identity gives

\[
\operatorname{Disc}_X(f)
=(-1)^{N(N-1)/2}\frac{\operatorname{Res}_X(f,\Phi)^2}{a_Na_0}
=D_m\left(
\frac{m\operatorname{Res}_X(f,\Phi)}{c^2d^mY^m}
\right)^2,
\tag{13}
\]

where

\[
D_m=(-1)^{(m-1)/2}m.
\]

In particular the square class is constant; there is no residual `Y` twist. Over an algebraically closed coefficient field, the ambient group is contained in `A_N` by (13) and contains the already computed `A_N`, so it equals `A_N`. Equation (9) proves the same assertion for the full polynomial Keller map.

Over `K=Q(i)`, the constant `D_m` is a square exactly when the positive integer `m` is a square. Since arithmetic monodromy contains the geometric `A_N`, it is `A_N` in that case and `S_N` otherwise. In the latter case (13) places `sqrt(D_m)` in the normal closure; the geometric group has index two, so the full constant extension is exactly `K(sqrt(D_m))=K(sqrt(m))`. For `m=3` this is `Q(i,sqrt(3))`, equivalently `Q(i,sqrt(-3))`. This statement concerns the map over `Q(i)` and does not claim that the lift has coefficients in `Q`.

This completes the theorem.

## 7. Exact degree-twelve check and scope

The retained [primary PDF](sources/kistner-shaska-graded-keller-2608.02863v1.pdf) has SHA-256 `3cd7e239d8b01c79239789a97c26e466208d9cb216d9773574562c468c915a2e`. Section 5.1, equations (21)–(24), gives the `m=3` formulas. Proposition 5.3 and Theorems 5.4–5.5 assert polynomiality, the Jacobian and monodromy; the argument above verifies the needed assertions directly and gives an elementary full-group proof.

The [model](alternating_shear_keller_model.py), [literal input](uniform_alternating_shear_keller_inputs.json) and [checker](verify_uniform_alternating_shear_keller.py) reproduce all twelve coefficient terms, expand `A,B` over `QQ(i)` and check (8) and both lift conditions exactly. The actual three-variable polynomials have `77,1484,2` monomials and degrees `33,159,3`; their full expanded Jacobian determinant is also differentiated directly and checked as `-1`.

The [default receipt](uniform_alternating_shear_keller_verification.json), [full Jacobian receipt](uniform_alternating_shear_keller_full_jacobian_verification.json) and [source ledger](uniform_alternating_shear_keller_source_ledger.json) separate those computations from the written all-odd-`m` proof. Bounded Laurent/critical-curve checks for `m=3,5,...,15` are supplementary; the lift-root identity is checked symbolically in `m`. Corrupted primary coefficients and even shear degrees are rejected.

```bash
.venv/bin/python research/jacobian/monodromy-theory/verify_uniform_alternating_shear_keller.py \
  --write-receipt research/jacobian/monodromy-theory/uniform_alternating_shear_keller_verification.json
.venv/bin/python research/jacobian/monodromy-theory/verify_uniform_alternating_shear_keller.py \
  --full-jacobian \
  --write-receipt research/jacobian/monodromy-theory/uniform_alternating_shear_keller_full_jacobian_verification.json
```

This architecture lies outside the quadratic-tilt compiler audited earlier. Its rational-cover parameters are the invariant products `P=F_2F_3,Y=F_1F_3`, not retained independent affine output coordinates. Also the Laurent pole coefficient is proportional to `Y^(2m)` and vanishes on `Y=0`, unlike the nonzero constant pole term in the fixed Laurent extension obstruction. The affine-base and unit hypotheses of that obstruction therefore do not apply to this different embedded extension.

The construction supplies generic geometric alternating monodromy for polynomial maps of affine three-space. It does not produce a plane Keller map or a universal obstruction to one, does not prove rational descent of these lifts, and does not imply that group-theoretic necessary conditions are sufficient. No formal verification, external review or novelty claim is made.