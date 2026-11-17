# Cusp-profile suspension theorem for homogeneous GVC(3)

## 1. Statement

Work over a characteristic-zero field.  Put

\[
 \rho=t^2+xy,\qquad A=\rho+x^2,
 \qquad
 C=\frac{\rho^3-t^2A^2}{x}.
 \tag{1.1}
\]

The quotient in (1.1) is polynomial.  Indeed,

\[
 \boxed{\rho^3-t^2A^2=xC},
 \qquad
 C=y\rho^2-2xt^2\rho-x^3t^2.
 \tag{1.2}
\]

Let \(r\ge1\), \(e,h\ge0\), and

\[
 S(z)=\sum_{j=0}^e s_jz^j\ne0.
\]

Its level-\(e\) homogeneous lift is

\[
 S^{\mathrm{hom}}
 =\sum_{j=0}^e
 s_j(t^2A^2)^j\rho^{3(e-j)}.
 \tag{1.3}
\]

Normally \(e=\deg S\).  Allowing trailing zero coefficients merely moves a
power of \(\rho^3\) from the profile into the radial suspension.

Define

\[
 \boxed{
 P_{r,S,h}=\rho^hA^rC^{2r}S^{\mathrm{hom}},
 \qquad
 N=6r+3e+h,
 \qquad
 \Lambda=\Delta^N,
 }
 \tag{1.4}
\]

where

\[
 \Delta=4\partial_x\partial_y+\partial_t^2.
\]

For a polynomial \(f(z)=\sum a_jz^j\), write

\[
 \mathcal I(f)=\sum_j\frac{a_j}{2j+1}.
 \tag{1.5}
\]

Over the reals this is \(\int_0^1f(v^2)\,dv\).  Put

\[
 c_m(S)
 =\mathcal I\left((1-z)^{2rm}S(z)^m\right)
 =\int_0^1(1-v^2)^{2rm}S(v^2)^m\,dv.
 \tag{1.6}
\]

> **Theorem 1.1 — cusp-profile suspension.**  The polynomial
> \(P_{r,S,h}\) is homogeneous of degree \(2N\).  If
> \(c_m(S)\ne0\) for every \(m\ge1\), then
>
> \[
>  \boxed{\Delta^{Nm}(P_{r,S,h}^m)=0}
>  \qquad(m\ge1),
>  \tag{1.7}
> \]
>
> and, for \(1\le\ell\le rm\),
>
> \[
> \boxed{
> \begin{aligned}
> \Delta^{Nm+\ell}
> \left(x^{2\ell}P_{r,S,h}^m\right)
> ={}&2^{Nm+\ell}(Nm+\ell)!
> (2Nm+2\ell+1)!!\\
> &\times
> \binom{rm-1}{\ell-1}c_m(S)\ne0.
> \end{aligned}}
> \tag{1.8}
> \]

Taking \(\ell=1\) shows that

\[
 \Lambda^m(x^2P_{r,S,h}^m)\ne0
 \qquad(m\ge1),
 \tag{1.9}
\]

because one further application of \(\Delta\) gives the nonzero scalar in
(1.8).  Thus \(Q=x^2\) is a strong GVC detector at every positive power.

A simple sufficient condition for (1.6) is

\[
 S(z)>0\qquad(0\le z\le1).
 \tag{1.10}
\]

For fixed \(r,e,h\), strict positivity is a nonempty open subset of the
\((e+1)\)-dimensional real coefficient space.  Hence (1.4) gives
full-dimensional real profile families, not isolated coefficient points.

## 2. Reynolds--apolar transfer

The differential passage in this theorem is an instance of one general
identity.

Let \(q(z)=z^{\mathsf T}Az\) be a nondegenerate quadratic form in \(n\)
variables, and define its dual Laplacian by

\[
 D_q=\partial^{\mathsf T}A^{-1}\partial.
 \tag{2.1}
\]

Let \(\mathcal R_q\) be the normalized orthogonal Reynolds functional on

\[
 k[z_1,\ldots,z_n]/(q-1).
\]

> **Lemma 2.1 — Reynolds--apolar transfer.**  If \(F\) is homogeneous of
> degree \(2k\), then
>
> \[
> \boxed{
> D_q^kF
> =4^kk!\left(\frac n2\right)_k
> \mathcal R_q(F|_{q=1}).}
> \tag{2.2}
> \]

Both sides of (2.2) are orthogonally invariant linear functionals on
\(\operatorname{Sym}^{2k}\), whose invariant-functional space is
one-dimensional.  It remains only to evaluate on \(F=q^k\).  Directly,

\[
 D_q(q^a)=4a\left(a+\frac n2-1\right)q^{a-1},
\]

and iteration gives \(4^kk!(n/2)_k\), while
\(\mathcal R_q(q^k|_{q=1})=1\).  This proves the lemma over every
characteristic-zero field after scalar extension and descent.

For \(q=\rho=t^2+xy\), formula (2.1) gives

\[
 D_\rho=4\partial_x\partial_y+\partial_t^2=\Delta.
\]

Since \(n=3\), (2.2) becomes

\[
 \boxed{
 \Delta^kF
 =2^kk!(2k+1)!!\,
 \mathcal R_\rho(F|_{\rho=1}).}
 \tag{2.3}
\]

Thus spherical moments, Gaussian moments, top Laplacian contractions, and
the algebraic quadric Reynolds functional are the same homogeneous
calculation with different normalizations.  The rank-three functional is
constructed directly from the \(\mathfrak {so}_3\)-derivations in the
[algebraic Haar note](ALGEBRAIC_HAAR_QUADRIC_AND_SL2.md).

## 3. The phase calculation

On \(\rho=1\), put

\[
 u=x^2,\qquad B=1+u.
\]

Equations (1.2)--(1.4) become

\[
 P_{r,S,h}\big|_{\rho=1}
 =u^{-r}B^r
   (1-t^2B^2)^{2r}S(t^2B^2).
 \tag{3.1}
\]

In the localization of the quadric ring obtained by substituting
\(y=(1-t^2)/x\), the functional \(\mathcal R_\rho\) first extracts the
constant phase in \(x\), kills odd powers of \(t\), and sends
\(t^{2j}\) to \(1/(2j+1)\).  This is the algebraic form of uniform phase
and height averaging.

After taking the \(m\)-th power of (3.1), phase extraction gives

\[
 \begin{aligned}
 \mathcal R_\rho(P_{r,S,h}^m)
 &=[u^{rm}]H_m(1+u),\\
 \mathcal R_\rho(x^{2\ell}P_{r,S,h}^m)
 &=[u^{rm-\ell}]H_m(1+u),
 \end{aligned}
 \tag{3.2}
\]

where

\[
 H_m(B)
 =B^{rm}\int_0^1
 \left((1-v^2B^2)^{2r}S(v^2B^2)\right)^m\,dv.
 \tag{3.3}
\]

After \(w=vB\),

\[
 H_m(B)=B^{rm-1}J_m(B),
 \tag{3.4}
\]

with

\[
 J_m'(B)=\left((1-B^2)^{2r}S(B^2)\right)^m.
 \tag{3.5}
\]

The derivative in (3.5) has a zero of order at least \(2rm\) at \(B=1\).
Consequently

\[
 J_m(1+u)=c_m(S)+O(u^{2rm+1}).
 \tag{3.6}
\]

Through degree \(rm\), equations (3.4)--(3.6) therefore give

\[
 H_m(1+u)=c_m(S)(1+u)^{rm-1}.
 \tag{3.7}
\]

It follows that

\[
 [u^{rm}]H_m(1+u)=0,
 \tag{3.8}
\]

while, for \(1\le\ell\le rm\),

\[
 \boxed{
 [u^{rm-\ell}]H_m(1+u)
 =\binom{rm-1}{\ell-1}c_m(S).}
 \tag{3.9}
\]

Applying (2.3) to (3.8)--(3.9) proves (1.7)--(1.8).  This is precisely the
full lower-jet ladder from the
[one-profile Hopf classification](HOPF_LIFT_CLASSIFICATION.md), now with
the endpoint contact and phase-square winding homogenized internally in
three variables.

## 4. Maximal shifted-power failure

The whole shifted-power tail also survives.  Put

\[
 \delta(X,Y,T)=4XY+T^2,
\]

the symbol of \(\Delta\).  Fischer adjointness gives

\[
 S(\partial)(x^{2d}F)|_0
 =(\partial_X^{2d}S)(\partial)F|_0,
 \tag{4.1}
\]

and

\[
 \partial_X^{2d}\delta^K
 =4^{2d}\frac{K!}{(K-2d)!}
 Y^{2d}\delta^{K-2d}.
 \tag{4.2}
\]

Take the ladder formula (1.8) at power \(m+d\) and level \(\ell=d\).
Equations (4.1)--(4.2) yield

\[
\boxed{
\begin{aligned}
&\Delta^{(N-1)d}\partial_y^{2d}
 \left(\Lambda^m(P_{r,S,h}^{m+d})\right)\\
&\quad=
 2^{Nm+(N-3)d}
 \bigl(Nm+(N-1)d\bigr)!
 \bigl(2Nm+2(N+1)d+1\bigr)!!\\
&\qquad\times
 \binom{r(m+d)-1}{d-1}c_{m+d}(S)\ne0.
\end{aligned}}
\tag{4.3}
\]

Therefore:

> **Corollary 4.1 — maximal power-tail failure.**  For every \(m,d\ge1\),
>
> \[
>  \boxed{\Lambda^m(P_{r,S,h}^{m+d})\ne0.}
> \tag{4.4}
> \]

Thus every member of the family violates every shifted-power formulation,
not merely the conclusion for one external multiplier.

## 5. Exact polyharmonic depth

The \(\ell=1\) case of (1.8), together with (4.2) for two \(X\)-derivatives,
gives

\[
\boxed{
 \partial_y^2\Delta^{Nm-1}(P_{r,S,h}^m)
 =2^{Nm-3}(Nm-1)!(2Nm+3)!!c_m(S)\ne0.}
\tag{5.1}
\]

Combining (5.1) with (1.7) proves:

> **Corollary 5.1 — exact depth.**  The power \(P_{r,S,h}^m\) has exact
> polyharmonic index \(Nm\):
>
> \[
>  \Delta^{Nm}(P_{r,S,h}^m)=0,
>  \qquad
>  \Delta^{Nm-1}(P_{r,S,h}^m)\ne0.
> \tag{5.2}
> \]

The pure cancellation always occurs at the last possible trace.

## 6. The three old families are specializations

The construction has three independent operations:

\[
 \boxed{
 \text{winding }r
 \quad+\quad
 \text{profile }S
 \quad+\quad
 \text{radial suspension }h.}
 \tag{6.1}
\]

They recover the previously separate families as follows.

### 6.1 The original witness

For

\[
 r=1,\qquad S=1,\qquad h=0,
\]

one obtains

\[
 P_{1,1,0}=AC^2,qquad N=6.
\]

This is the
[23-term homogeneous GVC(3) witness](THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md).

### 6.2 The endpoint-power family

For \(s\ge2\), take

\[
 r=1,\qquad S(z)=(1-z)^{s-2},\qquad h=0.
\]

Then \(e=s-2\) and

\[
 \begin{aligned}
 S^{\mathrm{hom}}
 &=(\rho^3-t^2A^2)^{s-2}\\
 &=x^{s-2}C^{s-2}.
 \end{aligned}
\]

Thus

\[
 P_{1,S,0}=Ax^{s-2}C^s=P_s,
 \qquad N=3s,
\]

which is exactly the infinite family in the original witness note.

### 6.3 Radial padding

For every \(k\ge6\), take

\[
 r=1,\qquad S=1,\qquad h=k-6.
\]

Then

\[
 P_{1,1,k-6}=\rho^{k-6}AC^2=P_k,
 \qquad N=k,
\]

which is the padded family used in the
[homogeneous spillover theorem](GVC3_HOMOGENEOUS_SPILLOVERS.md).

### 6.4 Rank and power transport

Let \(q\) be a quadratic form of rank at least three and let \(k\ge6\).
After scalar extension and a linear change of variables, choose a
nondegenerate ternary summand on which the associated second-order operator
is \(\Delta\).  Write the full operator as

\[
 D_q=\Delta+D_\perp.
\]

For every triple

\[
 r\ge1,\qquad e,h\ge0,\qquad 6r+3e+h=k,
 \tag{6.2}
\]

and every degree-\(e\) profile satisfying the moment hypothesis, extend
\(P_{r,S,h}\) by making it independent of the complementary variables.
Every term containing \(D_\perp\) then vanishes.  The pure cancellation,
multiplier ladder, shifted-power detectors, and exact depth are therefore
unchanged.

> **Corollary 6.1 — rank transport with profile strata.**  For every
> \(k\ge6\), a quadratic symbol of rank at least three carries all the
> cusp-profile counterexample strata indexed by (6.2).  Each stratum has
> the open real profile cone \(S>0\) on \([0,1]\).  Rank at most two carries
> no homogeneous GVC counterexample, because \(q^k\) is a split symbol
> after scalar extension.

Thus the existing rank-\(3\) threshold is strengthened from existence of
one padded witness to explicit profile-moduli strata with the maximal
power-tail and exact-depth conclusions.

The ladder, shifted-power failure, and exact depth in the companion notes
are precisely the corresponding special cases of (1.8), (4.3), and (5.1).

## 7. The transverse cusp defect

Put

\[
 U=\rho,\qquad V=tA.
\]

Then

\[
 \boxed{V^2-U^3=-xC.}
 \tag{7.1}
\]

On the divisor \(x=0\),

\[
 U=t^2,\qquad V=t^3.
 \tag{7.2}
\]

Thus \(x=0\) maps to the normalization

\[
 t\longmapsto(t^2,t^3)
\]

of the \(A_2\)-cusp \(V^2=U^3\).  Equation (7.1) has the following precise
interpretation:

- \(U^3-V^2\) is the cusp equation;
- \(x\) cuts out the transverse divisor;
- \(C=(U^3-V^2)/x\) is its first transverse defect;
- \(C^{2r}\) cancels the phase-square winding \(x^{-2r}\);
- \(A^r\) supplies the adjacent Taylor exponent forced by the Hopf
  classification.

The geometry also gives a converse to the construction.  On \(\rho=1\),
put \(B=1+x^2\) and consider an arbitrary normalized phase-square
one-profile expression

\[
 \Phi=x^{-2r}B^rR(t^2B^2).
 \tag{7.3}
\]

The endpoint factor \(1-t^2B^2=xC\) has exact generic \(x\)-adic order one
along the divisor \(x=0\) of the quadric.  Hence polynomiality of (7.3) is
equivalent to

\[
 \operatorname{ord}_{z=1}R(z)\ge2r.
 \tag{7.4}
\]

Write \(R(z)=(1-z)^{2r}S(z)\), with \(e=\deg S\).  Then

\[
 \Phi=B^rC^{2r}S(t^2B^2).
 \tag{7.5}
\]

Standard homogenization of (7.5) gives exactly (1.3)--(1.4).  More
explicitly, if a homogeneous lift has degree \(2N\), then

\[
 h=N-6r-3e\ge0
 \tag{7.6}
\]

and the lift is \(P_{r,S,h}\).  The lift is unique in its degree: a
homogeneous polynomial vanishing on \(\rho=1\) vanishes on the dense set of
scalar translates of that quadric and is therefore zero.

> **Proposition 7.1 — exhaustive normalized cusp lift.**  Up to the sphere
> torus and nonzero scaling, the polynomials \(P_{r,S,h}\) are precisely
> the homogeneous lifts in the same three variables of the normalized
> phase-square one-profile class (7.3).

Consequently every such lift satisfies

\[
 N=6r+3e+h\ge6r.
 \tag{7.7}
\]

For primitive winding \(r=1\), the first possibility is \(N=6\).  Equality
forces \(e=h=0\), so \(S\) is constant and the normalized polynomial is
\(AC^2\), up to scaling and the sphere torus.  This recovers the existing
minimum theorem and closes the arbitrary-profile and radial-suspension
loopholes inside that architecture.  It remains a scoped minimum, not a
global lower bound for all homogeneous ternary counterexamples.

## 8. The quadric quotient and the dimension resource

Let \(n\ge2\) and write a split nondegenerate quadratic form as

\[
 q_n=xy+q_0(t_1,\ldots,t_{n-2}).
\]

Let \(\mathbb G_m\) act by

\[
 x\mapsto\lambda x,
 \qquad
 y\mapsto\lambda^{-1}y,
 \qquad
 t_i\mapsto t_i.
 \tag{8.1}
\]

On the affine quadric \(q_n=1\), exactness of invariants for the torus gives

\[
\begin{aligned}
 \left(k[x,y,t_1,\ldots,t_{n-2}]/(q_n-1)\right)^{\mathbb G_m}
 &\cong
 k[xy,t_1,\ldots,t_{n-2}]/(q_n-1)\\
 &\cong k[t_1,\ldots,t_{n-2}].
\end{aligned}
\tag{8.2}
\]

> **Theorem 8.1 — quadric resource count.**  For \(n\ge2\), the angular
> quotient of a split rank-\(n\) quadric has dimension \(n-2\).

Call this dimension the **Beta rank**.  With the harmless convention that
the one-variable case has Beta rank zero,

\[
 \beta(n)=\max(n-2,0).
 \tag{8.3}
\]

For \(n=2\), the quadric is \(xy=1\), its invariant ring is \(k\), and the
angular quotient is a point.  For \(n=3\), the invariant ring is \(k[t]\):
\(t\) is the first genuine height/Beta parameter.  The three coordinate
resources used by (3.1) are therefore

\[
\begin{array}{c|c}
\text{resource}&\text{coordinate}\ \text{role}\\ \hline
\text{angular extraction}&x\\
\text{opposite-weight divisor}&y\\
\text{free quotient parameter}&t.
\end{array}
\tag{8.4}
\]

De Bondt's split-symbol theorem proves homogeneous GVC in dimensions one
and two.  The specialization \(P_{1,1,0}\) disproves it in dimension three,
and unused-variable padding handles all larger dimensions.  Combining those
independent theorems with (8.3) gives the exact numerical equivalence

\[
 \boxed{
 \operatorname{HGVC}(n)\text{ holds}
 \quad\Longleftrightarrow\quad
 \beta(n)=0
 \quad\Longleftrightarrow\quad
 n\le2.}
 \tag{8.5}
\]

Equation (8.5) packages the proved dimension classification; quotient
dimension alone is not being asserted as a general causal criterion for an
arbitrary moment problem.

## 9. What this says, and does not say, about GVC(2)

A nonhomogeneous binary operator

\[
 \Lambda=\Lambda_{d_1}+\cdots+\Lambda_{d_s}
\]

introduces selection counts

\[
 \alpha_1+\cdots+\alpha_s=m.
 \tag{9.1}
\]

After normalization, \(\alpha/m\) ranges over a simplex of dimension
\(s-1\).  Such moving frequencies can imitate a height parameter even
though the homogeneous binary quadric has Beta rank zero.  The imitation is
not measure-theoretically exact: its rows carry multinomial, contraction,
and radial factorial weights rather than the flat rank-three Reynolds law.

This gives a useful interpretation of the existing
[Beta--torus Hall problem](BINARY_GVC_UNIFORM_FACE_TERMINATION.md): the only
possible latitude in that Hall-shell programme was synthetic and factorially
weighted.  A cusp-oriented route target is therefore:

> **Synthetic-latitude obstruction target (open).**  Prove that a
> scale-compatible binary Hall packet cannot simultaneously reproduce the
> endpoint binomial jet, retain the opposite-weight divisor mark, and satisfy
> every factorially weighted pure cancellation unless it acquires a
> split-symbol separator or loses support.

This is not proved here, and it was never shown equivalent to unrestricted
GVC(2).  Such an equivalence would require a reduction showing that every
surviving binary packet must imitate the cusp mechanism.  The Hall-specific
promotion problem remains a well-defined target inside the parked route:
inherit a positive component or enough common Cartesian marks from the
prime-dependent signed shell over one common high quotient.  The cusp
resource count explains that target but does not eliminate any of its
factorial or carry hypotheses.  The independent
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md) proves unrestricted
binary GVC without this reduction.

## 10. Scope and literature boundary

The external binary positive input is Michiel de Bondt's
[*A few remarks on the Generalized Vanishing
Conjecture*](https://arxiv.org/abs/1206.2836), which proves GVC for products
of linear forms and hence for every homogeneous binary operator.  The
endpoint coefficient mechanism originates in Christopher D. Long's
[*Small Counterexamples to the Gaussian Moments
Conjecture*](https://arxiv.org/abs/2607.18186) and is classified in the
repository's Hopf note.

A search of the available primary literature on 2026-08-03 located neither
the homogeneous cusp-profile suspension (1.4) nor the geometric packaging
(7.1).  This note records the result as a repository theorem but makes no
priority claim before independent literature search, expert review, and
external refereeing.

Nothing here gives a counterexample for the ordinary second-order Laplacian
conjecture.  The operators are \(\Delta^N\) with \(N\ge6\).  The resource
interpretation is not the proof of unrestricted GVC(2); that proof is the
separate Hall-envelope theorem.

## 11. Reproduction

Run

```bash
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
```

The checker verifies the cusp identity, all homogeneous specializations,
five exact profile ladders through \(m=5\), top Laplacian contractions and
terminal traces for winding, profile, and radial examples, and two direct
shifted-power detectors for the non-power profile \(S=1+z\).  It also
enumerates the admissible rank/power strata through \(k=18\) and checks
that a trailing profile zero is exactly a three-step radial suspension.
It writes
`artifacts/generated-results/gvc3_cusp_profile_suspension.json`.

The bounded loops are regressions.  The all-order statements are the
Reynolds, Taylor-gap, and Fischer calculations above.
