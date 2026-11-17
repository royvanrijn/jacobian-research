# A growing selfadjoint square has no finite scalar mate

Let k have characteristic zero, let eta in k^times, and let d,p be
positive integers. Use the finite scalar differential algebra

\[
 \mathcal E=k(V)[F;D],\qquad
 D=\eta V^{1-d}\partial_V.
\]

Suppose U,W are rational functions of V, with

\[
 U=uV^p+\text{terms of smaller degree at infinity},\quad
 u\in k^\times,\qquad \deg_\infty W<2p.                  \tag{1}
\]

The case W=0 is included, with degree minus infinity. Put

\[
                         L=F^2+U,\qquad P=L^2+W.         \tag{2}
\]

There is no finite T over k(V), or any finite algebraic coefficient
cover, satisfying [P,T]=kappa for a nonzero scalar kappa in k.
The theorem has no partner-order bound and imposes no parameter
regularity assumption. It is an obstruction to scalar finite operators.
It is not a polynomial-centralizer theorem: when W=0, L commutes
with P=L^2 but is not polynomial in P.

The proof uses a nonzero infinity coefficient for every odd positive
fractional flow. It therefore treats all finite partner orders by one
argument, rather than by a finite list of string equations.

## Adjoint projection and the complete odd-flow normal form

The coefficient-fixing antiinvolution F^dagger=-F is well defined,
since it preserves [F,f]=D(f). Both L and P are selfadjoint.
If [P,T]=kappa, then [P,T^dagger]=-kappa, so
(T-T^dagger)/2 is a nonzero skewadjoint mate with the same commutator.
It suffices to exclude such a mate over k(V).

In the formal pseudodifferential algebra, multiplication is

\[
 F^n f=\sum_{j\ge0}{n\choose j}D^j(f)F^{n-j}.
\]

Operators are bounded above in F, and every specified coefficient is
a finite sum. There is a unique monic fourth root R=P^{1/4} with
leading term F. The coefficient recursion divides only by four.
Since (-R^dagger)^4=P and -R^dagger has leading term F, uniqueness
gives R^dagger=-R. Adjoint preserves positive and strictly negative
F-orders, so for m>=1

\[
                  B_m=(P^{m/4})_+,\qquad
                         B_m^\dagger=(-1)^m B_m.         \tag{3}
\]

The commutator [P,B_m] has order at most two, because its negative-part
expression has upper order 4-1-1=2.

For any skewadjoint finite operator T with [P,T] of order at most two,
its highest coefficient is constant: the coefficient at F^{M+3} is
4D(t_M). Its order M must be odd by skewadjointness. Subtract that
constant times B_M; the remainder stays skewadjoint, its commutator
still has order at most two, and its order strictly drops. An order-zero
skewadjoint coefficient is zero. Consequently every possible mate
has the exhaustive finite normal form

\[
                         T=\sum_{\substack{m\ge1\\ m\ {\rm odd}}}c_mB_m,
                                  \qquad c_m\in k.       \tag{4}
\]

No constant lower flow has been discarded. Finite repetition proves
(4) for arbitrary rational coefficient shapes and degrees.

## Weight bounds and the slow-potential gap

Give V weight one and F weight p/2. A nontrivial Ore contraction loses

\[
                              \tau=d+p/2>0               \tag{5}
\]

units of weight: D lowers V-degree by at least d, while a contraction
also removes one F. Thus for finite differential operators

\[
             \operatorname{wt}[A,B]
                 \le\operatorname{wt}A+\operatorname{wt}B-\tau.
                                                                  \tag{6}
\]

The zeroth-contraction terms cancel because position coefficients
commute. This is a bound on finite coefficients and does not evaluate
an unbounded spatial series.

Set B_m^L=(L^{m/2})_+. The square-root and fourth-root recursions imply

\[
 \operatorname{wt}B_m^L\le mp/2,\qquad
 \operatorname{wt}(B_m-B_m^L)\le mp/2-\epsilon,
 \quad\epsilon=2p-\deg_\infty W>0                         \tag{7}
\]

when W!=0; for W=0 their difference is zero. Here is the coefficientwise
justification of the second bound. The monic-root coefficients are
universal differential polynomials in U,W. At W=0 the fourth root is
L^{1/2}. Every coefficient of the difference contains at least one W
or a derivative of W. Assigning U weight p and W weight 2p gives the
usual upper weight for the root; using the actual W-degree loses at
least epsilon. Derivatives and nontrivial contractions lose further
positive weight. The recursion only divides by nonzero numerical
constants, so this weight induction holds at every negative F-order.
Raising the roots to the finite integer power m and then taking the
positive part proves (7).

## The nonzero quadratic-flow coefficient

For odd m>=1, put a=(m+1)/2. The quadratic positive flow has a scalar
commutator with L:

\[
 C_m=[L,B_m^L]=-2D\bigl(\operatorname{res}_F L^{m/2}\bigr).
                                                                  \tag{8}
\]

Indeed the coefficient of F^0 in [L,rF^{-1}+...] is 2D(r), and the
negative-part commutator has order at most zero. Its leading residue
at infinity is

\[
 \operatorname{res}_F L^{m/2}
      ={m/2\choose a}u^aV^{ap}
                          +\text{terms of smaller degree}.        \tag{9}
\]

The zero-contraction binomial expansion gives the displayed coefficient:
F^{m-2a}=F^{-1}. Any derivative term loses at least tau relative to
that weight, and every lower term of U also lowers it. The binomial
coefficient is nonzero for every positive odd m. Since ap>0, applying
D gives

\[
 C_m=-2\eta ap{m/2\choose a}u^aV^{ap-d}
                           +\text{terms of smaller degree}.       \tag{10}
\]

Now [L^2,B_m^L]=LC_m+C_mL, whose F^2 coefficient is 2C_m.
The other contributions to [P,B_m] are [W,B_m^L] and
[P,B_m-B_m^L]. By (6)--(7), both have F^2 coefficient of degree
at most ap-d-epsilon, strictly below (10). Therefore

\[
 [P,B_m]_{F^2}
    =-4\eta ap{m/2\choose a}u^aV^{ap-d}
                              +\text{terms of smaller degree}.    \tag{11}
\]

This coefficient never vanishes. In a nonempty finite sum (4), the
largest odd m has strictly greatest displayed V-degree: increasing
m by two increases ap-d by p>0. Lower flows cannot cancel it.
A nonzero coefficient of F^2 contradicts a scalar commutator.
This proves the theorem over k(V).

For a finite algebraic coefficient extension, first take the normalized
[coefficient trace](UNIVERSAL_P_POWER_OBSTRUCTION.md). The uniquely extended derivation commutes with trace,
so a finite cover mate would give a rational mate with the same nonzero
scalar commutator. No adjoint on the cover is needed.

The conditions in (1) are substantive. Outside the growing-potential
hypothesis, the ordinary Weyl pair F^4+V,-F has commutator eta when
D=eta partial_V. The proof also makes no assertion when W reaches
weight 2p. A finite prefix of flow coefficients is not being used to
cover those cases.

## The ordered resonant critical face

Consider [A,V]=eta/V and

\[
 q=A+\beta(V^{-4}A)^2,\qquad P_{\rm res}=gq^2+\mu V,
                    \qquad\beta,g,\eta\ne0.             \tag{12}
\]

Set E=V^{-4}A. Then [E,V]=eta V^{-5} and A=V^4E, so
q=beta E^2+V^4E. The rational coefficient-preserving shift
F=E+V^4/(2beta) gives exactly

\[
 q=\beta F^2-\frac{V^8}{4\beta}-2\eta V^{-2}.
\]

Thus P_res/(g beta^2) is of the form (2) with d=6, p=8,

\[
 U=-\frac{V^8}{4\beta^2}-\frac{2\eta}{\beta}V^{-2},
                       \qquad W=\frac{\mu V}{g\beta^2}.
                                                                  \tag{13}
\]

The degree of W is at most one, below 2p=16, so no finite scalar mate
exists. This remains true when the lowest divisor character vanishes
at both primitive twenty-fourth and primitive twelfth roots. That
vanishing only prevents the earlier local-divisor proof; it does not
contradict the infinity argument.

## The actual globally admitted candidate is excluded directly

Use the global A_0,C,S_b of
[the global admission note](GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md),
and put \(\mathsf L=X^4C\), \(W=\mathsf L^2\). These are actual finite operators
commuting with R. Take lambda to be a nonzero h-independent scalar.
Since [C,v]=rh/v,

\[
                 v\mathsf L=X^3C,\qquad
                 \mathsf L v=X^3C+rhX^5,
\]

both \(vW=(v\mathsf L)\mathsf L\) and \(Wv=\mathsf L(\mathsf L v)\) are globally
admitted. Hence so is E_0=vW+Wv and

\[
 S_\lambda=S_b-3\lambda h^3(RE_0+CW+WC)
                                  -3R\lambda^2h^6W^2,
                                      \quad\lambda\ne0. \tag{14}
\]

It has the prescribed Long classical reduction. With z=a-2/(3r) and
C=r(z-v), its exact generic expression is

\[
 S_\lambda-\frac4{3r}=-3r\bigl(z+\lambda h^3W\bigr)^2+2v.
                                                                  \tag{15}
\]

This full generic operator, not only its critical reduction, is covered
by the theorem. Set

\[
 E=v^{-4}(z-v),\qquad \beta=\lambda h^3r^2,\qquad
 F=E+\frac{v^4}{2\beta}.
\]

Then [F,v]=h v^{-5}, W=r^2E^2, and

\[
 z+\lambda h^3W
 =\beta E^2+v^4E+v
 =\beta F^2+v-\frac{v^8}{4\beta}-2hv^{-2}.               \tag{16}
\]

Dividing (15) by -3r beta^2 yields (2) with d=6, p=8,

\[
 U=\frac v\beta-\frac{v^8}{4\beta^2}
                         -\frac{2h}{\beta}v^{-2},\qquad
 W_{\rm slow}=-\frac{2v}{3r\beta^2}.                    \tag{17}
\]

The transformations between a and F are invertible affine first-order
Ore changes with rational position coefficients. They preserve finite
differential order. The nonexistence conclusion therefore applies in
the actual generic relative algebra, even allowing arbitrary finite
algebraic position covers. QEP1 localizes any globally admitted mate
commuting with exact R into that algebra. Consequently (14) has no
such finite mate, and cannot supply a bounded full Long lift.
No formal-parameter limiting argument or centralizer-finiteness
hypothesis is used in this application.

## Scope

The obstruction concerns a scalar completed square L^2 plus a strictly
slower rational position potential, after an actual rational first-order
change of momentum. It supplies no general polynomial-centralizer
assertion and no conclusion for arbitrary nonscalar positions, for
unrelated transport relations, or for added terms reaching the top
weight or changing the square-plus-position structure. No finite search,
formal verification, independent whole-proof replay, external review or
novelty claim is made. The proof is coefficientwise formal algebra used
to exclude actual finite mates at every order.
