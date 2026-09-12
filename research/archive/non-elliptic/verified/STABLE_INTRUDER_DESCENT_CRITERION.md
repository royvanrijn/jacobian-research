# Stable intruder descent and physical inertia

This note isolates the reusable argument behind the fixed-quintic and
all-rank quadratic-gauge stabilizer theorems.  It separates three logically
different gates:

1. descent from a stabilized polynomial ring to the physical polynomial
   ring;
2. identity on the physical target; and
3. identity on the physical source.

Only the first gate follows from the Newton polytope.  The second also needs
faithfulness of the normalized boundary decoration, and the third needs a
trivial deck group.  Keeping these hypotheses separate prevents the stable
intruder argument from being used as a global automorphism theorem.

Work over a characteristic-zero field `k`.  Put

\[
 R=k[x_1,\ldots,x_n],\qquad
 R^{[s]}=R[t_1,\ldots,t_s],
\]

and let `H in R` be nonzero.  A vertex of the Newton polytope of `H` is an
**intruder** if all of its coordinates are positive.

## 1. Stable base descent

> **Stable intruder descent criterion.**
>
> Suppose `H` has an intruder and
>
> \[
>  T\in\operatorname{Aut}_k(R^{[s]}),\qquad T(H)=\rho H
>  \quad(\rho\in k^\times).
> \]
>
> Then
>
> \[
>  T(R)=R.
> \]

For each `a`, let `tau_a` be the standard additive-group action translating
`t_a` and fixing the other variables.  Both conjugates

\[
 T\tau_aT^{-1},\qquad T^{-1}\tau_aT
\]

fix `H`.  If either conjugate acted nontrivially on `R`, then `H` would be a
stable additive-group invariant of `R`.  Kuroda's stable-invariant theorem
says that such an invariant has no intruder.  Hence both conjugates fix `R`
pointwise.  Equivalently, `T(x_i)` and `T^{-1}(x_i)` are independent of every
`t_a`.  Thus

\[
 T(R)\subseteq R,\qquad T^{-1}(R)\subseteq R,
\]

which proves equality.

This step allows arbitrary polynomial degree and any finite number of
identity stabilizations.  It does not say that `T|_R` is the identity.

## 2. The physical identity gate

The remaining unstabilized step has a similarly reusable form.

> **Boundary-faithful intruder criterion.**
>
> In addition to the hypotheses above, suppose that the induced
> automorphism of `R/(H)` is the identity.  Then
>
> \[
>  T|_R=\operatorname{id}_R.
> \]

By stable base descent, `T|_R` is an automorphism of `R`.  Identity modulo
`H` gives

\[
 T(x_i)=x_i+HV_i,\qquad V_i\in R.                     \tag{2.1}
\]

Assume that some `V_i` is nonzero.  Choose a weight in the open normal cone
of an intruder vertex `d` of `H`, perturbed so that it exposes one vertex
`v` of `V_i`.  Minkowski additivity of Newton polytopes makes `d+v` a vertex
of `HV_i`.  It has every coordinate positive and cannot be cancelled by the
monomial `x_i`.  Thus `T(x_i)` has an intruder.  This contradicts the
Derksen--Hadas--Makar-Limanov coordinate-polynomial theorem, which says that
a coordinate of a polynomial ring has no intruder.  Therefore every `V_i`
vanishes.

In applications, “identity modulo `H`” is not automatic.  It must be
supplied by an intrinsic boundary normalization whose residual symmetry
acts faithfully.  Primeness of `H` is often what makes preservation of the
boundary divisor imply `T(H)=rho H`; it is not an extra hypothesis in the
two algebraic implications once that equality is known.

## 3. The source and marked-orbit gates

Let `F:X->Y` be a dominant generically finite polynomial map, and suppose a
stable self-equivalence has already been shown to act identically on the
physical target.  Any resulting physical source restriction is then a deck
transformation of the generic cover.  Consequently,

\[
 \operatorname{Deck}(F_\eta)=1
 \quad\Longrightarrow\quad
 \text{the physical source restriction is the identity}.       \tag{3.1}
\]

For a transitive monodromy group `G` with point stabilizer `G_1`, the deck
group is

\[
 N_G(G_1)/G_1.                                        \tag{3.2}
\]

Thus the universal `S_N` cover has trivial deck group because
`N_{S_N}(S_{N-1})=S_{N-1}`.

Combining the three gates gives the following useful formulation.

> **Physical marked-rigidity corollary.**
>
> Suppose a stable marked self-equivalence preserves a prime boundary
> `(H)`, the polynomial `H` has an intruder, the normalized boundary
> decoration has trivial stabilizer, and the generic cover has trivial deck
> group.  Then the self-equivalence is the identity on every physical source
> and target coordinate after every identity stabilization.

Automorphisms involving only the added identity factors may remain.  If
`Aut_vert` denotes the kernel of restriction to the physical source and
target, the conclusion is the pointwise statement

\[
 \operatorname{Aut}_{\rm st}^{\rm marked}/
 \operatorname{Aut}_{\rm vert}=1.                     \tag{3.3}
\]

This is the precise sense in which the physical inertia is trivial modulo
vertical stabilization gauge.

## 4. Scope

The criterion proves a stabilizer statement for objects satisfying its
hypotheses.  By itself it does not:

- construct a stack of all stable polynomial maps;
- classify the vertical automorphism group;
- prove separatedness or finite type for a global quotient; or
- make a presentation-dependent compiler descend through Tschirnhaus
  equivalence.

The clean quadratic-gauge application, including its uniform intruder and
faithful boundary action, is in
[Stable moduli of the quadratic-gauge family](QUADRATIC_GAUGE_STABLE_MODULI.md).
The explicit quotient receiver and its global weight-one slice are in
[The clean quadratic-gauge decorated receiver](QUADRATIC_GAUGE_DECORATED_RECEIVER.md).

## External inputs

- Shigeru Kuroda,
  [*Initial forms of stable invariants for additive group actions*](https://arxiv.org/abs/1304.0313),
  Theorem 1.1.
- Harm Derksen, Ofer Hadas, and Leonid Makar-Limanov,
  [*Newton polytopes of invariants of additive group actions*](https://sites.lsa.umich.edu/hderksen/wp-content/uploads/sites/614/2018/05/A.I.a.8.pdf),
  Corollary 3.3.
