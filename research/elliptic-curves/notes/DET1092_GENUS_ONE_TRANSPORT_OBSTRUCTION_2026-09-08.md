# Genus-one carrier transport: geometric obstruction and a degree-six field

## Exact result and scope

**New verified application.** For the fixed first302 genus-one carrier
`C`, and for all nine unchanged generic-point controls, write `J=Jac(C)`
and `E` for the original elliptic fibre at the marked parameter. Then

\[
\operatorname{Hom}_{\overline{\mathbf Q}}(J,E)=0,
\qquad
\operatorname{Hom}_{G_{\mathbf Q}}(J[2],E[2])=0.
\]

Both statements hold in reverse. Thus lowering the carrier genus from two
to one does **not** enable a homomorphism-based transport of its class to
the original fibre. No nonconstant map `C -> E` exists over any number
field. The incidence map `C -> X` to the varying elliptic surface remains
valid; it is not a map to one fixed fibre.

**New exact construction.** A displayed zero-dimensional scheme of degree
six parametrizes all identifications of these two2-torsion modules. In
each pair it is the spectrum of a degree-six field, not six rational
points. A torsion identification exists over that field, but still no
elliptic isogeny. This is **not a seed cover** or a prospective seed rule.

**Verified input boundary.** Reuse only the carrier equations from the
[full Picard comparison](DET1092_GENUS_ONE_PICARD_SPECIFICITY_2026-09-08.md):
one historically calibrated first member and nine source-only members
`z=z_0(tau)`. No marked exceptional coordinate, subsequent direction or V3
artifact enters this calculation. The first member remains retrospectively
calibrated; stripping its point coordinates does not remove that provenance.

## Small independent nonisogeny certificate

**Verified application.** At43 the exact reduced monic models are

\[
J:\ y^2=x^3+37x^2+24x+12,
\qquad E_{302}:\ y^2=x^3+5x^2+26x+38.
\]

Their orders are56 and55. The traces are therefore `-12,-11`, and the
Frobenius discriminants are `-28,-51`. Both reductions are ordinary.
Their geometric endomorphism algebras are the distinct fields
`Q(sqrt(-7))` and `Q(sqrt(-51))`.

**Established literature, verified application.** An ordinary reduction has
geometric endomorphism algebra `Q(sqrt(a_p^2-4p))`; see
[Sutherland, Theorem14.5 and Corollary14.7](https://math.mit.edu/classes/18.783/2017/LectureNotes14.pdf).
Homomorphisms specialize injectively at good reduction; see
[Conrad, Proposition6.4](https://math.stanford.edu/~conrad/DarmonCM/2011Notes/SemistableReduction.pdf),
applied also after finite extension. Isogenous elliptic curves have
isomorphic rational endomorphism algebras by conjugating with an isogeny.
The two distinct quadratic fields therefore exclude a geometric isogeny
of **any** degree, not merely one over the prime field.

**Verified controls.** The same test gives the following exact witnesses.
The final two columns are Frobenius discriminants, not global curve
discriminants or maximal-order discriminants.

| Pair | Prime | Carrier trace | Original trace | Carrier discriminant | Original discriminant |
| --- | ---: | ---: | ---: | ---: | ---: |
| First302 | 43 | -12 | -11 | -28 | -51 |
| scale-0131232 | 17 | -6 | -8 | -32 | -4 |
| scale-0257585 | 67 | -14 | -15 | -72 | -43 |
| scale-0487239 | 19 | -4 | -8 | -60 | -12 |
| scale-0177036 | 17 | -5 | -8 | -43 | -4 |
| scale-0043332 | 17 | -6 | -8 | -32 | -4 |
| scale-0590501 | 43 | -11 | -8 | -51 | -108 |
| scale-0290097 | 17 | -5 | -8 | -43 | -4 |
| scale-0748009 | 17 | -6 | -8 | -32 | -4 |
| Generic-point302 | 17 | -6 | -8 | -32 | -4 |

## Independent torsion fields

**New verified deduction.** In every pair, irreducible cubic and `1+2`
reductions certify both cubic splitting groups as `S3`. A good prime
where the discriminant quadratic characters differ proves that their
quadratic subfields are distinct. All witnesses are retained, without
factoring any global discriminant or constructing a number field.

The intersection of two `S3` splitting fields is Galois. A nontrivial
common quotient of `S3` is `C2` or `S3`, and either would identify their
unique quadratic subfields. The intersection is therefore `Q`, giving

\[
G_{\mathbf Q}\longrightarrow
\operatorname{GL}_2(\mathbf F_2)\times\operatorname{GL}_2(\mathbf F_2)
\quad\text{with full image }S_3\times S_3.
\]

For first302, carrier irreducibility is witnessed at61, original
irreducibility at43, and the respective transpositions at43 and47.
At43 the factor types are `1+2` and `3`, giving the unequal quadratic
characters already at the nonisogeny prime.

**New finite-linear deduction.** On the16 linear maps between the two
two-dimensional modules, the independent left/right action has orbits of
sizes `1,9,6`, for ranks `0,1,2`. Only zero is fixed. The checker enumerates
all36 pairs of invertible matrices and all16 maps; it does not improperly
align separately selected Frobenius conjugacy representatives.

For any number field `L`, a rank-one equivariant map would have coefficient
field of degree9, and a rank-two map degree6. Consequently

\[
\text{rank-one map}\Rightarrow9\mid[L:\mathbf Q],\qquad
\text{rank-two map}\Rightarrow6\mid[L:\mathbf Q].
\]

In particular, extensions of degree at most five cannot create a nonzero
torsion-module transfer. No extension creates an elliptic isogeny, by the
separate geometric obstruction above.

## Explicit minimal torsion-identification scheme

**New equation construction.** Write the two monic cubic models as

\[
f_J(x)=x^3+ax^2+bx+c,\qquad f_E(x)=x^3+Ax^2+Bx+C,
\]

using the exact rational coefficient triples in the frozen equation files.
Set

\[
T=\begin{pmatrix}0&0&-c\\1&0&-b\\0&1&-a\end{pmatrix},
\quad N=uI+vT+wT^2,
\qquad \boxed{\det(sI-N)=s^3+As^2+Bs+C.}
\]

Equating the three coefficients is an explicit system of degrees one,
two and three in `u,v,w`. Its expanded equations are retained in the
[scheme certificate](../../artifacts/generated-results/elliptic-curves/det1092_genus1_transport_gate_v1/torsion-isomorphism-scheme.json).
The induced algebra map sends the original cubic root to
`u+v*theta+w*theta^2` in the carrier cubic algebra.

**New exact deduction.** Over the algebraic closure, evaluation at the
three distinct carrier roots is an invertible Vandermonde transformation.
The characteristic-polynomial condition says that their three images are
the three distinct original roots, in some order. Hence the scheme has
exactly six reduced geometric points. The full product `S3 x S3` acts
transitively on these bijections, so its coordinate algebra is a degree-six
field over `Q`. It has no rational point. Its field is a minimal-degree
extension admitting any nonzero torsion-module map, since the other
nonzero-map orbit has size9. No field arithmetic is required for this
degree and obstruction proof.

## Meaning for the seed question

**Verified conclusion.** The first carrier class is rational and outside
its inherited rank12 subgroup, but that fact cannot be pushed through a
curve homomorphism or a coefficient-module map to certify the original
MW17 quotient. The same obstruction holds for the nine original-generic
controls. Together with the completed specificity counterexamples, this
closes the proposed linear-transfer repair for these selected carriers.

**Unresolved.** A successful criterion must use additional marked nonlinear
incidence or the original elliptic Kummer class itself. These results do
not exclude either, all possible higher descents, different carriers, or
an empirical global arithmetic correlation. No full Selmer group, nonzero
Sha class, new elliptic direction or prospective302 producer is claimed.

## Reproduction and limits

**Verified computation.** The protocol freezes all ten pairs and the same
64 primes before outcomes. All640 exposures, including denominator and
bad-reduction skips, are retained. The constructor uses finite-field
factorization and point counting; the independent checker instead counts
roots and square residues with integer arithmetic and rederives the exact
models from the saved quartics and original equation. The only polynomial
matrix construction is the displayed3-by3 universal identity.

```sh
sage -python research/elliptic-curves/cas/construct_det1092_genus1_transport_gate.sage
sage -python research/elliptic-curves/cas/verify_det1092_genus1_transport_gate.sage
```

Both runs finish under one second, separately capped at25 seconds. The
[independent replay](../../artifacts/generated-results/elliptic-curves/det1092_genus1_transport_gate_v1/replay.json)
binds equations, outcomes, source files and checker hashes. No original
parameter sweep, point search, class group, full Selmer computation, paid
backend, pilot change or detached process was used.
