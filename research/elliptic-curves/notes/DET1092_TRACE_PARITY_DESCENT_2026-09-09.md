# Trace parity is the obstruction to descending a quadratic branch class

**All-prime extension.** The subsequent
[division and multisection theorem](DET1092_ALL_PRIME_DIVISION_AND_MULTISECTION_THEOREM_2026-09-09.md)
proves full saturation under genus below9, not only2-saturation. Every
genuine bisection, of any genus, is generically independent. These are
function-field statements; specialized seed incidence remains separate.

## Result: the RR branch and a generic Kummer class are different objects

**New deduction in this application, using established Kummer theory.**
Let `F=K(sqrt(d))` be a quadratic field extension, `Q in E(F)`, and
`Z=Q+sigma(Q) in E(K)`. Under the canonical identification of quadratic
twist2-torsion with `E[2]`, one has the exact criterion

\[
\boxed{\delta_F(Q)\in\operatorname{res}_{F/K}H^1(K,E[2])
       \quad\Longleftrightarrow\quad Z\in2E(K).}
\tag{1}
\]

The right side depends on the **trace parity**, not a ramification budget.
When it fails, allowing arbitrary ramification or larger coefficients
cannot produce a generic cubic norm-square class whose restriction is the
actual RR branch class. Translating `Q` by any generic point does not repair
it: the trace changes by twice that point.

On the determinant1092 parent, the historical norm10 centre, the norm8
first-seed carrier, and the orbit8044 conic all have nonzero trace parity.
Their generic branch classes therefore do **not** descend, even though
their split rational specializations can give either productive or inherited
points. This is an obstruction to a proposed *class continuation*, not to
the existing rank18 constructions.

**New geometric consequence of the completed halving argument.** For every
nonzero parity in the full generic `M17/2M17`, the curve of halves of a
representative section has degree4 and genus9. Therefore a nonconstant base
change from a curve of genus less than9 cannot make that trace even. In
particular neither a rational nor an elliptic base change can repair this
descent obstruction by adjoining half the trace. Constant number-field
extensions cannot repair it either, since the full geometric generic MW
group is already the displayed rational group.

## 1. Audit: what is already proved and what is added?

**Verified prior applications.** The
[trace/twist identity](DET1092_TRACE_TWIST_KUMMER_OBSTRUCTION_2026-09-08.md)
proves that the anti-trace's twist Kummer class equals the inherited trace
class. The
[first-centre singular-member theorem](DET1092_RR_NET_SINGULAR_MEMBER_GATE_2026-09-08.md)
already proves the genus-nine halving curve for the historical norm10 centre.
The [degree-two quotient](CURVE302_LOW_DEGREE_MULTISECTIONS_2026-09-07.md)
and its common chart-label interpretation already identify the finite
translation quotient with `M17/2M17`.

We do not rerun those enumerations, reconstruct the first exceptional
point, or re-prove a special first-centre genus. The additions are:

- the necessary-and-sufficient descent criterion (1), with a constructive
  even-trace direction;
- identification of the existing parity labels as the precise obstruction
  classes in (1), not labels of exceptional seed classes;
- the uniform genus-nine consequence for every nonzero parity and the
  resulting low-genus base-change obstruction;
- a precise explanation of why the recent generic norm-square lifting
  theorems do not characterize a moving odd-trace RR branch.

Novelty of the general Kummer identities is not claimed. No generic Selmer
dimension, class group, new point, or search score is computed.

## 2. Necessity: the branch has a nonzero norm

**Established literature.** Kummer maps commute with the point norm and
cohomological corestriction; restriction followed by corestriction is
multiplication by the extension degree. These are the standard norm maps
used, for example, in
[Morgan, section4, Lemma4.3](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12533).
Thus

\[
\operatorname{cor}_{F/K}(\delta_F(Q))=\delta_K(Z),\qquad
\operatorname{cor}_{F/K}\operatorname{res}_{F/K}=2=0
       \quad\text{on }H^1(K,E[2]).
\tag{2}
\]

If the branch class descends, its corestriction is zero. Injectivity of
`E(K)/2E(K) -> H1(K,E[2])` gives `Z in 2E(K)`.

**Explicit algebraic verification, independent of cohomological notation.**
Write `E:y^2=f(x)=x^3+a x^2+b x+c`, `s^2=d`, and

\[
 Q=(u+vs,w+zs),\quad \ell(x)=mx+n,\quad
 g(x)=(x-u)^2-v^2d,
\]

where `ell` is the common chord of `Q,sigma(Q)`. If
`Z=(c_Z,e_Z)`, then

\[
 f(x)-\ell(x)^2=(x-c_Z)g(x).
\]

In `A=K[theta]/f(theta)` and `B=A tensor_K F`, the branch element
`eta=u+vs-theta` satisfies

\[
\boxed{
N_{B/A}(\eta)=g(\theta),\qquad
g(\theta)(c_Z-\theta)=\ell(\theta)^2.}
\tag{3}
\]

Its norm squareclass is `delta(Z)`, not zero. If
`eta=alpha*gamma^2` for an element `alpha in A*`, taking this quadratic
norm would give a square. Thus no such `alpha` exists when `delta(Z)`
is nonzero, even before imposing a norm-square condition on `alpha` itself.
The displayed patch assumes finite nonbranch coordinates and unit
denominators; (2) covers the exceptional charts as well.

Adding an inherited class to the desired branch class changes its
corestriction by zero. Hence the obstruction also survives subtraction
of any generic point or multiplication by any restricted generic Kummer
class. Being conjugation-invariant after quotienting by inherited classes
must not be mistaken for possessing a descended representative.

## 3. Sufficiency is constructive when the trace is even

**New explicit application of the standard twist construction.** Suppose
`Z=2T` with `T in E(K)`. Put `P=Q-T`; then `sigma(P)=-P`.
On the nonbranch patch write

\[
 P=(r,s e),\qquad r,e\in K,\qquad f(r)=d e^2.
\]

The element

\[
\alpha_P=d(r-\theta),\qquad N_{A/K}(\alpha_P)=d^4e^2
\]

is a norm-square class over `K`. Since `d` becomes a square in `F`, its
restriction equals `delta_F(P)`. Multiplying by a representative of
`delta_K(T)` produces the desired descended class of `Q`. At the identity
or2-torsion charts, use the Kummer map on the quadratic twist: the same
construction is canonical through `E^(d)[2]=E[2]`. This proves the converse
of (1), not merely a necessary test.

Thus an anti-invariant point is exactly the even-trace situation. For an
odd trace, replacing the branch by its anti-trace instead gives
`R=2Q-Z`, whose twist class is inherited by the already completed identity.
That operation discards the branch information needed for a primitive seed.

## 4. The common parity universe has this arithmetic meaning

**New deduction from the completed divisor/translation formulas.** For a
degree-two divisor

\[
 C=2O+bF+\phi(w),
\]

restriction to the generic elliptic fibre identifies `C-2O` with the point
`P_w`. For a geometrically irreducible bisection with generic point `Q`,
this means `Q+sigma(Q)=P_w`. Translation by `S_x` changes this sum to
`P_w+2S_x`. Consequently

\[
\boxed{M_{17}/2M_{17}\ \xrightarrow{\delta}\
       \{\text{trace obstructions to branch-class descent}\}}
\tag{4}
\]

is precisely the same131072-element universe as the existing degree-two
translation quotient and half-lattice chart labels. The zero label is
the even-trace case. Every nonzero label obstructs descent of the generic
branch class of any bisection representing it.

In particular all40917 geometrically rational norm10 bisection orbits have
nonzero labels and this obstruction. This does **not** obstruct splitting
or new points. It says that their individual branch classes cannot be
packaged as ordinary generic Kummer classes before splitting. The complete
orbit count is reused, not recalculated, and genus-one lattice candidates
are not promoted to actual smooth bisections.

## 5. Why a low-genus base change cannot supply half the trace

**New uniform deduction from the existing first-centre proof.** Let
`Z in M17` have nonzero parity. Over `Qbar(t)` the curve
`T_Z=[2]^-1(Z)` is connected of degree4. A degree-one component would be
a generic half, contradicting the full geometric integral basis. A
degree-two component with conjugate halves `P,P'` would give the nonzero
rational2-torsion point `P+P'-Z=P'-P`, also impossible. Every nontrivial
partition of4 has a component of degree1 or2.

There is no ramification over a smooth fibre: multiplication by2 is
etale there, including when a slope chart has a pole. At an `I1` fibre,
work over a strictly henselian local trait in characteristic zero. The
section lies in the smooth Neron locus, whose special fibre is `G_m`.
Squaring is surjective on the algebraically closed residue field and
etale, so a half lifts locally by Hensel's lemma. Translation by that half
identifies the local halving torsor with `E[2]`. The `I1` inertia action
has one transposition on the three nonzero2-torsion points, and therefore
profile `(2,1,1)` on all four halves. There are24 such fibres and no others.

The same [Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
used in the completed special-centre theorem gives

\[
2g(T_Z)-2=-8+24=16,\qquad\boxed{g(T_Z)=9.}
\]

The required [etaleness of multiplication](https://stacks.math.columbia.edu/tag/0BFH)
and the semistable local argument are established tools; this is their
uniform application to the already proved full-MW17,24-I1 parent.
Replacing `Z` by `Z+2S` translates the halving cover, so this geometry is
intrinsic to the parity class.

If `C -> P1` is a nonconstant map from a smooth projective curve and a half
of `Z` exists over `Qbar(C)`, it defines a nonconstant map `C -> T_Z`.
Riemann--Hurwitz forces `g(C)>=9`. Consequently

\[
\boxed{M_{17}\cap2E(\overline{\mathbf Q}(C))=2M_{17}
             \quad\text{when }g(C)<9.}
\tag{5}
\]

This is a2-saturation statement about the inherited subgroup, not a bound
on the rank after base change. A conic or elliptic base change can and does
add independent points. It cannot do so by halving an odd generic trace.
The genus-nine curves describe halves of *centres*, never the exceptional
directions as literal halves of MW17.

For a quadratic carrier of genus below9, this also makes the conjugation
action on the branch class explicitly nontrivial:

\[
\sigma\delta_F(Q)=\delta_F(Q)+\operatorname{res}\delta_K(Z),
\qquad\operatorname{res}\delta_K(Z)\ne0.
\]

The shift vanishes only after taking the quotient by inherited classes;
that invariant quotient class still need not have a descended representative,
as the norm obstruction in section2 proves.

## 6. Apply the distinction to the completed evidence

**Verified generic-only application.** Three existing trace words are
checked in the full integral Gram lattice:

| Existing construction | Trace word / source | Norm | Branch class descends generically? |
|---|---|---:|---|
| Historical first-centre RR net | Pinned17-entry centre word |10| No |
| Norm8 first-seed carrier pencil | `P14-P15` |8| No |
| Orbit8044 conic | `-P1+P9-P13` |10| No |

Indices are zero-based. The norm8 conclusion applies to the entire
quadratic pencil wherever its member is geometrically irreducible, including
the calibrated first-seed member. No calibrated member coordinate is an
input to this proof. All three also have the elementary height obstruction
`height(Z)/4<4`; the new uniform argument is not confined to these norms.

**Verified prior controls, reinterpreted rather than recomputed.** At a
split rational specialization, the quadratic algebra becomes `Q x Q` and
the branch class is a pair

\[
\bigl(\delta(P),\delta(Z-P)\bigr),\qquad
\delta(P)+\delta(Z-P)=\delta(Z).
\tag{6}
\]

An odd trace means the two absolute classes differ. Modulo the inherited
image their classes agree, but this common residual class may be zero or
nonzero. The old productive conic fibres and the dependent conic all use
the same orbit8044 trace: its parity alone cannot distinguish them. The
completed halving-or-cycle certificates do distinguish the actual branches.
Likewise, the nine blinded M16 recoveries are new over their cores but old
over full M17; the reference subgroup, not the generic trace identity,
changes their classification. No later rank21 or V3 point is used.

## 7. Consequence for the next constructive step

**Exact obstruction, not a new predictor.** The recent unramified and
singular-fibre-only theorems concern classes already defined over `Q(t)`.
For the actual odd-trace RR branch there is a stronger problem: no such
descended class exists at all, with any ramification, if it is required to
restrict to that moving branch. A separately fitted class can still
specialize to the same point at302; the copied constant-abscissa continuation
is such a fit, not a generic continuation of the norm8 branch class.

The correct generic arithmetic object retains the quadratic extension and
the marked branch. In the notation of section2 it satisfies

\[
 [\eta]\in B^*/B^{*2},\qquad
 N_{B/F}(\eta)\in F^{*2},\qquad
 [N_{B/A}(\eta)]=\delta(Z).
\tag{7}
\]

The second norm lies in a prescribed **nonzero norm fibre**, not a norm
kernel. Equation(3) gives its explicit witness. At a split parameter one
must test a component of (6) against the specialized generic subgroup.
Its norm alone carries only the inherited trace and cannot do that job.
This is an exact reformulation of the remaining obstruction, not a claim
that solving (7) or detecting a split supplies independence automatically.

**Open endpoint.** Generic-input selection of a302 member whose split
branch escapes the inherited subgroup is still missing. The required
pointwise independence test is already available. The present theorem
rules out replacing that marked splitting problem by ordinary generic
Kummer-class descent or by a low-genus trace-halving trick.

## 8. Checkpoint and limits

**Verified computation.** The
[packet](../../artifacts/generated-results/elliptic-curves/det1092_trace_parity_descent_v1/)
pins only existing generic construction inputs for the theorem, and labels
the38 old control decisions as retrospective evaluation. The checker proves
the universal chord-norm identity, the even-trace norm construction and
the universal halving-quartic discriminant identity. It checks the three
trace parities/heights and the existing quotient and control summaries.
It does not repeat the131072-orbit census or reconstruct any exceptional point.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_trace_parity_descent.sage
```

The proof combines written standard cohomology/local geometry with exact
symbolic and input checks; it is not a formal proof-assistant result.
No point search, new parameter, class group, pilot change or detached job
was started.
