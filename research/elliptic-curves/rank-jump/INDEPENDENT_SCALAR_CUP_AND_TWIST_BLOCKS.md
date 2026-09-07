# Independent scalar cup products and simultaneous twist obstructions

Follow-up: [the nonscalar norm-lifting control](NONSCALAR_CUP_BLOCK_AND_SELF_GLUING.md)
independently computes the opposite maps for `1+theta` and `-(1+theta)`,
constructs two unramified norm lifts, and identifies the zero case's
elliptic self-gluing. It does not close the production nine-bit target.

The half-ideal Artin matrix contains an independently computable part of
the higher solubility obstruction. For strict classes, its failure to be
symmetric is exactly the scalar \(-1\) cup product:
\[
\boxed{\chi_\psi\bigl((-1)\cup\beta\bigr)
 =\chi_\psi([\mathfrak J_\beta])
  +\chi_\beta([\mathfrak J_\psi]).}
\]
Here \((\beta)=\mathfrak J_\beta^2\) outside the allowed set \(S\), and
values are written in \(\mathbf F_2\).

This yields new CT values for the \(-1\) quadratic twists of three retained
curves. No point search, cubic-field norm search at production scale, or
new elliptic CT backend calculation is needed.

| Original fibre | Known independent subgroup | Generic subgroup / observed quotient | Retained strict rational dimension \(k\) | CT rank on that space after twisting | Necessary rational subspace after twisting, dimension at most |
|---|---:|---:|---:|---:|---:|
| A1/MW16-05, \(307/206\) | 25 | \(16/9\) | 10 | 8 | 2 |
| Published R17, \(-2300/843\) | 24 | \(17/7\) | 8 | 6 | 2 |
| Published R17, \(-1561/3133\) | 17 | \(17/0\), censored | 6 | 6 | 0 |

The table's last column bounds only the transported **retained strict
space**. It is not a bound on the total rank of the twist. The generic
rational subgroup of the original family is not preserved by this twist,
and the last row's observed zero quotient is not a proved exact rank.

The two high-gain fibres and the low-gain fibre all show simultaneous
obstruction blocks. Their presence alone does not discriminate a high
jump from an ordinary specialization. The comparison instead isolates a
solubility switch while keeping the cubic algebra and strict incidence fixed.

## Precisely which cup product is computed

Let \(K\) be the cubic descent field, \(S\) contain 2, infinity and the bad
primes, and \(C=\operatorname{Cl}(\mathcal O_{K,S_K})\). Write \(U\) for the
strict class space: its elements are unramified outside \(S_K\) and
locally square at \(S_K\). Then \(U=\operatorname{Hom}(C,\mathbf F_2)\).

For \(\beta\in U\), the field Brauer class \((-1)\cup\beta\) vanishes:
\(\beta\) is square above 2 and at the other allowed primes, positive at
real places, and unramified at the remaining odd primes. Its restricted
ramification cup product can nevertheless survive in \(C/2\).

This scalar \(-1\) must not be confused with the deformation parameter
\(u=-1\). The earlier fixed-cubic obstruction uses
\(\gamma_{-1}=1+\theta\). The two cup matrices are different. For example,
on the retained five-character basis, the scalar matrix has entry
\((0,1)=1\), whereas the retained \(1+\theta\) matrix has entry \((0,1)=0\).

## Proof through the two Bockstein maps

Let \(B_{\mathrm{const}}\) be the connecting map for
\[
0\to\mathbf Z/2\to\mathbf Z/4\to\mathbf Z/2\to0
\]
and \(B_\mu\) the connecting map for
\[
1\to\mu_2\to\mu_4\to\mu_2\to1.
\]
These middle modules have different Galois actions when \(i\notin K\).
A lift of a binary cocycle to integers modulo four shows
\[
B_\mu(\beta)=B_{\mathrm{const}}(\beta)+(-1)\cup\beta.
\]
The Kummer sequences identify
\(B_\mu(\beta)=[\mathfrak J_\beta]\in C/2\).
This comparison is the \(p=2\) case of
[McCallum–Sharifi, Lemma 8.1](https://arxiv.org/pdf/math/0202161).

For completeness, the duality step needed here includes the strict local
conditions. A strict character \(\psi\) has a compact-support lift after
choosing its local trivializations at \(S_K\). The dual of the constant
\(\mathbf Z/4\) sequence is the \(\mu_4\) sequence. Naturality of
Poitou–Tate duality therefore gives
\[
\langle B_{\mathrm{const}}(\beta),\psi\rangle
 =\langle\beta,B_{\mu,c}(\psi)\rangle
 =\chi_\beta([\mathfrak J_\psi]).
\]
The second equality follows from the same Kummer diagram with compact
support: forgetting its boundary trivializations leaves the half ideal.
Changes in those trivializations contribute only local terms at \(S_K\),
which pair to zero because \(\beta\) is strict. Real-place terms also
vanish because both classes are locally trivial there.

The adjoint Bockstein and half-ideal description is the construction in
[Chung et al., Proposition 4.4 and Lemma 4.5](https://academic.oup.com/imrn/article/2019/18/5674/4656167).
The strict compact-support argument above supplies the boundary
justification for the present \(S\)-integer setting. We do not invoke
that paper's symmetry theorem under an unverified \(\mu_4\subset K\)
hypothesis.

Pairing the Bockstein comparison with \(\psi\) proves the boxed identity.
If
\[
A_{ij}=\chi_{\beta_i}([\mathfrak J_{\beta_j}]),
\]
the matrix is consequently
\[
\boxed{M_{-1}=A+A^{\mathsf T}.}
\]
This is alternating, as required, including when \(A\) itself is neither
symmetric nor alternating.

## Why this is an elliptic twist comparison

Let \(E^{(-1)}\) be the quadratic twist of \(E\), with the standard
identification of their two-torsion. Its four-torsion action differs
from that of \(E\) by the scalar quadratic character of \(-1\). In the
Baer difference of the two four-torsion extensions, the extension
cocycle is therefore \(\chi_{-1}I\).

The decorated CT comparison used in the
[earlier torsion analysis](TORSION_DIFFERENCE_AND_CT.md) applies here as
well. On two strict classes its local correction terms at \(S\) vanish;
outside \(S\) the local conditions are unramified. Thus
\[
\operatorname{CT}_{E}(\beta,\psi)
+\operatorname{CT}_{E^{(-1)}}(\beta,\psi)
=\chi_\psi\bigl((-1)\cup\beta\bigr).
\]
The norm-square projection introduces no extra scalar term, exactly as
in the [strict cup argument](CUP_IDEAL_AND_STRICT_LIFTING_OBSTRUCTION.md).

Twisting by \(-1\) adds no bad prime outside this \(S\). The entire strict
space \(U\) is consequently the same for both curves. All the retained
classes used in the production table are rational on \(E\), so the first
CT term is zero. The displayed \(A+A^{\mathsf T}\) is the actual
restriction of the twist's CT form, rather than merely a heuristic score.

A restricted CT rank \(r\) forces at least \(r\) independent images in
\(\Sha(E^{(-1)})[2]\) from that retained space. Its rational subspace has
dimension at most \(k-r\). Conversely, a radical vector need not be
rational, and additional classes may detect it.

## Independent positive norm control

Before using the formula for twist claims, the frozen experiment tested
it in the single small cubic field
\[
K_c=\mathbf Q(t),\qquad t^3-11t^2-14t-1=0,\qquad
S_c=\{2,163,\infty\}.
\]
This is an arithmetic identity control, not a curve candidate.
The bounded worker certified its class group \((\mathbf Z/2)^2\),
checked 31 squareclasses, and found two strict generators
\[
\beta_0=t^2-10t+1,\qquad
\beta_1=t^2-13t+12.
\]
Their half-ideal Artin matrix is
\[
A=\begin{pmatrix}1&0\\1&1\end{pmatrix}.
\]
Two independently solved norm equations in \(K_c(i)\) then produced
\[
a_j^2+b_j^2=\beta_j.
\]
For example,
\[
a_0=\frac{-139+261t-22t^2}{125},\qquad
b_0=\frac{48-2t+4t^2}{125}.
\]
The second witness and all ideal data are retained in the
[control certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_scalar_cup_control_v1.json).

The parity ideals of \(a_j+ib_j\), evaluated by the two unramified
characters, gave
\[
\begin{pmatrix}0&1\\1&0\end{pmatrix}=A+A^{\mathsf T}.
\]
This is a nonzero control, so a routine returning zero for every cup
product would fail it.

A separate verifier uses the retained rational coefficients, exact
polynomial reduction, and Hensel lifts of all three cubic roots and the
two roots of \(-1\) modulo \(5^{64}\). It recomputes the valuation parities
and Artin values without a class-group calculation, a norm solver, or
PARI ideal valuations. The coordinates are integral away from 5 and each
\(\beta_j\) has norm 625, proving that no other finite prime can have
been omitted from the parity ideal.

The global norm succeeds while the restricted cup remains nonzero.
That is the precise distinction the control was designed to test.

## The three retrospective solubility pairs

1. **A1/MW16-05 at \(307/206\), compared with its \(-1\) twist.**
   All ten retained strict classes are rational on the original fibre.
   The twist form has four nondegenerate alternating planes, leaving
   only a two-dimensional necessary rational subspace. The original
   nine quotient directions belong to an independent subgroup; the
   statement does not assign nine generic directions to the twist.

2. **Published R17 at \(-2300/843\), compared with its \(-1\) twist.**
   Eight retained strict classes are rational on the original fibre,
   whose certified subgroup gains seven directions above the generic
   17. Three alternating planes become obstructed on the twist.
   The two-dimensional remaining radical is not a point certificate.

3. **Published R17 at \(-1561/3133\), compared with its \(-1\) twist.**
   The observed subgroup has no recovered quotient gain. Nevertheless,
   all six of its strict generic combinations inject into the twist's
   Sha group: the restricted CT matrix is nonsingular. Strong
   simultaneous solubility changes therefore also occur in the
   low-gain control. They are not sufficient explanations of a large jump.

These pairs control the cubic algebra, coefficient scale up to signs,
and strict class-field incidence exactly. They do not control visibility,
which is not being measured.

Explicit alternating planes, in the retained strict-basis bitmask
coordinates, are:

| Fibre | Pairs with CT value one, mutually orthogonal | Radical basis |
|---|---|---|
| A1/MW16-05 \(307/206\) | \((1,2),(6,11),(30,132),(37,94)\) | \(392,539\) |
| R17 \(-2300/843\) | \((1,2),(7,8),(45,139)\) | \(18,76\) |
| R17 \(-1561/3133\) | \((1,2),(5,17),(14,34)\) | empty |

These are decompositions of a bilinear obstruction. They are not claimed
to be canonical Galois submodules or common rational constructions.

## Mechanisms, rejected shortcuts, and the next gate

Ranked by current evidential strength:

1. **Solubility:** labelled four-torsion extensions and their restricted
   cup products can obstruct many classes together. The scalar part is
   now computable independently from half-ideal Artin data, with a
   nontrivial independent norm control.
2. **Incidence:** strict class-field blocks survive this twist unchanged.
   They supply candidate directions but do not decide which orientation
   of the twist pair is rational.
3. **Weak explanation:** a large abstract class factor, or a large
   obstruction rank in isolation, does not predict the original jump.
   The low-gain control has a nonsingular six-dimensional obstruction.
4. **Missing arithmetic:** the non-scalar \(1+\theta\) cup values and the
   local corrections needed for the original nine CT bits still require
   independent evaluation.
5. **Missing positive implication:** even a vanishing full lifting
   obstruction must be followed by a rational-point construction or
   another exact solubility certificate.

For Agent 1, a half-ideal Artin matrix built from point-independent strict
Selmer generators would provide this relative solubility information
before exceptional points are supplied. Without a known soluble
reference, however, the difference of two CT forms does not identify
which curve should gain rank. None of the retrospective point masks used
here is a prospective selector input. No search policy was changed.

The next falsifiable gate is to recover at least one non-scalar cup
evaluation independently and distinguish it from the scalar formula.
The retained five-character pair \((0,1)\) is particularly diagnostic:
the two predicted bits differ. The existing bounded cubic-field norm
solvers have no demonstrated production-scale solution for that task;
a small-field control success is not evidence that they will scale.

## Replay and scope

The [protocol](SCALAR_CUP_PROTOCOL.json) freezes the four retained Artin
matrices, one small cubic control and a 30-second worker limit.
The [analysis certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_scalar_cup_v1.json)
contains the three production matrices, radical coordinates, obstruction
ranks and the separate scalar matrix on the fixed-cubic five-class space.

```sh
python3 elliptic-curves/rank-jump/scalar_cup.py check
sage -python elliptic-curves/rank-jump/verify_scalar_cup_control.py
sage -python -m unittest discover -s elliptic-curves/rank-jump -p test_scalar_cup.py
```

The small-field norm construction completed within its 30-second limit.
Its explicit witnesses are permanent replay inputs; replay does not
solve them again. All computations are retrospective or arithmetic
controls. No active-search file, candidate population, worker setting,
or mathematical-status entry was changed.
