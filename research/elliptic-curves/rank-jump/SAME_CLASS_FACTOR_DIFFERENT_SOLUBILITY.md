# The same elementary class factor can be rational or Sha

The retained fixed-cubic controls give a direct counterexample to a
class-field-only explanation of simultaneous rational solubility.

For each of the six deformations there is a nonzero block of anchor-rational
classes that is **strictly split at every retained bad place**, yet every
nonzero class in that block is CT-obstructed on the deformation. At \(u=-1\),
a three-dimensional elementary \(S\)-class direct factor has the same seven
nonzero dual characters on both curves: all seven are rational on the
anchor, and all seven represent nonzero Sha classes on the deformation.

The cubic field, specified \(S\), squareclasses, half ideals and Artin
evaluations are identical in this comparison. Their existence does not
determine whether the elliptic covers have rational points.

This advances the **solubility** analysis. It does not weaken the earlier
incidence certificates, identify the full Sha group, or establish a
prospective rank predictor.

## Fix the arithmetic data and change the covering map

Use the already retained Fermigier–Mestre anchor
\[
E_0:y^2=x^3+Ax+B,\qquad K=\mathbf Q(\theta),\quad\theta^3+A\theta+B=0,
\]
with
\[
\begin{aligned}
A&=-5750886029903523759416717668139307,\\
B&=167347710468055045100164888198438918505621536951206.
\end{aligned}
\]
Its twenty certified independent rational Kummer classes span \(W\).
They are anchor classes, not generic sections: this fixed-cubic pencil
has arithmetic generic rank zero.

The same labelled two-torsion algebra for \(E_u\) is obtained through
\(\alpha_u=\theta+u\theta^2\), where
\[
E_u:\quad y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]
Let \(S_u\) contain 2, infinity and the complete retained bad support of
both \(E_0\) and \(E_u\). Define
\[
V_u^0=\{\beta\in W:\beta\text{ is square in }K_{\mathfrak p}
       \text{ at every place above }S_u\}.
\]
This imposes zero localization, not just membership in both elliptic
local point images.

With this **same specified set \(S_u\)** on both curves,
\[
U_{S_u}(E_0)=U_{S_u}(E_u)
 \simeq\operatorname{Hom}
 \bigl(\operatorname{Cl}(\mathcal O_{K,S_{u,K}}),\mathbf F_2\bigr).
\]
The equality follows from the earlier
[strict Selmer identification](STRICT_SELMER_AND_ARTIN_BLOCKS.md).
The minimal bad-prime sets of the curves need not agree; this experiment
uses their union. It does not assert that every possible curve-specific
incidence feature is unchanged.

Every element of \(V_u^0\) is rational on \(E_0\), since \(V_u^0\subset W\).
The issue is its rationality after the same squareclass is viewed as an
\(E_u[2]\)-class.

## Strict blocks with injective solubility obstruction

Let \(W_u\subset W\) be the retained locally admissible space for \(E_u\).
It contains \(V_u^0\). Compute the rectangular cross-pairing
\[
T_u:V_u^0\longrightarrow W_u^\vee,\qquad
\beta\longmapsto
 \bigl(\gamma\mapsto\operatorname{CT}_{E_u}(\beta,\gamma)\bigr).
\]
All six cross-maps are injective.

| \(u\) | \(\dim W_u\) | \(\dim V_u^0\) | \(\operatorname{rank}T_u\) | Anchor-rational part of \(V_u^0\) | Deformation-rational part of \(V_u^0\) |
|---:|---:|---:|---:|---:|---:|
| \(-3\) | 17 | 2 | 2 | 2 | 0 |
| \(-2\) | 13 | 1 | 1 | 1 | 0 |
| \(-1\) | 18 | 5 | 5 | 5 | 0 |
| \(1\) | 13 | 3 | 3 | 3 | 0 |
| \(2\) | 13 | 2 | 2 | 2 | 0 |
| \(3\) | 15 | 4 | 4 | 4 | 0 |

These last zeros are exact statements **inside the displayed spaces**.
They are not full-curve rank bounds.

Rational Kummer classes pair trivially with every Selmer class. Therefore
an injective cross-map proves
\[
\boxed{V_u^0\cap\delta E_u(\mathbf Q)=0,\qquad
       V_u^0\hookrightarrow\Sha(E_u)[2].}
\]
This needs neither a complete Selmer basis nor a conjecture about the
finiteness of Sha. The standard Selmer exact sequence and CT interpretation
are recalled in [Fisher, *On binary quartics and the Cassels–Tate pairing*, Section 1](https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf).

Using only the self-pairing of \(V_u^0\) would lose this conclusion.
For example, at \(u=-2\) the strict space has dimension one, so its
alternating self-pairing is necessarily zero, while its cross-map has
rank one. Its single nonzero class is proved nonrational.

The strict space at the union of **all six** bad-place sets is zero.
That primary negative endpoint is retained. Among the fifteen pairs of
nonzero deformations, nine have a nonzero common strict space; each such
space is obstructed on both members. The solubility switch is proved by
comparison with the anchor, where those same classes are rational.

## An elementary direct factor does not repair the implication

The preceding actual-fibre experiment found
[large soluble elementary \(S\)-class factors](SOLUBLE_ELEMENTARY_S_CLASS_BLOCKS.md).
To test whether that stronger incidence structure could itself force
rationality, a second protocol fixes the largest strict block above:
\(u=-1\), with five anchor-coordinate masks
\[
17108,\quad34628,\quad65575,\quad404296,\quad528076.
\]
The choice is retrospective and made before its Artin matrix was evaluated.

For each corresponding class \(\beta_j\), exact ideal arithmetic constructs
\[
(\beta_j)=\mathfrak J_j^2.
\]
One reduction per ideal, with the principal multiplier retained, gives
the Artin matrix on the five half-ideal images in the common \(S\)-class
group:
\[
M_{ij}=\chi_{\beta_i}([\mathfrak J_j]),\qquad
M=
\begin{pmatrix}
1&1&0&1&0\\
0&0&0&0&1\\
1&1&1&1&0\\
1&1&0&1&0\\
1&1&0&1&1
\end{pmatrix}.
\]
It has rank three.

Taking character indices \(0,1,2\) and ideal indices \(0,2,4\), all
zero-based, gives the invertible submatrix
\[
\begin{pmatrix}1&0&0\\0&0&1\\1&1&0\end{pmatrix}.
\]
Thus
\[
C:=\operatorname{Cl}(\mathcal O_{K,S_{-1,K}})
 =H\oplus C',\qquad H\simeq(\mathbf Z/2)^3.
\]
The proof is the same exact retraction as in the previous factor
certificate: the selected characters restrict isomorphically on \(H\),
and the selected ideals have order dividing two. Dual ideal words for
the three selected characters are
\[
\mathfrak J_0\mathfrak J_2,\qquad
\mathfrak J_4,\qquad
\mathfrak J_2.
\]

The selected character space
\[
V'=\langle\beta_{17108},\beta_{34628},\beta_{65575}\rangle
 \simeq H^\vee
\]
is rational on \(E_0\) and injects into \(\Sha(E_{-1})[2]\).
For a compact CT witness, pair its basis against inherited anchor masks
\(1,6,128\). The matrix is
\[
\begin{pmatrix}
0&1&0\\
1&0&1\\
0&1&1
\end{pmatrix},
\]
also invertible. Hence all seven nonzero combinations in \(V'\) are
obstructed, not merely the three listed generators.

This is the concrete failed implication:
\[
\boxed{\text{an elementary strict \(S\)-class factor with its dual classes}
\quad\not\Longrightarrow\quad
\text{a simultaneously rational elliptic block}.}
\]
The factor and dual classes are exactly the same on the two curves.
Their elliptic cover maps, and the map from strict Selmer classes to
Sha, differ. No Artin evaluation changes when passing from one curve
to the other.

## What this changes in the mechanism ranking

1. **Solubility, strongest remaining mechanism:** the labelled
   4-torsion extension together with its local Kummer conditions.
   The [Jacobian comparison](JACOBIAN_LOCAL_CONDITIONS_AND_CT.md)
   identifies the CT difference with the isogeny-descent pairing on
   the glued genus-two Jacobian. This structure can distinguish the
   two rationality maps even on classes with zero localization.
2. **Incidence, still necessary and useful:** strict \(S\)-class characters
   and their elementary ideal factors. The high-fibre results remain
   valid, but their rationality requires additional information.
   This experiment rules out sufficiency, not possible statistical
   usefulness after controlling the generic subgroup.
3. **Disproved sufficient explanation:** fixed cubic field, norm-square
   representatives, zero bad-place localizations, half-ideal relations,
   and a nondegenerate Artin factor, even taken together. These data
   coincide on the displayed three-dimensional rational/Sha pair.
4. **Weak substitutes:** the abstract degree or group type of a shared
   descent field. The six retained torsion-difference fields already
   have the same degree-48 module type while their CT forms differ.
   The actual labelled extension and local lifting conditions matter.
5. **Visibility:** irrelevant to the nonzero classes excluded here.
   Their covers have no rational point; a larger point-search budget
   cannot expose one.

For Agent 1, this is a restriction on any future selector: an arithmetic
incidence score must be accompanied by a distinct solubility assessment.
None of these oracle-derived masks or ideal factors may be used to select
prospective candidates.

The next bounded mathematical target is the displayed \(3\times3\) CT
matrix. Derive it through the glued Jacobian's torsion-difference extension
and local lifting data, rather than only reading the retained elliptic
CT certificate. Its inputs are fixed: \(u=-1\), three strict masks and
three partner masks. Matching those nine bits would connect an
equation-defined higher-descent object to the simultaneous obstruction.
It would still not make CT-vanishing sufficient for a rational point.

## Scope, evidence and replay

The [strict-space protocol](STRICT_DEFORMATION_SOLUBILITY_PROTOCOL.json)
uses only the six existing parameters, twenty anchor classes and retained
CT matrices. The primary and secondary endpoints were frozen before the
zero-localization kernels were computed. It performs no new CT calculation.

The [strict-space report](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_deformation_solubility_v1.json)
retains every local constraint, all six kernels, all fifteen paired
comparisons, complete cross-pairing rows and explicit obstruction witnesses.
Independent replay verifies:

- complete support using the retained prime factorizations and
  \(\operatorname{disc}(f_u)=D(u)^2\operatorname{disc}(f)\);
- the strict kernels by a separate binary linear-algebra implementation;
- 537 finite local-square checks and exact positivity at the real places;
- injectivity of the cross-maps from the pinned CT matrices.

The CT values themselves are inherited from the repository's exact
certificates. This replay independently checks the new restrictions and
ranks, not a fresh computation of every underlying Hilbert symbol.

The [five-class Artin protocol](STRICT_SHA_ARTIN_PROTOCOL.json) permits
one reduction per ideal, 25 evaluations and a 30-second worker cap.
All 25 residues were units; no exceptional local branch or factorization
was needed. The
[half-ideal and residue inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_Sha_Artin_inputs_v1.json)
and [factor report](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_Sha_Artin_v1.json)
retain the exact evidence.

A separate verifier proves all five ideal-square identities by lattice
multiplication, checks their principal transports, verifies the cyclic
residue rings and recomputes every Jacobi symbol independently.
Five targeted tests also guard the full cross-pairing requirement,
the zero common endpoint, the elementary factor and all seven nonzero
CT-obstructed combinations.

```sh
python3 elliptic-curves/rank-jump/strict_deformation_solubility.py check
python3 elliptic-curves/rank-jump/strict_Sha_Artin.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_strict_deformation_solubility.py
sage -python elliptic-curves/rank-jump/verify_strict_deformation_solubility.py --all
sage -python elliptic-curves/rank-jump/verify_strict_Sha_Artin.py
```

The local replay launches six sequential workers with a 30-second cap each.
All completed. Arithmetic used Sage 10.9 and PARI 2.17.3.
No parameter sweep, point search, full class group, active-search file or
mathematical-status entry was changed.
