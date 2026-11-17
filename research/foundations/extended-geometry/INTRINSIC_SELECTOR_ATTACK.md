# Intrinsic-selector attack on the symmetric quintic

This note attacks the warning in stable normalization functoriality directly:
a root, completion point, or boundary order named by a presentation is not
automatically preserved by stable polynomial left--right equivalence.  The
smallest weighted test requested here is degree five, with one identity
variable allowed and with stable source transformations permitted to mix that
variable with a boundary coordinate.

The search finds the expected automorphism of the **coarse** normalized
incidence cover, but no counterexample to the full map-intrinsic
Zariski--Main package.  The failed realization supplies the missing reusable
uniqueness lemma: affine membership of a divisorial branch is invariant under
arbitrary polynomial stabilization and source automorphisms.  On the
boundary-clean weighted locus this singles out the distinguished root sheet
before any presentation coordinate is named.

## 1. The smallest symmetric seed

Take

\[
 H(W)=\frac{W^2-W^5}{3}
     =\frac{W^2(1-W^3)}3.                              \tag{1.1}
\]

It is a normalized exact-double weighted seed:

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1,
\]

and

\[
 H''(W)=\frac{2-20W^3}{3},\qquad
 \kappa=H''(1)=-6,\qquad a_0=-\frac54.                \tag{1.2}
\]

Let \(\zeta\) be a primitive cube root of unity.  Then

\[
 H(\zeta W)=\zeta^2H(W),\qquad H''(\zeta W)=H''(W).
                                                                    \tag{1.3}
\]

The full affine stabilizer of the effective Hessian divisor is exactly
\(\mu _3\).  Indeed, comparison in

\[
 \alpha^2H''(\alpha W+\beta)=\gamma H''(W)
\]

forces \(\beta=0\) and
\(\gamma=\alpha^2=\alpha^5\), hence \(\alpha^3=1\).
Thus this seed is the smallest exact-double degree-five test with a
nontrivial affine symmetry preserving both zero and infinity.

## 2. The coarse automorphism really permutes candidate marks

For

\[
 E(W;s,t)=H(W)-sW+t,
\]

equation (1.3) gives

\[
 E(\zeta W;\zeta s,\zeta^2t)=\zeta^2E(W;s,t).          \tag{2.1}
\]

Hence the normalized two-parameter incidence cover has a pencil-compatible
self-equivalence

\[
 (W,s,t)\longmapsto(\zeta W,\zeta s,\zeta^2t).         \tag{2.2}
\]

It preserves the effective Fitting divisor, its multiplicities, the node and
conductor constructions, zero, and the unique point at infinity.  In weighted
target coordinates it is compatible with

\[
 (A,B,C)\longmapsto(\zeta^2A,\zeta B,C),               \tag{2.3}
\]

so both reduced target boundary images are fixed.  But it cyclically permutes
the three nonzero primitive roots

\[
 \{1,\zeta,\zeta^2\}.                                  \tag{2.4}
\]

This is the sought automorphism of the unmarked finite-cover data.  It proves
that no rule using only the Fitting divisor, reduced target images, zero, and
infinity can select the root \(1\).  On that coarse object the candidates
must be retained as the finite \(\mu _3\)-set, equivalently as stacky or
unordered decorated data.

It is not a target-fixed deck transformation: the universal two-parameter
cover has none.  More importantly, (2.2) does not preserve the Zariski--Main
open immersion of the weighted polynomial map.

## 3. The regular-sheet calculation

Put \(c=-H'(1)=1\).  On the sheet of a nonzero primitive root
\(\rho^3=1\) over \(C=0\), the reconstruction equations give

\[
 \gamma_0=-H'(\rho)/c=\rho,\qquad
 u_0=W/\gamma=1,\qquad x=C/\gamma.                    \tag{3.1}
\]

The numerator of the remaining source coordinate is

\[
 \gamma-1-a_0(u-1).
\]

Its constant term on the \(\rho\)-sheet is

\[
 \rho-1.                                               \tag{3.2}
\]

For \(\rho=\zeta,\zeta^2\), this is nonzero.  Since
\(\operatorname{ord}_C(x)=1\), equation (3.2) gives

\[
 \operatorname{ord}_C(z)=-2.                          \tag{3.3}
\]

Those two sheets are boundary sheets.  For \(\rho=1\), write
\(W=1-BC+O(C^2)\).  Then

\[
 \gamma=1-5BC+O(C^2),\qquad
 u=1+4BC+O(C^2),
\]

and \(a_0=-5/4\) cancels the constant and linear terms.  Thus \(z\) is regular.
The root-one sheet is the unique degree-one unramified root sheet meeting the
affine source.

The coarse symmetry (2.2) sends this affine sheet to a pole-two boundary
sheet.  Therefore it is not an automorphism of

\[
 \mathcal B(F)=
 (\overline X_F\to Y,\ j_F:\mathbb A^3\hookrightarrow\overline X_F,
 \partial_F).
\]

This is exactly where an attempted realization by polynomial left--right
transformations fails.

## 4. Stable mixing cannot repair the failure

The preceding obstruction is independent of coordinates.

**Stable affine-membership lemma.**  Let \(v\) be a divisorial valuation of
\(k(x_1,\ldots,x_n)\).  Membership of its center in the affine source
\(\mathbb A^n\) is the condition

\[
 v(f)\geq0\quad\text{for every }f\in k[x_1,\ldots,x_n],
\]

equivalently \(v(x_i)\geq0\) for all coordinate generators.  If
\(\Phi\) is a polynomial automorphism, both \(\Phi\) and \(\Phi^{-1}\) express
each coordinate polynomially in the other set.  Hence this nonnegativity
condition is preserved in both directions.

After adjoining identity variables \(T_1,\ldots,T_q\), the boundary valuation
has its Gauss extension \(v(T_i)=0\).  The same argument applies to every
automorphism of
\(k[x_1,\ldots,x_n,T_1,\ldots,T_q]\).  Consequently no stable polynomial
source transformation can turn a boundary branch into an affine branch.
Mixing an identity variable with a pole merely transfers the negative
valuation to another coordinate; the inverse automorphism prevents all poles
from disappearing simultaneously.

For the requested one-variable test, the extra-root valuation contains

\[
 (v(x),v(y),v(z),v(T))=(1,0,-2,0).                    \tag{4.1}
\]

The elementary attacks \(T\mapsto T+z\), \(z\mapsto z+T\), coordinate swaps,
and \(z\mapsto z+P(x,y,T)\) all visibly retain a negative entry.  The lemma
covers wild as well as tame polynomial automorphisms, so bounded enumeration
of stable shears is unnecessary.

## 5. Downstream selector audit

| Claimed datum | Coarse \(\mu _3\) action | Intrinsic selector in the full package | Verdict |
|---|---|---|---|
| unordered Fitting/Hessian divisor | preserved, support permuted | full Fitting ideal | intrinsic only as unordered effective data |
| distinguished primitive root | \(1\mapsto\zeta\) | unique degree-one unramified sheet meeting \(j_F(\mathbb A^3)\) | intrinsic on the boundary-clean locus |
| zero-cluster center | fixed | unique center where the discriminant-normalization branch meets the second boundary | intrinsic under the exact-double uniqueness hypothesis |
| infinity | fixed by every affine scaling | unique completion point after affine-coordinate rigidity from at least two Fitting supports | intrinsic for the degree-five test |
| target boundary order | both components fixed | \(Z_\Delta\) uniquely receives a ramified boundary prime; \(Z_0\) does not | intrinsic |

Thus the attack produces a genuine warning example for the **coarse**
selector assumption but no pair of stably polynomially left--right-equivalent
weighted maps violating the current boundary-clean theorems.  Any statement
made before adjoining the regular-reconstruction incidence must weaken the
root choice to an unordered finite cover or quotient stack.  Statements that
use the full package may retain a point label only after proving the selector
listed in the third column.

## 6. Reusable uniqueness checklist

For each family, a presentation mark is stable only after the following
family-specific facts are established.

1. **Candidate completeness:** enumerate every generic branch over the
   relevant intrinsic target stratum by the degree sum.
2. **Stable branch status:** characterize candidates by membership in the
   Zariski--Main affine open; the stable affine-membership lemma makes this
   invariant under arbitrary identity-variable mixing.
3. **Unique signature:** prove exactly one candidate has the required tuple
   of affine/boundary status, ramification degree, residue degree, and
   incidence centers.
4. **Completion rigidity:** retain a compactification point only when the
   induced normalization-coordinate automorphism is forced to extend and fix
   it.
5. **Ordering rigidity:** order boundary images only by a unique intrinsic
   signature, such as the presence of ramification.

Failure of item 3 means the correct output is unordered or stacky decorated
data, not a selected point.

## 7. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_intrinsic_selector_attack.py
```

The checker enumerates the affine Hessian stabilizer, verifies the exact
pencil and target transformations, confirms the cyclic permutation of the
three candidate roots, computes the regular-versus-pole selector, and runs
representative one-variable stable mixing attacks.
