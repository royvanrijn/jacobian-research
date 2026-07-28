# One marked dyadic fiber in two stable Keller classes

The fixed common-fiber pair gives a compact completion of the unfinished
marked-\(\mathbb Q_2\) comparison.  Put

\[
 P(T)=P_{5,1}(T)=T^5+T^3-2T^2+T+1.                    \tag{1}
\]

This single polynomial is an inverse equation for two fixed
determinant-one Keller maps over \(\mathbb Q\).  Their complete fibers are
the same global finite etale scheme, while the ambient maps have different
stable ramified-stratum unit ranks.

## 1. The complete marked \(G_{\mathbb Q_2}\)-set

Modulo \(2\),

\[
 P(T)=(T+1)(T^4+T^3+1),                               \tag{2}
\]

and the quartic factor is irreducible over \(\mathbb F_2\).  The factors are
distinct, so

\[
 \mathbb Q_2[T]/(P)
 \simeq \mathbb Q_2\times U_{2,4},                    \tag{3}
\]

where \(U_{2,4}\) is the unramified quartic extension.  With labels chosen
along geometric Frobenius, its Roe--Turturean marking is

\[
 \boxed{
 \sigma=(1234)(5),\qquad \tau=x_0=x_1=1.
 }                                                     \tag{4}
\]

Both presentation relators hold, the wild normal closure is trivial, and
the exact image is \(C_4\).  Equation (2) identifies (4) with the complete
root action, not merely with its orbit partition.

The dependency-free certificate is
[`gq2_common_quintic_stable_pair.json`](../arithmetic/certificates/gq2_common_quintic_stable_pair.json).

## 2. The connected global fiber

The discriminant is

\[
 \operatorname{Disc}(P)=28085=5\cdot41\cdot137.        \tag{5}
\]

Reduction modulo \(17\) is irreducible.  Hence \(P\) is irreducible over
\(\mathbb Q\), and

\[
 A=\mathbb Q[T]/(P)                                   \tag{6}
\]

is a quintic field.  It has one real root.  Thus the global fiber is
connected even though its base change to \(\mathbb Q_2\) has two connected
components.

## 3. Two fixed determinant-one surroundings

Let

\[
 H(T)=T^5+T^3-2T^2,\qquad G(T)=H(T)+T.                 \tag{7}
\]

The fixed weighted map attached to \(H\) has determinant \(-4\).  At target

\[
 (-1/4,-1,1)
\]

its inverse polynomial is (1).  Scaling its first output by \(-1/4\) gives
a determinant-one map with common-fiber target

\[
 \boxed{q^{\rm wt}=(1/16,-1,1).}                      \tag{8}
\]

The fixed quadratic-gauge map attached to \(G\) has determinant \(-2\).
After the standard output scaling to determinant one, its target is

\[
 \boxed{q^{\rm quad}=(1,0,-2),}                        \tag{9}
\]

and its inverse polynomial is again (1).  The reconstruction theorems give
complete scheme-theoretic fibers

\[
 (F^{\rm wt})^{-1}(q^{\rm wt})
 \simeq \operatorname{Spec}A
 \simeq (F^{\rm quad})^{-1}(q^{\rm quad}).             \tag{10}
\]

The two canonical ramified target strata have normalized forms

\[
 \mathbb A^1\times\mathbb G_m
 \quad\text{and}\quad
 \mathbb G_m^2,                                        \tag{11}
\]

so their unit ranks are respectively \(1\) and \(2\).  These ranks survive
identity stabilization.  Therefore the two maps are not stably
polynomially left--right equivalent.

> **Marked dyadic stable-separation theorem.**  The complete marked
> \(G_{\mathbb Q_2}\)-set (4) occurs as the base change of one connected
> complete global fiber in two explicit determinant-one Keller maps whose
> stable ramified-stratum unit ranks are \(1\) and \(2\).  Consequently, the
> complete local Galois action of a Keller fiber does not determine the
> stable polynomial left--right class of its ambient map.

This conclusion is stronger than comparing two locally isomorphic fibers:
both maps in (10) contain the identical global algebra \(A\).

## 4. The fixed local ball

The discriminant in (5) is a \(2\)-adic unit, so the automatic local
stability exponent is

\[
 2v_2(\operatorname{Disc}P)+1=1.                      \tag{12}
\]

Thus every \(2\)-integral parameter \(u\equiv1\pmod2\) in the fixed pencil
\(P_{5,u}\) has the same algebra (3) and the same marked action (4), up to
relabeling.  The locally prescribed fixed-pair theorem, with an auxiliary
inert-prime condition or Hilbert irreducibility, therefore gives infinitely
many connected common fibers carrying this marked action.  Equation (1) is
the smallest displayed member, not a parameter-height minimality claim.

## Verification

Run

```bash
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_common_quintic_stable_pair.json --json
.venv/bin/python scripts/verify_marked_q2_stable_separation.py
```

The first command checks the marked presentation and the exact unramified
polynomial comparison without third-party dependencies.  The second checks
global irreducibility, discriminant, both inverse equations and targets,
both determinant-one Jacobians, and the stable-unit-rank ledger.
