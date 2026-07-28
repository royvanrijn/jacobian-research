# Universal unibranch spectator countermodels

> **Status: exact countermodels to a purely local exclusion.**  Finite
> flatness, normality, smooth integral source, clean boundary ramification,
> a singular unibranch target image, and a separate étale affine sheet are
> jointly consistent in every rank at least four.  The missing Keller input
> is global: the distinguished étale open must be \(\mathbb A^2\), with
> trivial units and the prescribed free boundary class group.

The exact replay is
[`cas/verify_unibranch_spectator_models.py`](cas/verify_unibranch_spectator_models.py).
Its pinned artifact is
[`../artifacts/generated-results/jc2_unibranch_spectator_countermodels.json`](../artifacts/generated-results/jc2_unibranch_spectator_countermodels.json).

## 1. The universal family

For every \(n\ge3\), define

\[
\pi_n:\mathbb A^2_{T,u}\longrightarrow\mathbb A^2_{u,v},
\qquad
v=T^{n+1}-T^n+uT. \tag{1.1}
\]

The source algebra has the monic presentation

\[
k[u,v,T]/(T^{n+1}-T^n+uT-v),
\]

so it is finite free of rank

\[
N=n+1
\]

over \(k[u,v]\), with basis \(1,T,\ldots,T^n\).  The source is the smooth
integral surface \(\mathbb A^2_{T,u}\).

Its Jacobian is

\[
J_n=u-nT^{n-1}+(n+1)T^n. \tag{1.2}
\]

Thus the complete ramification divisor is the smooth curve

\[
E_n=(J_n=0).
\]

On \(E_n\),

\[
u=nT^{n-1}-(n+1)T^n,\qquad
v=(n-1)T^n-nT^{n+1}. \tag{1.3}
\]

At \(T=0\), the two target coordinates have orders \(n-1\) and \(n\).
Since these are coprime, the image is a singular unibranch plane branch of
type \((n-1,n)\).  The generic transverse index of \(E_n\) is two.

## 2. The affine spectator

The fiber over the singular target value \((u,v)=(0,0)\) is

\[
T^{n+1}-T^n=T^n(T-1). \tag{2.1}
\]

It therefore has two points:

\[
\boxed{
k[T]/(T^n)\quad+\quad k\ \text{ at }T=1.
}
\]

The first is the length-\(n\) curvilinear boundary collision.  At the
spectator,

\[
J_n(1,0)=1,
\]

so \(T=1\) is a reduced étale affine sheet over the same singular target
value.  This is precisely the configuration that a purely local
finite-flat theorem was supposed to exclude.

No source singularity, nonnormality, nonflatness, or disconnected generic
cover is involved.

## 3. Euler budget

The ramified component has generic local multiplicity two and special local
multiplicity \(n\).  Its Orevkov cost is

\[
2+(n-2)=n=N-1. \tag{3.1}
\]

Hence the packet exactly saturates the global Euler-multiplicity budget; it
is not excluded.

This explains the role of the spectator.  Without it, the related rank-\(n\)
map

\[
(T,u)\longmapsto(u,T^n+uT)
\]

has the same cost \(n\) but budget \(n-1\), and is excluded.  One étale
spectator supplies exactly the missing unit of global degree.

## 4. The quartic extremal member

For \(n=3\), equation (1.1) is a quartic cover.  At the origin its fiber is

\[
T^3(T-1),
\]

the exact `3+1` cusp packet from
[`JC2_QUARTIC_PACKET_FRONTIER.md`](JC2_QUARTIC_PACKET_FRONTIER.md).

The same cover also has two critical points

\[
T_\pm=\frac{1\pm\sqrt3}{4}
\]

mapping to

\[
(u,v)=\left(\frac18,-\frac1{64}\right).
\]

The complete fiber factors as

\[
\left(T-T_+\right)^2\left(T-T_-\right)^2. \tag{4.1}
\]

Thus one finite polynomial cover realizes both surviving quartic packet
types:

\[
\boxed{3+1\quad\text{and}\quad2+2.}
\]

This is an exact model of why finite-flat fiber arithmetic and the Orevkov
budget do not close the quartic frontier by themselves.

## 5. Exact global failure

The family is not a Keller-normalization model.  Put \(W=J_n\).  Equation
(1.2) has the polynomial inverse

\[
u=W+nT^{n-1}-(n+1)T^n.
\]

Therefore `(T,W)` is a polynomial coordinate system on the source and

\[
\mathbb A^2_{T,u}\setminus E_n
\cong
\operatorname{Spec}k[T,W,W^{-1}]
=\mathbb A^1\times\mathbb G_m. \tag{5.1}
\]

In particular:

- \(W\) is a nonconstant unit on the étale open;
- \(E_n\) is principal;
- \(\operatorname{Cl}(\mathbb A^2)=0\).

An actual Keller normalization instead has distinguished open
\(\mathbb A^2\), unit group \(k^\times\), and freely independent classes of
its missing-boundary primes.  The universal family fails exactly these
global conditions and no local ones.

## 6. Consequence for the direct programme

The proposed local statement

> a finite-flat algebra over a singular unibranch plane branch cannot have
> all ramification in the deleted boundary while retaining an affine étale
> sheet

is false.

The correct theorem must include the global distinguished-open condition.
Class-group freeness alone is not immediately contradictory:
\(\operatorname{div}(g)=2E+A\) merely gives \([A]=-2[E]\).  The remaining
attack must combine:

1. \(U\simeq\mathbb A^2\) and \(U^\times=k^\times\);
2. the free boundary basis in \(\operatorname{Cl}(B)\);
3. connected monodromy of the finite cover; and
4. the extra jump-free collision required to connect the spectator to the
   unibranch block.

The maximal monodromy supported on the \(n\)-fold collision has orbit sizes
\((n,1)\) and fixes the spectator.  A connected degree-\(n+1\) cover needs
another event involving that sheet.  The quartic model shows that a
`2+2` boundary self-collision can provide such an event without exceeding
the Euler budget.

Thus local length, conductor, and Euler inequalities are exhausted.  Any
further direct progress must be global topology or the boundary
localization sequence, not a stronger local Artin-algebra enumeration.

## 7. Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_unibranch_spectator_models.py
```

Intentional artifact regeneration uses `--refresh`.  The replay checks
`n=3,...,10`, while every displayed identity proves the family for all
\(n\ge3\).
