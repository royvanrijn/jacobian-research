# Isotropic harmonic-channel obstruction below \(\Delta^6\)

## 1. Result and scope

Let \(V\) be a three-dimensional quadratic space over a characteristic-zero
field, let \(\rho\) be its quadratic form, and let \(\Delta\) be the dual
Laplacian.  After scalar extension, choose isotropic covectors
\(v_i\in V^*\), and write

\[
 L_i(z)=\langle v_i,z\rangle,
 \qquad \Delta L_i^{\ell_i}=0.
\tag{1.1}
\]

For an even balanced degree \(d=2k\), consider two or three distinct positive
even harmonic degrees \(\ell_i\le d\) and

\[
 \boxed{
 P=\sum_{i=1}^r a_i\rho^{(d-\ell_i)/2}L_i^{\ell_i},
 \qquad r\in\{2,3\}.
 }
\tag{1.2}
\]

Thus (1.2) uses one coherent state in each selected irreducible summand of

\[
 \mathbb C[V]_d
 =\bigoplus_{j=0}^{d/2}\rho^j\mathcal H_{d-2j}.
\tag{1.3}
\]

Under \(SO_3\simeq PGL_2\), the coherent states are the pure-power points of
the binary-form model \(\mathcal H_\ell\simeq\operatorname{Sym}^{2\ell}\).
This is therefore an equivariant multiple-harmonic-channel family, not a
one-profile cusp lift.  The scalar \(\mathcal H_0\) channel is omitted without
loss: its coefficient is the first Reynolds moment, so the pure-zero premise
sets it to zero immediately.

> **Theorem 1.1 — coherent harmonic-channel obstruction.**  For
> \(d=4,6,8,10\), the family (1.2) contains no homogeneous GVC witness for
> \(\Delta^{d/2}\).  More precisely, if its first nine pure Reynolds moments
> vanish, then all active \(v_i\) are proportional.  That remaining stratum is
> one-sided for a phase torus and satisfies the GVC conclusion for every fixed
> multiplier.

For degree four the third moment suffices; for degrees six, eight, and ten the
cutoffs supplied by the calculation are respectively five, seven, and nine.
The theorem treats every collision stratum of the isotropic directions and
every coefficient boundary.

This is not a classification of arbitrary vectors in
\(\mathcal H_{\ell_1}\oplus\mathcal H_{\ell_2}\oplus\mathcal H_{\ell_3}\),
of two profiles in the same irreducible summand, or of unbalanced
degree/order pairs.  In particular, it does not prove that \(\Delta^6\) is
globally minimal.

## 2. Invariant moment compiler

On \(\rho=1\), the radial powers in (1.2) disappear.  Put

\[
 g_{ij}=\langle v_i,v_j\rangle,
 \qquad D_i=\ell_i n_i,
 \qquad N=D_1+D_2+D_3.
\tag{2.1}
\]

Because every \(v_i\) is isotropic, Wick pairings have no loops.  For three
pairwise distinct isotropic directions the edge multiplicities are forced:

\[
 e_{12}=\frac{D_1+D_2-D_3}{2},\quad
 e_{13}=\frac{D_1+D_3-D_2}{2},\quad
 e_{23}=\frac{D_2+D_3-D_1}{2}.
\tag{2.2}
\]

The contribution is zero unless all three numbers in (2.2) are
nonnegative.  Otherwise normalized spherical Wick contraction gives

\[
 \mathcal R_\rho
 \left(L_1^{D_1}L_2^{D_2}L_3^{D_3}\right)
 =
 \frac{D_1!D_2!D_3!}
 {(N+1)!!e_{12}!e_{13}!e_{23}!}
 g_{12}^{e_{12}}g_{13}^{e_{13}}g_{23}^{e_{23}}.
\tag{2.3}
\]

All \(g_{ij}\) are nonzero on the distinct-direction chart.  Since the
\(\ell_i\) are even, their contribution can be absorbed without adjoining
roots by replacing

\[
\begin{aligned}
 a_1&\longmapsto
 a_1\left(\frac{g_{12}g_{13}}{g_{23}}\right)^{\ell_1/2},\\
 a_2&\longmapsto
 a_2\left(\frac{g_{12}g_{23}}{g_{13}}\right)^{\ell_2/2},\\
 a_3&\longmapsto
 a_3\left(\frac{g_{13}g_{23}}{g_{12}}\right)^{\ell_3/2}.
\end{aligned}
\tag{2.4}
\]

Consequently the complete invariant moment is the three-variable polynomial

\[
\boxed{
 \mu_m=
 \sum_{\substack{n_1+n_2+n_3=m\\2\max D_i\le N}}
 \binom{m}{n_1,n_2,n_3}
 \frac{\prod_iD_i!}
 {(N+1)!!\prod_{i<j}e_{ij}!}
 \prod_i a_i^{n_i}.
}
\tag{2.5}
\]

Formula (2.5), rather than an expansion in all coefficients of \(P^m\), is
the moment compiler used in the checker.  Each summand is hypergeometric in
the occupation numbers.  A coherent multiplier is compiled by adding its
stub degree before applying the same formula.

## 3. Two channels

For two nonproportional isotropic directions, a contraction can return to
the scalar channel only when

\[
 \ell_1n_1=\ell_2n_2.
\tag{3.1}
\]

Let \(g=\gcd(\ell_1,\ell_2)\).  At the first possible return,

\[
 n_1=\ell_2/g,\qquad n_2=\ell_1/g,
 \qquad m=(\ell_1+\ell_2)/g.
\tag{3.2}
\]

There is exactly one occupation vector, and its coefficient in (2.3) is
nonzero.  Hence two nonzero channels in distinct directions never have all
pure moments zero.

If the directions are proportional, both channels have the same positive
phase direction.  This is one of the terminal strata handled in Section 6.

## 4. Three distinct directions

Normalize \(a_1=1\) and saturate by \(a_2a_3\).  Modular Gröbner discovery at
\(p=101,103,107\) gives the same first unit cutoff at every prime.  Exact
Gröbner reduction over \(\mathbb Q\) then promotes each modular exclusion:

| harmonic degrees | exact unit cutoff |
|---|---:|
| \((2,4,6)\) | 5 |
| \((2,4,8)\) | 5 |
| \((2,4,10)\) | 3 |
| \((2,6,8)\) | 3 |
| \((2,6,10)\) | 6 |
| \((2,8,10)\) | 3 |
| \((4,6,8)\) | 4 |
| \((4,6,10)\) | 3 |
| \((4,8,10)\) | 4 |
| \((6,8,10)\) | 3 |

For example, saturation is imposed by adjoining

\[
 z a_2a_3-1.
\tag{4.1}
\]

The ideal generated by (4.1) and the moment numerators through the displayed
cutoff is the unit ideal over \(\mathbb Q\).  Thus the table is a collection
of characteristic-zero certificates, not a conclusion inferred from good
primes.

## 5. Two distinct directions

It remains to treat collisions of coherent profiles.  In a nondegenerate
three-dimensional quadratic space, two isotropic covectors are orthogonal if
and only if they are proportional.  Hence every boundary of the
distinct-direction chart is a partition with two or one distinct isotropic
directions.

Suppose channels \(i,j\) use one direction and channel \(r\) uses the other.
For an occupation vector \(n\), put

\[
 D=\ell_i n_i+\ell_jn_j.
\tag{5.1}
\]

The invariant contribution is zero unless \(D=\ell_rn_r\), and otherwise

\[
 \mathcal R_\rho(L^DM^D)=\frac{D!}{(2D+1)!!}\langle L,M\rangle^D.
\tag{5.2}
\]

This is a second, smaller invariant compiler.  For each of the ten harmonic
triples and each of its three collision pairs, the checker again discovers
the cutoff modulo \(101,103,107\) and certifies it over \(\mathbb Q\).  The
largest cutoff for each triple is:

| harmonic degrees | largest collision cutoff |
|---|---:|
| \((2,4,6)\) | 3 |
| \((2,4,8)\) | 6 |
| \((2,4,10)\) | 4 |
| \((2,6,8)\) | 4 |
| \((2,6,10)\) | 8 |
| \((2,8,10)\) | 5 |
| \((4,6,8)\) | 4 |
| \((4,6,10)\) | 5 |
| \((4,8,10)\) | 5 |
| \((6,8,10)\) | 5 |

The saturation covers three nonzero channel coefficients.  Its coefficient
boundaries reduce either to the two-channel calculation in Section 3 or to a
single direction.

## 6. The all-order terminal stratum

After Sections 3--5, vanishing of the first nine moments forces every active
\(v_i\) to be proportional.  Move their common isotropic covector to the
split coordinate \(x\).  Then

\[
 P=\sum_i a_i\rho^{(d-\ell_i)/2}x^{\ell_i}.
\tag{6.1}
\]

Every term has strictly positive phase weight
\(\operatorname{wt}(x)=1\), \(\operatorname{wt}(y)=-1\), and
\(\operatorname{wt}(t)=\operatorname{wt}(\rho)=0\).  The Laplacian preserves
this weight.  If \(q=\deg Q\) and
\(\ell_{\min}\) is the smallest active harmonic degree, then every term of
\(QP^m\) has weight at least \(m\ell_{\min}-q\), whereas
\(\Delta^{dm/2}(QP^m)\) has degree at most \(q\).  It therefore vanishes as
soon as

\[
 m\ell_{\min}>2q.
\tag{6.2}
\]

This proves the GVC conclusion on the entire residual stratum and completes
Theorem 1.1.

## 7. A genuine near survivor and the promotion gate

The calculation does find a nontrivial prefix survivor in the
\((2,4,6)\) chart.  With \(a_1=1\), put

\[
 a_3=-\frac{143}{60},
 \qquad 56a_2^2+272a_2+85=0.
\tag{7.1}
\]

The first two moments vanish by harmonic-degree constraints, while

\[
 \mu_3=\frac{8a_1a_2(143a_1+60a_3)}{15015},
\tag{7.2}
\]

\[
 \mu_4=
 \frac{64a_1a_3(85a_1^2+272a_1a_2+56a_2^2)}{255255}.
\tag{7.3}
\]

Thus \(\mu_1,\ldots,\mu_4\) vanish at (7.1).  The \(\mathcal H_2\)
multiplier obtained by inserting \(L_1^2\) already has

\[
 \nu_2=\mathcal R_\rho(L_1^2P^2)=\frac{8a_2}{315}\ne0.
\tag{7.4}
\]

So this point has exactly the desired finite-prefix shape: cancellation in
the scalar Reynolds channel and survival in a multiplier channel.  It is not
promoted.  Its fifth pure moment is

\[
 \mu_5=
 \frac{64a_2(3146a_2^2-96600a_2-163875)}{58503375},
\tag{7.5}
\]

and the numerator in (7.5) is coprime to the quadratic in (7.1).  Hence both
algebraic points die at order five.

The general moment formula (2.5) is a finite hypergeometric sum, and the same
compiler supplies every multiplier moment.  Nevertheless the required pure
all-order premise fails before any closed-form, rational-diagonal, or
P-recursive promotion question arises.  No candidate from this search is
reported as an all-order survivor.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_gvc3_isotropic_harmonic_channels.py
```

The checker compiles the invariant moments, performs the three-prime modular
discovery, repeats all ten distinct-direction and thirty collision-chart
unit calculations exactly over \(\mathbb Q\), verifies the two-channel
return formula, and checks the near survivor and its multiplier defect.  It
writes
[`gvc3_isotropic_harmonic_channels.json`](../artifacts/generated-results/gvc3_isotropic_harmonic_channels.json).
