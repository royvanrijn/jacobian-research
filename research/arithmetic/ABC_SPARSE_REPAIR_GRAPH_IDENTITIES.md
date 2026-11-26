# Sparse repair identities on finite additive graphs

A fixed additive configuration can force identities in the
original complete labels while forcing no sparse identity that
distinguishes their high-power repair factors. This distinction
can be made exact for every finite configuration that has a
positive primitive realization. The proof uses the same actual
least repair as the primitive quotient count.

## Which identities survive every primitive realization

Let a finite tripartite hypergraph have \(m\) complete labels
\(\ell_1,\ldots,\ell_m\). On every edge require \(a+b=c\);
labels in the same part are distinct. Assume at least one
realization in positive integers with \(\gcd(a,b)=1\) on
each edge. Write \(\mathcal V\) for the rational linear
space cut out by the edge equations, with coordinates the
complete labels. The functions \(\ell_i\) restricted to
\(\mathcal V\) are pairwise nonproportional. In fact uniform
translation by \((u,v,u+v)\) on the three parts is allowed.
For same-part proportional functions the translation would
force proportionality constant one, inconsistent with the
given realization. The three translation directions as
functions of \((u,v)\) are pairwise nonproportional across
the parts.

For \(M\ge2\), let \(\mu_M(t)\) be the least multiplier
that makes \(t\) \(M\)-full and \(H_M(t)=t\mu_M(t)\).
Fix \(s\ge2\) and rational functions \(A_\nu\) on
\(\mathcal V\) of degree at most \(D\) in both numerator
and denominator, all nonzero. Fix integer arrays
\(\alpha_{\nu,i},\beta_{\nu,i},\gamma_{\nu,i}\), allowing
Laurent exponents and taking \(E\ge\max_{\nu,i}
|\gamma_{\nu,i}-\beta_{\nu,i}|\).
Consider the identity
\[
 \sum_{\nu=1}^s A_\nu(\boldsymbol\ell)
    \prod_{i=1}^m H_M(\ell_i)^{\alpha_{\nu,i}}
                 \mu_M(\ell_i)^{\beta_{\nu,i}}
                 \ell_i^{\gamma_{\nu,i}}=0 .                 \tag{GI1}
\]
Group the indices \(\nu\) by the vector
\(\boldsymbol u_\nu=(\alpha_{\nu,i}+\beta_{\nu,i})_i\).
If
\[
          M>2E+2D+\binom{s-1}{2}(m+2sD+1),                  \tag{GI2}
\]
then (GI1) holds on every positive primitive realization
where its coefficients are defined if and only if, for
each occurring vector \(\boldsymbol u\),
\[
       \sum_{\nu:\,\boldsymbol u_\nu=\boldsymbol u}
          A_\nu(\boldsymbol\ell)
          \prod_{i=1}^m\ell_i^{\gamma_{\nu,i}+\alpha_{\nu,i}}
                       =0 \quad\hbox{on }\mathcal V .         \tag{GI3}
\]
The sufficiency uses \(H_M(\ell_i)=\ell_i\mu_M(\ell_i)\)
and holds at every \(M\). Thus the universal identities
in this range come from the original identities after
factoring a common repair monomial within each group.
In particular a universal nondegenerate relation cannot
involve two distinct vectors \(\boldsymbol u_\nu\).
The threshold reflects the number of distinct complete
labels and the degrees of the coefficient functions,
with no bound on the magnitudes of their coefficients.

## Primitive pencils and actual least repairs

Fix a positive primitive realization \(\boldsymbol w\).
Choose a positive integral direction \(\boldsymbol v\)
in \(\mathcal V\) so that the affine forms
\(L_i(t)=\ell_i(\boldsymbol w+t\boldsymbol v)\) have
nonzero slopes and pairwise distinct zeros. These
conditions exclude finitely many proper linear subsets
of \(\mathcal V\) for the slope choice. To ensure that
finitely many specified nonzero rational functions on
\(\mathcal V\) restrict to nonzero defined functions,
avoid their numerator and denominator zero loci on
the line. For example, for each nonzero polynomial
\(P\) it suffices to choose \(\boldsymbol v\) with
\(P(\boldsymbol w+\boldsymbol v)\ne0\). The required
conditions are proper algebraic exclusions in a nonempty
open positive cone. One can first choose a rational
direction and scale it to an integer, since scaling
keeps the line and the distinct-zero condition.
For common dyadic height shells, choose the direction
near the translation \((1,1,2)\) by roles, before scaling.
All targets then have slopes whose maximum is less
than twice their minimum.

Let \(\mathcal P\) be a fixed finite prime pool containing
all prime divisors of the slopes, all nonzero two-form
resultants, and all primes up to a fixed \(P_0=P_0(m)\).
Choose a fixed integer \(D_0\) supported on \(\mathcal P\)
such that at \(t=D_0 z\) each \(L_i(D_0 z)\) has its
\(t=0\) valuations exactly at every \(p\in\mathcal P\).
For example use one higher power of each such prime
than the maximum order at \(t=0\) among the labels.
All edges are primitive throughout this progression.
Outside \(\mathcal P\) the resultant rules out a common
prime divisor of the two summands. Inside \(\mathcal P\)
their orders agree with the given primitive realization.

There is a positive proportion of integers \(z\) in
each interval \([Z,(1+\xi)Z]\), for any fixed \(\xi>0\)
and all sufficiently large \(Z\), for which every
\(L_i(D_0 z)\) is squarefree outside \(\mathcal P\).
For a prime outside the pool, the condition
\(p^2\mid L_i(D_0 z)\) eliminates one residue class
modulo \(p^2\). The interval count is at most its length
divided by \(p^2\), plus \(O(1)\). Only primes
\(p\ll_{\mathcal V,\boldsymbol v,\boldsymbol w,D_0,\xi}
\sqrt Z\) can occur. The sum of the \(O(1)\) errors is
\(O_m(\sqrt Z)\) with constants allowed to depend on
the pencil. Choose \(P_0\) so the sum of the main terms
is less than half the interval length. Discard also the finitely many same-part label
collisions. On the remaining values of \(z\), the bad-prime part \(B_i\) of each label is
fixed, and hence
\[
  \mu_M(L_i)=\mu_M(B_i)(L_i/B_i)^{M-1},\qquad
  H_M(L_i)=H_M(B_i)(L_i/B_i)^M .                             \tag{GI4}
\]
This uses the actual valuations and the unique least
multiplier. It works for every \(M\ge2\), including
orders below a fixed exceptional valuation.

## Proof of the universal identity classification

Use a pencil on which every nonzero left side of (GI3)
and all coefficient numerators and denominators stay
nonzero as rational functions. On the infinite set of
primitive parameters from (GI4), the terms of (GI1)
are values of rational functions of the form
\[
       f_\nu(t)=a_\nu(t)
             \prod_{i=1}^m L_i(t)^{
                          M u_{\nu,i}+\gamma_{\nu,i}-\beta_{\nu,i}},
\]
where \(a_\nu\) is the restriction of \(A_\nu\) multiplied
by a nonzero constant. Its numerator and denominator
have degrees at most \(D\). If (GI1) is universal,
the rational function \(\sum_\nu f_\nu\) vanishes
identically. Partition into minimal vanishing subsums.

For a minimal block, divide by any of its terms. All
ratios are units off the zeros of the \(m\) linear
forms, infinity, and the zeros and poles of the
coefficient functions in the block, so the support has
size at most \(m+1+2sD\). Brownawell–Masser on the
rational function field bounds the degrees of these
ratios by \(\binom{s-1}{2}(m+2sD+1)\). The conclusion
also holds trivially for a block of size two or a ratio
that is constant. At the zero of \(L_i\), the order of
the ratio \(f_\nu/f_{\nu'}\) differs from
\(M(u_{\nu,i}-u_{\nu',i})\) by at most \(2E+2D\).
Under (GI2) the valuation bound therefore forces equality
of the high exponent vectors throughout the block.
[Brownawell–Masser in the formulation of Fuchs–Zannier, Lemma 2.2](https://ems.press/content/serial-article-files/31810#page=8)

Every group of equal vectors is a union of such blocks
and has zero sum on the pencil. Factor out its common
repair monomial using (GI4), to see that the corresponding
left side of (GI3) is identically zero on the pencil.
By its choice this is possible only if all of (GI3)
holds on \(\mathcal V\). This proves necessity.

## Finite dictionaries on the same numerical scale

On the squarefree-outside-the-pool subset of the
primitive pencil, fix the integer exponent
arrays and a dictionary of \(J\) coefficient lists
of numerator and denominator degrees at most \(D\)
in the pencil parameter. The lists may vary with
the height. Suppose the high exponent vectors are
not all the same and \(M\) obeys (GI2). For any list,
if the resulting rational function vanishes identically,
the minimal-subsums argument makes the relation
degenerate at every parameter where it is defined.
Otherwise clearing denominators gives a nonzero
polynomial of degree bounded in \(M,s,m,D\) and the
fixed exponent arrays. Hence there are
\(O_{M,s,m,D,\alpha,\beta,\gamma}(J)\) parameters
where at least one list gives a nondegenerate relation.
The bound does not depend on the heights of
the coefficient entries.

After the height limit at fixed budget and order, one
can take any subset of these primitive realizations of
size \(n^\tau\) with \(0<\tau L<1\) and \(n=X^L\).
A coefficient dictionary of size \(n^{\tau-\eta}\)
cannot extract a nondegenerate relation on a
full-exponent subset in this model. The positive
primitive configurations have genuine exact additive
equations, including every trade in the prescribed
pattern. Their radicals are large: the variable parts
are squarefree and their numerical heights are
comparable to \(X\). Thus the example tests forcing
by the additive pattern and coefficient complexity.
The original small radical conditions would need to
supply further arithmetic that the pencil lacks.
[Least repair and primitive quotient count](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
· [Polynomial relations on repaired heights](ABC_REPAIRED_HEIGHT_POLYNOMIAL_RELATIONS.md)
