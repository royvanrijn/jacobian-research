# A computable absolute eventual cutoff for Multivariate SFC(2)

Let \(a,b\) be distinct nonnegative integer exponent columns of
any finite dimension and
\[
 Q_j(t)=\sum_{k=0}^j\binom jk
        \prod_i((j-k)a_i+kb_i)!\,t^k .
\]

**Theorem.** An absolute integer \(J\) is computable such that
\(\gcd(Q_{N-1},Q_N)=1\) for every \(N\ge J\), uniformly in
\(a,b\) and the coefficient field of characteristic zero.
The computation is specified by the cutoff constructions
below. No numerical value of \(J\) is evaluated here.

The [effective global reduction](MULTIVARIATE_EFFECTIVE_GLOBAL_REDUCTION.md)
excludes all private profiles and confines the remaining
candidate common roots to fully shared curvature
\(W\in[W_0,128]\), where \(W_0>0\) is an absolute computable
rational. For any prescribed rational \(\kappa>0\) it gives
a computable order above which the remaining profiles satisfy
\(|\delta|<\kappa\), where \(\delta=\max\rho-1\).
The [effective quartic entry](MULTIVARIATE_EFFECTIVE_SHARED_CRITICAL_ENTRY.md)
gives a computable joint entry to any prescribed physical
quartic field radius if \(|\delta|\) is sufficiently small.
The remaining work is therefore to cover all scales within
that neighborhood by a single computable cutoff hierarchy.

## Linked field scales

Use \(q_*,w_*,a_2,g\) from the
[effective bounded quartic theorem](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md)
on \(W\in[W_0,128]\), and \(\beta_j=w_*-w_j\).
Put
\[
      \ell=(-4g)^{-1/4},\quad \ell_0=\ell^{4/3},\quad
      T_j=-2a_2\ell^2\sqrt j,\qquad H_j=\ell\,j^{3/4}\beta_j .
\]
The effective bounds for \(q_*-q_m\) and \(\rho''\) give
\(\rho(q_*)-1=\delta+O(\delta^2)\) with computable constants.
After decreasing the allowable \(|\delta|\), \(a_2\) has the
sign of \(\delta\) and has absolute value between computable
positive multiples of \(|\delta|\). In particular it is zero
at \(\delta=0\).

For \(\beta_j\ne0\), the two criteria in the
[escaping-field theorem](MULTIVARIATE_EFFECTIVE_SHARED_ESCAPING_QUARTIC_FIELDS.md)
use exactly
\[
   j\ell_0|\beta_j|^{4/3}=|H_j|^{4/3},\qquad
   |a_2|\ell_0/|\beta_j|^{2/3}
              =|T_j|/(2|H_j|^{2/3}).
\tag{1}
\]
For \(a_2\ne0\), set \(r=\sqrt{2|a_2|\ell_0}\),
\(c=\ell_0r^4\), \(M_j=jc\), \(b_j=\beta_j/r^3\).
The thermal field and mass at the two orders satisfy
\[
                     M_j=T_j^2,\qquad H_j=b_j M_j^{3/4}.
\tag{2}
\]
On the subcritical side the physical length in the fold
proof is \(\ell_0 r\). Its field \(b_j\) and mass \(M_j\)
are those in (2). The mass coefficient \(c\) is bounded
above and below by computable positive multiples of
\(\delta^2\) on both signed sides.

## A finite hierarchy of cutoffs

First obtain the escaping-field constants
\(d_e,b_e,\tau_e,M_e,n_e\) on this \(W\) band by the
escaping-field theorem. Both physical fields will be
kept in \(|\beta_j|\le b_e\) by quartic entry at the
end of the construction.

Choose a rational \(B\ge1\) with \(B^{2/3}\ge1/\tau_e\).
Obtain the [effective subcritical fold](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_FOLD.md)
field radius \(b_f\), physical length cap and mass threshold.
Apply the
[subcritical regular thermal theorem](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_REGULAR.md)
with this \(B\) and some rational \(0<\varepsilon<b_f\).
It gives a computable cap on \(|\delta|\) and lower bound
on \(n\delta^2\) for localization to the same fold disk.
Since \(M_n\le C n\delta^2\) for a computable \(C\), a
computable lower bound on \(M_n\) suffices for that localization.
Apply the
[supercritical bounded thermal theorem](MULTIVARIATE_EFFECTIVE_SHARED_BOUNDED_THERMAL_CUSP.md)
with the same \(B\) to obtain its curvature cap, mass threshold,
and order threshold.

Choose a rational \(M_{\rm th}\) at least \(8M_e\), at least
one, and at least all three thermal mass thresholds just
obtained. Consider a prospective common root with \(M_n\ge
M_{\rm th}\). Here \(a_2\ne0\).
By (1) and (2), if either \(|b_j|\ge B\), then at that order
\[
       j\ell_0|\beta_j|^{4/3}=M_j |b_j|^{4/3}\ge8M_e,
       \qquad
       \frac{|a_2|\ell_0}{|\beta_j|^{2/3}}
                  =\frac1{2|b_j|^{2/3}}\le\frac{\tau_e}{2}.
\]
The escaping theorem excludes it. Thus both fields lie
in \(|b_j|\le B\). On the supercritical side the bounded
thermal theorem then excludes the common root.
On the subcritical side the regular theorem puts both
fields in the same prescribed fold disk. Shrinking the
admitted curvature margin computably also gives the fold
physical length cap, since \(\ell_0 r=O(\sqrt{|\delta|})\).
The effective fold theorem excludes a common root there.

In the other case \(T_n^2<M_{\rm th}\), including \(a_2=0\),
both \(|T_j|\) are bounded by some rational \(B_T\)
with \(B_T^2\ge2M_{\rm th}\).
Choose a rational \(B_H\ge1\) so that
\[
             B_H^{4/3}\ge8M_e,\qquad
             B_H^{2/3}\ge B_T/\tau_e.
\]
By (1), the escaping theorem excludes a root with either
\(|H_j|\ge B_H\). The remaining \(T_j,H_j\) therefore lie in
the bounded quartic window at both orders.
Use its computable curvature cap and order threshold for
these \(B_T,B_H\).

Finally apply the critical-entry theorem with physical radius
\(b_e\) to obtain a curvature cap and an order threshold.
Take a positive rational \(\kappa<1\) at most \(d_e\) and
all the specified curvature caps for the center, thermal,
critical-entry, and bounded quartic estimates. It is also
small enough for the physical length condition just used.
Take the global reduction cutoff at this \(\kappa\) and the
maximum of the order bounds in this hierarchy, increasing by
one where a lower-order bound \(n\) rather than a terminal
bound was used. Above that maximum the global reduction
and entry guarantee the hypotheses at the start of the
calculation. The two scale cases exhaust the possibilities.
Since the coefficients are integers, complex coprimality
gives a rational Bézout identity and transfers to every
characteristic-zero field.

The [terminal-eight theorem](MULTIVARIATE_TERMINAL_EIGHT_ALL_PROFILES.md)
and the earlier fixed-terminal results independently cover
all terminals through eight. Completing the whole conjecture
still requires covering the fixed terminals from nine below
the absolute cutoff, or a direct proof independent of the
cutoff. The threshold calculation may be large even though
each step terminates.
