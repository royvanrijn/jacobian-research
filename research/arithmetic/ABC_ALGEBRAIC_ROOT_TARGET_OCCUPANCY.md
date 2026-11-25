# One-prime-power loss in the terminal reconstruction of affordable fibres

For a fixed affordable pair lock, a grouped shape on
the *actual third letter* gives a distinct-target
occupancy bound. The terminal p-adic step loses the
Jacobian valuation once instead of twice. It uses
algebraic roots over the local field closure, so it
requires neither an integral local lift nor a
Hensel inequality \(e>2v_p(J)\).

## Recovery through a nearest algebraic root

Let \(P\) be a nonzero univariate polynomial over
\(\mathbb Q_p\), and let \(z\in\mathbb Q_p\) satisfy
\(v_p(P(z))\ge e\). If \(P(z)\ne0\) and
\(d=v_p(P'(z))<e\), there is an algebraic root
\(\zeta\) of \(P\) with
\[
                         v_p(z-\zeta)\ge e-d .             \tag{AR1}
\]
Indeed, factor \(P\) over \(\overline{\mathbb Q}_p\)
with roots counted with multiplicity. The logarithmic
derivative identity gives
\[
 d\ \ge\ v_p(P(z))-\max_{\zeta:P(\zeta)=0}v_p(z-\zeta).
\]
If \(P(z)=0\), use the root \(z\) itself. Roots need
only be labeled inside a fixed algebraic closure at
each prime. Two rational numbers assigned to the
same labeled root at precision \(e-d\) agree to that
precision by the ultrametric inequality.

## Improvement to the one-root congruence

Keep the hypotheses and notation of (QD1) and the
first three conditions of (QD2) from the one-grouped-root
count. In particular \(q,R=\lfloor r/L\rfloor,M\)
are fixed before the inner height limit and \(q>R\);
\(b_Q\) is the actual normalized shallow-modulus weight.
The determinant first step is unchanged. For
clarity, write its remaining parameters as:
\[
 B=hr+EM,\qquad
 r<\liminf qL,\quad
 \alpha(M+1)r<\kappa,\quad
 b\alpha(M+1)r/2>B .
\]
The last requirement in (QD2) can now be replaced by
\[
                       \frac{B}{1-\alpha}<w .             \tag{AR2}
\]
The conclusion is unchanged: at most
\(X^{(b_Q+o(1))L}\) pairs, with the original fixed-budget
epsilon allowances. The fixed product amplification
also works: substitute \((\ell h,\ell E)\) in the
conditions, then obtain normalized cost \(b_Q/\ell\).
[One-root congruence (QD1)–(QD2)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-quantitative-congruence-estimate-for-one-grouped-root)

On an auxiliary factor \(F(z,v)=0\) of bidegree at
most \((R,M)\), use
\[
                          J=zF_z-qvF_v .
\]
By the argument of (QD1), the samples with \(J=0\)
are bounded in number: an irreducible torus curve
with identically zero logarithmic derivative would
have \(vz^q\) constant, impossible for \(z\)-degree
\(R<q\). The exact same coefficient bound gives
\(\log_X|\operatorname{num}J|\le B/(1-\alpha)+o(1)\)
at every remaining point.

For each \(p^e\Vert n\), take an integral representative
of \(K\) and form the univariate polynomial
\[
                        P(Z)=Z^{qM}F(Z,KZ^{-q}) .
\]
It is nonzero, even if the auxiliary coefficients
are arbitrary, because \(i+q(M-j)\) are distinct as
\(i\le R\), and \(K\) is a unit. At the original
point \(P(z)\equiv0\pmod {p^e}\), and the chain rule
gives
\[
                 P'(z)\equiv z^{qM-1}J(z,v)\pmod {p^e}.
\]
Thus whenever \(d_p=v_p(J)<e\), choose a root in
(AR1) and record its label and \(d_p\). For
\(d_p\ge e\), mark the prime as discarded. There
are boundedly many root labels per prime at fixed
parameters, and the depths \(d_p<e\) cost at most
the divisor count of \(n\) in all. Record count is
subpower in the inner height limit.

Two original \(z\)'s in one record agree modulo
the product of \(p^{e-d_p}\) over the recorded
primes. The total log weight lost is
\[
            \sum_{p\mid n}\min(e_p,d_p)\log_X p
               \ \le\ \log_X|\operatorname{num}J| .
\]
The denominator of \(J\) is a \(p\)-adic unit at
every prime dividing \(n\). Under (AR2) the
remaining precision has positive fixed numerical
weight, and cross multiplication recovers the
short rational \(z\) uniquely. Then \(v\) is
prescribed modulo the whole modulus.
For rational \(v\), its uniqueness requires
\(2E<w\), already implied by
\(E\le B<b\alpha(M+1)r/2\le w/2\).
For an integral wide coordinate the reconstruction
only needs \(E<w\).
The color argument and divisor multiplicity for
amplification are exactly those of (QD1).

## Distinct targets after two locks

Fix an affordable full-divisor pair \(F\mid U\),
\(G\mid V\) at total normalized root cost \(t\).
On the fibre write \(U=Fu\), let \(u\le X^{E_U}\),
and keep actual targets of a grouped form
\[
                  T=D W z^q,\qquad
                  z\le X^{hL+o(L)},\qquad W\le X^{E_W+o(1)} .
\]
Here \(D\) is fixed within a subrecord, \(W,z\) are positive
integers, and the cost of all additional \(D\)-subrecords
can be bounded separately. These restrictions can be
supplied by grouping target valuation slots; the target
radical cap is kept. Primitivity gives
\(\gcd(W,u)=1\) and unit residues of \(u,W,z,D\)
at every prime of \(G\).
Modulo \(G\) the additive equation gives, up to its
fixed sign,
\[
                  (W/u)z^q\equiv F/D\pmod G .
\]
Set \(v=W/u\) and \(E=\max(E_U,E_W)\). Exact values
\((z,v)\) determine the actual target \(T\):
the reduced numerator and denominator of \(v\)
recover \(W,u\), respectively. Thus (AR2) and the
three original determinant inequalities directly bound
the number of distinct targets in this fixed lock
fibre by \(n^{b_Q/\ell+o(1)}\) per subrecord, whether
or not the target fibre is a matching in the complete
letters. If the extra subrecord cost is \(s_0\), this
gives a strict saving below coefficient \(1/3\) as soon
as \(t+s_0+b_Q/\ell<1/3\) with a fixed margin.

For instance consider the following admissible strict
parameters in the revised calculation:
\[
 \begin{gathered}
 w=0.9,\quad \kappa=22.5,\quad b=b_Q=0.04,\quad
 h=0.12,\quad E=0.0081,\quad qL>1.7,\\
 t=0.065,\quad s_0=1/3-0.12,\quad
 \ell=1,\quad \alpha=1/2,\quad M=23,\quad r=1.7 .
 \end{gathered}
\]
They give \(B=0.3903\), slope precision \(20.4<22.5\),
Taylor lower bound \(0.408>B\), and
\(B/(1-\alpha)=0.7806<w\).
The chart count yields
\(t+s_0+b_Q=1/3-0.015\).
For this same \(h,E,b,w\), no choice of the old
(QD2) parameters works. In fact the Taylor inequality
and the old terminal inequality imply, with
\(x=\alpha(M+1)>2h/b=6\) and \(M\ge6\),
\[
 \frac Ew
   <\frac{1-\alpha}{2M}\!\left(1-\frac{6}{x}\right)
   \le\frac7{12}\frac{\alpha(1-\alpha)}x
                     \left(1-\frac6x\right)
   \le\frac7{1152}
   <0.009 .
\]
In the last inequality use \(x^{-1}(1-6/x)\le1/24\).
This separates the estimates for a prescribed packet,
not every alternative choice of grouping or records.

The fixed-lock argument uses the true target factorization
to supply \(z,W\). It requires a specific bound \(E\) on
both the short free remainder and the group's cofactor.
Affordability alone supplies a shorter remainder than
the opposite modulus but no such quantitative grouping
bound for an arbitrary third-letter valuation profile.
Therefore this improvement does not close the global
occupancy problem.

## With the primitive kernel estimate, the terminal condition disappears

The same proof applies to (PK1). Put
\[
 \Lambda=\frac{b(M+1)r}{2},\qquad
 G_*=(\sqrt{\Lambda}-\sqrt{\Lambda-B})^2.
\]
The height bound for the logarithmic derivative is
now \(B+G_*\). The root proof thus replaces the terminal
condition of (PK1) by \(B+G_*<w\). This is automatic:
\[
 G_*\le B,\qquad
 2B<2\alpha\Lambda=b\alpha(M+1)r<b\kappa\le w.
\]
The last inequality is implied by the defining
condition on the shallow divisor \(Q\mid n\).
The wide rational reconstruction is automatic as
well since \(E\le B\). Therefore a rank-forcing
rectangle with \(R<q\) gives whole-chart occupancy
at bounded multiplicity up to divisor records.
This remains true for a relative shallow divisor in
an arbitrary full modulus. By moving \(\alpha\) close
to one and, if needed, shortening \(r\) while holding
\(\alpha r\) fixed, the exact strict test for this
rank forcing is that some integer \(M\ge1\) satisfies
\[
 E<
 \frac{b(M+1)/2-h}{M}\min\!\left(a,\frac{\kappa}{M+1}\right),
 \qquad b(M+1)>2h,\quad a=\liminf qL.
 \tag{AR3}
\]
The converse here concerns the displayed rank-forcing
inequalities only.

A parameter separation from even the optimized old (PK1)
takes the same fixed-lock costs and parameters \(h,b,w,\kappa\)
above, but \(E=0.0144\), \(M=13,r=1.5,\alpha=0.95\).
Here \(B=.3672<.399=\alpha\Lambda\), and
\(\alpha(M+1)r=19.95<22.5\). Thus the direct occupancy
saving is still \(0.015\).
For every old (PK1) choice of \(M,r,\alpha\) on this
prescribed packet, its terminal inequality requires
\[
 \frac Ew <
 \frac{1-\sqrt{6/(M+1)}}{2M}
 \le\frac2{135}<.016 .
\]
For the first inequality, the unconstrained minimum
of \(B+G_*\) over \(r\) is
\(EM/(1-\sqrt{h/A_0})\), where \(A_0=b(M+1)/2>h\).
This is the primitive-kernel minimum already derived
in (PK1). For the second inequality put
\(x=\sqrt{6/(M+1)}\). Since
\(1/M=x^2/(6-x^2)\le x^2/5\), the expression is at most
\(x^2(1-x)/10\le2/135\).
The same necessary bound holds for a smaller lower
bound than \(b=.04\).

As a direct profile application, keep both opposite
hypotheses of the fixed-difference sector (PK3) and
allow the entire target support in \([2.8,3.2]\) with
\(\int|r-3|\,d\nu_T(r)\le.0155\). This excludes the
minimum nonatomic endpoint in that wider sector:
use the same fixed rows, shallow divisor of weight
\(b_Q\le.1401\), \(\kappa=4\), \(qL\to3\), and same
unit ratio after removing the \(D\)-supported parts.
Even ignoring any common-gcd gain, the bounds are
\(h\le1/3+o(1)\), \(E\le.0155+o(1)\).
Take \(M=9,r=.395,\alpha=.99,b=.1399\). Then
\[
 \alpha(M+1)r=3.9105<4,\qquad
 B\le.2711667+o(1)<.273539475=\alpha\Lambda .
\]
The chart count is \(b_Q\le.1401<1/6\), and the
fixed difference recovers the column as in (PK3).
No upper bound on the numerical size of the removed
common gcd is needed for this extension. In particular
it excludes the same formal opposite profiles and
target \(\frac16\delta_3+\frac16\operatorname{Unif}[2.814,3.186]\).

## Obstruction for the symmetric diffuse profile persists

For the symmetric test profile from (SB1)–(SB4),
\(\nu_U=\nu_V=\nu_T=\frac16\delta_3+
\frac16\operatorname{Unif}[2,4]\), the strengthened
one-root rank condition still gives no strict saving
when counting the frozen records in the same way.
This statement allows arbitrary retained
submeasures on both free letters, integer group
multipliers, modulus and shallow divisors and
fixed amplification. It is a formal obstruction
to that counting criterion.

Use the existing notation: \(1/3+g\) pooled retained
mass, \(t\ge g\) diffuse part of it, \(\mu\) lock mass,
\(b'=b/\ell\) effective chart cost after amplification,
\(a\) grid spacing, and \(h,E\) unamplified widths.
Success needs \(\mu+b'<g\) and \(w\le4\mu\).
The first three determinant conditions with
\((\ell h,\ell E)\) give
\(EM<w/(2\ell)\le2\mu\).
Let \(t_0\) be retained mass assigned multiplier zero.
Every degree in the exact profile is at least two, so
\(t_0\le E\), since the two signs together have wider
weight at most \(2E\). The nonzero multipliers each
have absolute value at least one, hence
\(2h\ge1/3+g-t_0\).
If \(M=1\), rank forcing gives \(2h<2b'\), and then
\[
 E>1/3+g-2b'>1/3-g+2\mu\ge2\mu,
\]
a contradiction. Thus \(M\ge2\), \(E<\mu\), and
\(2h-b'>1/3\).

Set \(x=b'(M+1)/(2h)>1\). Rank and precision give
\[
 E<\frac{w}{2\ell M}(1-1/x)
 \le\frac{b'w}{8\ell(2h-b')}
 \le\frac{g^2}{8(2h-b')}.
 \tag{AR4}
\]
For the first inequality combine
\(EM< (b'\alpha(M+1)/2-h)r\) with
\(r<\kappa/(\alpha(M+1))\), \(b\kappa\le w\),
and \(\alpha\le1\).
For the second use \(1/M\le
b'/((2h-b')x)\) and \((x-1)/x^2\le1/4\).
The final inequality follows from \(b'w\le4b'\mu<g^2\).

If \(a\ge2\), projection of the pooled diffuse density
\(1/6\) on \([2,4]\) modulo \(a\) has density at
most \(1/6\) a.e. in the original degree units.
Folding to distance from \(a\mathbb Z\) gives density
at most \(1/3\). The depth integral among points of
total mass \(t\) is at least \(3t^2/2\); each
actual integer multiplier has correction at least
nearest-grid distance. Therefore \(E\ge3t^2/4\),
contrary to (AR4), which gives \(E<3g^2/8\).

If \(0<a<2\), (SB1) gives
\[
 2h-b'\ge\frac{11/12+3g-2E-ab'}{a}
            >\frac{11/12+g}{a}
\]
using \(E<\mu\) and \(\mu+b'<g\).
Now (AR4) gives \(E<(3/22)ag^2\). But the general
folded-density bound (SB2) gives
\(E\ge3a t^2/(4(2+a))>(3/16)ag^2\), another
contradiction. Hence there is no fixed-spacing
one-group closure even when the final determinant
reconstruction conditions have been eliminated.
