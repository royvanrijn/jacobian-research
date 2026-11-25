# Three-target stars with a small prime-order defect

Three actual targets with prescribed relative odds can be counted
jointly by a genus-at-least-two curve. A prime cyclic cover uses
three target factors, so neither summand needs a power condition.
This gives a new bound on common-neighbor columns of the actual
outside-factor graph.

Fix an odd prime \(p\). For \(v\ge1\), write
\(\operatorname{rad}_p(v)=\prod_{\ell:\ p\nmid v_\ell(v)}\ell\).
For any prime set \(S\), a subscript \(\ell\notin S\) on a radical
means retain only the factors outside \(S\).

## Count for three actual targets at prescribed odds quotients

Consider three ordered primitive positive triples
\[
 T_i=\sigma U_i+\tau V_i>0,\qquad
             U_i,V_i,T_i\le X,\qquad i=1,2,3,
\]
with a common choice of signs \(\sigma,\tau\in\{1,-1\}\).
Prescribe the relative complete odds
\[
               r_i=\frac{U_i V_1}{U_1 V_i},\qquad r_1=1,
\]
with \(r_1,r_2,r_3\) pairwise distinct. The rationals may vary
with \(X\). If the prescribed data can occur, their reduced
numerators and denominators are at most \(X^2\). Fix any prime
set \(S\) whose product is at most \(X^A\) for fixed \(A\).
Uniformly in these data, the number of ordered three-tuples
of triples with
\[
   \left.\operatorname{rad}_p(T_1T_2T_3^{p-2})\right|_{\ell\notin S}
                                                       \le B
                                                               \tag{PT1}
\]
is \(O_{p,A,\varepsilon}(B X^\varepsilon)\) for \(B\ge1\).
There is no condition on the odd or other valuations on
\(U_i,V_i\). One can use the product of the three individual
prime-order defect radicals in place of the mixed one in (PT1).

To prove this, put \(x=\gcd_i U_i,\ y=\gcd_i V_i\).
They are coprime. At a prime \(\ell\) outside the support of
the numerators and denominators of the \(r_i\),
\[
                v_\ell(U_i)-v_\ell(V_i)
                   =v_\ell(U_1)-v_\ell(V_1).
\]
Within each primitive edge, at most one of the two valuations
is positive, so the common difference determines both of them.
Thus the quotients \(d_i=U_i/x,\ e_i=V_i/y\) have support on a
fixed integer of polynomial height. The uniform fixed-support
divisor estimate gives \(X^{o_X(1)}\) choices of their values.
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Fix these values. The actual targets are the three binary forms
\(L_i(x,y)=\sigma d_i x+\tau e_i y\), and the prescribed distinct
odds ensure that the forms are pairwise nonproportional. Express
\(L_3=aL_1+bL_2\) with \(a,b\in\mathbb Q^\times\). Their numerators
and denominators are polynomial in height by Cramer's rule.
Write \(T_1T_2T_3^{p-2}=\Delta w^p\), where \(\Delta\) is
positive and \(p\)-th-power-free. There are at most
\(B X^{o_X(1)}\) choices of \(\Delta\): its outside radical
has no more than \(B\) possibilities, its support inside \(S\)
has at most \(2^{\omega(\prod_{\ell\in S}\ell)}\) choices, and
each occupied prime has at most \(p-1\) choices of exponent.
All relevant integers have polynomial height for the fixed \(p,A\).

Now set \(s=aL_1(x,y)/L_3(x,y)\) and \(z=w/L_3(x,y)\).
The complete additive relations give
\[
                     s(1-s)=ab\Delta z^p.                   \tag{PT2}
\]
The variable \(s\) is a nonconstant fractional linear
function of \(x/y\), and hence determines \((x,y)\) by
primitivity and positivity. The counting problem after fixing
the data therefore reduces to rational points on (PT2).

## Uniform subpower bound for the cyclic twists

For fixed prime \(p\ge5\), let \(h\in\mathbb Q^\times\) have
numerator and denominator of size at most \(X^C\) for fixed \(C\).
The smooth projective curve \(C_h\) of \(s(1-s)=h z^p\) has
\[
                       \# C_h(\mathbb Q)
                               \ll_{p,C,\varepsilon}X^\varepsilon .
                                                               \tag{PT3}
\]
Here the input is uniform Mordell–Lang for genus at least two,
together with rational isogeny torsion over a fixed cyclotomic
field. The argument does not use point bounds in varying
number fields of controlled discriminant.

Over \(K=\mathbb Q(\zeta_p)\) let \(\theta\) act on the curve by
\((s,z)\mapsto(s,\zeta_p z)\) and on its Jacobian \(J\).
The map with coordinate \(s\) has degree \(p\), with total
ramification at \(0,1,\infty\), by the orders \(1,1,-2\) of
\(s(1-s)\). Riemann–Hurwitz gives genus \((p-1)/2\).
The quotient by \(\theta\) is the projective \(s\)-line,
so the trace operator \(1+\theta+\cdots+\theta^{p-1}\) on
\(J\) is zero. Let \(\pi=1-\theta\). The cyclotomic identity
\((1-\zeta_p)^{p-1}=p\) times a unit shows
\(\pi^{p-1}=[p]\) times an automorphism. Indeed the trace
identity permits substitution in \(\mathbb Z[\zeta_p]\).
This proves that \(\pi\) is an isogeny of degree \(p\),
since \([p]\) has degree \(p^{2g}=p^{p-1}\).

The totally ramified points \(P_0,P_1\) at \(s=0,1\) are
rational, fixed by \(\theta\), and distinct. Their divisor
difference gives a nonzero \(K\)-point of \(J[\pi]\).
It is nonzero since a function with divisor \(P_1-P_0\)
would have degree one on a positive-genus curve.
It generates \(J[\pi]\), which has order \(p\). This is the
three-branch case of cyclic Jacobian descent.
[Poonen–Schaefer, cyclic covers and explicit descent](https://math.mit.edu/~poonen/papers/descent.pdf)

Write \(h=m/n\) with coprime signed numerator and positive
denominator. In the hyperelliptic model of (PT2), \(u=1-2s\)
gives \(u^2=1-4h z^p\). At places of \(K\) outside the finite
set \(\Sigma\) over \(2pmn\), the polynomial on the right has
unit leading coefficient and separable reduction. On the
chart at infinity \(z'=1/z,\ u'=u(z')^{(p+1)/2}\), its
equation is \((u')^2=(z')^{p+1}-4hz'\), with a simple root
at \(z'=0\). The two charts give good reduction of the
curve and its Jacobian. The isogeny extends
to the abelian scheme and is étale there, since its degree
is the residue-unit \(p\). Thus the isogeny Kummer image of
\(J(K)/\pi J(K)\) is unramified outside \(\Sigma\). Since
\(J[\pi]\) is constant and \(K\) contains \(\mu_p\), it lies
in the Kummer classes in \(K^\times/K^{\times p}\) with
valuations divisible by \(p\) outside \(\Sigma\).
The usual unit and class group sequence bounds the
\(\mathbb F_p\)-dimension of these classes by
\(|\Sigma|+O_p(1)\). In detail, a class gives the
\(p\)-torsion class of its fractional ideal after ignoring
\(\Sigma\), with kernel given by \(\Sigma\)-units modulo
\(p\)-th powers. The \(\Sigma\)-class group is a quotient of
the class group of the fixed \(K\), and the unit rank goes
up by at most \(|\Sigma|\).

On the free quotient of \(J(K)\), the cokernel of \(\pi\)
has order \(p^{r/(p-1)}\) if that quotient has rank \(r\).
This follows by taking determinants in
\(\pi^{p-1}=[p]\) times an automorphism. It is a quotient
of \(J(K)/\pi J(K)\), so
\(r\le(p-1)(|\Sigma|+O_p(1))\).
Use \(P_0\) for the Abel–Jacobi embedding of \(C_h\) in \(J\).
Uniform Mordell–Lang gives at most \(C_p^{1+r}\) points in
\(J(K)\) on that curve. Since \(|\Sigma|=O_p(1+\omega(mn))\)
and \(mn\) in absolute value is polynomial in \(X\), this is
\(X^{o_X(1)}\) uniformly in \(h\). This proves (PT3), and,
summed over the coefficient and mixed-defect records, proves
(PT1).
[Gao–Ge–Kühne, Theorem 1.1](https://arxiv.org/pdf/2105.15085v4)

## The cubic variant from the linear relation of the targets

The statement (PT1) also holds for \(p=3\), with
mixed product \(T_1T_2T_3\). After fixing \(d_i,e_i\),
put into the fixed prime support all primes dividing
the three nonzero pair determinants of the \(L_i\).
Any prime dividing two \(L_i(x,y)\) must be on that
support, since it divides the corresponding determinant
times both \(x\) and \(y\), and \(\gcd(x,y)=1\).
Thus off this support, a prime occurs on only one target.
Its nonzero valuation modulo three is detected in the
mixed defect. Recording the three cube-free parts \(k_i\)
of the targets costs at most \(B X^{o_X(1)}\) there,
by allocating each outside prime to its target and taking
one of the two nonzero exponents. On the fixed support
the records cost only \(X^{o_X(1)}\). This applies when
the pool \(S\) in (PT1) is already present.

The three \(L_i\) satisfy a linear relation with fixed
nonzero integer coefficients and polynomial height.
Writing \(T_i=k_i v_i^3\) gives a signed diagonal cubic
in the three roots. Their common divisor cubed divides
the nonzero determinant of any two \(L_i\), again by
\(\gcd(x,y)=1\), so it has divisor-many possibilities.
After dividing out this common divisor, apply Heath-Brown's
uniform primitive diagonal cubic bound at root height
at most \(X\). The two independent \(L_i\) recover
\((x,y)\) from the target values. The completion count is
\(X^{o_X(1)}\) for each cube-free record. This proves
(PT1) for the cubic variant, and the graph and endpoint
arguments below then apply verbatim.
[Heath-Brown, Theorem 1.3](https://arxiv.org/abs/2601.07817)

## The unweighted product for any fixed order at least three

For any fixed integer \(q\ge3\), in the statement of (PT1)
one can replace the left side by
\[
                \left.\operatorname{rad}_q(T_1T_2T_3)
                                                \right|_{\ell\notin S}.
                                                               \tag{PT6}
\]
Here \(\operatorname{rad}_q\) uses \(q\nmid v_\ell\).
If \(q\) has an odd prime divisor \(p\), use (PT1) at \(p\).
Off the fixed pair-determinant support a prime occurs on
at most one of the three targets. The three weights
\(1,1,p-2\) are units modulo \(p\). Thus the outside
\(p\)-defect in (PT1) is bounded by the outside
\(q\)-defect in (PT6), after adding the pair-determinant
primes to the fixed pool. This pool may be chosen before
enumerating the cofactors: for reduced \(r_i=a_i/b_i\),
every determinant prime divides the product of the
\(a_i b_i\) and the nonzero differences \(a_i b_j-a_j b_i\),
as in the mixed squarefree count.

If \(q\) is a power of two, use \(4\mid q\). In the proof
for the cubic variant replace cube-free records by
fourth-power-free records, with the same counting cost.
The linear relation of the three targets now gives a
diagonal quartic in three positive roots. After fixing
and dividing out their gcd as before, view it as
\(A v_1^4+B v_2^4+C v_3^4=0\), with all coefficients
nonzero and of polynomial height. Put
\(t=v_2/v_3,\ w=(v_1/v_3)^2,\ a=-B/A,\ b=-C/A\).
Then \(w^2=a t^4+b\), so the rational point
\[
                         (\xi,\eta)=(a t^2,a t w)
\]
lies on \(\eta^2=\xi^3+ab\,\xi\). This elliptic curve has
the rational point of order two \((0,0)\), and its
coefficients and the \(x\)-coordinates of the displayed
points have polynomial height. Dujella's bounded-height
elliptic estimate gives \(X^{o_X(1)}\) such points,
uniformly in these data. From \(\xi,\eta\) the ratios
of the primitive positive roots have bounded fibres
by the displayed formulas, and a primitive positive
root vector is determined by its projective point.
The gcd had only divisor-many possibilities because
its fourth power divides a fixed pair determinant.
This proves (PT6) for \(4\mid q\).
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/pdf/2312.03655v1)
· [Pair-determinant support](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#fixed-cross-ratios-and-the-common-factors)

## The star bound and a strict limiting exclusion

Consider a simple bipartite graph with one actual primitive
witness per edge as above. Take \(R\) rows, \(C\) columns and
\(E\) edges. Suppose for every ordered triple of distinct rows,
the relative-odds vector of their three edges at any shared
column has at most \(M\) possible values. For any such three
edges the witness triples are distinct, so their complete
odds are distinct by primitivity. Assume that the complete triple determines the two edge
labels, and that (PT1) holds for every such star with the same cap \(B\) and a
fixed prime pool, or with a pool fixed by the row tuple and
the relative-odds vector. Distinct shared columns give
distinct ordered three-tuples of actual witnesses. Hence
\[
    E\ll C+(R^3 C^2 M B)^{1/3}X^\varepsilon.                \tag{PT4}
\]
For \(E\ge4C\), convexity or the power mean gives
\(\sum_v\binom{\deg v}{3}\gg E^3/C^2\) on the column side.
The upper bound is \(O(R^3 M B X^\varepsilon)\) by (PT1).
The branch \(E<4C\) proves the remaining case.

For the endpoint application, use the original outside
factor graph of \(U,V\), with limiting nonatomic masses
\(c_U,c_V,c_T\) in base \(n=X^L\). Assume \(c_T>0\) and
that every selected target satisfies
\(\operatorname{rad}_q(T)\le n^{f+o(1)}\) for a fixed \(f\)
and a fixed \(q\ge3\), using (PT6) in the graph bound.
Here the primes in the three radical caps still come from
the original primitive triples. Choose a fixed full block
on \(T\) of arbitrarily small normalized cost, using a
positive compact band of its nonatomic measure. Fix the
atomic windows and use (TR2)–(TR3) to count complete odds
on each outside pair and relative odds for each pair of
rows at shared columns. Both counts cost arbitrarily
small normalized exponent after amplification. Choose
one witness per outside pair. Then \(R\le n^{c_U+o(1)}\),
\(C\le n^{c_V+o(1)}\); two relative-odds records cost
arbitrarily small exponent; and \(B\le n^{3f+o(1)}\)
by the individual target caps. Taking the inner height
limits first, then taking outer allowances to zero,
(PT4) gives
\[
                   d\le\max(c_V,\ c_U+\tfrac23c_V+f).
                                                               \tag{PT5}
\]
Exchanging the roles of \(U,V\) gives the reverse
orientation. This uses the bounded full modulus count
only to specify a vector of actual ratios. The power
saving in the star comes from the third-letter radical
constraint through the cyclic curve.
[Complete cross-ratio count (TR2)–(TR3)](ABC_DIFFUSE_FACTOR_GRAPH.md#counting-complete-cross-ratios-by-a-full-modulus)

For example suppose \(d=1/3\), the total nonatomic mass
is at most \(1/2+1/60\), and only the target has
\(f\le1/72\). At the balanced endpoint the pair
nonatomic masses are at least \(1/3\) when all three
are positive. If one vanished the other two would
both be \(1/3\), excluded in this range. Thus each
\(c_i\le1/6+1/60\), and the right side in (PT5) is at
most \(23/72<1/3\). This excludes all balanced endpoint
sequences with the indicated total nonatomic mass
and target defect for any fixed integer \(q\ge3\).
The two free summand spectra may be arbitrary subject
to the existing endpoint constraints.
[Pair nonatomic constraints (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

The order \(q\) can be chosen separately for each target before
forming the graph. Indeed
\[
 q\bigl(\log_X\operatorname{rad}(T)
                -\log_X\operatorname{rad}_q(T)\bigr)
                                    \le \log_X T .
\]
The difference in parentheses is at least
\((1/3-f-o(1))L\) at the balanced endpoint. For
\(f<1/3\), the possible \(q\)'s are consequently bounded
in terms of the fixed budget before the height limit.
A pigeonhole choice of one \(q\) preserves the normalized
count. The proof constants may depend on this order in
the inner limit. This also allows \(q=2\) if one uses
the four-target estimate at that order. For instance
\(\Gamma\le 1/2+1/100\) and \(f\le1/100\) satisfy the
strict inequalities for both star sizes, so an adaptive
order \(q\ge2\) is sufficient in that sector.
[Four-target parity bound](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#a-saving-near-the-minimum-nonatomic-mass-with-only-a-target-parity-cap)

On general endpoint profiles \(c_U+c_V\) can be much larger,
or the \(p\)-th-power defects of the three targets in a star
can cost too much. The statement does not give a uniform
bound on the mixed defect dictionary in such stars. That
would be a separate arithmetic input to (PT4).
