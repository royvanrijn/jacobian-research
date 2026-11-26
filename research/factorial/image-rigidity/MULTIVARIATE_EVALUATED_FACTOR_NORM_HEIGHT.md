# Evaluated norm bounds for an isolated common quadratic

Fix a nonreal algebraic reference root \(\alpha\) of degree \(D>2\).
Use a common rational argument normalization of the two moments,
as in the [simultaneous height comparison](MULTIVARIATE_SIMULTANEOUS_CONTACT_HEIGHT.md).
In this variable multiply the reference root by a fixed positive
integer if needed to make it integral, and apply the same argument
normalization to both actual parents. Write \(T\in\mathbb Z[z]\)
for the monic minimal polynomial of the resulting \(\alpha\).
Put \(K=\mathbb Q(\alpha)\). Let \(S\) be any fixed finite set of
rational primes. For a positive integer \(m\) put
\(m_{S^c}=\prod_{p\notin S}p^{v_p(m)}\).

Take primitive integer parent polynomials \(F,G\) in this chart.
Suppose their rational gcd has primitive integer form
\[
             H(z)=x_0+x_1z+x_2z^2
\]
with a nonreal root \(\beta\) approaching \(\alpha\). Write
\(X=\max_i|x_i|\) and
\[
       C=\gcd\bigl(|\operatorname{Res}(T,F)|,\,
                    |\operatorname{Res}(T,G)|\bigr).
\]
The convention \(\gcd(0,m)=m\) applies. Both resultants cannot
vanish because the common gcd has degree two and \(D>2\).

**Theorem.** Put \(\kappa=4\) if \(D\ge5\). One may put
\(\kappa=3\) if \(D\ge4\) and either \(K\) has no subfield of
index two or the coordinates of
\((|\alpha|^2,-2\Re\alpha,1)\) are linearly independent over
\(\mathbb Q\). Then for every \(0<\delta<D-\kappa\),
\[
 \begin{split}
 C_{S^c}&\gg X^{D-\kappa-\delta},\\
 C_{S^c}&\gg |\beta-\alpha|^{-
                    2(D-\kappa-\delta)/(\kappa+\delta)}.
                                                               \tag{1}
 \end{split}
\]
The constants may depend on the reference root, \(S\), and
\(\delta\). The second condition for \(\kappa=3\) uses a
neighborhood of the reference quadratic. In particular, if
the two evaluated parent resultants have a gcd supported on
fixed primes, this excludes accumulating common quadratics.
More generally, an approach rate \(|\beta-\alpha|=O(\eta^\gamma)\)
requires prime-to-\(S\) gcd growth at least the power of
\(1/\eta\) given by multiplying the exponent in the second
line by \(\gamma\).

The test concerns primes dividing the two evaluated norms.
It can be used when the parent coefficient arrays have
unrestricted prime support. A uniform upper bound on the
evaluated gcd for arbitrary remaining factorial tails is
a separate obligation.

## A norm-form inequality in three variables

Let
\[
            U=|\operatorname{Norm}_{K/\mathbb Q}(H(\alpha))|.
\]
This is a positive integer because \(\alpha\) is integral
and \(2<D\). The proof gives \(U_{S^c}\gg X^{D-\kappa-\delta}\)
for every primitive coefficient vector in the universal
\(\kappa=4\) case and in the subfield case of \(\kappa=3\).
In the rational independence case it gives the bound in a
neighborhood of the target line.

Work in a splitting field \(E\) of \(T\), and let \(S_E\)
consist of the Archimedean places and all places over \(S\).
For the \(D\) distinct roots \(\alpha_i\) of \(T\) use the
linear forms
\[
                    L_i(x)=x_0+x_1\alpha_i+x_2\alpha_i^2.
\]
At each place the absolute values are normalized for the
product formula. Since \(\prod_i L_i(x)\) is the rational
norm of \(H(\alpha)\), one has
\[
          \prod_{w\in S_E}\prod_{i=1}^D |L_i(x)|_w
                         =U_{S^c}.
\]
The primitive rational vector has norm one at every finite
place, and \(\prod_{w\mid\infty}\|x\|_w=X\).

Any three distinct evaluations are independent by their
Vandermonde determinant. At \(w\) choose three of smallest
absolute value at \(x\). The third of these is bounded below
by \(c_w\|x\|_w\), since the inverse of the matrix of the
three forms is bounded and there are only finitely many
choices. The product of the three selected evaluations over
\(S_E\) is consequently at most a constant times
\(U_{S^c}/X^{D-3}\). If \(U_{S^c}\le X^{D-3-\delta}\)
holds along an unbounded sequence, these products are at
most \(X^{-\delta/2}\) eventually. There are only finitely
many choices of systems, one at each place. The projective
[Subspace Theorem, equation (1.1)](https://arxiv.org/pdf/math/0408381)
puts the points in finitely many proper linear subspaces.
The rational points in each such subspace lie in a proper
rational subspace, since independent rational vectors stay
independent over \(E\).

If the limiting coefficient vector has rational rank three,
it belongs to no proper rational subspace. Primitive vectors
projectively approaching it therefore leave all the resulting
subspaces. This proves the rank-three target case.

## Restriction to rational planes

Fix one rational two-dimensional subspace \(W\) from the first
application. A basis of \(W\) gives independent polynomials
\(f,g\in\mathbb Q[z]\) of degree at most two. Neither can vanish
at a conjugate \(\alpha_i\), since \(D>2\). The restricted
evaluations at \(\alpha_i,\alpha_j\) are proportional exactly
when
\[
                      \frac{f(\alpha_i)}{g(\alpha_i)}
                  =   \frac{f(\alpha_j)}{g(\alpha_j)}.
\]
For each value of the ratio the equation \(f-z_0g=0\) is a
nonzero polynomial of degree at most two in its argument.
Hence a proportionality class has at most two members.
Moreover, a repeated ratio at distinct conjugates implies
that \(K\) has an index-two subfield. In the embedding
corresponding to \(\alpha_i\), a conjugation taking \(\alpha_i\)
to \(\alpha_j\) fixes their ratio. The field generated by
this ratio has extension degree at least two in
\(\mathbb Q(\alpha_i)\) and at most two from \(f-z_0g=0\).

At each \(w\), select the smallest evaluation \(L_i\) and an
independent one \(L_j\) among the first three smallest.
Their maximum controls \(\|x\|_w\). Thus all but the first
two evaluations are bounded below by \(c_w\|x\|_w\),
and both of the first two are at least \(|L_i(x)|_w\).
Since \(|L_j(x)|_w\ll_w\|x\|_w\), this gives
\[
        \prod_{k=1}^D |L_k(x)|_w
        \gg_w \|x\|_w^{D-4}\,
                    |L_i(x)L_j(x)|_w^2 .
\]
The estimates are uniform over the finitely many selections
on the fixed plane. Multiplying over \(S_E\), a violation
\(U_{S^c}\le X^{D-4-\delta}\) forces the product of the two
selected forms to be \(O(X^{-\delta/2})\).
Express \(x\) in an integral lattice basis of the saturated
lattice \(W\cap\mathbb Z^3\). Its two coordinates are
primitive, have Archimedean norm comparable to \(X\),
and have norm one at each finite place. The two-form
Subspace Theorem therefore reduces the violating points
to finitely many rational lines. On each rational line
there are only two primitive integer vectors.

If \(K\) has no index-two subfield, choose the first two
evaluations. They are independent. All but the first one
are bounded below by \(c_w\|x\|_w\), giving instead
\[
           \prod_{k=1}^D |L_k(x)|_w
             \gg_w \|x\|_w^{D-2}|L_i(x)L_j(x)|_w.
\]
This excludes violations of the \(D-3-\delta\) bound on the
same planes. Finitely many exceptional primitive vectors
can be absorbed by lowering the positive implicit constant.
This proves all the norm estimates. The proof applies
directly to primitive vectors, including when \(H(\alpha)\)
is a unit or a root of unity.

## Passing back to the factorial parents

Multiplicativity of the integer resultant gives
\[
        U=|\operatorname{Res}(T,H)|\ \mid\ C .
\]
For the Archimedean comparison \(H(z)=x_2(z-\beta)(z-\bar\beta)\),
so the evaluations at \(\alpha,\bar\alpha\) are both
\(O(X|\beta-\alpha|)\). The remaining conjugate evaluations
are \(O(X)\). Hence
\[
            U_{S^c}\le U\ll X^D|\beta-\alpha|^2 .
\]
Together with the norm lower bound this implies
\(X^{\kappa+\delta}|\beta-\alpha|^2\gg1\), and then (1).
For an algebraic positive diffuse contact of a fixed ordinary
core, the [rational chart and algebraic factor selections](MULTIVARIATE_SIMULTANEOUS_CONTACT_HEIGHT.md)
supply the setting above. The norm lower bound uses a
qualitative Subspace Theorem.

## A joint prime-ideal budget

The same conclusion holds with a smaller budget than \(C\).
Enlarge \(S\) to contain the prime divisors of
\(T(0)\operatorname{disc}(T)\). At a rational prime \(p\notin S\)
write \(\mathfrak p\) for the primes of \(\mathcal O_K\) over
\(p\), \(f_{\mathfrak p}\) for their residue degrees, and
\[
  m_{\mathfrak p}
     =\min\bigl(v_{\mathfrak p}(F(\alpha)),v_{\mathfrak p}(G(\alpha))\bigr).
\]
Here the valuations are normalized to take integer values on
principal ideals. If one evaluation vanishes, use infinite
valuation for it. The norm of the common evaluation ideal
\((F(\alpha),G(\alpha))\) gives a budget with \(p\)-exponent
\(\sum_{\mathfrak p\mid p}f_{\mathfrak p}m_{\mathfrak p}\).
This can be smaller than the scalar gcd of norms, which takes
the minimum only after summing over \(\mathfrak p\).

Let \(r_p\le2\) be any upper bound on the number of unit roots
in the proposed common factor at \(p\). One can take two or
use the [Newton and truncated face width tests](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md).
Define
\[
 B_p=\max_{\substack{I\subseteq\{\mathfrak p\mid p\}\\
                    \sum_{\mathfrak p\in I}f_{\mathfrak p}\le r_p}}
                   \ \sum_{\mathfrak p\in I}f_{\mathfrak p}m_{\mathfrak p}.
\]
Both inequalities in (1) hold with
\(\prod_{p\notin S}p^{B_p}\) in place of \(C_{S^c}\).

Indeed \(H(\alpha)\) divides both evaluations in \(\mathcal O_K\).
The primes at which it is nonunit correspond to distinct
irreducible factors of \(\bar T\) dividing \(\bar H\).
The correspondence follows from \(p\nmid\operatorname{disc}(T)\):
the order \(\mathbb Z[\alpha]\) is maximal at \(p\), because its
discriminant is the field discriminant times the square of its
index, and the distinct factors of \(T\bmod p\) give the primes
and their residue degrees. As \(p\nmid T(0)\), these factors have
nonzero constant coefficients. Their degrees count nonzero
residual roots of the primitive factor, bounded by the number of
unit roots using the horizontal Newton width. Thus their total
degree is at most \(r_p\). Adding the valuation upper bounds at
these primes gives \(v_p(U)\le B_p\). This proves the refined
budget.

For a strict local comparison, take \(p=31\) and
\[
       T(z)=\prod_{j=1}^5(z-2j)+62\cdot1000001 .
\]
The polynomial is Eisenstein at two and has five distinct
linear residual factors at \(31\). With
\(F(z)=(z-2)(z-4)\) and \(G(z)=(z-6)(z-8)(z-10)\), the two
evaluations are divisible at disjoint primes above \(31\).
The scalar gcd of their norms has valuation two there, while
the common ideal gives zero. With \(F=T+31\) and
\(G=zT+31^2\), the ideal exponent is five, while the budget with
\(r_p=2\) is two. These examples use polynomial evaluations to
compare the budgets.

## The index-two loss is real

Let \(D=2m\) with \(m\ge2\), put \(\rho=2^{1/m}\), and take
\(\alpha=i\sqrt\rho\). Its polynomial
\(z^{2m}-(-1)^m2\) is Eisenstein at two, and
\(\mathbb Q(\rho)\) has index two in \(K\).
For coprime positive rational convergents \(p/q\) to \(\rho\),
put \(H_q(z)=qz^2+p\). Its nonreal roots approach those of
\(T\) at the chosen embedding and its norm is exactly
\[
        |\operatorname{Norm}_{K/\mathbb Q}(H_q(\alpha))|
                          =(p^m-2q^m)^2 .
\]
One conjugate factor \(p-q\rho\) is \(O(q^{-1})\), while
the other \(m-1\) factors in the norm from
\(\mathbb Q(\rho)\) are \(O(q)\). Thus this norm is
\(O(q^{D-4})\), showing the general norm exponent is
sharp. In degree four a constant-norm sequence is explicit.
Write \(a_q+b_q\sqrt2=(1+\sqrt2)^{2q}\) and take
\(H_q(z)=b_qz^2+a_q\). The coefficients are primitive by
their quadratic norm equation, and the norm of
\(H_q(\alpha)\) is one while its root tends to \(\alpha\).
These fixtures concern the norm condition on quadratics.
