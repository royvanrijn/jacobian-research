# Approximate polynomial zeroes and Smith exponents

Let \(v(p)=1\) on \(\mathbb C_p\), with valuation ring \(\mathcal O\).
Let \(F=(F_1,\ldots,F_m)\) be restricted power series of Gauss norm
at most one on the integral unit polydisc, and let \(a\in\mathcal O^m\)
satisfy \(v(F_i(a))\ge e\) for every \(i\). Suppose the Jacobian
\(J=DF(a)\) is invertible over \(\mathbb C_p\). Write
\(s_1\le\cdots\le s_m\) for its Smith valuations over \(\mathcal O\).
They can be obtained by successively choosing an entry of minimum
valuation as pivot and clearing its row and column. Set \(s_0=0\)
when \(m=1\).

If \(e>s_{m-1}+s_m\), there is a common zero \(\zeta\) in the
polydisc with
\[
                      \min_i v(\zeta_i-a_i)\ge e-s_m .                 \tag{AZ1}
\]
In particular \(e>v(\det J)\) suffices, and the precision loss is
at most \(v(\det J)\). The zero may belong to a ramified extension.
For example with \(p\) odd, \(F(x,y)=(x,y^2-p^3)\) at \((0,p^2)\)
has \(e=3\) and Smith valuations \(0,2\). There is a nearby zero
over the algebraic closure and no zero over \(\mathbb Q_p\).
The distinction matters for congruence reconstruction, where a
fixed finite set of local centers can be used even if some centers
are not rational.

## Proof over the algebraic closure

Translate \(a\) to zero. In one variable apply Weierstrass
preparation on the closed unit disc. In the nonzero case \(F(0)\ne0\),
write the series as a scalar times a polynomial with roots in the
disc times a Gauss-norm-one unit whose inverse also has norm one.
The unit contributes an integral logarithmic derivative at zero.
Thus \(v(F'(0)/F(0))\le s_1-e<0\) implies that a reciprocal of
some root has valuation at most \(s_1-e\), by the logarithmic
derivative of the polynomial. That root has valuation at least
\(e-s_1\). If \(F(0)=0\), take the center itself.
[Weierstrass preparation for the Tate algebra, Bosch, §2.2](https://www.math.purdue.edu/~tongliu/seminar/rigid/Bosch.pdf)

Suppose \(m\ge2\) and let \(s=s_1\). Permute rows and columns to
put an entry of valuation \(s\) in the first pivot. Since
\(e>s_{m-1}+s_m\ge2s\), one-variable Hensel applied just to the
first equation, with all remaining coordinates fixed, moves the
first coordinate by a displacement of valuation at least \(e-s\)
to an exact zero of that equation. At the new center every other
residual still has valuation at least \(e\). The linear shift
has valuation at least \(s+(e-s)=e\), and its square has valuation
greater than \(e\). The Jacobian changes by a matrix \(\Delta\)
with entries of valuation at least \(e-s>s_m\). Since the entries
of \(J^{-1}\) have valuations at least \(-s_m\),
\(1+\Delta J^{-1}\) is an integral invertible matrix. Hence the
Smith valuations are preserved and the first pivot still has
valuation \(s\).

Take \(t\in\mathbb C_p\) with \(v(t)=s+\delta\), for a small
positive rational \(\delta\). Rescale all variables about the
new center by \(t\), and divide each equation by \(t\) times
the first pivot. The new first equation has constant zero and
unit first derivative equal to one. Its other linear
coefficients are integral. Every nonlinear coefficient in each
equation has valuation at least \(\delta\), since the original
series and its translations have integral coefficients and
each term of degree at least two gains \(2v(t)\). The remaining
constant terms have valuation at least \(e-2s-\delta\).

The first equation can be solved for the new first coordinate
\(X_1=\Phi(Y_2,\ldots,Y_m)\) in the integral Tate algebra.
Write \(X_1\) as the negative of the other linear terms and the
nonlinear part, and iterate this map in the unit ball of the
Tate algebra in the \(Y\)'s. The nonlinear part is a strict
contraction by the coefficient bound, so this produces an
integral Tate series with \(\Phi(0)=0\). The remaining equations
after substitution still have residuals of valuation at least
\(e'=e-2s-\delta\). Their Jacobian is the Schur complement of
the unit pivot. The Smith valuations of that complement are
\(s_2-s,\ldots,s_m-s\): normalization first subtracts \(s\) from
all valuations, then clearing the unit pivot leaves the other
ones. Choose \(\delta\) so small that the induction threshold
holds, namely \(e'>s_{m-1}+s_m-2s\) when \(m>2\), and
\(e'>s_m-s\) when \(m=2\).

The reduced solution has every coordinate of valuation at least
\(e'-(s_m-s)=e-s_m-s-\delta>0\). The same holds for
\(\Phi(Y)\), since it has integral coefficients and vanishes
at zero. Undoing the scaling gives precision at least \(e-s_m\).
The first shift had valuation at least \(e-s\ge e-s_m\), so
the same precision holds relative to the original center.
This proves (AZ1).

Strictness at the sum of the two largest Smith valuations can be
necessary for the original system. For integers \(1\le b_1\le b_2\)
set \(e=b_1+b_2\) and take
\[
 F_1=p^{b_2}x+p^e,\qquad
 F_2=(x+p^{b_1})y+p^e .
\]
At \((0,0)\) the residual depth is \(e\) and the Smith valuations
are \(b_1,b_2\). The equation \(F_1=0\) forces \(x=-p^{b_1}\),
on which \(F_2\) is the nonzero constant \(p^e\). Thus \(e>s_m\)
alone would fail when \(m=2\). Add identity equations to extend
the example to larger \(m\).

## Uniform recovery in bounded-degree congruence systems

Let \(F_1,\ldots,F_m\in\mathbb Z[\mathbf X]\) have degree at most
\(D\), with \(D,m\) fixed in the inner height limit. At every
sampled rational tuple \(\mathbf x\), assume its denominators
are units modulo \(N\), all \(F_i(\mathbf x)\equiv0\bmod N\),
and \(\det DF(\mathbf x)\ne0\). Suppose the coordinate rational
heights are at most \(H\), and
\[
                           |\operatorname{num}\det DF(\mathbf x)|\le J_0,
                           \qquad 2H^2J_0<N .                       \tag{AZ2}
\]
Then the number of distinct tuples is
\(O_{m,D,\epsilon}(N^\epsilon)\), uniformly in the coefficients.
The bound also holds with \(J_0\) replaced by a cap on the
product of prime losses described below.

The polynomial map \(F\) is dominant and generically finite
because its Jacobian is not identically zero in characteristic
zero. Pick a single \(\mathbf c\in N\mathbb Z^m\) such that
the fiber \(F^{-1}(\mathbf c)\) is finite. Such \(\mathbf c\)
exists by Zariski density and generic finiteness. Its height
plays no role, and the degrees and the derivatives are
unchanged. Bézout bounds its geometric cardinality by \(D^m\).
At each prime \(p^e\Vert N\) with admissible Smith valuations
take a labeled root of \(F-\mathbf c\) from (AZ1), and record
\(l_p=s_m\). At the other primes record \(l_p=e\) and discard
their precision. All these integers are nonnegative and bounded
by \(e\), since the original evaluation is over \(\mathbb Q_p\).
Two tuples with the same retained root label agree to depth \(e-l_p\). At a discarded prime
\(l_p=e\le\sum_i s_i\). At a retained prime
\(l_p\le\sum_i s_i\) as well. Thus the product of the losses
\(\prod p^{l_p}\) divides \(\operatorname{num}\det DF(\mathbf x)\)
in absolute value. One may use a better bound for this product
if one has the Smith valuations separately.

The depth and root labels cost at most
\((1+D^m)^{\omega(N)}\tau(N)=O_{m,D,\epsilon}(N^\epsilon)\).
Within one record every reduced coordinate of height at most \(H\)
has at most one value: its cross difference for any two samples
is an integer of absolute value at most \(2H^2\) divisible by
\(N/\prod p^{l_p}\). This proves the asserted count. For polynomial
systems each center is algebraic since the chosen fiber is
zero-dimensional. Positive-dimensional components of \(F=0\)
therefore pose no problem for this prime-power count.

The proof improves the terminal step in which a determinant
method has already found finitely many auxiliary equations of
bounded degree and an invertible evaluation Jacobian. It does
not supply those equations or that Jacobian at an unrestricted
balanced \(abc\) fibre.
[One grouped root and actual target recovery](ABC_ALGEBRAIC_ROOT_TARGET_OCCUPANCY.md)
· [Classical multivariate Hensel comparison](https://kconrad.math.uconn.edu/blurbs/gradnumthy/multivarhensel.pdf)
