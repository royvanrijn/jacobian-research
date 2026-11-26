# Directional heights for transverse quadratic factors

Fix a finite ordinary two-monomial core, adjacent orders
\(n,N=n+1\), \(n\ge2\), and the signed integer gap \(h\) of the
appended rows. Consider a simple positive physical Gaussian
contact with parameters \(W_*,\lambda_*\), and suppose the
derivative of the phase difference in \(W\) is nonzero.
Use the [exact signed jet notation](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md).
In particular the tail rows have midpoint
\(m_i\), signed integer gap \(d_i\), and moments
\[
 W=\sum_i d_i^2/m_i,\qquad K=\sum_i d_i^3/m_i^2,\qquad
 Z=\sum_i d_i^2/m_i^2 .
\]
Put \(\eta=\max_i |d_i|/m_i\). All asymptotic assertions concern
ordinary shared tails approaching the reference as \(\eta\to0\).
No bound is imposed on the tail row count, gaps or heights.

Suppose the two original moments have a quadratic common factor
whose upper-half-plane root \(\alpha\) continues the indicated
contact. The conjugate root is \(\bar\alpha\).
Define the rational pair
\[
                      T=\frac{\alpha}{\bar\alpha}
                          +\frac{\bar\alpha}{\alpha},\qquad
                      M=S\,\alpha\bar\alpha ,
\]
with a rational squared argument scale \(S\) as follows. Write
\(R_{j,k}\) for the factorial coefficient multiplier of the tail
after its constant factor is removed and set
\(r_{j,k}=R_{j,k}/R_{j,k-1}\). Choose fixed signed integers
\(e_{j,k}\) with \(\sum_{j,k} e_{j,k}=2\) and put
\(S=\prod_{j,k} r_{j,k}^{e_{j,k}}\).
This includes the square of any adjacent tail ratio and
the product of two adjacent ratios in either or both orders.
The extra choice of \(S\) affects the norm comparison
and its arithmetic but not the trace.

## The exact first-order direction

Let \(\zeta_j(W)\) be the Gaussian root of the raw core polynomial
in the tail endpoint-normalized variable. At \(W_*\) write
\(\alpha_j=\partial_W\log\zeta_j\) and \(\kappa_j\) for its cubic
response, as in the signed jet formulas. A star on an evaluation
is suppressed in the following constants. Define
\[
       t=\Im(\alpha_N-\alpha_n)\ne0,\qquad
       \nu=\frac{\Im(\kappa_N-\kappa_n)}{t},\qquad
       L_j=\kappa_j-\nu\alpha_j .
\]
In particular \(\Im L_N=\Im L_n\).
At an exact contact the two clock displacements satisfy
\[
 W-W_*=-\nu K+O(Z),\qquad
 \lambda-\lambda_*=-2nN\,\Re(L_N-L_n)K+O(Z).                   \tag{1}
\]
For the first equation, the imaginary part of the exact contact
equation has terms \(t(W-W_*)+\Im(\kappa_N-\kappa_n)K\) and
error \(O((W-W_*)^2+|W-W_*|\,|K|+Z)\), since
\(B_4\le WZ\) and \(K^2\le WB_4\le W^2 Z\).
Absorbing the small multiple of \(|W-W_*|\) first gives
\(|W-W_*|=O(|K|+Z)\). Substitution makes the error \(O(Z)\).
Taking the real part gives the second equation, using
\(|L_2|\le Z\) for the center correction in the exact formula.
The upper and lower endpoint-normalized roots therefore have
logarithmic changes \(K L_j+O(Z)\) from \(\zeta_j(W_*)\).

For the scale \(S\), use the zero-endpoint log shapes \(F_j,H_j\)
of the signed jet formulas and put
\[
 f=\sum_{j,k}e_{j,k}\bigl(F_j(k)-F_j(k-1)\bigr),\quad
 g=\sum_{j,k}e_{j,k}\bigl(H_j(k)-H_j(k-1)\bigr),\quad
 m=\sum_k e_{N,k}.
\]
Let \(D_j\) be the endpoint scale of the added rows in order \(j\),
so \(\alpha D_j\) is the associated normalized root. Summing
the exact coefficient log expansions gives
\[
 \log(S/D_n^2)
   = Wf+Kg+m\bigl(h\log(N/n)-\lambda/(2nN)\bigr)+O(Z).
\]
Here the clock-difference error consists of \(O(Z)\), since
the signed smaller terms in the center expansion are bounded
by a constant times \(Z\). Consequently
\[
 \begin{aligned}
 T&=T_*+t_1 K+O(Z),&
 t_1&=-4\sin(2\vartheta_*)\,\Im L_n,\\
 M&=M_*+m_1 K+O(Z),&
 m_1&=M_*\bigl[(2-m)\Re L_n+m\Re L_N+g-\nu f\bigr],             \tag{2}
 \end{aligned}
\]
where \(\vartheta_*=\arg\zeta_n(W_*)\),
\(T_*=2\cos(2\vartheta_*)\) and
\[
 M_*=|\zeta_n(W_*)|^2
          q_*^{-2nN f}(N/n)^{mh}v_*^m,\qquad
 q_*=e^{-W_*/(2nN)},\quad v_*=e^{-\lambda_* /(2nN)} .
\]
The reference values and the response coefficients in (2) are
algebraic. The [finite physical exception theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
gives algebraic \(q_*,v_*\) for the finite core. The exponent
\(2nNf\) is integral, since each \(F_j(k)\) has denominator
dividing \(2j\). The remaining constants come from simple
algebraic roots and fixed polynomial derivatives.
All errors in (1) and (2) are uniform using the isolating
disks for the two chosen reference roots.

## Height bounds and rational alignment

Let \(D\) be the least common denominator of \(T,M\),
and put \(\rho=|K|+Z\). The primitive vector \(D(1,T,M)\)
has projective height comparable to \(D\) near the reference.
Write \(x_*=(1,T_*,M_*)\) and \(v=(0,t_1,m_1)\).

**Theorem.** Suppose \(1,T_*,M_*\) are linearly independent
over \(\mathbb Q\). For every \(\epsilon>0\), every sufficiently
close quadratic contact as above satisfies
\[
       D^{\,3+\epsilon}
       \begin{cases}
          Z\rho,&v\ne0,\\
          Z^2,&v=0
       \end{cases}
                    \ \ge\ c_\epsilon>0 .                    \tag{3}
\]
If the rational rank of \(x_*\) is two, let \(\mathcal L\)
be its unique rational two-dimensional linear space.
The same bound holds for points outside \(\mathcal L\),
after making the neighborhood smaller. For points in
\(\mathcal L\), one has \(D^{2+\epsilon}\rho\gg_\epsilon1\).
If \(v\notin\mathcal L\), or if \(v=0\), the latter bound
improves to \(D^{2+\epsilon}Z\gg_\epsilon1\).

For the proof in projective dimension two, take a linear form
with algebraic coefficients annihilating \(x_*\) and \(v\),
and a second independent one annihilating \(x_*\). Use the
first coordinate for the third form. If \(v=0\), choose
any independent basis of the annihilator of \(x_*\) for
the first two. On the primitive vector their normalized
magnitudes are \(O(Z),O(\rho)\), or \(O(Z),O(Z)\) if
\(v=0\). Apply the [projective Subspace Theorem, equation (1.1)](https://arxiv.org/pdf/math/0408381)
to these three independent forms at the real place.
The violating rational points lie in finitely many proper
rational subspaces, taking the rational span of the rational
points of any exceptional subspace defined over an extension.
Every subspace whose real extension misses \(x_*\) can be
discarded in a sufficiently small neighborhood. Rational rank
three excludes all such exceptions. At rank two the only one
which can contain \(x_*\) is \(\mathcal L\).

On \(\mathcal L\), a rational linear form normal to it
evaluated on (2) gives \(K=O(Z)\) when its value on \(v\)
is nonzero. Thus the effective distance to the reference is
\(O(Z)\) when \(v\notin\mathcal L\) or \(v=0\), and
\(O(\rho)\) otherwise. The Subspace Theorem in projective
dimension one gives the corresponding lower bound. Its
proper rational subspaces miss \(x_*\) because the rational
rank is two. A rational identification of \(\mathcal L\)
with \(\mathbb Q^2\) changes the projective height by only
a bounded factor. This proves all cases of (3).

If \(x_*\) is rational, away from exact attainment the
elementary numerator estimate gives \(D\rho\gg1\), improved
to \(DZ\gg1\) when \(v=0\). If \(v\ne0\) and
\(\operatorname{span}(x_*,v)\) is rational, either the
actual point lies exactly in this span or \(DZ\gg1\).
Indeed its rational normal has integer coefficients after
a fixed clearing of denominators and evaluates to \(O(DZ)\)
on \(D(1,T,M)\). These elementary alternatives are effective.

## Direct local budgets for the denominators

At a prime, let \(s_1,s_2\) be the physical root valuations
of the possible quadratic, in the original common variable.
The exact denominator of \(T\) has local exponent
\(|s_1-s_2|\), by the [reciprocal trace rule](MULTIVARIATE_TRACE_HERMITE_JET_HEIGHT.md).
One has \(v_p(M)=s_1+s_2+v_p(S)\).
In particular the exponent of the simultaneous denominator
\(D\) is exactly
\[
             \max\!\left\{|s_1-s_2|,
                                (-s_1-s_2-v_p(S))_+\right\}.  \tag{4}
\]
This is valid also for \(T=0\), where the roots have equal
valuation. The p-adic budgets for \(S\) are the exact adjacent
factorial-block ratios used in its definition. The two
[parent face tests](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md)
let one maximize (4) over only the admissible slope allocations
at each prime. The product of these local caps can upper-bound
\(D\) in (3), with a separate accounting for points on
\(\mathcal L\) if the rank is two.

The norm alone also furnishes a distinct test.
When \(M_*\) is irrational, its reduced denominator \(D_M\)
satisfies
\[
 D_M^{\,2+\epsilon}
  \begin{cases}
     \rho,&m_1\ne0,\\
     Z,&m_1=0
  \end{cases}
                  \ \gg_\epsilon 1,                          \tag{5}
\]
by the same projective-dimension-one argument. Its prime
exponent is only \((-s_1-s_2-v_p(S))_+\).
Thus local allocations with \(s_1+s_2+v_p(S)\ge0\)
contribute nothing to this norm denominator, even if they
give arbitrarily large trace-denominator contributions.
Changing \(S\) changes this test while keeping the physical
common-factor problem fixed. When \(T_*\ne0\), one can also
project to \((1,T/M)\), whose denominator exponent is
\((s_1+s_2+v_p(S)-v_p(T))_+\). For \(T/M\) the first-order
coefficient is \((t_1 M_*-T_*m_1)/M_*^2\); its vanishing
gives the same \(Z\) improvement in (5) whenever the limit
is irrational.

## Several rational norm scales

Choose \(b\ge1\) fixed scales \(S_1,\ldots,S_b\) and form the
rational vector \((1,T,M_1,\ldots,M_b)\), where \(M_i=S_i
\alpha\bar\alpha\). Each coordinate has expansion (2)
with the same \(K\). If the \(b+2\) reference coordinates
are linearly independent over \(\mathbb Q\), repeat the
projective argument with \(b\) independent forms annihilating
the reference and its direction, one further form annihilating
the reference, and the first coordinate. For a nonzero
direction this gives
\[
                    D^{b+2+\epsilon}\rho Z^b\gg_\epsilon1 .
\]
For zero direction use \(D^{b+2+\epsilon}Z^{b+1}\gg_\epsilon1\).
The prime exponent of \(D\) is the maximum of
\(|s_1-s_2|\) and all \(b\) norm contributions from (4).
Thus the several norms can share a denominator budget
instead of requiring a sum of separate budgets. This
assertion uses full rational rank of the enlarged reference
vector.

## Rational angle coordinates for bounded quadratic fields

If the possible root fields belong to a fixed finite list,
select one of them and let \(\Delta<0\) be its discriminant.
A rational signed angle chart valid at every nonreal root is
\[
 J=\sqrt{-\Delta}\cot(\arg\alpha).
\]
Indeed \((\alpha-\bar\alpha)/\sqrt{\Delta}\) and
\(\alpha+\bar\alpha\) are rational, and their quotient in the
opposite order is \(J\) with the upper-half-plane square root
convention. Its first-order response is
\[
                  J_1=-\sqrt{-\Delta}\,
                    \csc^2(\vartheta_*)\,\Im L_n .
\]
The derivative is algebraic. Thus all the directional
height and alignment tests hold for \((1,J,M)\) with
\(v=(0,J_1,m_1)\), using the rational rank and primitive
denominator of this new reference selection.

The relation \(J^2=(-\Delta)(T+2)/(2-T)\) gives the exact
denominator exponent, for \(J\ne0\),
\[
       v_p(D_J)=
         \tfrac12\bigl(v_p(T-2)-v_p(T+2)-v_p(\Delta)\bigr)_+ .  \tag{6}
\]
When \(J=0\) its denominator is one. Every prime with
\(v_p(T)<0\) gives zero in (6), because both \(T-2,T+2\)
then have valuation \(v_p(T)\). Such a prime also splits
in the root field, since the two roots have unequal
valuations. At odd primes the numerator in
\(2-T=(2D_T-A_T)/D_T\), with \(T=A_T/D_T\) in lowest
terms, is the only possible contribution to (6).
The inverse chart \(1/J\), usable when \(T_*\ne-2\),
has first-order coefficient \(-J_1/J_*^2\), and its
denominator is controlled by the numerator of \(T+2\)
and by \(\Delta\). Again the trace denominator primes
give no contribution, since \(\Delta\) is a unit at
split primes. These charts combine the scaled norm budget
with numerator or discriminant budgets rather than paying
for all the trace spread. The
[joint parity and discriminant budget](MULTIVARIATE_QUADRATIC_JOINT_PARITY_DISCRIMINANT_BUDGET.md)
supplies these numerator caps from the two parent arrays.
The charts are compatible with the
[fixed-field conic parameter](MULTIVARIATE_TAME_FIELD_DISCRIMINANT_SQUARECLASSES.md#quadratic-contact-denominator-growth-with-bounded-fields).

These bounds turn the angular and modulus equations of a
transverse quadratic contact into an anisotropic approximation
requirement with two exact factorial-dependent scales
\(\rho\) and \(Z\). The universal bounds are
\(Z\le W\eta\) and \(|K|\le W\eta\).
Under the additional condition \(Z=O(\eta^2)\), the
full-rank bound for one norm and \(v\ne0\) gives
\(D\gg_\epsilon\eta^{-3/(3+\epsilon)}\), using the different
orders of the two errors. The extra condition holds, for
example, for a bounded number of tail rows.
Profiles with \(|K|=O(Z)\) use two errors of order \(Z\).
The Subspace Theorem bounds are qualitative. Uniform upper
bounds for the local budgets, or a substitute separating
prime, are still needed for all ordinary diffuse tails.

