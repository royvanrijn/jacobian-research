# Persistence of an interior Gaussian contact under fixed ordinary rows

Fix an ordinary two-monomial core \(C\) and adjacent orders
\(n,N=n+1\), with \(n\ge2\). Use the physical reference of the
[finite diffuse exception theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts):
\[
 P_j(\tau,z)=
 \sum_{k=0}^j [z^k]Q_j^C(z)
             e^{-\tau k(j-k)/(2j)}z^k,\qquad j=n,N .
\]
At signed tail gap \(h\in\mathbb Z\), compare \(P_n(\tau,z)\)
and \(P_N(\tau,\Lambda^h e^{-\lambda/(2nN)}z)\), where
\(\Lambda=N/n\) and \(|\lambda|\le\tau\). Here \(Q_j^C\) is the
exact factorial moment of the fixed core in the original variable.

Suppose these reference polynomials have a common root at
\((\tau_*,\lambda_*)\) with \(0<\tau_*\) and
\(|\lambda_*|<\tau_*\). Suppose the contacting roots are simple.
Take the lower and upper nonreal root branches \(z(\tau),w(\tau)\)
near that contact before the positive upper rescaling. Choose
the continuous angular difference
\(\theta(\tau)=\arg w(\tau)-\arg z(\tau)\) with value zero at
\(\tau_*\). Assume this difference changes sign there. The
condition includes a transverse crossing and also an odd-order
contact of higher order. It is an effective test at the
algebraic exceptional clock, by the local calculation in the
[finite diffuse exception theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#a-local-bound-when-the-contacting-roots-are-simple).

**Theorem.** Fix an open neighborhood of \((\tau_*,\lambda_*)\)
and a threshold \(\rho>0\). Every sufficiently small-curvature
finite collection \(A\) of ordinary shared rows with relative gaps
less than \(\rho\) can be retained exactly while adding other
ordinary shared rows such that:

- the total added signed gap is \(h\), and the two diffuse
  parameters of the full added profile are in the prescribed
  neighborhood;
- the additional rows beyond \(A\) can have arbitrarily small
  relative gaps;
- the actual two moments have roots in the same physical variable
  whose normalized distance is arbitrarily small.

For the last property the common normalization is by the lower
tail endpoint scale, and both roots stay bounded away from zero.
The normalization in the first two properties uses the height and
signed gap of each ordinary row, as in the finite diffuse exception
theorem. If \(A\) contains an unequal row, the maximum relative gap
of the full tail can be kept exactly equal to that of \(A\).

The conclusion still holds after prescribing arbitrary valuations
and unit residues of both actual constant-normalized coefficient
arrays at finitely many primes, to any finite precision.
Thus a sign-changing interior reference contact can obstruct a
uniform positive separation margin for ordinary tails, even at a
fixed nonzero maximum relative gap. The theorem gives arbitrarily
close pairs of roots. Their exact common-root status remains open.

## Keeping the fixed ordinary rows exact

Write \(S_A=\sum_A d_i^2/M_i\), where \(M_i=\max(a_i,b_i)\)
and \(d_i=b_i-a_i\), and put \(h_A=\sum_A d_i\).
For \(j=n,N\), let \(T_j^A\) be the positive endpoint scale
of the tail \(A\) alone. Strip its constant coefficient and
normalize its leading coefficient to one. Let \(a_{j,k}\) be
the resulting coefficient multiplier after the binomial factor
is removed. These are exact real algebraic numbers.

The
[factorial chord estimate](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size)
and the elementary signed center bound give
\[
 a_{j,k}=1+O_{n,N}(S_A),\qquad
 E_A:=\log(T_N^A/T_n^A)-h_A\log\Lambda=O_{n,N}(S_A).
\]
These bounds are uniform in the number of added rows and their
signed gap. For the chord bound one may take \(S_A<1/2\),
which implies \(d_i/M_i\) in absolute value is at most \(S_A\)
and hence \(M_i/2\le\min(a_i,b_i)\).

Now replace \(P_j(\tau,z)\) by the exact anchored reference
\(P_j^A(\tau,z)\) with \(k\)-th coefficient multiplied by
\(a_{j,k}\). Choose a short closed interval
\([\tau_-,\tau_+]\) around \(\tau_*\) on which the contacting
root branches stay simple and nonreal and the endpoint angular
differences have opposite signs. The branches and signs persist
for \(S_A\) sufficiently small. This follows by isolating the
two simple branches in small root disks and applying the implicit
function theorem over the compact interval.
For some \(\tau_A\) in the interior of that interval, the
perturbed branches \(z_A,w_A\) have equal arguments.

The remaining signed parameter needed for a physical contact is
\[
                 \lambda_A=2nN\left(h\log\Lambda+E_A
                              -\log|w_A/z_A|\right).
\]
The sign follows from
\(w_A=\Lambda^h e^{E_A}e^{-\lambda_A/(2nN)}z_A\).
By shortening the interval before fixing the smallness bound
on \(A\), the pair \((\tau_A,\lambda_A)\) stays close to
\((\tau_*,\lambda_*)\) and still satisfies
\(|\lambda_A|<\tau_A\). This permits replacing the remaining
Gaussian tail by ordinary rows without changing the fixed
rows in \(A\).

## Ordinary approximation and finite local data

Use the
[attainability construction](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#attainability-of-the-limiting-parameters)
for \((\tau_A,\lambda_A)\) with prescribed gap \(h-h_A\).
Concretely, for a large integer \(L\) take
\(\lfloor (L+1)(\tau_A+\lambda_A)/2\rfloor\)
increasing and \(\lfloor (L+1)(\tau_A-\lambda_A)/2\rfloor\)
decreasing unit-gap rows with smaller exponent \(L\).
Correct their signed gap by rows of the required orientation
at a larger height that makes the correction's curvature
and signed parameter tend to zero. This gives exact ordinary
tails with the required signed gap and relative gaps tending
to zero. Including \(A\) adds \(S_A\) to their variance
parameter and adds at most \(S_A\) to the magnitude of the
signed parameter. The chosen neighborhood has room for both.

The
[uniform physical approximation](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#uniform-approximation-of-the-physical-pair)
for these ordinary rows makes their normalized arrays and center
converge to the anchored reference at
\((\tau_A,\lambda_A,h)\). The limiting roots are equal in the
shared physical scale and have nonzero modulus. Continuity of
polynomial roots now gives the claimed small distance.
The maximum relative gap remains the largest one of \(A\)
once the new rows are sufficiently diffuse.

For the local assertions, append the balanced rows furnished by the
[independent real and finite-prime approximation theorem](MULTIVARIATE_REAL_LOCAL_COEFFICIENT_APPROXIMATION.md).
For each finite ordinary profile just constructed those rows can
be given any prescribed local coefficient multipliers while their
real multipliers tend to one. Their signed total gap is zero
and both diffuse parameters tend to zero. They can also be made
so diffuse that the retained rows and the prescribed maximum
relative gap are unchanged. The common lower endpoint normalization
is continuous in the positive endpoints, so the approximating
root distances remain arbitrarily small.

This mechanism applies at an interior sign-changing Gaussian
contact whenever such a contact exists for the exact core.
All its hypotheses can be checked separately at the finitely many
algebraic exceptional clocks of a fixed core.

## Exact persistence with real bases and integer row counts

Under the same sign-changing hypothesis, every sufficiently
large integer \(t\) gives an exact simple-root contact if the
two new unit-gap row bases may be real. The row multiplicities,
the signed gaps of individual rows, and their total signed gap
remain integers. The two real bases tend to infinity and the
diffuse parameters tend to the specified interior contact.
This is also true with any fixed small collection \(A\) from
the theorem retained exactly.

First take \(A\) empty. For prospective parameters
\((\tau,\lambda)\) near the interior contact put
\(S_\pm=(\tau\pm\lambda)/2>0\), and set
\[
 M_+=t+\max(h,0),\qquad M_-=t+\max(-h,0),\qquad
                 m_\pm=M_\pm/S_\pm .
\]
Use \(M_+\) increasing unit rows with midpoint \(m_+\)
and \(M_-\) decreasing unit rows with midpoint \(m_-\).
Their smaller exponents are \(m_\pm-1/2\), and their total
signed gap is exactly \(h\). A real base here means that the
moment of a real exponent \(x\ge0\) is \(\Gamma(x+1)\).
At unit gap the needed ratios are finite rising products,
so the exact coefficients use the same expression as for
ordinary unit rows.

The dependence on \(\upsilon=1/t\) extends analytically through
\(\upsilon=0\) after endpoint normalization of each tail.
For an increasing row, its normalized log multiplier at index
\(k\) is
\[
 \sum_{l=1}^k\log\!\left(1+\frac{l-j/2}{jm}\right)
       -\frac{k}{j}\sum_{l=1}^j
                     \log\!\left(1+\frac{l-j/2}{jm}\right).
\]
For a decreasing row reflect \(k\) to \(j-k\).
The linear term in \(1/m\) in either case is
\(-k(j-k)/(2jm)\).
For each group, \(m_\pm^{-1}=S_\pm\upsilon/
(1+\max(\pm h,0)\upsilon)\) is analytic and zero at
\(\upsilon=0\). Multiplication of the normalized log by
\(M_\pm\) therefore preserves analyticity, with total limit
\(-\tau k(j-k)/(2j)\).

The same property holds for the log relative argument scale.
For one increasing row, subtract \(\log\Lambda\) from the
difference of its two endpoint logs. The result is the
difference at \(j=N,n\) of
\[
              \frac1j\sum_{l=1}^j
                    \log\!\left(1+\frac{l-j/2}{jm}\right).
\]
This difference is analytic in \(1/m\) with linear term
\(-1/(2nN m)\). The decreasing row contributes its negative.
Adding the row groups gives an analytic log relative scale
with limit \(h\log\Lambda-\lambda/(2nN)\) at \(\upsilon=0\).
The term \(h\log\Lambda\) uses the exact difference of the
two multiplicities.

Continue the two simple roots analytically in
\((\tau,\lambda,\upsilon)\) near the reference, and form the
difference between the logarithm of their upper-to-lower
ratio and the log relative argument scale. Choose the log
branch of the root ratio real at the reference.
An exact physical contact is the vanishing of this complex
difference. At \(\upsilon=0\), the derivative of its real
part in \(\lambda\) is \(1/(2nN)\). The real implicit function
theorem therefore solves the modulus equation with an analytic
\(\lambda=\lambda(\tau,\upsilon)\) on a fixed small rectangle.
This function stays strictly inside the bias cone.
At \(\upsilon=0\) the imaginary equation is precisely the
sign-changing reference angle difference. Choose two endpoints
of a small strength interval with opposite signs and with
the reference contact as its only zero. The signs persist
for sufficiently small \(\upsilon\). The intermediate value
theorem then gives an exact phase equality on the modulus
slice, for every large integer \(t=1/\upsilon\).
Every such solution in the interval tends to the reference
strength, by uniform convergence and uniqueness of its zero
there. The bias parameter tends to \(\lambda_*\).
If the crossing is transverse, the full two-dimensional
implicit function theorem makes the two parameters analytic
in \(\upsilon\).

For a fixed small collection \(A\), apply the same construction
to its anchored reference, using the new signed gap \(h-h_A\)
in the two multiplicities and retaining \(T_N^A/T_n^A\) in the
log relative scale. The opposite endpoint signs in the proof
above give a sign-changing zero of the anchored angular
difference. Its roots are still simple, and the nearby bias
is interior. The same modulus and phase argument gives an
exact contact with \(A\) retained.

Thus continuous estimates that allow real bases while keeping
the exact unit-gap factor blocks, integer multiplicities and
integer signed gap must admit contacts near such an interior
exception. An exclusion for ordinary profiles there needs
information about the integral smaller exponents or an
independent exclusion of the reference contact.
