# Effective slope selection and cutoffs from actual Bernstein certificates

This proof establishes the finite extraction step for the
[graph certificate construction](OAKU_GRAPH_CERTIFICATES_AND_ACTUAL_POLYNOMIAL_DECISION.md).
It uses the supplied identities for the actual polynomial powers. The
[selector implementation](bernstein_slope_producer.py) performs the resulting
rational algebra and finite original-power evaluations. The pointwise Gaussian
construction is one way to obtain the inputs to this proof.

## 1. Certificate contract

Let \(K\) have characteristic zero, \(P\) have positive degree \(D\), and
\(T=\lambda(\partial)\), where the nonzero formal symbol \(\lambda\) has
positive order \(r\). Put
\[
  \tau_k=\min\{j\ge0:T^k(P^j)\ne0\},
\]
with value \(\infty\) for an empty set. Suppose a nonzero \(p\in\mathbf Q[X,Y]\)
and integers \(k_0\ge1\), \(J\ge0\) certify the following properties.

1. \(p(k,\tau_k)=0\) at every finite \(\tau_k\) with \(k\ge k_0\).
2. For \(k\ge k_0\), if \(\tau_k\) is finite then \(\tau_{k+1}\) is finite
   and \(0\le\tau_{k+1}-\tau_k\le J\).

If \(\tau_k=\infty\), every subsequent order is infinite because \(T^{k+1}\)
factors through \(T^k\). If \(\tau_k<\infty\), the degree drop gives
\(\tau_k\ge rk/D\). The finite extraction algorithm takes the all-index
properties as proved inputs.

## 2. Rational slopes and the finite set of affine lines

Suppose the sequence remains finite after \(k_0\). Remove repeated and vertical
factors of \(p\). Beyond the finitely many zeros of the leading coefficient and
discriminant in \(Y\), the surviving real branches are disjoint smooth graphs
\(Y=f_i(X)\) on a final real interval. Every such algebraic branch, and every
difference between two branches, is eventually monotone or constant. Indeed,
its derivative is algebraic, and a nonzero algebraic function on the interval
has only finitely many zeros there. To see this last assertion, use a polynomial
equation irreducible over \(\mathbb R(X)\) for the function and clear
denominators. Its specialization at zero in the function variable is a
nonzero polynomial in \(X\), unless the function is identically zero.

Since \(\tau_k=O(k)\), discard branches for which \(|f_i(X)|/X\) tends to
infinity. Every remaining branch has \(f_i(X)=O(X)\), since \(f_i(X)/X\) is also
algebraic and eventually monotone. Its derivative is eventually monotone or
constant by the same reasoning applied to the second derivative. An
unbounded eventually monotone derivative would force the ratio \(f_i(X)/X\)
to be unbounded, by integrating along the final interval. Therefore each
remaining branch has bounded derivative.

Group the remaining branches by bounded difference. Each difference outside
a group tends to infinity in absolute value. The bound on the derivatives
and the step bound on \(\tau_k\) prevent a change between distinct groups at
consecutive large integers. Hence a single group supplies \(\tau_k\) at every
sufficiently large integer \(k\).

Each nonlinear branch is eventually strictly convex or concave because its
nonzero algebraic second derivative has finitely many zeros. Its arc length
through \(X=N\) is \(O(N)\). Consecutive lattice chord vectors in increasing
\(X\) order are distinct because their slopes are strictly ordered. There are
\(O(R^2)\) integer vectors of length at most \(R\), so the total length of
\(L\) distinct nonzero integer vectors is \(\gg L^{3/2}\). The chord lengths
sum to at most the arc length. Thus a nonlinear branch contributes only
\(O(1+N^{2/3})\) lattice points through \(X=N\). An affine line with irrational
slope or irrational intercept contains at most one lattice point. The selected
group supplies \(N-O(1)\) points, so it contains a rational affine line of
slope \(\rho=a/b\), written in lowest terms. The bounded group differences,
degree drop and step bound give
\[
  \tau_k=\rho k+O(1), \qquad r/D\le \rho\le J.
\]
The integers \(b\tau_k-ak\) now belong to a finite set after a fixed index.
For an offset \(c\) in that set, either \(p(X,(aX+c)/b)\) vanishes identically
or it has finitely many roots. Every sufficiently late offset therefore comes
from a line factor \(bY-aX-c\) of \(p\).

## 3. Finite extraction by actual power tests

Factor \(p\) over \(\mathbb Q\) and collect the distinct rational line slopes
in \((0,J]\). For each one \(\rho=a/b\), write
\[
 q(X,Z)=p(X,(aX+Z)/b)=\sum_{i=0}^d X^i q_i(Z),\qquad
 f(Z)=q_d(Z)\ne0.
\]
Let \(\alpha\) be the leading coefficient of \(f\). Take a rational Cauchy
bound \(R\) on the moduli of its roots, or \(R=0\) if \(f\) is constant. Define
\[
 B=\lceil R\rceil+1,\qquad L=a+bJ,\qquad A=B+L+1,\qquad
 E=\sum_{i<d,j\ge0}|[X^i Z^j]q|A^j,\qquad
 K_\rho=\max(1,\lfloor E/|\alpha|\rfloor+1).
\]
The root bound gives \(|f(Z)|\ge|\alpha|\) for \(|Z|\ge B\). For \(X\ge1\)
and \(|Z|\le A\), the lower-\(X\)-degree terms of \(X^{-d}q\) have absolute
sum at most \(E/X\). Consequently no real root of \(q(X,Z)\) lies in either
closed interval \(B\le Z\le A\) or \(-A\le Z\le-B\) at \(X\ge K_\rho\).
For the actual selected slope, \(z_k=b\tau_k-ak\) has jumps of absolute
size at most \(L\). It eventually lies in \(|Z|<B\), since every retained line
offset is a root of \(f\). Steps cannot jump across either of the two
root-free intervals. Thus
\[
                      |b\tau_k-ak|<B
                      \quad(k\ge\max(k_0,K_\rho)).
\]

Choose a common integer \(K\ge k_0,K_\rho\) beyond the roots of the leading
\(Y\)-coefficient of \(p(X,Y)\), large enough that the intervals
\(|y-\rho K|<B/b\) are pairwise disjoint. These are explicit rational
inequalities. List every nonnegative integer root \(j\) of the nonzero
univariate polynomial \(p(K,Y)\), and evaluate the original \(T^K(P^j)\) for
each of them. A test at \(j\) requires only the symbol coefficients through
order \(Dj\). If none survives then \(\tau_K=\infty\). Otherwise the first
survivor is \(\tau_K\). The propagation premise gives an all-finite tail,
and the disjoint intervals identify its unique actual slope.

## 4. Safety, sufficient cutoffs, and critical markers

Write \(h_m=\max\{k:T^k(P^m)\ne0\}\). For a polynomial marker \(Q\), commuting
\(T\) through multiplication by \(Q\) replaces its polynomial coefficients
by derivatives of those coefficients. Constant-coefficient operators commute
with \(T\). Thus the iterated commutators of \(T\) and multiplication by \(Q\)
vanish after \(\deg(Q)+1\) commutators. In the binomial commutator expansion
of \(T^m(QP^m)\), every term contains \(T^{m-i}(P^m)\) with
\(i\le\deg(Q)\). In particular \(m-\deg(Q)>h_m\) implies the marked output
vanishes.

If \(\tau_K=\infty\) then \(h_m\le K-1\) for every \(m\). This gives the
sufficient cutoff \(m>\deg(Q)+K-1\). When the selected slope is
\(\rho=a/b>1\), put
\[
 H_0=\max(K-1,\lceil B/a\rceil),\quad
 C=\lceil a/(a-b)\rceil,\quad E_0=\lceil aH_0/(a-b)\rceil.
\]
For \(h_m\ge K\), \(T^{h_m}(P^m)\ne0\) implies \(\tau_{h_m}\le m\).
The strip bound shows \(h_m\le m/\rho+B/a\). For \(h_m<K\) we have
\(h_m\le m/\rho+H_0\) as well. Thus \(m>C\deg(Q)+E_0\) implies
\(m-\deg(Q)>h_m\) and the marked output vanishes. If all positive pure
outputs vanish, constant markers are already covered. For positive degrees,
a sufficient linear slope is \(K\) in the eventual-annihilation branch
and \(C+E_0\) in the slope-\(\rho>1\) branch.

Conversely, suppose every marker has eventual vanishing. Testing \(Q=1\)
gives \(T^j(P^j)=0\) for all large \(j\), hence \(T^k(P^j)=0\) whenever
\(k\ge j\) and \(j\) is large. All the finitely many smaller \(j\) also die
for large \(k\) by degree drop. For a fixed \(d_0\), testing the markers
\(1,P,\ldots,P^{d_0}\) then gives \(\tau_k>k+d_0\) for all sufficiently
large \(k\). As this holds for every \(d_0\), a finite slope must exceed
\(1\). Hence the finite selector decides qualitative all-marker safety.

If every positive pure output vanishes, \(\tau_k>k\) for \(k\ge1\).
When \(\rho=1\), the eventual line offsets are positive integers. At least
one offset \(d\) occurs infinitely often and the fixed original marker
\(P^d\) fails eventual vanishing. The positive offsets of the line factors
give a finite family containing such a marker.

## 5. A coefficient-query limitation

For a bare oracle returning coefficients of an arbitrary formal symbol, no
algorithm can always halt and distinguish an eventually infinite order
sequence from a finite rational-slope sequence. This already holds for
\(P=x\) in two variables, even when the input is pure and safe. Suppose the
algorithm halts on \(\lambda=\xi_y\), for which every positive order is
infinite. Choose \(M\ge2\) so that it did not query the coefficient of
\(\xi_x^M\). On \(\lambda=\xi_y+\xi_x^M\), the algorithm receives the same
answers, but on \(K[x]\) the operator is \(\partial_x^M\), and
\(\tau_k=Mk\). This input is pure and safe by the degree bound above. The
obstruction concerns exact selection of the alternative under this input
representation.
