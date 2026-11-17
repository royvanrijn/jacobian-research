# A parity-sensitive Newton obstruction for mixed finite operators

**QPN1.** The noncentral coefficient in [QPC1](UNIVERSAL_P_POWER_OBSTRUCTION.md)
survives an infinite class of genuinely mixed polynomial corrections. The
hypotheses permit arbitrarily high generic momentum order. They do not permit
arbitrary mixed coefficients, negative powers of the position parameter, or a
normalization of every Long lift to the stated form.

## The exact polynomial support region

Let K be a characteristic-zero field, h a nonzero element of K, and
D=K(v)[a;h delta], where delta=v^-1 partial_v. All displayed operators are
finite normal-ordered sums with coefficients on the left. Suppose n>=2 and
\[
 S=g a^n+\mu v+\sum_{(j,e)\in A}c_{j,e}v^e a^j,\qquad
 g\mu\ne0,\quad c_{j,e}\in K,\quad j,e\ge0,                    \tag{1}
\]
where A is finite and excludes the two displayed monomials. Impose:
\[
 \begin{array}{ll}
 e=2k:&j-k<n,\\
 e=2k+1:&j<k.
 \end{array}                                                \tag{2}
\]
Then S has no finite scalar T with [S,T]=kappa for a nonzero
kappa in K, on any finite algebraic coefficient cover of K(v).

The constants may themselves be Laurent series in a deformation parameter.
There is no hypothesis that g is a valuation unit. The theorem therefore
includes h-small highest momentum coefficients and all finite partner orders.

For n=2 the allowed corrections include v^2 a^2, v^3, v^5 a,
and v^(2k) a^(k+1) for every k>=1. The latter family has no
common upper bound on momentum order as k varies. Lower constant-coefficient
powers a^j, j<n, are also allowed. The odd-v restriction excludes v^3 a:
it lies exactly on the equality boundary.

## The same coefficient cannot be changed inside this region

Reduce the finitely many constants and rational coefficients of a
hypothetical mate at any sufficiently large good prime p>2n+1, as in
QPC1. In the reduced finite Ore ring, put
\[
 J=(n-1)p-(2n+1).
\]
The exact coefficient remains
\[
 [v^{-2p}a^J]S^p
       =\frac{n^2-1}{24n}\,g^{p-2}\mu^2h^{p+1}\ne0.          \tag{3}
\]
Here and below the fraction is evaluated in the good residue field.
The source and all arithmetic exclusions of QPC1's coefficient formula
are retained; the new step is the complete selection of words in the
presence of (2).

For a normal monomial v^e a^j define
\[
             w(j,e)=j-e/2,\qquad d(j,e)=n-w(j,e).
\]
An Ore contraction lowers a order by one and v exponent by two.
Thus w is unchanged by a contraction and is additive under multiplication.
Every contribution to the target monomial in (3) has
\[
       \sum_{\text{p factors}}d(j,e)=
                np-(J+p)=2n+1.                            \tag{4}
\]
The top factor a^n has deficit zero. Every other even-v factor
has strictly positive integer deficit. Every odd-v factor has
deficit at least n+1/2, with equality only for the displayed mu*v;
the strict odd inequality in (2) gives deficit greater than n+1/2
for every other odd-v factor.

A word containing no odd-v factors belongs to K[v^2][a;h delta].
That subring is closed under multiplication, since
delta((v^2)^k)=2k(v^2)^(k-1) for k>=0. Its normal order has no
negative v powers, so it cannot contribute to v^-2p.

Every nonzero contraction preserves the parity of the total v exponent.
Since the target exponent -2p is even, the word must contain an even
number of odd-v factors. It consequently has at least two.
Their combined deficit is at least 2n+1. Equation (4), strictness
in (2), and positivity of every remaining non-top deficit show
that there are exactly two odd-v factors, both mu*v, and all p-2
other factors are g*a^n.

Thus precisely the same two-v words as in QPC1 contribute. Its
finite-field summation proves (3), without a new combinatorial
identity or a bound on the total number or degree of permitted corrections.
No cancellation from another support point is possible.

## Finite-characteristic necessity and descent

Over a good residue field of characteristic p, the generic central
localization of K(v)[a;h delta] has central degree p. In a splitting
representation a relation [S,T]=kappa becomes a p-by-p matrix
commutator. Conjugacy by I+epsilon*T over epsilon^2=0 makes the
characteristic polynomial of S translation invariant to first order.
It is therefore z^p+c, and Cayley--Hamilton forces S^p to be central.

The target order J in (3) is positive and is not divisible by p.
The center is K(v^p)[a^p], so (3) contradicts that necessity.
All facts used in this paragraph and the existence of arbitrarily
large good reductions are justified in QPC1. Only finitely many
constant-field elements occur, even if they are scalar Laurent series;
the argument does not evaluate a formal series at h=1.

For a finite algebraic coefficient extension E/K(v), normalized field
trace, applied to the normal coefficients of T, is a bimodule retraction
D(E)->D(K(v)). It commutes with delta and fixes S and kappa.
A mate on the cover would consequently give a mate over K(v).
This proves the asserted arbitrary finite-cover result. The same
normalized matrix trace excludes mates for the fixed scalar image
S*I_m in every finite matrix size. A genuinely nonscalar first image
does not admit that reduction.

## The strict boundary and the Long scope

The selection argument requires strict inequalities in (2).
At n=2 and p=7, the finite operator
\[
                    S=a^2+v+3v^3a
\]
has zero coefficient at v^-14*a^2 in S^7, whereas a^2+v has
coefficient 4. The new term has exactly the same deficit as v.
This is a countercontrol to preservation of this one witness after
weakening the odd inequality; it supplies neither a mate nor a
characteristic-zero counterexample.

For contrast, P=a^2+v^2/2, Q=-a is an exact finite canonical
pair. It has no nonzero v term and is outside (1).
These controls keep both the nonzero mu condition and the strict
support region visible.

A bounded exact-R scalar Long proposal whose localized generic S
satisfies (1)--(2) is excluded, assuming its actual global admission
and full partner relations. No theorem puts every Long representative
into this support region. In particular, a generic correction with
an odd-v momentum monomial on or above the boundary, a higher
even-v weight, a negative-v Laurent term, a deformed first image,
or a genuinely nonscalar carrier remains outside the result.

[Exact controls](check_mixed_p_power_obstruction.py) recompute full
sparse Ore powers for representative mixed operators, including one
whose generic momentum order exceeds n, compare the target coefficient
with (3), check the support gate rejects the equality case, and verify
the stated equality-boundary and positive-pair controls. Finite checks
support implementation; the all-support, all-order and arbitrary-cover
claims use the written word-selection proof and QPC1.

No finite-left Long carrier, bounded full Long lift, DC1/DC2 solution,
all-representative nonlifting, formal verification, external review,
independent whole-proof replay or novelty priority is claimed.
