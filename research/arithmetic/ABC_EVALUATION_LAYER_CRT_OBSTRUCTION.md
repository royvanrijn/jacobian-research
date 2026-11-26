# A robust CRT obstruction for varying evaluation layers

The finite-field evaluation layers have no bounded-height integer
lift on almost all marked incidences, even with arbitrary label maps.
A finite difference on several levels turns the simultaneous congruences
into an exact integer identity. A stability argument on the cyclic
digit then contradicts injectivity.

## The marked template

Let \(k\ge1\) be fixed and let \(q\) tend to infinity through primes.
Use the digit box \(\mathcal A\), cardinality \(N\), and evaluation maps
\(f_a(t)\) of (DI1). The constant digit ranges over \(\mathbb F_q\).
The linear digit belongs to \(\{0,\ldots,s_1-1\}\) with \(s_1\to\infty\)
and \(s_1\le q\). Any further digits have the prescribed ranges
of (DI1). In particular that linear digit can be truncated.
The two sides are copies of the same digit box, and before excluding
\(f_a+f_b=0\), the level at \(t\) consists of the \(q\) complete
bipartite blocks \(f_a(t)=r,\ f_b(t)=-r\).
There are \(N^2\) such edge-level incidences over all levels.
Pairs with \(f_a+f_b=0\) remove at most \(qN=o(N^2)\) incidences.

Give each level a positive integer modulus \(M_t\). Assume
\[
          \operatorname{lcm}_{t\in I} M_t>2^k X
                         \qquad (|I|=k).
                                                               \tag{EL1}
\]
Then there is a fixed \(\epsilon_k>0\) such that for large \(q\)
no injective map \(A:\mathcal A\to [1,X]\cap\mathbb Z\) and map
\(B:\mathcal A\to\mathbb Z\) realize all but \(\epsilon_k N^2\)
of the marked incidences by the congruence
\(A(a)+B(b)\equiv0\pmod{M_t}\). The same assertion holds
with a specified right side constant in each block.

For the numerical parameters in (DI1), take \(M_t=p_t^e\).
The primes are distinct and at least \(q\), so for
\(L\le\kappa/2\),
\(\log_X(p_t^e)\ge\beta\kappa/2\). It suffices to fix
\(k>2/(\beta\kappa)\) before the height and budget limits.
There is no polynomial, continuity, or affine constraint on
either map.
[The varying evaluation layers and their pair-capacity test](ABC_DEEP_TARGET_LAYER_INCIDENCE_OBSTRUCTION.md#varying-partitions-with-bounded-pair-agreement)

## Stability of integer finite differences

Let \(G\) be a finite abelian group and let \(f:G\to\mathbb Z\).
Write \(\Delta_h f(x)=f(x+h)-f(x)\).
For any \(j\ge1\) and \(0<\alpha<1/2\) there is a
\(\delta_j(\alpha)>0\) independent of \(G\) such that
\[
  \Pr_{x,h_1,\ldots,h_j}
      (\Delta_{h_1}\cdots\Delta_{h_j}f(x)=0)
                    \ge 1-\delta_j(\alpha)
                                                               \tag{EL2}
\]
implies that \(f\) takes one value on at least
\(1-\alpha\) of its domain. The variables in the probability
are independent and uniform.

First recall the elementary majority correction for approximate
homomorphisms. If \(c:G\to\mathbb Z\) fails
\(c(x+y)=c(x)+c(y)\) on at most a fraction \(\rho<1/6\)
of pairs, the distribution of \(D_x(y)=c(x+y)-c(y)\)
has a value \(g(x)\) with frequency at least \(1-2\rho\)
for every \(x\). To see this, for independent uniform \(y,z\)
use the addition tests at \((x+y,z)\) and \((y,z)\).
Except on at most \(2\rho\) of choices they give
\(D_x(y+z)=D_x(y)\). The two arguments \(y+z,y\) are
independent and uniform, so the collision probability of
\(D_x\) bounds its maximum frequency from below.
For fixed \(x,z\), the identity
\(D_{x+z}(y)=D_x(z+y)+D_z(y)\) gives
\(g(x+z)=g(x)+g(z)\) by choosing \(y\) for which all three
values are modal. The union bound allows this since \(6\rho<1\).
When \(c(x)\ne g(x)\), the original addition test fails at
that \(x\) for at least \(1-2\rho\) of \(y\).
Thus \(c\) differs from the exact homomorphism \(g\)
on at most \(\rho/(1-2\rho)\) of \(G\).
This is the majority argument of the homomorphism test.
[Blum, Luby and Rubinfeld, Self-testing/correcting with applications to numerical problems](https://people.csail.mit.edu/ronitt/papers/selftestcorrect.pdf)

For \(j=1\), (EL2) is the collision probability of \(f\)
and one can take \(\delta_1(\alpha)=\alpha\).
Inductively take
\(\delta_j(\alpha)=(\alpha/100)\delta_{j-1}(\alpha/4)\).
At most a fraction \(\alpha/100\) of shifts \(h\) have
\((j-1)\)-difference failure rate above
\(\delta_{j-1}(\alpha/4)\) for \(\Delta_hf\).
For each other shift the induction gives a constant \(c(h)\)
agreeing with \(\Delta_h f\) on at least \(1-\alpha/4\) of \(G\).
Assign arbitrary values of \(c\) at the exceptional shifts.
If \(h,h',h+h'\) are all good, choose an \(x\) where the
three constant predictions in
\[
    \Delta_{h+h'} f(x)=\Delta_h f(x)+\Delta_{h'} f(x+h)
\]
hold together. This gives \(c(h+h')=c(h)+c(h')\).
The fraction of failures in the addition test for \(c\)
is at most \(3\alpha/100\). Majority-correct \(c\) to a
homomorphism \(G\to\mathbb Z\), which is zero because \(G\)
is finite. On all but at most
\(\alpha/100+(3\alpha/100)/(1-6\alpha/100)\) of shifts,
\(\Delta_hf\) is therefore zero on at least \(1-\alpha/4\)
of \(G\). The resulting collision probability of \(f\)
is at least \(1-\alpha\), which proves the induction.

## Applying the stability argument to the layers

Suppose only \(\epsilon N^2\) marked incidences fail.
Including the deleted \(f_a+f_b=0\) pairs, at most
\(2\epsilon N^2\) incidences are absent from the complete
level blocks for large \(q\). By Markov's inequality
all but \(O(\epsilon^{1/2}q)\) levels lose at most
\(O(\epsilon^{1/2}qH^2)\) incidences, where \(H=N/q\)
is each block's side size. At each of these levels,
all but \(O(\epsilon^{1/4}q)\) blocks lose at most
\(O(\epsilon^{1/4}H^2)\) edges. In a remaining block
all but \(O(\epsilon^{1/4}H)\) first-side vertices lie
in its largest first-side connected component.
Indeed, if its first-side size is \(a_*\), the total
edge count is at most \(a_*H\). Let \(G_t\subseteq\mathcal A\)
be the union of the selected first-side components at such
a level \(t\). Its complement has size \(O(\epsilon^{1/4}N)\).
Any two points in \(G_t\) with the same \(f_a(t)\) have
equal first-side labels modulo \(M_t\), by the congruences
along a path in the block.

Keep all digits above the linear digit unchanged. For \(k\)
nonzero pure constant-digit shifts \(h_i\in\mathbb F_q\),
choose \(2k\) distinct good levels \(t_i,s_i\) with
\(s_i-t_i=h_i\) in \(\mathbb F_q\).
For each \(i\) there are at least \(q-O(\epsilon^{1/2}q)\)
choices before removing the \(O(k)\) already used levels.
Choose these levels as functions of the shifts, independent of
the base point. In the first two digits decompose
\[
           (h_i,0)=u_i+v_i,\qquad
           u_i=(-t_i,1),\quad v_i=(s_i,-1).
\]
Translation by \(u_i\) preserves evaluation at \(t_i\),
and translation by \(v_i\) preserves evaluation at \(s_i\).
The constant digit is always read modulo \(q\).
The linear digit moves by ordinary integers.

Extend the labels arbitrarily to integer linear digits
outside the digit box to define translation operators
\(\mathrm T_u\). Expand each operator by
\(\Delta_{u_i+v_i}=\Delta_{u_i}+\mathrm T_{u_i}\Delta_{v_i}\).
The product of the \(k\) operators is a sum of \(2^k\)
translated mixed differences, each using \(k\) distinct
kernel levels. Use base points with linear digit in
\([k,s_1-1-k]\), so all evaluated vertices are inside the
original box. For a fixed shift tuple, each evaluated vertex
is a translate of the base and the translation is injective
on this restricted box. Thus for a uniformly random base the
probability that any required vertex is outside the pertinent
\(G_t\) is at most \(O_k(\epsilon^{1/4})\), by a union bound
over the fixed number of terms, levels and vertices.

For every term for which all its vertices pass, take one of
its kernel directions at level \(t\). Pairing the vertices
along that direction cancels the term modulo \(M_t\).
The term is therefore divisible by the least common multiple
of its \(k\) moduli. Its absolute value is at most \(2^kX\)
and (EL1) makes it zero over \(\mathbb Z\).
Consequently the \(k\)-fold difference in the pure constant-digit
shifts vanishes with probability at least \(1-O_k(\epsilon^{1/4})\)
over the restricted base and independent uniform shifts.
If some shift is zero the pure difference vanishes directly.

Average over every digit of the base except the constant
digit. One fixed slice gives an injective
\(f:\mathbb F_q\to\mathbb Z\) satisfying (EL2) for \(j=k\)
when \(\epsilon=\epsilon_k\) is small enough, say with
\(\alpha=1/4\). This contradicts injectivity.
The real bound on the labels entered only through the
integer bound \(2^kX\). The proof applies to composite
moduli if they obey (EL1).

This obstruction concerns near-complete retention of a specified
evaluation template. It does not assert that every balanced abc
graph contains such a template, and it gives no power bound
on a subgraph retaining a vanishing fraction of its incidences.
