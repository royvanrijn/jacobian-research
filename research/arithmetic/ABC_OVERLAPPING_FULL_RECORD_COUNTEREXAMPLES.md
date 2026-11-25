# Simultaneous affordable-record obstruction with complete integer labels

This is an arithmetic obstruction to an inference based only on the two
free letters, primitive complete labels, and entropy of a list of locks.
The third letter constructed here violates the small-radical bound. The
construction can use arbitrarily small fixed budgets and a wholly diffuse
profile on each free letter. It attains the residual counting exponent
after *every* simultaneous choice of fixed individual degree coordinates.

## Fixed parameters and complete integers

Take an integer \(m\) large, set \(L=1/m\), \(n=X^L\), and put
\[
d=\tfrac13,\qquad J=2m-1,\qquad
\beta=d/J,\qquad
\mathcal E=\{2m+1,\dots,4m-1\}.
\]
All asymptotics in \(X\) take \(m\) fixed. For each of the \(2J\) distinct
root slots, select a disjoint pool of \(P\asymp_m n^\beta/\log X\) primes
in \([n^\beta/4,n^\beta/2]\), with exactly the same cardinality \(P\) for
each slot. The prime number theorem gives enough primes for these pools.
Independently choose a root in every slot and form
\[
U=\prod_{e\in\mathcal E}p_e^e,\qquad
V=\prod_{e\in\mathcal E}q_e^e,\qquad T=U+V .
\]
Thus \(U,V,T\) are pairwise coprime. The \(U,V\) dictionaries each have
size \(M_0=P^J=n^{d-o_X(1)}\) and their entire Cartesian product is
available. Each of \(U,V,T\) is \(\asymp_m X\), and they are all at
most \(X\) for sufficiently large \(X\): indeed
\(\beta L\sum_e e=1\) while the factor \(2^{-\sum_e e}\) gives
\(U+V\le X\). Also
\(\operatorname{rad}(U),\operatorname{rad}(V)\le n^d\).
The scaled root mass measure of each of \(U,V\) is, in the \(X\)-limit,
\(\beta\sum_{e\in\mathcal E}\delta_{eL}\). As \(m\to\infty\) it tends
to \(d\,\mathrm{Unif}[2,4]\).

## A uniform fixed-target bound

Select all degree slots in both \(U,V\) satisfying \(eL\ge79/25\).
Let \(t_m\) denote their normalized root mass and \(w_m\) their normalized
ordinary logarithmic weight. Then
\[
t_m\longrightarrow 7/25<d,\qquad
w_m\longrightarrow \tfrac13\int_{79/25}^{4}r\,dr
=\tfrac{1253}{1250}>1.
\]
For large fixed \(m\), each assignment of these root slots gives coprime
divisors \(F\mid U\) and \(G\mid V\) with \(FG>X\), uniformly for
sufficiently large \(X\). At a fixed \(T\), the congruences
\(U\equiv0\bmod F\), \(U\equiv T\bmod G\) determine \(U\) uniquely in
\([1,X]\). Consequently a fixed \(T\) has at most
\(P^{t_m/\beta}\le n^{t_m}\) representations in the Cartesian product.

Uniformly for \(1\le R\le X\), the number of integers \(T\le X\) with
\(\operatorname{rad}(T)=R\) is \(X^{o(1)}\). For example, for every
\(\epsilon>0\) it is at most
\[
 X^\epsilon\prod_{p\mid R}(1-p^{-\epsilon})^{-1}
 \le X^\epsilon \exp(o_\epsilon(\log X)).
\]
For the last bound, separate the finitely many primes below a fixed
cutoff \(Y\), and bound the logs of the remaining factors by
\(O_\epsilon(Y^{-\epsilon}\log X)\), then let \(Y\to\infty\). Since
\(m\) is fixed, the number of candidate edges with
\(\operatorname{rad}(T)\le n^d\) is at most
\(n^{d+t_m+o_X(1)}\).

## Simultaneous entropy saturation

Sample each Cartesian edge independently with probability
\(\pi=(M_0\log X)^{-1}\), so the mean sample size is
\(S=M_0/\log X=n^{d-o_X(1)}\). The expected number of colliding
unordered pairs in either \(U\) or \(V\) is \(O(S/\log X)\).
The fixed-target bound makes the expected number of colliding pairs in
\(T\) at most \(\pi^2 M_0^2 n^{t_m}=o(S)\). The expected number
of selected edges with \(\operatorname{rad}(T)\le n^d\) is also \(o(S)\).

For any of the finitely many subsets \(I\) of the \(2J\) slots, write
\(Y_I\) for the vector of the selected prime roots and \(t_I=\beta |I|\).
The candidate fibre of each assignment to \(Y_I\) has size
\(P^{2J-|I|}\). Its sampled size has mean
\[
\mu_I=\pi P^{2J-|I|}
      =n^{d-t_I+o_X(1)} .
\]
For all \(I\) and all assignments simultaneously, all sampled fibres have
size \(O_m(\max(\log X,\mu_I))\) with probability \(1-o(1)\). Apply the
binomial upper-tail bound at threshold
\(\max(2\mu_I,K_m\log X)\), with \(K_m\) large enough for the union bound
over at most \(2^{2J}\max(1,P^{2J})\) assignments. This property remains
valid under arbitrary deletions.

There is a sample with at least \(S/2\) edges, with this fibre property
and with \(o(S)\) colliding pairs and low-radical edges. This follows by
the binomial lower-tail bound and Markov's inequality. Delete all
low-radical edges and at most one edge per collision. The remaining
\(\mathcal H\) has size \(n^{d-o_X(1)}\), is a matching in **all three
complete integer labels**, and has \(\operatorname{rad}(T)>n^d\) on
all edges.

Use natural-log entropy for the uniform distribution on \(\mathcal H\).
As \(T\) determines the edge, the dictionary size of \(Y_I\) and the
simultaneous fibre caps give
\[
\begin{aligned}
H(Y_I)&=\min(d,t_I)\log n+o_X(\log n),\\
H(T\mid Y_I)&=(d-t_I)_+\log n+o_X(\log n).
\end{aligned}
\]
Indeed \(H(Y_I)\le\min(\log|\mathcal H|,\,|I|\log P)\) and
\(H(T\mid Y_I)\le\log\max_y\#\{T:Y_I=y\}\), while their sum is
\(\log|\mathcal H|\). For \(t_I<d\), there is in particular a single
assignment admitting at least \(n^{d-t_I-o_X(1)}\) distinct targets.

The formula applies to any finite collection of records, using the
union \(I\) of the slots they name. A product of their fixed prime-power
coordinates and the vector of those coordinates determine one another
because the pools are disjoint. For two explicit distinct affordable
records, start with the common slots of scaled degree at least \(79/25\)
and add the slots in \([313/100,79/25)\) on the first letter in the
first record, on the second letter in the second record. Each addition
has limiting mass \(1/200\). Each record has limiting cost \(57/200\)
and the union has limiting cost \(29/100<d\). Both records already
have numerical weight greater than one before the additions. The common
fibre count can still attain the remaining exponent \(d-29/100\) in
the two limits.

## A relative-residue statement that does survive

Suppose two records prescribe \(F\mid U_i\) and \(G_i\mid V_i\) in
\(U_i+V_i=T_i\), for \(i=1,2\), and
\(\gcd(F,G_1G_2)=1\). Put \(G_\cap=\gcd(G_1,G_2)\). If
\(T_1\equiv T_2\pmod {G_\cap}\), then
\[
 F G_\cap\mid U_1-U_2 .
\]
In particular, when \(FG_\cap\) exceeds the diameter of the possible
\(U\)-interval, \(U_1=U_2\). For exact equal target values, \(V_1=V_2\)
also, and the intersection of the two target dictionaries equals the
target dictionary of the union lock
\((F,\operatorname{lcm}(G_1,G_2))\), assuming the records restrict
the same ambient family. The same statement holds with a common
\(G\) and varying \(F_i\), with the corresponding roles exchanged.

This relative congruence can be applied to bounded affine relations
as well. For instances with common lock \(F\mid U_i,G\mid V_i\),
if \(\sum_i c_i T_i\equiv0\bmod G\) and
\(|\sum_i c_iU_i|<FG\), then \(\sum_i c_iU_i=0\). When
\(\sum_i c_iT_i=0\) exactly, the relation then holds in the \(V_i\)
as well. These lifts recover relations *inside the prescribed common
overlap*. The construction above shows that their existence, complete
matching, and the nested fixed-slot entropy conditions alone permit
the full residual exponent.


## Simultaneous removal of bounded-order additive relations

The entropy saturation can hold while every affordable record, including
every affordable union of records, has no nontrivial relations among at
most eight distinct targets with any prescribed bounded coefficients.
This gives an obstruction even when the target dictionaries have
minimal additive energy. It still prescribes no small target radical.

Use the same random sample. The combined scaled mass of \(U,V\) is
uniformly approximated, as \(m\to\infty\), by the measure with density
\(1/3\) on \([2,4]\). At mass \(t\), the limiting maximum weight of
a set of slots is the top-mass integral \(4t-\tfrac32t^2\).
This can also be seen directly by taking the largest degrees first,
and the discrete Riemann sums converge at each fixed \(t\).
At \(t_*=139/500\) the maximum is \(498037/500000<1\).
Thus for sufficiently large fixed \(m\), every slot subset of ordinary
weight \(w_I>1\) has \(t_I>139/500\). Because there are only finitely
many \(I\) at fixed \(m\), at sufficiently large \(X\) every actual
assignment to any \(I\) with \(w_I>1\) has \(FG>X\).
Therefore its Cartesian fibre contains at most one edge for each target.

Fix \(2\le k\le8\) and nonzero integers \(c_1,\ldots,c_k\) from a
fixed finite coefficient range containing ({-4,ldots,4}). In each fibre of \(I\) with \(w_I>1\),
there are at most \(P^{(k-1)(2J-|I|)}\) ordered tuples of distinct
target edges satisfying \(\sum_{i=1}^k c_iT_i=0\), since the
first \(k-1\) target values determine the last one. The expected number
of selected tuples over all its \(P^{|I|}\) fibres is at most
\[
 \pi^k P^{|I|+(k-1)(2J-|I|)}
 = n^{(k-2)(d-t_I)+o_X(1)} .
\]
The exponent is at most
\(6(1/3-139/500)=83/250<1/3\).
Summing over all fixed coefficient patterns, \(k\), and \(I\)
there are \(o(S)\) such tuples in expectation. Use Markov's inequality
again when selecting the sample and delete one edge per surviving
relation. This preserves the matching size and all the entropy estimates.
By combining equal entries first, every equality between sums of up
to four targets is now a permutation of the summands, for every
affordable fixed record and every affordable fixed union. In particular
every such target dictionary \(\mathcal D\) has exactly
\(2|\mathcal D|^2-|\mathcal D|\) ordered equal-sum quadruples.


## Prime-power congruences alone need not give sparse unit reductions

Let \(e_1,\ldots,e_k\) be the active degrees on a free letter, let
\(D\ge2\), and temporarily let its roots vary through all units modulo
\(D\). Write \(A=(\mathbf Z/D\mathbf Z)^\times\) and
\(g=\gcd(e_1,\ldots,e_k)\). The image of
\((z_1,\ldots,z_k)\mapsto\prod_j z_j^{e_j}\) is exactly
\(A^g=\{z^g:z\in A\}\). One inclusion is immediate, and the other
follows by applying a fixed integer Bezout relation for the degrees to
an arbitrary \(z\in A\). In particular if two consecutive degrees
are active, the image is all of \(A\) for every \(D\), even a product of
many prime powers from two or more records. For an odd prime power
\(D=p^a\), the same image has index \(\gcd(g,\varphi(D))\), at most
\(g\), by cyclicity. At a fixed budget and for bounded scaled degrees,
this index is independent of the size of \(p\).

Thus reducing at primes assigned to the other letter gives no small
unit image solely from the unrestricted power-map description. For
example, reducing \(U+V=T\) modulo a divisor of \(U\) gives \(T\equiv
V\), and the same is true simultaneously modulo any list of divisors
of \(U\) by passing to their least common multiple. Actual root interval
cuts, root dictionaries, or interaction with the target factorization
can reduce the image. Those features are essential to any sparse-image
inference from this calculation.


## Pairwise overlap forces relative residue spread in a matching

Take any list of edges \(U_i+V_i=T_i\) with assigned records
\((F_i,G_i)\), where \(F_i\mid U_i,G_i\mid V_i\). Assume that the
\(U_i\) are all distinct and range in an interval of diameter \(B\).
For an arbitrary integer modulus \(M\), every pair with
\(T_i\equiv T_j\bmod M\) satisfies
\[
D_{ij}(M)=\gcd(F_i,F_j)\gcd(G_i,G_j,M)\mid U_i-U_j .
\]
Indeed the second factor divides both \(V_i-V_j\) and \(T_i-T_j\).
The factors are coprime by primitivity of either edge. Hence if
\(D_{ij}(M)\ge Q\) for all distinct pairs in a subfamily, each residue
class of the targets modulo \(M\) contains at most
\(1+\lfloor B/Q\rfloor\) members of the subfamily. This follows by
sorting the corresponding \(U_i\), whose adjacent gaps are each at
least \(Q\). In particular if \(Q>B\), target reduction modulo \(M\)
is injective on that subfamily. One may use the same proof with the
opposite letter in the interval or with a signed equation.

This is a necessary relative-residue condition across records with
possibly different locks. It uses neither the factorization of the
targets nor their height. It does require a lower bound for the shared
divisors, and the magnitude of \(M\) can greatly exceed the residual
dictionary size when the condition gives injectivity.

## Failure of bounded affine covers from the two free radicals

The [affine rank test](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#sparsity-for-finitely-many-small-integer-parameters)
needs no target radical bound if the two free spectra are
wholly diffuse. Its minor calculation is unchanged. In the
rank reduction, use witness differences only on the moving
free letters and strip their supported parts. If no moving
free letter has positive residual diffuse mass, its remaining
mass is zero by domination by the original diffuse measure.
Enumerate both residual dictionaries at arbitrarily small
cost. Otherwise a positive band on a moving free letter
forces another rank drop, as before.

Consequently every large fixed-record fibre of the small-budget
construction above requires at least \(n^{d-t-\epsilon-o(1)}\)
bounded-feature systems for any fixed \(\epsilon>0\) with
width chosen by that test, if \(t\) is its prescribed mass.
This obstruction uses actual additive heights and can hold
simultaneously with minimal additive energy in every
affordable target dictionary. It places the missing cover
input in the target arithmetic, beyond the primitive
relation and the two free radical caps.

## Three radical caps and two affordable records at budget \(L=3/2\)

There is a primitive arithmetic model with an exact target radical cap
in which two distinct affordable records and their union attain the
residual exponent. Its fixed budget is \(L=3/2\), outside the small-budget
regime. This supplements the arbitrarily small-budget example with
unrestricted target radical.

For large \(N\), let \(\mathcal P_0,\mathcal Q_0\) be disjoint prime
pools in \( [N^{13/20},2N^{13/20}]\), both of size
\(\asymp N^{13/20}/\log N\). Let \(\mathcal P_1,\mathcal Q_1\)
be disjoint prime pools in \([N^{3/20},2N^{3/20}]\), both of size
\(\asymp N^{3/20}/\log N\). The four pools are mutually disjoint.
For each choice \(p_i\in\mathcal P_i,q_i\in\mathcal Q_i\), set
\(P=p_0p_1\) and \(Q=q_0q_1\). Use parameters \(x,y\) satisfying
\[
3N\le x\le4N,\quad N\le y\le2N,\quad
x\equiv1\pmod {72},\quad y\equiv12\pmod {72},\quad
x=P u,\quad y=Q v,\quad
\gcd(u,v)=1,\quad \gcd(uv,PQ)=1 .
\]
There are \(\gg N^{2/5}\) choices for every prime tuple, uniformly.
In fact \(u,v\) range in specified classes modulo 72 in intervals of
lengths \(H_0=N/P,H_1=N/Q\asymp N^{1/5}\). Neither 2 nor 3 divides
\(u\). For a prime \(\ell\ge5\) the number of pairs with \(\ell\mid
\gcd(u,v)\) is at most
\((H_0/(72\ell)+1)(H_1/(72\ell)+1)\).
Sum through \(O(N^{1/5})\), using
\(\sum_{\ell\ge5}\ell^{-2}\le\sum_{k\ge5}k^{-2}<1/4\);
the cross terms cost \(O(N^{1/5}\log N)\).
The four additional primes of \(PQ\) exclude \(o(H_0H_1)\) pairs
by the same single-coordinate progression count.

The total number of admissible \(x,y\), allowing the primes to vary,
is \(\gg N^2/(\log N)^4\). To see this, multiply the fixed-tuple
count by the four pool cardinalities, and note that each admissible
\(x,y=O(N)\) has only boundedly many divisors in the four prime
ranges. Each such prime has size at least \(N^{3/20}\).
The constraints also imply \(\gcd(x,y)=1\).

Now form
\[
 U=x^3(x+2y),\qquad V=y^3(2x+y),\qquad
 T=(x+y)^3(x-y),
 \qquad U=V+T .
\]
All complete values are \(\asymp N^4\). The six linear forms
\(x,y,x+2y,2x+y,x+y,x-y\) have two-by-two coefficient determinants
nonzero and of absolute value at most 3. Thus any prime at least 5
divides at most one of them. The imposed residue classes give that
2 and 3 divide neither \(U\) nor \(T\). Therefore the three labels
are pairwise coprime. This argument also shows that each selected
\(p_i,q_i\) has valuation exactly 3 in its complete label, so the
locks below are complete prime-power divisors.

All six linear forms have radical \(N^{1-o(1)}\) after discarding
\(o(N^2/(\log N)^4)\) parameter pairs. Indeed, for each fixed
\(\epsilon>0\) there are at most \(N^{1-\epsilon+o(1)}\) values
\(z\ll N\) with radical at most \(N^{1-\epsilon}\) by the
fixed-radical support estimate, and each value of any one form has
\(O(N)\) parameter pairs. One can take \(\epsilon=\epsilon_N\to0\)
sufficiently slowly after this estimate. The two forms in each complete
letter have greatest common divisor dividing 6, so all three complete
radicals are \(N^{2-o(1)}\). Each also is at most a fixed constant
times \(N^2\).

The surviving parameter pairs number \(N^{2-o(1)}\).
For each such pair choose one admissible prime tuple. Greedily retain
a matching in the three complete labels, of size \(N^{2-o(1)}\).
This costs only a factor \(N^{o(1)}\) since a fixed \(U\), \(V\) or
\(T\) has at most \(N^{o(1)}\) parameter pairs. In fact \(x^3\mid U\),
\(y^3\mid V\) and \((x+y)^3\mid T\), respectively, give divisor-many
possibilities for \(x\), \(y\) or \(x+y\). Each determines the other
parameter via the respective complete value and defining formula. So
each selected edge removes at most \(N^{o(1)}\) collision neighbours.

Take \(X=(100N)^4,\ n=(100N)^6\) and \(d=1/3\).
The complete radical bounds are all at most \(n^d\), and the complete
values are at most \(X\). On this matching consider the two records
on \(U,V\) and their union
\[
\begin{aligned}
Y_1&=((p_0p_1)^3,q_0^3),\\
Y_2&=(p_0^3,(q_0q_1)^3),\\
Y_1\vee Y_2&=((p_0p_1)^3,(q_0q_1)^3).
\end{aligned}
\]
Each record has normalized cost \(29/120+o(1)<d\) and ordinary
weight \(87/80+o(1)>1\), so each is affordable. The union has cost
\(4/15+o(1)\). The named prime coordinates can be recovered from each
record value by the disjoint prime ranges.

More generally, let \(Y_I\) name any subset \(I\) of the four prime
coordinates, and \(a_I\) be the sum of their log-\(N\) weights, taking
weights \(13/20\) for \(p_0,q_0\) and \(3/20\) for \(p_1,q_1\).
There are at most \(O(N^{a_I})\) possible assignments to \(Y_I\).
For a fixed assignment, \(x\) must be a multiple of the product of
the selected \(p\)-primes, and \(y\) a multiple of the product of
the selected \(q\)-primes. The sum of the selected weights on either
letter is at most \(4/5<1\). Thus the two box ranges give at most
\(O(N^{2-a_I})\) parameter pairs in each fibre, uniformly. In the
uniform distribution on the complete matching, the same entropy
argument as in the small-budget construction gives
\[
 H(Y_I)=a_I\log N+o(\log N),\qquad
 H(T\mid Y_I)=(2-a_I)\log N+o(\log N).
\]
There is one matching edge at which **every** associated fibre also
attains its residual exponent. Indeed, for any fixed
\(\epsilon>0\), the \(Y_I\) assignments with fibre cardinality smaller
than \(N^{2-a_I-\epsilon}\) account for only \(O(N^{2-\epsilon})\)
edges. Sum this estimate over the finitely many \(I\), then take
\(\epsilon\to0\) sufficiently slowly.
In particular this edge has \(N^{11/20-o(1)}\) distinct targets in
each affordable parent fibre and \(N^{2/5-o(1)}\) in their common
union fibre, exactly the respective normalized residual exponents.
