# Exact simultaneous overlap bounds for two shifts

The two-difference Motzkin proof constructs a character from a signed dual-grid orientation. In that proof the occupied set has no forbidden pairs. Counting the marked edges of an arbitrary set gives separate budgets for the two directional overlaps. Filling the fibres of the resulting character proves the converse. [Pinned primary proof, The upper bound for finite groups](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/two-difference-motzkin/PROOF.md)

## The exact finite feasible region

Let \(H=\langle a,b\rangle\) be a finite abelian group of order \(N\), with \(a,b\ne0\). For \(A\subseteq H\) set
\[
 e_d(A)=|A\cap(A+d)|,\qquad d\in\{a,b\}.
\]
For an integer \(0\le s\le N/2\) and real budgets \(B_a,B_b\ge0\), there is a set \(A\) of size \(s\) satisfying \(e_d(A)\le B_d\) in both directions if and only if there is a character \(\psi:H\to\mathbb R/\mathbb Z\) satisfying
\[
                       N\|\psi(d)\|\ \ge\ s-B_d
                       \qquad(d=a,b).                                  \tag{1}
\]
The notation \(\|\cdot\|\) means distance to the nearest integer. For \(s>N/2\), pass to the complement, which has size \(N-s\) and budgets \(B_d-(2s-N)\). A negative resulting budget precludes feasibility. For a larger ambient finite group, apply the criterion separately in each coset of \(H\) and add sizes and budgets.

For necessity, write \(x\) for the indicator of \(A\). On the positive edge of type \(d\) from \(z\), mark it if its endpoints have the same occupation. There are precisely
\[
                   N-2s+2e_d(A)
\]
marked edges in direction \(d\). This follows by summing \(|1-x(z)-x(z+d)|\), which is \(1-x(z)-x(z+d)+2x(z)x(z+d)\). The parity of a mark is \(1+x(z)+x(z+d)\bmod2\). The sum around any unit square is even, including squares in the quotient with identified corners or parallel labeled sides. Orient the marked dual edges by Euler circuits, with each labeled edge used once. As in the primary proof, crossing a positive grid edge defines a signed increment \(w_d\in\{-1,0,1\}\), supported on marked edges. The sum of the increments around each face is zero because every dual vertex has equal in-degree and out-degree.

Lift to the square grid and sum along paths. The resulting integer potential has period increment \(L(u,v)\) for each relation \(ua+vb=0\). The increment is independent of the starting point because the lifted increments are periodic, and averaging starting points gives \(L(u,v)=u m_a+v m_b\), where \(m_d=N^{-1}\sum_z w_d(z)\). The mark parity along a path gives \(L(u,v)\equiv u+v\bmod2\), since only its endpoint occupation terms survive, and those endpoints project to the same vertex. Thus the rule
\[
                  (u,v)\longmapsto
                  \tfrac12 u(1-m_a)+\tfrac12 v(1-m_b)\pmod 1
\]
descends to a character \(\psi\) on \(H\). The marked-edge count gives
\[
 |m_d|\le 1-2s/N+2e_d(A)/N,\qquad
 N\|\psi(d)\|\ge s-e_d(A).
\]
This proves the two necessary inequalities in (1) with a single character.

For sufficiency, let the character image have order \(m\), with fibre size \(g=N/m\). Fill the first \(\lfloor s/g\rfloor\) cyclic fibres and an arbitrary part of the next fibre to get exactly \(s\) points. Write \(h\) for the distance of \(\psi(d)\) from zero in image steps, so \(N\|\psi(d)\|=gh\). If \(h>0\), the two filled intervals of fibres meet in exactly the amount \((s-gh)_+=\max(0,s-gh)\). To see this it suffices by symmetry to translate by \(+h\), \(h\le m/2\). The overlap starts \(h\) fibres into the interval. If it reaches the final partial fibre, that fibre is intersected with an entirely full fibre of the translated interval. Wraparound can occur only at the halfway endpoint when \(m\) is even and both lengths reach \(m/2\), and the partial fibre is then empty. For \(h=0\), the overlap is at most \(s=(s-gh)_+\). Thus this single set satisfies \(e_d(A)\le (s-N\|\psi(d)\|)_+\le B_d\) for both shifts.

Writing \(\alpha(H;a,b)\) for the unrestricted two-difference independence number, (1) gives the exact symmetric budget
\[
               \min_{|A|=s}\max(e_a(A),e_b(A))
                         =\max(0,s-\alpha(H;a,b))
               \qquad(0\le s\le N/2).                                  \tag{2}
\]
This follows by taking the same budget in both directions and maximizing the lesser character distance. For \(H=\mathbb Z/24\mathbb Z\), \(a=3,b=4,s=11\), the distinct undominated overlap-budget pairs are \((0,3),(2,0)\). They are obtained by running through the characters \(z\mapsto tz/24\) in (1). In particular the two shifts separately permit zero overlap at size \(11\), while a common bound of one is infeasible. The counts \(e_d\) use labeled positive-direction pairs, so a jump of order two counts both orientations through its two possible start sites.

## Consequences for shifted arithmetic gap statistics

Consider a subset of the integers in growing intervals. For fixed positive shifts \(a,b\), let \(\rho\) and \(c_a,c_b\) be simultaneous subsequential limits of the normalized occupied count and the two normalized overlap counts. Periodize each interval, keeping its subset as the cyclic fibre. Choosing the interval length within \(O_{a,b}(1)\) of an integer coprime to \(\gcd(a,b)\) ensures the two jumps generate the cyclic group. The periodization changes the counts at the boundary by only \(O_{a,b}(1)\). Any limiting pair of characters from (1) obeys
\[
  B\,\psi(a)-A\,\psi(b)=0\pmod 1,\qquad A=a/\gcd(a,b),\quad B=b/\gcd(a,b).
\]
If \(A+B\) is odd, choose representatives
\(\psi(a)=(1+\eta_a)/2\), \(\psi(b)=(1+\eta_b)/2\) in \([0,1]\). The relation implies \(B\eta_a-A\eta_b\) is an odd integer. With
\(|\eta_d|\le1-2\rho+2c_d\), this gives, for \(0\le\rho\le1/2\),
\[
                  Bc_a+Ac_b\ \ge\
                  \max\{0,(A+B)\rho-(A+B-1)/2\}.                       \tag{3}
\]
Equality is realizable for every such \(\rho\). Use a rational character with evaluations \((1+1/(A+B))/2\) and \((1-1/(A+B))/2\). These satisfy the relation, so they are obtained from one character of the integer line on rescaling by \(\gcd(a,b)\). On each residue modulo that gcd, choose periods that are increasing multiples of the character's order, and fill its fibres as above with density approaching \(\rho\). Both shift distances are positive unless \(A+B=1\), excluded here. Thus each overlap tends to \((\rho-(A+B-1)/(2(A+B)))_+\), giving equality. If \(A+B\) is even, both \(A,B\) are odd and the character with both evaluations \(1/2\) realizes zero overlap for every \(0\le\rho\le1/2\).

For a concrete application take the subset to be the gaps of the actual critical random-prime \(k\)-fold sumset inside a finite window. Equations (1)–(3) are unconditional consistency tests for any proposed gap density and two proposed shifted gap-overlap limits. Complementation also permits density above \(1/2\), replacing \((\rho,c_d)\) by \((1-\rho,1-2\rho+c_d)\). The tests concern overlaps of the gap indicators at two target shifts. The [critical intensity route](RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md#a-sufficient-condition-for-the-critical-density) still needs estimates for representation events and shared unexposed primes before it identifies those limits.

Exact subset enumeration matched the character frontier in 1,913 cases, over generating shift pairs and sizes on cyclic groups of orders 3 through 13 and on the products of orders \((2,2),(2,4),(3,3),(2,6)\).
