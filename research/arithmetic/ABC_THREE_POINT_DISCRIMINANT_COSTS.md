# Discriminant costs of three-point covers on diffuse prime degrees

A tame specialization of a three-point cover sees the original
valuation through powers of a branch permutation. The gain in
its normalized discriminant exponent is bounded by the largest
single-prime divisibility mass of the degree distribution. This
bound is uniform in the degree and the monodromy group of the
cover. It applies even when the cover is chosen after the
primitive valuations are known.

## An inequality for every branch permutation

Let \(\sigma\) be a permutation of \(d\) letters. Write
\(\operatorname{ind}(\sigma)=d-\#\operatorname{cycles}(\sigma)\).
For nonnegative, finitely supported weights \(\beta_e\) on positive
integers put
\[
   b=\sum_e\beta_e,\qquad
   v_*=\sup_{\ell\ {\rm prime}}\sum_{\ell\mid e}\beta_e.
\]
Then
\[
        \sum_e\beta_e\operatorname{ind}(\sigma^e)
                  \ge (b-v_*)\operatorname{ind}(\sigma).
                                                            \tag{DI1}
\]
Indeed, a cycle of length \(m\) splits into \(\gcd(m,e)\)
cycles in \(\sigma^e\). Its decrease in index is
\[
            \gcd(m,e)-1
               =\sum_{\substack{u\mid m,\ u\mid e\\u>1}}\varphi(u).
\]
For each \(u>1\), its valuation weight \(\sum_{u\mid e}\beta_e\)
is at most \(v_*\) by using any prime factor of \(u\).
Summing the decreases gives at most \(v_*(m-1)\) for this
cycle, and then (DI1) for the permutation. A permutation
with one nontrivial cycle of a prime length \(\ell\) attains
the corresponding bound with \(v_*\) replaced by the mass
of degrees divisible by \(\ell\). Thus the best bound in
terms of the single scalar \(v_*\) is attained in a single
branch calculation.

For pairwise coprime degree labels \(e\), the parameter \(v_*\)
is at most the largest individual weight. In particular if
the weights are uniform on \(s\) distinct prime degrees then
\[
         \frac1s\sum_e\operatorname{ind}(\sigma^e)
              \ge (1-1/s)\operatorname{ind}(\sigma).
                                                            \tag{DI2}
\]
No bound on the cycle lengths is needed.

There is also an Artin conductor version at a tame place.
For a finite-order linear transformation \(\sigma\) in
characteristic zero, put \(a(\sigma)=\dim V-\dim V^\sigma\).
An eigenvalue other than one becomes one in \(\sigma^e\)
only when its order, which is greater than one, divides \(e\).
Applying the same prime-divisibility bound separately to each
eigenvalue gives
\[
                 \sum_e\beta_e a(\sigma^e)
                          \ge (b-v_*)a(\sigma).
\]
For the permutation representation this is (DI1).
[Serre, Local Fields, Chapter VI, §2](https://link.springer.com/book/10.1007/978-1-4612-2241-9)

## Normalized discriminant of a three-point specialization

Let \(t=a/c\) for an actual positive primitive triple \(a+b=c\),
and let a finite cover \(\pi:Y\to\mathbf P^1\) of degree \(d\)
over a number field \(K\) be geometrically connected and branched
only above \(0,1,\infty\). Assume \(g(Y)\ge1\). Let \(E_t\)
be the degree-\(d\) finite étale algebra of its fibre at \(t\).
For each role \(i=a,b,c\), choose a set \(P_i\) of primes dividing
that role. Let \(R>1\) be the radical normalization, and let
\(\mathcal B\) contain every selected prime at which the following
specialization hypotheses are unavailable: \(p\) is unramified
in \(K\), and at every prime of \(K\) over \(p\) the tame
specialization inertia is the cycle action of the relevant branch
generator raised to the contact order. For any fixed cover and base field the excluded primes belong
to a fixed finite set. For a non-Galois cover one can apply the
specialization theorem to its Galois closure and the action on
the sheets. A fixed constant extension of the base is unramified
away from a fixed set of primes.
For a moving cover or base field the following estimate charges
the excluded primes on the selected support.

Define the original and exceptional masses, and the maximal
divisibility mass at good primes, by
\[
 \begin{split}
    b_i&=\sum_{p\in P_i}\log_Rp,\qquad
    u_i=\sum_{p\in P_i\cap\mathcal B}\log_Rp,\\
    v_i&=\sup_{\ell\ {\rm prime}}
            \sum_{\substack{p\in P_i\setminus\mathcal B\\
                            \ell\mid v_p(i)}}\log_Rp .
 \end{split}
\]
Let \(\sigma_0,\sigma_1,\sigma_\infty\) be the geometric
branch permutations on the \(d\) sheets, assigned in order to
\(a,b,c\). Set \(I_i=\operatorname{ind}(\sigma_i)/d\).
Then the relative discriminant ideal of the finite étale algebra satisfies
\[
 \frac{\log_R N_{K/\mathbf Q}\mathfrak d(E_t/K)}{d[K:\mathbf Q]}
                       \ \ge\ \sum_i (b_i-u_i-v_i)I_i .
                                                            \tag{DI3}
\]
In particular, if \(b_i\ge b_0\) and \(v_i\le v_0\) for
every \(i\), with \(b_0\ge v_0\), then
\[
 \frac{\log_R N_{K/\mathbf Q}\mathfrak d(E_t/K)}{d[K:\mathbf Q]}
                     +\sum_i u_i
                              \ \ge\ 2(b_0-v_0).
                                                            \tag{DI4}
\]
This estimate is uniform in \(d\), in the geometric monodromy
group and in the inertia orders.

For the proof, primitivity identifies the contact order at
\(p\in P_i\) with the original valuation \(v_p(i)\).
For example at \(p\mid b\) use \(1-t=b/c\). At each good place
over \(p\), the inertia orbits on the sheets give the tame
discriminant exponent
\(\operatorname{ind}(\sigma_i^{v_p(i)})\). This follows from
the tame local formula \(\sum_j f_j(e_j-1)\) for the
discriminant exponent of an étale algebra. An inertia orbit has
length \(e_j\), and each local factor contributes \(f_j\) such
orbits. Sum the residue degrees of the places above the
unramified rational prime to get \([K:\mathbf Q]\).
Now sum (DI1) on the good primes of each role, dropping the
nonnegative discriminant contributions at other primes.
This proves (DI3). By Riemann–Hurwitz
\(\sum_i I_i=2+(2g(Y)-2)/d\ge2\), and each \(I_i\le1\).
Adding \(\sum_i u_i\) proves (DI4).
[Tame specialization on fibre factors, Conrad, Theorem 4.1 and Corollary 4.3](https://math.stanford.edu/~conrad/papers/inertiapaper.pdf)
· [Darmon–Granville, Propositions 3.1–3.2](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)

For a disconnected non-Galois fibre the bound is for the full
étale algebra. Individual factors may have smaller normalized
discriminants. If the cover is regular Galois with constant
deck group over \(K\), each fibre factor is a copy of the same
Galois specialization field as a \(K\)-algebra. Indeed the
finite fibre is a torsor for the constant deck group, and its
Galois orbits are translates of a subgroup orbit, all with the
same point field. In this regular situation the same normalized
lower bound applies to each factor individually.

## Prime valuation profiles with actual quadratic realizations

The construction of (CV1–CV3) works with \(I_J\) replaced by
the set \(E_J\) of prime integers in \([J,2J]\). Fix \(J\)
large, let \(s=|E_J|\), \(S=\sum_{e\in E_J}e\), and set
\[
                    X_0=(2Y)^{WS},\qquad L=3s/S .
\]
Here \(W\) is the least common multiple of the orders in
any prescribed finite list of power signatures, including
\((2,2,2)\), as in (CV1). For fixed \(J,W\) take the height
limit \(Y\to\infty\) first. The Ramsey construction assigns
blocks of \(W\) support primes in \([Y,2Y]\) to the \(s\) degrees
on each of the three letters. It gives \(X_0^{L+o_Y(1)}\)
possible support sets. Its local-solvability proof uses only
the length \(W\) of each block and the exact degree assigned
to it. The ternary quadratic reflection and CRT realization
also use only the prescribed valuations and the caps
\(Y^{WS}\le P_i\le X_0\). Thus the same primitive sums exist
with \(X_0^{1-o_Y(1)}\le c\ll X_0^{6+17L/3}\). Fresh primes
in those realizations can increase the radical.
[Simultaneous local power tests and quadratic realization (CV1–CV3)](ABC_SIMULTANEOUS_LOCAL_POWER_PROFILE_OBSTRUCTION.md)

Use \(R=X_0^L\) in (DI3)–(DI4). On the planted supports,
each role has mass at least
\(\log Y/(3\log(2Y))=1/3-o_Y(1)\). Since all planted degrees are
distinct primes and each degree block has mass at most \(1/(3s)\),
one has \(v_i\le1/(3s)\) even before removing bad primes.
Consequently every primitive realization of the planted profile
satisfies, for every cover in (DI4),
\[
 \frac{\log_{X_0^L} N_{K/\mathbf Q}\mathfrak d(E_t/K)}
                         {d[K:\mathbf Q]}
     +\sum_i u_i
                      \ \ge\ \frac23(1-1/s)-o_Y(1).
                                                            \tag{DI5}
\]
The error comes just from the prime interval and is independent
of the cover. For \(J\to\infty\), the prime number theorem
gives \(s\to\infty\), \(L\sim2/J\), and uniform distribution of
the scaled degrees \(Le\) on \([2,4]\) after normalization by
the planted radical on each letter. This is a diffuse
balanced profile. Any integer divisibility order greater than
one captures at most one planted degree on a letter.

In particular choosing a new three-point cover or large inertia
orders for a particular prime-degree profile cannot substantially
reduce the generic ramification cost on that profile. Formula
(DI5) charges every excluded planted prime separately, even if
the bad-prime set changes with the cover. The result is a lower
bound on this ramification cost, not a bound on the number of
specialization fields of that discriminant. It gives no
primitive realization with the unrestricted exceptional-set
radical cap, and it leaves open arguments using the actual
height and simultaneous original radical constraints.
