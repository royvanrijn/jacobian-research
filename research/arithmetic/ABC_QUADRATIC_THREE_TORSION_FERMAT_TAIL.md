# Prime paths and a sharper quadratic 3-torsion tail

For a positive squarefree integer \(d\), put \(h_3(\pm d)=|\mathrm{Cl}(\mathbb Q(\sqrt{\pm d}))[3]|\). For either sign and \(X\ge2,H\ge1\),
\[
 \mathcal N_\pm(H;X):=\#\{d\in[X,2X):d\ {\rm squarefree},\ h_3(\pm d)>H\}
 \ll_\varepsilon X^\varepsilon\left(XH^{-2}+\min\{X^{11/6}H^{-5}+X^{5/2}H^{-8},\,X^3H^{-9}\}\right).                     \tag{1}
\]
The witnesses are paths of three split primes in a single ideal-class coset. A Fermat cubic obstruction guarantees different cube-free parts on many adjacent edges. The coherent principal generators put a full prime cube into the difference-of-squares argument. This strengthens the upper tail in [Heath-Brown–Pierce, Proposition 6.1](https://arxiv.org/pdf/1409.3177) on its useful upper-tail range.

## Witness selection

It suffices to treat \(h_3(-d)\). Scholz reflection gives \(h_3(d)\le h_3(-3d)\) and \(d\mapsto\operatorname{sf}(3d)\) is an involution on positive squarefree integers, mapping one dyadic box into a bounded union of boxes of comparable scale. For \(H\le X^{3/14}\) both decaying bounds in (1) dominate the classical first-moment tail \(X/H\). For \(H>C X^\kappa\), with \(\kappa=1597/5000\), the count vanishes by [Chan–Koymans, Theorem 1.1](https://arxiv.org/pdf/2505.00611). The older \(X^{1/3+o(1)}\) pointwise bound also suffices for the prime-range step below. Choose \(0<\eta\ll1\) fixed and \(Z=X^{1/2+4\eta}/H\). In the remaining range \(X^{1/6}\le Z\le X^{1/3}\) for large \(X\). Equations (2.5)–(2.7) in [Heath-Brown–Pierce](https://arxiv.org/pdf/1409.3177) give \(\gg Z/\log Z\) unramified split primes in \([Z,2Z)\) outside \(O_\eta(X^\eta)\) exceptional \(d\). Select one prime ideal over each split prime and partition them by the quotient \(\mathrm{Cl}(K)/\mathrm{Cl}(K)[3]\) with \(K=\mathbb Q(\sqrt{-d})\). The number of parts is \(O_\eta(X^{1/2+\eta}/H)\) by \(h(K)\ll_\eta X^{1/2+\eta}\). Cosets with fewer than ten selected primes account for \(o(Z/\log Z)\) of them.

In a remaining coset, fix an anchor \(\mathfrak p_0\). Choose a generator \(\beta_i=(a_i+b_i\sqrt{-d})/2\) of the principal ideal \((\overline{\mathfrak p_0}\mathfrak p_i)^3\). Define
\[
 \alpha_{ij}=\overline{\beta_i}\beta_j/p_0^3
              =(u_{ij}+v_{ij}\sqrt{-d})/2.
\]
Its principal ideal is \((\overline{\mathfrak p_i}\mathfrak p_j)^3\), so \(\alpha_{ij}\) is integral, and for \(i\ne j\),
\[
                  4(p_ip_j)^3=u_{ij}^2+d v_{ij}^2,\qquad (p_ip_j,v_{ij})=1.
\]
If an odd prime \(p\mid p_ip_j\) also divided \(v_{ij}\), it would divide \(u_{ij}\), and \(\alpha_{ij}\) would belong to \(p\mathcal O_K\) despite its ideal being supported at only one prime above \(p\). Also \(u_{ij}\) and \(v_{ij}\) are nonzero: a rational or purely imaginary generator gives an ideal invariant under conjugation, whereas this ideal is not.

The signed numbers \(v_{ij}=(a_i b_j-b_i a_j)/(2p_0^3)\) obey the rank-two determinant identity
\[
                     v_{12}v_{34}+v_{14}v_{23}=v_{13}v_{24}.
\]
Color each unordered pair \(ij\) by the positive cube-free part of \(|v_{ij}|\). A monochromatic \(K_4\) would give a nonzero solution of \(s^3+t^3=u^3\) in signed integers, contrary to Fermat's theorem for cubes. Every nonmonochromatic \(K_4\) contains two adjacent edges with different colors, since its line graph is connected. Each adjacent pair is contained in exactly \(s-3\) complete graphs on four vertices in a coset of size \(s\). Therefore the number of different-color adjacent pairs in a coset of size \(s\ge10\) is \(\gg s^3\). By Hölder over the large cosets, each nonexceptional \(d\) in \(\mathcal N_-(H;X)\) supplies
\[
                    \gg_\eta Z^3 H^2 X^{-1-4\eta}                    \tag{2}
\]
ordered paths \(p_1,p_2,p_3\) of distinct primes with different colors on the edges \(12,23\). Indeed \(\sum s^3\ge(\sum s)^3/(\#\text{cosets})^2\); the three factors of \(\log Z\) are absorbed by \(X^{2\eta}\). For each unordered adjacent edge pair one can take a deterministic ordering of its two outer primes. All counted paths are distinct as prime triples in a fixed \(d\).

## The common prime cube and the three small congruence parameters

For a counted path abbreviate \((u_1,v_1)=(u_{12},v_{12})\) and \((u_2,v_2)=(u_{23},v_{23})\). Let \(\delta=\gcd(|v_1|,|v_2|)>0\), \(y_i=v_i/\delta\), and \(k=y_2u_1+y_1u_2\), retaining signs. The two colors differ, so \(y_1/y_2\) is not a rational cube. Normalization of the generators gives \(\alpha_{12}\alpha_{23}=p_2^3\alpha_{13}\). Comparing coefficients of \(\sqrt{-d}\) shows \(k=2 p_2^3 v_{13}/\delta\). Since \(p_2\nmid\delta\), we have
\[
     j:=k/p_2^3=2v_{13}/\delta\in\mathbb Z\setminus\{0\},\qquad
     |y_1|,\ |y_2|,\ |j|\ll V,\qquad V=Z^3 X^{-1/2}.
\]
A fixed constant in each bound is harmless. This uses the actual third edge and bounds all three parameters at the same scale, instead of bounding \(k\) at scale \(Z^3V\). Subtracting the two norm equations after multiplication by the opposite \(y_i^2\) gives
\[
                  j(y_2u_1-y_1u_2)=4(y_2^2p_1^3-y_1^2p_3^3).         \tag{3}
\]
Thus \(p_1,p_2,p_3,y_1,y_2,j\) determine the sum and difference \(y_2u_1\pm y_1u_2\), hence the two signed \(u_i\). Each positive integer \(4(p_i p_{i+1})^3-u_i^2\) determines \(d\) and \(|v_i|\) uniquely by squarefree decomposition. Signed \(y_i\) determine the signs of \(v_i\). In particular the greatest common divisor \(\delta\) is never charged.

Fix \(p_2\). Put \(q_i=\operatorname{rad}_{\rm odd}(|y_i|)\) and \(r=\operatorname{rad}_{\rm odd}(|j|)\). The norm equations and (3) imply
\[
 \begin{split}
 4y_2^2p_1^3&\equiv p_2^3j^2\pmod {q_1},\\
 4y_1^2p_3^3&\equiv p_2^3j^2\pmod {q_2},\\
 y_2^2p_1^3&\equiv y_1^2p_3^3\pmod r.
 \end{split}                                                       \tag{4}
\]
For the first two congruences use \(y_2u_1\equiv k\bmod q_1\), \(y_1u_2\equiv k\bmod q_2\), and cancel the unit \(p_2^3\). The moduli are pairwise coprime. In fact \((y_1,y_2)=1\); an odd prime dividing both \(y_i\) and \(j\) would divide \(k\) and then \(u_i\), forcing it to divide \(p_i p_{i+1}\) despite norm primitivity. The third congruence and \(p_1\ne p_3\) show that both \(p_1,p_3\) are units modulo \(r\).

CRT turns (4) into \(O_\varepsilon(X^\varepsilon)\) lattice cosets in \((p_1,p_3)\) of determinant \(q_1 q_2 r\). At each prime in \(q_i\) one coordinate belongs to one of at most three residue classes; at each prime of \(r\), the ratio \(p_1/p_3\) belongs to one of at most three unit classes. In a homogeneous lattice the conditions are \(q_1\mid\mu_1,q_2\mid\mu_3\) and the corresponding fixed ratio congruence modulo \(r\). Let \(\lambda_1\) denote the length of a shortest nonzero homogeneous vector. As in [Heath-Brown–Pierce, Sections 4.2–4.4](https://arxiv.org/pdf/1409.3177), each coset has at most \(O(1+Z^2/(q_1 q_2 r)+Z/\lambda_1)\) points in the size-\(O(Z)\) box, by successive minima after translating by one point if occupied.

The first term sums to \(O_\varepsilon(X^\varepsilon V^3)\). The second sums to \(O_\varepsilon(X^\varepsilon Z^2)\), by the prescribed odd radical count in Heath-Brown–Pierce Lemma 4.4. For the shortest-vector term, take \(\mu=(\mu_1,\mu_3)\ne0\) of length \(\lambda_1\ll\sqrt{q_1q_2r}\ll V^{3/2}\). If both components are nonzero, divisor choices and the specified-radical count recover \(q_1,q_2,y_1,y_2\) with subpower multiplicity. The odd radical \(r\) divides \(y_2^2\mu_1^3-y_1^2\mu_3^3\). This integer cannot vanish, because \((y_1,y_2)=1\) and their ratio is not a cube. Thus \(r\), then \(j\), have divisor/radical-many choices as well. Summing \(Z/\|\mu\|\) contributes \(O_\varepsilon(X^\varepsilon ZV^{3/2})\). If \(\mu_1=0\), then \(q_2,r\mid\operatorname{rad}_{\rm odd}(\mu_3)\), so only \(y_1\) has \(O(V)\) choices and the contribution is \(O_\varepsilon(X^\varepsilon ZV)\) by the harmonic sum. The other case is symmetric. Both shortest-vector contributions are bounded by \(O_\varepsilon(X^\varepsilon(V^3+Z^2))\), since \(ZV\le Z^2+V^3\) for \(V\ge1\) and \(ZV^{3/2}\) is the geometric mean of \(Z^2\) and \(V^3\).

Summing over \(p_2\) gives \(O_\varepsilon(X^\varepsilon(ZV^3+Z^3))\) witness paths. Dividing by (2) yields
\[
 \mathcal N_-(H;X)\ll_{\eta,\varepsilon}X^{O(\eta)+\varepsilon}
    \left(\frac{X}{H^2}+\frac{Z^7 X^{-1/2}}{H^2}\right)+O_\eta(X^\eta).
\]
Substitute \(Z=X^{1/2+4\eta}/H\) and choose the auxiliary allowances in terms of the final \(\varepsilon\). The exceptional set is absorbed by \(X/H^2\) in the pointwise range. This proves the bound with the second decaying term of (1).

## Full moduli, including the factor 2

As before, \(p_1,p_2,p_3\) are distinct odd primes in \([Z,2Z)\), \(v_i=\delta y_i\), \((y_1,y_2)=1\), and \(p_2^3j=y_2u_1+y_1u_2\), where the signed nonzero \(y_i,j\) are \(O(V)\). Also
\[
          j(y_2u_1-y_1u_2)=4(y_2^2p_1^3-y_1^2p_3^3).
\]
Define
\[
       q_i=\frac{|y_i|}{\gcd(|y_i|,4)},\qquad
       r=\frac{|j|}{\gcd(|j|,4)},\qquad Q=q_1q_2r .
\]
The exact norm and path equations give the congruences
\[
 \begin{split}
 4y_2^2p_1^3&\equiv p_2^3j^2\pmod {|y_1|},\\
 4y_1^2p_3^3&\equiv p_2^3j^2\pmod {|y_2|},\\
 y_2^2p_1^3&\equiv y_1^2p_3^3\pmod r .
 \end{split}                                                      \tag{5}
\]
For each of the first two relations, square the congruence \(y_2u_1\equiv p_2^3 j\bmod |y_1|\) or its analogue, use the norm equation, and cancel the unit \(p_2^3\) by norm primitivity. The last relation follows from \(j\mid 4(y_2^2p_1^3-y_1^2p_3^3)\). No odd prime can divide both \(y_i\) and \(j\): it would divide \(u_i\) and force a prime factor of \(p_i p_{i+1}\) to divide \(v_i\), a contradiction. If \(q_i\) is even, \(v_2(y_i)\ge3\), so the corresponding first or second relation in (5) and oddness of \(y_{3-i}\) imply \(v_2(j)=1\). In this case divide that relation by 4 at 2 to obtain a cube congruence between units modulo \(2^{v_2(y_i)-2}\). If \(r\) is even, the last relation makes \(y_1,y_2\) both odd. Therefore \(q_1,q_2,r\) are pairwise coprime, including their 2-parts.

At odd prime powers, the ratios are units by \((p_i p_{i+1},v_i)=1\) in the first two relations and by \(p_1\ne p_3\) and coprimality of \(r\) to \(y_1,y_2\) in the last. The same unit assertions hold at 2 by the preceding paragraph. A unit has at most three cube roots modulo any odd prime power (the unit group is cyclic) and at most one modulo a power of 2 (cubing is an automorphism of its unit group). Thus CRT covers prime triples by \(O_\varepsilon(X^\varepsilon)\) homogeneous lattices of index \(Q\) in \(\mathbb Z^3\). Each lattice imposes \(p_1=c p_2\bmod q_1\), \(p_3=c'p_2\bmod q_2\), \(p_1=c''p_3\bmod r\), with unit coefficients. They satisfy \(4y_2^2c^3\equiv j^2\bmod q_1\), \(4y_1^2(c')^3\equiv j^2\bmod q_2\), \(y_2^2(c'')^3\equiv y_1^2\bmod r\). At 2 the first two written congruences follow from the stronger divided-unit condition. Each choice of the three moduli has only a bounded number of lifts to \((y_1,y_2,j)\), since each ratio of an original signed variable to its modulus lies in \(\{\pm1,\pm2,\pm4\}\).

## Rank two and the primitive coefficient plane

Each lattice is counted in a ball of radius \(O(Z)\) in \(\mathbb R^3\). Rank one permits at most one positive primitive prime triple per lattice, so contributes \(O_\varepsilon(X^\varepsilon V^3)\). In rank three the bound per lattice is \(O(Z^3/Q)\) by successive minima and the index; summing over the full moduli up to \(O(V)\) costs only the product of three harmonic sums. Its total is \(O_\varepsilon(X^\varepsilon Z^3)\).

In rank two let \(P=\lambda_1\lambda_2\). Use a basis of the entire rank-two intersection lattice in the span of small vectors, and its nonzero integer cofactor \(n=(n_1,n_2,n_3)\). Its norm is comparable to \(P\) by Minkowski's second theorem. The basis itself need not lie in the box. The cofactor is proportional to each local congruence normal modulo the corresponding full modulus. This is a direct identity from substituting either basis vector's linear relation into its minors; it does not require a rank statement over \(\mathbb Z/q\mathbb Z\). For example modulo \(q_1\), the two vectors have first coordinate \(c\) times the second, so \(n_3=0\bmod q_1\) and \(n_2=-c n_1\bmod q_1\). Thus \(q_1\mid n_3,q_2\mid n_1,r\mid n_2\). Each lattice contributes \(O(Z^2/P)\) prime triples by the successive minima estimate.

Dyadically fix the three moduli, with positive scales \(Q_1,Q_2,Q_3\), product \(Q_*\asymp Q\), and smallest and largest scales \(u,w\), respectively. Also dyadically fix \(P\). If one normal coordinate is zero, the two moduli whose defining relations involve that coordinate divide both remaining coordinates, because both local coefficients are units. Two zero coordinates leave no actual prime triples in the plane. With just one zero, \(P\) is at least a constant times the product of those two moduli. There are \(O_\varepsilon(X^\varepsilon Q_*)\) possible lattice data in the boxes, so this branch contributes \(O_\varepsilon(X^\varepsilon Z^2 w)\).

Now all normal coordinates are nonzero. In particular \(P\gg w\). Write
\[
 (n_3,n_1,n_2)=(a_1q_1,a_2q_2,a_3r),\qquad
 (y_1,y_2,j)=(b_1q_1,b_2q_2,b_3r).
\]
The calculations with the three normal directions give
\[
            q_1q_2r\mid S:=b_2^2b_3^2a_1^3q_1+b_1^2b_3^2a_2^3q_2
                              +4b_1^2b_2^2a_3^3r.                    \tag{6}
\]
For instance modulo \(q_1\), cube \(c n_1=-n_2\), use \(4y_2^2c^3=j^2\), substitute the displayed definitions, cancel the unit \(q_2^2r^2\), and multiply by \(b_1^2\). Modulo \(q_2\) perform the symmetric operation. Modulo \(r\) cube \(c'' n_1=-n_3\), use \(y_2^2(c'')^3=y_1^2\), cancel \(q_1^2q_2^2\), and multiply by \(b_3^2\). These operations require only the unit cancellations furnished by pairwise coprimality, so they apply to complete prime powers, also at 2.

The \(b_i\) are bounded, and \(|S|\ll P^3/u^2\). If \(P\) exceeds a sufficiently small fixed multiple of \(Q_*^{1/3}u^{2/3}\), the trivial count of lattice data, divided by \(P\), costs at most \(O_\varepsilon(X^\varepsilon Q_*^{2/3}/u^{2/3})\ll X^\varepsilon V^{4/3}\), before the common factor \(Z^2\). This is the required bound in that range.

If \(P\) is at most that fixed multiple, (6) gives the exact equation \(S=0\). Write \(a_i=s a_i'\) with \(\gcd(a_1',a_2',a_3')=1\) and \(s>0\). The common factor \(s^3\) cancels from the equation. Dyadically pin \(P'\asymp\max_i |a_i'| Q_i\); it satisfies \(P'\ll P\). There are \(O(P'^3/Q_*)\) choices for the nonzero \(a_i'\) in this box. The coefficients \(L_i\) of the linear equation in \((q_1,q_2,r)\) are fixed by \(a_i',b_i\), and their gcd \(g\) is bounded by a fixed constant. Indeed each coefficient is a cube of \(a_i'\) times a nonzero fixed multiplier determined by \(b\); at each prime one of the \(a_i'\) is a unit. Also \(\max_i |L_i| Q_i\gg P'^3/w^2\), by taking an index with \(|a_i'| Q_i\asymp P'\).

To count the pairwise-coprime modulus triples for these fixed coefficients, scale the three coordinates by \(Q_i^{-1}\). The lattice of integer solutions of \(L_1q_1+L_2q_2+L_3r=0\), after scaling, has two-dimensional covolume
\[
     \frac{\|\operatorname{diag}(Q_1,Q_2,Q_3)L\|_2}{g\,Q_*}.
\]
Indeed the unscaled lattice has primitive integer normal \(L/g\), hence covolume \(\|L\|_2/g\). The scaling on its plane multiplies area by \(\det(D)\|D^{-T}L\|_2/\|L\|_2\) for \(D=\operatorname{diag}(Q_i^{-1})\). If the candidate vectors span rank at most one, they contribute at most two primitive integral vectors before scaling, by pairwise coprimality. If rank two, successive minima in a fixed ball give an upper bound by a constant times reciprocal covolume. Thus in all cases the number is
\[
                         O(1+Q_* w^2/P'^3).
\]
Multiplying by the \(a'\) count and summing \(P'\ll P\) gives \(O_\varepsilon(X^\varepsilon(P^3/Q_*+w^2))\) possible data. The scale \(s\) is deliberately not counted: fixed \(a',b\) and moduli determine the original signed data and the possible local root classes, regardless of which scale led to the normal. Dividing by \(P\) and multiplying \(Z^2\), use
\[
           P^2/Q_* \ll u^{4/3}/Q_*^{1/3} \ll V^{1/3},
           \qquad w^2/P\ll w\ll V .
\]
This finishes the small-\(P\) case. Summing rank two over the dyadic scales gives \(O_\varepsilon(X^\varepsilon Z^2V^{4/3})\).

The [rank-two divisor bound](ABC_RANK_TWO_NORM_PLANE_DIVISORS.md#counting-the-rank-two-congruences) also gives \(O_\varepsilon(X^\varepsilon(Z^2V+Z^3V^{1/2}))\), by retaining the nonzero integer quotient \(S/Q\) in (6) and factoring after one modulus is fixed.

Therefore the number of selected prime paths is \(O_\varepsilon(X^\varepsilon(Z^3+V^3+Z^2V^{4/3}))\). The six-variable recovery of the full witness, including \(\delta,d\), is unchanged. Divide by the per-discriminant lower bound \(\gg_\eta Z^3 H^2 X^{-1-4\eta}\) from (2) in the witness selection and substitute \(Z=X^{1/2+4\eta}/H\) to get the two new decaying terms of (1) on \(H\ge X^{3/14}\). For \(H<X^{3/14}\), its \(X^{5/2}/H^8\) term already dominates the classical first-moment bound \(X/H\). Absorb the split-prime exceptions and handle real signs as in the witness selection.

## Moment consequences and remaining terms

For \(k\ge1\), the tail (1), with the sum as its decaying bound, the ordinary first moment, and the pointwise cap \(h_3(d)\ll |d|^\kappa\), \(\kappa=1597/5000\), give
\[
 \sum_{0<|d|\le X}^{*} h_3(d)^k
 \ll_{k,\varepsilon}
 X^{\max(1+3(k-1)/14,\ 1+\kappa(k-2))+\varepsilon}.
\]
The rank-one arm \(X^{5/2}/H^8\) crosses the ordinary first moment at \(H=X^{3/14}\); the other decaying arm \(X^{11/6}/H^5\) crosses earlier at \(H=X^{5/24}\). For \(k\) beyond the respective two decay degrees (8 and 5), each decaying arm reaches its maximum at the pointwise cap and is dominated there by \(X/H^2\) times \(H^k\). In particular the second and fourth moment exponents are \(17/14\) and \(23/14\).

For balanced kernels at scale \(N\), the three marginal tails have the forms \(N^2/H^2+N^{11/3}/H^5+N^5/H^8\), with subpower losses. For \(H\ge N^{5/11}\), the last term is dominated by the middle one, since they cross at \(H=N^{4/9}\). The Frobenius triangle bound then gives \(O_\varepsilon(N^\varepsilon(N^3/H^3+N^{11/2}/H^{15/2}))\). The decaying term crosses the CM2 \(N^3/H^2\) bound at \(H=N^{5/11}\). Thus the balanced \(2\le k\le3\) minimum moment is \(O_{k,\varepsilon}(N^{3+(5/11)(k-2)+\varepsilon})\), and the high sector \(H\ge N^{5/9}\) admits \(k=3\) at no price. In arbitrary dyadic shapes the same price-free cubic sector holds for the common threshold \(H\ge\max(AB,AC,BC)^{5/18}\).

A direct improvement of the rank-one count \(V^3\), or a separate bound on which congruence lattices support an actual prime direction, would further improve the marginal moments. The intermediate and low polynomial torsion ranges remain open for a price-free \(k>2\) moment.

### A limitation of bounds using only rows and edge tails

The marginal tails and all shared-center second-moment bounds allow
a formal \(k>2\) obstruction. Fix \(0<a<4/9\). Partition each of
three \(N\)-element vertex sets into \(M\asymp N^a\) aligned classes
of comparable size. Give the edges joining like classes weight
\(M\), and give all others weight one. All row sums are \(O(N)\)
and all shared-center second moments are \(O(N^3)\). At a threshold
near \(M\), each edge set has \(O(N^2/M)\) exceptional entries.
Both the \(N^6/H^9\) term and the sum \(N^{11/3}/H^5+N^5/H^8\)
cover this, and monotonicity gives the bounds down to threshold
one. At a same-class triple all three weights are \(M\), contributing
\(\gg N^3 M^{k-2}\) to the minimum moment. This is a formal weight
model. An arithmetic \(k>2\) bound requires additional information
about how the three edge weights align.
