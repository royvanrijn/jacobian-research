# A count saving for small square kernels on all three letters

For a positive integer \(T\), let
\(\operatorname{sf}_2(T)=\prod_{q:\ v_q(T)\ {\rm odd}}q\)
be its squarefree kernel. There are absolute
\(f_*,\delta,\lambda_*>0\) such that for every fixed
\(0<\lambda<\lambda_*\) the number of primitive \(a+b=c\le X\)
satisfying
\[
 \operatorname{rad}(abc)<c^\lambda,\qquad
 \prod_{T\in\{a,b,c\}}\operatorname{sf}_2(T)\le c^{f_*\lambda}
 \tag{K1}
\]
is \(O_{\lambda,\varepsilon}(X^{(1/3-\delta)\lambda+\varepsilon})\).
There is no short-letter condition and no exact valuation divisibility
requirement beyond the prescribed parity defect. In particular even
actual valuation masses may be distributed over many different
degrees. The argument uses the actual three-letter height and radical
constraints. Larger parity defects require additional weight, as in the selector below.

## A conic of controlled coefficient height

Work on a dyadic shell \(c\asymp X\). Write
\[
            a=A x^2,\qquad b=B y^2,\qquad c=C z^2
\]
with \(A,B,C\) squarefree. There are
\(X^{f_*\lambda+o_X(1)}\) triples of kernels in (K1), since
their product and its ordered prime allocation determine them.
Fix one such triple. The roots have no common divisor, and
their rational projective point lies on the nonsingular conic
\(Q(x,y,z)=A x^2+B y^2-C z^2=0\).

Set \(M=\max(A,B,C)\). Cassels' small-zero theorem gives
a nonzero integral zero \(\mathbf u\) of \(Q\) with height
\(O(M)\). Complete \(\mathbf u\) to a rational basis by two
standard vectors \(\mathbf v,\mathbf w\), and let \(B_Q\) be
the symmetric bilinear form for which \(Q(\mathbf v)=B_Q(\mathbf v,\mathbf v)\).
For \(\mathbf r=s\mathbf v+t\mathbf w\), define
\[
        (F_1,F_2,F_3)(s,t)
                = Q(\mathbf r)\mathbf u-2B_Q(\mathbf u,\mathbf r)\mathbf r .
\]
These are integral quadratic forms of coefficient height \(O(M^2)\).
They parametrize the conic by secants through the point of
\(\mathbf u\). They have no common geometric zero:
such a zero would give an isotropic two-dimensional space
spanned by \(\mathbf u,\mathbf r\) in a nondegenerate
three-dimensional quadratic space over the algebraic closure,
which is impossible. The parametrization is birational, with
inverse given off the base point by the line through \(\mathbf u\).
It is therefore an isomorphism of the projective line onto the
conic. The three coordinate lines meet the conic transversely
at distinct pairs of points, and no two such lines meet on
the conic. Hence each \(F_i\) has two simple geometric roots
and the pairs are disjoint.
[Cassels, Theorem on small zeros of an isotropic integral quadratic form](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/bounds-for-the-least-solutions-of-homogeneous-quadratic-equations/9C3C57BC5549819DF3DBF5F3BC42DCE8)

## Numerical and finite-place bounds on the conic parameter

For primitive integer parameters \((s,t)\), put
\(H=\max(|s|,|t|)\) and \(g=\gcd(F_1(s,t),F_2(s,t),F_3(s,t))\).
The resultant identities in both affine charts show
\[
 g\le M^{O(1)},\qquad
 \max_i|F_i(s,t)|=H^2 M^{O(1)}
\]
where the second equality means upper and lower bounds with
fixed polynomial factors in \(M\). The formulas for the same
root choice are \(|F_1|=g x,|F_2|=g y,|F_3|=g z\), since the
original roots are primitive. There are two complete values
of height comparable to \(X\), so \(H=X^{1/4+O(f_*\lambda)+o_X(1)}\)
and all three coordinate logarithmic heights are
\[
       \log_X |F_i(s,t)|
         =\tfrac12\log_X T_i+O(f_*\lambda)+o_X(1).
 \tag{K2}
\]
Here and below the suppressed constants are absolute.

For an outer sequence of small-budget boxes, let \(h\)
be its limiting smallest complete height exponent.
One may partition the real parameter into finitely many
charts in which every cross determinant between two
primitive parameters in a chart satisfies
\[
        |s t'-s't|\le X^{h/2+O(f_*\lambda)+o(1)}.
 \tag{K3}
\]
The \(o(1)\) includes convergence of the box height exponents.
If \(h=1\), use the usual \(2H^2\). If \(h<1\), choose one
of two real affine charts in which the denominator has
absolute value \(H\). The short coordinate polynomial is
quadratic, or linear in an affine chart with a root at
infinity. Its distinct finite roots are separated by
\(M^{-O(1)}\), from the coefficient and discriminant bounds.
Equation (K2), divided by \(H^2\), forces the chart coordinate
to lie within \(X^{-(1-h)/2+O(f_*\lambda)+o(1)}\) of one
real root. Indeed a first root proximity follows by taking
the square root of the bound on the polynomial value;
after this is much smaller than the separation, only that
root is close, and the full bound applies after dividing by
the other factor. A nonreal root and its conjugate would
contradict the separation. Record that root in the chart,
and multiply the interval length by the two denominators
to prove (K3).
[Height and real correction for a projective carrier](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#height-local-branch-weights-and-the-real-correction)

Fix a bad-prime product \(J\ne0\) containing \(2ABC\), the
leading coefficients and each nonzero homogeneous pairwise
resultant of the \(F_i\), and the derivative resultants of
their nonconstant dehomogenizations in both affine charts.
It is \(M^{O(1)}\) up to an absolute factor. If a prescribed
whole prime power \(q^e\Vert T_i\) is in a lock, work in a
\(q\)-adic projective chart whose denominator is a unit.
For \(q\nmid J\), the coordinate polynomial has a simple
root and its value has valuation \(e/2\). Recording the chart
and the branch therefore fixes the projective congruence
to depth \(e/2\).
For \(q\mid J\) the value valuation is at least
\((e-1)/2\) even accounting for the evaluated gcd. If its
valuation exceeds \(2v_q(J)\), the derivative resultant
bounds the derivative valuation by \(v_q(J)\) and Hensel
gives a root class with depth at least
\((e-1)/2-v_q(J)\). At lesser valuations discard the
precision at cost at most \(O(v_q(J)+1)\). The total lost
log weight across the selected bad support is \(O(\log J)\).
A chart and one of at most two roots per prime cost
\(X^{o_X(1)}\) combined records, by the divisor estimate.
Thus locks on any collection of letters of combined
numerical weight \(w\) give a divisibility of the cross
determinant of any two parameters in the same record by a
modulus of size
\[
                         X^{w/2-O(f_*\lambda)-o_X(1)} .
 \tag{K4}
\]
The selected prime sets on the three letters are disjoint
by primitivity. The weights are the valuations of the primitive quotient throughout.
[Exceptional support in projective residue records](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#projective-residue-records-including-exceptional-primes)

## Using all actual endpoint locks

The affordable whole-coordinate pair-lock result gives an
absolute \(\eta_0>0\) usable at all balanced endpoint
sequences. It supplies full divisors of some two letters
of combined normalized radical cost at most \(1/3-\eta_0\)
and numerical weight at least the shorter height of that
pair plus \(\eta_0\). Its height threshold is at least the
global shortest height \(h\), and it permits macroscopic
individual primes. Use these locks on a hypothetical
sequence of boxes near saturation before conditioning on
the square kernels. After conditioning, the old bound on
the possible lock values still holds, and the kernel
record adds at most \(f_*\lambda\) to the exponent.
[Affordable whole-coordinate pair locks at every balanced endpoint](ABC_WIDE_COFACTOR_CONGRUENCES.md#affordable-whole-coordinate-pair-locks-at-every-balanced-endpoint)

Choose \(f_*<\eta_0/4\). If the claimed bound failed for
arbitrarily small fixed budgets and arbitrarily small absolute
gaps, form a sequence tending to the balanced \(1/3\)
endpoint as in the existing small-budget box reduction.
The total radical exponent \(L\) satisfies \(L/\lambda\to1\).
The available divisor-record count, including all kernel
choices, has exponent at most
\((1/3-\eta_0+f_*+o(1))\lambda\).
For such a record (K4) beats (K3) by an absolute numerical
margin (multiplied by \(1/2\)), since \(f_*\lambda\to0\).
It fixes the primitive parameter, and hence the actual
triple. This contradicts saturation and yields the absolute
\(\delta,\lambda_*\). Every kernel, resultant, branch chart,
and finite truncation has been paid at its fixed budget
before the outer limit.

In the balanced minimally repaired sector with
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\) and repair
height \(H_N(c)\le Y\) for a fixed order \(N\), the height
cutoff implies \(c\le Y^{1/(N(1/3-\xi)\lambda)}\).
Applying the original count gives repaired exponent at most
\((1/3-\delta)/(N(1/3-\xi))+\varepsilon\), which is below
\(1/N\) for \(\xi<\delta\). This formulation takes the
original total radical budget to be \(\lambda\), with fixed
small enlargements absorbed when using the alternative
central-band conventions.
[One-target minimal repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## A quantitative charged-kernel selector

The fixed-record proof also yields a direct sector at
non-negligible kernel mass. Let \(L=\log_X\operatorname{rad}(abc)\).
Define the normalized kernel mass \(f=\log_X(ABC)/L\) and the
numeric weight
\[
          W=\sum_{q\mid ABC} v_q(abc)\log_X q .
\]
The sum uses entire odd-valuation prime powers, not just their
squarefree parts. Once \(A,B,C\) are paid, all exact prime
powers on their support have only subpower many possibilities,
by the fixed-support divisor estimate. Include all of them
in (K4), gaining weight \(W/2\) with no additional radical cost.
If \(W>h\) by a fixed margin, this already proves a count
exponent \(fL+o(L)\) whenever \(f<1/3\) by a fixed margin.

More generally, for the remaining even-valuation support,
a whole-prime selector of normalized mass \(s\) and numerical
weight \(w\) yields exponent \((f+s)L+o(L)\) whenever
\(W+w>h\) by a fixed margin. If its largest individual
normalized prime mass is at most \(\tau\), one can use the
explicit bound
\[
       f+\frac{(h-W)_+(1-f)}{2+h-W}+\tau
 \tag{K5}
\]
for the coefficient, up to arbitrary fixed allowances and
vanishing small-budget losses, in the range \(W\le h\).
Indeed the remaining support has mass \(1-f\) and numerical
weight \(2+h-W+o(1)\), since two complete letters are tall.
Sort its primes by their scaled valuation \(L v_q(abc)\).
To achieve a desired numerical weight \(h-W+\gamma\)
with an arbitrarily small fixed margin \(\gamma>0\), the
fractional top tail costs at most that target times
\((1-f)/(2+h-W+o(1))\). Round up the last prime to obtain
an actual full divisor, increasing the cost by at most
\(\tau\). Enumerating the selected product by its radical
cap costs exactly this normalized amount up to subpower
losses. Neither the order nor the subset need be common
across the triples. The proof applies in boxes with
uniform values of these parameters up to arbitrarily
small slack; dyadic shells and fixed parameter bins
handle their variation uniformly.

Thus with \(\tau=o(1)\) and \(W\le h\), the explicit
sufficient condition is \(W>h-1+3f\) by a fixed margin,
since this places (K5) strictly below \(1/3\).
In particular a positive ordinary short-height gap
permits a fixed positive kernel cap without any lower
bound on \(W\). With no height gap, a positive proportion
of odd radical at sufficiently high degrees can supply
the needed weight even when \(f\) is a non-negligible
fraction of the total budget. This calculation does
not guarantee a saving on every balanced parity profile.

## Two-square-kernel version and its rate obstruction

If only two letters \(U,V\) have fixed squarefree kernels,
write them as \(A x^2,B y^2\). The third letter \(T\) is
\(A x^2\pm B y^2\) up to an overall sign, so the primitive
parameter is simply \([x:y]\). It has height
\(X^{1/2+O((\log_X A)+(\log_X B))+o(1)}\).
There are now one contact of multiplicity two on
each kernel letter and two simple contacts on \(T\).
Consequently a full power divisor contributes half its
numerical weight on \(U,V\), and all its weight on \(T\).
The same finite-place precision proof applies with
coefficient error \(O(\log_X(AB))\). The cross determinant
threshold in ordinary numerical units is
\[
 \begin{cases}
   1, & \text{all three complete heights tend to one},\\
   (1+h)/2, & \text{one of \(U,V\) is short of height \(h<1\)},\\
   h, & \text{\(T\) is short of height \(h<1\)} .
 \end{cases}
\]
Indeed the real chart widths have respective exponents
zero, \(-(1-h)/2\) near the double zero, and \(-(1-h)\)
near a simple zero. Thus if a kernel letter has a
fixed positive short-height gap and the tall remaining
letter \(T\) has a prescribed positive-mass bounded
scaled-degree diffuse block, a strict saving follows
for sufficiently small joint kernel cost on **only**
\(U,V\). Start with the full lock on \(T\) of mass
\(1/3+o(1)\) and numerical weight \(1+o(1)\), then omit
a small positive block from that diffuse band: its
numerical loss can be made smaller than the fixed
determinant margin \((1-h)/2\). Pay the two kernel
radicals from the omitted radical margin. Equivalently
use actual bounded-degree-band primes small enough to
round the omitted block. If \(T\) is short by a fixed
margin below height \(1/2\), the analogous count starts
with one tall kernel letter of numerical branch
weight \(1/2\) and omits a diffuse block there if
available.

This specific two-kernel argument has a height-one
residual obstruction at the level of profiles. In the
balanced case with all heights tending to one and
each scaled degree supported on a narrow diffuse
interval about \(3\), all \(U,V\) half-rates are less
than every \(T\) full rate. The entire radical budget
\(1/3\) on \(T\) gives exactly the threshold weight one,
and replacing any of it by equally cheap locks from
\(U,V\) decreases weight. Thus no such threshold lock
gives a saving on this profile, even with trivial
kernels on \(U,V\). This does not assert that the
profile is realized by primitive triples, or rule out
a different use of the conic and the third radical cap.

## A small squarefree kernel suffices for abc finiteness

Fix any \(\eta>0\). The full \(abc\) conjecture is equivalent
to the following restriction: for every fixed \(0<\mu<1\)
there are only finitely many positive triples
\[
  x^2+D y^2=z^2,\qquad
  \gcd(x^2,Dy^2,z^2)=1,\qquad
  D\ \text{squarefree},\quad D\le (z^2)^\eta,\qquad
  \operatorname{rad}(Dxyz)<(z^2)^\mu.
  \tag{K2}
\]
The same fixed \(\eta\) works for all \(\mu\). This is a
finiteness reduction with a changed radical threshold.

To see the quantitative transfer, take a positive primitive
triple \(a+b=c\), put \(R=\operatorname{rad}(abc)\), and set
\[
 g_0=\gcd(a-b,c)\in\{1,2\},\quad s=(a-b)/g_0,\quad
 t=c/g_0,\quad E=t^2-s^2=4ab/g_0^2 .
\]
Primitivity gives the asserted gcd and \(\gcd(s,t)=1\).
For a fixed integer \(k\ge1\), write
\((s+\sqrt{-E})^k=P+Q\sqrt{-E}\) with \(P,Q\in\mathbb Z\).
The norm and the difference of the two conjugate powers give
\[
       P^2+EQ^2=t^{2k},\qquad
       |P|\le t^k,\qquad |Q|\le kt^{k-1}.
\]
For the last inequality, divide the difference of the two
kth powers by the difference of their bases, then bound
each of the \(k\) summands using \(|s\pm\sqrt{-E}|=t\).
After setting \(t=1\), the expressions \(P,Q\) are nonzero
polynomials in \(s\), with respective values \(1,k\) at
\(s=1\). Homogeneity shows that \(PQ=0\) occurs at finitely
many ratios \(s/t\) for fixed \(k\). Each ratio recovers the
single primitive ordered triple from \(a/c=(1+s/t)/2\).

The even binomial terms give
\(P\equiv2^{k-1}s^k\bmod t\), since \(-E\equiv s^2\bmod t\).
Hence \(g=\gcd(P,t^k)\le2^{k(k-1)}\), using
\(\gcd(P,t^k)\mid\gcd(P,t)^k\).
When \(PQ\ne0\), divide the norm relation by \(g^2\).
Its two square letters have coprime roots \(x=|P|/g\) and
\(z=t^k/g\). The other letter is a positive integer with
squarefree kernel \(D=\operatorname{sf}_2(E)=\operatorname{sf}_2(ab)\),
since multiplication by the rational square \(Q^2/g^2\)
preserves valuation parity. This yields the primitive triple
in (K2), with \(C=z^2\) satisfying
\[
  C\asymp_k c^{2k},\qquad D\le R,\qquad
  \operatorname{rad}(Dxyz)\le R|P Q|
                      \ll_k R c^{2k-1}.
\]
For the radical bound, primes in \(tE\) already belong to
\(abc\), including 2, and division by \(g^2\) adds no prime
to the three integer letters.

If \(abc\) fails, choose \(0<\lambda<1\) and an unbounded
sequence with \(R<c^\lambda\). Fix \(k\) with \(1/(2k)<\eta\).
The images have unbounded \(C\) and \(D\le C^\eta\) for all
sufficiently large \(c\). Their radicals are
\(O_k(C^{1-(1-\lambda)/(2k)})\). Any fixed
\(\mu\) strictly between \(1-(1-\lambda)/(2k)\) and \(1\)
therefore gives infinitely many exceptions to (K2).
The other implication is restriction to the triples in (K2).
Ellenberg's congruence reduction also uses a height
amplification built from \((a-b)/c\). The norm power here
keeps the third letter in the initial squareclass of \(ab\)
throughout the amplification.
[Ellenberg, Lemma 4 and Theorem 3](https://people.math.wisc.edu/~ellenberg/CongruenceABC.pdf).

The restriction has an exact binary sextic form. Take coprime
positive integers \(u,v\) with \(u^2>Dv^2\) and put
\[
  J=u^2-Dv^2,\qquad K=2uv,\qquad S=u^2+Dv^2,\qquad
  b=\gcd(J,K,S).
\]
Then \((x,y,z)=(J,K,S)/b\) gives precisely the primitive
positive conic points. The inverse parameter is
\([u:v]=[z+x:y]\). At an odd prime dividing \(b\), the
relations and \(\gcd(u,v)=1\) force the prime to divide \(u,D\)
and be prime to \(v\). Its valuation on \(S\) is one because
\(D\) is squarefree. Thus \(b\) is supported on the primes
of \(2D\). The prime 2 always appears in the primitive
additive triple, so every prime removed by division by \(b\)
still lies in \(\operatorname{rad}(Dxyz)\). In particular
\[
 \operatorname{rad}(Dxyz)
   =\operatorname{rad}\bigl(2Duv(u^2-Dv^2)(u^2+Dv^2)\bigr),
 \qquad z=S/b .
\]
The remaining finiteness estimate can be posed directly
for those sextic values, with the exact moving-kernel
constraint \(D\le(S/b)^{2\eta}\) and threshold \((S/b)^{2\mu}\).

## Primitive Pythagorean triples and a fixed exact-power hypotenuse

Fix \(k\ge1\). The full \(abc\) conjecture is equivalent to
finiteness at every strict radical threshold for primitive
positive triples \(X^2+Y^2=Z^{2k}\). Here the threshold is
\(\operatorname{rad}(XYZ)<(Z^{2k})^\mu\) for fixed
\(0<\mu<1\). The case \(k=1\) fixes the kernel in (K2) to \(D=1\).

Given a primitive \(a+b=c\), choose \(u,v\) of opposite
parity from the three letters. If \(c\) is odd, take the two
summands. If \(c\) is even, take \(c\) and one summand.
Thus the remaining letter is \(u+v\) or \(|u-v|\),
\(\gcd(u,v)=1\), and \(\max(u,v)\asymp c\).
Put \(S=u^2+v^2\) and write
\[
   (u+i v)^{2k}=P+iQ,\qquad P^2+Q^2=S^{2k}.
\]
The algebras \(\mathbb F_p[i]\) for odd \(p\) are reduced
algebras. If \(p\) divides both \(P,Q\), the reduction of
\(u+i v\) is nilpotent and hence zero, contrary to
coprimality of \(u,v\). Since \(S\) is odd, 2 cannot divide
both \(P,Q\) either. Outside the finitely many rational
ratios for which \(PQ=0\), the displayed triple is primitive
and positive after taking absolute values. Finiteness of
those ratios follows from the nonzero homogeneous real and
imaginary part polynomials, and a ratio determines the
coprime \(u,v\).

Each of \(u,v,u-v,u+v\) divides the homogeneous polynomial
\(P Q\). At the corresponding points, \((u+i v)^2\) is
either real or purely imaginary. The product of the four
distinct linear factors therefore divides \(P Q\) over
\(\mathbb Z\), by Gauss's lemma. Its degree is \(4k\), so after
removing the three original letters, the remaining factor
has degree \(4k-3\). Since \(S=O(c^2)\), this proves
\[
  C=S^{2k}\asymp_k c^{4k},\qquad
  \operatorname{rad}(PQS)\ll_k
     \operatorname{rad}(abc)c^{4k-1}.
\]
An unbounded original sequence below \(c^\lambda\),
with fixed \(\lambda<1\), yields unbounded primitive
exceptions at any fixed \(\mu\) strictly between
\(1-(1-\lambda)/(4k)\) and \(1\). The reverse implication
is restriction of full \(abc\).

The degree \(4k\) is minimal for a rational map whose
three specified fibres have local degrees divisible by
\(2,2,2k\), respectively. With degree \(d\), those
fibres alone contribute at least
\((2-1/(2k))d\) to the ramification divisor, so
Riemann–Hurwitz gives \(d\ge4k\). The map
\(P^2/S^{2k}\) has degree \(4k\), with difference from
one equal to \(-Q^2/S^{2k}\), and attains the bound.

The exponent conversion
between the two full conjecture statements is exact:
\(\operatorname{rad}(abc)<c^\lambda\) is
\(c>\operatorname{rad}(abc)^{1+\varepsilon}\) for
\(\lambda=1/(1+\varepsilon)\), which ranges over \((0,1)\).

## A fixed icosahedral power signature suffices for abc

The full \(abc\) conjecture is equivalent to finiteness at
every strict radical threshold for primitive signed relations
\[
       X^2-Y^3=1728 Z^5,\qquad
       XYZ\ne0,\qquad \gcd(X^2,Y^3,1728 Z^5)=1 .
       \tag{K3}
\]
The threshold is \(\operatorname{rad}(6XYZ)<C^\mu\) for
\(0<\mu<1\), where \(C=\max(X^2,|Y|^3,1728|Z|^5)\).
The signs allow the three absolute letters to be reordered
as a positive primitive additive triple.

Use the following forms from Klein's icosahedral identity.
In the display, \(s=U^5\) and \(t=V^5\).
\[
\begin{aligned}
 F(U,V)&=UV(s^2+11st-t^2),\\
 G(U,V)&=s^4-228s^3t+494s^2t^2+228st^3+t^4,\\
 T(U,V)&=s^6+522s^5t-10005s^4t^2-10005s^2t^4-522st^5+t^6,\\
 T^2-G^3&=1728F^5 .
\end{aligned}
\]
The quotient map \(1728F^5/T^2\) has degree 60. Its fibres
over \(0,1,\infty\) have uniform orders \(5,3,2\)
respectively. The forms have no common geometric zero in
pairs, and their separate roots in \(\mathbf P^1\) are simple.
[Darmon–Granville, §7.4 and Klein's displayed identity](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf).

Choose \(u,v\) from the input as in the Pythagorean reduction,
then put \(U=u^2-v^2\) and \(V=2uv\). Swap \(u,v\) if
\(u\equiv 2v\pmod 5\). This ensures \(u\not\equiv2v\pmod5\)
since \(u/v=2\) changes to \(3\) in that case. It preserves
the three input linear factors up to sign. Now \(U,V\)
are coprime, with \(U\) odd and \(V\) even. Also
\(U\not\equiv2V\pmod5\), since
\(u^2-4uv-v^2\equiv (u-2v)^2\pmod5\).

The evaluated relation is primitive. At 2, \(G\) is odd.
At 3, \(G\equiv(s^2+t^2)^2\) and is nonzero for coprime
\(U,V\). At any other prime dividing both \(F\) and \(G\),
neither \(U\) nor \(V\) can vanish, because \(G\) on either
axis is a pure power with coefficient one. With \(w=s/t\),
the polynomial \(H(w)=w^2+11w-1\) must vanish. The
remainder of the displayed quartic in \(s,t\) at
\((s,t)=(w,1)\) modulo \(H\) is
\(3125(1-11w)=3125w^2\bmod H\).
Thus a common prime of \(F,G\) would have to be 5.
There \(H(w)=(w-2)^2\) and \(w=U/V\bmod5\), already
excluded. Identity (K3) then gives the required gcd.
Finite exceptional rational input ratios at zeros of \(F,G,T\)
can be omitted, since the composition with
\((u^2-v^2)/(2uv)\) is nonconstant.

The output height \(C\) is comparable to \(c^{120}\).
Indeed \(\max(|U|,|V|)\asymp c^2\) by
\(U^2+V^2=(u^2+v^2)^2\). On normalizing that maximum
to one, the maximum of the three absolute output forms of
common degree 60 is bounded below, by compactness and
absence of a common projective zero.

The radical estimate uses the ramification of the quadratic
substitution. The polynomial \(T(U,V)\) is divisible by
\(U^2+V^2\), as is seen by setting \(U=\pm iV\) in the
displayed sextic in \(s,t\). This factor becomes \(S^2\).
Consequently \(\operatorname{rad}(T(U,V))\ll c^{58}\),
using \(S=O(c^2)\) and the remaining degree-28 form
in \(U,V\). For \(F\), the factor \(UV=2uv(u-v)(u+v)\)
contains the three original letters. The remaining linear
factor has size \(O(c)\), and the residual degree-ten form
has size \(O(c^{20})\). Together with \(G=O(c^{40})\),
these yield
\[
       \operatorname{rad}(6TGF)\ll
         \operatorname{rad}(abc)c^{119}.
\]
If \(\operatorname{rad}(abc)<c^\lambda\) along an
unbounded sequence with fixed \(\lambda<1\), the output
gives infinitely many strict exceptions in (K3) at any
\(\mu\) between \(1-(1-\lambda)/120\) and \(1\).
Conversely, the absolute letters of a solution to (K3)
form a primitive positive additive triple with radical
\(\operatorname{rad}(6XYZ)\). Each triple determines
only finitely many assignments, roots and signs.

## Minimum degree for anchoring three rational points

Consider nonconstant rational maps from \(\mathbf P^1\)
over \(\mathbb Q\) for which the local degrees over
\(0,1,\infty\) are divisible by \(5,3,2\), respectively.
Among maps sending all three input points \(0,1,\infty\)
into that set, the minimum degree is \(120\).

Every such map \(f\) has degree divisible by \(60\).
Over \(\overline{\mathbb Q}\), pull back the Klein quotient
\(f_0=1728F^5/T^2\) along \(f\) and normalize. The
projection to the source of \(f\) is unramified,
since locally the quotient is a power of order \(e\),
and the pullback local order is a multiple of \(e\).
In local coordinates the relation is \(w^e=t^m\) with
\(e\mid m\), and its normalized branches project with
degree one locally. Away from \(0,1,\infty\) the quotient
was already unramified. Riemann–Hurwitz excludes an
unramified cover of \(\mathbf P^1\) of degree greater than
one, so each component gives a lift of \(f\) through
\(f_0\). In degree \(60\), the lift is a Möbius map.

In any degree-60 model over \(\mathbb Q\), at most two
points above \(0,1,\infty\) are rational. First suppose
there is a rational order-five point. The cyclic
stabilizer is Galois stable and its other fixed point
is rational as well. Put these two points at \(0,\infty\)
by a rational coordinate change. Using a deck
transformation of the Klein map to send the first image
to 0, the lift is now scaling by \(\rho\). The remaining
order-five divisor in these coordinates is given by
\(\rho^{10}z^{10}+11\rho^5 z^5-1\). This polynomial
has rational coefficients because its divisor is rational
and its constant coefficient is normalized. In
particular \(\rho^5\in\mathbb Q^\times\). A further
rational ramification point would give a rational root
\(w=\rho^5z^5\) of one of the residual polynomials
in the display in \(s,t\), with \(t=1\).
The quadratic for \(F\) has nonsquare discriminant \(125\).
For \(G,T\) the rational root test leaves only \(\pm1\),
where both polynomials are nonzero. Thus this case has
exactly two rational points above the three fibres.

Suppose next that there is no rational order-five point.
Identify \(\operatorname{Aut}(A_5)\) with \(S_5\) acting
by conjugation. This identification follows from the
action on the five Sylow two-subgroups, which are the
Klein four groups fixing the respective letters. The
action on those groups by \(A_5\) is faithful, so an
automorphism fixing all five induces the identity on
\(A_5\). Conjugation by \(S_5\) realizes all permutations.

Choose a Galois automorphism \(\sigma\) acting by
squaring on fifth roots of unity and fixing third roots
of unity. A Möbius element of order five has a pair of
inverse multipliers at its two fixed points. Under
\(\sigma\), this pair is squared, giving a different
pair. Thus \(\sigma\) acts by an outer automorphism on
the geometric deck group, represented by an odd
permutation \(h\in S_5\).
Complex conjugation preserves the multiplier pair and
acts by an inner automorphism, because the two
\(A_5\)-classes of five-cycles are exactly the two pairs
of possible multipliers. The implementing element
\(\kappa\) is a nontrivial involution. A trivial action
would realize \(A_5\) in \(\operatorname{PGL}_2(\mathbb R)\).
Such an action preserves orientation on the real
projective circle, by simplicity, and a finite group
of orientation-preserving circle homeomorphisms is
cyclic. Indeed a finite-order homeomorphism fixing a
point is the identity, and a free orbit gives a faithful
action by cyclic permutations.

The involution \(\kappa\) has no real fixed point.
Otherwise, in real coordinates sending its two fixed
points to \(0,\infty\), it acts by \(z\mapsto-z\).
The antiholomorphic map obtained by composing
\(\kappa\) with complex conjugation centralizes all
deck transformations. It is now \(z\mapsto-\bar z\),
which becomes ordinary conjugation after a complex
scaling. Again this would realize \(A_5\) over \(\mathbb R\).

The centralizer of \(\kappa\) in \(A_5\) is a Klein
four group. The centralizer of \(h\) in \(A_5\) is
\(S_3,C_2\), or \(C_3\), according as \(h\) is a
transposition, a four-cycle, or a disjoint transposition
and three-cycle. Their intersection has at most one
nontrivial involution. Each rational order-two point
defines a rational involution by its stabilizer, so
there are at most two rational order-two points.
Each rational order-three point defines a cyclic
subgroup fixed pointwise by \(h\), since \(\sigma\)
fixes the multipliers there. Each displayed centralizer
has at most one order-three subgroup, whose two
fixed points bound the number of rational order-three
points by two.

If both types of rational point existed, \(h\) would
have to be a transposition. The rational involution
would then invert the distinguished order-three
subgroup in its \(S_3\) centralizer. Complex conjugation
also inverts that subgroup, by its multipliers at a
rational point, so \(\kappa\) and the rational involution
belong to its normalizer \(S_3\) in \(A_5\). They commute
because the latter involution is rational. Distinct
involutions in \(S_3\) do not commute, so the rational
involution is \(\kappa\). Its lack of real fixed points
contradicts a rational order-two point. This completes
the two-point bound in degree 60.

Finally \(f_0((z^2-1)/(2z))\) has degree 120 and
sends \(0,1,\infty\) to the order-five fibre. The two
branch values of the quadratic map are \(\pm i\),
which are in the order-two fibre of \(f_0\). Thus this
composition has no branch value outside \(0,1,\infty\),
and meets all three local divisibility conditions.

## Counting primitive Pythagorean exceptions

For fixed \(0<\mu\le1\), the number of primitive positive
\(x^2+y^2=z^2\) with \(z^2\le X\) and
\(\operatorname{rad}(xyz)\le (z^2)^\mu\) is
\[
 O_{\mu,\varepsilon}\!\left(X^{
   \min\{\mu/3,\ \max(1/4,(3\mu-1)/7)\}+\varepsilon}\right).
 \tag{K4}
\]
At the unit budget the exponent is \(2/7\).

Take \(x\) odd by ordering the two legs. The coprime
integers \((z+x)/2,(z-x)/2\) have square product
\((y/2)^2\). Writing them as \(u^2,v^2\) gives the
usual coprime opposite-parity parameters \(u>v>0\).
On a shell \(u\asymp H\) the five values
\[
     u,\quad v,\quad u+v,\quad u-v,\quad S=u^2+v^2
\]
are pairwise coprime. The two sums or differences are
odd, and a common odd prime in two of the values would
divide both \(u,v\). The product of their radicals is
exactly \(\operatorname{rad}(xyz)\), since one of
\(u,v\) is even. The shell therefore has total radical
budget \(O(H^{4\mu})\).

In polynomial-size ranges, the number of integers with
any specified radical is \(O_\varepsilon(H^\varepsilon)\).
For instance, the Rankin bound with exponent
\(\varepsilon/2\) has Euler product over the prescribed
primes, of size \(H^{o(1)}\). The latter follows by
ordering those primes and using \(p_j\ge j+1\) and
\(\omega=O(\log H)\). Also \(S\) has
\(O_\varepsilon(H^\varepsilon)\) sum-of-two-squares
representations by the usual Gaussian divisor bound.
Thus \(\operatorname{rad}(S)\le H^r\) gives
\(O_\varepsilon(H^{r+\varepsilon})\) pairs.
Set \(r=4\mu/3\). In the complementary case the two
cheapest radicals among the four linear values have
product \(O(H^{4\mu/3})\). Enumerating these two
radicals and then their values costs
\(O_\varepsilon(H^{4\mu/3+\varepsilon})\), and two
distinct linear values determine \(u,v\).

For the second bound, take
\(r=\max(1,(12\mu-4)/7)\). When
\(\operatorname{rad}(S)>H^r\), choose the three
cheapest linear radicals. Their product is
\(O(H^{3\mu-3r/4})\). Each of the four possible
choices gives a primitive positive additive triple
of height comparable to \(H\). Put \(L_\pm=u\pm v\).
The relations are \(L_-+v=u\), \(u+v=L_+\),
\(L_-+2v=L_+\), and \(L_-+L_+=2u\).
Doubling a letter costs at most a factor two in the
radical. Knowing the relation and the triple recovers
\(u,v\) up to a bounded number of choices.

For primitive additive triples with maximum letter
\(\asymp H\), letter product \(\asymp V\) and radical
product at most \(R\), Kane's upper-bound proof gives
\(O_\varepsilon(H^\varepsilon(R/H+V^{1/3}))\).
In the proof of Proposition 7, before
replacing the three radicals by powers of their
individual letters, the dyadic bound is
\(H^\varepsilon(R/H+\prod_i v_i)\), with
\(v_i=\prod_{p^2\mid A_i}p\) on the three letters.
Use the split in the proof of Theorem 2 with the
actual letter product. Write \(L\) for the product of
the square roots of the three largest square divisors.
When \(L\le V^{1/3}\), the product of the \(v_i\) is
at most \(V^{1/3}\). Above \(V^{1/3}\), the two terms
in Proposition 10, apart from \(H^\varepsilon\), are
\(O(V/L^2)\) and \(O(\sqrt{V/L})\), each bounded
by \(O(V^{1/3})\). There are only logarithmically
many dyadic intervals for \(V\).
[Kane, Proposition 7, Proposition 10 and proof of Theorem 2](https://cseweb.ucsd.edu/~dakane/ABC2.pdf).

Here \(V\ll H^3\) and \(3\mu-3r/4-1\le r\), so this bound gives
\(O_\varepsilon(H^{r+\varepsilon})\) in the second
case as well. Summing the shells with \(H\ll X^{1/4}\)
proves (K4).

If one linear value is at most \(u^\theta\) for fixed
\(0\le\theta\le1\), the second count also has exponent
\[
 \frac14\max\{(2+\theta)/3,\ (8\mu+\theta-3)/5\}.
\]
For this bound let \(r\) be the displayed maximum
without the factor \(1/4\). Keep the short value
and the two cheapest of the other three linear
radicals. Their combined exponent is at most
\((8\mu-2r+\theta)/3\) in the shell \(H\) when
\(\operatorname{rad}(S)>H^r\). Now \(V\ll H^{2+\theta}\),
and the height-sensitive bound costs \(H^{r+\varepsilon}\).
At \(\mu=1\) this improves the \(2/7\) exponent when
\(\theta<5/7\).

## Four linear factors suffice for abc finiteness

For coprime integers \(u>v>0\) of opposite parity put
\[
 F(u,v)=uv(u-v)(u+v),\qquad R(u,v)=\operatorname{rad}(F(u,v)).
\]
The four factors are pairwise coprime. In particular \(u-v\) and
\(u+v\) are odd, and their gcd divides both \(2u\) and \(2v\).
The abc conjecture is equivalent to finiteness of the pairs with
\(R(u,v)<u^{2-\delta}\) for every fixed \(0<\delta<2\).

In one direction, apply abc to the primitive triple
\(v^2+(u^2-v^2)=u^2\), whose radical is exactly \(R(u,v)\).
Conversely, suppose infinitely many primitive \(a+b=c\) satisfy
\(\operatorname{rad}(abc)<c^\lambda\) for some fixed \(0<\lambda<1\).
If \(c\) is odd, take \(u=\max(a,b)\), \(v=\min(a,b)\). If \(c\)
is even, take \(u=c\), \(v=a\). In both cases \(u\asymp c\),
with the required coprimality and parity. The product \(F(u,v)\)
is \(abc\) times one positive integer at most \(2c\), so
\(R(u,v)\ll u^{1+\lambda}\). These unbounded pairs violate the
four-factor threshold for every fixed \(0<\delta<1-\lambda\).
The four points \(0,1,-1,\infty\) and the squaring map also occur
in Smirnov's arithmetic Hurwitz formulation.
[Smirnov, Hurwitz inequalities for number fields](https://www.mathnet.ru/eng/aa316).

## Counting small radicals of the four linear factors

For fixed \(0<\beta\le2\), the number of these pairs with \(u\le X\)
and \(R(u,v)\le u^\beta\) is
\[
 O_{\beta,\varepsilon}\!\left(X^{e(\beta)+\varepsilon}\right),
 \qquad
 e(\beta)=\min\left\{\frac{\beta}{2},
 \max\left(\frac{4\beta}{9},\frac{\beta}{4}+\frac7{18},\frac45\right)\right\}.
 \tag{K5}
\]
In particular \(e(2)=8/9\). For comparison, using only the three
cheapest radicals and applying Theorem 1.2 to their additive
relation gives \(15/16\) at \(\beta=2\), since their combined
budget is \(3\beta/4\) and doubling costs at most a factor two.
The proof of (K5) uses the lattice of all four divisibilities and
the Fourier bounds from all four additive relations.
[Bernert, Browning, Lichtman and Teräväinen, Theorem 1.2 and Proposition 3.1](https://arxiv.org/html/2410.12234v2).

### Valuation slots and a four-factor lattice bound

Use dyadic ranges \(u\asymp H\) and \(\min(v,u-v)\asymp W\).
The product of the four linear values \(L_i\) is \(\asymp H^3W\).
After fixing which value is short, use \((u,v)\) or \((u,u-v)\)
as integral coordinates, with box dimensions \(O(H),O(W)\).

Fix a large integer \(m\ge4\). Write
\(L_i=q_i\prod_{j=1}^m x_{ij}^j\), where \(x_{ij}\) is the
product of primes of exact valuation \(j\) and \(q_i\) contains
the full prime powers of valuation above \(m\). Since
\(\operatorname{rad}(q_i)\ll H^{1/(m+1)}\), the fixed-radical
bound used for (K4) gives \(O_{m,\eta}(H^{4/(m+1)+\eta})\)
possibilities for the ordered \(q_i\). Fix these and the dyadic
scales \(H^{a_{ij}}\) of all \(x_{ij}\). Write
\(a_i=\sum_j a_{ij}\), \(p_j=\sum_i a_{ij}\), \(B=\sum_j p_j\).
Then \(B\le\beta+o(1)\), \(\sum_j jp_j\le4+o(1)\), and
\(3a_{i3}\le1+o(1)\).

Choose a subset \(I\) of slots to leave free and fix all other
slots. Put \(S=\sum_{(i,j)\in I}a_{ij}\) and
\(T=\sum_{(i,j)\in I}(j-1)a_{ij}\). For each \(i\) the fixed
variables prescribe a divisor \(M_i=q_i\prod_{j:(i,j)\notin I}x_{ij}^j\)
of \(L_i\). The \(M_i\) are pairwise coprime whenever the pattern
has a solution. The Chinese remainder theorem gives a lattice of
index \(\prod_i M_i\), since each linear form is primitive modulo
its own modulus. The index is \(\asymp_m H^3W H^{-S-T}\).

The number of primitive points for this fixed pattern is
\(O(1+HW/\prod_i M_i)\). To see this, scale the two box dimensions
to bounded size. If the lattice points in the box are collinear
through the origin, there are at most two primitive integral
preimages. Otherwise two independent points of bounded size
generate a sublattice of bounded covering radius. The full lattice
has at most that covering radius, so the disjoint Voronoi cells of
the counted points lie in a box of bounded area. Each cell has area
\(\prod_i M_i/(HW)\). Multiplying by \(O_m(H^{B-S})\) choices of
fixed slots gives a count \(O_m(H^{B-S}+H^{B-2+T})\).

### Fourier bounds and exponent elimination

Each of the four triples of linear forms determines \((u,v)\)
and participates in one of the relations
\[
 (u-v)+v=u,\qquad u+v=(u+v),\qquad
 (u-v)+2v=(u+v),\qquad (u-v)+(u+v)=2u.
\]
Use the second and fourth moment estimates in the proof of
Proposition 3.1, equations (3.1)–(3.2), on each letter of a relation.
Its coordinate product is \(H^{a_i}\). For every \(2\le j\le m\),
the moments are at most \(H^{a_i+\eta}\) and
\(H^{3a_i-a_{ij}+\eta}\), up to constants. The fourth moment saves
the product over all multiples of \(j\), which includes slot \(j\).
The fixed coefficient \(q_i\) and any doubling only dilate the
integral argument by an integer. Orthogonality, Cauchy–Schwarz and
Hölder give the exponent \((2/3)\sum_i a_i-(1/6)\sum_i a_{ij}+\eta\)
for the three letters in the relation. Dropping the constraints on
the omitted form can only increase the count. Averaging the four
bounds gives \(B/2-p_j/8+\eta\).
[Bernert, Browning, Lichtman and Teräväinen, proof of Proposition 3.1](https://arxiv.org/html/2410.12234v2#S3).

It suffices to rule out limiting exponents \(D>e(\beta)\) for
fixed boxes. Indeed, if the stated estimate fails for some
\(\varepsilon>0\), take \(m\) large and \(\eta\) small so that
the loss from \(q_i\) and the logarithmically many scales and shells
is less than \(\varepsilon/2\). Along a sequence of offending
boxes take a convergent subsequence of the bounded log sizes.
The geometric and Fourier estimates yield
\[
 D\le\max(B-S,B-2+T),\qquad p_j\le4B-8D\quad(j\ge2),
\]
for every choice of free slots, with \(B\le\beta\),
\(\sum_j jp_j\le4\), and \(a_{i3}\le1/3\). The arbitrary
\(H^\eta\) loss in the moment bound disappears in this limit.

If \(p_2\ge D\), Fourier gives \(D\le4B/9\). Otherwise leave all
degree-one and degree-two slots free. Since \(B\le2\), the second
geometric exponent is at most \(p_2<D\), giving
\(p_1+p_2\le B-D\). The height inequality then implies
\[
 4B-4\le3p_1+2p_2+p_3,
 \qquad p_2\le4-B-3D+p_3\le4+3B-11D.
\]
If \(p_2+2/3\ge D\), this gives \(D\le B/4+7/18\).
In the remaining case also leave the largest degree-three slot
free, of size \(a\in[p_3/4,1/3]\). The second geometric exponent
is at most \(p_2+2a<D\), giving \(p_1+p_2+p_3/4\le B-D\).
Therefore
\[
 4B-4\le3(p_1+p_2)+p_3-p_2
 \le3(B-D)+p_3/4\le4B-5D,
\]
so \(D\le4/5\). Fourier also gives \(D\le B/2\). This proves
(K5), uniformly including shells with one small linear value.

## The lattice of repeated-prime divisors

For a fixed admissible pair, take \(E_i=L_i(u,v)/\operatorname{rad}(L_i(u,v))\)
and \(E=\prod_i E_i\). Let \(\Lambda\) consist of integral \((r,t)\)
with \(E_i\mid L_i(r,t)\) for each of the four forms. Then
\[
 \Lambda=\{(r,t)\in\mathbb Z^2:E\mid ut-vr\}
         =\mathbb Z(u,v)+E\mathbb Z(r_0,t_0),\qquad ut_0-vr_0=1.
 \tag{K6}
\]
Modulo \(E_u\), the coefficient \(v\) is a unit and \(u=0\).
Modulo \(E_v\), the coefficient \(u\) is a unit and \(v=0\).
Modulo each of \(E_+\) and \(E_-\), the coefficient \(u\) is a
unit and \(v\equiv\mp u\), respectively. Thus each divisibility
condition is equivalent to \(E_i\mid ut-vr\). Their product is
permitted by pairwise coprimality. The final equality in (K6)
uses the determinant-one integral basis given by Bézout's identity.

Take \(s=\min(v,u-v)\) and let \(\ell\) be the corresponding
linear form. In the real coordinates \((r/u,\ell(r,t)/s)\),
the pair is \((1,1)\) and the lattice has determinant
\[
 \Delta=\frac{E}{us}
       =\frac{\max(v,u-v)(u+v)}{R(u,v)}\asymp \frac{u^2}{R(u,v)}.
\]
The least sup norm \(\rho\) of a lattice vector independent of
\((1,1)\) satisfies
\[
 \frac{\Delta}{2}\le\rho\le\frac{\Delta+1}{2}.
\]
Indeed the difference of the two coordinates of an independent
vector is a nonzero integral multiple of \(\Delta\) in absolute
value. The second basis vector in (K6) has coordinate difference
of absolute value \(\Delta\). Translate by a multiple of \((1,1)\)
to put the coordinate mean within \(1/2\) of zero.

Full abc is therefore also equivalent to a bound
\(\rho\ll_\varepsilon u^\varepsilon\) for each \(\varepsilon>0\).
In a hypothetical strict exception with \(R(u,v)<u^{2-\delta}\),
the lattice has only one direction in every bounded normalized box
once \(u\) is large. Bounding the volume contribution in (K5)
still leaves such primitive points. A pointwise argument using
these congruences must constrain which repeated-prime divisors can
occur at the input pair. Ordinary lattice counting also permits
the isolated primitive directions described by (K6).

### Continued fractions for a prescribed divisor pattern

Prescribe four pairwise coprime positive divisors \(M_i\) and fix
which linear value is \(s=\min(v,u-v)\). In the coordinates
\((u,s)\), write \(u=M_0 x\), \(s=M_1 y\), where \(M_0,M_1\)
correspond to \(u,s\). Put \(Q=M_2M_3\) for the other two moduli.
If \(s=v\), the other forms give \(s\equiv u,-u\) modulo their
respective moduli \(M_{u-v},M_{u+v}\). If \(s=u-v\), they give
\(s\equiv u,2u\) modulo \(M_v,M_{u+v}\). Hence they prescribe
\(y\equiv kx\pmod Q\) for a uniquely determined \(0\le k<Q\).
The residue is obtained by multiplying the indicated slopes by
\(M_0 M_1^{-1}\) modulo each remaining modulus.

Suppose \(M=\prod_i M_i>2us\). Then \(2xy<Q\), so for
\(n=(kx-y)/Q\) we have
\[
 0<\frac{k}{Q}-\frac{n}{x}=\frac{y}{Qx}<\frac{1}{2x^2}.
\]
The fraction \(n/x\) is reduced because \(\gcd(x,y)=1\).
By Legendre's criterion it is a convergent of \(k/Q\). To include
the terminating rational convention, one can use both finite
continued-fraction expansions. The criterion follows by extending
the two expansions of \(n/x\). Their intervals end at mediants
with the respective preceding fractions, on opposite sides and
at distances \(1/[x(x+b)]\) with preceding denominator \(b<x\).
For an integer approximation the upper interval starts with that
integer. Here the exact final convergent cannot occur since \(y>0\).

Thus the Euclidean algorithm lists every admissible pair in this
range for the prescribed pattern. Recover \(y=kx-nQ\) and check
positivity, the choice of short value, coprimality and parity.
There are \(O(\log(2Q))\) convergents. In the expansion giving
\(n/x\), the next partial quotient is \(M/(us)+O(1)\). This follows
from the exact convergent error, since its inverse after
multiplication by \(x^2\) is the next complete quotient plus the
ratio of the preceding denominator to \(x\).

For \(M_i=E_i\) the next partial quotient is therefore
\(\Delta+O(1)\) and captures every sufficiently large strict
exception. Conversely, a family of prescribed patterns and
reconstructed admissible pairs with unbounded \(u\) disproves abc
if \(M/\operatorname{rad}(M)\gg u^{1+\delta}s\) for some fixed
\(\delta>0\). Indeed \(M/\operatorname{rad}(M)\) divides \(E\),
so \(R(u,v)\ll u^{2-\delta}\). Short lattice vectors for prescribed
prime powers are also used in Dokchitser's search on the
three-letter equation. The congruence here uses all four factors
and depends on the actual short value.
[Dokchitser, LLL & ABC](https://arxiv.org/abs/math/0307322).

### A fixed set of primes cannot supply the excess

Fix any finite prime set \(S_0\). Let \(F_{S_0}\) denote the
full \(S_0\)-part of \(F(u,v)\). The Subspace Theorem gives
\[
 F_{S_0}\ll_{S_0,\varepsilon} us\,u^\varepsilon
 \qquad(\varepsilon>0).
 \tag{K7}
\]
Use the two independent real forms \(u,\ell(u,v)\) and, at each
\(p\in S_0\), the single form \(L_i\) divisible by \(p\). If none
is divisible by \(p\), use \(u\). There are finitely many choices
of these form systems and the product in Theorem 4.1 is exactly
\(us/F_{S_0}\). For a fixed system, a product at most
\(u^{-\varepsilon}\) confines the pairs to finitely many proper
rational subspaces of \(\mathbb Q^2\). Each line contains at most
two primitive integer vectors.
[Schlickewei, Theorem 4.1, with two real forms and one form per finite place](https://matwbn.icm.edu.pl/ksiazki/aa/aa31/aa3148.pdf).

In particular, if \(R(u,v)<u^{2-\delta}\), the part of \(E\)
supported outside any fixed \(S_0\) is
\(\gg_{S_0,\varepsilon}u^{\delta-\varepsilon}\). This follows
from \(E\gg us\,u^\delta\) and \(E_{S_0}\le F_{S_0}\).
An unbounded disproof family through the divisor reconstruction
must obtain a power-sized contribution from primes outside every
fixed set. The pointwise gap is control of the normalized excess
from such varying prime supports.

## Fixed polynomial prescriptions for the repeated parts

Index the four forms by \(u,v,u-v,u+v\), and write
\(s=\min(v,u-v)\). Fix four nonzero polynomials
\(M_i(T)\in\mathbb Z[T]\). For every integer \(n\ge2\) with all
\(M_i(n)\ne0\), and every admissible pair satisfying
\[
 |M_i(n)|\mid \frac{L_i(u,v)}{\operatorname{rad}(L_i(u,v))}
 \qquad\text{for all }i,
\]
one has the effective bound
\[
                    \prod_i |M_i(n)|\ll_{(M_i)} us . \tag{K8}
\]
In particular, a fixed polynomial prescription for all four
exact repeated parts gives \(R(u,v)\gg_{(M_i)}u^2\). The polynomials
can have values supported on primes that grow with \(n\). A finite
list of polynomial prescriptions, including a fixed set of residue
classes with one prescription on each class, gives the same bound.
The order \(us\) is attained by \(u=2n^k,v=1\) for \(n\) a power
of two, with \(M_u(T)=T^k\) and the other moduli equal to one.

### Finite polynomial directions for the specialized lattices

First assume only \(|M_i(n)|\mid L_i(u,v)\) and put
\(M(T)=\prod_i M_i(T)\). If there is a nonconstant common polynomial
factor of two moduli, its primitive integral representative divides
both integral polynomials by Gauss's lemma. For large \(n\) its
absolute value exceeds one, contrary to coprimality of distinct
linear values. Thus, for the large-\(n\) argument, the moduli are
pairwise coprime in \(\mathbb Q[T]\). Finitely many remaining \(n\)
give bounded \(|M(n)|/(us)\).

In the chart \((u,s)\), write the forms as \(u,s,u-s,u+s\) or
\(u,s,u-s,2u-s\), with the moduli reordered accordingly. Over
\(\mathbb Q[T]\), the simultaneous divisibilities define a free
rank-two module with a basis of determinant \(M(T)\). Explicitly,
write the first two coordinates as \(M_0 x,M_1 y\). The other two
congruences give \(y\equiv kx\pmod{M_2M_3}\) by the polynomial Chinese
remainder theorem. The columns \((M_0,M_1k)\) and
\((0,M_1M_2M_3)\) are the desired basis. A modulus that is constant
imposes the trivial congruence over \(\mathbb Q[T]\).

There is a constant \(C\) and a finite list of nonzero polynomial
directions, depending only on the moduli and chart, such that
\(|M(n)|>Cus\) implies \((u,s)\) is proportional to the evaluation
of a listed direction. For the proof, let
\(\theta=\log_n(u/s)\). In the range \(|M(n)|>us\), this lies in a
fixed bounded interval because \(u\le |M(n)|\). Order the coordinate
monomials by degree on the first coordinate and degree plus
\(\theta\) on the second, resolving ties in favor of the first.
The order is constant on each open interval between successive
integers and on each integer singleton. In each of these finitely
many cases, reduce the polynomial basis by subtracting a monomial
multiple of one column from the other whenever their leading
positions agree. Choose the multiple to cancel the leading term
of the column with the larger degree in that position. Each step
strictly lowers a leading monomial in the fixed order, so the
process terminates with distinct leading positions. This is the
two-column form of reduction to weak Popov form.
[Mulders and Storjohann, §2](https://cs.uwaterloo.ca/~astorjoh/jscpopov.pdf).

Multiply both columns by a fixed positive integer \(D\) to clear
denominators of their entries and all four polynomial quotients.
The resulting columns \(b_1,b_2\), evaluated at \(n\), span a
sublattice of index \(D^2\) in the integer congruence lattice. Indeed
the evaluated determinant has absolute value \(D^2|M(n)|\), and
the full lattice has determinant \(|M(n)|\) by coprimality of the
evaluated moduli and the Chinese remainder theorem.

Use \(\|(x,y)\|_\theta=\max(|x|,n^\theta|y|)\) after specialization.
The opposite leading positions give
\[
 n^\theta D^2|M(n)|\asymp
             \|b_1(n)\|_\theta\,\|b_2(n)\|_\theta .
\]
In the polynomial determinant, the product of the leading entries
has strictly higher degree than the crossed product. If one column
has equal weighted degrees, the other has a strict inequality by
the tie rule. The specialization estimate is uniform within each
of the finitely many orders because the absolute value of every
nonzero fixed polynomial is comparable to \(n^{\deg}\) for large
\(n\). Writing \((u,s)=\alpha_1b_1(n)+\alpha_2b_2(n)\), one has
\(\alpha_j\in D^{-2}\mathbb Z\). Cramer's rule and
\(\|(u,s)\|_\theta=u\) give \(|\alpha_j|\ll u/\|b_j(n)\|_\theta\).
If both are nonzero, the last two estimates imply \(|M(n)|\ll us\).
Taking the basis columns over the finitely many orders proves the
claim for arbitrary reconstructed integer pairs. General parametric
LLL addresses a fixed metric. The argument here explicitly allows
the aspect ratio \(u/s\) to depend on the reconstructed pair.
[Bogart, Goodrick and Woods, Theorems 1.2 and 1.4](https://arxiv.org/abs/1909.04762).

### Specialization and repeated polynomial factors

Normalize a listed direction to coprime \(U,V\in\mathbb Z[T]\) in
the \((u,s)\) chart, writing \(V(n)\) for the short coordinate before
primitive normalization. The gcd of the evaluated coordinates divides
a fixed nonzero integer by the polynomial Bézout identity. Thus
\(us\asymp |U(n)V(n)|\), and \(|M_i(n)|\) divides the corresponding
form value \(A_i(n)=L_i(U(n),V(n))\) in this chart. Directions with a
vanishing form have no admissible pairs. A constant direction gives a
bounded primitive pair and contributes bounded prescribed divisors.

For all sufficiently large \(n\), the evaluation divisibilities force
\(M_i(T)\mid A_i(T)\) over \(\mathbb Q[T]\). In polynomial division,
clear denominators and compare the remainder of smaller degree
with \(|M_i(n)|\), which divides its integer value. A nonzero
remainder is impossible for large \(n\). Write \(Y_i=A_i/M_i\).
For a primitive irreducible factor \(f\) of \(M_i\) with
\(f\nmid Y_i\), divide a constant multiple of \(Y_i\) by \(f\)
to obtain a nonzero integral remainder of degree below \(\deg f\).
Under the repeated-part divisibilities in (K8), \(|M_i(n)|
\operatorname{rad}|M_i(n)|\) divides \(A_i(n)\), so every prime
divisor of \(f(n)\) divides the remainder at \(n\). This uses Gauss's
lemma to pass from \(f\mid M_i\) over \(\mathbb Q[T]\) to integral
divisibility. Since \(f\) and the remainder are coprime over
\(\mathbb Q[T]\), another polynomial Bézout identity confines all
prime divisors of \(f(n)\) to a fixed finite set.

If one such \(f\) has degree at least two, it has at least two
distinct roots. If there are two distinct such linear factors,
their product has two distinct roots and values supported in the
union of the two fixed sets. In either case the effective
\(S\)-part estimate for a polynomial with at least two distinct
roots bounds \(n\). Thus for large \(n\) at most one linear factor
is missing from the corresponding \(Y_i\).
[Bugeaud, Evertse and Győry, Theorem 2.2](https://arxiv.org/abs/1708.08290).

For a nonconstant direction, the Wronskian \(W=UV'-U'V\) is
nonzero and has degree at most \(\deg U+\deg V-1\). For any form
\(A_i\), complete it by an independent linear form in \(U,V\).
Their Wronskian is a nonzero constant multiple of \(W\). Hence a
factor of multiplicity \(e\) in \(A_i\) divides \(W\) to at least
multiplicity \(e-1\). Each factor of \(M_i\) has gained a
multiplicity in \(A_i=M_iY_i\), except for the possible single
missing linear factor. Since the \(M_i\) are pairwise coprime,
\(M\mid W\) or \(M\mid fW\) over \(\mathbb Q[T]\), with \(f\)
linear. In either case \(\deg M\le\deg U+\deg V\), giving (K8)
on the finite list of directions. Denominator clearing, polynomial
comparisons and the effective \(S\)-part bound make the constant
computable from the prescriptions.

There is also a bound for prescriptions that use ordinary
divisors. Let \(Q(T)\) be the product of the distinct primitive
irreducible factors of \(M(T)\). If \(|M(n)|>Cus\), use the same
list of directions and \(M/Q\mid W\) to obtain, for large \(n\),
\[
 \frac{|M(n)|}{\operatorname{rad}|M(n)|\,us}
 \ll_{(M_i)}\frac{|Q(n)|}{n\operatorname{rad}|Q(n)|}. \tag{K9}
\]
Here \(\operatorname{rad}|Q(n)|\mid\operatorname{rad}|M(n)|\).
A disproof that uses only the excess guaranteed by a fixed
polynomial prescription of ordinary divisors would therefore
need a separate power-sized repeated part in values of its
squarefree polynomial \(Q\). The repeated-part prescription in
(K8) cannot produce a fixed positive deficit. These bounds leave
open prescriptions with varying polynomials or growing exponents,
as well as repeated parts outside the prescribed moduli.

## Uniform bounds for moving polynomial divisors

Let \(P_i(T)\in\mathbb Z[T]\) be monic, pairwise coprime over
\(\mathbb Q[T]\), and put \(P=\prod_i P_i\). Write \(d=\deg P\ge1\),
\(Q=\operatorname{rad}_{\mathbb Q[T]}P\) in monic normalization,
and \(B=P/Q\). Suppose all complex roots of \(P\) have modulus at
most \(H\ge1\). If \(n\ge 2H\) is an integer and \(P_i(n)\mid L_i(u,v)\)
on an admissible pair, then
\[
 \frac{B(n)}{us}\le \frac{(16(1+H))^{32d^2}}{n},
 \qquad s=\min(v,u-v).                                      \tag{K10}
\]
The estimate permits the polynomials, their coefficients and
their multiplicities to vary with \(n\). In particular,
\[
 \frac{P(n)}{us\,\operatorname{rad}P(n)}
 \le (16(1+H))^{32d^2}\,\frac{Q(n)}{n\operatorname{rad}Q(n)}.
\]
The loss depends on the total degree and root bound. This makes
the specialization in (K9) uniform without assuming the normalized
congruence lattice has only one direction.

### Polynomial bases with bounded coefficients

Work in the appropriate \((u,s)\) chart from (K8), so the four
forms are \(u,s,u-s,u+s\) or \(u,s,u-s,2u-s\). Set
\(\theta=\log_n(u/s)\) and \(D=\deg B\le d-1\). If \(\theta>2d\),
then \(us=s^2n^\theta>n^{2d}\) and \(B(n)\le(3n/2)^D\) give
(K10) directly. In the remaining range use the leading-position
order with weights \(0,\theta\) and ties favoring the first
coordinate, as in the fixed-prescription argument.

For each of the two leading positions, choose a vector in the
polynomial congruence module with the smallest possible leading
monomial at that position. Its degree there is at most \(d\),
since the corresponding coordinate vector times \(P\) belongs to
the module. The degree in the other coordinate is at most \(3d\).
The two vectors form a basis over \(\mathbb Q[T]\). Cancel the
leading monomial of any module vector by a monomial multiple of
the basis vector in the same position and iterate. These are
elementary shifted polynomial basis arguments of the type used
for modular polynomial equations.
[Neiger, Problems 1 and 2](https://arxiv.org/abs/1602.00710).

The basis vectors \(b_1,b_2\) can be chosen in \(\mathbb Z[T]^2\)
with coefficients of absolute value at most
\(A=(16(1+H))^{12d^2}\). To see this, set the leading coefficient
at each selected position to one and use the remaining allowed
coefficients as unknowns. The congruences give \(d\) linear
equations from their remainders, together with the leading
coefficient equation. Their integer entries have absolute value
at most \((16(1+H))^{4d}\). For a monic modulus of degree \(q\ge1\),
the quotient of \(T^k\), \(k\le3d\), is the polynomial part of
\(T^{k-q}\prod_{j=1}^q(1-\alpha_j/T)^{-1}\). Its coefficients are
bounded by \(2^{4d}(1+H)^{3d}\), using the expansion in the roots.
The modulus has coefficient sum at most \((1+H)^q\). Subtraction
then bounds each remainder coefficient, with the linear form
coefficients contributing at most a factor two. A nonsingular
minor of rank at most \(d+1\), Cramer's rule and denominator
clearing give the stated \(A\). Scaling a column by a nonzero
integer does not change its leading monomial or its status as
a basis over \(\mathbb Q[T]\).

Write \(w_j\) for the weighted leading degree of \(b_j\). The
determinant is \(cP(T)\) for some nonzero integer \(c\) by the
polynomial Chinese remainder theorem. Distinct leading positions
give \(w_1+w_2=d+\theta\). The basis columns satisfy the integral
congruences after evaluation, since division by a monic polynomial
preserves integral quotients. Their span has index \(|c|\) in the
evaluated congruence lattice. The latter has determinant \(P(n)\)
because its four moduli are coprime as divisors of the four
pairwise coprime linear values.

### The Wronskian bounds every admissible short point

In the norm \(\|(x,y)\|_\theta=\max(|x|,n^\theta|y|)\), each
evaluated column has norm at most \(2A n^{w_j}\). Also
\(P(n)\ge(n/2)^d\) and \(\|(u,s)\|_\theta=u\). In the expansion
\((u,s)=\lambda_1 b_1(n)+\lambda_2 b_2(n)\), the coefficients belong
to \(c^{-1}\mathbb Z\). Cramer's rule gives, for every nonzero
\(\lambda_j\),
\[
 \frac1{|c|}\le |\lambda_j|
 \le \frac{2^{d+2} A u}{|c| n^{w_j}}.
\]
Thus a nonzero coefficient at a weight \(w_j\) forces
\(u\ge n^{w_j}/(2^{d+2}A)\).

If both coefficients are nonzero, choose the larger weight.
It satisfies \(2w_j-\theta\ge d\ge D+1\). If only one coefficient
is nonzero, write its column as \((U,V)\). For a nonzero
Wronskian \(UV'-U'V\), each factor of \(P_i\) of multiplicity
\(m\) occurs in that Wronskian to order at least \(m-1\).
Indeed complete \(L_i(U,V)\) by an independent form in \(U,V\)
and differentiate. The two Wronskians differ by a nonzero scalar.
The four polynomial moduli are coprime, so \(B\) divides the
Wronskian and \(D\le\deg U+\deg V-1\le 2w_j-\theta-1\).
If the Wronskian vanishes, the column has constant direction.
Every form on it is a nonzero constant multiple of \(V\),
because its evaluation gives an admissible pair. Thus all \(P_i\)
divide \(V\), and \(w_j\ge d+\theta\) again suffices.

In all cases, \(us=u^2/n^\theta\ge n^{D+1}/(2^{2d+4}A^2)\).
Combining this with \(B(n)\le(3n/2)^D\) proves (K10).
The radical form follows from \(P=BQ\) and
\(\operatorname{rad}Q(n)\mid\operatorname{rad}P(n)\).

### Nearby prime bases with growing valuation orders

For each form, suppose a product \(G_i\) of powers \(q^k\),
\(k\ge1\), divides its repeated part \(E_i\). Assume the bases
are distinct primes or pairwise coprime squarefree integers
\(q=n+a\) with \(|a|\le H\). Replace a selected power \(q^k\) by
the polynomial \((T+a)^{k+1}\) in the corresponding \(P_i\).
The evaluated product divides \(L_i\), since
\(G_i\operatorname{rad}G_i\mid L_i\). In (K10) the value \(B(n)\)
is exactly \(\prod_i G_i\). The total polynomial degree is at most
twice the sum of the selected excess exponents \(k\).

In an unbounded sequence with \(n\to\infty\), if
\(d\log(2+H)=o(\log n)\), then
\[
                       \prod_i G_i \le \frac{us}{n}\,u^{o(1)} .
\]
Indeed one modulus has degree at least \(d/4\), so
\(u\ge (n/2)^{d/4}/2\). The logarithm of the degree and root
loss in (K10) is therefore \(o(\log u)\). This regime includes
fixed offsets with unbounded total excess order \(o(\log n)\).
At any fixed strict deficit \(R(u,v)<u^{2-\delta}\), the ratio
\(E/\prod_i G_i\) must then be at least \(n u^{\delta-o(1)}\).
The possible excess at other primes or in higher valuations
remains a pointwise obligation.

## Fixed valuation orders and the short factor

Fix positive integers \(A_i,k_i\) on the four forms and consider
admissible pairs with \(L_i(u,v)=A_i x_i^{k_i}\), \(x_i\ge1\). Put
\[
 \sigma=\sum_i\frac1{k_i},\qquad
 J=\prod_i\frac{x_i}{\operatorname{rad}(x_i)},\qquad
 s=\min(v,u-v)=L_j(u,v).
\]
The same results apply with fixed \(k_i\) and all valuations of
\(L_i\) divisible by \(k_i\) outside a fixed finite set of primes.
The \(k_i\)-power-free coefficients then range over finitely many
choices. There is no restriction on where primes in the \(x_i\)
lie or how far apart they are.

* If \(\sigma<2\), there are only finitely many admissible pairs.
* In all other cases,
  \[
             R(u,v)J\gg_{(A_i),(k_i),\varepsilon} u^{2-\varepsilon}
             \qquad(\varepsilon>0),                         \tag{K11}
  \]
  with one residual case. This is the case of exactly two
  distinct orders greater than one, both on long forms, with
  \(k_j=1\) on the short form.

When at least three orders exceed one, or exactly two do and
they are equal, (K11) holds with \(J\) replaced by \(J_{\overline S}\),
the part of \(J\) outside any fixed prime set \(S\). The same
strengthening holds when all orders are one.

Consequently a strict radical deficit in the covered profiles
requires a power-sized repeated part in the roots \(x_i\). In
particular, (K11) gives full finiteness at every fixed strict
threshold on these profiles if \(J\le u^{o(1)}\). In the strengthened
cases it suffices for each prime valuation outside \(S\) to be
either zero or the order of its letter. Valuations at the fixed
primes can be arbitrary.

### The fixed cover and the real proximity estimate

Write \(t=v/u\). The four marked points on \(\mathbb P^1\) are
\(\infty,0,1,-1\). Give the corresponding points orders \(k_i\),
omitting the order-one points from the ramification data. If
at least three orders exceed one, or if exactly two do and they
are equal, the resulting orbifold has a finite cover by a smooth
projective curve. In complex terms it is a quotient of a Riemann
surface by a finite group, with local orders exactly \(k_i\).
Riemann existence and descent give a cover
\(\pi:C\to\mathbb P^1\) over a number field, ramified only as
prescribed. For degree \(d=\deg\pi\), Riemann–Hurwitz gives
\[
                       2-2g(C)=d(\sigma-2).
\]
[Behrend and Noohi, Proposition 5.7](https://webspace.maths.qmul.ac.uk/b.noohi/papers/Uniformization.pdf)
· [Darmon and Granville, §3.1](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf).

All admissible \(t\) in this fixed profile lift to \(C(K)\) for
one fixed number field \(K\). For the descent, take the root stack
on the marked projective line, so the cover is finite étale over
this stack. Outside a fixed finite set of primes the marked
sections are disjoint and the cover and its finite étale property
spread out. Enlarge this set by the primes dividing the \(A_i\).
The point \([v:u]\) has intersection order \(v_p(L_i(u,v))\)
with the corresponding section, by primitivity. These divisors
are divisible by their orders, so their quotient divisors give
the required root line bundles with sections on the base ring.
After base extension the intersection orders remain divisible.
Pulling back the cover gives a finite étale scheme of degree \(d\)
over the ring of integers with the fixed primes inverted. Its
field factors have bounded degree and ramification in that fixed
set. Hermite–Minkowski gives one finite extension containing all
the lifts. This is the specialization argument in the
Darmon–Granville method.
[Darmon and Granville, Proposition 3.2 and proof of Theorem 2](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf).

For \(\sigma<2\), there are at least three nontrivial orders and
\(C\) has genus greater than one. Faltings gives finiteness of
the lifts, and \(t\) determines the primitive pair.
[Faltings, Satz 7](https://doi.org/10.1007/BF01388432).

For the other genera, the primitive height of \(t\) is \(u\).
At each archimedean place, a lift of \(t\) close to the short
mark approaches some preimage with order \(k_j\). At a place
above a prime \(p\in S\), use the mark whose form is divisible
by \(p\), if there is one. Coprimality ensures uniqueness.
Near each preimage, the morphism to that mark has local order
\(k_i\), including order one at an unramified mark. Thus the
weighted sum of local logarithmic proximities is
\[
 \frac{\log(u/s)}{k_j}+
     \sum_{p\in S}\sum_i\frac{v_p(L_i)}{k_i}\log p+O(1).
\]
The weights are the local degrees divided by \([K:\mathbb Q]\)
and sum to one above each place of \(\mathbb Q\). The indicated
proximities differ by bounded constants from the local expressions
also at the finitely many bad primes, by the local form of a
map of curves and compactness. There are finitely many choices
of target preimages. The identity cover gives the same argument
when all orders are one.

When \(g(C)=0\), enlarge \(K\) once to obtain a parameter on
\(C\). Its logarithmic height is \((\log u)/d+O(1)\). Apply the
number-field Roth bound with the chosen target at each of the
fixed places. In the dimension-one Subspace Theorem, take
one linear form vanishing at the target and one that is a unit
in its local projective chart. This gives the upper bound
\((2+\eta)\log u/d+O_{\eta}(1)\) for the weighted sum. Taking
all archimedean places preserves the coefficient \(2\) with
absolute logarithmic heights. Here \(2/d=\sigma-2\).
[Evertse and Ferretti, equation (1.1) and §1.4](https://arxiv.org/abs/math/0408381).

When \(g(C)=1\), enlarge \(K\) once to give the curve an origin
and to make the preimages rational over \(K\). The fixed-place
approximation estimate on an elliptic curve gives upper bound
\(\eta\log u+O_\eta(1)\) for that sum for every \(\eta>0\). To
apply it at a preimage, translate that point to the origin.
The local canonical height is the logarithmic proximity to the
origin up to a bounded term. The global canonical height is
\(O(1+h(t))\) by comparison of ample heights on the fixed curve.
Apply the fixed-place height estimate to each translate, using
that the local heights at the other fixed places are bounded
below.
Here \(\sigma=2\).
[Gross and Silverman, Proposition 9](https://msp.org/pjm/1995/167-2/pjm-v167-n2-p05-s.pdf).

Since the four letters are pairwise coprime,
\(RJ\asymp_{(A_i)}\prod_i x_i\). On the other hand
\(\prod_i x_i\asymp_{(A_i)}u^\sigma(s/u)^{1/k_j}\), as the
three long forms are comparable to \(u\). Also
\(\log(J/J_{\overline S})\le\sum_{p\in S,i}v_p(L_i)\log p/k_i\),
so the proximity bounds give the strengthened (K11). If exactly
one order exceeds one, or two
orders do and the short form has one of them, the elementary
bound \(\sigma-1/k_j\ge2\) already suffices.

### Two unequal powers and an unpowered short form

Suppose the two powered long forms have distinct orders
\(k,\ell\ge2\), and the short form is unpowered. In the \((u,s)\)
charts the long forms are \(u,u-s,u+s\) or \(u,u-s,2u-s\).
The two powered forms therefore give a difference
\(|A x^k-B y^\ell|\asymp s\), after multiplying one coefficient
by two if necessary. The compared terms are \(\asymp u\) and
have gcd at most two. In this case the exact radical accounting
above gives
\[
                  RJ\asymp u^{1+1/k+1/\ell}s .               \tag{K12}
\]
Thus the remaining instance of (K11) is the generalized Hall
lower exponent \(1-1/k-1/\ell\) for the gap between the selected
powers. For \((k,\ell)=(2,3)\) this is the weak Hall exponent
\(1/6\) in the height of the compared terms.
[Laniewski and Müller, Weak Hall's conjecture, approximation gains, and continued fractions](https://arxiv.org/abs/2609.35252).

Conversely, a family of coprime positive powers with fixed
coefficients and unbounded maximum \(u\) would disprove abc if
\(0<h=|A x^k-B y^\ell|<u^{1-1/k-1/\ell-\delta}\) for one
fixed \(\delta>0\). For opposite parities, use the larger and
smaller terms as \(u,v\). For two odd terms, use \(u\) and \(v=h\).
Both constructions are admissible, with the smaller term
appearing as \(v\) or \(u-v\) and with \(h\) as the short form
for large \(u\). The radical is at most a fixed constant times
\(u^{1+1/k+1/\ell}h\), which produces a strict four-factor
deficit. This identifies a genuine different-powers route for
the residual profile, alongside the remaining repeated parts
inside roots and the profiles without fixed orders.
