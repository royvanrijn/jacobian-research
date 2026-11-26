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
