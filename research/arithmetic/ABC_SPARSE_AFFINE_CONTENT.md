# An affine count for sparse fibres in short target intervals

Fix two full divisor records \(F\mid U,G\mid V\) in an
anatomic box at \(n=X^L\), with primitive signed equation
\(T=Fu+Gv\) and positive residual values \(u,v\).
Here \(F,G\) include the fixed signs. Write \(w_F,w_G\) for
the numerical weights of the divisors and \(H_u,H_v\) for
the numerical height caps of the residuals in base \(X\).
Let \(c_u,c_v\) be the residual nonatomic scaled-degree
masses. Suppose the sample in the fixed pair record is a
matching in \(u,v\) and its targets lie in one interval of width
\(X^{h+o(1)}\). A target can supply a full divisor \(M\) of
positive numerical weight at least \(w\), with at most
\(n^{\rho+o(1)}\) choices of \(M\) and radical at most
\(n^{B_0+o(1)}\) for a fixed \(B_0\).
Assume
\[
              w>\min(H_u+h-w_G,H_v+h-w_F).
\]
Then in the limiting box conventions the exponent in the fixed
pair record is at most \(\rho+\min(c_u,c_v)\). The matching
condition may be omitted if \(w+\min(w_F,w_G)>h\).

## A fixed divisor supplies a primitive shifted line

Freeze \(M\). All differences from any base point of the sample satisfy \(F\Delta u+G\Delta v\equiv0\bmod M\). Thus the determinant of any two difference vectors is divisible by \(M\), because \(F,G\) are units. Multiplying by \(G\) rewrites it as the determinant of the two \((\Delta u,\Delta T)\) vectors. The analogous multiplication by \(F\) uses \((\Delta T,\Delta v)\). Its absolute value is consequently at most \(2X^{\min(H_u+h-w_G,H_v+h-w_F)+o(1)}\). The affine threshold gives collinearity.

A horizontal or vertical line meets the matching in at most one point. A line through the origin contains at most one positive primitive pair. Otherwise use a primitive equation
\[
 a u+b v=c,\qquad \gcd(a,b)=1,\qquad abc\ne0 .
\]
Its fixed integer coefficients have polynomial height by taking the primitive line through two polynomial-height points. Set \(D=Fb-Ga\). The identities
\[
        bT=D u+Gc,\qquad aT=-D v+Fc
\]
give the exact equality
\[
              \gcd(D,M)=\gcd(c,M).
 \tag{SL2}
\]
Indeed for \(p^e\Vert M\), all \(u,v,F,G\) are units, so the first identity makes the truncated valuations \(\min(v_p(D),e)\) and \(\min(v_p(c),e)\) equal. This includes \(D=0\). Taking differences of the identities shows that all \(u\) in the matching, and also all \(v\), occupy one residue class modulo \(M/\gcd(M,c)\).

To strip coefficients on a nondegenerate line, fix the parts
of \(u,v\) supported on \(abc\) at subpower cost, uniformly
in the line's polynomial height. Include the stripped factors
in the coefficients and divide by their fixed common gcd.
For \(p\mid M\) this gcd has valuation zero, since \(u,v\)
are units and \(\gcd(a,b)=1\). Thus the part of the intercept
supported on \(M\) survives. The remaining variable roots
are coprime both to the two new coefficients and to the
remaining fixed shift.
[Uniform fixed-support count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## A stronger diffuse cap without a content-weight threshold

The count in a matching has normalized limiting exponent at most
\[
                  \rho+\min(c_u,c_v).
 \tag{SL4}
\]
Even the matching assumption can be discarded if \(w+\min(w_F,w_G)>h\). In that case a horizontal or vertical line yields a singleton by the width of the target and its fixed divisor: for example if \(u\) is fixed, two values of \(v\) differ by a nonzero multiple of \(M\), hence their targets differ by at least \(|G|M\). The line through the origin is still a singleton by primitivity.

Here is the proof for a fixed \(M\) after collinearity on a nondegenerate line. In the usual iterated box limits, pass to a subsequence on which the weight of \(M/\gcd(c,M)\) converges. If it has positive limit, that integer has radical at most \(n^{B_0+o(1)}\). At a suitably small fixed positive normalized cutoff one can discard low prime orders, preserving fixed positive numerical weight. The resulting \(k\)-full divisor (with \(kL\) bounded away from zero) separately puts \(u\) and \(v\) in unit residue classes, by (SL2) and its difference consequences. For each of these one-letter congruences, fix the roots outside shrinking windows around the zero atom and finitely many positive atoms. The cost approaches the nonatomic mass of the selected letter when the finite list is enlarged. Group the retained roots and use the amplified unit congruence (GW), which releases the exact value at arbitrarily small additional normalized cost. The line determines the other free value. This bounds the fixed-record exponent both by \(c_u\) and by \(c_v\).

If the weight of \(M/\gcd(c,M)\) tends to zero, the intercept \(c\) retains essentially the entire original weight of \(M\) on its support. Strip and fix the root parts supported on the line coefficients as above, then divide the fixed coefficient gcd. Since \(u,v\) are units at the primes of \(M\) and \(a,b\) are coprime, the large-weight divisor of the fixed shift survives. Its high-valuation part is again a positive-weight \(k\)-full modulus. This case does not require close residual caps. If both remaining free degree measures have positive nonatomic mass, either can supply an additional full modulus of arbitrarily small normalized cost from a compact positive band. Freeze such a modulus on one letter; the other satisfies a fixed unit congruence since all variable roots are prime to the line coefficients after support stripping. Apply the same one-letter grouped counting argument and then swap the letters. This bounds the fixed-record exponent by the smaller nonatomic mass. If just one remaining free measure has positive nonatomic mass, use its cheap full modulus to release the wholly atomic other letter, at arbitrarily small cost. If neither does, use the fixed shift modulus and the amplified two-letter grouped congruence to release their exact primitive ratio and therefore both values. The support stripping only decreases each nonatomic mass, by domination of weak measures. This proves the fixed-record bound, uniformly by the subsequence argument, then summation over at most \(n^{\rho+o(1)}\) target-lock records proves (SL4).

The conditional statement only needs that a sample in the original fixed pair record lives in a single target interval of the specified width. That interval may be chosen differently for each original fixed pair record. It need not have small absolute height or a prescribed centre.

### A fixed-cost sector beyond the elementary interval and projective bounds

Use limiting root masses \(1/3\) on all three letters. The shorter free letter \(U\) has full height \(63/100\); \(V,T\) both have full height one. Restrict \(T\) to an interval of numerical width exponent \(h=7/20\) within each prescribed pair record. Use the following genuine degree profiles, described here as formal box data, without asserting the existence of saturating triples.

- On \(U\), fix a full block at scaled order four, root cap \(59/400\), weight \(59/100\). The remaining numerical weight is \(1/25\). Make its residual degree law the sum of mass \(169/2000\) uniform on \([9/100,11/100]\) and an atom of mass \(38/375\) at \(1893/6080\). The diffuse part has weight \(169/20000\) and the atom weight \(631/20000\), which sum to \(1/25\). Every residual degree is below \(8/25\).
- On \(V\), take mass \(9/100\) uniformly on \([399/100,401/100]\), all of it fixed as the factor \(G\). It has weight \(9/25=36/100\). The remaining mass \(73/300\) is uniform on an interval of half-width \(1/8\) centred at \(192/73\). Its numerical weight is \(16/25\), and every residual degree is below \(69/25=2.76\). Thus the whole \(V\)-law is nonatomic.
- On \(T\), take mass \(1/4\) uniformly on \([61/20,67/20]\) and an atom of mass \(1/12\) at \(12/5\). These parts have weights \(4/5\) and \(1/5\) respectively. Choose from the diffuse band the target lock of root cap \(\rho=1/100\), of weight in \([61/2000,67/2000]\) up to vanishing errors.

The pair record cost is \(f+g=19/80\), leaving unsaved allowance \(r=23/240\). The affine threshold is \(H_u+h-w_G=1/25+7/20-9/25=3/100\), strictly below the specified target modulus weight. That target weight is strictly below \(H_u=1/25\), so one cannot infer positive intercept content from that height comparison alone. The horizontal and vertical cases are singleton even with \(M=1\), since \(w_F,w_G>h\). The after-record exponent from (SL4) is at most \(169/2000\), smaller than \(r-\rho=103/1200\). Restoring all three records gives
\[
          d\le 19/80+1/100+169/2000=83/250=1/3-1/750 .
\]

For comparison, any target lock with cost at most \(r\) weighs at most \((67/20)r<7/20\), so direct reconstruction in the target interval fails. The full \(V\) values lie in an interval of width at most \(X^{63/100+o(1)}\) by the short \(U\) height and the specified \(T\) interval. With its initial fixed \(G\), a further affordable residual \(V\) lock weighs at most \((69/25)r<27/100=63/100-9/25\), so direct interval reconstruction on this free letter also fails. A further residual \(U\) lock weighs at most \((8/25)r<1/25=63/100-59/100\), so cannot reconstruct that letter from its ordinary height interval. A projective lock of total affordable root cap \(1/3\) has numerical weight at most \((401/100)/3<163/100\), the needed shortest-pair exponent. These tests allow arbitrary whole prime subfactors; the example does not depend on an indivisible target atom. The three pair nonatomic sums are at least \(1/3\), since \(c_u+c_T=169/2000+1/4>1/3\) and \(V\) has nonatomic mass \(1/3\). The sector distinguishes the new affine diffuse-cap count from these elementary reconstruction tests. No assertion is made that it escapes every other specialized profile estimate.
