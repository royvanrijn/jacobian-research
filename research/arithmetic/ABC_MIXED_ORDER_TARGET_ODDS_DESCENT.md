# Two-target parity with an odd power class of primitive odds

Fix an odd prime \(p\), a height parameter \(X\ge2\), and \(H>0\).
Use positive integer triples
\[
 T_i=\sigma U_i+\tau V_i,\qquad
 \gcd(U_i,V_i)=1,\qquad
 U_i,V_i,T_i\le X^H
\]
with common fixed signs \(\sigma,\tau\in\{-1,1\}\). All three original radical caps may be imposed throughout. For positive rationals let \([x]_q\) denote the multiplicative class modulo rational \(q\)-th powers. Let \(\operatorname{sf}_2(x)\) be the positive squarefree representative of \([x]_2\), including primes of odd negative valuation.

## A point count with two quadratic characters

**MO1.** Fix positive rational \(r\ne1\), the relative odds
\[
 r=\frac{U_2/V_2}{U_1/V_1},
\]
fix a class \([U_1/V_1]_p\) of representative \(c>0\), and fix
\[
 d=\operatorname{sf}_2\big((T_1/V_1)(T_2/V_2)\big),\qquad
 d'=\operatorname{sf}_2(dr).
\]
Use the representative \(c=\prod_{\ell}\ell^{v_\ell(U_1/V_1)\bmod p}\), with residues in \([0,p-1]\). Set \(K=\mathbb Q(\zeta_p)\) and
\[
 h_p(u)=\big|\mathrm{Cl}\big(K(\sqrt u)\big)[p]\big|
 \quad (u\text{ positive squarefree}).
\]
For some constant \(C_p>0\), the number of ordered pairs realizing the three records is
\[
 \ll_{p,H,\varepsilon} X^\varepsilon
          \big(h_p(d)h_p(d')\big)^{C_p}.
 \tag{1}
\]
The statement is uniform in those records. If both fields lie in a fixed finite list, the right side is \(O(X^\varepsilon)\). More generally the exponent loss is \(O_p(\zeta)\) in units of \(\log n\) when \(d d'\le n^\zeta\), uniformly for \(n=X^L\) and fixed \(L>0\). The exact bound (1) can be better, since it depends on torsion in the two class groups.

Write \(U_1/V_1=c t^p\) with \(t\in\mathbb Q_{>0}\). The two triples give a point on
\[
 \mathcal C:\quad
 w^2=d(\sigma c t^p+\tau)(\sigma rc t^p+\tau).
 \tag{2}
\]
The two factors have disjoint simple roots since \(c,r\ne0\) and \(r\ne1\). The smooth projective curve has genus \(p-1\). The coordinate \(t\) determines both primitive odds and therefore both triples. In particular there is no free common scale to count.

Over \(K\) the action \(\theta(t,w)=(\zeta_p t,w)\) gives a degree-\(p\) map to the smooth projective conic
\[
 D:\quad w^2=d(\sigma c s+\tau)(\sigma rc s+\tau),\qquad s=t^p .
\]
The conic has genus zero, including when it has no \(K\)-point. Over an algebraic closure, \(s\) has simple zeros at \(s=0,w=\pm\sqrt d\), and simple poles at the two points over infinity, separated by \(w/s=\pm c\sqrt{dr}\) up to sign. Thus the cover ramifies exactly at these four points, fully, with respectively the Kummer weights \(1,1,-1,-1\). Denote their unique inverse images by \(P_0^+,P_0^-,P_\infty^+,P_\infty^-\).

For the Jacobian \(J\) of \(\mathcal C\), the trace of the deck action factors as pullback after pushforward through the Jacobian of \(D\), so it vanishes. With \(\phi=1-\theta\) the cyclotomic identity gives \(\phi^{p-1}=[p]\) times an automorphism and \(\deg\phi=p^2\). As a \(G_K\)-module,
\[
 J[\phi]=
 \{(a_+,a_-,b_+,b_-)\in\mathbb F_p^4:
                       a_++a_-+b_++b_-=0\}
/\langle(1,1,-1,-1)\rangle
 \simeq \mathbb F_p(\chi_d)\oplus\mathbb F_p(\chi_{d'}).
 \tag{3}
\]
Here \(\chi_u\) is the quadratic character of \(K(\sqrt u)/K\), with the trivial character allowed. To verify (3), differences of branch points are fixed by the deck action and killed by \(p\), since \(p\) times such a difference is pulled back from a principal divisor on the geometric conic. In any principal divisor supported on the four points, the defining function is a deck eigenfunction because its divisor is invariant. It is \(t^jR\) for a function \(R\) of the conic, so its four orders modulo \(p\) are \(j(1,1,-1,-1)\). The vector is realized by \(\operatorname{div}(t)\), and multiples of \(p\) of branch differences are principal. This proves the displayed quotient, whose order \(p^2\) exhausts the kernel. For odd \(p\) the antisymmetric difference of the two zero points and the antisymmetric difference of the two pole points form a basis. Their Galois actions are the two stated characters. This describes the field loss without passing to the compositum of the two quadratic extensions.

For \(c,r\in\mathbb Q^\times\) let \(B\) be the positive integer obtained by multiplying \(2pd\) by the absolute reduced numerators and denominators of \(c,r,r-1\). Away from the places above \(B\), (2) gives a smooth proper model by gluing the two charts with parameters \(t\) and \(t^{-1}\), rescaling \(w\) by \(t^{-p}\) on the second chart. The leading coefficient is a unit and the two factors have disjoint simple roots on reduction. At a simple root the derivative in \(t\) is a unit, and off the roots the degree-two cover is étale. Hence \(J\) has good reduction there. The endomorphism \(\phi\) extends and is finite étale outside \(B\), since it divides \([p]\). The isogeny Kummer image of \(J(K)/\phi J(K)\) is unramified there by extending a \(K\)-point to the abelian scheme and pulling back the finite étale isogeny.

For each summand \(\mathbb F_p(\chi_u)\) in (3), restriction to \(K_u=K(\sqrt u)\) is injective on \(H^1\) by inflation–restriction, since \([K_u:K]\) divides two and \(p\) is odd. The target is \(K_u^\times/K_u^{\times p}\) because \(K_u\) contains \(\mu_p\). Restriction of the local unramified classes remains unramified. Thus the image in these two summands is bounded by the square of an \(S\)-unit factor of size \(\exp O_p(1+\omega(B))\), times \(h_p(d)h_p(d')\). Indeed, a Kummer class in either field can have nonzero valuations modulo \(p\) only above \(B\). Divide its outside ideal exponents by \(p\) to obtain a \(p\)-torsion class in the class group of the \(S\)-integers. The kernel consists of \(S\)-units modulo \(p\)-th powers, with rank bounded by \(|S|+O_p(1)\). The \(p\)-torsion of the quotient class group has cardinality at most the cardinality of the original \(p\)-torsion. This gives
\[
 \operatorname{rank}J(K)
 \ll_p 1+\omega(B)+\log h_p(d)+\log h_p(d').
\]
For this implication, on the free Mordell–Weil group the determinant identity for \(\phi^{p-1}\) gives index \(p^{\operatorname{rank}J(K)/(p-1)}\), and the free quotient index is no larger than the index on \(J(K)\).

If there are any candidate pairs then \(\mathcal C(\mathbb Q)\) has a point for the Abel–Jacobi embedding. The polarized degree of the image is determined by its genus. An Abel–Jacobi curve of genus at least two contains no positive-dimensional abelian translate. Uniform Mordell–Lang gives a point count exponential in \(1+\operatorname{rank}J(K)\). Since \(B\le X^{O_{p,H}(1)}\) for any occupied record, the \(1+\omega(B)\) term is uniformly \(o_X(\log X)\). This proves (1). The discriminant bound \(|\operatorname{disc} K_u|\ll_p u^{p-1}\) follows from the discriminants of the cyclotomic and quadratic fields. Minkowski's class-number bound in fixed degree gives \(\log h_p(u)\ll_p1+\log u\), proving the claimed bound for small squarefree representatives.

The torsion and point-count steps use [Poonen–Schaefer, Proposition 6.2](https://math.mit.edu/~poonen/papers/descent.pdf) and [Gao–Ge–Kühne, Theorem 1.1](https://arxiv.org/html/2105.15085v4). The proof above describes which character appears at each branch pair.

## The exact original parity records

For each original triple define
\[
 A_i=[T_iV_i]_2,\qquad B_i=[T_iU_i]_2 .
\]
The two character classes in (1) are exactly \(A_1 A_2\) and \(B_1 B_2\). The second identity follows from \([r]_2=[U_1U_2V_1V_2]_2\). As the three letters are pairwise coprime, the pair \((A_i,B_i)\) determines the three exact squarefree parity parts. At any prime the three active roles give respectively the vectors \((0,1)\), \((1,0)\), and \((1,1)\), from \(U_i,V_i,T_i\). Thus the signature does not lower the unconditional price of the full parity record. Its graph use is the price of this joint record *conditional on the chosen outside column*.

Another way to price the first character is to fix the original target quotient class \([T_1/T_2]_2\). For fixed \(r=u/v\) in lowest terms the quotient \(V_1/V_2\) has only divisor-many choices. To check this, from the reduced fraction \(U_1/V_1\) the numerator and denominator of \(r U_1/V_1\) are obtained by cancelling
\[
 g=\gcd(u,V_1)\gcd(v,U_1).
\]
Primewise this follows from \((u,v)=(U_1,V_1)=1\), including when the exponents in the two crossing terms are unequal. Then \(V_2=v V_1/g\). There are at most \(\tau(u)\tau(v)=X^{o_X(1)}\) possibilities for \(g\), so pricing the original target quotient class also prices \(d\), up to a subpower. It does not make \(d'\) fixed as \(r\) varies. Equality of original target parity alone does not give trivial characters. For example the primitive rows \((1,3,4)\) and \((2,7,9)\) under \(T=U+V\) give \(d=21\) and \(d'=2\).

For a support-based estimate, write
\[
 D_q(x)=\prod_{\ell:\ q\nmid v_\ell(x)}\ell
\]
for positive rational \(x\). Outside any prescribed prime pool of polynomial prime product, the price of the anchor class can be bounded by the product of the outside \(D_p\) for its two original summands. The two summands have disjoint prime supports. The price of the target quotient class can be bounded by its outside mixed \(D_2\), which is at most the product of the two actual target parity parts. For each support cap the number of squarefree products is bounded by the cap, and recording the bounded residue exponent at each known prime costs \(X^{o_X(1)}\). Fixing the parts on a polynomial-height pool also costs \(X^{o_X(1)}\). Removing a pool for the catalogue estimate does not remove its possible contribution to the fields in (1).

## Outside graph consequences and the complete-label comparison

**MO2.** Use a simple bipartite graph with one genuine primitive triple on each edge. Its complete triple determines its outside labels. Let \(R,C,E\) be its row, column, and edge counts. For fixed ordered distinct rows, assume the common columns have at most \(K_0\) relative odds records and, conditional on each such record, at most \(F\) first-odds \(p\)-classes. Distinct row labels ensure that their odds cannot be equal, since the primitive fraction determines the complete triple and both labels. Assume the relative odds records use the fixed sign convention of the point count.

Suppose each column uses at most \(S\) distinct signatures \((A_i,B_i)\). Split the column by its signature. There are at most \(CS\) new columns. A wedge in a split column has \(d=d'=1\), so each ordered pair of rows, relative odds record, and first-odds class comes from at most \(O_{p,H,\varepsilon}(X^\varepsilon)\) such columns by (1). The ordered tuple determines the original outside column and its signature. If \(E\ge2CS\), Cauchy–Schwarz gives at least \(E^2/(2CS)\) wedges. Hence
\[
 E\ll_{p,H,\varepsilon} X^\varepsilon
           \left(CS+R(CS K_0 F)^{1/2}\right).
 \tag{4}
\]
One can split on the two signatures modulo a subgroup \(G\) of the rational squareclass group instead. Then \(d,d'\in G\). Accounting for their possible values gives the same estimate with \(F\) multiplied by the maximum, for one fixed \(r\), of
\[
 \sum_{\substack{d\in G\\ \operatorname{sf}_2(dr)\in G}}
       (h_p(d)h_p(\operatorname{sf}_2(dr)))^{C_p}.
\]
Each \(d\) in the sum is interpreted as its positive squarefree representative. For a fixed finite group this factor is absorbed in the constant. If \(G\) varies, the displayed factor is required. Bounding only the cardinality of the signature quotient does not control it.

Alternatively leave columns intact. For fixed ordered row labels and relative odds, suppose there are at most \(F\) first-odds classes and at most \(G_0\) original target quotient classes, and let \(P\) bound their uniform field factor in (1). Then
\[
 E\ll_{p,H,\varepsilon} X^\varepsilon
              \left(C+R(C K_0 F G_0 P)^{1/2}\right).
 \tag{5}
\]
In place of \(FG_0P\) one can use the sum of the factors in (1) over the distinct actually occupied joint records \((c,d)\) conditional on \(r\) and the row pair. The conversion of target quotient classes to \(d\) costs a uniform subpower as above. The same wedge count proves (5).

For instance, if the row and column exponents both equal \(1/6\) in units \(\log n\), (4) is strict below \(1/3\) when \(s+k+f<1/6\) with \(S\le n^s,K_0\le n^k,F\le n^f\), by a fixed margin. The exponent at general part exponents \(a,b\) is \(\max(b+s,a+(b+s+k+f)/2)\). At equal part exponent \(a=b\ge2/9\), this expression gives no strict improvement even with zero catalogue charges. The original caps remain imposed, and the graph theorem does not assert that the required conditional catalogues exist in an unrestricted endpoint family. For fixed \(p,L\), the inner \(X^\varepsilon\) allowances can be chosen after a strict margin is specified.

In the affine complete-label graph, with \(U_1,U_2\) fixed and the complete value \(V\) common, the original target difference is a fixed nonzero integer \(D=\sigma(U_2-U_1)\). If \(T_1,T_2\) have equal \(q\)-free coefficients \(a\), then \(T_i=a y_i^q\), and \(a(y_2-y_1)\) divides \(D\). There are divisor-many possibilities for \(a,y_2-y_1\), and for each the difference of the two positive \(q\)-th powers is strictly monotone in the smaller positive root. Hence the pair count is \(X^{o_X(1)}\) without any odd-order condition at \(U,V\).

More generally, for a fixed exact squareclass \(d\) of \(T_1T_2\), the same affine graph has \(X^{o_X(1)}\) completions without any odd-order condition. Indeed \(T_1T_2=d z^2\) with \(z\) integral, and
\[
 (T_1+T_2)^2-4d z^2=D^2 .
\]
If \(d=1\) factor the left side over the integers. Otherwise take the algebraic integer \(T_1+T_2+2z\sqrt d\) in the real quadratic field. Its principal ideal divides \((D^2)\) and has only divisor-many possibilities uniformly in \(d\). Each ideal gives one orbit under units. Both real embeddings are positive, their product is \(D^2\), and the larger one is at most \(4X^H\). After choosing the logarithm of the first embedding, the solutions in an orbit in the specified height range lie in an interval of length \(O_H(\log X)\). Unit spacing is bounded below absolutely because any real quadratic unit \(u>1\) has integer trace and norm \(1\) or \(-1\), giving \(u\ge (1+\sqrt5)/2\). If both embeddings are positive the unit ratios have norm one, which is more restrictive. The ideals dividing \((D^2)\) are \(X^{o_X(1)}\) since over each rational prime there are at most two prime ideals and the exponents are bounded by \(4v_\ell(D)\). This proves the claim. The nonaffine outside graph allows \(V_1\ne V_2\) in the same column and has no such fixed target difference.

The same-prime variant in which \(m\) is a power of \(\ell\), the odds are fixed modulo \(m\)-th powers and both target classes are fixed modulo \(\ell\)-th powers, does give a cyclic-cover count after substituting \(t^m\). Existing statements are stronger for that treatment: in the parity case MTK1 uses the single product class of the two common coordinate forms and the two targets. For odd \(\ell\), the proof of PT1 applies to the two targets and one common coordinate form, and outside the fixed relative-odds support the odds class already prices that coordinate's \(\ell\)-free class. Formula (1) addresses the mixed case with odd \(p\) and target parity.

## Necessity of a further arithmetic condition in unrestricted height

Split characters alone do not give a subpower pair estimate. Fix \(r=4\) and the positive sign convention. Take reduced rational \(a\) in \([9/4,7/3]\) and let
\[
 s=\frac{5-2a}{a^2-4}>0,\qquad w=1+a s.
\]
Then \(w^2=(s+1)(4s+1)\). The two positive primitive fractions of odds \(s,4s\) give a pair with \(d=d'=1\). If \(a\) has reduced numerator and denominator bounded by \(c\sqrt X\) for a sufficiently small constant \(c>0\), both triples have height at most \(X\). There are \(\gg X\) such reduced \(a\) by the standard positive density of coprime pairs in a fixed planar cone. The map \(a\mapsto s\) is injective on these parameters: \(w>0\) is uniquely determined by \(s\), and then \(a=(w-1)/s\). These are genuine primitive triples. This construction does not impose the small radical caps and does not imply an obstruction under them.

## Projective form and an orientation choice

The point count applies to four pairwise nonproportional homogeneous rational linear forms \(L_0,L_1,L_2,L_3\) of polynomial coefficient height, evaluated at a positive coprime integer pair \((x,y)\) of polynomial height, with all four values positive. Fix a class \([L_0(x,y)/L_1(x,y)]_p\) of representative \(c>0\) and a squarefree class \(e=[L_2(x,y)L_3(x,y)]_2\). Express
\[
 L_i=a_iL_0+b_iL_1\qquad(i=2,3).
\]
The coefficients are nonzero, of polynomial height, and \(a_2b_3-a_3b_2\ne0\). Put \(L_0/L_1=c t^p\) and divide the parity equation by \(L_1^2\). The resulting curve is
\[
 w^2=e(a_2c t^p+b_2)(a_3c t^p+b_3).
\]
The same proof gives (1) with the two fields of signed squareclasses \(e b_2b_3\) and \(e a_2a_3\). A negative squareclass is represented by its signed squarefree integer in \(h_p\), with the same argument. The primitive pair is recovered from \(t\), since \(L_0/L_1\) is a nonconstant fractional linear map of \(x/y\). The character computation can be done directly from the projective configuration: the two classes are those of
\[
 e\,\det(L_0,L_2)\det(L_0,L_3),\qquad
 e\,\det(L_2,L_1)\det(L_3,L_1).
\]
The denominator in Cramer's rule occurs squared. In particular their ratio is the squareclass of a cross ratio of the four zeroes, independent of \(e\) and of scaling the forms. Both character fields are trivial over \(\mathbb Q\) precisely when the relevant cross ratio is a rational square and one of the two displayed determinant classes is a rational square.

For the pair of actual triples at fixed relative odds, one may use \(x=\gcd(U_1,U_2)\), \(y=\gcd(V_1,V_2)\) and the two actual target forms \(L_i(x,y)=T_i\). These are four distinct forms for \(r\ne1\). Fixing their coefficients costs a uniform subpower because the prime support of the common-factor quotients is contained in the fixed numerator and denominator of \(r\). Thus one may choose any two of these four forms for the odd-order quotient and use the other two in the parity product. In particular this allows an odd-order condition on the quotient of the two actual targets and a parity condition on the two original summand directions. Any one of the six choices can be used in (5) with its own conditional records and actual field factor.

There is a smaller signature split in the original orientation if all the prescribed relative odds in the wedges being counted are squares. Then \(d\) and \(d'\) have the same class, so splitting columns by \([T_iV_i]_2\) alone suffices for (4). More generally, if relative odds classes lie in a fixed finite squareclass group \(G\), split on \([T_iV_i]_2\) modulo \(G\); both fields then lie in the finite list from \(G\). The all-wedge count (4) requires the relative-odds condition for every wedge retained by this split.

## Two independent descents for distinct orders

**MO3.** Fix a prime \(p\) and an integer \(q\ge2\)
with \(p\nmid q\). Use the four nonproportional linear forms above and fix two positive rational classes
\[
             L_0/L_1 = c\,t^p,\qquad L_2/L_3=e\,w^q.
\]
Here \(c,e\) fix the indicated classes, \(t,w\) are positive rational, and their sign restriction is inessential to the point bound. Assume the four form coefficients and the two representatives have heights at most \(X^{O(1)}\), with a fixed height exponent. Write \(L_i=a_iL_0+b_iL_1\) for \(i=2,3\), and set \(A_i=a_i c,\ B_i=b_i\). Thus all four coefficients are nonzero and \(A_2 B_3-A_3 B_2\ne0\). The resulting curve is
\[
 e\,w^q=\frac{A_2 t^p+B_2}{A_3 t^p+B_3},
 \qquad
 t^p=\frac{B_2-e B_3 w^q}{e A_3w^q-A_2}.
 \tag{6}
\]
It is geometrically integral and has genus \((p-1)(q-1)\), at least two. Put \(K=\mathbb Q(\zeta_p)\) and form the following two finite étale \(K\)-algebras:
\[
 E_{p,0}=K[W]/(W^q-B_2/(eB_3)),\qquad
 E_{p,\infty}=K[W]/(W^q-A_2/(eA_3)).
\]
For an étale algebra \(E\) with field factors \(E_j\), write
\[
               \mathcal H_{\ell}(E)=\prod_j|\mathrm{Cl}(E_j)[\ell]|.
\]
Set \(H_p=\mathcal H_p(E_{p,0})\mathcal H_p(E_{p,\infty})\).
Then the number of rational points on the smooth complete model of (6) has bound
\[
             O_{p,q,\varepsilon}(X^\varepsilon H_p^{C_{p,q}}).
 \tag{7}
\]
The constant also depends on the uniform coefficient height exponent and \(C_{p,q}\) only depends on the orders. No splitting fields or composita are used.

Indeed the degree-\(p\) cover of the \(w\)-line in (6) ramifies with index \(p\) at the \(2q\) distinct simple zeros and poles displayed by the two binomials, and is unramified elsewhere, in particular at infinity. The zeros already have valuation one in the geometric function field of the \(w\)-line, so its Kummer class is nontrivial and the geometric cover is connected. Riemann–Hurwitz gives the displayed genus. The cyclotomic isogeny \(\phi=1-\theta\) for \(\theta(t,w)=(\zeta_p t,w)\) has degree \(p^{2q-2}\), since the trace vanishes and \(\phi^{p-1}\) is \([p]\) times an automorphism.

Let \(\Omega_0,\Omega_\infty\) denote the two sets of branch points on the \(w\)-line, and let \(M_i=\mathbb F_p^{\Omega_i}\) be their permutation modules, with augmentation \(\operatorname{aug}\) the sum of the coordinates. The branch divisor description gives
\[
 J[\phi]
 \cong
 \frac{\ker(\operatorname{aug}:M_0\oplus M_\infty\longrightarrow\mathbb F_p)}
 {\langle(\boldsymbol 1_0,-\boldsymbol 1_\infty)\rangle}.
 \tag{8}
\]
For a direct proof, branch differences are fixed by \(\theta\) and their \(p\)-multiples are principal by pullback from the line. A defining function for any principal branch divisor is a deck eigenfunction \(t^j R(w)\), so the branch orders mod \(p\) lie in the span of the vector of zeros \(1\) and poles \(-1\). The divisor of \(t\) realizes this vector. Therefore the degree-zero classes give precisely the indicated quotient, whose cardinality is \(p^{2q-2}\) and exhausts \(J[\phi]\). This is a special case of [Poonen–Schaefer, Proposition 6.2](https://math.mit.edu/~poonen/papers/descent.pdf).

As \(q\) is a unit in \(\mathbb F_p\), each \(M_i\) splits into the constants and its augmentation-zero submodule. In (8) the constants contribute only the one line being quotiented out. Thus \(J[\phi]\) is the direct sum of the two augmentation-zero submodules, a direct summand of \(M_0\oplus M_\infty\). Its \(H^1\) injects into the cohomology of the two permutation modules. Shapiro's lemma applied to the stabilizer of each branch root identifies the latter cohomology with the product of \(E_j^\times/E_j^{\times p}\) for the field factors of \(E_{p,0}\) and \(E_{p,\infty}\). More explicitly, the isomorphism for one root orbit restricts a cocycle to the stabilizer of the chosen root and reads its coordinate at that root. The coefficient module over \(K\) is identified with \(\mu_p\).

Take the finite set of rational primes dividing \(p q\) or the numerator or denominator of a nonzero coefficient \(e,A_i,B_i,A_2 B_3-A_3 B_2\). Away from these primes, the two sets of branch points are finite étale and remain disjoint. Normalize the cyclic cover of the \(w\)-line using the equation for \(t\) near zeros, that for \(t^{-1}\) near poles, and the étale cover elsewhere, including infinity. These charts give smooth proper reduction, since the branch orders are one and \(p\) is a unit. The isogeny is finite étale there. Its Kummer image is unramified outside the specified primes. The explicit Shapiro map just given preserves this local condition, by restricting to inertia in each factor field and then taking the root coordinate. Kummer classes in every factor field therefore have valuations divisible by \(p\) away from the designated primes. Apply the \(S\)-unit and \(p\)-class bound to each factor. The sum of the factor degrees is bounded in terms of \(p,q\), so the unit factor has logarithm \(O_{p,q}(1+\omega(B))\), where \(B\) is the product of the rational primes in the finite set. Since \(B\) has polynomial height this is \(o_X(\log X)\). The free-group determinant identity for \(\phi^{p-1}\) then bounds the rank in terms of this logarithm and \(\log H_p\). Uniform Mordell–Lang gives (7), with a rational base point for the Abel–Jacobi embedding chosen only if there is a rational point to count. The tame Kummer and étale-algebra framework also appears in [Bruin–Poonen–Stoll, Sections 3, 6, and 7](https://math.mit.edu/~poonen/papers/genus3.pdf).

When \(q\) is another prime, reverse the two projections. Use the \(q\)-cover of the \(t\)-line. It has \(p\) zeros and \(p\) poles, separate and simple, from the two binomials
\[
                 T^p+B_2/A_2,\qquad T^p+B_3/A_3 .
\]
Over \(K'=\mathbb Q(\zeta_q)\) form the separate factor algebras \(E_{q,0},E_{q,\infty}\) of these binomials and set \(H_q=\mathcal H_q(E_{q,0})\mathcal H_q(E_{q,\infty})\). Since \(q\nmid p\), the same proof gives the alternative to (7) with \(H_q\). The point count thus uses \(\min(H_p,H_q)^{C_{p,q}}\) after taking the larger of the absolute-in-the-record exponents and constants. The second factor \(H_q\) is independent of the class \(e\), while the first factor \(H_p\) is independent of the class \(c\) because \(A_2/A_3=a_2/a_3\). Each unwanted class still affects the bad prime support, already covered by \(X^\varepsilon\).

In particular, at \(q=2\) the \(p\)-descent gives precisely the two quadratic characters in (3) up to bounded factors from the fixed base field. The other descent uses only \(2\)-torsion of the separate \(p\)-root factor fields over \(\mathbb Q\). It is uniformly subpower if the two coefficient ratios \(-B_i/A_i\) range over a fixed finite list modulo \(p\)-th powers, regardless of the parity twist \(e\). In the original odds orientation these ratios are \(-\tau/(\sigma c)\) and \(-\tau/(\sigma r c)\). Thus it suffices that the relative odds and first odds range over fixed finite lists of \(p\)-classes. For a varying collection the bound uses the actual factor-field torsion or the finite-height discriminant bound. Representatives of branch classes supported on squarefree products at most \(n^\rho\) have heights at most \(n^{O_{p,q}(\rho)}\), and each resulting field discriminant is at most \(n^{O_{p,q}(\rho)}\) by clearing denominators in the binomial. This gives field loss \(n^{O_{p,q}(\rho)}\) by Minkowski. A fixed finite list can therefore be replaced by genuinely small branch-class height; cardinality alone does not control this cost.

The graph argument (5) uses the joint records of the two classes at exponents \(p,q\) and either of these field factors. There is no condition on numerical divisor depth. All original caps can still be retained. In the balanced outside box with part exponents \(1/6\), subpower relative odds count and uniform subpower field cost, the resulting sector is strict when the joint conditional class price is below \(1/6\) by a fixed margin. No estimate here supplies this record-price hypothesis from the three caps in general.

There is an explicit obstruction to omitting the target class from a point bound based only on two original caps. Fix odd prime \(p\), positive integer \(m\), and \(L=3/(mp)\). Choose coprime \(a,b\) with \(a\) even and \(b\) odd, in fixed intervals of size asymptotic to a positive constant times \(X^{1/(mp)}\), chosen so that \(a^{mp}\) lies between \(X/20\) and \(X/16\) and \(b^{mp}\) between \(3X/5\) and \(5X/8\). Set
\[
 (U_1,V_1,T_1)=(a^{mp},b^{mp},a^{mp}+b^{mp}),\qquad
 (U_2,V_2,T_2)=(4a^{mp},b^{mp},4a^{mp}+b^{mp}).
\]
Both original targets are in the shell \(X/2<T_i\le X\). The two original summand radical caps \(n^{1/3}\) hold in each row, since the extra factor \(4\) introduces no new prime. The relative odds are \(4\), and the first odds have trivial \(p\)-class. The fixed intervals and congruences have \(\gg_{m,p}X^{2/(mp)}=\gg_{m,p}n^{2/3}\) coprime pairs, by elementary Möbius inversion restricted to odd common divisors. The example imposes no target class restriction.
At sufficiently large fixed \(m\), each of its two row
projections is injective in the parameter pair. Deleting
the rows with \(\operatorname{rad}(T_i)\le n^{1/3}\)
loses only \(O(n^{1/3+\eta}X^\varepsilon)\) choices for
each projection by the effective original count, since
all three caps would then hold. For small fixed
\(\eta,\varepsilon/L\) the remaining set still has
\(\gg n^{2/3}\) pairs. Both complete targets violate
their caps. The preceding rational conic family keeps
the parity fields split while leaving the odd-order
record free.
[Effective original count (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
