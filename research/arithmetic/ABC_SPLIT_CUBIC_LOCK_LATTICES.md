# Split cubic forms and the index of common-lock lattices

Three full divisors on a primitive additive family give a rank-two lattice. Its three residual letters are primitive linear forms with explicitly prescribed pair determinants. A further common projective congruence has a cyclic index. The number of possible residue classes for this further index can be computed exactly, and their realization by actual primitive triples includes one forced parity correction.

## Lattice normalization

Let \(\Lambda\subseteq\mathbb Z^2\) be a rank-two lattice containing a vector \((a,b)\) with \(\gcd(a,b)=1\). Denote by \(F,G,H\) the positive generators of the images of \(\Lambda\) under \(a,b,a+b\), respectively. They are pairwise coprime by the assumed primitive vector. The lattice of just these divisibilities is
\[
        \Lambda_0=\{(a,b)\in\mathbb Z^2:F\mid a,\ G\mid b,\ H\mid a+b\}.
\]
Putting \(a=Fx,\ b=Gz\) gives \(z\equiv -F G^{-1}x\bmod H\). Hence \([\mathbb Z^2:\Lambda_0]=FGH\). Set \(J=[\Lambda_0:\Lambda]\).

Choose an oriented basis of \(\Lambda\), and let \(L_1,L_2,L_3\) be the cofactor forms \(a/F,b/G,(a+b)/H\) in its integer parameters. The definitions of the three contents make the forms primitive, and
\[
  FL_1+GL_2=HL_3,\qquad
  \det(L_1,L_2)=JH,\qquad
  \det(L_1,L_3)=JG,\qquad
  \det(L_2,L_3)=-JF.                                         \tag{SL1}
\]
Indeed the determinant of \((a,b)\) in the oriented basis is \(FGHJ\), and the signed relation replaces the corresponding form in the other two determinants. In particular the split binary cubic \(L_1L_2L_3\) has discriminant \((FGHJ^3)^2\), the product of the squares of the three pair determinants. This uses the ordinary binary cubic discriminant, including any root at infinity. For \(J=1\) the pair determinants are precisely the three common divisors up to the displayed signs.

An oriented unimodular change sets \(L_1=x\), since it is primitive. After an integral shear in \(y\), the unique normal form is
\[
 L_1=x,\qquad L_2=ux+JH y,\qquad L_3=vx+JG y,\qquad
                 Gu-Hv=-F,\quad 0\le u<JH.                 \tag{SL2}
\]
The shear changes \(u\) by a multiple of \(JH\). Two such oriented bases with first cofactor \(x\) differ by just such a shear, so the range for \(u\) is canonical. Conversely these forms define the lattice through \((a,b)=(Fx,G(ux+JH y))\). If \(L_2,L_3\) are primitive it has the prescribed contents and index. It contains a primitive vector since no prime can divide all entries of its basis matrix. To see this directly, such a prime would divide \(F\) and both \(Gu,GJH\). It is coprime to \(G,H\), contradicting \(\gcd(u,JH)=1\). Equivalently the first invariant in Smith normal form is one.

## Exact count of oblique congruences and unit constraints

Let \(u_0\) be the unique class modulo \(H\) with \(Gu_0\equiv-F\bmod H\), taken in \([0,H)\), and let \(v_0=(F+Gu_0)/H\). Then the candidates in (SL2) are
\[
                       u=u_0+kH,\quad v=v_0+kG,\quad
                                  0\le k<J.
\]
Primitivity of the forms is exactly \(\gcd(u,JH)=\gcd(v,JG)=1\). Primes dividing \(H\) already make \(u\) a unit, and primes dividing \(G\) make \(v\) a unit. At \(p\mid J\) and \(p\nmid FGH\), the vanishing classes of \(u\) and \(v\) in \(k\bmod p\) are distinct, since \(Gu-Hv=-F\) is a unit. At \(p\mid J\) and \(p\mid FGH\) there is just one forbidden class. When \(p\mid F\) the two classes coincide, and when \(p\mid G\) or \(p\mid H\) one of the two forms is already a unit. Thus the exact number of lattices with these contents and index is
\[
          J\prod_{\substack{p\mid J\\p\nmid FGH}}\left(1-\frac2p\right)
           \prod_{\substack{p\mid J\\p\mid FGH}}\left(1-\frac1p\right).
                                                               \tag{SL3}
\]
These are distinct sublattices of the original \((a,b)\)-plane. At \(p\mid J\) the valid congruence records information beyond the three common divisors. Its degree of refinement is expressed by \(b/G\equiv u\,a/F\bmod JH\).

Local primitive representation in these coordinates is exact. Outside \(FGHJ\), the maps \((x,y)\mapsto(a,b)\) are invertible modulo \(p\), and the three cofactor lines have distinct reductions. Thus a primitive parameter is required, and at most one cofactor is divisible by \(p\). At \(p\mid J\), the three forms are nonzero multiples of \(x\) modulo \(p\). A primitive original row must have \(x\) a unit, and then all three cofactors are units. At a prime of \(FGH\) outside \(J\), permute the roles for the local test so the prime belongs to \(F\). Primitivity of the original row is equivalent to \(L_2\) being a unit. This also makes \(L_3\) a unit by (SL1). The full divisibility condition on \(F\) at this prime requires \(L_1\) to be a unit as well. The excluded lines are distinct since \(\det(L_1,L_2)=JH\) is a unit. For the signed permutations of the additive equation the unit statements are unchanged.

Consequently on the full-lock primitive locus the density of allowed residues modulo \(p\), among all \(p^2\) parameter pairs, equals \(1-p^{-2}\) off \(FGHJ\), \(1-1/p\) on \(J\), and \((1-1/p)^2\) on \(FGH\) outside \(J\). The three original residual radical caps correspond exactly to separate radical caps on \(L_1,L_2,L_3\), with the fixed lock primes stripped. More generally all valuations at a prescribed polynomial-height support can be recorded before imposing the caps outside it. The original real height region is carried to the same three linear inequalities in (SL2). This yields a split-cubic formulation without dropping the separate caps.

## The exact span of the full-lock primitive locus

For a row in \(\Lambda\), call its locks exact if \(L_1\) is a unit at every prime of \(F\), \(L_2\) is a unit at every prime of \(G\), and \(L_3\) is a unit at every prime of \(H\). Define
\[
                    E=\begin{cases}
                       2,&2\mid FGH,\ 2\nmid J,\\
                       1,&\text{otherwise}.
                       \end{cases}
\]
The lattice generated by all positive primitive rows with these exact locks has index \(E\) in \(\Lambda\). Moreover two such rows with all three complete projections distinct generate a lattice of that index. Their largest letters can lie in any prescribed relative range \(c_2/c_1\in(1,1+\xi)\), for fixed \(\xi>0\). Thus
\[
 \gcd_{\mathbf a,\mathbf b}
          |a_{\mathbf a} b_{\mathbf b}-a_{\mathbf b} b_{\mathbf a}|
                          = FGHJ E ,                         \tag{SL4}
\]
where the gcd ranges over that primitive exact-lock locus.

The lower bound follows from membership in \(\Lambda\), with the extra mod-two class when \(E=2\). In that case the two distinct excluded parameter lines modulo two leave one projective class, so every cross determinant of parameters is even.

For attainment, put \(Q=\operatorname{rad}(FGHJ)\) and take coprime positive integers \(x=jQ+1,\ x'=(j+1)Q+1\), with \(j\ge1\) large enough that \(1<x'/x<1+\xi\). Choose integer \(d,d'\) with \(xd'-x'd=E\). Use the two parameter pairs
\[
                         (x,xy+d),\qquad (x',x'y+d').
\]
Their determinant is \(E\), and their individual coordinate gcds divide \(E\). For \(E=2\) both first entries are odd, so the pairs are primitive. At primes on \(J\) the unit conditions are automatic since \(x,x'\equiv1\bmod Q\). At each \(p\mid FGH\) outside \(J\), these first entries are both one modulo \(p\). There is just one excluded residue of \(y\) for each row. This follows in (SL2) by testing \(L_2\) when \(p\mid F\) or \(p\mid G\), and \(L_3\) when \(p\mid H\). The tested slope is a unit in each case. For odd \(p\) at most two excluded classes leave a choice. The only case with \(p=2\) here has \(E=2\), so \(d'\equiv d\bmod2\) and the two excluded classes coincide. Choose \(y\) by the Chinese remainder theorem and make it sufficiently large and positive. All constraints hold and the largest-letter ratio tends to \(x'/x\) from values greater than one. The three complete projections are distinct for large \(y\), as \(x'\ne x\). These rows have cross determinant of absolute value \(FGHJ E\). They prove (SL4) and generate the asserted sublattice. Their three common contents are precisely \(F,G,H\): a further common prime outside \(FGH\) would have to divide their determinant, but the cofactors are units on \(J\).

In particular, if \(\Lambda\) is already generated by a family of positive primitive rows whose contents \(F,G,H\) are exact on every row, then \(2\mid J\) if and only if \(2\mid FGH\). One direction follows from (SL3) and the other from the forced parity in (SL4). Conversely every lattice counted in (SL3) with this parity condition is generated by two such rows. Thus the parity constraint on the cyclic index is exact, with no additional analogous constraint at an odd prime.

If \(N\) distinct positive primitive exact-lock rows of this lattice have \(c\le X\) and all \(a/c\) lie in an interval of width \(\Delta\), then
\[
                         N\le 1+\frac{\Delta X^2}{FGHJ E}.  \tag{SL5}
\]
Order the distinct ratios \(a/c\). The difference of two consecutive ratios is their projective determinant divided by the product of their targets. Its numerator is the same cross determinant as in (SL4) up to sign, hence at least \(FGHJ E\). Summing the consecutive spacings proves (SL5).

## The split order and its index form

There is an arithmetic interpretation of the extra index \(J\). Let \(\Lambda^\flat\) and \(\Lambda_0^\flat\) be the images of \(\Lambda\) and \(\Lambda_0\) under \((a,b)\mapsto(a,-b)\). The latter is closed under coordinatewise multiplication. Its coordinates are multiples of \(F,G\), and they agree modulo \(H\). Since \(J\Lambda_0^\flat\subseteq\Lambda^\flat\), the lattice \(J\Lambda^\flat\) is also closed under coordinatewise multiplication. In fact for any positive integer \(d\),
\[
       d\Lambda^\flat\text{ is closed under multiplication}
                            \quad\Longleftrightarrow\quad J\mid d .
                                                               \tag{SL6}
\]
For necessity, apply closure to the square of \(d(F,-Gu)\). This requires \((dF^2,-dG^2u^2)\in\Lambda\). The original congruence in (SL2) then requires \(JH\mid d u(Gu+F)=dHuv\). Both \(u\) and \(v\) are units modulo \(J\), so \(J\mid d\). For sufficiency, the product of two vectors of \(\Lambda^\flat\) lies in \(\Lambda_0^\flat\), hence \(d\) times that product lies in \(\Lambda^\flat\) when \(J\mid d\).

The resulting order in the split cubic algebra \(\mathbb Q^3\) is
\[
            \mathcal R=\{(t+Ja,t-Jb,t):(a,b)\in\Lambda,\ t\in\mathbb Z\}
                              \ \subseteq\ \mathbb Z^3 .
\]
The diagonal copy of \(\mathbb Z\) and multiplication closure in (SL6) show directly that this is a subring. Its index in \(\mathbb Z^3\) is \(FGHJ^3\), from the lattice of its two coordinate differences. For \(\gamma=(t+Ja,t-Jb,t)\) with all three coordinates distinct, its monogenic suborder is spanned by \(1,\gamma,\gamma^2\), since the characteristic polynomial of \(\gamma\) is monic of degree three and is minimal. The absolute Vandermonde determinant of these vectors is \(J^3|ab(a+b)|\). Therefore
\[
                     [\mathcal R:\mathbb Z[\gamma]]
                                 =|L_1 L_2 L_3|.
\]
This identifies the split cubic in (SL1) with an explicit index form. Translation of \(\gamma\) by a diagonal integer does not affect the correspondence. The three original radical caps restrict the three factors of its index separately, and the numerical box imposes three separate bounds on differences of the coordinates of \(\gamma\). The proof uses the explicit split order rather than a count of arbitrary cubic orders of the same discriminant.

Its conductor in \(\mathbb Z^3\) has coordinate ideals \(FHJ^2\mathbb Z,\ GHJ^2\mathbb Z,\ FGJ^2\mathbb Z\). A vector on the first axis has \(b=t=0\); from \(u(a/F)\equiv0\bmod JH\) and \(\gcd(u,JH)=1\), the least nonzero absolute first coordinate is \(FHJ^2\). The second axis has \(a=t=0\), so the least absolute second coordinate is \(GHJ^2\). On the third axis \(a=-b=FGw\) for an integer \(w\), up to sign. The congruence for the lattice is equivalent to \(J\mid wv\) by \(Gu+F=Hv\), so \(J\mid w\). The least absolute third coordinate is \(FGJ^2\). Since any ideal in \(\mathbb Z^3\) is the product of its coordinate ideals, these three minima give the conductor. Their common gcd is \(J^2\).

## Many actual residue classes at fixed radical cost

Formula (SL3) gives an explicit obstruction to pricing the extra congruence by its radical. Fix distinct odd primes \(p,q,r,s\) and a real \(\theta\ge0\). Set
\[
 J=p^k,\qquad
 G=r^{\lfloor\theta\log_r J\rfloor},\qquad
 H=s^{\lfloor\theta\log_s J\rfloor},\qquad
 F=q^{\lfloor\log_q(GHJ)\rfloor},
 \qquad X=GHJ .
\]
Then \(X\asymp J^{1+2\theta}\) and \(X/q<F\le X\). The lock primes and the radical of \(J\) are taken from the same fixed finite set. There are \(J(1-2/p)\) different lattices, and every one is generated by two positive primitive exact-lock triples with all complete letters \(\asymp X\), with constants uniform in the lattice and \(k\). Indeed one may use the preceding construction with the fixed common multiple \(pqrs\) in place of \(Q\), fixed \(j,d,d'\), and choose \(y\) from a fixed block of length \(pqrs\) sufficiently far to the right. Then \(u<JH\) and \(v=(F+Gu)/H\) give those height bounds. By making the starting \(y\) sufficiently large, the largest letters can all be put in one common dyadic shell.

A positive proportion of these classes can be realized by such pairs for which both unlocked cofactors \(b/G,c/H\) on both rows have radical greater than \(X^{1/[2(1+2\theta)]}\). To see this, pigeonhole a single allowed \(y\) for a fixed positive proportion of the classes. Along the remaining lifts \(u=u_0+jH\), each of the four cofactors is a nonconstant affine function of \(j\), so each is injective. All are \(O(X)\). The number of positive integers of size \(O(X)\) and radical at most \(X^{1/[2(1+2\theta)]}\) is \(O_\varepsilon(X^{1/[2(1+2\theta)]+\varepsilon})=o(J)\), by enumerating the radicals and their supported powers. Discard these cases for the four functions. [Uniform supported-power count](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

In this construction the sum of the three lock weights \(\log_X(FGH)\) tends to \(2-1/(1+2\theta)\). The determinant refinement fills the remaining \(1/(1+2\theta)\) of numerical weight and has fixed radical, yet it still has \(X^{1/(1+2\theta)+o(1)}\) available classes. The construction preserves the full cap on the original first letter, since the additional factor there is fixed. It does not preserve the two additional small radical caps. The classification therefore gives the exact obstruction to treating the residual determinant modulus as determined by its radical and complete divisors. The joint occupancy bound for the three residual split-cubic factors remains an arithmetic requirement.

The unimodular normalization and the cubic discriminant identity are classical lattice and form calculations. The present use keeps the three ordered cofactor contents, the further cyclic index, and the primitive full-lock unit conditions as separate data. [Bhargava–Shnidman, binary cubic discriminants and orders](https://arxiv.org/pdf/1206.4746#page=7) · [Common three-role divisor records](ABC_THREE_ROLE_DIVISOR_RECORD_OBSTRUCTION.md)
