# Refined spectral splittings in the degree-forty-two residual

**Theorem PolynomialWeylFortyTwoSecondSpectralSplitting.**
For the actual degree-forty-two residual over \(\mathbb C\), take
the centered \(P\) and the mate \(Q^\ast\) from the
[normalized cusp relation](POLYNOMIAL_WEYL_FORTY_TWO_NORMALIZED_CUSP_RELATION.md),
so that \(S=(Q^\ast)^2-F(P)\) has lowest Euler grade \(-51\).
For every spectral value \(\lambda\), each size of Jordan block
at each eigenvalue of the finite monodromy of \(P-\lambda\)
occurs with multiplicity divisible by six.
For the finite monodromy of \(Q-\lambda\), each size at
each eigenvalue occurs with multiplicity divisible by
twenty-one.

For a noninteger \(f\)-root class modulo \(\mathbb Z\),
count its six-root packets from the cusp relation in each
occupied \(3\mathbb Z\)-suborbit. In the integer class
include one base packet in each of the three integer
suborbits. The decreasing list of positive counts is
\(r=(r_1,\ldots,r_d)\).
The normalized actual partitions for \(P\) and \(Q\),
with multiplicities divided by six and twenty-one
respectively, both dominate \(r\).
If the class is noninteger and occupies one suborbit,
both partitions are exactly \(r\). Both integer partitions
after this division have exactly three parts, and equal
\((1+e,1,1)\) if all \(e\) excess integer packets
occupy one suborbit.

## The scaled spectral matrix to the next nonscalar order

Use the flat spectral construction of the
[regular spectral theorem](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_SPECTRAL_MONODROMY.md)
for \(P-\lambda\), with \(Q^\ast\) as the mate.
Set \(\lambda=s^{18}\), \(\tau=s^{-1}\), \(X=z/s\).
The operators
\[
 \mathcal P_s=s^{-18}P(z/s,sD_z),\qquad
 \mathcal Q_s=s^{-63}Q^\ast(z/s,sD_z),\qquad
 \mathcal S_s=s^{-51}S(z/s,sD_z)
\]
have coefficients polynomial in \(\tau\), with respective
limits \(P_0,Q_0,S_0\) in the \(z\) coordinate.
The scalar relation and cusp relation give
\[
 [\mathcal Q_s,\mathcal P_s]=\tau^{81},\qquad
 \mathcal Q_s^2=s^{-126}F(s^{18}\mathcal P_s)
                                 +\tau^{75}\mathcal S_s .       \tag{1}
\]
Take jets of solutions of \((\mathcal P_s-1)y=0\)
at \(z=1\). Let \(\mathsf Q(s)\) denote the jet
matrix of \(\mathcal Q_s y\), obtained by reducing
derivatives in the scalar equation.
The leading differential coefficient is the
\(\tau\)-independent \(\alpha z^{12}\), so all
reductions can be done over \(\mathbb C[[\tau]]\)
for germs at \(z=1\). Modulo \(\tau^{81}\),
\(\mathcal Q_s\) commutes with \(\mathcal P_s\),
so it acts as an endomorphism of the solution
space over that truncated ring.
Consequently the matrix of \(\mathcal Q_s^2\)
on solution jets is \(\mathsf Q(s)^2\bmod\tau^{81}\).
In (1) the term \(s^{-126}F(s^{18}\mathcal P_s)\)
acts exactly as the scalar \(a(s)=s^{-126}F(s^{18})\).
The coefficient of \(\tau^{75}\mathcal S_s\) at
its leading order gives the endomorphism
\(\mathsf S_0\) induced by \(S_0\) on the solution
space of \(P_0-1\).

Let \(\mathcal B(s)\) denote the horizontal spectral
matrix for the scaled jets at \(z=1\).
The chain rule at fixed \(z\) gives
\(\partial_s y=-18s^{80}\mathcal Q_s y-(z/s)D_zy\).
Thus \(K(s)=-\mathcal B(s)/(18s^{80})\) satisfies
\[
                             K(s)^2
               = a(s)I+\tau^{75}\mathsf S_0+O(\tau^{76}).        \tag{2}
\]
Indeed \(K(s)=\mathsf Q(s)+O(\tau^{81})\).
The leading matrix at \(\tau=0\) is the endomorphism
\(\mathsf Q_0\) induced by \(Q_0\).
Its square is \(\kappa I\) from the earlier
spectral theorem, and by the cusp relation
\[
                       \mathsf S_0^3=c\,\mathsf Q_0,
                       \qquad c\ne0.                            \tag{3}
\]
In particular \(\mathsf S_0\) commutes with
\(\mathsf Q_0\).

## The second formal splitting and its deck orbit

Use the formal block diagonalization in the earlier
spectral theorem for the two leading signs
\(u_\pm(\infty)=\pm\sqrt{\kappa}\).
Write \(u_\pm(s)=\pm\sqrt{a(s)}\), as formal
scalar units in \(\mathbb C[[\tau^{18}]]\).
A gauge \(T=I+O(\tau)\) changes \(K\) by conjugation
and a derivative term of order at least \(\tau^{82}\),
since \(s^{-80} T^{-1}T'=O(\tau^{82})\).
Hence (2) still holds in the block frame, now with
\(\mathsf S_0\) restricted to each sign.
The sum \(K_\pm+u_\pm(s)I\) is invertible on each
block. Factoring \(K_\pm^2-u_\pm(s)^2I\) in (2)
then gives the spectral block matrices in the form
\[
 \mathcal B_\pm(s)
   =-18s^{80}u_\pm(s)I
       -\frac{9}{u_\pm(\infty)}s^5\mathsf S_{0,\pm}
       +O(s^4).                                                \tag{4}
\]
By (3), the leading nonscalar matrix in each sign is
diagonalizable, with possible eigenvalues the three
different cube roots of \(-729c/\kappa\).
Indeed the cube of \((-9/u_\pm(\infty))\mathsf S_{0,\pm}\)
is \((-729/u_\pm(\infty)^3)c\,u_\pm(\infty)I\).

Within each sign the scalar term in (4) commutes
with every gauge. Thus the same off-diagonal
commutator recursion as before, now with leading
nonscalar order \(s^5\), gives the formal splitting
by these secondary scalars. The spectral monodromy
matrix is horizontal, so preserves every joint
block. For two distinct signs the off-diagonal
horizontal equation is ruled out by its leading
order \(s^{80}\), and for equal signs and distinct
secondary scalars it is ruled out by its leading
order \(s^5\), exactly as in the first splitting.

Let \(\zeta\) be a primitive eighteenth root of
unity. The analytic deck gauge obtained by the
scaled basepoint transport \(1\to\zeta^{-1}\)
still intertwines monodromy and the pullback
connection, as in the earlier theorem. In (4)
\(u_\pm(\zeta s)=u_\pm(s)\), so the pullback
multiplier \(\zeta^{81}=-1\) exchanges the signs
of the full scalar term.
The secondary multiplier is \(\zeta^6\), a
primitive cube root of unity. After the formal
gauges, the same off-diagonal valuation argument
forces the deck identification to take each joint
block to the corresponding permuted block.
The action on the six possible labels is
transitive, since simultaneous cycles of lengths
two and three have period six.
Thus all six blocks occur with identical
monodromy Jordan partitions over the formal
Laurent field. Analytic specialization at large
\(s\) and then flat spectral transport give the
stated Jordan multiplicities for every finite
\(\lambda\), as in the first splitting.


## The mate's twenty-one spectral blocks

For \(Q^\ast-\lambda\), take \(\lambda=s^{63}\) and
the same \(X=z/s\) and rescaled operators.
Spectral transport uses \(+P\), so the normalized
connection matrix at \(z=1\) is
\(K_Q(s)=\mathcal B_Q(s)/(63s^{80})\).
The chain term contributes \(O(\tau^{81})\) to
\(K_Q\), and \(\mathcal P_s\) acts on the
solution jets modulo \(\tau^{81}\), as before.
Since \(\mathcal Q_s=1\) on solutions, (1) gives
\[
        \Phi_s(K_Q)=I-\tau^{75}\mathsf S_0+O(\tau^{76}),
        \qquad \Phi_s(T)=s^{-126}F(s^{18}T),                 \tag{5}
\]
where \(\mathsf S_0\) now acts on solutions of \(Q_0-1\).

The seven possible leading eigenvalues \(v_i(\infty)\)
are the roots of \(\kappa v^7=1\).
For each there is a unique formal scalar solution
\(v_i(s)\) of \(\Phi_s(v_i(s))=1\) with that
leading value, in \(\mathbb C[[\tau^{18}]]\).
Indeed the constant-coefficient derivative with
respect to \(v\) is \(7\kappa v_i(\infty)^6\ne0\),
so successive coefficients are uniquely determined.
Block diagonalize by these seven leading values.
The gauge derivative only contributes at
\(O(\tau^{82})\) in (5), and factoring the polynomial
difference at \(v_i(s)I\) yields the connection
in block \(i\) as
\[
    \mathcal B_{Q,i}(s)
       =63s^{80}v_i(s)I
        -\frac{9}{\kappa v_i(\infty)^6}s^5\mathsf S_{0,i}
        + O(s^4).                                         \tag{6}
\]
On that block \(\mathsf S_{0,i}^3=c\,v_i(\infty)^5I\),
by (2) of the cusp relation with \(Q_0=1\).
Thus the secondary matrix in (6) is diagonalizable
with possible eigenvalues the cube roots of
\(-729c\,v_i(\infty)/\kappa\), using
\(\kappa v_i(\infty)^7=1\).

For a primitive sixty-third root \(\zeta\), the same
analytic deck gauge identifies the scaled companion
equations at \(s,\zeta s\). Write \(\omega=\zeta^{18}\).
Because \(\Phi_{\zeta s}(T)=\Phi_s(\omega T)\),
the scalar series in (6) is permuted exactly by
\(v_i(s)\mapsto\omega v_i(\zeta s)\), with
\(\zeta^{81}=\omega\).
The secondary eigenvalues are multiplied by \(\zeta^6\).
Cubing gives \(\omega\), in agreement with their
defining scalars. The primary labels have a cycle of
length seven. After seven deck steps, the secondary
multiplier \(\zeta^{42}\) has order three, so all
twenty-one joint labels form one orbit.
The same formal block and horizontal intertwiner
argument as for \(P\) gives matching Jordan
partitions on these twenty-one blocks and hence
divisibility by twenty-one for \(Q^\ast-\lambda\).

The mate change preserves the finite monodromy
of the scalar spectral equation. To see this, let
\(L_b=Q-bD(P)/2\) for \(b\in\mathbb C\), and choose a
polynomial \(G\) with \(G'=D/2\). Then
\([L_b,G(P)]=D(P)/2\). Along solutions of
\((L_b-\lambda)y=0\) one has
\((L_b-\lambda)(\partial_b-G(P))y=0\).
Since \(\deg D\le3\), \(D(P)\) has differential
order at most \(90<105\). Thus \(L_b\) has the same
leading differential coefficient for every \(b\),
and right reduction gives a flat analytic jet
transport in \(b\) on the punctured \(X\)-line.
By the same argument as for spectral flatness, this
conjugates finite monodromy from \(L_0=Q\) to
\(L_1=Q^\ast\).

## Discrete monodromy consequences

Each six-root packet of \(f\) has one root in each
residue modulo eighteen in its \(3\mathbb Z\)
progression. Thus the \(P_0-1\) model has each
block size \(r_i\) repeated six times, by the
homogeneous-model Jordan calculation of the regular
spectral theorem. Its three \(9\mathbb Z\)-suborbits
give \(Q_0-1\) each block size \(r_i\) repeated
twenty-one times, since each suborbit contributes
seven copies. Rank specialization of the
powers of \(M-\chi I\) at each eigenvalue \(\chi\)
therefore gives dominance of \(r\) by both actual
partitions with multiplicities divided by six and
twenty-one respectively.
If \(r\) has just one part, this forces equality.
For the integer eigenvalue the exact fixed dimensions
\(18,63\) give exactly three parts in the corresponding
normalized partitions. When \(r=(1+e,1,1)\), dominance forces
the largest part to be at least \(1+e\), and the
total size and the two other positive parts force
equality.

For \(P\), a noninteger class of size six has
scalar monodromy. A noninteger class of size twelve
has either twelve singleton blocks or six size-two
blocks, the latter forced if both packets occupy the
same \(3\mathbb Z\)-progression.
The corresponding \(Q\) classes are scalar in size
twenty-one, and have either forty-two singleton
blocks or twenty-one size-two blocks in size
forty-two, with the same forced case.
At the eigenvalue one, zero, six or twelve excess
integer roots of \(f\) give the normalized partitions
\((1,1,1)\), \((2,1,1)\), or one of
\((3,1,1),(2,2,1)\), respectively.
The first of the last two is forced if both
integer packets occupy the same \(3\mathbb Z\)-progression.
