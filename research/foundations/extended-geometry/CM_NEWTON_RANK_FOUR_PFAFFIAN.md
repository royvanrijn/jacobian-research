# A Pfaffian test for \(\mu\)-ordinary reduction in relative rank four

**Theorem.** Let \(A/K\) be a polarized abelian variety over a number
field with a Rosati-compatible geometric CM coefficient field \(E\)
preserved by \(G_K\). Put \([E:\mathbb Q]=2e\) and \(\dim A=4e\).
Suppose exactly one conjugate Hodge pair has signature \((2,2)\)
and all the others have signature \((0,4)\) or \((4,0)\).
For some auxiliary prime \(\ell\), assume the connected algebraic
monodromy has projective image \(\mathrm{PGL}_4\) on each
embedding block of \(E\). Then \(A/K\) has \(\mu\)-ordinary
reduction at a set of places of natural density one. Here
\(\mu\)-ordinary means that the relative Newton polygon attains
the average of the sorted Hodge vectors on each coefficient
orbit.

Let \(L/K\) be the field defining the \(E\)-action. When \(e>1\),
the ordinary density over \(K\) is \(1/[L:K]\). When \(e=1\)
the ordinary density is one. In the latter case the stated
monodromy hypothesis holds if \(E\) is the full geometric
rational endomorphism algebra, by the argument below.

Fité's Theorem 4 proves positive-density ordinary reduction for
fourfolds with a geometric imaginary quadratic action of
signature \((2,2)\). It permits extra geometric endomorphisms.
For full geometric imaginary quadratic endomorphism algebra,
the present result gives density one over the original
field. The higher-dimensional statement allows any CM
coefficient field under the given monodromy hypothesis. The
simple-signature results of Cantoral Farfán–Li–Mantovan–Pries–Tang
require a mixed pair with smaller multiplicity one.
[Fité, Theorem 4 and §8](https://arxiv.org/pdf/2005.10185v3),
[Cantoral Farfán et al., Assumption 2.5 and Theorems 4.18–4.19](https://arxiv.org/html/2508.11174v1).

The strict coefficient-norm inequalities already give the
desired midpoint comparison when the cyclic field action of
Frobenius excludes complex conjugation. On the remaining
components, a half power of Frobenius exchanges each block with
its symplectic dual. The determinant of one plus its normalized
square is a Pfaffian square divided by a determinant of even
valuation. This parity detects midpoint failures which the
coefficient-norm inequalities by themselves can allow.

## Midpoint bounds on the two types of field action

Work at a good absolute degree-one place \(v\) of residue
characteristic \(p\), excluding ramification in the coefficient
and endomorphism fields, primes dividing the order index or
polarization degree, and \(p=\ell\).
Write \(\gamma\) for the arithmetic Frobenius action on \(E\),
\(t=\operatorname{ord}(\gamma)\), and \(c\) for complex conjugation
on \(E\). Write \(\theta\) for coefficient Frobenius on the
embeddings, acting by postcomposition. The maps \(\theta\) and
precomposition by \(\gamma\) commute. The special-fiber
Hodge multiplicities \(m_\sigma\in\{0,2,4\}\) specialize from
characteristic zero and are invariant under
\(\sigma\mapsto\theta\circ\sigma\circ\gamma^{-1}\). In fact,
on the Hodge space the residue field automorphism raises
scalar eigenvalues to the \(p\)-th power and acts on the
endomorphisms by \(\gamma\). Applying it to the coefficient
eigenvector equation gives the stated permutation.
[Averaged Hodge polygons](CM_NEWTON_DISCONNECTED_SPECTRUM_INEQUALITIES.md#averaged-hodge-polygons-and-compression).

Let \(\lambda\) be a coefficient place, \(f_\lambda\) the length
of its \(\theta\)-orbit, and \(S_{\lambda,j}\) and \(B_{\lambda,j}\)
the prefixes of length \(j\) of its actual and reference slopes.
Polarization identifies the conjugate slopes with the reversed
complements. Actual polygons lie above the reference and share
their endpoints. Denote \(\Delta_{\lambda,j}=S_{\lambda,j}-B_{\lambda,j}\).
The exact-power argument in the
[simple-signature note](CM_NEWTON_SIMPLE_SIGNATURE_NONABELIAN.md#the-frobenius-power-equals-the-field-action-order)
uses only the individual \(\mathrm{PGL}_4\) images once these
have been established. It therefore applies here and gives,
outside a density-zero set, the midpoint inequality
\[
    \sum_{\lambda\mid p} f_\lambda\,\lceil t\Delta_{\lambda,2}\rceil
           \le 2t-1.                                             \tag{1}
\]
Here \(tB_{\lambda,2}\) is integral because \(\theta^t\)
preserves the original multiplicities, and the allowance
is \(2t-1\) because the only mixed pair has smaller multiplicity
two.

The permutation \(\theta\gamma^{-1}\) either fixes both mixed
labels \(\tau,\bar\tau\) or exchanges them. Let \(\epsilon=0\)
or \(1\) in the two cases and put \(d=c^\epsilon\gamma\).
Then \(\theta^i\circ\tau=\tau\circ d^i\) for every \(i\).
The mixed orbit therefore has length \(\operatorname{ord}(d)\)
and is self-conjugate exactly when \(c\in\langle d\rangle\).
If \(c\notin\langle\gamma\rangle\), either the two mixed labels
are in distinct conjugate orbits of length \(t\), or there is
one orbit of length \(2t\) containing both. In detail,
\(\epsilon=0\) always gives the first type. For \(\epsilon=1\),
the second type occurs when \(t\) is odd, since then
\((c\gamma)^t=c\), and the first type occurs when \(t\) is even.
Every reference polygon in this signature is flat or has its
only kink at the midpoint. Any departure from the reference
must increase the midpoint of a mixed orbit, by convexity and
the matching endpoints. Polarization gives the same midpoint
excess on the conjugate orbit, since
\(\Delta_{\bar\lambda,j}=\Delta_{\lambda,4-j}\).
Equation (1) excludes such a failure in either type, where its
cost is at least \(2t\).

Suppose henceforth \(c\in\langle\gamma\rangle\), so that \(t\)
is even and \(\gamma^{t/2}=c\). If \(\epsilon=0\) then \(d=\gamma\)
has order \(t\) and its orbit contains both mixed labels.
If \(\epsilon=1\), write \(d=\gamma^{t/2+1}\).
Its order is \(t\) when \(4\mid t\), with both labels in
one orbit, and \(t/2\) when \(t\equiv2\pmod 4\), with one
label in each of two orbits. In all cases the fraction of
mixed labels on a mixed orbit is \(2/t\).

Let \(\pi\) denote Frobenius at \(v\). Conjugation by
\(\pi^{t/2}\) exchanges each embedding with its conjugate
and commutes with \(\pi^t\), so conjugate blocks have
identical Newton slopes. With polarization, their polygons
are self-dual and each has slope sum two. The reference
endpoints match. Therefore every nonmixed block is flat of
slope \(1/2\), and the mixed blocks have reference slopes
\[
    \left(\tfrac12-\tfrac1t,\tfrac12-\tfrac1t,
          \tfrac12+\tfrac1t,\tfrac12+\tfrac1t\right).
                                                                    \tag{2}
\]
Indeed the jump between the second and third reference slopes is
the mixed fraction \(2/t\). The perturbation obtained by
averaging these two slopes is numerically allowed by (1).
It raises the midpoint by \(1/t\) and costs \(t\) including all
conjugates of the perturbed block. The next argument rules out
any positive midpoint excess at a density-one set of places.

## Even valuation from crystalline Frobenius

Let \(M\) be the first crystalline lattice of \(A_v/\mathbb F_p\)
and base change it to a finite unramified Witt ring \(W\) where
the coefficient algebra splits and the geometric action is
defined. Primes dividing the degree of the polarization or the
order index have been excluded. Thus \(M_W\) has a perfect
alternating pairing and is the direct sum of integral
coefficient lattices \(M_\sigma\) of rank four. The Rosati
condition pairs \(M_\sigma\) perfectly with \(M_{\bar\sigma}\).
On \(M_W\) the crystalline semilinear Frobenius is the linear
action of \(\pi\) followed by Witt Frobenius on scalars with
respect to the original \(\mathbb Z_p\)-lattice. On first
cohomology the linear action is by pullback, which reverses
the action on endomorphisms. Either direction for the induced
\(\gamma\)-permutation gives the same exchange by
\(\gamma^{t/2}=c\).

Semilinear Frobenius induces maps of the coefficient
lattices around \(\theta\)-orbits with elementary divisors
in \(\{1,p\}\). The number of \(p\) divisors is the
corresponding Hodge multiplicity, up to the covariant
versus contravariant convention which exchanges complementary
multiplicities. For example, the quotient of the target
lattice by the Frobenius image is killed by \(p\) by the
Frobenius–Verschiebung relation. Its dimensions in the
coefficient blocks are those of the corresponding Hodge
quotient in integral Dieudonné theory. All these dimensions
and their complements are even in the present signature.
The scalar Witt Frobenius is a semilinear lattice isomorphism
that permutes the coefficient idempotents, so it changes no
determinant valuations. Each one-step block map of the
linear operator \(\pi\) thus has even determinant valuation.
So does each block map of \(u=\pi^{t/2}\), by multiplication
along its \(t/2\) steps.
[Moonen, §§1.1.3 and 1.3.1 for the integral Dieudonné type](https://www.numdam.org/article/ASENS_2004_4_37_2_223_0.pdf),
[CFLMPT, §3.1 and Proposition 3.6 for the Hodge-orbit comparison](https://arxiv.org/html/2508.11174v1),
[Moonen, §§1.1–1.3 for the Dieudonné type and reference polygon](https://www.numdam.org/article/ASENS_2004_4_37_2_223_0.pdf).

For any \(\sigma\), take dual integral bases for
\(M_\sigma,M_{\bar\sigma}\) under the alternating pairing.
On this sum write
\[
 u=\begin{pmatrix}0&C\\ B&0\end{pmatrix}.
\]
Its similitude is \(p^{t/2}\), whence
\(C=-p^{t/2}B^{-\mathsf T}\). Consequently
\[
 Q_\sigma
  :=\det\!\left(1+\frac{\pi^t|_{M_\sigma}}{p^{t/2}}\right)
  =\frac{\det(B^{\mathsf T}-B)}{\det B}
  =\frac{\operatorname{Pf}(B^{\mathsf T}-B)^2}{\det B}.
                                                                    \tag{3}
\]
Since \(W\) is unramified and \(v_p(\det B)\) is even,
\(v_p(Q_\sigma)\) is even or infinite.

For a mixed block its normalized slope list
\(t(s_{\lambda,i}-1/2)\) has the form \((-a,-b,b,a)\) with
\(a\ge b\ge0\). A nonzero midpoint excess in (2) gives
\(a+b<2\). The valuations of the factors in \(Q_\sigma\)
are bounded below by \(\min(0,t(s_{\lambda,i}-1/2))\).
Their sum is at least \(-a-b>-2\). With (3) this implies
\(v_p(Q_\sigma)\ge0\). The same conclusion holds on the
conjugate block at a distinct coefficient place because its
slopes coincide. Nonmixed blocks already have normalized
slopes zero. Thus when the unique mixed pair fails to attain
the reference, every coefficient embedding contributes a
\(p\)-integral \(Q_\sigma\).

## A bounded rational function of Frobenius

On each algebraic monodromy component whose field action
\(\gamma\) has \(c\in\langle\gamma\rangle\), use the regular
function
\[
   J(g)=\prod_\sigma
        \det\!\left(1+\frac{g^t|_{V_\sigma}}{\chi(g)^{t/2}}\right),
                                                                    \tag{4}
\]
where \(\chi\) is the similitude character. The union of
these components is invariant under conjugation, and (4)
is a class function there. At Frobenius its value is the
norm from \(E\) of the relative determinant of
\(1+\pi^t/p^{t/2}\). The relative characteristic polynomial
of \(\pi^t\) has integral coefficients in \(E\),
by the Tate and crystalline compatibility of the relative
polynomials. Hence \(J(\pi)\) is rational with denominator
a power of \(p\). Its archimedean absolute value is at
most \(2^{8e}\), since every normalized Weil root in (4)
has modulus one. If the reference fails, (3) shows that
\(J(\pi)\) is an integer in this fixed finite interval.
[Relative-polynomial compatibility](CM_NEWTON_DISCONNECTED_SPECTRUM_INEQUALITIES.md#valuations-and-coefficient-places).

To see that these levels have density zero, fix a component
and write \(D=(G^\circ)^{\mathrm{der}}\). The projective
hypothesis implies every simple adjoint factor of \(D\) is
of type \(A_3\), and each block sees one such factor.
The exact-power argument picks a representative \(g_0\)
whose action preserves a pinning and whose power
\(g_0^t\) acts by a scalar \(b_\sigma\) on each block.
Let \(\nu\) be the cocharacter given by the sum
of all positive coroots of \(D\). It is fixed by the
pinned action, and its weights in each block are
\(3,1,-1,-3\), since the block is standard or dual for
its factor. Thus \((g_0\nu(z))^t=g_0^t\nu(z^t)\) and
\(\chi(g_0\nu(z))=\chi(g_0)\). Along this curve (4) is
a product of factors
\(1+b_\sigma \chi(g_0)^{-t/2}z^{tw}\) with
\(w\in\{3,1,-1,-3\}\). All are nonzero Laurent
polynomials, and their product has positive highest degree
and negative lowest degree. It is therefore nonconstant
on every indicated component. Algebraic Chebotarev
excludes the finitely many integer levels on a set of
density zero. Together with (1) this proves the
\(\mu\)-ordinary assertion.
[Rajan, Theorem 3](https://repository.ias.ac.in/38317/1/2-PUB.pdf).

## Ordinary density and imaginary quadratic coefficients

If \(\gamma=1\), the oriented Hodge multiplicities are
constant on the \(\theta\)-orbits, so their reference polygons
are ordinary. If \(c\notin\langle\gamma\rangle\) and \(t>1\),
the jump in each mixed block of the reference is \(1/t\).
That jump is incompatible with a polygon with slopes zero
and one. If \(c\in\langle\gamma\rangle\) and \(t>2\), the mixed
reference slopes in (2) are fractional. For \(t=2\), the mixed
label fraction \(2/t\) is one, so the mixed orbits exhaust only
the two mixed labels. When \(e>1\) there is then a nonmixed
coefficient place with slope \(1/2\).
Thus for \(e>1\) ordinary reduction occurs at the identity
class of \(\operatorname{Gal}(L/K)\) outside a set
of density zero. For \(e=1\), \(t\) is one or two and
the reference is ordinary in both cases, using (2) for
the nontrivial action.

Finally, if \(e=1\) and \(E\) is the full geometric rational
endomorphism algebra, choose \(\ell\) split in a finite
extension that defines the action and contains \(E\).
Faltings gives two nonisomorphic irreducible oriented
representations of \(G^\circ\), both of dimension four.
As the signature is \((2,2)\), the relative determinant
characters have infinity type \((2,2)\), so on \(G^\circ\)
both equal \(\chi^2\): after division by the square
of the cyclotomic character their algebraic Hecke
characters have trivial infinity type and hence finite
image. If \(z,w\) are the scalar characters of the
connected center on the two blocks, polarization gives
\(zw=\chi\) and the determinant gives \(z^4=\chi^2\).
Connectedness forces \(z=w\). If either derived block
were self-dual, polarization would then identify the
two complete \(G^\circ\)-representations, contrary to
Faltings. Among faithful irreducible four-dimensional
representations of semisimple Lie algebras over an
algebraically closed field of characteristic zero, only
the standard and dual of \(\mathfrak{sl}_4\) are
nonself-dual. The others are the standard of
\(\mathfrak{sp}_4\), the tensor product of two standards of
\(\mathfrak{sl}_2\oplus\mathfrak{sl}_2\), and
\(\operatorname{Sym}^3\) of the standard of
\(\mathfrak{sl}_2\). Thus the required projective image
is \(\mathrm{PGL}_4\).
[Fité, Proposition 14](https://arxiv.org/pdf/2005.10185v3),
[Faltings, Theorems 3 and 4](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf),
[Hui, §3.1](https://arxiv.org/pdf/2208.04004).

**Corollary.** Suppose the full geometric rational endomorphism
algebra of an abelian fourfold \(A/K\) is imaginary quadratic. Its
ordinary density is one in signature \((2,2)\). In unordered
signature \((1,3)\) the ordinary density is \(1/[L:K]\), where \(L\)
defines the geometric endomorphisms. The only possible densities
are one and one half, and both occur.

The only other unordered signature would be \((0,4)\). In that
signature each coefficient block in rational Betti cohomology has
a single Hodge type, so every \(E\)-linear rational endomorphism
preserves the Hodge decomposition. Clearing denominators gives
complex-torus endomorphisms, which are algebraic and descend to
the algebraic closure of \(K\). This would give an endomorphism
algebra containing \(M_4(E)\), contrary to the hypothesis. The
density for \((1,3)\) and examples with quadratic \(L/K\) follow
from the simple-signature theorem. Base extension to \(L\) then
also realizes density one in that signature.
[Milne, *Complex Multiplication*, Propositions 6.1 and 7.9](https://www.jmilne.org/math/CourseNotes/CM.pdf),
[Simple-signature theorem and examples](CM_NEWTON_SIMPLE_SIGNATURE_NONABELIAN.md#examples-with-ordinary-density-12).

For any CM field \(E\), the real PEL construction in the
[simple-signature examples](CM_NEWTON_SIMPLE_SIGNATURE_NONABELIAN.md#examples-with-ordinary-density-12)
also realizes the hypotheses here. Use a diagonal
hermitian space of rank four with signature \((2,2)\)
at one place and definite signature everywhere else.
The same conjugation of the lattice, action labels and
level gives a real point of the descended fine moduli
space. The monodromy argument uses the unique
noncompact real factor, now \(\mathrm{SU}(2,2)\).
It gives all the independent rank-four unitary factors,
and the same real Hilbert specialization yields an
absolutely simple fiber with full geometric algebra \(E\)
and a quadratic endomorphism field. For \(e>1\) this
produces ordinary density \(1/2\), with slope \(1/2\) on
the nonmixed factors at the nontrivial Frobenius class.

## An explicit genus-four Jacobian over \(\mathbb Q\)

Let \(C\) be the smooth projective model of
\[
   y^3=A(x)B(x)^2,\qquad
   A=x^3+x+1,\qquad B=x^3-2x^2-2x+1 .
\]
Its Jacobian \(J\) has full geometric endomorphism ring
\(\mathbb Z[\zeta_3]\), with endomorphism field
\(\mathbb Q(\zeta_3)\), and has ordinary reduction at a set
of rational primes of natural density one.

The branch discriminants are \(-31\) and \(125\), and the
resultant is \(31\). Thus the six branch points are distinct.
Infinity is unramified since the degree is nine. Riemann–Hurwitz
gives genus four. Away from \(3,5,31\) the normalized cover has
good reduction. After an étale extension at a branch point,
a unit can be absorbed into \(y\) and the cover has the
local form \(x-a=s^3,\ y=s\) or \(s^2\). Away from the
branches the cover is étale, including the chart at infinity
with coordinate \(y/x^3\).

The deck transformation acts through \(E=\mathbb Q(\zeta_3)\)
on \(H^1\), since its invariant cohomology is the cohomology
of the quotient \(\mathbb P^1\). It preserves the canonical
polarization of \(J\). The forms \(x^j\,dx/y\) and
\(x^j B(x)\,dx/y^2\) for \(j=0,1\) are holomorphic. At a
branch point \(dx\) has order two and \(y\) has order one or
two. At a zero of \(B\) the extra factor \(B(x)\) has order
three. At infinity \(x\) has a simple pole and \(y\) a pole
of order three. The forms are independent in their two
deck eigenspaces and give signature \((2,2)\).

For \(p=7,13\), fix the primitive cube root \(2,3\) in
\(\mathbb F_p\), respectively. For \(q=p^r\), let \(n_i\)
count \(x\in\mathbb F_q\) with \(A(x)B(x)^2\ne0\) for which
\((A(x)B(x)^2)^{(q-1)/3}\) equals the \(i\)-th power of
that root. With the deck action labelled by
\(\zeta=\zeta_3\), the \(r\)-th Frobenius power trace on the
\(\zeta\)-block is \(-1-\sum_{i=0}^2n_i\zeta^i\).
Indeed, for the composite of \(q\)-power Frobenius with
the \(j\)-th deck power there are three fixed points over
each counted \(x\) of character index \(-j\) and one over
each rational branch point. The three points at infinity
are fixed just when \(j=0\). Fourier projection of the
Lefschetz trace gives the formula, with the term \(-1\)
from infinity.

The exact character counts and Newton identities in the
[PARI/GP verifier](../scripts/verify_rank_four_pfaffian_curve.gp)
give the relative Frobenius polynomials
\[
\begin{aligned}
 P_7(X)&=X^4+(2-\zeta)X^3+(3-3\zeta)X^2
                  +(7-14\zeta)X-49\zeta,\\
 P_{13}(X)&=X^4+(5+\zeta)X^3+16X^2
                  +(52-13\zeta)X+169 .
\end{aligned}
\]
The other relative polynomial is the conjugate \(\bar P_p\).
The absolute middle coefficients of \(P_p\bar P_p\) are
\(167\) and \(1101\), prime to the respective residue
characteristics. The reductions are therefore ordinary.
Each relative determinant has \(p\)-adic valuation two, so
the four slopes in each block are \(0,0,1,1\).

To certify the geometric endomorphism algebra, the
following squarefree reductions give the indicated
factor degrees over the residue field \(\mathbb F_\ell\).

| \(\ell\) | Image of \(\zeta\) | \(P_7\) | \(P_{13}\) |
| ---: | ---: | :---: | :---: |
| 73 | 8 | \(1,3\) | \(4\) |
| 79 | 23 | \(4\) | \(1,3\) |

Thus each quartic is irreducible over \(E\), and its Galois
group contains a four-cycle and a three-cycle. A subgroup
of \(S_4\) containing both has order divisible by 12 and
contains an odd permutation. It is therefore \(S_4\).
Each \(\bar P_p\) is distinct from \(P_p\), so
\(P_p\bar P_p\) is irreducible over \(\mathbb Q\).
For any root \(\alpha\) of \(P_p\), its field
\(K_p=\mathbb Q(\alpha)\) contains \(E\) and has no proper
intermediate extension of \(E\), since the point stabilizer
\(S_3\) is maximal in \(S_4\).

Let \(\alpha_i\) be the four roots of \(P_p\). The roots of
\(\bar P_p\) are \(p/\alpha_i\) by polarization. Fix a
\(p\)-adic place of the splitting field. Two of the \(\alpha_i\)
have valuation zero and two have valuation one. If
\(\alpha_i/\alpha_j\) for \(i\ne j\) were a root of unity,
its valuations after all \(S_4\) permutations would
vanish. A permutation sending \(i,j\) to roots of
different valuations contradicts this. The same argument
applies to pairs from \(\bar P_p\). A ratio between the two
quartics would be \(\alpha_i\alpha_j/p\). If \(i=j\) its
valuation is already nonzero, and if \(i\ne j\) one can
send both roots to valuation zero. Thus the eight roots
stay distinct after every positive power. Tate's
endomorphism theorem over each finite extension of
\(\mathbb F_p\) gives geometric rational endomorphism
algebra \(K_p\).
[Tate, Theorem 2(c)](https://pazuki.perso.math.cnrs.fr/index_fichiers/Tate66.pdf).

The absolute polynomials \(P_7\bar P_7\) and
\(P_{13}\bar P_{13}\) are squarefree modulo \(2\), with
factor degrees \(4,4\) and \(2,6\), respectively. By
Dedekind factorization the fields \(K_7,K_{13}\) are
therefore nonisomorphic. Good-reduction specialization
embeds the characteristic-zero geometric endomorphism
algebra of \(J\) into both \(K_p\), carrying the deck action
to the corresponding copy of \(E\). Each \(K_p/E\) has no
intermediate fields, so the only possible common endomorphism
algebra is \(E\). Its integral order contains
\(\mathbb Z[\zeta_3]\), which is already maximal. The deck
action is defined over \(\mathbb Q(\zeta_3)\) and conjugated
to its inverse by the nontrivial Galois element, so that
quadratic field is the endomorphism field. The theorem now
gives density-one ordinary reduction over \(\mathbb Q\).
[Costa–Mascot–Sijsling–Voight, equation (7.3.2)](https://arxiv.org/pdf/1705.09248).

The verifier also checks ordinary reduction at the inert
prime \(11\), where Frobenius acts nontrivially on \(E\).
The point counts are \(12,132,1332,14172\) in extension
degrees one through four. The middle coefficient is
\(-105\equiv5\bmod 11\).
From the repository root, run the verifier.

```sh
gp -fq research/foundations/scripts/verify_rank_four_pfaffian_curve.gp
```

On Windows with PARI/GP in WSL, prefix the command by
`wsl`. The reference run used PARI/GP 2.15.4.
