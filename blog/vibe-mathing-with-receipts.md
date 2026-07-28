---
title: "Vibe mathing, with receipts"
subtitle: "What happened when I treated mathematical research like an unfamiliar codebase—and refused to trust the output"
author: "Roy van Rijn"
date: "2026-07-28"
status: "Draft for publication after the Zenodo records are live"
description: "A developer's account of using AI, exact computation, Lean, and repeated falsification to move from one checkable counterexample to two constructive papers about Keller fibers."
---

> **Paper I:** *Over Characteristic Zero, Every Finite Étale Algebra of Rank at Least Three Is a Full Keller Fiber*  
> Zenodo record: [{{ZENODO_FIBERS_RECORD_URL}}]({{ZENODO_FIBERS_RECORD_URL}})  
> DOI: `{{ZENODO_FIBERS_DOI}}`  
> Direct PDF: [{{ZENODO_FIBERS_PDF_URL}}]({{ZENODO_FIBERS_PDF_URL}})
>
> **Paper II:** *Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map*  
> Zenodo record: [{{ZENODO_HASSE_RECORD_URL}}]({{ZENODO_HASSE_RECORD_URL}})  
> DOI: `{{ZENODO_HASSE_DOI}}`  
> Direct PDF: [{{ZENODO_HASSE_PDF_URL}}]({{ZENODO_HASSE_PDF_URL}})
>
> **Code and formalization:**  
> Repository: [{{REPOSITORY_URL}}]({{REPOSITORY_URL}})  
> Archived release: [{{ZENODO_CODE_RECORD_URL}}]({{ZENODO_CODE_RECORD_URL}})  
> DOI: `{{ZENODO_CODE_DOI}}`  
> Commit: `{{PUBLICATION_COMMIT_SHA}}`

# Vibe mathing, with receipts

I am a developer.

Right now it is difficult to ignore the argument around “vibe coding.”
People who could not have written an application in the traditional way are
describing what they want to an AI, running the result, reporting what broke,
and iterating. Some of the resulting software is fragile. Some of it is
surprisingly capable. Either way, the phenomenon is real.

The uncomfortable part for experienced developers is not only that code is
being generated. It is that the people producing it do not necessarily read
or understand every line. They may never have been able to write those lines
themselves.

As a developer, this feels less alien to me than it might sound. Large
software systems have always exceeded the understanding of any one person.
We navigate them through interfaces, types, tests, logs, debuggers,
reproducible builds, and carefully chosen experiments. None of these makes a
program correct. Together, however, they let us ask much better questions
than “does this code look convincing?”

I began to wonder whether something similar was possible in mathematics.

Could I enter a field in which I did not have the expected training, use AI
to navigate definitions and possible constructions, generate exact data,
loop through failed ideas, and gradually replace intuition with artifacts
that could be checked? Could I begin exploring before I understood the whole
territory in the traditional order?

Or, stated more provocatively: could I vibe-math my way to a genuine
mathematical frontier?

The experiment was never supposed to be “can AI produce a paper that looks
like mathematics?” We already know that it can. That is precisely the
problem. A plausible-looking proof is now almost free, while the time needed
to find the hidden gap still belongs to a human reader.

The real question was:

> Can AI-assisted exploration produce a result that survives exact
> computation, independent implementations, formal proof, explicit statement
> boundaries, and eventually human mathematical review?

I now have two papers that I believe give a meaningful partial answer. They
are not externally refereed yet, and their novelty claims remain open to
specialist correction. But they have passed a much stronger internal process
than “the model said it was true.”

This is the story of how I got there.

## Why I started with counterexamples

My first instinct was to search for counterexamples.

This was not accidental. Counterexamples have a useful asymmetry: they can be
hard to find but easy to challenge. If someone claims that a polynomial map
has constant Jacobian and sends two distinct points to the same target, I do
not need to trust the story that produced it. I can differentiate the map,
expand the determinant, substitute the points, and compare the outputs.

The witness is concrete. It can be attacked from several directions. A
second computer algebra system can repeat the calculation. A short
dependency-free script can check the integers. A proof assistant can verify
the identities from definitions.

That is an unusually good match for the current AI moment. AI can generate
many candidates, including many bad ones. Exact witnesses let us discard the
bad candidates without debating whether their explanations sound
intelligent.

Proofs and constructions are different. A twenty-page argument is not one
object with one obvious test. Its correctness is spread across definitions,
quantifiers, genericity assumptions, field extensions, literature
interfaces, and sentences that may quietly claim more than the lemmas
support. AI can make such an argument look smooth long before it is sound.

This difference became the central difficulty of the project. Searching for
a counterexample was the easy-to-check beginning. Moving from a witness to a
general construction, and from a construction to a precise paper, was much
harder.

## The seed I did not discover

The starting point was an explicit three-dimensional polynomial Keller map.
A Keller map is a polynomial map whose Jacobian determinant is a nonzero
constant. For generations, the Jacobian conjecture asked whether every such
map in characteristic zero must be polynomially invertible.

In July 2026, Levent Alpöge publicly announced an explicit counterexample,
crediting Akhil Mathew with suggesting the problem and Fable with producing
the example. The determinant and collision were subsequently formalized
independently in Isabelle/HOL by André Ramos, Davi Hulak, and Ruy de
Queiroz.

- [Alpöge's announcement](https://x.com/__alpoge__/status/2079028340955197566)
- [Independent Isabelle/HOL formalization](https://isa-afp.org/entries/Jacobian_Counterexample.html)
- [Terence Tao's explanation of the construction](https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/)

I did not discover that cubic map. My work began with a different question:

> What is this map the first example of?

The original witness was already exact. Its Jacobian and collision could be
checked directly. I reproduced those calculations, normalized the
determinant to one, and began changing coefficients and coordinates to see
which parts of the mechanism survived.

Most changes failed. That was useful. A failed symbolic identity is much
more informative than an AI explanation of why a construction “should”
work.

The first real milestone was recognizing that the complicated
three-variable formula was hiding a one-variable inverse equation.

## Milestone 1: stop staring at the map

The map looked complicated in source coordinates. The inverse problem did
not.

After changing perspective, a source point could be represented by a marked
root \(S\) of a scalar equation

\[
E_{\Pi,B,C}(S)=0.
\]

The crucial quantity was its derivative

\[
D=\frac{\partial E_{\Pi,B,C}}{\partial S}.
\]

At first the derivative looked like an auxiliary expression. It turned out
to be the organizing principle of the construction.

It did two jobs at once:

1. \(D\) was the unit needed to reconstruct the original source coordinates
   from the marked root.
2. The same \(D\) was the factor that cancelled from the Jacobian
   calculation.

That was the moment the example stopped looking like a miraculous cubic
formula. The inverse equation was the real object. The large polynomial map
was its pullback into coordinates where the denominators cancelled.

This was also the first place where the developer's way of working helped
me. I did not initially understand the entire geometry. I could still track
the same quantity through several representations, generate expansions,
compare identities, and ask which invariant appeared in every successful
version.

The derivative kept appearing. Eventually I understood why.

## Milestone 2: choose the fiber before the map

Once the inverse equation was visible, the direction of the problem could be
reversed.

Normally one chooses a polynomial map and then studies its fibers. I asked
whether I could choose a finite fiber first and build the map around it.

Suppose \(P(T)\) is a separable polynomial of degree \(N\). Choose a
translation \(a\) and form

\[
G(S)=P(a+S)-P(a).
\]

The translated polynomial has a marked root at \(S=0\). Feeding \(G\) into
the inverse construction makes the distinguished inverse equation become
the polynomial \(P(a+S)\). A root of \(P\) then determines a source point,
and the reconstruction sends that source point back to the root.

This sounds almost too simple when written afterward. It was not simple to
make precise.

It is easy to arrange for the desired roots to appear among the points of a
fiber. That is not enough. I needed to prove all of the following:

- the displayed formulas are polynomial rather than merely rational;
- the Jacobian determinant is exactly one;
- every root reconstructs to a source point;
- every source point in the fiber comes from a root;
- this equivalence works over every commutative test algebra, not only over
  an algebraic closure;
- the literal fiber scheme is \(\operatorname{Spec} K[T]/(P)\);
- the actual function-field degree of the map is \(N\);
- the fiber has rank \(N\), so it contains every generic inverse sheet.

That last equality is what I call fullness:

\[
\operatorname{rank}_K F^{-1}(y)
=N
=\operatorname{gdeg}(F).
\]

The distinction matters. A list of \(N\) points is a computation. A natural
identification of the complete fiber algebra, together with the actual
generic degree of the map, is a theorem.

The resulting statement is:

> Over every characteristic-zero field \(K\), every separable polynomial
> \(P\in K[T]\) of degree \(N\geq 3\) produces a polynomial map
> \(F_P:\mathbb A^3_K\to\mathbb A^3_K\) with Jacobian determinant one,
> geometric degree \(N\), and a distinguished full fiber naturally
> isomorphic to \(\operatorname{Spec} K[T]/(P)\).

Every finite étale algebra over an infinite field is generated by one
element. In characteristic zero, that turns the polynomial statement into
the paper's abstract conclusion:

> Every finite étale algebra of rank at least three over a
> characteristic-zero field occurs as a full fiber of a Jacobian-one
> polynomial self-map of affine three-space.

This became the first paper.

## Milestone 3: the proof assistant stopped being decoration

At the beginning, my Lean knowledge was basic. I used AI to help translate
definitions, locate library lemmas, interpret compiler errors, generate
proof attempts, and reorganize statements when the first formulation was
too vague to formalize.

This was looped prompting in the most literal sense:

1. ask for a formal statement;
2. compile it;
3. read the error;
4. reduce the claim;
5. expose the missing hypothesis;
6. try again;
7. independently compare the result with the paper.

Lean did not tell me that the theorem was interesting or new. It did
something narrower and indispensable: it made it difficult to slide between
similar-looking statements.

For example, these are not interchangeable:

- exhibiting roots and representing the whole fiber;
- bounding the inverse polynomial degree and proving the actual
  function-field degree;
- choosing a polynomial presentation and choosing one naturally;
- working in characteristic zero and working in every characteristic other
  than two;
- proving an equality over the base field and proving it naturally after
  scalar extension.

The formalization forced these boundaries into public declarations.

The characteristic-zero polynomial-presentation theorem is now formalized
end to end in Lean, including the literal natural fiber, finite étaleness,
rank, the function-field comparison, geometric degree, admissible
translation, monogenicity, and the abstract finite-étale corollary. The
publication certificate contains no `sorry`, no `admit`, and no
project-specific axioms.

The broader supplied-presentation theorem in characteristic different from
two is proved in the paper but is not yet formalized end to end. The
rank-two exclusion uses classical results of Campbell, Razar, and Wright and
also lies outside the Lean certificate.

Those are not footnotes to hide. They are the boundary of what I am asking a
reader to trust.

## Milestone 4: one arithmetic fiber

The construction can carry arbitrary finite étale arithmetic into a Keller
fiber. A particularly striking example comes from the classical polynomial

\[
(X^3-19)(X^2+X+1).
\]

It has no rational root, but it has a root over the real numbers and over
every \(p\)-adic field \(\mathbb Q_p\). The two factors divide the local work
between them. The cubic supplies roots at one collection of primes; the
quadratic cyclotomic factor supplies them at the others.

The first paper turns this algebra into a complete degree-five Keller fiber.
But the map produced by the general construction depends on the polynomial.
Changing \(19\) to another parameter generally changes the map.

That left a more demanding question:

> Can one fixed polynomial Keller map contain infinitely many complete
> fibers that are everywhere locally soluble but have no rational point?

This is a small change in quantifiers and a large change in difficulty.

The first theorem says:

\[
\text{for every fiber algebra, there exists a map.}
\]

The new question asks:

\[
\text{there exists one map with infinitely many such fiber algebras.}
\]

## Milestone 5: finding the line

The natural arithmetic family was

\[
f_a(X)=(X^3-a)(X^2+X+1).
\]

Expanded in \(X\), the parameter \(a\) moves three coefficients:

\[
f_a(X)=X^5+X^4+X^3-aX^2-aX-a.
\]

That did not fit the target pencil I had. Repeated searches for a fixed map
kept producing formulas that were almost right and geometrically wrong.

Then came the smallest observation in the project:

\[
X=S-\frac12.
\]

Under this translation,

\[
X^2+X+1=S^2+\frac34.
\]

Now the moving part is

\[
-a\left(S^2+\frac34\right).
\]

Only the quadratic and constant coefficients move. Those are exactly the two
coefficients controlled by a line in the target of one quadratic-gauge map.

The isolated example had become a rational target line.

This is the kind of moment I had hoped AI-assisted exploration might make
possible. No individual computation was deep. The value came from keeping
many equivalent forms alive, testing them, and noticing which tiny change
made the coefficient motion match the geometry.

For the resulting fixed determinant-one map \(\Phi\), the targets are

\[
y_a=\left(-1,\frac{32a}{9},\frac{8a+1}{3}\right),
\]

and the complete fiber is

\[
\Phi^{-1}(y_a)
\simeq
\operatorname{Spec}\mathbb Q[X]/
\bigl((X^3-a)(X^2+X+1)\bigr).
\]

## Milestone 6: the prime family was not the real family

At first I concentrated on primes

\[
a=\ell\equiv 1\pmod 9.
\]

That already gives infinitely many Hasse failures. But when I rewrote the
local proof place by place, I discovered that primality was not doing the
work I had assigned to it.

The actual sufficient conditions are:

\[
a>1,\qquad
a\equiv1\pmod9,\qquad
a\notin\mathbb Q^3,\qquad
p\mid a\Longrightarrow p\equiv1\pmod3.
\]

The congruence-and-prime-support core is multiplicatively closed. The
noncube condition is imposed afterward, because perfect cubes would give a
rational root. This distinction mattered enough to formalize explicitly:
the core is a semigroup; the final admissible set is not claimed to be one.

For every admissible \(a\):

- the fiber is reduced, finite étale, and has rank five;
- it has a real point;
- it has a point over every \(\mathbb Q_p\);
- it has no rational point;
- its primitive projective target coordinates are
  \([9:-9:32a:24a+3]\);
- its height is exactly \(32a\);
- different parameters give different targets.

The first parameter is \(a=19\). The old example did not disappear. It
became the first visible point on the line.

The multiplicative family is much larger than the prime progression. A
character-twisted Selberg--Delange argument gives the asymptotic

\[
\#\{a:H(y_a)\leq B\}
\sim
\frac{G_3(1)}{96\sqrt{\pi}}\,
\frac{B}{\sqrt{\log B}},
\]

for an explicit positive Euler-product constant \(G_3(1)\).

A bounded dependency-free enumeration through \(a=10^6\) reproduces the
pinned artifact and finds 26,846 members of the full sufficient family.
That computation is a regression test, not the proof of the asymptotic.

Finally, the classical low-degree theory of intersective polynomials shows
that degree five is the first degree in which such a zero-dimensional Hasse
failure can occur. The fixed map therefore has the smallest possible
geometric degree.

This became the second paper.

## Milestone 7: formalizing the fixed-map family

The second Lean certificate proves considerably more than a few sample
values of \(a\).

For every admissible parameter, it verifies:

- the exact displayed fixed map and its determinant-one normalization;
- the centered inverse factorization;
- the function-field degree-five calculation;
- the literal quotient-algebra representation of the fiber;
- separability, finite étaleness, and rank five;
- the rational-point obstruction;
- existence of a real point;
- existence of a root over every \(\mathbb Q_p\);
- multiplicative closure of the congruence-and-support core;
- the prime subfamily;
- primitive target coordinates, exact height, and target distinctness.

Again, the certificate has no `sorry`, no `admit`, and no project-specific
axioms.

The Selberg--Delange asymptotic is still an ordinary mathematical proof. So
is the degree-minimality argument using the Berend--Bilu criterion. I have
made both interfaces explicit in the paper, but neither should be described
as Lean-verified.

This is where the difference between finding counterexamples and writing
papers became impossible to ignore.

The original collision could be checked by differentiating and
substitution. The two new papers required me to control entire families,
generic degrees, naturality, base change, local fields, analytic counting,
external literature, and the exact relationship between formal and informal
statements. There was no single green test that made all of that true.

The verification had to become layered.

## What “verified” means here

I do not use “verified” as one undifferentiated label.

| Claim layer | How it is checked | What the check does not establish |
|---|---|---|
| Foundational cubic determinant and collision | Direct exact calculation; independent Isabelle/HOL formalization | That every later generalization is correct |
| Generated polynomial identities and explicit examples | SymPy, Sage, Singular, dependency-free Python, and cross-format hash checks | Generic theorems outside the checked identities |
| Characteristic-zero prescribed-fiber theorem | Isolated Lean publication certificate | Novelty, importance, or the broader characteristic-\(\ne2\) statement |
| Fixed-map Hasse family, local points, and target height | Isolated Lean publication certificate plus exact audits | Selberg--Delange or degree minimality |
| Asymptotic count | Written Euler-product and character-twisted Selberg--Delange proof; bounded enumeration as regression | A machine-checked analytic-number-theory proof |
| Minimal geometric degree | Written reduction to classical intersective-polynomial theory | A new formalization of the external theorem |
| Novelty | Literature audit and explicit comparison with nearby work | Exhaustive priority or peer review |

The repository also has a machine-readable mathematical status file. It
distinguishes theorem, conditional result, computation, experiment, and open
problem. A successful bounded search is not promoted into a theorem.

This sounds administrative. In an AI-assisted project, I think it is part of
the mathematics.

## What the machines did not prove

Lean did not prove that I chose the most meaningful question.

Singular did not prove that the construction is new.

Passing tests did not prove that my definitions are the definitions the
community will care about.

An AI-assisted literature search did not prove that I found every relevant
paper.

Formalizing a statement did not prove that the prose surrounding it is
incapable of misleading a reader.

And none of these tools can perform the social and mathematical function of
an expert referee who understands the surrounding field.

The papers therefore disclose the use of generative AI. I used it during
exploration, symbolic-checker development, proof review, literature search,
Lean formalization, and manuscript editing. AI-generated suggestions were
not treated as evidence. References were checked against their original
sources, and claims were checked through written arguments, exact
computations, or the Lean declarations identified in the verification
appendices. I reviewed the manuscripts and take responsibility for them.

That is not a claim that I manually derived every line before allowing a
machine to touch it. I did not. It is a claim that the route by which an idea
was generated is separate from the evidence offered for its correctness.

## The receipts

The two papers are available as fixed Zenodo records:

1. [*Over Characteristic Zero, Every Finite Étale Algebra of Rank at Least
   Three Is a Full Keller Fiber*]({{ZENODO_FIBERS_RECORD_URL}}), DOI
   `{{ZENODO_FIBERS_DOI}}`.
2. [*Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller
   Map*]({{ZENODO_HASSE_RECORD_URL}}), DOI `{{ZENODO_HASSE_DOI}}`.

The archived code and formalization release is available at
[{{ZENODO_CODE_RECORD_URL}}]({{ZENODO_CODE_RECORD_URL}}), DOI
`{{ZENODO_CODE_DOI}}`. The corresponding public repository commit is
`{{PUBLICATION_COMMIT_SHA}}`.

The narrow reproduction commands are:

```bash
# Paper I: generated paper/Lean correspondence
.venv/bin/python scripts/verify_common_arithmetic_fibers_correspondence.py

# Paper II: fixed map and uniform Hasse family
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py

# Lean publication certificates
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.PaperCertificate
lake build FiniteEtaleKeller.FixedHassePaperCertificate
```

The repository contains more exploratory material than these two papers.
That material has mixed status and should not be read as one enormous claim.
The two publication certificates and their verification documents are the
intended narrow entry points.

## What I think I learned

The strongest version of “vibe mathing” is not mathematics without
understanding.

It is mathematics in which understanding does not always have to arrive
first.

I could begin by navigating representations, generating examples, asking AI
for possible bridges, and trying to falsify every bridge that appeared.
Exact computation let me reject bad directions. Formalization forced vague
claims to split into precise ones. Writing the papers forced me to explain
why the surviving path worked.

My understanding followed the artifacts. It did not precede all of them.

That changes what may be possible for people outside the normal research
pipeline. It does not abolish expertise. In fact, it makes the remaining
expert work easier to see. Once calculations and formal statements are
machine-checkable, the scarce human questions become:

- Is this the right theorem?
- Is the interface with the literature correct?
- Is an apparently small hypothesis carrying hidden mathematical content?
- Is the ordinary proof outside the formal boundary sound?
- Has the result already appeared in another language?

Those are exactly the questions on which I now want expert criticism.

I began with counterexamples because they were cheap to distrust and easy to
check. I ended with constructions because the first witness contained more
structure than I expected. Crossing from one to the other was far harder
than generating the original candidates, and it required more than prompting:
it required a growing system of definitions, tests, exact certificates,
failed searches, formal proofs, literature boundaries, and prose.

I do not know yet how the mathematical community will judge the two results.
That judgment cannot be automated.

But I do think the experiment has reached the point where the right response
is no longer “do we trust the AI?” Nobody should.

The right response is to inspect the statement, run the checks, read the
proof, and try to break it.

The vibes helped find the path. Verification decided which parts survived.

---

## Citation placeholders

Replace every `{{...}}` token before publication.

### Paper I

```text
Roy van Rijn, "Over Characteristic Zero, Every Finite Étale Algebra of
Rank at Least Three Is a Full Keller Fiber", Zenodo, 2026.
DOI: {{ZENODO_FIBERS_DOI}}
Record: {{ZENODO_FIBERS_RECORD_URL}}
```

### Paper II

```text
Roy van Rijn, "Quantitative Hasse-Principle Failures in the Fibers of a
Fixed Keller Map", Zenodo, 2026.
DOI: {{ZENODO_HASSE_DOI}}
Record: {{ZENODO_HASSE_RECORD_URL}}
```

### Code and formalization archive

```text
Roy van Rijn, "Jacobian Research: publication certificates for prescribed
Keller fibers and fixed-map Hasse failures", Zenodo, 2026.
Version: {{ZENODO_CODE_VERSION}}
DOI: {{ZENODO_CODE_DOI}}
Record: {{ZENODO_CODE_RECORD_URL}}
Commit: {{PUBLICATION_COMMIT_SHA}}
```
