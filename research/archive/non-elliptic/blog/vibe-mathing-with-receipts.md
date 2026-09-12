---
title: "Vibe mathing, with receipts"
subtitle: "How AI-assisted exploration found the exact dimensional boundary of the Generalized Vanishing Conjecture"
author: "Roy van Rijn"
date: "2026-08-04"
status: "Draft article accompanying the published GVC preprint"
description: "A developer's account of using AI, exact computation, Lean, and repeated falsification to prove the two-variable Generalized Vanishing Conjecture and find its first counterexample in dimension three."
---

> **Paper:** *The Generalized Vanishing Conjecture: The Two-Variable Theorem and the First Failing Dimension*<br>
> Zenodo record: [zenodo.org/records/21782342](https://zenodo.org/records/21782342)<br>
> DOI: [10.5281/zenodo.21782342](https://doi.org/10.5281/zenodo.21782342)<br>
> Repository source: [`papers/generalized-vanishing-two-variables/main.tex`](../papers/generalized-vanishing-two-variables/main.tex)<br>
> Lean proof: [github.com/royvanrijn/jacobian-research/tree/main/formal/gvc](https://github.com/royvanrijn/jacobian-research/tree/main/formal/gvc)<br>
> Formalization: the counterexample and failure in every dimension at least three are fully Lean-verified; nearly all supporting algebraic machinery and much of the binary proof infrastructure are checked as well

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

Vibe coders sometimes even seem to have a major advantage. Their
unwillingness to read the code makes it possible to work with several agents
in parallel, producing incredible speedups. However, I am not going to settle the
questions of viability or quality here.

Instead, I decided to take that instinct somewhere else: to vibe my way into
a field about which I had very limited knowledge, not burdened by knowledge,
deep understanding, or preconceptions.

As a developer, this feels less alien to me than it might sound. I have
spent a lot of time using tools like Codex, refining prompts, and running
parallel experiments. Large software systems have always exceeded the
understanding of any one person. We navigate them through interfaces, types,
tests, logs, debuggers, reproducible builds, and carefully chosen
experiments. None of these makes a program correct. Together, however, they
let us ask much better questions than “does this code look convincing?”

I began to wonder whether something similar was possible in mathematics.

Could I enter a field in which I did not have the expected training, use AI
to navigate definitions and possible constructions, generate exact data,
loop through failed ideas, and gradually replace intuition with artifacts
that could be checked? Could I begin exploring before I understood the whole
territory in the traditional order?

Or, stated more provocatively: could I vibe-math my way to a genuine
mathematical frontier?

The real spark was Levent posting the [JC(3) counterexample](https://x.com/__alpoge__/status/2079028340955197566).

That construction changed the question for me. Until then, I had been
thinking mostly about whether AI could help me navigate unfamiliar
mathematics. Suddenly there was a concrete mathematical object to inspect:
a three-variable construction connected to the Jacobian Conjecture, with
the kind of exact algebraic structure that could be tested, generalized,
formalized, and attacked from several directions.

It also suggested a sharper experiment. Instead of asking whether AI could
produce mathematics that looked plausible, I could ask whether an
AI-assisted process could locate the exact boundary between a theorem and
its first counterexample.

The experiment was never supposed to be “can AI produce a paper that looks
like mathematics?” We already know that it can. That is precisely the
problem. A plausible-looking proof is now almost free, while the time needed
to find the hidden gap still belongs to a human reader.

The real question was:

> Can AI-assisted exploration turn a concrete spark into a result that
> survives exact computation, independent implementations, formal proof,
> explicit statement boundaries, and eventually human mathematical review?

The Generalized Vanishing Conjecture became my test case. The result is a
sharp boundary: it holds for every constant-coefficient differential
operator in one or two variables, and it fails from dimension three onward.

This is the story of how Levent's post became a theorem, a counterexample,
and a collection of receipts.

## A conjecture with an infinite premise

Let \(\Lambda\) be a constant-coefficient differential operator and let
\(P\) be a polynomial. The Generalized Vanishing Conjecture asks whether

\[
\Lambda^m(P^m)=0\quad\text{for every }m\geq1
\]

forces

\[
\Lambda^m(QP^m)=0\quad\text{for all sufficiently large }m
\]

for every fixed polynomial multiplier \(Q\).

The statement grew out of work around the Vanishing Conjecture and the
Jacobian Conjecture. But what made it especially interesting to me was its
shape. The premise is an infinite sequence of exact identities, and the
conclusion is eventual vanishing after inserting any fixed multiplier.

This is exactly the kind of statement on which vibe mathing can go wrong.
A search can verify the first ten powers, or the first hundred, and still
prove nothing about the next one. An AI can explain why a visible pattern
“should” continue and bury the missing quantifier inside a polished proof.

So the first rule was severe:

> A bounded search was never a proof of GVC.

That did not make bounded searches useless. It changed what I allowed them
to mean.

## The long route through exact computation

I began classifying low-degree cases in two variables. The calculations
separated operator symbols by root type, exposed exceptional coefficient
branches, and tracked which weighted faces could survive the first several
pure identities.

Again and again, the same patterns appeared. Surviving supports became
one-sided. Moment ideals developed staircase-shaped radicals. A complicated
coefficient system would collapse onto a face where a simple weighted
degree deficit forced every mixed expression to vanish eventually.

Degree by degree, exact calculations closed polynomial degrees four, five,
six, and seven for arbitrary binary constant-coefficient operators. These
were real theorems in their stated ranges. They were not a license to draw
an all-degree curve through four data points.

The attempted general proof grew into a much larger machine: Hall
matchings, factorial packets, prime-power carries, characters, Graver
bases, and finite traces. AI was useful here in the same way it is useful in
an unfamiliar codebase. It could suggest representations, translate between
formulations, generate experiments, and help turn a failed idea into a
smaller exact question.

Many pieces of that machine became valid mathematics. The machine as a
whole did not close the conjecture. One promotion step kept moving away: I
could understand a packet after it had been isolated, but I could not prove
that the original infinite vanishing premise always exposed one fixed
packet in the required way.

This was an important part of the experiment. The repository accumulated
more positive results without quietly changing the status of the missing
step. The route was productive, but it remained incomplete.

Then the proof became shorter.

## Two moving Newton faces

Take the symbol of \(\Lambda\) and the support of \(P\). For every positive
weight

\[
w_s=(s,1),
\]

look at two faces:

- the lowest-weight face of the operator symbol;
- the highest-weight face of the polynomial.

The pure vanishing identity passes to these extremal faces. Hall's marriage
theorem orders their horizontal Newton intervals just after ordinary
degree. A separate prime-dilation argument shows that, while the two faces
have unequal weights, their intervals cannot overlap.

In two variables, this has a decisive consequence. The projection of each
face is an interval on a line. Two disjoint intervals cannot exchange order
without meeting.

Now define the lower operator envelope \(L(s)\) and the upper polynomial
envelope \(U(s)\). Their gap

\[
\Delta(s)=U(s)-L(s)
\]

is piecewise linear because both supports are finite. Just after \(s=1\),
Hall localization puts the operator interval strictly to the right of the
polynomial interval. Shifted-ray rigidity prevents the intervals from
overlapping, so they cannot pass one another while \(\Delta(s)>0\).

The gap cannot stay positive forever. After the last breakpoint its slope
is forced to be negative, so it must reach zero. At the first meeting, every
operator monomial lies on or above one common weighted threshold and every
polynomial monomial lies on or below it. A fixed multiplier contributes
only a bounded defect; it cannot repair a deficit that grows linearly with
the power \(m\).

That proves:

> Over a characteristic-zero field, every constant-coefficient
> differential operator in two variables satisfies the Generalized
> Vanishing Conjecture.

The short proof also explains the staircase radicals from the earlier
computer algebra. They were not mysterious survivors. They were finite
shadows of a geometric fact: the moving Newton intervals were forbidden to
overlap.

Years of imagined mathematical labor had collapsed into following two
piecewise-linear envelopes until they met.

## Then the conjecture broke

The binary theorem raised the obvious question: was two variables merely
the first case of a higher-dimensional principle?

No.

Put

\[
\rho=t^2+xy,\qquad A=\rho+x^2,
\]

\[
C=y\rho^2-2xt^2\rho-x^3t^2,
\qquad P=AC^2,
\]

and let

\[
\Delta=4\partial_x\partial_y+\partial_t^2,
\qquad \Lambda=\Delta^6.
\]

Then, for every \(m\geq1\),

\[
\Lambda^m(P^m)=0,
\qquad
\Lambda^m(x^2P^m)\ne0.
\]

The first identity supplies the complete pure premise. The second uses the
single fixed multiplier \(Q=x^2\) and fails at every power, not merely
infinitely often. This is a homogeneous counterexample in three variables.

Adding unused variables propagates the same failure to every larger finite
dimension. Combined with the binary theorem, the exact dimensional
classification is:

\[
\boxed{\operatorname{GVC}(n)\text{ holds if and only if }n\leq2.}
\]

The same boundary already holds when the differential operator is required
to be homogeneous.

The contrast is sharp. In two variables, Newton faces project to intervals
and inherit a total order. In three variables, higher-dimensional faces can
move around one another. The mechanism that makes the binary proof work is
precisely what disappears at the first failing dimension.

## From a witness to an all-order Lean theorem

The [public Lean development](https://github.com/royvanrijn/jacobian-research/tree/main/formal/gvc)
is not a small certificate attached to the end of the project. It currently
contains more than 3,100 lines across sixteen modules and 164 theorem or
lemma declarations, with no `sorry`, `admit`, or explicit `axiom`
declarations. We have now checked almost every part of the argument that can
be cleanly isolated into an algebraic interface.

Most importantly, the three-variable counterexample is formalized end to
end.

The development does not check a list of sample powers. It defines the
literal polynomials \(\rho,A,C,P\), the operator \(\Delta\), its sixth
power \(\Lambda\), and the multiplier \(Q\). It then proves the algebraic
machinery connecting them:

- coefficientwise semantics for constant-coefficient differential
  operators;
- composition of operator symbols and apolar top contraction;
- the cusp identity behind the construction;
- the Reynolds--phase identity by exact coefficient expansion;
- the endpoint coefficient extraction for every power;
- the zero pure coefficient and exact nonzero neighboring coefficient;
- base change from the rationals to every characteristic-zero field;
- padding from three variables to every larger finite dimension.

The resulting Lean theorem is unconditional: GVC fails over every
characteristic-zero field in every finite dimension at least three.

That is the entire negative half of the dimensional classification, fully
checked for all powers, fields, and finite dimensions in its stated scope.
Nothing in the counterexample is left as an assumed computational bridge.

The formalization goes considerably further. It checks the coefficientwise
definition of GVC, operator composition, top contraction, the full endpoint
coefficient ladder, formal beta evaluation, exact factorial valuations,
Reynolds expansion, Laurent phase extraction, coefficient-ring base change,
unused-variable padding, the winding--profile--radial degree formulas, and
the final piecewise-linear envelope crossing used by the positive theorem.

What remains outside Lean is now concentrated rather than diffuse: the
global Hall-localization and shifted-ray separation bridge that drives the
short binary proof, together with its no-reversal and common-threshold
termination argument. Those steps are proved in the paper. Everything
around that bridge has been reduced to explicit definitions, exact lemmas,
or independent degree-four through degree-seven regressions.

This combination matters to me. Exploration found the construction and the
proof architecture. Exact computation kept killing false branches.
Formalization then converted the entire negative half of the dimensional
classification, and much of the positive machinery, into declarations that
can be checked from definitions.

## What I think I learned

The strongest version of “vibe mathing” is not mathematics without
understanding.

It is mathematics in which understanding does not always have to arrive
first.

I could begin by navigating representations, generating examples, asking AI
for possible bridges, and trying to falsify every bridge that appeared.
Exact computation let me reject bad directions. Formalization forced vague
claims to split into precise ones. Writing the proof forced me to explain
why the surviving path worked.

My understanding followed the artifacts. It did not precede all of them.

That changes what may be possible for people outside the normal research
pipeline. It does not abolish expertise, and I would not say that it makes
the work easier for them. It creates a new experience: people can begin
exploring unfamiliar mathematics before they have acquired the usual
training, but they will also produce a great deal of low-quality material
that needs to be filtered.

AI may be useful for that first pass, helping rate, sort, and identify which
artifacts deserve closer human attention. But that is only a pre-check, not
a substitute for expert judgment. The remaining human questions become:

- Is this the right theorem?
- Is the interface with the literature correct?
- Is an apparently small hypothesis carrying hidden mathematical content?
- Is the ordinary proof outside the formal boundary sound?
- Has the result already appeared in another language?


Those are exactly the questions on which I now want expert criticism.

The GVC project forced the process in both directions at once: a general
positive theorem in two variables and a concrete failure in three. Crossing
between those modes was far harder than generating plausible candidates. It
required a growing system of definitions, tests, exact certificates, failed
searches, formal proofs, literature boundaries, and prose.

I do not know yet how the mathematical community will judge the result.
That judgment cannot be automated.

But I do think the experiment has reached the point where the right response
is no longer “do we trust the AI?” Nobody should.

The right response is to inspect the statement, run the checks, read the
proof, and try to break it.

The vibes helped find the path. Verification decided which parts survived.

## The receipts

- [Published GVC preprint on Zenodo](https://zenodo.org/records/21782342), DOI [10.5281/zenodo.21782342](https://doi.org/10.5281/zenodo.21782342)
- [Public Lean proof and formalization](https://github.com/royvanrijn/jacobian-research/tree/main/formal/gvc)
- [Repository source](../papers/generalized-vanishing-two-variables/main.tex)
- [Unrestricted binary GVC by Hall-envelope separation](../extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md)
- [Homogeneous three-variable counterexample](../extended-geometry/THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md)
- [Lean coverage and build instructions](../formal/gvc/README.md)
- [Reproduction catalogue](../REPRODUCE.md)

The formal development is pinned to Lean and Mathlib versions recorded in
the repository. From the repository root, build the GVC verification with:

```bash
make verify-gvc-lean
```
