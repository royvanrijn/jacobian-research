# Erdős #404: factorial sums and unbounded 5-adic cancellation

For an integer \(a\ge1\) and a prime \(p\), let \(f(a,p)\) be the supremum of \(v_p(\sum_{i\in I}i!)\) over finite sets of distinct positive indices \(I\) with \(\min I=a\); allow the value \(\infty\). [Erdős #404](https://www.erdosproblems.com/404) also asks for an infinite increasing sequence of indices whose partial factorial sums have valuations tending to infinity. The concrete open target here is
\[
                         f(1,5)=\infty. \tag{404.1}
\]
Reported finite witnesses reach valuation at least 20,000, but they have not been imported or replayed in this repository. They do not prove (404.1).

## Uniform lifting obligation

For each precision \(j\), track the reachable residues of sums of **distinct** factorials modulo \(5^j\), with index 1 retained. Find an invariant that (i) permits cancellation of the target residue, (ii) survives the passage from \(j\) to \(j+1\), and (iii) accounts for every newly exposed carry, including at indices around \(5^j\). A finite automaton constructed after fixing a maximum precision is not a uniform invariant. The [existing digit-orbit](image-rigidity/MULTIVARIATE_FACTORIAL_DIGIT_ORBITS.md) and [5-adic automaton](image-rigidity/MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md) notes concern factorial coefficients in polynomial rigidity; their finite residue mechanisms are method analogies, not a proof of (404.1).

If (404.1) is proved, compactness of \(\{0,1\}^{\mathbb N}\) gives an infinite set of indices containing 1 whose 5-adic factorial sum is zero: choose finite cancelling sets at unbounded precisions, take a coordinatewise convergent subsequence, and use that \(v_5(i!)\to\infty\). The limit set cannot be finite, because a finite sum of positive factorials is nonzero. Its increasing enumeration has partial-sum valuations tending to infinity. This implication does not supply the missing cancellations.

## Small completed special case

For every \(a\ge1\) with \(a\equiv6\pmod9\), one has
\[
                         f(a,3)=v_3(a!)+1. \tag{404.2}
\]
After division by \(a!\), the terms at \(a+1,a+2\) are \(7,2\pmod9\); every later term is \(0\pmod9\). The four subset residues are \(1,8,3,1\pmod9\), so valuation one is attained and valuation two is impossible. This is an elementary internal observation; novelty and external review have not been assessed.
