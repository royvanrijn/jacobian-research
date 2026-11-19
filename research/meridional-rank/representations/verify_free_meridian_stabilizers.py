#!/usr/bin/env python3
"""Finite word controls for the exact free-meridian stabilizer theorem.

The unrestricted theorem uses its separate arc-cutting proof. These
standard-library calculations verify its literal band conventions,
framing, relations and strand-deletion left inverse on bounded inputs.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("This certificate requires Python assertions enabled.")

import argparse
import hashlib
import json
from pathlib import Path
import platform


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def reduced(word):
    out = []
    for letter in word:
        if out and out[-1] == -letter:
            out.pop()
        else:
            out.append(letter)
    return tuple(out)


def inverse(word):
    return tuple(-letter for letter in reversed(word))


def multiply(*words):
    return reduced(letter for word in words for letter in word)


def hurwitz(values, braid):
    values = list(values)
    for letter in braid:
        j = abs(letter) - 1
        a, b = values[j:j+2]
        values[j:j+2] = (
            [b, multiply(inverse(b), a, b)] if letter > 0
            else [multiply(a, b, inverse(a)), a]
        )
    return values


def band(i, j):
    assert 1 <= i < j
    return tuple(range(i, j)) + tuple(range(-(j-2), -i+1))


def lift(word, mobile):
    out = []
    for letter in word:
        j = abs(letter) - 1
        generator = band(mobile[j], mobile[j+1])
        out.extend(generator if letter > 0 else inverse(generator))
    return tuple(out)


def delete_protected(braid, n, protected):
    """Delete strand trajectories, tracking labels through every crossing."""
    labels = list(range(1, n+1))
    mobile = set(labels) - set(protected)
    out = []
    for letter in braid:
        j = abs(letter) - 1
        if labels[j] in mobile and labels[j+1] in mobile:
            k = sum(label in mobile for label in labels[:j]) + 1
            out.append(k if letter > 0 else -k)
        labels[j], labels[j+1] = labels[j+1], labels[j]
    assert all(labels[i-1] == i for i in protected)
    return tuple(out)


def framed(values, mobile, protected):
    if not mobile:
        return []
    out = []
    for j in mobile:
        prefix = tuple(i for i in range(mobile[0], j) if i in protected)
        out.append(multiply(prefix, values[j-1], inverse(prefix)))
    return out


def band_formula_controls():
    count = 0
    for n in range(2, 13):
        values = [(i,) for i in range(1, n+1)]
        for i in range(1, n):
            for j in range(i+1, n+1):
                W = tuple(range(i+1, j))
                expected = values.copy()
                expected[i-1] = multiply(W, (j,), inverse(W))
                expected[j-1] = multiply((-j,), inverse(W), (i,), W, (j,))
                got = hurwitz(values, band(i, j))
                assert got == expected
                assert hurwitz(got, inverse(band(i, j))) == values
                count += 1
    return count


def subset_controls():
    counts = {"protected_subsets": 0, "signed_generator_controls": 0,
              "framed_braid_relations": 0, "framed_commuting_relations": 0,
              "lifted_word_controls": 0}
    examples = []
    for n in range(1, 9):
        values = [(i,) for i in range(1, n+1)]
        for mask in range(1 << n):
            protected = {i for i in range(1, n+1) if mask & (1 << (i-1))}
            mobile = [i for i in range(1, n+1) if i not in protected]
            s = len(mobile)
            initial = framed(values, mobile, protected)
            counts["protected_subsets"] += 1
            for j in range(1, s):
                for sign in [-1, 1]:
                    word = (sign*j,)
                    literal = lift(word, mobile)
                    after = hurwitz(values, literal)
                    assert all(after[i-1] == (i,) for i in protected)
                    assert framed(after, mobile, protected) == hurwitz(initial, word)
                    assert delete_protected(literal, n, protected) == word
                    counts["signed_generator_controls"] += 1
            for j in range(1, s-1):
                left = lift((j, j+1, j), mobile)
                right = lift((j+1, j, j+1), mobile)
                assert hurwitz(values, left) == hurwitz(values, right)
                counts["framed_braid_relations"] += 1
            for i in range(1, s):
                for j in range(i+2, s):
                    assert hurwitz(values, lift((i, j), mobile)) == hurwitz(
                        values, lift((j, i), mobile))
                    counts["framed_commuting_relations"] += 1
            if s >= 2:
                words = [tuple(range(1, s)), tuple(range(s-1, 0, -1)),
                         tuple(j if j % 2 else -j for j in range(1, s)),
                         tuple(range(1, s))*3]
                if s >= 3:
                    words += [(1, 2, -1, -2, 1, -2), (1, 2, 1, -2, -1, -2)]
                for word in words:
                    literal = lift(word, mobile)
                    after = hurwitz(values, literal)
                    assert all(after[i-1] == (i,) for i in protected)
                    assert framed(after, mobile, protected) == hurwitz(initial, word)
                    assert delete_protected(literal, n, protected) == word
                    counts["lifted_word_controls"] += 1
                if n == 8 and protected == {1, 4, 5, 8}:
                    word = (1, 2, -1, -2, 1, -2)
                    examples.append({"strands": n, "protected": sorted(protected),
                                     "mobile": mobile, "mobile_word": word,
                                     "lifted_band_word": lift(word, mobile),
                                     "deleted_word": delete_protected(lift(word, mobile), n, protected)})
    return counts, examples


def cyclically_reduced(word):
    word = list(reduced(word))
    while len(word) > 1 and word[0] == -word[-1]:
        word = word[1:-1]
    return tuple(word)


def negative_controls():
    basis = [(1,), (2,), (3,)]
    # The opposite route between the same endpoint strands does not fix
    # the intermediate based meridian, despite the same strand permutation.
    assert hurwitz(basis, band(1, 3))[1] == (2,)
    assert hurwitz(basis, (-1, 2, 1))[1] != (2,)
    # A pure full twist fixes peripheral conjugacy classes, not free words.
    after = hurwitz([(1,), (2,)], (1, 1))
    assert after != [(1,), (2,)]
    assert [cyclically_reduced(w) for w in after] == [(1,), (2,)]
    return {"opposite_band_rejected": True,
            "fixed_conjugacy_is_not_fixed_word": True,
            "pure_full_twist_output": after}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    bands = band_formula_controls()
    counts, examples = subset_controls()
    negative = negative_controls()
    inputs = [HERE / "FREE_MERIDIAN_STABILIZERS.md", Path(__file__),
              HERE / "free_meridian_stabilizers_source_ledger.json",
              HERE / "sources/farb-margalit-primer-published.pdf",
              HERE / "sources/birman-ko-lee-band-generators-1998.pdf"]
    result = {"status": "passed", "python": platform.python_version(),
              "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in inputs},
              "exact_band_formula_controls": bands,
              **counts, "examples": examples, "negative_controls": negative,
              "scope": "Bounded exact Artin free-word formula, protected subsets, framed Hurwitz actions, braid relations, and strand-deletion left inverses. The unrestricted full stabilizer theorem depends on the written based-arc and cutting proof, not on this finite computation. No formal or external review is asserted."}
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {"sources", "examples", "negative_controls", "scope"}}, indent=2))


if __name__ == "__main__":
    main()
