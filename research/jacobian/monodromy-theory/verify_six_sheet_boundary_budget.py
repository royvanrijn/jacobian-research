"""Check the finite monodromy profiles after the geometric support bound."""
from itertools import combinations_with_replacement

from sage.all import libgap


def profile_census(sheet_count, branch_support_cap):
    symmetric_group = libgap.SymmetricGroup(sheet_count)
    transitive_count = 0
    survivors = set()
    subgroup_classes = libgap.ConjugacyClassesSubgroups(symmetric_group)
    for subgroup_class in subgroup_classes:
        group = libgap.Representative(subgroup_class)
        if not libgap.IsTransitive(group, list(range(1, sheet_count + 1))):
            continue
        transitive_count += 1
        candidates = []
        for conjugacy_class in libgap.ConjugacyClasses(group):
            representative = libgap.Representative(conjugacy_class)
            support = int(libgap.Size(libgap.MovedPoints(representative)))
            if 2 <= support <= branch_support_cap:
                normal = libgap.NormalClosure(group, libgap.Group([representative]))
                candidates.append((support, normal))
        for number in range(1, (sheet_count - 1) // 2 + 1):
            for branches in combinations_with_replacement(candidates, number):
                profile = tuple(sorted(support for support, _ in branches))
                if sum(profile) > sheet_count - 1:
                    continue
                normal_generators = [
                    generator
                    for _, normal in branches
                    for generator in libgap.GeneratorsOfGroup(normal)
                ]
                generated = libgap.Group(normal_generators)
                if generated == group:
                    survivors.add((int(libgap.Size(group)), profile))
    return transitive_count, survivors


def main():
    assert profile_census(5, 2) == (5, {(120, (2,)), (120, (2, 2))})
    assert profile_census(6, 3) == (
        16, {(720, (2,)), (720, (2, 2)), (720, (2, 3)), (360, (3,))},
    )
    print("Five sheets: S5 with support (2,) or (2,2).")
    print("Six sheets: S6 with support (2,), (2,2), or (2,3); A6 with (3,).")


if __name__ == "__main__":
    main()
