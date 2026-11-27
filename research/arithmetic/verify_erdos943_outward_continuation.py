"""Exact rational optimization for the outward-continuation gcd cutoff."""

from itertools import product

from sage.all import MixedIntegerLinearProgram, QQ

CYCLE = range(6)
CLIP = QQ(33) / 2
PAIR_BOUND = 26
TRIPLE_BOUND = QQ(669) / 20
F_BOUND = QQ(43) / 2


def maximum_for_branch(forward_high, reverse_high, f_high):
    program = MixedIntegerLinearProgram(solver="PPL")
    f = program.new_variable(nonnegative=True)
    k = program.new_variable(nonnegative=True)
    smaller_matching_sum = program.new_variable(nonnegative=True)[0]

    def clipped(expression, high):
        if high:
            program.add_constraint(expression >= CLIP)
            return CLIP
        program.add_constraint(expression <= CLIP)
        return expression

    u, forward, reverse = [], [], []
    for i in CYCLE:
        program.add_constraint(f[i] <= F_BOUND)
        forward_size = f[i] + k[(i - 2) % 6]
        reverse_size = f[i] + k[(i - 1) % 6]
        program.add_constraint(forward_size <= 33)
        program.add_constraint(reverse_size <= 33)
        u.append(clipped(f[i], i in f_high))
        forward.append(clipped(forward_size, i in forward_high))
        reverse.append(clipped(reverse_size, i in reverse_high))

    total_f = sum(f[i] for i in CYCLE)
    for parity in (0, 1):
        matching_sum = sum(k[i] for i in range(parity, 6, 2))
        program.add_constraint(smaller_matching_sum <= matching_sum)
        program.add_constraint(total_f + 2 * matching_sum <= 198)

    for i in CYCLE:
        next_vertex, second_next = (i + 1) % 6, (i + 2) % 6
        previous, second_previous = (i - 1) % 6, (i - 2) % 6
        program.add_constraint(forward[i] + forward[next_vertex] <= PAIR_BOUND)
        program.add_constraint(reverse[i] + reverse[previous] <= PAIR_BOUND)
        program.add_constraint(reverse[i] + forward[second_next] <= PAIR_BOUND)
        program.add_constraint(
            u[i] + forward[next_vertex] + forward[second_next] <= TRIPLE_BOUND
        )
        program.add_constraint(
            u[i] + reverse[previous] + reverse[second_previous] <= TRIPLE_BOUND
        )

    program.set_objective(2 * (total_f + smaller_matching_sum))
    return program.solve()


def main():
    masks = [
        frozenset(i for i in CYCLE if bits & (1 << i))
        for bits in range(1 << 6)
    ]
    independent = [
        mask for mask in masks if all((i + 1) % 6 not in mask for i in mask)
    ]
    expected = {0: QQ(764) / 5, 1: QQ(719) / 5, 2: QQ(769) / 5}
    maxima = {size: QQ(0) for size in expected}
    examined = 0
    for forward, reverse in product(independent, repeat=2):
        if any((i + 2) % 6 in forward for i in reverse):
            continue
        for f_high in masks:
            if not f_high <= forward & reverse:
                continue
            value = maximum_for_branch(forward, reverse, f_high)
            assert value <= expected[len(f_high)], (forward, reverse, f_high, value)
            maxima[len(f_high)] = max(maxima[len(f_high)], value)
            examined += 1
    assert examined == 250
    assert maxima == expected
    print(f"{examined} rational branches. Maxima by number of high F positions: {maxima}")


if __name__ == "__main__":
    main()
