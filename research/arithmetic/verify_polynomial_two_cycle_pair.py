#!/usr/bin/env python3
"""Exact forest normalization and pure-braid components for two cycle pairs."""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import factorial, lcm
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_duplicated_cycle_paths as x

if not __debug__:
    raise RuntimeError("Two-cycle-pair checker requires active assertions")
HERE = Path(__file__).resolve().parent
c, m, l, h, u = x.c, x.m, x.l, x.h, x.u


def forest(state):
    return x.forest(state, repeated=False)


def path_data(state):
    d = forest(state)
    assert d["kind"] == "path"
    al, ai = (d["da"].index(k) for k in (1, 2))
    bl, bi = (d["db"].index(k) for k in (1, 2))
    p, = d["ac"][al] & d["bc"][bi]
    q, = d["ac"][ai] & d["bc"][bi]
    r, = d["ac"][ai] & d["bc"][bl]
    ai_len, al_len = len(d["ac"][ai]), len(d["ac"][al])
    bi_len, bl_len = len(d["bc"][bi]), len(d["bc"][bl])
    alpha, beta = x.distance(state[0], q, r), x.distance(state[1], p, q)
    assert 1 <= alpha < ai_len and 1 <= beta < bi_len
    return dict(d, al=al, ai=ai, bl=bl, bi=bi, p=p, q=q, r=r,
                ai_len=ai_len, al_len=al_len, bi_len=bi_len, bl_len=bl_len,
                alpha=alpha, beta=beta)


def path_seed(ai_len, al_len, bi_len, bl_len, alpha=1, beta=1):
    assert min(ai_len, al_len, bi_len, bl_len) >= 2
    assert 1 <= alpha < ai_len and 1 <= beta < bi_len
    ai = tuple(range(ai_len)); nxt = ai_len
    bi = (0,)+tuple(range(nxt, nxt+bi_len-1)); nxt += bi_len-1
    p, r = bi[bi_len-beta], ai[alpha]
    al = (p,)+tuple(range(nxt, nxt+al_len-1)); nxt += al_len-1
    bl = (r,)+tuple(range(nxt, nxt+bl_len-1)); nxt += bl_len-1
    N = nxt+1
    state = (m.perm(N, [ai, al]), m.perm(N, [bi, bl]), m.perm(N, [(nxt, 0)]))
    d = path_data(state)
    assert (d["ai_len"], d["al_len"], d["bi_len"], d["bl_len"],
            d["alpha"], d["beta"]) == (ai_len, al_len, bi_len, bl_len, alpha, beta)
    return state


def star_to_path(state):
    d = forest(state)
    if d["kind"] == "B-star":
        changed = star_to_path(m.braid(state, (1,)))
        result = m.braid(changed, (-1,))
        assert forest(result)["kind"] == "path"
        return result
    assert d["kind"] == "A-star"
    isolated = d["ac"][d["da"].index(0)]
    ends, = x.supports(state[2]); start, = ends-isolated
    Bonly = set.union(*d["bc"])-set.union(*d["ac"])
    result = x.capture(x.rotate_to(state, start, min(Bonly)))
    assert forest(result)["kind"] == "path"
    return result


def reattach(state, other_b, target):
    """Cut the A-leaf edge and reattach to either B cycle at a chosen distance."""
    old = path_data(state)
    length = old["bl_len"] if other_b else old["bi_len"]
    assert 1 <= target < length
    ends, = x.supports(state[2]); start, = ends-old["outside"]
    cut = x.capture(x.rotate_to(state, start, old["p"]))
    d = forest(cut)
    assert d["kind"] == "A-star" and cut[0] == state[0]
    isolated = old["ac"][old["al"]]
    ends, = x.supports(cut[2]); start, = ends-isolated
    retained = old["r"] if other_b else old["q"]
    y = x.power(m.inverse(cut[1]), target)[retained]
    assert y not in set.union(*d["ac"])
    result = x.capture(x.rotate_to(cut, start, y))
    new = path_data(result)
    assert result[0] == state[0]
    assert (new["ai_len"], new["al_len"]) == (old["ai_len"], old["al_len"])
    expected_alpha = old["ai_len"]-old["alpha"] if other_b else old["alpha"]
    assert (new["alpha"], new["beta"]) == (expected_alpha, target)
    assert (new["bi_len"], new["bl_len"]) == (
        (old["bl_len"], old["bi_len"]) if other_b else (old["bi_len"], old["bl_len"]))
    if not other_b:
        leaf = old["bc"][old["bl"]]
        assert all(result[1][z] == state[1][z] for z in leaf)
    return result


def change_alpha(state, target):
    old = path_data(state)
    assert 1 <= target < old["ai_len"]
    swapped = m.braid(state, (1,))
    d = path_data(swapped)
    assert (d["alpha"], d["beta"]) == (
        old["bi_len"]-old["beta"], old["ai_len"]-old["alpha"])
    changed = reattach(swapped, False, old["ai_len"]-target)
    result = m.braid(changed, (-1,))
    new = path_data(result)
    assert (new["alpha"], new["beta"]) == (target, old["beta"])
    return result


def switch_a(state):
    old = path_data(state)
    changed = reattach(m.braid(state, (1,)), True, 1)
    result = m.braid(changed, (-1,))
    new = path_data(result)
    assert result[1] == state[1]
    assert (new["ai_len"], new["al_len"], new["bi_len"], new["bl_len"]) == (
        old["al_len"], old["ai_len"], old["bi_len"], old["bl_len"])
    return result


def normalize_path(state):
    d = path_data(state)
    if d["bi_len"] > d["bl_len"]:
        state = reattach(state, True, 1)
    d = path_data(state)
    if d["ai_len"] > d["al_len"]:
        state = switch_a(state)
    state = change_alpha(reattach(state, False, 1), 1)
    d = path_data(state)
    ends, = x.supports(state[2]); start, = ends-d["outside"]
    state = x.rotate_to(state, start, d["q"])
    expected = path_seed(min(d["lengths_a"]), max(d["lengths_a"]),
                         min(d["lengths_b"]), max(d["lengths_b"]))
    absolute = lambda z: c.cyclic_label(c.v.one.fixed_product(z))
    assert absolute(state) == absolute(expected)
    return state


def matching_seed(a1, a2, b1, b2, phase):
    assert min(a1, a2, b1, b2) >= 2
    lengths = (a1+b1-1, a2+b2-1)
    assert 0 <= phase < lengths[1]
    ac, bc = [], []; root = 0
    for a, b in ((a1, b1), (a2, b2)):
        ac.append(tuple(range(root, root+a)))
        bc.append((root,)+tuple(range(root+a, root+a+b-1)))
        root += a+b-1
    A, B = m.perm(root, ac), m.perm(root, bc)
    y = x.power(m.compose(A, B), phase)[lengths[0]]
    state = A, B, m.perm(root, [(0, y)])
    assert forest(state)["kind"] == "matching"
    return state


def matching_to_path(state):
    d = forest(state)
    assert d["kind"] == "matching"
    ua, ub = set.union(*d["ac"]), set.union(*d["bc"])
    kind = lambda z: "O" if z in ua & ub else "A" if z in ua else "B"
    period = lcm(*map(len, d["components"]))
    for k in range(period):
        # Locate a differing support type without repeatedly expanding long twists.
        ends, = x.supports(state[2])
        rk = x.power(d["R"], k)
        pair = {kind(rk[z]) for z in ends}
        if len(pair) == 1:
            continue
        trial = x.rotate(state, k)
        if pair == {"O", "B"}:
            # Q1^-1 leaves A in the second heavy slot and makes the endpoints exclusive.
            swapped = m.braid(trial, (-1,))
            result = m.braid(x.capture(swapped), (1,))
        else:
            result = x.capture(trial)
            if pair == {"O", "A"}:
                assert forest(result)["kind"] == "A-star"
                result = star_to_path(result)
        assert forest(result)["kind"] == "path"
        return result, "-".join(sorted(pair))
    # Equality of the periodic three-letter words forces equal lengths and phase.
    assert len(set(d["lengths_a"])) == len(set(d["lengths_b"])) == 1
    return None, "aligned"


def pure_words():
    words = []
    for i in range(3):
        for end in range(i+1, 3):
            path = tuple(range(end, i+1, -1))
            word = path+(i+1, i+1)+tuple(-z for z in path[::-1])
            words.extend((word, tuple(-z for z in word[::-1])))
    for word in words+[(1, 1), (2, 2), (-1, 2, 2, 1), (1, 2, 2, -1)]:
        strands = list(range(3))
        for z in word:
            i = abs(z)-1; strands[i], strands[i+1] = strands[i+1], strands[i]
        assert strands == [0, 1, 2]
    return words


def pure_closure(seed):
    words = pure_words()
    start = c.cyclic_label(seed); seen = {start}; queue = deque([start])
    pattern = tuple(tuple(m.cycle_lengths(g)) for g in seed)
    while queue:
        state = queue.popleft()
        for word in words:
            changed = c.cyclic_label(m.braid(state, word))
            assert tuple(tuple(m.cycle_lengths(g)) for g in changed) == pattern
            if changed not in seen:
                seen.add(changed); queue.append(changed)
    return seen


def normalization_controls():
    domains = list(product(range(2, 5), repeat=4))+[
        (2, 11, 3, 17), (7, 2, 11, 3), (5, 13, 2, 7), (9, 4, 8, 3)]
    rows = []
    for ai, al, bi, bl in domains:
        modes = Counter(); aligned = 0
        for phase in range(al+bl-1):
            state = matching_seed(ai, al, bi, bl, phase)
            changed, mode = matching_to_path(state)
            if changed is None:
                assert ai == al and bi == bl and phase == 0
                aligned += 1
            else:
                normalize_path(changed); modes[mode] += 1
        small = max(ai, al, bi, bl) <= 4
        distances = list(product(range(1, ai), range(1, bi))) if small else sorted(
            set(product((1, ai-1), (1, bi-1))))
        beta_controls = alpha_controls = switches = star_controls = 0
        for alpha, beta in distances:
            state = path_seed(ai, al, bi, bl, alpha, beta)
            for target in (range(1, bi) if small else (1, bi-1)):
                reattach(state, False, target); beta_controls += 1
            for target in (range(1, ai) if small else (1, ai-1)):
                change_alpha(state, target); alpha_controls += 1
            reattach(state, True, 1); switch_a(state); switches += 2
            normalize_path(state)
        state = path_seed(ai, al, bi, bl)
        d = path_data(state)
        star = x.capture(x.rotate_to(state, d["q"], d["p"]))
        sd = forest(star)
        isolated = sd["ac"][sd["da"].index(0)]
        main = next(z for z in sd["components"] if z != isolated)
        ends, = x.supports(star[2]); start, = ends-isolated
        for target in sorted(main):
            trial = x.rotate_to(star, start, target)
            normalize_path(star_to_path(trial))
            normalize_path(star_to_path(m.braid(trial, (1,))))
            star_controls += 2
        rows.append(dict(lengths=[ai, al, bi, bl], degree=ai+al+bi+bl-2,
                         matching_escape_modes=dict(modes), aligned_matchings=aligned,
                         path_distance_pairs=len(distances), beta_controls=beta_controls,
                         alpha_controls=alpha_controls, role_switches=switches,
                         star_controls=star_controls))
    return rows


def primitive_group(seed):
    N = len(seed[0]); base, blocks, edges = l.group_from_transposition(seed)
    assert (base, blocks, edges) == (N, 1, N*(N-1)//2)
    assert c.block_divisors(seed) == [1, N] and c.stabilizer_size(seed) == 1
    return dict(actual_group_order=factorial(N), normal_transposition_support_edges=edges,
                all_block_count_divisors=[1, N], deck=1)


def complete_components():
    profiles = list(combinations_with_replacement(range(2, 5), 2))
    domains = list(combinations_with_replacement(profiles, 2))+[
        ((2, 5), (3, 7)), ((2, 6), (4, 5)), ((3, 8), (3, 8)), ((4, 9), (4, 9))]
    rows = []
    for aa, bb in domains:
        N, _, _, passport = u.first_deficit((aa, bb), 1)
        denominator = (1+(aa[0] == aa[1]))*(1+(bb[0] == bb[1]))
        total = Fraction(N*N*(sum(aa)-1)*(sum(bb)-1), denominator)
        assert total.denominator == 1
        total = int(total)
        assert total == c.character_total(passport) == c.comparison_total(passport)
        placements = 3*(2-(aa == bb))
        assert c.type_placements(passport) == placements
        repeated = aa[0] == aa[1] and bb[0] == bb[1]
        primitive_seed = c.v.one.fixed_product(path_seed(*aa, *bb))
        orbit = pure_closure(primitive_seed)
        primitive_literal = total-((N//2)**2 if repeated else 0)
        assert len(orbit)*N == primitive_literal
        components = [dict(primitive_group(primitive_seed), component="primitive",
                           literal_seed=primitive_seed, fixed_order_pure_absolute_size=len(orbit),
                           fixed_order_pure_literal_size=primitive_literal,
                           full_absolute_size=placements*len(orbit),
                           full_literal_size=placements*primitive_literal,
                           orbit_sha256=c.digest(orbit))]
        if repeated:
            imp_seed = c.v.one.fixed_product(matching_seed(*aa, *bb, 0))
            imp_orbit = pure_closure(imp_seed)
            control = x.group_control(imp_seed, False, aa[0], bb[0])
            assert len(imp_orbit) == N//2 and not orbit & imp_orbit
            assert len(orbit)*N+len(imp_orbit)*(N//2) == total
            components.append(dict(control, component="quadratic",
                                   literal_seed=imp_seed, fixed_order_pure_absolute_size=len(imp_orbit),
                                   fixed_order_pure_literal_size=(N//2)**2,
                                   full_absolute_size=placements*len(imp_orbit),
                                   full_literal_size=placements*(N//2)**2,
                                   orbit_sha256=c.digest(imp_orbit)))
        rows.append(dict(profiles=[aa, bb], degree=N, ordered_hook_total=total,
                         full_type_placements=placements, full_literal_domain=placements*total,
                         components=components))
        print(json.dumps(dict(profiles=[aa, bb], degree=N,
                              pure_absolute_sizes=[z["fixed_order_pure_absolute_size"] for z in components])),
              flush=True)
    return rows


def three_plus_one_controls():
    domains = [((2, 2, 2), 3), ((2, 2, 3), 3), ((2, 3, 3), 4),
               ((2, 3, 4), 3), ((3, 3, 3), 5), ((2, 4, 6), 7)]
    rows = []
    for heavy, q in domains:
        N, _, _, passport = u.first_deficit((heavy, (q,)), 1)
        seed = h.tree_seed((heavy, (q,)), 1)
        orbit = pure_closure(seed)
        total = c.character_total(passport)
        assert total == c.comparison_total(passport) == N*len(orbit)
        multiplicity = 1
        for count in Counter(heavy).values():
            multiplicity *= factorial(count)
        assert total == Fraction(N*N*q*(q-1), multiplicity)
        group = primitive_group(seed)
        placements = c.type_placements(passport)
        assert placements == 6
        rows.append(dict(heavy=heavy, pure_cycle=q, degree=N, group=group,
                         literal_seed=seed, fixed_order_pure_absolute_size=len(orbit),
                         ordered_hook_total=total, full_literal_domain=placements*total,
                         orbit_sha256=c.digest(orbit)))
    return rows


def negative_controls():
    rejected = 0
    for fn, args in ((path_seed, (1, 2, 3, 4)), (path_seed, (3, 4, 2, 5, 3, 1)),
                     (matching_seed, (3, 4, 2, 3, 6))):
        try:
            fn(*args)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError("invalid domain accepted")
    # An unequal aligned-phase matching must not be mistaken for a quadratic lift.
    state, mode = matching_to_path(matching_seed(3, 2, 2, 3, 0))
    assert state is not None and mode != "aligned"
    # This has no immediate opposite A/B capture but does escape via a shared sheet.
    witness = matching_seed(3, 2, 2, 3, 1)
    d = forest(witness)
    ua, ub = set.union(*d["ac"]), set.union(*d["bc"])
    ends, = x.supports(witness[2])
    for k in range(lcm(*map(len, d["components"]))):
        rk = x.power(d["R"], k)
        changed = {rk[z] for z in ends}
        assert not (len(changed & (ua-ub)) == len(changed & (ub-ua)) == 1)
    state, mode = matching_to_path(witness)
    assert state is not None and "O" in mode
    return dict(rejected=rejected, unequal_phase_zero_not_imprimitive=True,
                shared_endpoint_escape_required=True, no_opposite_color_capture=True,
                shared_escape_literal_witness=witness)


def main():
    inputs = [
        "POLYNOMIAL_TWO_CYCLE_PAIR_COMPONENTS.md",
        "POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md",
        "POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md",
        "POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md",
        "POLYNOMIAL_PASSPORT_COMPONENTS.md", "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
        "verify_polynomial_duplicated_cycle_paths.py",
        "verify_polynomial_first_deficit.py", "verify_polynomial_stable_hermite.py",
        "verify_polynomial_phase_interpolation.py", "verify_polynomial_involution_components.py",
        "verify_polynomial_lift_transport.py", "verify_polynomial_passport_components.py",
        "sources/jones-zvonkin-braid-cacti-2002.pdf",
        "sources/goulden-jackson-cacti-1992.pdf",
        "sources/wajnryb-two-special-fibers-1995.pdf",
        "sources/chen-gendron-prado-tahar-isoresidual-curves-2026.pdf"]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    controls = normalization_controls()
    components = complete_components()
    tail = three_plus_one_controls()
    result = dict(
        schema_version=1, kind="two_cycle_pair_pure_braid_components", result="PASS",
        python=platform.python_version(), checker_sha256=sha(Path(__file__).resolve()),
        input_sha256={path: sha(HERE/path) for path in dict.fromkeys(inputs)},
        normalization_controls=controls, complete_component_controls=components,
        three_plus_one_controls=tail,
        complete_passport_count=len(components)+len(tail),
        complete_literal_domains=sum(z["full_literal_domain"] for z in components+tail),
        negative_controls=negative_controls(), pure_generators=pure_words(),
        scope="Written Y1-Y5 classify every four-cycle two-heavy-factor plus one-T passport "
              "under both full and fixed-order pure braids. Exactly one primitive component, "
              "plus one quadratic component iff both heavy factors duplicate a cycle length. "
              "The three-plus-one case follows from the cited single-cycle-tail theorem, closing "
              "the complete first-deficit k2,s1 sector. Finite domains control the uniform proof, "
              "not an extrapolation. Higher-p, k>=3, deeper and zero-simple profiles remain open.",
        independent_replay=False, formal_verification=False, external_review=False)
    (HERE/"polynomial_two_cycle_pair_verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print("PASS:"+str(len(controls))+" normalization domains;"+
          str(result["complete_passport_count"])+" complete pure-braid passports;"+
          str(result["complete_literal_domains"])+" literal tuples represented")


if __name__ == "__main__":
    main()
