#!/usr/bin/env sage -python
"""Certify the three QQ halves, their MW names, and a saturated 13-basis."""

import hashlib
import itertools
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
HALVINGS = LOCAL / "q12o5867-three-target-halvings-qq.json"
CROSS = LOCAL / "q12o5867-two-primary-target-support-cross-prime.json"
QUOTIENT = LOCAL / "q12o5867-two-primary-cosets-mod89.json"
GATE = GENERATED / "elkies-k3-p1229-mw13-degree1-recovery-gate.json"
CLASSIFIER = LOCAL / "q12o5867-p0-shell-lattice-classification-all-profiles-mod89.json"
OUTPUT = GENERATED / "elkies-k3-q12o5867-two-primary-boundary.json"
Q3 = LOCAL / "q12o5867-degree1-compiler-branch-qq.json"
ABEL = LOCAL / "q12o5867-abel-trace-named-seeds-qq.json"
FRAME = GENERATED / "elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt"
SUPPORT_ARTIFACTS = (
    LOCAL / "q12o5867-saturation-completion-class21-qq.json",
    LOCAL / "q12o5867-replacement-word-seeds-qq.json",
    LOCAL / "q12o5867-two-primary-target-support-sections-qq.json",
    LOCAL / "q12o5867-two-primary-target-support-sections-v2-qq.json",
    LOCAL / "q12o5867-support-class170-shell32-qq.json",
)
PROFILE_SHELLS = {
    prime: LOCAL / f"q12o5867-p0-shell-all-profiles-mod{prime}.json"
    for prime in (83, 89, 137)
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
halvings = json.loads(HALVINGS.read_text())
cross = json.loads(CROSS.read_text())
quotient = json.loads(QUOTIENT.read_text())
gate = json.loads(GATE.read_text())
classifier = json.loads(CLASSIFIER.read_text())
profile_shells = {
    prime: json.loads(path.read_text()) for prime, path in PROFILE_SHELLS.items()
}
classes = classifier["lattice_shell"]["classes"]
height_gram = matrix(QQ, gate["parent_lattice"]["saturated_shioda_height_gram"])
assert height_gram.det() == QQ(237)/4


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


frame = load_matrix(FRAME)
ns_gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)

R = PolynomialRing(QQ, "u")
u = R.gen()
A = R(q8["child"]["minimal_A_coefficients_low_to_high"])
B = R(q8["child"]["minimal_B_coefficients_low_to_high"])
discriminant = 4*A**3+27*B**2


def polynomial_section(x_values, y_values):
    x = R(x_values)
    y = R(y_values)
    assert y**2 == x**3+A*x+B
    return x, y


def smooth_intersection(left, right):
    x_difference = left[0]-right[0]
    y_difference = left[1]-right[1]
    if not x_difference and not y_difference:
        return None
    if not x_difference:
        common = y_difference.monic()
    elif not y_difference:
        common = x_difference.monic()
    else:
        common = x_difference.gcd(y_difference).monic()
    if common.gcd(discriminant).degree() > 0:
        return None
    infinity_orders = []
    if x_difference:
        infinity_orders.append(4-int(x_difference.degree()))
    if y_difference:
        infinity_orders.append(6-int(y_difference.degree()))
    return int(common.degree()+max(0, min(infinity_orders)))


exact_sections = {}
for path in SUPPORT_ARTIFACTS:
    data = json.loads(path.read_text())
    for key, record in data["sections"].items():
        identity = (str(path.relative_to(ROOT)), key)
        section = record["section"]
        exact_sections[identity] = {
            "kind": "target_support",
            "point": polynomial_section(
                section["x_coefficients_low_to_high"],
                section["y_coefficients_low_to_high"],
            ),
        }
q3 = json.loads(Q3.read_text())
q3_identity = (str(Q3.relative_to(ROOT)), "exact_Q3")
exact_sections[q3_identity] = {
    "kind": "existing_exact_Q3",
    "point": polynomial_section(
        q3["section"]["x_coefficients_low_to_high"],
        q3["section"]["y_coefficients_low_to_high"],
    ),
}
abel = json.loads(ABEL.read_text())
q4_key = "Q4_candidate1_shell220_rank12_seed"
q4_identity = (str(ABEL.relative_to(ROOT)), q4_key)
q4_section = abel["sections"][q4_key]["section"]
exact_sections[q4_identity] = {
    "kind": "existing_exact_Q4",
    "point": polynomial_section(
        q4_section["x_coefficients_low_to_high"],
        q4_section["y_coefficients_low_to_high"],
    ),
}
for attempt in halvings["attempts"]:
    half = attempt["verified_rational_halves"][0]
    key = f"declared_class{attempt['declared_mod89_half_target_class_index']}"
    identity = (str(HALVINGS.relative_to(ROOT)), key)
    exact_sections[identity] = {
        "kind": "recovered_half",
        "point": polynomial_section(
            half["x"]["numerator_coefficients_low_to_high"],
            half["y"]["numerator_coefficients_low_to_high"],
        ),
    }

audit_by_identity = {
    (row["source_artifact"], row["source_key"]): row
    for row in cross["sections"]
}
assert set(exact_sections) == set(audit_by_identity)
for identity, record in exact_sections.items():
    audit = audit_by_identity[identity]
    profiles = {tuple(reduction["equation_component_profile"])
                for reduction in audit["reductions"]}
    assert len(profiles) == 1
    record["equation_component_profile"] = next(iter(profiles))
    record["supplied_class_index"] = int(audit["supplied_class_index"])


def class_vector(class_index):
    return vector(ZZ, classes[class_index]["current_4A1_mw"])


class_section_cache = {}
def class_section_vector(class_index):
    if class_index in class_section_cache:
        return class_section_cache[class_index]
    item = classes[class_index]
    mw = class_vector(class_index)
    component_pairings = vector(ZZ, item["current_component_pairings"])
    rhs = component_pairings-mw*frame[4:, :4]
    root_tail = frame[:4, :4].transpose().solve_right(rhs)
    assert all(value in ZZ for value in root_tail)
    tail = vector(ZZ, list(root_tail)+list(mw))
    assert tail*frame[:, :4] == component_pairings
    answer = vector(ZZ, [1, 1]+list(tail))
    class_section_cache[class_index] = answer
    return answer


def expected_intersection(left_index, right_index):
    return int(class_section_vector(left_index)*ns_gram*class_section_vector(right_index))


# The all-identity-profile labels are independent of any profile transport;
# Q3 and Q4 additionally have pinned exact Abel names.  Verify their literal
# intersections against the full NS vectors before using them as anchors.
anchor_identities = []
for identity, record in exact_sections.items():
    audit = audit_by_identity[identity]
    supplied = record["supplied_class_index"]
    all_identity_profile = record["equation_component_profile"] == (0, 0, 0, 0)
    unique_at_all_primes = all(
        reduction.get("profile_compatible_lattice_class_indices") == [supplied]
        for reduction in audit["reductions"]
    )
    if identity == q3_identity or (
        all_identity_profile and unique_at_all_primes
    ):
        anchor_identities.append(identity)
assert len(anchor_identities) >= 8
for left_identity, right_identity in itertools.combinations(anchor_identities, 2):
    observed = smooth_intersection(
        exact_sections[left_identity]["point"], exact_sections[right_identity]["point"]
    )
    if observed is None:
        continue
    expected = expected_intersection(
        exact_sections[left_identity]["supplied_class_index"],
        exact_sections[right_identity]["supplied_class_index"],
    )
    assert observed == expected, (
        left_identity, exact_sections[left_identity]["supplied_class_index"],
        right_identity, exact_sections[right_identity]["supplied_class_index"],
        observed, expected,
    )

for identity, record in exact_sections.items():
    comparisons = []
    audit = audit_by_identity[identity]
    good_degree = next(
        reduction["inverse_parent_degree"] for reduction in audit["reductions"]
        if reduction["prime"] == 137
    )
    candidates = [
        int(item["class_index"]) for item in classes
        if int(item["q4o164_parent_degree"]) == good_degree
    ]
    for anchor_identity in anchor_identities:
        if anchor_identity == identity:
            continue
        observed = smooth_intersection(
            record["point"], exact_sections[anchor_identity]["point"]
        )
        if observed is None:
            continue
        anchor_index = exact_sections[anchor_identity]["supplied_class_index"]
        candidates = [
            class_index for class_index in candidates
            if expected_intersection(class_index, anchor_index) == observed
        ]
        comparisons.append({
            "anchor_source_artifact": anchor_identity[0],
            "anchor_source_key": anchor_identity[1],
            "anchor_class_index": anchor_index,
            "observed_smooth_intersection": observed,
            "surviving_class_indices": candidates,
        })
        if len(candidates) == 1:
            break
    assert len(candidates) == 1, (identity, candidates)
    record["certified_class_index"] = candidates[0]
    record["intersection_fingerprint"] = comparisons
    record["mw_vector"] = class_vector(candidates[0])
    record["exact_shioda_height"] = (
        record["mw_vector"]*height_gram*record["mw_vector"]
    )

# Verify all three exact group-law identities again at the MW-vector level,
# now using the intersection-certified names of every support section.
identity_by_supplied_class = {}
for identity, record in exact_sections.items():
    if record["kind"] == "target_support":
        identity_by_supplied_class[record["supplied_class_index"]] = identity

halving_verifications = []
for attempt in halvings["attempts"]:
    declared = int(attempt["declared_mod89_half_target_class_index"])
    half_identity = (str(HALVINGS.relative_to(ROOT)), f"declared_class{declared}")
    half_record = exact_sections[half_identity]
    target_vector = vector(ZZ, 13)
    support_terms = []
    for term in attempt["exact_QQ_target_support_word"]:
        supplied = int(term["supplied_physical_class_index"])
        support_identity = identity_by_supplied_class[supplied]
        support_record = exact_sections[support_identity]
        coefficient = ZZ(term["coefficient"])
        target_vector += coefficient*support_record["mw_vector"]
        support_terms.append({
            "coefficient": int(coefficient),
            "supplied_class_index": supplied,
            "certified_class_index": support_record["certified_class_index"],
        })
    assert 2*half_record["mw_vector"] == target_vector
    halving_verifications.append({
        "supplied_half_class_index": declared,
        "certified_half_class_index": half_record["certified_class_index"],
        "quotient_key": attempt["declared_quotient_key"],
        "equation_component_profile": list(half_record["equation_component_profile"]),
        "exact_shioda_height": str(half_record["exact_shioda_height"]),
        "support_word_with_certified_names": support_terms,
        "literal_curve_substitution": True,
        "literal_QQ_doubling": True,
        "exact_MW_vector_doubling": True,
    })

current_hnf = matrix(ZZ, quotient["current_subgroup"]["row_hnf"])
half_vectors = matrix(ZZ, [
    exact_sections[(str(HALVINGS.relative_to(ROOT)), f"declared_class{row['supplied_half_class_index']}")][
        "mw_vector"
    ]
    for row in halving_verifications
])
recovered_module = current_hnf.stack(half_vectors).row_module(ZZ)
assert abs(recovered_module.basis_matrix().det()) == 1
assert recovered_module.basis_matrix().elementary_divisors() == [1]*13


def current_quotient_key(mw_vector):
    remainder = vector(ZZ, mw_vector)
    for row_index in range(10):
        assert current_hnf[row_index, row_index] == 1
        remainder -= remainder[row_index]*current_hnf[row_index]
    assert all(remainder[index] == 0 for index in range(10))
    return [int(remainder[index] % 2) for index in range(10, 13)]


for row in halving_verifications:
    identity = (
        str(HALVINGS.relative_to(ROOT)),
        f"declared_class{row['supplied_half_class_index']}",
    )
    actual_key = current_quotient_key(exact_sections[identity]["mw_vector"])
    assert actual_key == row["quotient_key"]
    row["certified_quotient_key"] = actual_key

section_identities = list(exact_sections)
exact_section_matrix = matrix(ZZ, [
    exact_sections[identity]["mw_vector"] for identity in section_identities
])
exact_section_module = exact_section_matrix.row_module(ZZ)
print("name corrections", [
    (identity[1], exact_sections[identity]["supplied_class_index"],
     exact_sections[identity]["certified_class_index"])
    for identity in section_identities
    if exact_sections[identity]["supplied_class_index"]
       != exact_sections[identity]["certified_class_index"]
])
print("exact section span rank/index", exact_section_matrix.rank(),
      abs(exact_section_module.basis_matrix().det()))
print("current classes outside exact span", [
    class_index for class_index in quotient["current_subgroup"]["physical_class_indices"]
    if class_vector(class_index) not in exact_section_module
][:10])
assert exact_section_matrix.rank() == 13
assert abs(exact_section_module.basis_matrix().det()) == 2
exact_pool_height_gram = (
    exact_section_module.basis_matrix()*height_gram
    *exact_section_module.basis_matrix().transpose()
)
assert exact_pool_height_gram.det() == 237
outside_exact_span = [
    class_index for class_index in quotient["current_subgroup"]["physical_class_indices"]
    if class_vector(class_index) not in exact_section_module
]
assert outside_exact_span == [21, 48, 216, 516, 543, 701]

payload = {
    "schema": "q12o5867-two-primary-boundary-v1",
    "status": "CLOSED_EXACT_TWO_PRIMARY_BOUNDARY_EQUATION_POOL_REMAINS_INDEX_2",
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (
            (Q8, HALVINGS, CROSS, QUOTIENT, GATE, CLASSIFIER, Q3, ABEL, FRAME)
            +SUPPORT_ARTIFACTS+tuple(PROFILE_SHELLS.values())
        )
    },
    "complete_component_profile_pass": {
        str(prime): {
            "signed_section_count": data["direct_shell"]["unique_signed_section_count"],
            "component_profile_scope": data["direct_shell"]["component_profile_scope"],
            "component_profile_histogram": data["direct_shell"]["component_profile_histogram"],
            "all_sixteen_profiles_nonempty": (
                len(data["direct_shell"]["component_profile_histogram"]) == 16
                and min(data["direct_shell"]["component_profile_histogram"].values()) > 0
            ),
        }
        for prime, data in profile_shells.items()
    },
    "three_halvings": halving_verifications,
    "smith_index_reduction": {
        "before_smith_diagonal": quotient["current_subgroup"]["smith_diagonal"],
        "before_index": quotient["current_subgroup"]["index"],
        "after_smith_diagonal": [1]*13,
        "after_index": 1,
    },
    "equation_level_basis_gate": {
        "exact_section_count_audited": len(exact_sections),
        "rank": 13,
        "index_in_saturated_MW13": 2,
        "smith_diagonal": list(map(
            int, exact_section_module.basis_matrix().elementary_divisors()
        )),
        "height_gram_determinant_of_exact_pool_lattice": "237",
        "saturated_height_gram_determinant": "237/4",
        "saturated_thirteen_section_basis_recovered": False,
        "current_class_indices_in_missing_exact_pool_coset": outside_exact_span,
        "bounded_completion_lifts": {
            "nonlifting_mod89_branches": [216, 701, 543, 48],
            "lifted_but_intersection_relabelled": {"21": 26, "516": 521},
        },
    },
    "all_exact_sections_name_audit": [
        {
            "source_artifact": identity[0],
            "source_key": identity[1],
            "kind": record["kind"],
            "supplied_class_index": record["supplied_class_index"],
            "certified_class_index": record["certified_class_index"],
            "equation_component_profile": list(record["equation_component_profile"]),
            "exact_shioda_height": str(record["exact_shioda_height"]),
            "intersection_fingerprint": record["intersection_fingerprint"],
        }
        for identity, record in exact_sections.items()
    ],
    "proof_boundary": (
        "One complete sixteen-profile pass at 83/89/137 and exactly three targeted "
        "QQ(u) halvings were performed. No unrestricted inverse-degree shell, fourth "
        "halving, parent-fibre search, or ratpoints-bound increase is authorized."
    ),
    "point_factory_decision": {
        "decision": "CLOSE_Q12_ORBIT5867_AS_RANK32_POINT_PRODUCTION_ROUTE",
        "reason": (
            "The abstract quotient is killed, but the bounded exact equation pool remains "
            "index two and therefore cannot support certified MW13 words for the controls."
        ),
        "retain_exact_birational_map_theorem": True,
        "public_control_MW13_words_computed": False,
        "bounded_parent_MW_enumeration_authorized": False,
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(f"wrote {OUTPUT}")
print("halves", [(row["certified_half_class_index"], row["quotient_key"])
                 for row in halving_verifications])
print("exact equation pool index", abs(exact_section_module.basis_matrix().det()))
print("decision", payload["point_factory_decision"]["decision"])
