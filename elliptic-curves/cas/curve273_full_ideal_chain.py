"""Exact data for the first full-ideal descent chain on ICARM curve 273.

The data are kept separate from both the verifier and the relation-pool
analyser so those two independent consumers reconstruct the same algebraic
integers rather than maintaining parallel copies of the certificate.
"""

from __future__ import annotations


T0 = (
    (28691731813798755604363789, 17957201189903465826327159),
    (
        7159638381133483906634203654283170391,
        12780381281373253031851035100853459,
    ),
)

T1 = (
    (8108401645961, 3241943091229),
    (13757921887007, 10290585023712),
    (1998874339580503775477, 182072033289697848802),
)

T2 = (
    (1998874339580503775477, 182072033289697848802),
    (
        210549181078644643738293358684427,
        12055357087325822035141506894195,
    ),
)

T3 = (
    (6757889, 3837804),
    (9828720251573, 3497176632472),
    (48922801174561, 47639959748846),
    (50777509904197, 8431789759716),
)

T4 = (
    (50777509904197, 8431789759716),
    (
        28667248277432793199773894520960218073,
        12725678201649134529227295218218488637,
    ),
)

T5 = (
    (110693621, 3604894),
    (118153954113061, 26623425927990),
    (222635866003829, 168812229950784),
    (241667315076287, 75987985552523),
)


T6 = (
    (
        110082763798456567967312309869978698968187793,
        12357613814041997795032595200900142542684867,
    ),
)


T7 = (
    (3069949, 1911346),
    (63217513, 22678367),
    (168510611, 8392763),
    (7754996948873, 5321926736053),
    (38169724803043, 26317327110158),
)


T8 = (
    (19767179, 6801961),
    (2024658179, 1400637434),
    (12736098241512293, 9617972020917348),
    (231573912732876911, 136784656449821585),
)


T9 = (
    (7636709, 6045380),
    (56740549, 8493651),
    (104379329, 22584173),
    (413582930291, 366762253550),
    (555442705507, 488627836116),
)


PAIR_49_74 = (
    (296876930090401, 101399158692471),
    (9466946362153418565457, 5491186297291468700542),
)


PAIR_25_96 = (
    (29455103, 17783085),
    (
        73146185081283799302921629257,
        25214402471548825187489833593,
    ),
)


SUPPORTS = (T0, T1, T2, T3, T4, T5, T6, T7, T8, T9)


RELATION_SPECS = (
    {
        "name": "I1",
        "kind": "ideal",
        "before": T0,
        "after": T1,
        "declared": T0 + T1,
        "source": T0,
        "coordinates": (
            40452015641036084505034,
            -633786344511290342857518,
            188325461970930895,
        ),
    },
    {
        "name": "CRT2",
        "kind": "linear",
        "before": T1,
        "after": T2,
        "declared": T1[:2] + (T2[1],),
        "m": -619592186148608109179738569159086,
    },
    {
        "name": "I3",
        "kind": "ideal",
        "before": T2,
        "after": T3,
        "declared": T2 + T3,
        "source": T2,
        "coordinates": (
            8738173789610100,
            -22339703817949299,
            6638096711,
        ),
    },
    {
        "name": "CRT4",
        "kind": "linear",
        "before": T3,
        "after": T4,
        "declared": T3[:3] + (T4[1],),
        "m": -273946804572678797468511649991383692093,
    },
    {
        "name": "I5",
        "kind": "ideal",
        "before": T4,
        "after": T5,
        "declared": T4 + T5,
        "source": T4,
        "coordinates": (
            23936444895462507,
            -25909320534488312,
            -201409657931,
        ),
    },
    {
        "name": "I6",
        "kind": "ideal",
        "before": T5,
        "after": T6,
        "declared": T5 + T6,
        "source": T5,
        "coordinates": (
            1876028045451,
            -1977200499447,
            -8769733,
        ),
    },
    {
        "name": "I7",
        "kind": "ideal",
        "before": T6,
        "after": T7,
        "declared": T6 + T7,
        "source": T6,
        "coordinates": (
            4081147165751,
            -41375042746495,
            12294323,
        ),
    },
    {
        "name": "I8",
        "kind": "ideal",
        "before": T7,
        "after": T8,
        "declared": T7 + T8,
        "source": T7,
        "coordinates": (
            291950561784079,
            -6928772233568803,
            2058839297,
        ),
    },
    {
        "name": "I9",
        "kind": "ideal",
        "before": T8,
        "after": T9,
        "declared": T8 + T9,
        "source": T8,
        "coordinates": (
            10612874054972376,
            -12569391345381019,
            3734912328,
        ),
    },
    {
        "name": "I7_pair_49_74",
        "kind": "ideal",
        "before": T6,
        "after": PAIR_49_74,
        "declared": T6 + PAIR_49_74,
        "source": T6,
        "coordinates": (
            17624888393018,
            -178680918274026,
            -275292813,
        ),
    },
    {
        "name": "I8_pair_25_96",
        "kind": "ideal",
        "before": PAIR_49_74,
        "after": PAIR_25_96,
        "declared": PAIR_49_74 + PAIR_25_96,
        "source": PAIR_49_74,
        "coordinates": (
            5380042323,
            -7292970586,
            2051100,
        ),
    },
)


def prime_ideal(field, theta, label):
    """Return the degree-one prime ideal encoded by ``(q, theta mod q)``."""

    q, residue = label
    return field.ideal(q, theta - residue)


def ideal_element(field, theta, target, coordinates):
    """Reconstruct an integral element from coordinates in an ideal basis."""

    ideal = field.ideal(1)

    for label in target:
        ideal *= prime_ideal(field, theta, label)

    basis = tuple(ideal.basis())

    if len(basis) != 3:
        raise RuntimeError(f"unexpected ideal-basis length {len(basis)}")

    return sum(
        (
            coordinates[index] * basis[index]
            for index in range(3)
        ),
        field(0),
    )


def build_relations(field, theta):
    """Reconstruct all certified principal-ideal relations."""

    relations = []

    for spec in RELATION_SPECS:
        if spec["kind"] == "linear":
            alpha = field(spec["m"]) - theta
        else:
            alpha = ideal_element(
                field,
                theta,
                spec["source"],
                spec["coordinates"],
            )

        relation = dict(spec)
        relation["alpha"] = alpha
        relations.append(relation)

    return tuple(relations)
