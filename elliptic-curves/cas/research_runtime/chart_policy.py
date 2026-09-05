"""Independent search representation stages and calibrated metric policies.

A policy is scheduling information. A score, ordering, calibration win or
bounded miss never supplies a rank bound or a mathematical pruning rule.
"""
from dataclasses import asdict, dataclass
from fractions import Fraction
from time import monotonic

from .cvp import VoronoiIterator, ldl
from .store import digest


@dataclass(frozen=True)
class ChartPolicy:
    model_normalization: str = "raw"
    chart_parameterization: str = "pointed-quartic"
    enumeration_backend: str = "gmp-pointed-sieve"
    metric_name: str = "weighted-gram"
    quotient_weight: str = "1"
    chart_metric_kind: str = "metric"
    chart_metric_weight: str = "1"
    diversity_window: int = 1

    def __post_init__(self):
        for name in ('model_normalization', 'chart_parameterization', 'enumeration_backend', 'metric_name'):
            if not isinstance(getattr(self, name), str) or not getattr(self, name):
                raise ValueError("each search stage requires an explicit policy")
        chart_weight=Fraction(str(self.chart_metric_weight))
        if self.chart_metric_kind not in ("metric","gauss","raw") or chart_weight<=0:
            raise ValueError("invalid pointed-chart metric")
        object.__setattr__(self,"chart_metric_weight",str(chart_weight))
        weight = Fraction(str(self.quotient_weight))
        if weight <= 0 or type(self.diversity_window) is not int or self.diversity_window < 1:
            raise ValueError("metric weight and diversity window must be positive")
        object.__setattr__(self, 'quotient_weight', str(weight))

    @property
    def key(self):
        return digest(self)

    def metric(self, generic_gram, quotient_gram):
        """Sum orthogonal/projection component Grams with a declared weight.

        The caller supplies the component forms; this does not guess which
        coordinates belong to a generic subgroup. Positive definiteness of
        the resulting form is checked before enumeration.
        """
        n = len(generic_gram)
        if len(quotient_gram) != n or any(len(row) != n for row in (*generic_gram, *quotient_gram)):
            raise ValueError("metric component dimensions differ")
        weight = Fraction(self.quotient_weight)
        gram = [[Fraction(str(generic_gram[i][j]))+weight*Fraction(str(quotient_gram[i][j]))
                 for j in range(n)] for i in range(n)]
        return ldl(gram)[0]

    def holes(self, state, generic_gram, quotient_gram, *, checkpoint=None, target=None):
        binding = {"state": state.key, "policy": self.key}
        if checkpoint is not None:
            return VoronoiIterator.resume(checkpoint, binding=binding)
        return VoronoiIterator(self.metric(generic_gram, quotient_gram), target=target,
                               seen=state.parity.seen_holes, binding=binding)


class RepresentationPipeline:
    """Stage registries let a bounded benchmark vary one choice at a time."""
    def __init__(self, *, normalizers, parameterizations, enumerators):
        self.normalizers = dict(normalizers)
        self.parameterizations = dict(parameterizations)
        self.enumerators = dict(enumerators)

    def run(self, state, policy, *, centre, limits):
        stages = (("model_normalization", self.normalizers),
                  ("chart_parameterization", self.parameterizations),
                  ("enumeration_backend", self.enumerators))
        # Check the full combination before expensive arithmetic starts.
        for name, registry in stages:
            if getattr(policy, name) not in registry:
                raise ValueError(f"unsupported {name}: {getattr(policy, name)}")
        measurements = {}
        start = monotonic()
        normalized = self.normalizers[policy.model_normalization](state, limits)
        measurements['model_normalization_seconds'] = monotonic()-start
        start = monotonic()
        chart = self.parameterizations[policy.chart_parameterization](state, normalized, centre, policy, limits)
        measurements['chart_parameterization_seconds'] = monotonic()-start
        start = monotonic()
        updated, witnesses = self.enumerators[policy.enumeration_backend](state, chart, limits)
        measurements['enumeration_backend_seconds'] = monotonic()-start
        if updated.arithmetic != state.arithmetic or any(point not in updated.basis for point in state.basis):
            raise ArithmeticError("enumerator changed the curve or discarded a known point")
        return updated, {"policy": asdict(policy), "state_before": state.key, "state_after": updated.key,
                         "measurements": measurements, "witnesses": witnesses,
                         "mathematical_exclusion": False}


def calibration_protocol(*, panel, policies, limits, outcome_commitment, controls):
    """Freeze blind identities/configuration before opening recovery outcomes.

    Callers persist this record before dispatch. It deliberately excludes
    labels, target points and outcomes from policy inputs.
    """
    if not panel or len(set(panel)) != len(panel) or not policies or not limits or not outcome_commitment:
        raise ValueError("a blinded panel, policy sweep, budget and outcome commitment are required")
    if set(panel) & set(controls):
        raise ValueError("calibration and held-out controls must be disjoint")
    row = {"schema": "elliptic-curves.chart-policy-calibration.v1", "panel": list(panel),
           "policies": [asdict(p) for p in policies], "limits": limits,
           "outcome_commitment": outcome_commitment, "held_out_controls": list(controls),
           "ranking_rule": ["certified_independent_recoveries_descending", "total_wall_seconds_ascending", "policy_key"],
           "used_as_mathematical_exclusion": False}
    return {**row, "protocol_hash": digest(row)}


def rank_calibration(protocol, measurements):
    """Compare complete blind cells. Missing cells never become zero recovery."""
    row = dict(protocol); expected = row.pop('protocol_hash')
    if digest(row) != expected:
        raise ValueError("calibration protocol changed after commitment")
    policies = {digest(p): p for p in protocol['policies']}
    expected_cells = {(case, policy) for case in protocol['panel'] for policy in policies}
    cells = {(r['case'], r['policy']): r for r in measurements}
    if len(cells) != len(measurements) or set(cells) != expected_cells:
        raise ValueError("incomplete or duplicate blind policy sweep")
    result = []
    for key, policy in policies.items():
        entries = [cells[case, key] for case in protocol['panel']]
        for entry in entries:
            if (type(entry['certified_independent_recoveries']) is not int or entry['certified_independent_recoveries'] < 0
                or Fraction(str(entry['wall_seconds'])) < 0 or entry.get('protocol_hash') != expected):
                raise ValueError("invalid calibration measurement or changed protocol")
        result.append({"policy": policy, "policy_key": key,
            "certified_independent_recoveries": sum(e['certified_independent_recoveries'] for e in entries),
            "total_wall_seconds": str(sum(Fraction(str(e['wall_seconds'])) for e in entries))})
    return sorted(result, key=lambda r: (-r['certified_independent_recoveries'], Fraction(r['total_wall_seconds']), r['policy_key']))
