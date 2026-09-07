"""Theorem-bound pruning; missing, stale and heuristic evidence never excludes."""

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

from .binary import BinaryBasis
from .store import digest


def theorem_fingerprint(entry):
    # Bind the exact mathematical statement, not unrelated navigation changes.
    return digest({key:entry.get(key) for key in
        ("id","state","scope","canonical_source","artifact_hash","dependencies")})


@dataclass(frozen=True)
class SearchRequest:
    kind: str
    scope: tuple[tuple[str,str], ...]
    target_rank: int = 1
    class_mask: int | None = None
    class_dimension: int | None = None
    height: str | None = None

    def __post_init__(self):
        object.__setattr__(self,"scope",tuple(tuple(row) for row in self.scope))
        if len(dict(self.scope))!=len(self.scope):raise ValueError("duplicate search scope keys")
        if self.kind not in ("nontorsion_section","rank_target","kummer_class","birational_chart"):
            raise ValueError("search requests need explicit mathematical semantics")
        if type(self.target_rank) is not int or self.target_rank<1:raise ValueError("positive target rank required")
        if self.kind=="kummer_class":
            if type(self.class_dimension) is not int or self.class_dimension<0:
                raise ValueError("labelled squareclass dimension required")
            if type(self.class_mask) is not int or not 0<=self.class_mask<1<<self.class_dimension:
                raise ValueError("class mask outside its declared subspace")


class PruningRegistry:
    def __init__(self,root:Path,registry_path=None):
        self.root=Path(root)
        self.registry_path=Path(registry_path) if registry_path else self.root/"elliptic-curves/data/search_constraints_v1.json"
        self._stamp=None
        self._rules=()
        self.inactive=[]
        self._watched=(self.root/"MATH_STATUS.json",self.registry_path)

    def _path(self,name):
        path=(self.root/name).resolve()
        if not path.is_relative_to(self.root.resolve()):raise ValueError("proof path outside repository")
        return path

    def _load(self):
        def stamp():
            result=[]
            for path in self._watched:
                try:
                    stat=path.stat();result.append((str(path),stat.st_mtime_ns,stat.st_ctime_ns,stat.st_size,stat.st_ino))
                except FileNotFoundError:result.append((str(path),None))
            return tuple(result)
        if self._stamp==stamp():return
        status_path=self.root/"MATH_STATUS.json"
        status=json.loads(status_path.read_text())
        entries={entry["id"]:entry for entry in status["entries"]}
        records=json.loads(self.registry_path.read_text())["constraints"] if self.registry_path.exists() else []
        # New theorem entries may carry constraints directly in the authority.
        # They become visible on the very next request without rebuilding a queue.
        for entry in entries.values():
            for rule in entry.get("search_constraints",[]):
                records.append({**rule,"theorem_id":entry["id"],"theorem_fingerprint":theorem_fingerprint(entry)})
        rules=[];self.inactive=[]
        for rule in records:
            entry=entries.get(rule.get("theorem_id"))
            reason=None
            if entry is None or entry.get("state")!="proved":reason="theorem is not proved in the authority"
            elif theorem_fingerprint(entry)!=rule.get("theorem_fingerprint"):reason="mathematical statement changed"
            elif not rule.get("evidence"):reason="no pinned proof witnesses"
            else:
                for name,expected in rule["evidence"].items():
                    try:actual=sha256(self._path(name).read_bytes()).hexdigest()
                    except FileNotFoundError:actual=None
                    if actual!=expected:
                        reason="proof witness missing or changed";break
            if reason:
                self.inactive.append({"theorem_id":rule.get("theorem_id"),"reason":reason})
            else:rules.append(rule)
        self._rules=tuple(rules)
        self._watched=tuple(sorted({status_path,self.registry_path,
            *(self._path(name) for rule in records for name in rule.get("evidence",{}))}))
        self._stamp=stamp()

    def matching_rules(self,scope):
        self._load()
        scope=dict(scope)
        return tuple(rule for rule in self._rules if all(scope.get(key)==value for key,value in rule["scope"].items()))

    def decision(self,request:SearchRequest):
        if not isinstance(request,SearchRequest):raise TypeError("typed search request required")
        # Authority reads and witness hashes are cheap compared with a worker.
        # Deliberately do not cache a proved exclusion across file changes.
        self._load()
        scope=dict(request.scope)
        exclusions=[]
        for rule in self._rules:
            if any(scope.get(key)!=value for key,value in rule["scope"].items()):continue
            kind=rule["kind"]
            if kind=="rank_upper" and request.kind in ("nontorsion_section","rank_target"):
                bound=rule["upper"]
                if type(bound) is not int or bound<0:raise ValueError("invalid theorem rank bound")
                if request.target_rank>bound:exclusions.append(rule["theorem_id"])
            elif kind=="allowed_squareclass_subspace" and request.kind=="kummer_class":
                if request.class_dimension!=rule["dimension"]:continue
                basis=BinaryBasis(rule["dimension"])
                for mask in rule["basis_masks"]:basis,_=basis.append(mask)
                if basis.reduce(request.class_mask)[0]:exclusions.append(rule["theorem_id"])
            elif kind not in ("rank_upper","allowed_squareclass_subspace","scheduling_only"):
                raise ValueError("unknown pruning-rule semantics")
        return {"search_allowed":not exclusions,"status":"EXCLUDED_BY_THEOREM" if exclusions else "UNKNOWN",
                "theorems":exclusions,"inactive_constraints":list(self.inactive)}

    def filter(self,requests):
        """Evaluate the mathematical gate before a caller constructs any chart."""
        for request in requests:
            if self.decision(request)["search_allowed"]:yield request


def known_constraints(root):
    """Import existing exact conclusions; never infer a theorem from a score."""
    from .arithmetic import CurveModel,TwoTorsionContext
    from .regulator import Surface
    from .subspace import local_intersection,restricted_radical
    root=Path(root)
    entries={r["id"]:r for r in json.loads((root/"MATH_STATUS.json").read_text())["entries"]}
    results=root/"artifacts/generated-results"
    rules=[]
    def add(theorem,scope,kind,files,**data):
        entry=entries[theorem]
        files=[*files,root/entry["canonical_source"]]
        rules.append({"theorem_id":theorem,"theorem_fingerprint":theorem_fingerprint(entry),
            "scope":scope,"kind":kind,**data,
            "evidence":{str(p.relative_to(root)):sha256(p.read_bytes()).hexdigest() for p in files}})
    direct_path=results/"elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
    pairs_path=results/"elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
    direct=json.loads(direct_path.read_text())["weierstrass_model"]
    pairs={r["pair_key"]:r for r in json.loads(pairs_path.read_text())["pairs"]}
    original=results/"elkies-k3-r17-product-19bad-083ad-rank-zero-v1.json"
    sweep=results/"elkies-k3-r17-product-regulator-sweep-v1.json"
    old=json.loads(original.read_text());new=json.loads(sweep.read_text())
    targets=[(old["pair_key"],"EC-K3-R17-PRODUCT-19BAD-083AD-ARITHMETIC-RANK-ZERO",original,old["rank_over_QQ_u"])]
    targets.extend((row["pair_key"],"EC-K3-R17-PRODUCT-REGULATOR-OBSTRUCTION-SWEEP",sweep,row["rank_over_QQ_u"]["upper"]) for row in new["targets"])
    for key,theorem,path,upper in targets:
        surface=Surface(tuple(direct["A_coefficients_low_to_high"]),tuple(direct["B_coefficients_low_to_high"]),
                        tuple(pairs[key]["product_quartic_coefficients_low_to_high"]))
        add(theorem,{"surface":surface.key},"rank_upper",[path,direct_path,pairs_path],upper=upper)
    source=results/"elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json"
    ct=results/"elliptic-curves/fixed_cubic_u_minus1_cassels_tate_v1.json"
    local=json.loads(source.read_text());pairing=json.loads(ct.read_text())
    row=next(r for r in local["runs"] if r["parameter_u"]=="-1")
    basis=local["anchor"]["known_kummer_basis_beta_power_coordinates"]
    algebra=TwoTorsionContext(tuple(local["anchor"]["base_polynomial_ascending"]))
    space=digest({"algebra":algebra.key,"basis":basis})
    maps=[r["known_span_quotient_rows"] for r in row["finite_local_conditions"]]+[row["real_local_condition"]["known_span_quotient_rows"]]
    masks=local_intersection(len(basis),maps)
    radical=restricted_radical(masks,pairing["arithmetic"]["matrix"],global_dimension=len(basis))
    if radical["obstructed_class_count"]!=pairing["arithmetic"]["obstructed_class_count"]:
        raise ArithmeticError("CT summary does not match its exact matrix")
    add("EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE",
        {"curve":CurveModel(tuple(row["raw_curve_ainvariants"])).key,"squareclass_space":space},
        "allowed_squareclass_subspace",[source,ct,root/pairing["evidence"]],
        dimension=len(basis),basis_masks=radical["radical_global_masks"])
    return {"schema":"elliptic-curves.search-constraints.v1","constraints":rules}
