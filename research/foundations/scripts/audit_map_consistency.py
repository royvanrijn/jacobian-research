#!/usr/bin/env python3
"""Mechanically audit every direct restatement of the announced 3D map.

The manifest is intentional: adding a new consumer requires deciding whether
it restates, normalizes, specializes, or derives the map. Python expressions
are read from their AST without importing the consumer modules.
"""

from __future__ import annotations

import ast
import contextlib
import importlib.util
import io
import operator
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
x, y, z = sp.symbols("x y z")
u = 1 + x*y
CANONICAL = (
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
)
ENVIRONMENT = {"x": x, "y": y, "z": z, "u": u, "w": u}


def assignment(path: str, name: str) -> ast.expr:
    tree = ast.parse((ROOT / path).read_text(), filename=path)
    matches = [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == name
                for target in node.targets)
    ]
    assert len(matches) == 1, (path, name, len(matches))
    return matches[0]


def component_nodes(node: ast.expr) -> list[ast.expr]:
    if isinstance(node, ast.Call):
        assert len(node.args) == 1
        node = node.args[0]
    assert isinstance(node, (ast.List, ast.Tuple)), ast.dump(node)
    return list(node.elts)


def expression(node: ast.expr, extra: dict[str, object] | None = None) -> sp.Expr:
    allowed = (
        ast.Expression, ast.BinOp, ast.UnaryOp, ast.Name, ast.Constant,
        ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd,
        ast.Load,
    )
    assert all(isinstance(part, allowed) for part in ast.walk(node)), ast.dump(node)
    environment = dict(ENVIRONMENT)
    if extra:
        environment.update(extra)
    return sp.sympify(eval(
        compile(ast.Expression(node), "<map-audit>", "eval"),
        {"__builtins__": {}},
        environment,
    ))


def extracted_map(path: str, name: str = "F") -> tuple[sp.Expr, ...]:
    return tuple(expression(node) for node in component_nodes(assignment(path, name)))


def assert_same(path: str, actual: tuple[sp.Expr, ...]) -> None:
    differences = tuple(sp.expand(got - want)
                        for got, want in zip(actual, CANONICAL))
    assert len(actual) == 3 and differences == (0, 0, 0), (path, differences)
    print("PASS map formula:", path)


DIRECT_MAPS = (
    "archive/tooling/analyze_gradient_infinity.py",
    "archive/tooling/commuting_flows.py",
    "scripts/image_nonproperness.py",
    "archive/tooling/keller_tangent.py",
    "archive/tooling/nonproper_fiber_benchmark.py",
    "scripts/verify_counterexample.py",
    "scripts/verify_marked_root_model.py",
)
for source in DIRECT_MAPS:
    assert_same(source, extracted_map(source))

assert_same(
    "scripts/cubic_homogeneous_reduction.py:announced",
    extracted_map("scripts/cubic_homogeneous_reduction.py", "announced"),
)
assert_same(
    "scripts/cubic_model.py:(a,b,c)",
    tuple(expression(assignment("scripts/cubic_model.py", name))
          for name in ("a", "b", "c")),
)

# The finite-field evaluators express the same formula through modular or
# abstract field operations. Compare them on a grid in several characteristics.
def load_functions(path: str):
    spec = importlib.util.spec_from_file_location("map_audit_consumer", ROOT / path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_value(point: tuple[int, int, int]) -> tuple[int, int, int]:
    substitutions = dict(zip((x, y, z), point))
    return tuple(int(component.subs(substitutions)) for component in CANONICAL)


distribution = load_functions("archive/tooling/finite_field_distribution.py")
refinements = load_functions("archive/tooling/finite_field_refinements.py")
for point in ((a, b, c) for a in range(-2, 3)
              for b in range(-2, 3) for c in range(-2, 3)):
    expected = canonical_value(point)
    abstract = distribution.image_point(
        *point, operator.add, operator.neg, operator.mul, 1
    )
    assert abstract == expected, (point, abstract, expected)
    for prime in (2, 3, 5, 1009):
        reduced = tuple(value % prime for value in expected)
        assert refinements.evaluate(point, prime) == reduced
print("PASS map formula: finite-field evaluators")

# The dependency-free verifier stores the expanded integer coefficient maps.
# Execute it in an isolated namespace and compare those coefficient dictionaries
# to SymPy's independent expansion of the canonical formula.
independent_source = ROOT / "scripts/verify_counterexample_independent.py"
namespace = {"__name__": "map_consistency_audit"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(independent_source.read_text(), str(independent_source), "exec"), namespace)
independent_mapping = namespace["mapping"]
for got, want in zip(independent_mapping, CANONICAL):
    assert got == sp.Poly(sp.expand(want), x, y, z).as_dict()
print("PASS map formula: dependency-free verifier")

# Julia's homotopy system subtracts a moving target from the first coordinate.
julia = "".join((ROOT / "archive/tooling/nonproper_fiber_homotopy.jl").read_text().split())
for fragment in (
    "u=1+x*y",
    "u^3*z+y^2*u*(4+3*x*y)-(-1//4+s)",
    "y+3*x*u^2*z+3*x*y^2*(4+3*x*y)",
    "2*x-3*x^2*y-x^3*z",
):
    assert fragment in julia, fragment
print("PASS map formula: archive/tooling/nonproper_fiber_homotopy.jl")

# Human-facing full displays are checked textually; allow harmless prose changes
# in the sentence introducing the fixed formula.
display = "".join([
    "(1+xy)^3z+y^2(1+xy)(4+3xy),",
    "y+3x(1+xy)^2z+3xy^2(4+3xy),",
    "2x-3x^2y-x^3z",
])
compact = "".join((ROOT / "archive/legacy-notes/FINITE_FIELD_VALUE_DISTRIBUTION.md").read_text().split())
assert display in compact, "archive/legacy-notes/FINITE_FIELD_VALUE_DISTRIBUTION.md"
readme = "".join((ROOT / "README.md").read_text().split())
assert any(
    fragment in readme
    for fragment in (
        "Let`u=1+xy`anddefine`F:A^3->A^3`by",
        "Put`u=1+xy`anddefine",
    )
), "README formula introduction"
for fragment in (
    "u^3z+y^2u(4+3xy)",
    "y+3xu^2z+3xy^2(4+3xy)",
    "2x-3x^2y-x^3z",
):
    assert fragment in readme, fragment
facts = "".join((ROOT / "archive/core-support/FACTS.md").read_text().split())
for fragment in (
    "Put`u=1+xy`andwrite`F=(a,b,c)`",
    "a=u^3z+y^2u(4+3xy)",
    "b=y+3xu^2z+3xy^2(4+3xy)",
    "c=2x-3x^2y-x^3z",
):
    assert fragment in facts, fragment
print("PASS map formula: README and facts notes")

print("PASS: all inventoried map restatements agree in signs and coordinate order")
