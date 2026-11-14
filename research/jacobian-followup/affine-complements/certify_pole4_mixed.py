"""Exact localized identities for all mixed double/simple-root charts."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib
HERE=Path(__file__).resolve().parent
n=S.Symbol("n1");inv=S.Symbol("inverse_n1")
start=monotonic()
data=json.loads((HERE/"pole4_symbol_atlas.json").read_text())
out=[]
for row in data["cases"]:
    if row["kind"]!="double_simple":continue
    equations=[S.sympify(value) for degree,value in row["coefficient_equations"]]
    if row["name"]=="double_at_marked":
        weights=[S.Integer(degree==4) for degree,value in row["coefficient_equations"]]
        coefficient=-S.Rational(20,3)
    elif row["name"]=="simple_at_marked":
        weights=[S.Integer(degree==9) for degree,value in row["coefficient_equations"]]
        coefficient=S.Rational(20,3)
    else:
        weights=[S.Integer(degree) for degree,value in row["coefficient_equations"]]
        coefficient=-S.Rational(20,3)
    assert S.expand(sum(a*b for a,b in zip(weights,equations))-coefficient*n**4)==0
    multipliers=[S.expand(inv**4*weight/coefficient) for weight in weights]
    equations.append(1-inv*n)
    multipliers.append(1+inv*n+(inv*n)**2+(inv*n)**3)
    assert S.expand(sum(a*b for a,b in zip(equations,multipliers))-1)==0
    out.append(dict(name=row["name"],equations=list(map(str,equations)),
                    multipliers=list(map(str,multipliers)),
                    coefficient_functional="coefficient w^4" if row["name"]=="double_at_marked"
                      else "coefficient w^9" if row["name"]=="simple_at_marked"
                      else "derivative at w=1",
                    extracted_condition=str(coefficient*n**4)))
report=dict(as_of="2026-11-14",outcome="PASS",certificates=out,
            input_sha256=hashlib.sha256((HERE/"pole4_symbol_atlas.json").read_bytes()).hexdigest(),
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            elapsed_seconds=round(monotonic()-start,6))
(HERE/"pole4_mixed_certificates.json").write_text(json.dumps(report,indent=2)+"\n")
print("PASS exact localized identities for all three mixed-root charts")

