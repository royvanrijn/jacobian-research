"""Remove two redundant carrier directions from the published 19-variable map."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    path=HERE/"sources/degree3-map-external.json";d=json.loads(path.read_text())
    old_names=["x","y","z"]+[f"w{i}" for i in range(1,17)]
    old=sp.symbols(old_names);env=dict(zip(old_names,old))
    original=sp.Matrix([sp.expand(sum(sp.Rational(c)*sp.sympify(m.replace("^","**"),locals=env) for m,c in row.items())) for row in d["components"]])
    x,y,z=old[:3];u,v=sp.symbols("u v")
    # u=6*w11-w12 and v=3*w14-w15. All remaining outputs
    # are independent of the two complementary carrier directions.
    rows=list(original[:3])+list(original[3:13])+[6*original[13]-original[14],original[15],3*original[16]-original[17],original[18]]
    sub={old[13]:0,old[14]:-u,old[16]:0,old[17]:-v}
    rows=[sp.expand(q.subs(sub,simultaneous=True)) for q in rows]
    names=old_names[:13]+["u","w13","v","w16"]
    syms=sp.symbols(names);assert len(rows)==len(syms)==17
    # Normalize the original 3x3 tangent block.
    rows=[rows[2]/2,rows[1],rows[0]]+rows[3:]
    J=sp.Matrix(rows).jacobian(syms);assert J.subs(dict.fromkeys(syms,0))==sp.eye(17)
    A=J[:3,:3];B=J[:3,3:];C=J[3:,:3];D=J[3:,3:];L=D-sp.eye(14)
    Di=sp.eye(14);term=sp.eye(14);nilpotence=None
    for exponent in range(1,15):
        term=(term*(-L)).applyfunc(sp.expand)
        if term==sp.zeros(14):nilpotence=exponent;break
        Di=(Di+term).applyfunc(sp.expand)
    assert nilpotence and (D*Di).applyfunc(sp.expand)==sp.eye(14)
    S=(A-B*Di*C).applyfunc(sp.expand)
    det=sp.expand(S.det(method="berkowitz"));assert det==1
    points=[]
    for p in d["points"]:
        vec={name:sp.Rational(c) for name,c in p.items()}
        vec["u"]=6*vec["w11"]-vec["w12"];vec["v"]=3*vec["w14"]-vec["w15"]
        points.append([vec[n] for n in names])
    images=[[sp.expand(f.subs(dict(zip(syms,p)))) for f in rows] for p in points]
    assert images[0]==images[1]==images[2] and len({tuple(p) for p in points})==3
    degrees=[sp.Poly(f,*syms).total_degree() for f in rows]
    coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in sp.Poly(f,*syms).terms()] for f in rows]
    counts=[sum(sum(t["exponent"])>=2 for t in row) for row in coefficients]
    witness=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),external_data_sha256=sha(path),
      dimension=17,degree=max(degrees),coordinate_degrees=degrees,nonlinear_counts=counts,nonlinear_total=sum(counts),
      variables=names,coefficients=coefficients,points=[[str(q) for q in p] for p in points],
      common_image=list(map(str,images[0])),determinant=str(det),identity_derivative=True,
      carrier_nilpotence_exponent=nilpotence,schur_complement=[[str(q) for q in row] for row in S.tolist()],
      transformations=dict(retained_u="6*w11-w12",retained_v="3*w14-w15",
        removed_complementary_directions=["w11","w14"],tangent_normalization="(G3/2,G2,G1,carriers)"),
      scope="Explicit cubic witness in dimension17; no minimum dimension or support claim.")
    (HERE/"cubic17-witness.json").write_text(json.dumps(witness,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in witness.items() if k not in ("coefficients","schur_complement","points")}))
if __name__=="__main__":main()

