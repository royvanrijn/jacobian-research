
#!/usr/bin/env sage -python
import contextlib, io, json, sys
from pathlib import Path
from sage.all import QQ, ZZ, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
CLOSE=SCRIPTS/"close_h92_q8_q24_by_q6_translation.sage"
OUTFRAME=LOCAL/"route-scout/q32-d12-mw5-root-adapted-frame.txt"
OUTMETA=LOCAL/"route-scout/q32-d12-mw5-frame.json"

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

def run_scope(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    buf=io.StringIO()
    try:
        sys.argv=[str(path)]
        with contextlib.redirect_stdout(buf):
            exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope

cl=run_scope(CLOSE)
need=("ns","F8eq","O8","Badapt","adapted","H13")
missing=[k for k in need if k not in cl]
if missing:
    raise SystemExit("close scope missing "+",".join(missing))

ns=cl["ns"]
F=vector(ZZ,cl["F8eq"])
O=vector(ZZ,cl["O8"])
Badapt=cl["Badapt"]
G=cl["adapted"]
H=cl["H13"]

# Reconstruct the exact q32 class from its intrinsic current-D13 data.
R=G[:13,:13]
C=G[:13,13:]
graph={i:set() for i in range(13)}
for i in range(13):
    assert R[i,i]==2
    for j in range(i+1,13):
        assert R[i,j] in (0,-1)
        if R[i,j]==-1:
            graph[i].add(j); graph[j].add(i)
branch=[i for i in graph if len(graph[i])==3]
assert len(branch)==1
spin=[j for j in graph[branch[0]] if len(graph[j])==1]
assert len(spin)==2

z=vector(ZZ,(-2,1,-1,1))
assert z*H*z==52
labels=vector(ZZ,[ZZ(i in spin) for i in range(13)])

rhs=vector(QQ,labels)-vector(QQ,z)*C.transpose()
r=rhs*R.inverse()
assert all(v in ZZ for v in r)
r=vector(ZZ,r)
w=vector(ZZ,list(r)+list(z))
norm=ZZ(w*G*w)
assert norm==64
Dcoords=vector(ZZ,[16,2]+list(w))
Gns=cl["Gadapt"]
assert Dcoords*Gns*Dcoords==0
D32=vector(ZZ,Dcoords*Badapt)
assert D32*ns*D32==0
assert D32*ns*F==2
assert D32*ns*O==14

# Exact U-neighbor child and root-adapted minimization.
split=primitive_hyperbolic_split(ns,D32)
raw_child=matrix(ZZ,split["child_frame"])
minimized=minimize_child_frame(raw_child)
frame=matrix(ZZ,minimized["frame"])
root_data=tuple(map(int,minimized["root_data"]))
assert root_data==(12,264,4),root_data
assert frame.det()==948

root_rank=12
root=frame[:root_rank,:root_rank]
coupling=frame[:root_rank,root_rank:]
tail=frame[root_rank:,root_rank:]
height=tail-coupling.transpose()*root.inverse()*coupling
assert height.det()==QQ(237)

OUTFRAME.parent.mkdir(parents=True,exist_ok=True)
OUTFRAME.write_text(
    "# H3 q32 geometrically clean D12/MW5 child\n"
    "# exact lattice class in current canonical D13 equation frame\n"
    "# q=32, factor order (16,2), old-fibre degree 2\n"
    + "\n".join(" ".join(map(str,row)) for row in frame.rows()) + "\n"
)

# Compare with historical orbit85 D12 frame if available.
Q24=GEN/"elkies-k3-h3-q6-q8-d13-q24-degree2.json"
hist_equal=False
hist_height=None
hist_isometric_height=None
if Q24.exists():
    data=json.loads(Q24.read_text())
    hit=next(r for r in data["neighbors"] if int(r["orbit_index"])==85)
    hist=matrix(ZZ,hit["child_root_adapted_frame"])
    hist_equal=(hist==frame)
    Hr=hist[:12,:12]
    Hc=hist[:12,12:]
    Ht=hist[12:,12:]
    hist_height=Ht-Hc.transpose()*Hr.inverse()*Hc
    # rank-5 qfisom is cheap enough at scaled integral level.
    from sage.all import lcm, pari
    scale=ZZ(1)
    for v in list(height.list())+list(hist_height.list()):
        scale=lcm(scale,ZZ(QQ(v).denominator()))
    try:
        iso=pari((scale*height).change_ring(ZZ)).qfisom(
            pari((scale*hist_height).change_ring(ZZ))
        )
        hist_isometric_height=(str(iso)!="0")
    except Exception:
        hist_isometric_height=None

payload={
    "schema":"elkies-k3.h3-q32-d12-mw5-exact-frame.v1",
    "status":"PASS_EXACT_Q32_D12_MW5_FRAME",
    "q":32,
    "factor_order":[16,2],
    "old_fibre_degree":2,
    "D_dot_O":14,
    "spinor_nodes_1_based":[i+1 for i in spin],
    "mw_coordinates":[-2,1,-1,1],
    "root_coordinates":list(map(int,r)),
    "root_data":list(root_data),
    "mw_height":[[str(v) for v in row] for row in height.rows()],
    "frame":str(OUTFRAME.relative_to(ROOT)),
    "historical_orbit85_frame_equal":hist_equal,
    "historical_orbit85_MW_height_isometric":hist_isometric_height,
}
OUTMETA.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(
    "Q32D12FRAME|"
    f"spinors={','.join(str(i+1) for i in spin)}|"
    f"DdotO=14|root_data=12,264,4|MW=5|"
    f"height_diag={','.join(str(height[i,i]) for i in range(5))}|"
    f"hist_frame_equal={int(hist_equal)}|"
    f"hist_height_isometric={hist_isometric_height}|"
    "status=PASS_EXACT_Q32_D12_MW5_FRAME",
    flush=True,
)
print(f"OUTPUT|{OUTFRAME}",flush=True)
print(f"OUTPUT|{OUTMETA}",flush=True)
