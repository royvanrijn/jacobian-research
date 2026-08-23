#!/usr/bin/env sage -python
"""
Bridge H92 native-q24 and current equation-q8 through the ABSTRACT adapted D13
frame, never through the mismatched ambient source-H3 coordinates.

If the two adapted positive-definite 17x17 Gram matrices coincide, coordinates
in that abstract frame are directly portable.  Then replay the native q24
divisor in the equation ambient marking and certify D12.

If they do not coincide, emit exact Gram/height diagnostics and stop before any
unsafe coordinate reuse.
"""

import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    h = Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
        h/"git"/"jacobian-research",
        h/"projects"/"jacobian-research",
    ]
    seen=set()
    for c in candidates:
        try: c=c.resolve()
        except Exception: continue
        if c in seen: continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")


ROOT=locate_repo()
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
EQ=SCRIPTS/"certify_h92_q8_equation_ns_divisor.sage"
Q24=SCRIPTS/"audit_h92_q8_q24_effective_zero_choices.sage"
OUT=LOCAL/"q8-q24-abstract-d13-bridge.json"
TMP=LOCAL/"q8-q24-abstract-d13-native-temp.json"


def run(path, argv):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    try:
        sys.argv=[str(path)]+list(argv)
        exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


print("Q8Q24ABS|stage=equation",flush=True)
eq=run(EQ,[])
print("Q8Q24ABS|stage=native_q24",flush=True)
qp=run(Q24,["--output",str(TMP)])

# Helpers from q24 producer, used symmetrically on both child lattices.
child_frame_with_zero=qp["child_frame_with_zero"]
d13_root_adaptation=qp["d13_root_adaptation"]
roots_and_data=qp["roots_and_data"]

# ---------------- equation frame ----------------
nsE=eq["ns"]
FE=vector(ZZ,eq["F8eq"])
OE=vector(ZZ,eq["selected"]["E8"]["simple_root_vectors_in_source_h3_ns"][0])
S3E=vector(ZZ,eq["S3std"])
assert OE*nsE*FE==1

childE,BzeroE=child_frame_with_zero(nsE,FE,OE)
AE,adaptE,HE=d13_root_adaptation(childE)
BadaptE=block_diagonal_matrix(identity_matrix(ZZ,2),AE)*BzeroE
GE=block_diagonal_matrix(matrix(ZZ,((0,1),(1,0))),-adaptE)
assert BadaptE*nsE*BadaptE.transpose()==GE

def ecoords(C):
    c=vector(QQ,C)*BadaptE.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ,c)

cS3=ecoords(S3E)
assert cS3[1]==52
zS3=vector(ZZ,cS3[-4:])

# ---------------- native q24 frame at II*_E8_1 ----------------
nsN=qp["source_ns"]
FN=vector(ZZ,qp["F_actual"])
DN=vector(ZZ,qp["D_actual"])
native_profile=next(p for p in qp["profiles"] if p["zero"]=="II*_E8_1")
ON=vector(ZZ,native_profile["zero_source_h3_ns"])
assert ON*nsN*FN==1

childN,BzeroN=child_frame_with_zero(nsN,FN,ON)
AN,adaptN,HN=d13_root_adaptation(childN)
BadaptN=block_diagonal_matrix(identity_matrix(ZZ,2),AN)*BzeroN
GN=block_diagonal_matrix(matrix(ZZ,((0,1),(1,0))),-adaptN)
assert BadaptN*nsN*BadaptN.transpose()==GN

def ncoords(C):
    c=vector(QQ,C)*BadaptN.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ,c)

cDN=ncoords(DN)
assert cDN[1]==2
zN=vector(ZZ,cDN[-4:])

adapt_equal=(adaptE==adaptN)
height_equal=(HE==HN)
mw_equal=(zS3==zN)

diff_count=0
diff_max=0
if not adapt_equal:
    for a,b in zip(adaptE.list(),adaptN.list()):
        if a!=b:
            diff_count+=1
            diff_max=max(diff_max,abs(int(a-b)))

print(
    "Q8Q24ABS_FRAME|"
    f"adapted_equal={int(adapt_equal)}|height_equal={int(height_equal)}|"
    f"gram_diff_entries={diff_count}|gram_diff_max={diff_max}|"
    f"AJ_mw={','.join(map(str,zS3))}|"
    f"native_q24_mw={','.join(map(str,zN))}|mw_equal={int(mw_equal)}|"
    f"native_Dcoords={','.join(map(str,cDN))}|"
    "status=PASS_DIAGNOSTIC",
    flush=True,
)

payload={
    "schema":"elkies-k3.h92-q8-q24-abstract-d13-bridge.v1",
    "adapted_equal":adapt_equal,
    "height_equal":height_equal,
    "AJ_S3_mw":list(map(int,zS3)),
    "native_q24_mw":list(map(int,zN)),
    "mw_equal":mw_equal,
    "native_q24_adapted_coordinates":list(map(int,cDN)),
}

if not adapt_equal:
    payload["status"]="ABSTRACT_D13_FRAMES_DIFFER"
    payload["gram_diff_entries"]=diff_count
    payload["gram_diff_max"]=diff_max
    payload["equation_height_gram"]=[[str(x) for x in row] for row in HE.rows()]
    payload["native_height_gram"]=[[str(x) for x in row] for row in HN.rows()]
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUT}",flush=True)
    print(
        "Q8Q24ABS_RESULT|"
        f"adapted_equal=0|height_equal={int(height_equal)}|mw_equal={int(mw_equal)}|"
        "status=NEED_ABSTRACT_D13_ISOMETRY",
        flush=True,
    )
    raise SystemExit(2)

# The full adapted Gram is identical, so the same integer coordinate row is a
# legitimate lattice class in the equation frame.
DE = cDN * BadaptE

assert DE*nsE*DE==0
assert DE*nsE*FE==2

mate=eq["isotropic_mate"](nsE,DE)
orth=matrix(ZZ,[list(DE*nsE),list(mate*nsE)]).right_kernel_matrix()
child=-(orth*nsE*orth.transpose())
rd=eq["roots_and_data"](child)[2]

print(
    "Q8Q24ABS_REPLAY|"
    f"D_square={DE*nsE*DE}|D_degree={DE*nsE*FE}|"
    f"root_data={rd[0]},{rd[1]},{rd[2]}|"
    f"MW_rank_if_rho19={17-rd[0]}|"
    f"status={'PASS_D12' if rd==(12,264,4) else 'UNEXPECTED_CHILD'}",
    flush=True,
)

# Compare native horizontal section and equation AJ more strongly.
PN=vector(ZZ,native_profile["section_source_h3_ns"])
cPN=ncoords(PN)
assert cPN[1]==1
zPN=vector(ZZ,cPN[-4:])
horizontal_equal=(zPN==zS3)

# Native vertical adapted coordinates are now also portable.
VE_native=DN-ON-PN
cV_native=vector(QQ,VE_native)*BadaptN.inverse()
assert all(x in ZZ for x in cV_native)
cV_native=vector(ZZ,cV_native)

# Replay native q24 section and vertical class in equation abstract frame.
PE=cPN*BadaptE
VE=cV_native*BadaptE
assert DE==OE+PE+VE

print(
    "Q8Q24ABS_HORIZONTAL|"
    f"AJ_mw={','.join(map(str,zS3))}|"
    f"q24_section_mw={','.join(map(str,zPN))}|"
    f"same={int(horizontal_equal)}|"
    f"native_height={native_profile['height']}|"
    f"native_corr={native_profile['D13_correction']}|"
    f"native_PdotO={native_profile['P_dot_O']}|status=PASS",
    flush=True,
)

print(
    "Q8Q24ABS_VERTICAL|"
    f"adapted_coords={','.join(map(str,cV_native))}|"
    f"V2={VE*nsE*VE}|VO={VE*nsE*OE}|VP={VE*nsE*PE}|"
    "status=PASS_EXACT_VERTICAL_REPLAY",
    flush=True,
)

payload.update({
    "status":"PASS_EXACT_Q24_ABSTRACT_D13_BRIDGE" if rd==(12,264,4) else "REPLAY_NOT_D12",
    "equation_q24_source_h3_ns":list(map(int,DE)),
    "equation_q24_section_source_h3_ns":list(map(int,PE)),
    "native_vertical_adapted_coordinates":list(map(int,cV_native)),
    "equation_vertical_source_h3_ns":list(map(int,VE)),
    "horizontal_equal":horizontal_equal,
    "child_root_data":list(map(int,rd)),
    "child_root_lattice":"D12" if rd==(12,264,4) else None,
    "MW_rank_if_rho19":17-int(rd[0]),
})

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)

if rd!=(12,264,4):
    print("Q8Q24ABS_RESULT|status=ABSTRACT_REPLAY_NOT_D12",flush=True)
    raise SystemExit(3)

print(
    "Q8Q24ABS_RESULT|"
    f"adapted_equal=1|height_equal={int(height_equal)}|"
    f"horizontal_same={int(horizontal_equal)}|"
    "child=D12|MW=5|"
    "status=PASS_EXACT_Q24_ABSTRACT_D13_BRIDGE",
    flush=True,
)
