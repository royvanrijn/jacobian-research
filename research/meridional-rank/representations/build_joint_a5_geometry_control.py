#!/usr/bin/env python3
"""Second construction of the known P(2,3,3,5) diagram, with Spherogram."""

import hashlib
from importlib.metadata import version
import json
from pathlib import Path

from spherogram import RationalTangle

HERE=Path(__file__).resolve().parent


def main():
    assert version("spherogram")=="2.4.1"
    tangle=RationalTangle(0)
    for denominator in (2,3,3,5):
        tangle=tangle+RationalTangle(1,denominator)
    link=tangle.numerator_closure()
    assert len(link.link_components)==1
    pd=[list(c) for c in link.PD_code()]
    assert len(pd)==13
    result={"status":"CONSTRUCTED","definition":"Numerator closure of 1/2+1/3+1/3+1/5",
            "spherogram_version":version("spherogram"),"pd":pd,
            "input_sha256":{"build_joint_a5_geometry_control.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/"joint_a5_geometry_control.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"CONSTRUCTED","crossings":13,"components":1}))


if __name__=="__main__":
    main()
