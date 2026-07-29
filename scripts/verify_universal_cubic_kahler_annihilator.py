#!/usr/bin/env python3
"""Universal locally free cubic Kähler-different/annihilator equality.

For the Deligne--Faddeev cubic algebra

    w^2 = -a*c + b*w - a*t,
    w*t = -a*d,
    t^2 = -b*d + d*w - c*t,

compute Omega over Q[a,b,c,d].  This checker verifies exactly that its
zeroth Fitting ideal (the Kähler different) equals Ann(Omega).

Applied on the punctured spectrum of the Koszul trace-free module, this
identifies the canonical different with the actual annihilator.  Extending
the equality across the collision origin still uses relative cotangent
saturation; this checker does not assert that separate universal theorem.
"""

from __future__ import annotations

import shutil
import subprocess


PROGRAM = r"""
LIB "homolog.lib";
ring polynomial_ring=0,(a,b,c,d,w,t),dp;
ideal algebra_relations=
  w2+a*c-b*w+a*t,
  w*t+a*d,
  t2+b*d-d*w+c*t;
qring cubic_algebra=std(algebra_relations);

module differentials=
  [2*w-b,a],
  [t,w],
  [-d,2*t+c];
differentials=std(differentials);

module annihilator_action=[1,0,0,1];
module annihilator_relations=
  [differentials[1][1],differentials[1][2],0,0],
  [differentials[2][1],differentials[2][2],0,0],
  [differentials[3][1],differentials[3][2],0,0],
  [0,0,differentials[1][1],differentials[1][2]],
  [0,0,differentials[2][1],differentials[2][2]],
  [0,0,differentials[3][1],differentials[3][2]];
module annihilator_preimage=std(
  modulo(annihilator_action,annihilator_relations)
);
ideal annihilator_ideal=std(ideal(annihilator_preimage));
ideal kahler_different=std(fitting(differentials,0));
ideal first_difference=simplify(
  reduce(annihilator_ideal,kahler_different),2
);
ideal second_difference=simplify(
  reduce(kahler_different,annihilator_ideal),2
);

print(
  "ANNIHILATOR_GENERATORS="
  +string(size(annihilator_ideal))
);
print(
  "KAHLER_DIFFERENT_GENERATORS="
  +string(size(kahler_different))
);
print(
  "IDEAL_DIFFERENCE="
  +string(size(first_difference)+size(second_difference))
);
quit;
"""


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "-q"],
        input=PROGRAM,
        text=True,
        capture_output=True,
        check=True,
    )
    values: dict[str, int] = {}
    wanted = {
        "ANNIHILATOR_GENERATORS",
        "KAHLER_DIFFERENT_GENERATORS",
        "IDEAL_DIFFERENCE",
    }
    for line in result.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in wanted:
            values[key] = int(value)
    assert values == {
        "ANNIHILATOR_GENERATORS": 3,
        "KAHLER_DIFFERENT_GENERATORS": 3,
        "IDEAL_DIFFERENCE": 0,
    }, (values, result.stdout, result.stderr)
    print(
        "PASS: the universal locally free cubic algebra has "
        "Fitt_0(Omega)=Ann(Omega)"
    )
    print(
        "PASS: the canonical different equals the actual annihilator "
        "on the punctured Koszul base"
    )
    print(
        "CONDITIONAL: relative cotangent saturation extends this equality "
        "across the collision origin"
    )


if __name__ == "__main__":
    main()
