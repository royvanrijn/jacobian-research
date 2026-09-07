"""Small, auditable wrapper around the msolve F4/F4SAT executable."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import shutil, subprocess, tempfile
import sympy as sp

@dataclass
class MSolveResult:
    characteristic:int; returncode:int; empty:bool; positive_dimensional:bool
    contains_one:bool; output:str; stderr:str; command:list[str]

def available(): return shutil.which("msolve") is not None

def _integral_string(f,variables):
    _,poly=sp.Poly(f,*variables,domain=sp.QQ).clear_denoms()
    return sp.sstr(poly.as_expr()).replace("**","^")

def input_text(equations,variables,characteristic=0):
    strings=[_integral_string(f,variables) for f in equations]
    return ",".join(map(str,variables))+"\n"+str(characteristic)+"\n"+",\n".join(strings)+"\n"

def run(equations,variables,prime=1000003,threads=4,saturate=None,
        groebner=True,timeout=None,eliminate=None):
    if not available(): raise RuntimeError("msolve executable not found")
    eqs=list(equations)
    if saturate is not None: eqs.append(saturate)
    with tempfile.TemporaryDirectory(prefix="jc-msolve-") as tmp:
        inp,out=Path(tmp)/"input.ms",Path(tmp)/"output.ms"
        inp.write_text(input_text(eqs,variables,prime))
        command=[shutil.which("msolve"),"-f",str(inp),"-o",str(out),"-t",str(threads),"-v","0"]
        if groebner: command += ["-g","2"]
        if eliminate is not None:
            if not 0 < int(eliminate) < len(variables):
                raise ValueError("eliminate must be between 1 and nvars-1")
            command += ["-e",str(int(eliminate))]
        if saturate is not None: command += ["-S"]
        proc=subprocess.run(command,text=True,capture_output=True,timeout=timeout)
        output=out.read_text() if out.exists() else ""
    stripped=output.strip()
    contains_one="[1]" in output or stripped=="[-1]"
    return MSolveResult(prime,proc.returncode,contains_one,stripped.startswith("[1,"),
                        contains_one,output,proc.stderr,command)
