#!/usr/bin/env python3
"""Re-run every stored characteristic-zero msolve proof certificate."""
import subprocess,tempfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
certdir=root/"results"/"certificates"
inputs=sorted(certdir.glob("*_Q.input"))
assert inputs
for source in inputs:
  expected=source.with_suffix(".msolve")
  with tempfile.NamedTemporaryFile(suffix=".out") as output:
    command=["msolve","-f",str(source),"-o",output.name,"-g","2","-v","0"]
    if "elimination" in source.stem:
      command += ["-e","8"]
    subprocess.run(command,check=True)
    text=Path(output.name).read_text()
  stored=expected.read_text()
  if "elimination" in source.stem:
    relation="8*C3^69*Fm4^3+18*d1*dm1^6*C3^23*Fm4+27*d0*dm1^9"
    assert relation in text and relation in stored,source
    print("PASS",source.name,"contains exact published elimination relation (5.9)")
  else:
    assert "[1]" in text and "[1]" in stored,source
    print("PASS",source.name,"has reduced Groebner basis [1]")
