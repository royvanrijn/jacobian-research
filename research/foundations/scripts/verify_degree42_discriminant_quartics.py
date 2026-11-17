#!/usr/bin/env python3
"""Certify the quartic closure of the degenerate discriminant strata."""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def main() -> None:
    normals, bases, residuals, _defect = transformed_problem()
    e1, e2, translation, w0, _w1, w2 = bases
    a = sp.Symbol("sync42p_a")
    variables = normals + bases + (a,)
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"

    program = f"""
ring q=0,({",".join(map(str, variables))}),(dp(5),dp(7));
ideal I={",".join(serialize(item) for item in residuals)};

proc evt(poly p)
{{
  return(subst(p,{zero_normals},{w0},0,{e1},{a}^2,
    {e2},2*{a},{translation},0));
}}
proc Dt(poly p)
{{
  return({a}*diff(p,{normals[0]})+diff(p,{normals[1]}));
}}
poly ta34=evt(diff(I[5],{normals[3]}));
poly ta35=evt(diff(I[5],{normals[4]}));
poly ta45=evt(diff(I[11],{normals[4]}));
poly th3=evt(Dt(Dt(I[5])))/2;
poly th4=evt(Dt(Dt(I[11])))/2;
poly th5=evt(Dt(Dt(I[17])))/2;
poly tz25=th5;
poly tz24=th4+ta45*tz25;
poly tz23=th3+ta34*tz24+ta35*tz25;
proc tcross(poly p)
{{
  return(tz23*evt(Dt(diff(p,{normals[2]})))
    +tz24*evt(Dt(diff(p,{normals[3]})))
    +tz25*evt(Dt(diff(p,{normals[4]}))));
}}
poly tk3=evt(Dt(Dt(Dt(I[5]))))/6+tcross(I[5]);
poly tk4=evt(Dt(Dt(Dt(I[11]))))/6+tcross(I[11]);
poly tk5=evt(Dt(Dt(Dt(I[17]))))/6+tcross(I[17]);
poly tz35=tk5;
poly tz34=tk4+ta45*tz35;
poly tz33=tk3+ta34*tz34+ta35*tz35;
proc tquartic(poly p)
{{
  poly ans=evt(Dt(Dt(Dt(Dt(p)))))/24;
  ans=ans+tz23*evt(Dt(Dt(diff(p,{normals[2]})))/2)
    +tz24*evt(Dt(Dt(diff(p,{normals[3]})))/2)
    +tz25*evt(Dt(Dt(diff(p,{normals[4]})))/2);
  ans=ans+tz33*evt(Dt(diff(p,{normals[2]})))
    +tz34*evt(Dt(diff(p,{normals[3]})))
    +tz35*evt(Dt(diff(p,{normals[4]})));
  ans=ans+tz23^2*evt(diff(diff(p,{normals[2]}),{normals[2]}))/2
    +tz24^2*evt(diff(diff(p,{normals[3]}),{normals[3]}))/2
    +tz25^2*evt(diff(diff(p,{normals[4]}),{normals[4]}))/2;
  ans=ans+tz23*tz24*evt(diff(diff(p,{normals[2]}),{normals[3]}))
    +tz23*tz25*evt(diff(diff(p,{normals[2]}),{normals[4]}))
    +tz24*tz25*evt(diff(diff(p,{normals[3]}),{normals[4]}));
  return(ans);
}}
poly tq18=tquartic(I[18]);
poly tq19=tquartic(I[19]);
int t_zero_quartic=(
  tq18==15/16*{w2}*{a} and tq19==5/64*{w2}
);

proc evc(poly p)
{{
  return(subst(p,{zero_normals},{w0},0,{e1},0,{translation},0));
}}
proc Dc(poly p)
{{
  return({e2}*diff(p,{normals[0]})+diff(p,{normals[1]}));
}}
poly ca34=evc(diff(I[5],{normals[3]}));
poly ca35=evc(diff(I[5],{normals[4]}));
poly ca45=evc(diff(I[11],{normals[4]}));
poly ch3=evc(Dc(Dc(I[5])))/2;
poly ch4=evc(Dc(Dc(I[11])))/2;
poly ch5=evc(Dc(Dc(I[17])))/2;
poly cz25=ch5;
poly cz24=ch4+ca45*cz25;
poly cz23=ch3+ca34*cz24+ca35*cz25;
proc ccross(poly p)
{{
  return(cz23*evc(Dc(diff(p,{normals[2]})))
    +cz24*evc(Dc(diff(p,{normals[3]})))
    +cz25*evc(Dc(diff(p,{normals[4]}))));
}}
poly ck3=evc(Dc(Dc(Dc(I[5]))))/6+ccross(I[5]);
poly ck4=evc(Dc(Dc(Dc(I[11]))))/6+ccross(I[11]);
poly ck5=evc(Dc(Dc(Dc(I[17]))))/6+ccross(I[17]);
poly cz35=ck5;
poly cz34=ck4+ca45*cz35;
poly cz33=ck3+ca34*cz34+ca35*cz35;
proc cquartic(poly p)
{{
  poly ans=evc(Dc(Dc(Dc(Dc(p)))))/24;
  ans=ans+cz23*evc(Dc(Dc(diff(p,{normals[2]})))/2)
    +cz24*evc(Dc(Dc(diff(p,{normals[3]})))/2)
    +cz25*evc(Dc(Dc(diff(p,{normals[4]})))/2);
  ans=ans+cz33*evc(Dc(diff(p,{normals[2]})))
    +cz34*evc(Dc(diff(p,{normals[3]})))
    +cz35*evc(Dc(diff(p,{normals[4]})));
  ans=ans+cz23^2*evc(diff(diff(p,{normals[2]}),{normals[2]}))/2
    +cz24^2*evc(diff(diff(p,{normals[3]}),{normals[3]}))/2
    +cz25^2*evc(diff(diff(p,{normals[4]}),{normals[4]}))/2;
  ans=ans+cz23*cz24*evc(diff(diff(p,{normals[2]}),{normals[3]}))
    +cz23*cz25*evc(diff(diff(p,{normals[2]}),{normals[4]}))
    +cz24*cz25*evc(diff(diff(p,{normals[3]}),{normals[4]}));
  return(ans);
}}
poly cq18=cquartic(I[18]);
poly cq19=cquartic(I[19]);
int cusp_quartic=(
  cq18==5/8*{e2}*{w2} and cq19==5/64*{w2}
);

print("DEGREE42_DISCRIMINANT_QUARTICS");
print(t_zero_quartic);
print(cusp_quartic);
"""
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    compact = process.stdout.split()
    marker = compact.index("DEGREE42_DISCRIMINANT_QUARTICS")
    assert compact[marker + 1 : marker + 3] == ["1", "1"], (
        process.stdout + process.stderr
    )
    print(
        "PASS: the t=0 and cusp discriminant strata have terminal "
        "quartic coefficient 5*w2/64 and synchronize on D(w1*w2)"
    )


if __name__ == "__main__":
    main()
