#!/usr/bin/env python3
"""Generate and invariant-deduplicate short one/two-blow-up chart words."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.charts import enumerate_words,x,y
from jcsearch.canonical import chart_signature

parser=argparse.ArgumentParser();parser.add_argument("--depth",type=int,default=2)
parser.add_argument("--output",default="results/chart_catalog.json");args=parser.parse_args()
classes={}
failures=Path("results/scan_2d_5_20_collision.json")
penalty=1 if failures.exists() else 0
for chart in enumerate_words(args.depth):
  try: signature=chart_signature(chart.expressions,(x,y))
  except (sp.PolynomialError,ValueError,ZeroDivisionError): continue
  key=json.dumps(signature,sort_keys=True)
  item={"name":chart.name,"word":chart.word,"expressions":tuple(map(sp.sstr,chart.expressions)),
        "jacobian":sp.sstr(chart.jacobian),"signature":signature,
        "score":len(chart.word)+penalty}
  classes.setdefault(key,item)
records=sorted(classes.values(),key=lambda r:(r["score"],r["name"]))
path=Path(args.output);path.write_text(json.dumps(records,indent=2));print("raw words deduplicated to",len(records),"invariant classes; wrote",path)

