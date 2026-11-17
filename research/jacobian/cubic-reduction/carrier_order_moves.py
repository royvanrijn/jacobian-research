#!/usr/bin/env python3
"""Five changed-order constructions moving a zero/positive weight carrier."""
from pathlib import Path
import json,hashlib
from carrier_order_search import run as model
from carrier_rank_search import run as rank
HERE=Path(__file__).resolve().parent
ORDERS=[[5,6,7,8,9,10,11,12,4],[12,4,5,6,7,8,9,10,11],[4,5,7,8,9,10,11,12,6],[12,11,10,9,8,7,6,5,4],[6,7,8,5,9,10,11,4,12]]
def main():
 cases=[]
 for order in ORDERS:
  folder=HERE/'carrier-order-search'/'_'.join(map(str,order))
  if (folder/'checkpoint.json').exists():raise RuntimeError('Existing exact handle/checkpoint must be inspected, not restarted: '+str(folder))
  model([j-1 for j in order]);rank(folder);cases.append(dict(order=order,folder=str(folder.relative_to(HERE))))
 (HERE/'carrier-order-moves-checkpoint.json').write_text(json.dumps(dict(status='complete',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=cases,scope='Five explicit changed-order models, each20-second primary identity/rank coefficient budget. Each per-case receipt is authoritative about its actual status. Not exhaustive orders or a global dimension bound.'),indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
