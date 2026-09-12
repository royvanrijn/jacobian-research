#!/usr/bin/env python3
"""Logged access to the public Magma calculator, with its standard limits.

Reads only a self-generated request. Stores exact request/response and timing.
This wrapper does not manufacture, discover, or import any rational point.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import resource
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'research/artifacts/local/elliptic-curves/constructed-class-blind-cover-v1'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('request')
    ap.add_argument('--kind', choices=['preparation', 'search'], default='preparation')
    ns = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    request_path = Path(ns.request).resolve()
    assert request_path.parent == OUT
    code = request_path.read_text()
    assert len(code.encode()) <= 50000
    assert 'assert GetMemoryLimit() gt 0 and GetMemoryLimit() le 2147483648;' in code
    rows = [json.loads(s) for s in (OUT / 'arithmetic-ledger.jsonl').read_text().splitlines()]
    used = sum(float(r['arithmetic_wall_seconds']) for r in rows)
    assert used + 65 <= 1200, used
    record = dict(step=request_path.stem, kind=ns.kind, request_sha256=hashlib.sha256(code.encode()).hexdigest(),
                  started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  public_url='https://magma.maths.usyd.edu.au/xml/calculator.xml',
                  service_limit_seconds=60, service_input_limit_bytes=50000,
                  memory_limit_ceiling_bytes=2147483648,
                  memory_enforcement='assert existing remote memory limit is positive and <= 2 GiB',
                  input_bytes=len(code.encode()))
    request = urllib.request.Request(record['public_url'],
                data=urllib.parse.urlencode({'input': code}).encode(),
                headers={'Content-Type':'application/x-www-form-urlencoded',
                         'Referer':'https://magma.maths.usyd.edu.au/calc/'})
    start = time.monotonic()
    try:
        with urllib.request.urlopen(request, timeout=65) as response:
            raw = response.read()
            record['http_status'] = response.status
    except Exception as exc:
        raw = (type(exc).__name__ + ': ' + str(exc)).encode()
        record['error'] = raw.decode()
    record['request_wall_seconds'] = time.monotonic() - start
    request_path.with_suffix('.response.xml').write_bytes(raw)
    try:
        root = ET.fromstring(raw)
        result = '\n'.join(''.join(line.itertext()) for line in root.findall('.//results/line'))
        headers = {e.tag: e.text for e in root.findall('.//headers/*')}
        record['headers'] = headers
        record['reported_arithmetic_cpu_seconds'] = float(headers.get('time', 60))
    except Exception:
        result = raw.decode(errors='replace')
    # Charge full request elapsed time conservatively, including network wait.
    record['arithmetic_wall_seconds'] = record['request_wall_seconds']
    request_path.with_suffix('.result.txt').write_text(result)
    request_path.with_suffix('.metadata.json').write_text(json.dumps(record, indent=2) + '\n')
    with (OUT / 'arithmetic-ledger.jsonl').open('a') as log:
        log.write(json.dumps(record) + '\n')
    print(json.dumps(record, indent=2), flush=True)
    print(result, flush=True)


if __name__ == '__main__':
    main()
