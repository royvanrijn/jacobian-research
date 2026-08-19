from pathlib import Path
import json
import re
import subprocess
import time
import urllib.request
import urllib.parse

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "archive-hunt"
OUT.mkdir(parents=True, exist_ok=True)

HOSTS = [
    "math.harvard.edu/~elkies/*",
    "www.math.harvard.edu/~elkies/*",
    "people.math.harvard.edu/~elkies/*",
]

TERMS = [
    "6,79",
    "6*79",
    "474",
    "1311",
    "rank 17",
    "rank17",
    "rk17",
    "mordell-weil rank 17",
    "shimura",
    "k3",
    "28 independent",
]

def get(url):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0"}
    )
    with urllib.request.urlopen(req, timeout=30) as f:
        return f.read()

rows = []

for host in HOSTS:
    q = urllib.parse.urlencode({
        "url": host,
        "output": "json",
        "filter": "statuscode:200",
        "collapse": "urlkey",
    })

    url = "https://web.archive.org/cdx/search/cdx?" + q
    print("CDX", host, flush=True)

    try:
        data = json.loads(get(url))
    except Exception as e:
        print("FAILED", host, repr(e), flush=True)
        continue

    if not data:
        continue

    header = data[0]

    for r in data[1:]:
        d = dict(zip(header, r))
        rows.append(d)

# Dedup by original URL.
uniq = {}
for d in rows:
    uniq[d.get("original", "")] = d

rows = list(uniq.values())

print("unique URLs =", len(rows), flush=True)

(OUT / "cdx_inventory.json").write_text(
    json.dumps(rows, indent=2) + "\n"
)

# First filename ranking.
interesting = []

for d in rows:
    u = d.get("original", "")
    low = u.lower()

    score = sum(
        1
        for t in [
            "k3", "shimura", "rank", "elliptic",
            "mw", "mordell", "79", "474",
            "highrank", "record"
        ]
        if t in low
    )

    if score:
        interesting.append((score, d))

interesting.sort(
    key=lambda x: (-x[0], x[1].get("original", ""))
)

with (OUT / "interesting_urls.txt").open("w") as f:
    for score, d in interesting:
        f.write(
            f"{score}\t{d.get('timestamp','')}\t"
            f"{d.get('original','')}\n"
        )

print("filename-interesting =", len(interesting), flush=True)

# ------------------------------------------------------------
# Fetch archived textual files.
#
# Focus first on HTML/TXT/PAR/GP/MAGMA/SAGE/C/C++/MAPLE-ish things,
# and only a few thousand candidates.
# ------------------------------------------------------------

text_ext = re.compile(
    r'\.(html?|txt|tex|dat|data|gp|pari|m|mag|magma|sage|py|c|cc|cpp|h|maple)$',
    re.I
)

candidates = []

for d in rows:
    u = d.get("original", "")

    if text_ext.search(u):
        candidates.append(d)
        continue

    low = u.lower()

    if any(
        x in low
        for x in [
            "k3", "shimura", "rank", "elliptic",
            "mw", "mordell", "record"
        ]
    ):
        candidates.append(d)

# Most promising filenames first.
def url_score(d):
    u = d.get("original", "").lower()
    return sum(
        3 if x in u else 0
        for x in ["k3", "rank17", "rk17", "shimura"]
    ) + sum(
        1 if x in u else 0
        for x in ["rank", "elliptic", "mw", "79", "474"]
    )

candidates.sort(key=lambda d: -url_score(d))
candidates = candidates[:5000]

print("content candidates =", len(candidates), flush=True)

hits = []

for ix, d in enumerate(candidates, 1):
    ts = d.get("timestamp", "")
    original = d.get("original", "")

    if not ts or not original:
        continue

    archive_url = (
        "https://web.archive.org/web/"
        + ts
        + "id_/"
        + original
    )

    try:
        raw = get(archive_url)

        # Skip obvious binaries.
        if b"\x00" in raw[:5000]:
            continue

        text = raw.decode("utf-8", errors="ignore")
    except Exception:
        continue

    low = text.lower()

    found = [t for t in TERMS if t in low]

    if found:
        print(
            "HIT",
            found,
            original,
            flush=True
        )

        hits.append({
            "timestamp": ts,
            "original": original,
            "archive_url": archive_url,
            "terms": found,
        })

        safe = re.sub(r'[^A-Za-z0-9_.-]+', '_', original)[-180:]
        (OUT / ("hit-" + safe + ".txt")).write_text(text)

    if ix % 100 == 0:
        print("scanned", ix, flush=True)

    time.sleep(0.03)

(OUT / "content_hits.json").write_text(
    json.dumps(hits, indent=2) + "\n"
)

print("TOTAL HITS =", len(hits))
