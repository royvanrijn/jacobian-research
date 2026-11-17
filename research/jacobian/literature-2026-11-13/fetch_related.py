"""Read-only arXiv metadata query; never executes downloaded material."""
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
QUERY = 'all:"Generalized Vanishing Conjecture"'
url = 'https://export.arxiv.org/api/query?' + urlencode(dict(search_query=QUERY, max_results=100, sortBy='lastUpdatedDate', sortOrder='descending'))
raw = urlopen(url, timeout=45).read()
(HERE / 'related.xml').write_bytes(raw)
ns = {'a': 'http://www.w3.org/2005/Atom'}
rows = []
for entry in ET.fromstring(raw).findall('a:entry', ns):
    row = {key: ' '.join(entry.findtext('a:' + key, '', ns).split()) for key in ('id', 'title', 'published', 'updated', 'summary')}
    row['abstract'] = row.pop('summary')
    row['authors'] = [author.findtext('a:name', '', ns) for author in entry.findall('a:author', ns)]
    rows.append(row)
(HERE / 'related-records.json').write_text(json.dumps(rows, indent=2, ensure_ascii=False) + '\n')
(HERE / 'related-query.json').write_text(json.dumps({'query': QUERY, 'url': url, 'results': len(rows)}, indent=2) + '\n')
for row in rows:
    print(row['id'].rsplit('/', 1)[-1], row['updated'][:10], row['title'])
