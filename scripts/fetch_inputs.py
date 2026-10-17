"""Download the pinned original inputs and verify their SHA-256 digests."""
from pathlib import Path
import hashlib, json, urllib.request
root = Path(__file__).resolve().parents[1]
manifest = json.loads((root/'inputs/source-manifest.json').read_text())
out = root/'inputs/original'
out.mkdir(parents=True, exist_ok=True)
for row in manifest['files']:
    name = row['name']
    assert Path(name).name == name and '/' not in name and '\\' not in name
    target = out/name
    if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() == row['sha256']:
        continue
    with urllib.request.urlopen(row['url'], timeout=60) as response:
        data = response.read()
    if hashlib.sha256(data).hexdigest() != row['sha256']:
        raise RuntimeError('Source hash mismatch: '+name)
    target.write_bytes(data)
print('Verified', len(manifest['files']), 'original input files.')
