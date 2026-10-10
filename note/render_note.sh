#!/bin/bash
# Markdown -> HTML (pandoc) -> PDF (headless Chromium); run from note/
set -e
pandoc thermal_lever_note.md -f markdown -t html5 -s --css note.css --metadata pagetitle="The thermal lever against the duplicate-cell floor" -o thermal_lever_note.html
python3 - <<'PY'
import re
p = 'thermal_lever_note.html'
h = open(p).read()
h = re.sub(r'<header id="title-block-header">[\s\S]*?</header>', '', h)
unit = r'(?:mOhm|mW/K|W/mK|W/K|kJ/mol|mm|UTC|C|K|s|W|%)(?![\w/])'
parts = re.split(r'(<[^>]+>)', h)
for i, t in enumerate(parts):
    if not t.startswith('<'):
        parts[i] = re.sub(r'(\d) (?=' + unit + ')', '\\1\u00a0', t)
        parts[i] = re.sub(r'\b([hfsc]) = ', '\\1\u00a0=\u00a0', parts[i])
open(p, 'w').write(''.join(parts))
PY
CHROME=${CHROME:-$(ls -d /opt/pw-browsers/chromium-*/chrome-linux/chrome 2>/dev/null | head -1)}
"$CHROME" --headless --no-sandbox --disable-gpu --no-pdf-header-footer --print-to-pdf=thermal_lever_note.pdf --allow-file-access-from-files "file://$(pwd)/thermal_lever_note.html" 2>/dev/null
pdfinfo thermal_lever_note.pdf | grep -E 'Title|Pages'
