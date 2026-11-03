#!/usr/bin/env python3
"""Check a rendered docx against the PI formatting rule of 2026-10-21 (comment in the Word file): all text Arial
and black (figures excluded); text that is not a heading at 11 pt. Resolves, for every paragraph and run in
word/document.xml (and footnotes), the effective font, colour and size through run properties, character style,
paragraph style, basedOn chain and document defaults. Prints the departures; exit 1 if any.
Usage: s7_check_docx_format.py file.docx ..."""
import re, sys, zipfile, collections
HEAD = re.compile(r"(Heading\d|Title|Subtitle)$")
def props(xml):
    f = re.search(r'<w:rFonts [^>]*w:ascii="([^"]+)"', xml); ft = re.search(r'<w:rFonts [^>]*w:asciiTheme="([^"]+)"', xml)
    c = re.search(r'<w:color [^>]*w:val="([^"]+)"', xml); ct = re.search(r'<w:color [^>]*w:themeColor="([^"]+)"', xml)
    s = re.search(r'<w:sz w:val="(\d+)"', xml)
    return {"font": f.group(1) if f else ("theme:" + ft.group(1) if ft else None),
            "color": ("theme:" + ct.group(1)) if ct else (c.group(1).lower() if c else None), "sz": s.group(1) if s else None}
bad = 0
for f in sys.argv[1:]:
    z = zipfile.ZipFile(f); st = z.read("word/styles.xml").decode()
    parts = [z.read(n).decode() for n in ("word/document.xml", "word/footnotes.xml") if n in z.namelist()]
    dd = re.search(r"<w:rPrDefault>(.*?)</w:rPrDefault>", st, re.S); base = props(dd.group(1) if dd else "")
    base = {k: (v if v is not None else {"font": "default", "color": "000000", "sz": "20"}[k]) for k, v in base.items()}
    styles = {}
    for m in re.finditer(r'<w:style [^>]*w:styleId="([^"]+)"[^>]*>(.*?)</w:style>', st, re.S):
        b = re.search(r'<w:basedOn w:val="([^"]+)"', m.group(2)); rp = re.search(r"<w:rPr>(.*?)</w:rPr>", m.group(2), re.S)
        styles[m.group(1)] = (b.group(1) if b else None, props(rp.group(1)) if rp else {"font": None, "color": None, "sz": None})
    def chain(sid, k):
        seen = 0
        while sid and sid in styles and seen < 20:
            v = styles[sid][1][k]
            if v is not None: return v
            sid = styles[sid][0]; seen += 1
        return None
    dep = collections.Counter()
    for xml in parts:
        for p in re.finditer(r"<w:p[ >].*?</w:p>", xml, re.S):
            px = p.group(0); ps = re.search(r'<w:pStyle w:val="([^"]+)"', px); ps = ps.group(1) if ps else "Normal"
            for r in re.finditer(r"<w:r[ >].*?</w:r>", px, re.S):
                rx = r.group(0)
                if "<w:t" not in rx: continue
                rpr = re.search(r"<w:rPr>(.*?)</w:rPr>", rx, re.S); rp = props(rpr.group(1)) if rpr else {"font": None, "color": None, "sz": None}
                cs = re.search(r'<w:rStyle w:val="([^"]+)"', rx); cs = cs.group(1) if cs else None
                eff = {k: rp[k] or (chain(cs, k) if cs else None) or chain(ps, k) or base[k] for k in ("font", "color", "sz")}
                if eff["font"] != "Arial": dep[("font", ps, cs, eff["font"])] += 1
                if eff["color"] not in ("000000", "auto"): dep[("color", ps, cs, eff["color"])] += 1
                if not HEAD.match(ps) and eff["sz"] != "22": dep[("size", ps, cs, eff["sz"])] += 1
    n = sum(dep.values()); bad += n
    print(f, "| text runs departing:", n)
    for k, v in dep.most_common(12): print("   ", k, v)
sys.exit(1 if bad else 0)
