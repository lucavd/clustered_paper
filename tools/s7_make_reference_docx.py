#!/usr/bin/env python3
"""Build the Word reference document for manuscript v4 (PI comment of 2026-10-21 in the Word file:
all text Arial and black, figures excluded; body text, not headings, 11 pt).
Starts from the pandoc default reference.docx and rewrites styles.xml and theme1.xml.
Usage: s7_make_reference_docx.py default_reference.docx out_reference.docx"""
import re, sys, zipfile
src, dst = sys.argv[1], sys.argv[2]
HEAD = re.compile(r'w:styleId="(Heading\d|Title|Subtitle|Heading\dChar|TitleChar|SubtitleChar)"')
ARIAL = '<w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="Arial" w:cs="Arial"/>'
def fix_rpr(rpr, body_size):
    rpr = re.sub(r'<w:rFonts[^>]*/>', '', rpr)
    rpr = re.sub(r'<w:color [^>]*/>', '', rpr)
    if body_size:
        rpr = re.sub(r'<w:sz [^>]*/>', '', rpr); rpr = re.sub(r'<w:szCs [^>]*/>', '', rpr)
    inner = rpr[len('<w:rPr>'):-len('</w:rPr>')]
    add = ARIAL + '<w:color w:val="000000"/>' + ('<w:sz w:val="22"/><w:szCs w:val="22"/>' if body_size else '')
    return '<w:rPr>' + add + inner + '</w:rPr>'
def fix_style(m):
    st = m.group(0); is_head = bool(HEAD.search(st)); is_char_table = 'w:type="table"' in st or 'w:type="numbering"' in st
    if is_char_table: return st
    if '<w:rPr>' in st:
        st = re.sub(r'<w:rPr>.*?</w:rPr>', lambda r: fix_rpr(r.group(0), not is_head), st, count=1, flags=re.S)
    elif '<w:rPr/>' in st:
        st = st.replace('<w:rPr/>', fix_rpr('<w:rPr></w:rPr>', not is_head))
    else:
        st = st.replace('</w:style>', fix_rpr('<w:rPr></w:rPr>', not is_head) + '</w:style>')
    return st
zin = zipfile.ZipFile(src); zout = zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED)
for item in zin.infolist():
    data = zin.read(item.filename)
    if item.filename == 'word/styles.xml':
        s = data.decode('utf-8')
        s = re.sub(r'<w:rPrDefault>.*?</w:rPrDefault>', '<w:rPrDefault><w:rPr>' + ARIAL +
                   '<w:color w:val="000000"/><w:sz w:val="22"/><w:szCs w:val="22"/><w:lang w:val="en-GB"/></w:rPr></w:rPrDefault>', s, flags=re.S)
        s = re.sub(r'<w:style [^>]*>.*?</w:style>', fix_style, s, flags=re.S)
        data = s.encode('utf-8')
    elif item.filename == 'word/theme/theme1.xml':
        s = data.decode('utf-8')
        s = re.sub(r'(<a:(?:major|minor)Font>\s*<a:latin typeface=")[^"]*"', r'\1Arial"', s)
        data = s.encode('utf-8')
    zout.writestr(item, data)
zout.close()
