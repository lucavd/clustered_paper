
import zipfile, re, shutil, sys
p = sys.argv[1]; tmp = p + ".tmp"
zin = zipfile.ZipFile(p); zout = zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED)
SZ = {"Title": "36", "Heading1": "32", "Heading2": "28", "Heading3": "24"}
for it in zin.infolist():
    b = zin.read(it.filename)
    if it.filename == "word/styles.xml":
        s = b.decode()
        for sid, v in SZ.items():
            m = re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>' % sid, s, re.S); assert m, sid
            blk = m.group(0)
            new = re.sub(r'<w:sz w:val="\d+"\s*/>', '<w:sz w:val="%s" />' % v, blk)
            new = re.sub(r'<w:szCs w:val="\d+"\s*/>', '<w:szCs w:val="%s" />' % v, new)
            if sid == "Title": new = re.sub(r'<w:jc w:val="center"\s*/>', '<w:jc w:val="left" />', new)
            s = s.replace(blk, new)
        b = s.encode()
    zout.writestr(it, b)
zout.close(); shutil.move(tmp, p)
s = zipfile.ZipFile(p).read("word/styles.xml").decode()
for sid in SZ:
    blk = re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>' % sid, s, re.S).group(0)
    print(sid, re.findall(r'<w:sz w:val="(\d+)"', blk), re.findall(r'<w:jc w:val="(\w+)"', blk))
