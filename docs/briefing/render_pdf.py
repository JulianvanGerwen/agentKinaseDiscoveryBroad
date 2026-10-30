"""Render BRIEFING.md to BRIEFING.html and BRIEFING.pdf.

Usage (from the repository root):
    pip install markdown xhtml2pdf
    python docs/briefing/render_pdf.py
"""
import re
from pathlib import Path

import markdown
from xhtml2pdf import pisa

HERE = Path(__file__).resolve().parent
SRC = HERE / "BRIEFING.md"

CSS = """
@font-face { font-family: DejaVu; src: url("fonts/DejaVuSans.ttf"); }
@font-face { font-family: DejaVu; src: url("fonts/DejaVuSans-Bold.ttf"); font-weight: bold; }
@font-face { font-family: DejaVu; src: url("fonts/DejaVuSans-Oblique.ttf"); font-style: italic; }
@font-face { font-family: DejaVu; src: url("fonts/DejaVuSans-BoldOblique.ttf"); font-weight: bold; font-style: italic; }
@page { size: A4; margin: 2cm 1.9cm 2cm 1.9cm;
        @frame footer { -pdf-frame-content: footer; bottom: 1cm; margin-left: 1.9cm; margin-right: 1.9cm; height: 0.6cm; } }
body { font-family: DejaVu; font-size: 9.5pt; line-height: 1.4; color: #222; }
h1 { font-size: 19pt; color: #1b4f72; margin-bottom: 2pt; }
h2 { font-size: 14pt; color: #1b4f72; border-bottom: 1px solid #aab7c4; padding-bottom: 2pt; margin-top: 16pt; }
h3 { font-size: 11.5pt; color: #1b4f72; margin-top: 12pt; }
p { margin: 4pt 0 6pt 0; text-align: justify; }
li { margin-bottom: 3pt; }
em { color: #333; }

table { border: 0.5px solid #999; margin: 6pt 0 8pt 0; }
th { background-color: #e8eef4; font-size: 8.5pt; padding: 3pt; border: 0.5px solid #999; text-align: left; }
td { font-size: 8.5pt; padding: 3pt; border: 0.5px solid #ccc; vertical-align: top; }
hr { color: #aab7c4; }
code { font-family: Courier; font-size: 8.5pt; }
#footer { font-size: 8pt; color: #777; text-align: center; }
"""


LIST_RE = re.compile(r"^(\s*)([-*]|\d+\.)\s")


def fix_lists(text):
    """GitHub-style markdown -> Python-Markdown: blank line before lists, 4-space nesting."""
    out, prev = [], ""
    for line in text.split("\n"):
        m = LIST_RE.match(line)
        if m:
            indent = len(m.group(1))
            line = " " * (indent * 2) + line.lstrip()  # 2-space nesting -> 4-space
            if prev.strip() and not LIST_RE.match(prev) and not prev.startswith(" "):
                out.append("")
        out.append(line)
        prev = line
    return "\n".join(out)


def wrap_figures(text):
    """Keep each figure image and its caption together on one page."""
    pat = re.compile(r"!\[[^\]]*\]\(([^)]+)\)\n\n(\*Figure .*?\*)\n", re.S)

    def rep(m):
        cap = markdown.markdown(m.group(2))
        return (f'<div style="page-break-inside: avoid; margin: 8pt 0;">'
                f'<img src="{(HERE / m.group(1)).as_posix()}" style="width: 17cm;"/>'
                f'<div style="font-size: 8.5pt; color: #333;">{cap}</div></div>\n')
    return pat.sub(rep, text)


def main():
    text = SRC.read_text(encoding="utf-8")
    text = re.sub(r"\A---\n.*?\n---\n", "", text, flags=re.S)  # drop YAML front matter
    text = fix_lists(text)
    text = wrap_figures(text)
    body = markdown.markdown(text, extensions=["tables", "sane_lists", "md_in_html"])
    html = (
        "<html><head><meta charset='utf-8'><style>" + CSS + "</style></head><body>"
        + body
        + "<div id='footer'>AI-generated research briefing; not verified by a human. Page <pdf:pagenumber></div>"
        + "</body></html>"
    )
    (HERE / "BRIEFING.html").write_text(html.replace(HERE.as_posix() + "/", ""), encoding="utf-8")
    with open(HERE / "BRIEFING.pdf", "wb") as fh:
        res = pisa.CreatePDF(html, dest=fh, path=str(HERE) + "/")
    print("errors:", res.err)


if __name__ == "__main__":
    main()
