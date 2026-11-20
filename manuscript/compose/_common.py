"""Shared by the compose scripts: panel lookup by symlink_figures.sh name, vector placement, crops of the published
manuscript pages (raster) and of our own panels, grey placeholders for panels not built, the panel-letter overlay, and
compose/order/<figure>.order.yml written from the layout itself."""
import argparse
import io
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore", message=".*fitz.*")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import fitz  # noqa: E402

MANUSCRIPT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# published manuscript (figures are raster images on these pages); MANUSCRIPT_PDF=<pdf> to use another copy
PUBLISHED_PDF = os.environ.get("MANUSCRIPT_PDF", os.path.join(os.path.dirname(MANUSCRIPT), "Watson et al MOSAIC 2026 v4.docx.pdf"))

PT = 72
A4 = (8.27, 11.69)
BASE_RC = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "pdf.fonttype": 42}
LETTER_PT = 9
LETTER_W = 0.17  # inches kept left of a panel for its letter
MIN_SCALE = 0.4  # below this, the panels' 10 pt text prints under 4 pt: warned
PNG_DPI = 300


def cli(fig):
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-dir", default=os.path.join(MANUSCRIPT, "figures", fig), help="linked panels (symlink_figures.sh)")
    ap.add_argument("--out-dir", default=os.path.join(MANUSCRIPT, "composed"))
    ap.add_argument("--final", action="store_true", help="fail on placeholders, missing panels, crops or text below 4 pt")
    return ap.parse_args()


def linked(panel_dir, path, max_len=250):
    """Name symlink_figures.sh gives a BASE_DIR path: '/' -> '_', long names cut to their last 200 characters."""
    name = path.replace("/", "_").lstrip("_")
    if len(name) > max_len:
        stem, ext = os.path.splitext(name)
        name = stem[-200:] + ext
    return os.path.join(panel_dir, name)


def size(path):
    """Page size of a PDF panel, inches."""
    with fitz.open(path) as d:
        r = d[0].rect
    return r.width / PT, r.height / PT


def place(page, path, x, y, w, h, clip=None):
    with fitz.open(path) as src:
        page.show_pdf_page(fitz.Rect(x, y, x + w, y + h) * PT, src, 0, clip=None if clip is None else fitz.Rect(*clip) * PT)


def _clip(item):
    """'crop:<pdf>[@page]#x,y,w,h[!x,y,w,h]' or 'part:<panel>#x,y,w,h' -> (path, page index, (x0, y0, x1, y1) inches)."""
    path, box = item[5:].split("#")
    box = box.split("!")[0]
    pno = 0
    if "@" in path:
        path, p = path.rsplit("@", 1)
        pno = int(p) - 1
    cx, cy, cw, ch = map(float, box.split(","))
    return path, pno, (cx, cy, cx + cw, cy + ch)


def crop(page_no, x, y, w, h, letter=None, pdf=None):
    """Item for a region (inches, top-left origin) of a page of the published manuscript (1-based page number).
    letter: (x, y, w, h), or a list of them, of the published letter(s) to erase when not in the clip's top-left corner."""
    rects = [] if letter is None else [letter] if isinstance(letter[0], (int, float)) else letter
    erase = "".join("!" + ",".join(f"{v:.3f}" for v in r) for r in rects)
    return f"crop:{pdf or PUBLISHED_PDF}@{page_no}#{x:.3f},{y:.3f},{w:.3f},{h:.3f}{erase}"


_SOURCES = {}


def _source(path):
    """Crop sources opened once per run, so a composed PDF holds one copy of each published page image."""
    return _SOURCES.setdefault(path, fitz.open(path))


def erase_letter(item):
    """Paint white the published panel letter of a 'crop:' item (redrawn by the letter overlay): the regions after '!'
    in the item, else the clip's top-left LETTER_W x 0.2 in. Done on the source page before any crop of it is placed."""
    path, pno, (x0, y0, _, _) = _clip(item)
    rects = [fitz.Rect(x0, y0, x0 + LETTER_W, y0 + 0.2)]
    if "!" in item:
        rects = [fitz.Rect(ex, ey, ex + ew, ey + eh) for ex, ey, ew, eh in (map(float, r.split(",")) for r in item.split("!")[1:])]
    for r in rects:
        _source(path)[pno].draw_rect(r * PT, color=None, fill=(1, 1, 1), overlay=True)


def place_crop(page, item, x, y, w, h):
    """item 'crop:<pdf>[@page]#x,y,w,h[!x,y,w,h...]' (inches, top-left origin): that region of another PDF page (vector
    or raster), scaled to fit the box (up or down), centred and top-aligned. Returns the scale and the covered rect."""
    path, pno, (x0, y0, x1, y1) = _clip(item)
    src = _source(path)
    assert fitz.Rect(x0, y0, x1, y1) * PT in src[pno].rect, f"{item}: clip outside the source page"
    s = min(w / (x1 - x0), h / (y1 - y0))
    cw, ch = (x1 - x0) * s, (y1 - y0) * s
    page.show_pdf_page(fitz.Rect(x + (w - cw) / 2, y, x + (w + cw) / 2, y + ch) * PT, src, pno, clip=fitz.Rect(x0, y0, x1, y1) * PT)
    return s, fitz.Rect(x + (w - cw) / 2, y, x + (w + cw) / 2, y + ch)


def placeholder(page, x, y, w, h, text):
    r = fitz.Rect(x * PT, y * PT, (x + w) * PT, (y + h) * PT)
    page.draw_rect(r, color=(0.75, 0.75, 0.75), fill=(0.93, 0.93, 0.93), width=0.5)
    page.insert_textbox(r + (6, 6, -6, -6), text, fontsize=7, fontname="helv", color=(0.35, 0.35, 0.35),
                        align=fitz.TEXT_ALIGN_CENTER)


def overlay_letters(page, letters):
    """letters: [(letter, x, y)] in inches from the top-left; sans-serif bold LETTER_PT."""
    w, h = page.rect.width / PT, page.rect.height / PT
    plt.rcParams.update(BASE_RC)
    ov = plt.figure(figsize=(w, h))
    ov.patch.set_alpha(0)
    for letter, x, y in letters:
        ov.text(x / w, 1 - y / h, letter, fontsize=LETTER_PT, fontweight="bold", ha="left", va="top")
    buf = io.BytesIO()
    ov.savefig(buf, format="pdf", transparent=True)
    plt.close(ov)
    with fitz.open("pdf", buf.getvalue()) as ovd:
        page.show_pdf_page(page.rect, ovd, 0, overlay=True)


def write_order(fig, order, out_dir=None):
    """compose/order/<fig>.order.yml: {letter: [panel path, 'crop: ...' or 'pending: ...']} in compose order.
    A non-default --out-dir (scratch composes) gets <out_dir>/<fig>.order.yml instead, so the repo file is not touched."""
    default = os.path.join(MANUSCRIPT, "composed")
    if out_dir and os.path.abspath(out_dir) != default:
        os.makedirs(out_dir, exist_ok=True)
        out = os.path.join(out_dir, f"{fig}.order.yml")
    else:
        out = os.path.join(MANUSCRIPT, "compose", "order", f"{fig}.order.yml")
    with open(out, "w") as f:
        f.write(f"# written by compose/{fig}.py from its layout; do not edit (rerun the compose script)\n")
        f.writelines(f"{k}: {json.dumps(v)}\n" for k, v in order.items())  # JSON lists are YAML flow sequences
    return out


def save(doc, out_dir, name):
    """<name>.pdf (all pages) and <name>.png, or <name>_p<N>.png per page for a multi-page figure."""
    os.makedirs(out_dir, exist_ok=True)
    pdf = os.path.join(out_dir, f"{name}.pdf")
    doc.save(pdf, deflate=True, garbage=3)
    for i, page in enumerate(doc):
        png = f"{name}.png" if len(doc) == 1 else f"{name}_p{i + 1}.png"
        page.get_pixmap(dpi=PNG_DPI).save(os.path.join(out_dir, png))
    print(f"Saved: {pdf} (+ .png, {len(doc)} page{'s' if len(doc) > 1 else ''})")


def check_inside(layout, page_wh):
    """Every box inside the page."""
    page = fitz.Rect(0, 0, *page_wh)
    outside = [letter for letter, (x, y, w, h), _ in layout if fitz.Rect(x, y, x + w, y + h) not in page]
    assert not outside, f"boxes outside the {page_wh[0]} x {page_wh[1]} in page: {outside}"


def check_no_overlap(drawn):
    """drawn: [(letter, fitz.Rect inches)] areas actually covered by panels (boxes may overlap where panels are narrower)."""
    overlaps = [(a, b) for i, (a, rect_a) in enumerate(drawn) for b, rect_b in drawn[i + 1:]
                if (rect_a & rect_b).width > 0.01 and (rect_a & rect_b).height > 0.01]
    assert not overlaps, f"overlapping panels: {overlaps}"


def grid(letter, box, items, ncol, gap=0.05):
    """Layout entries for one lettered panel made of several files, row-major in ncol columns inside box.
    The first cell is widened by LETTER_W on the left so every cell has the same panel width."""
    x, y, w, h = box
    nrow = -(-len(items) // ncol)
    cw = (w - LETTER_W - gap * (ncol - 1)) / ncol
    ch = (h - gap * (nrow - 1)) / nrow
    out = []
    for i, item in enumerate(items):
        r, c = divmod(i, ncol)
        cx = x + LETTER_W + c * (cw + gap)
        out.append((letter, (cx - LETTER_W, y + r * (ch + gap), cw + LETTER_W, ch) if i == 0 else (cx, y + r * (ch + gap), cw, ch), item))
    return out


class Published:
    """Maps inches on a published figure page to inches on our A4 page: the region `content` (x0, y0, x1, y1) of
    manuscript page `page_no` is scaled to the page width minus `margin` on each side and moved to `top`."""

    def __init__(self, page_no, content, margin=0.2, top=0.25):
        self.page_no, self.content, self.top, self.margin = page_no, content, top, margin
        self.k = (A4[0] - 2 * margin) / (content[2] - content[0])

    def box(self, x, y, w, h):
        return (self.margin + (x - self.content[0]) * self.k, self.top + (y - self.content[1]) * self.k, w * self.k, h * self.k)

    def crop(self, x, y, w, h, letter=None):
        return crop(self.page_no, x, y, w, h, letter)


def _aspect(panel_dir, item):
    if item.startswith(("part:", "crop:")):
        _, _, (x0, y0, x1, y1) = _clip(item)
        return (x1 - x0) / (y1 - y0)
    src = linked(panel_dir, item)
    if os.path.exists(src):
        w, h = size(src)
        return w / h
    return 1.5


def flow(rows, panel_dir, top=0.25, bottom=11.45, left=0.2, right=8.07, gap=0.12, max_h=4.0):
    """Supplementary-figure pages from rows of (letter, item): each row gets the height at which its panels fill the
    width between left and right (at most max_h, or the row's own cap given as a third element (letter, item, cap) of
    its first entry), centred; a row that does not fit, or a None row, starts a new page. A letter column is kept for
    the first box of each letter. Returns [layout per page]."""
    pages, layout, y, seen = [], [], top, set()
    for row in rows:
        if row is None:
            pages.append(layout)
            layout, y = [], top
            continue
        cap = row[0][2] if len(row[0]) > 2 else max_h
        aspects = [_aspect(panel_dir, e[1]) for e in row]
        lws = []
        for e in row:
            lws.append(LETTER_W if e[0] and e[0] not in seen else 0)
            seen.add(e[0])
        avail = right - left - gap * (len(row) - 1) - sum(lws)
        h = min(cap, avail / sum(aspects))
        if y + h > bottom and layout:
            pages.append(layout)
            layout, y = [], top
        x = left + (avail - h * sum(aspects)) / 2
        for e, a, lw in zip(row, aspects, lws):
            layout.append((e[0], (x, y, h * a + lw, h), e[1]))
            x += h * a + lw + gap
        y += h + gap
    pages.append(layout)
    return pages


def compose(fig, layout, args, page_wh=A4, extra_letters=()):
    """layout: [(letter, (x, y, w, h) inches, item)] or a list of such layouts (one per page). item = a BASE_DIR path
    (as in symlink_config.yml), 'part:<BASE_DIR path>#x,y,w,h' (that region of our panel, vector), 'crop:...' (crop(), a
    region of a published page) or 'pending: <what is needed>'. The letter goes at the box's top-left (first box of each
    letter only; later boxes of the same letter are continuations and use their full width). A letter ending in '*'
    (e.g. 'c*') keeps the published letter inside its crop: not erased, not redrawn. Panels are scaled down, never up,
    to fit the box, centred horizontally and top-aligned; crops are scaled to the box. A panel not linked in
    --panel-dir becomes a placeholder (logged). extra_letters: [(letter, x, y)] drawn on the first page."""
    pages = layout if layout and isinstance(layout[0], list) else [layout]
    seen = set()
    for letter, _, item in (e for lay in pages for e in lay):
        if item.startswith("crop:") and letter and not letter.endswith("*") and letter not in seen:
            erase_letter(item)
        seen.add(letter.rstrip("*"))
    doc = fitz.open()
    order, log = {}, []
    for pno, lay in enumerate(pages):
        check_inside(lay, page_wh)
        page = doc.new_page(width=page_wh[0] * PT, height=page_wh[1] * PT)
        letters, drawn = [], []
        if len(pages) > 1:
            log.append(f" page {pno + 1}")
        for letter, (x, y, w, h), item in lay:
            keep = letter.endswith("*")
            letter = letter.rstrip("*")
            first = bool(letter) and letter not in order
            kind = item.split(":", 1)[0] if item.startswith(("crop:", "part:", "pending:")) else "panel"
            if letter:
                order.setdefault(letter, []).append(f"crop: page {_clip(item)[1] + 1} of the manuscript PDF" if kind == "crop" else item)
            if first and not keep:
                letters.append((letter, x, y))
            if kind == "crop":
                s, r = place_crop(page, item, x, y, w, h)
                drawn.append((letter, r))
                log.append(f"  {letter}: crop page {_clip(item)[1] + 1} #{item.split('#')[1].split('!')[0]} scale {s:.2f}")
                continue
            if first and not keep:
                x, w = x + LETTER_W, w - LETTER_W
            if kind == "pending":
                log.append(f"  {letter}: placeholder, {item}")
                placeholder(page, x, y, w, h, f"{letter}  {item}")
                drawn.append((letter, fitz.Rect(x, y, x + w, y + h)))
                continue
            path, clip = item, None
            if kind == "part":
                path, _, clip = _clip(item)
            src = linked(args.panel_dir, path)
            if not os.path.exists(src):
                log.append(f"  {letter}: MISSING {path} (run its analysis, then symlink_figures.sh {fig})")
                placeholder(page, x, y, w, h, f"{letter}  missing: {path}")
                drawn.append((letter, fitz.Rect(x, y, x + w, y + h)))
                continue
            pw, ph = (clip[2] - clip[0], clip[3] - clip[1]) if clip else size(src)
            s = min(1.0, w / pw, h / ph)
            place(page, src, x + (w - pw * s) / 2, y, pw * s, ph * s, clip)
            drawn.append((letter, fitz.Rect(x + (w - pw * s) / 2, y, x + (w + pw * s) / 2, y + ph * s)))
            log.append(f"  {letter}: {item} scale {s:.2f}" + ("  [WARN: text < 4 pt]" if s < MIN_SCALE else ""))
        check_no_overlap(drawn)
        overlay_letters(page, letters + (list(extra_letters) if pno == 0 else []))
    print(f"{fig}:\n" + "\n".join(log))
    not_final = [line for line in log if any(tag in line for tag in ("MISSING", "placeholder", "crop page", "WARN"))]
    if args.final and not_final:
        sys.exit(f"{fig} is not final (--final):\n" + "\n".join(not_final))
    print("order ->", write_order(fig, order, args.out_dir))
    save(doc, args.out_dir, fig)
    doc.close()
