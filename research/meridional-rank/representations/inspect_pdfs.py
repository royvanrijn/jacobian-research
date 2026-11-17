#!/usr/bin/env python3
"""Read text and assemble Poppler previews; needs bundled pypdf and Pillow."""

import json
import logging
import sys
from pathlib import Path
from pypdf import PdfReader
from PIL import Image, ImageOps, ImageDraw

HERE = Path(__file__).resolve().parent


def main():
    logging.getLogger("pypdf").setLevel(logging.ERROR)
    report = {}
    for title in (sys.argv[1:] or ["PROOF", "LITERATURE"]):
        path = HERE / "output/pdf" / (title + ".pdf")
        reader = PdfReader(path)
        text = "\n".join(page.extract_text() for page in reader.pages)
        (HERE / "research/.scratch/pdfs" / (title + ".txt")).write_text(text, encoding="utf-8")
        previews = [p for p in (HERE / "research/.scratch/pdfs").glob(title + "-*.png")
                    if p.stem[len(title) + 1:].isdigit()]
        assert len(previews) == len(reader.pages)
        report[title] = {"pages": len(reader.pages), "preview_count": len(previews)}
        contact = Image.new("RGB", (1200, 570 * ((len(previews) + 2) // 3)), "#dddddd")
        draw = ImageDraw.Draw(contact)
        for i, preview in enumerate(sorted(previews)):
            with Image.open(preview) as original:
                thumb = ImageOps.contain(original.convert("RGB"), (390, 535))
            x, y = (i % 3) * 400, (i // 3) * 570
            contact.paste(thumb, (x, y + 22))
            draw.text((x + 8, y + 4), preview.name, fill="black")
        contact.save(HERE / "research/.scratch/pdfs" / (title + "-contact.png"))
    source = PdfReader(HERE / "sources/bbk21.pdf")
    paragraphs = []
    for i, page in enumerate(source.pages):
        text = page.extract_text()
        if "pretzel" in text.lower() and "cycle" in text.lower():
            paragraphs.append({"pdf_page": i + 1, "excerpt": text[:1900]})
    print(json.dumps({"outputs": report, "primary_pretzel_quotient_passages": paragraphs}, indent=2))


if __name__ == "__main__":
    main()
