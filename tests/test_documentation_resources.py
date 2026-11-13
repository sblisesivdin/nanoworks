# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Standard-library checks for local Sphinx documentation resources.

Run with:
    python -m unittest tests/test_documentation_resources.py

These tests do not build Sphinx or execute simulation engines. They check
that locally referenced documents, images and downloadable files exist.
External web links and generated HTML anchor IDs remain out of scope.
"""

from pathlib import Path
from typing import Optional
import re
import unittest


_SOURCE = Path(__file__).resolve().parents[1] / "docs" / "source"
_RST = tuple(sorted(_SOURCE.rglob("*.rst")))

_DOC = re.compile(r":doc:`(?:[^`<>]*<)?([^`<>]+)>?`", re.DOTALL)
_DOWNLOAD = re.compile(r":download:`[^`]*?<([^<>]+)>`", re.DOTALL)
_IMAGE = re.compile(r"^\s*\.\. (?:figure|image)::\s+(\S+)", re.MULTILINE)
_HTML_PAGE = re.compile(r'\bhref="([^"]+\.html(?:#[^"]*)?)"')


def _resolve(page: Path, target: str) -> Optional[Path]:
    """Return a local resource path; leave external and special refs alone."""
    target = target.strip()
    if (not target or target == "self" or "://" in target
            or target.startswith(("#", "mailto:", "/"))):
        return None
    return (page.parent / target).resolve()


def _documents():
    for page in _RST:
        yield page, page.read_text(encoding="utf-8")


def _toc_entries(text: str):
    """Read ordinary indented entries following a toctree directive."""
    inside = False
    for line in text.splitlines():
        if re.match(r"^\.\. toctree::", line):
            inside = True
            continue
        if not inside:
            continue
        if line.strip() and not line.startswith(("   ", "\t")):
            inside = False
            continue
        entry = line.strip()
        if entry and not entry.startswith(":"):
            yield entry


class DocumentationResourceTests(unittest.TestCase):
    def test_local_document_roles(self):
        for page, text in _documents():
            for match in _DOC.finditer(text):
                target = match.group(1).strip().rstrip(">")
                if target.startswith("~"):
                    target = target.lstrip("~")
                resolved = _resolve(page, target)
                if resolved is None:
                    continue
                with self.subTest(page=page.relative_to(_SOURCE), target=target):
                    self.assertTrue(
                        resolved.with_suffix(".rst").is_file()
                        or resolved.with_suffix(".md").is_file(),
                        f"Missing local :doc: target: {target}",
                    )

    def test_local_figures_and_images(self):
        for page, text in _documents():
            for match in _IMAGE.finditer(text):
                target = match.group(1)
                resolved = _resolve(page, target)
                if resolved is None or "*" in target:
                    continue
                with self.subTest(page=page.relative_to(_SOURCE), target=target):
                    self.assertTrue(resolved.is_file(), f"Missing image: {target}")

    def test_downloadable_reference_files(self):
        for page, text in _documents():
            for match in _DOWNLOAD.finditer(text):
                target = match.group(1)
                resolved = _resolve(page, target)
                if resolved is None:
                    continue
                with self.subTest(page=page.relative_to(_SOURCE), target=target):
                    self.assertTrue(resolved.is_file(), f"Missing download: {target}")

    def test_toctree_targets(self):
        for page, text in _documents():
            for target in _toc_entries(text):
                # Sphinx can use glob patterns, external references or captions.
                if target.startswith(("http://", "https://", ":")) or "*" in target:
                    continue
                if "<" in target and target.endswith(">"):
                    target = target.rsplit("<", 1)[1][:-1]
                resolved = _resolve(page, target)
                if resolved is None:
                    continue
                with self.subTest(page=page.relative_to(_SOURCE), target=target):
                    self.assertTrue(
                        resolved.with_suffix(".rst").is_file()
                        or resolved.with_suffix(".md").is_file(),
                        f"Missing toctree target: {target}",
                    )

    def test_embedded_html_document_links(self):
        # Only check page existence, not Sphinx-generated fragment identifiers.
        for page, text in _documents():
            for match in _HTML_PAGE.finditer(text):
                target = match.group(1).split("#", 1)[0]
                resolved = _resolve(page, target)
                if resolved is None:
                    continue
                with self.subTest(page=page.relative_to(_SOURCE), target=target):
                    self.assertTrue(
                        resolved.with_suffix(".rst").is_file()
                        or resolved.with_suffix(".md").is_file(),
                        f"Missing local HTML page source: {target}",
                    )


if __name__ == "__main__":
    unittest.main()
