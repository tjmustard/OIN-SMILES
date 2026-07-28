#!/usr/bin/env python3
"""Wrap the retrospective page source into a standalone HTML document for GitHub Pages.

ONE source, TWO destinations. ``docs/agentic-notes/<release>/retrospective.page.html``
holds only the page *body* — no ``<!doctype>``, ``<html>``, ``<head>`` or ``<body>`` tags —
because that is exactly what the Artifact publisher expects to be handed. GitHub Pages needs
a complete document, so this script adds the shell: doctype, meta, a CSS reset, an inline
emoji favicon, and a theme toggle that stamps ``data-theme`` on the root element (the same
attribute the Artifact runtime's own toggle uses, so the page's tokens work unchanged in
both places).

Keeping the wrapper here rather than in the source means the two outputs cannot drift: edit
the body, re-run this, and both the hosted page and the Artifact carry the same content.

Usage
-----
    python3 tools/build_retrospective_page.py \
        --source docs/agentic-notes/v0.4.12/retrospective.page.html \
        --out    /path/to/gh-pages/index.html

Writes nothing else and reads nothing else. No third-party dependencies, deliberately: this
runs from a gh-pages worktree that has no virtualenv.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

# The <head> the body source deliberately does not carry. Everything is inline — GitHub
# Pages serves this from a bare branch, and the Artifact CSP blocks external hosts outright,
# so a linked stylesheet or font would fail silently in one destination and hard in the other.
HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Accuracy and speed across every OIN-SMILES point release since v0.4.4 — including the release where the headline number fell ten points on purpose.">
<meta name="color-scheme" content="light dark">
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'%3E%3Ctext y='.9em' font-size='90'%3E%E2%9A%97%EF%B8%8F%3C/text%3E%3C/svg%3E">
<style>
/* minimal reset - the Artifact wrapper supplies its own; a bare branch does not */
*, *::before, *::after { box-sizing: border-box; }
body { margin: 0; }
img, svg { max-width: 100%; }
table { border-collapse: collapse; }

/* theme toggle - Artifact hosts provide their own chrome for this, GitHub Pages does not */
#theme-toggle {
  position: fixed; top: 14px; right: 14px; z-index: 50;
  width: 38px; height: 38px; border-radius: 50%;
  display: grid; place-items: center; cursor: pointer;
  background: var(--surface, #fff); color: var(--ink-2, #444);
  border: 1px solid var(--rule, #ddd);
  font: 15px/1 ui-sans-serif, system-ui, sans-serif;
  box-shadow: 0 1px 2px rgba(0,0,0,.06), 0 6px 18px -12px rgba(0,0,0,.35);
}
#theme-toggle:hover { color: var(--ink, #000); }
#theme-toggle:focus-visible { outline: 2px solid var(--accent, #008C9E); outline-offset: 3px; }
@media print { #theme-toggle { display: none; } }
</style>
</head>
<body>
<button id="theme-toggle" type="button" aria-label="Toggle light and dark theme" title="Toggle theme">◐</button>
"""

# Reads the stored preference before first paint where possible, then lets the button flip it.
# Stamping `data-theme` is what the page's own token blocks key off, so nothing else changes.
TAIL = """
<script>
(function () {
  var root = document.documentElement;
  var KEY = "oin-retrospective-theme";
  try {
    var saved = localStorage.getItem(KEY);
    if (saved === "light" || saved === "dark") root.setAttribute("data-theme", saved);
  } catch (e) { /* private mode - fall back to prefers-color-scheme */ }

  document.getElementById("theme-toggle").addEventListener("click", function () {
    var explicit = root.getAttribute("data-theme");
    var current = explicit
      || (window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
    var next = current === "dark" ? "light" : "dark";
    root.setAttribute("data-theme", next);
    try { localStorage.setItem(KEY, next); } catch (e) { /* ignore */ }
  });
})();
</script>
</body>
</html>
"""

# Match the OPENING TAG, not the substring. A naive ``"<head" in body`` check fires on
# ``<thead>`` and on a class named ``sec-head`` -- which it did, on the first run, and a guard
# that refuses valid input is worse than no guard because the fix is to delete it.
BANNED_TAG = re.compile(r"<\s*(!doctype|html|head|body)\s*[\s>]", re.IGNORECASE)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    args = ap.parse_args(argv)

    if not args.source.is_file():
        print(f"error: source not found: {args.source}", file=sys.stderr)
        return 2

    body = args.source.read_text(encoding="utf-8")

    # The source is shared with the Artifact publisher, which wraps it in its own document
    # shell. A stray <html> or <body> here would nest documents there and be invisible until
    # someone opened the published page, so refuse rather than emit something subtly broken.
    found = sorted({m.group(1).lower() for m in BANNED_TAG.finditer(body)})
    if found:
        print(
            f"error: {args.source} contains document-level tag(s) {found}.\n"
            "       This file is the page BODY only - it is also handed to the Artifact\n"
            "       publisher, which supplies its own <head>/<body>. Remove them.",
            file=sys.stderr,
        )
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(HEAD + body + TAIL, encoding="utf-8")
    print(f"wrote {args.out} ({args.out.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
