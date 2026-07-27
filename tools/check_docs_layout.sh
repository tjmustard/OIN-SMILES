#!/usr/bin/env bash
# Guard the docs/ root against agentic session notes.
#
# docs/ is split by audience: the root carries a small, fixed set of product docs; every
# session artifact (measurement reports, lane write-ups, A/B results, refuted hypotheses,
# status snapshots) belongs under docs/agentic-notes/<release>/.
#
# Invoked from .git/hooks/pre-commit. The logic lives here, tracked, so the rule survives a
# fresh clone even though .git/hooks/ does not.
#
# Full rule: .agents/rules/docs-layout.md
set -euo pipefail

# The complete product-doc allowlist. Adding a fifth entry requires maintainer sign-off and
# a matching update to .agents/rules/docs-layout.md, AGENTS.md and docs/README.md.
ALLOWED=(
    "docs/README.md"
    "docs/OPTIMIZERS.md"
    "docs/GENERATION_PIPELINE.md"
    "docs/KNOWN_LIMITATIONS.md"
)

# Only ADDED files are checked. Editing an existing root doc is always fine, and this keeps
# the guard from blocking a rename that moves a file *out* of the root.
added=$(git diff --cached --name-only --diff-filter=A || true)

violations=()
for file in $added; do
    # docs/ root only — one path segment after "docs/". Subdirectories are unrestricted.
    case "$file" in
        docs/*/*) continue ;;
        docs/*) ;;
        *) continue ;;
    esac

    allowed=0
    for ok in "${ALLOWED[@]}"; do
        [[ $file == "$ok" ]] && allowed=1 && break
    done
    [[ $allowed -eq 1 ]] || violations+=("$file")
done

if [[ ${#violations[@]} -gt 0 ]]; then
    echo "❌ COMMIT BLOCKED — new file(s) at the docs/ root:"
    for v in "${violations[@]}"; do
        echo "  -> $v"
    done
    echo ""
    echo "The docs/ root is closed. It holds four product docs and nothing else:"
    printf '  %s\n' "${ALLOWED[@]}"
    echo ""
    echo "Session notes — anything recording what you measured, tried, or refuted — go in"
    echo "docs/agentic-notes/<release>/, where <release> is the release the work is FOR."
    echo ""
    echo "  git restore --staged ${violations[0]}"
    echo "  mkdir -p docs/agentic-notes/<release>"
    echo "  git mv ${violations[0]} docs/agentic-notes/<release>/"
    echo ""
    echo "If this really is product documentation, get maintainer sign-off and add it to"
    echo "ALLOWED in tools/check_docs_layout.sh in the same commit, with a reason."
    echo "Do NOT use --no-verify: it is banned for routine commits and also skips the"
    echo "commit-msg trailer rewrite. See .agents/rules/docs-layout.md"
    exit 1
fi

exit 0
