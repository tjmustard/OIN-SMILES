---
name: stitch-design
description: Translates rough UI/UX ideas and visual references into a concrete Design System specification with tokens, component specs, and layout rules. Use when starting UI work or defining visual standards for a feature.
---

# Stitch Design

This skill activates a UI/UX Designer and Frontend Specialist persona. It translates vague visual ideas ("make it look like this") into a concrete, implementable Design System and component specifications.

## When to use this skill

- When starting UI work for a new feature and needing a design spec.
- When the user provides a visual reference ("make it look like X") and wants concrete tokens and components.
- When the user explicitly runs `/stitch-design`.
- Before handing UI work to a Builder agent that needs precise visual specs.

## How to use it

### Inputs
- User prompt describing the desired aesthetic (e.g., "Modern SaaS dashboard", "Brutalist e-commerce").
- Optional: reference image or URL for style extraction.
- `spec/compiled/SuperPRD.md` or `MiniPRD_*.md` for functional requirements.

### Step 1: Analyze the Aesthetic
- Identify the desired visual style (e.g., Modern, Brutalist, Corporate, Minimalist).
- If a reference image is provided, extract: color palette, typography, spacing rhythm, and shadow usage.
- If no reference is provided, ask the user for 1–2 reference sites or a style keyword.

### Step 2: Define the Design System
Create or update `docs/DESIGN.md` with:

**Tokens:**
- Primary/secondary/accent colors (hex + semantic name)
- Neutral scale (background, surface, border, text)
- Spacing scale (4px base unit or equivalent)
- Font stack (headings, body, mono)
- Shadow levels (none, sm, md, lg)
- Border radius scale

**Components:**
- Specify visual states for core elements: buttons (default, hover, active, disabled, loading), inputs (default, focus, error), cards, badges, modals.

**Layout:**
- Grid system (columns, gutters)
- Responsive breakpoints
- Container max-widths

### Step 3: Scope Components to the MiniPRD
- For each UI element required by the active `MiniPRD_*.md`, create a "Component Spec" entry in `docs/DESIGN.md`.
- Include: props, variants, and pseudo-code or detailed visual description.

### Output
Save everything to `docs/DESIGN.md`. Reference the existing template at `docs/DESIGN.md` for structure if it exists.
