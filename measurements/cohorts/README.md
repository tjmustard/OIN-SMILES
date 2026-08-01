# Cohort manifests — the sample definitions every published rate rests on

**These were untracked until v0.4.16.** `tmCAT-tmPHOTO_xyz_dataset/` is gitignored in its entirety,
so the file naming which 5000 molecules the corpus figure is computed over lived one `rm -rf` — or
one `Kulik_TMC_Dataset` branch switch, which silently deletes 26,232 files while `git status` stays
clean — from oblivion.

This project already applies the rule *"a rate without its sample is not reproducible"* to lane
populations (`measurements/v0.4.15/pop_*.txt`). It had not applied it to the **corpus**, which is
the sample under every headline it has ever published.

| manifest | n | used by |
|---|---:|---|
| `cohort-v0.4.5-5k_manifest.json` | **5000** | **the corpus of record** — v0.4.6, v0.4.8-honest, v0.4.14 sweeps, every lane A/B since |
| `cohort-v047-slow100_manifest.json` | 100 | v0.4.7 slow-molecule profiling |
| `cohort-v048-live300_manifest.json` | 300 | v0.4.8's live 300-molecule confirmation |
| `cohort-v049-strata_manifest.json` | 328 | the frozen runtime benchmark (reproduces to 0.28%) |

Each records `seed`, `n`, `commit_id`, `built_at`, the dedup priority between the `cat/` and
`photo/` subdirectories, and **every molecule name**. The dedup priority is load-bearing:
**1033 basenames exist in BOTH subdirs** and the harness keys reports by basename, so a cohort
rebuilt with a different priority is a different cohort that looks identical.

## Rebuilding vs restoring

`tools/build_sweep_cohort.py` rebuilds a cohort from the dataset — but it needs the dataset, and it
reproduces *the same draw* only at the same seed and dedup priority. **These manifests are the
authority on what was actually drawn**, so prefer restoring the symlink farm from a manifest over
re-drawing.

⚠ Local paths are scrubbed to `<HOME>/` or made repo-relative — `measurements/` is public. That is
the harvester's own `scrub()`, applied here because these were copied in by hand rather than
harvested, and the hand path skips the guard.

## Checking a cohort is intact

```bash
find tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k -xtype l | wc -l   # MUST be 0
ls tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k | wc -l              # MUST equal the manifest's n
```

A dangling symlink is the silent failure mode: the sweep runs, produces a plausible short table, and
exits 0.
