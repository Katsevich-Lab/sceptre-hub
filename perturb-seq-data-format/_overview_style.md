```{=html}
<style>
.tree { border: 1px solid rgba(120, 130, 150, 0.35); border-radius: 8px; padding: 0.5em 0.7em 0.7em; margin: 0.8em 0 1.4em; font-size: 0.93em; }
.tree .bar { display: flex; gap: 0.5em; justify-content: flex-end; margin-bottom: 0.2em; }
.tree .bar button { font-size: 0.8em; padding: 0.1em 0.6em; border: 1px solid rgba(120, 130, 150, 0.45); border-radius: 4px; background: none; color: inherit; cursor: pointer; }
.tree .bar button:hover { background: rgba(120, 130, 150, 0.15); }
.tree ul { list-style: none; margin: 0; padding-left: 0; }
.tree ul ul { margin-left: 0.55em; padding-left: 0.85em; border-left: 1px solid rgba(120, 130, 150, 0.35); }
.tree summary.row { display: grid; grid-template-columns: calc(12.5em - var(--d) * 1.4em) 1fr auto; gap: 0 1em; align-items: baseline;
  padding: 0.22em 0.4em; border-radius: 4px; cursor: pointer; list-style: none; }
.tree summary.row::-webkit-details-marker { display: none; }
.tree summary.row:hover { background: rgba(120, 130, 150, 0.14); }
.tree .name::before { content: "▸"; display: inline-block; width: 1.1em; opacity: 0.6; }
.tree details[open] > summary .name::before { content: "▾"; }
.tree .name code { background: none; padding: 0; color: inherit; font-size: 0.95em; }
.tree .folder .name code { font-weight: 700; }
.tree .what { opacity: 0.92; }
.tree .size { opacity: 0.6; font-size: 0.88em; white-space: nowrap; font-variant-numeric: tabular-nums; }
.tree .body { margin: 0.3em 0 0.9em 1.5em; padding: 0.5em 0.7em; overflow-x: auto; font-size: 0.88em; background: rgba(120, 130, 150, 0.08); border-radius: 6px; }
.tree .body pre { margin: 0; padding: 0; background: none; border: none; font-size: 0.95em; }
.tree table.rows { width: auto; margin: 0.2em 0; border-collapse: collapse; font-variant-numeric: tabular-nums; }
.tree table.rows th, .tree table.rows td { padding: 0.15em 0.8em 0.15em 0; text-align: left; border-bottom: 1px solid rgba(120, 130, 150, 0.25); white-space: nowrap; }
.tree table.rows th { font-weight: 600; }
.tree .more, .tree ul.columns { margin: 0.2em 0; opacity: 0.8; }
.tree ul.columns { border: none; padding-left: 1.2em; list-style: disc; margin-bottom: 0.5em; }
@media (max-width: 700px) { .tree summary.row { grid-template-columns: 1fr auto; } .tree .what { grid-column: 1 / -1; grid-row: 2; padding-left: 1.1em; } }
details.fold { border-bottom: 1px solid rgba(120, 130, 150, 0.3); padding: 0.1em 0; }
details.fold > summary { cursor: pointer; padding: 0.45em 0.2em; list-style: none; }
details.fold > summary::-webkit-details-marker { display: none; }
details.fold > summary::before { content: "▸"; display: inline-block; width: 1.2em; opacity: 0.55; }
details.fold[open] > summary::before { content: "▾"; }
details.fold > summary:hover { background: rgba(120, 130, 150, 0.12); border-radius: 4px; }
details.fold > :not(summary) { margin-left: 1.4em; }
details.fold .where { opacity: 0.6; font-size: 0.88em; white-space: nowrap; }
table.status td, table.status th { padding: 0.35em 0.8em 0.35em 0; vertical-align: top; text-align: left; border-bottom: 1px solid rgba(120, 130, 150, 0.25); }
table.status td:nth-child(n+3) { white-space: nowrap; font-variant-numeric: tabular-nums; }
table.defs { font-size: 0.93em; }
table.defs td, table.defs th { padding: 0.3em 0.9em 0.3em 0; vertical-align: top; text-align: left; border-bottom: 1px solid rgba(120, 130, 150, 0.25); }
table.defs td:nth-child(2) { white-space: nowrap; opacity: 0.75; }
pre.layout { font-size: 0.88em; line-height: 1.35; }
</style>
<script>
function treeAll(button, open) {
  const tree = button.closest(".tree");
  // folders only: a file's contents open on a click of its own
  tree.querySelectorAll("details").forEach(d => { if (d.querySelector(":scope > ul")) d.open = open; else if (!open) d.open = false; });
  if (!open) tree.querySelector("details").open = true;
}
</script>
```
