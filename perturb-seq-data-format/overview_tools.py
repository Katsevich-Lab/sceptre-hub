"""What an overview page is built from: the tree of a data file, the table of check
results and the list of principles. The page sets SCREEN and screen_spec before
calling these."""
import html
import json
import re
from pathlib import Path

import anndata as ad
import h5py
import pandas as pd

# Set by the page that uses this module: the example file it shows, and the
# specification's own file, from which every description on the page is read.
SCREEN = screen_spec = None
ROWS = 5
pd.set_option("display.width", 120)
pd.set_option("display.max_columns", 12)


def text(value) -> str:
    """Plain text as HTML, with `code` spans."""
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", html.escape(str(value)))


def n(count, noun) -> str:
    return f"{count:,} {noun}{'' if count == 1 else 's'}"


def cell(value) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    if hasattr(value, "__len__") and not isinstance(value, str):
        shown = ", ".join(f"{v:.3g}" for v in list(value)[:4])
        return f"[{shown}{', …' if len(value) > 4 else ''}]"
    return str(value)


def rows(frame: pd.DataFrame, index_name=None) -> str:
    """The first rows of a table."""
    head = frame.head(ROWS)
    if index_name is not None:
        head = head.rename_axis(index_name).reset_index()
    body = head.map(cell).to_html(index=False, border=0, classes="rows", escape=True)
    more = f"<p class='more'>… {n(len(frame) - len(head), 'more row')}</p>" if len(frame) > len(head) else ""
    return body + more


def columns(described: dict, present) -> str:
    """What the specification says about the columns that are present."""
    items = [f"<li><code>{html.escape(c)}</code>: {text(described[c])}</li>" for c in present if c in described]
    return f"<ul class='columns'>{''.join(items)}</ul>" if items else ""


def as_json(value) -> str:
    return f"<pre>{html.escape(json.dumps(value, indent=2, ensure_ascii=False))}</pre>"


def row(name, what, size, depth, folder=False) -> str:
    tag = "summary"
    folder = " folder" if folder else ""
    return (f"<{tag} class='row{folder}' style='--d:{depth}'><span class='name'><code>{html.escape(name)}</code></span>"
            f"<span class='what'>{text(what)}</span><span class='size'>{html.escape(size)}</span></{tag}>")


def leaf(name, what, size, body, depth) -> str:
    return f"<li><details>{row(name, what, size, depth)}<div class='body'>{body}</div></details></li>"


def branch(name, what, size, children, depth) -> str:
    """A folder or a group: closed at first, except the top one."""
    return f"<li><details{' open' if depth == 0 else ''}>{row(name, what, size, depth, True)}<ul>{''.join(children)}</ul></details></li>"


def kilobytes(path) -> str:
    files = [path] if path.is_file() else [f for f in path.rglob("*") if f.is_file()]
    return f"{sum(f.stat().st_size for f in files) / 1024:,.0f} KB"


def as_tree(top) -> str:
    return ("<div class='tree'><div class='bar'><button onclick=\"treeAll(this, true)\">Open all folders</button>"
            f"<button onclick=\"treeAll(this, false)\">Close all</button></div><ul>{top}</ul></div>")


def screen_tree() -> str:
    says = screen_spec["descriptions"]
    in_every_modality = {c["name"]: c["description"] for c in screen_spec["cell_tables"]["every_modality"]}
    in_grna_table = {c["name"]: c["description"] for c in screen_spec["grna_table"]}
    children = []
    with h5py.File(SCREEN) as f:
        for name in sorted(f["mod"], key=lambda m: m == screen_spec["modalities"]["grna"]):
            grna = name == screen_spec["modalities"]["grna"]
            counts = ad.io.read_elem(f[f"mod/{name}/X"])
            obs, var = ad.io.read_elem(f[f"mod/{name}/obs"]), ad.io.read_elem(f[f"mod/{name}/var"])
            corner = pd.DataFrame(counts[:ROWS, :6].toarray(), index=obs.index[:ROWS], columns=var.index[:6])
            shape = f"{n(counts.shape[0], 'cell')} × {n(counts.shape[1], 'gRNA' if grna else 'feature')}"
            children.append(branch(f"mod/{name}", says["grna_modality" if grna else "modality"], shape, [
                leaf("X", says["counts"], n(counts.nnz, "nonzero count"),
                     "<p class='more'>The first cells and features:</p>" + rows(corner, "cell") + "<p class='more'>…</p>", 2),
                leaf("obs", says["cell_table"], f"{n(len(obs), 'row')} × {n(obs.shape[1], 'column')}",
                     columns(in_every_modality, obs.columns) + rows(obs, "cell"), 2),
                leaf("var", says["grna_table" if grna else "feature_table"], f"{n(len(var), 'row')} × {n(var.shape[1], 'column')}",
                     (columns(in_grna_table, var.columns) if grna else "") + rows(var, "gRNA" if grna else "feature"), 2),
            ], 1))
        shared = ad.io.read_elem(f["obs"])
        children.append(leaf("obs", says["shared_cell_table"], f"{n(len(shared), 'row')} × {n(shared.shape[1], 'column')}", rows(shared, "cell"), 1))
        entries = screen_spec["metadata"]["entries"]
        stamp = {k: ad.io.read_elem(f[f"uns/{k}"]) for k in entries if k in f["uns"]}
        values = pd.DataFrame([{"entry": k, "value": str(v), "meaning": entries.get(k, {}).get("description", "")} for k, v in stamp.items()])
        values["meaning"] = values["meaning"].str.replace("`", "")
        children.append(leaf(screen_spec["metadata"]["path"], says["metadata"], f"{len(values)} entries", values.to_html(index=False, border=0, classes="rows", escape=True), 1))
    return as_tree(branch(SCREEN.name, says["file"], kilobytes(SCREEN), children, 0))


def status(shows: dict) -> str:
    """What each check of the verification suite found when the page was rendered.
    `shows` says, for each group of checks, what it establishes."""
    found = {}
    for path in sorted(Path("verify/out/results").glob("*.json")):
        name, language = path.stem.rsplit("-", 1)
        found.setdefault(name, {})[language] = [c["result"] for c in json.loads(path.read_text())]
    cases = {c["name"]: c for c in json.loads(Path("conformance/cases.json").read_text())}
    for language in ("r", "python"):
        results = json.loads(Path(f"conformance/results-{language}.json").read_text())
        found.setdefault("conformance", {})[language] = [
            "pass" if sorted(r["rules"]) == cases[name]["expect"] and r["conforms"] == cases[name]["conforms"]
            and False not in [v for k, v in r.items() if k != "conforms"] else "fail"
            for name, r in results.items()]
    if set(found) != set(shows) or any(set(v) != {"r", "python"} for v in found.values()):
        raise SystemExit("Run verify/run.sh before rendering: this page reports its results.")

    def passed(results) -> str:
        fails = results.count("fail")
        return f"{results.count('pass'):,} of {len(results):,}" + ("" if fails == 0 else " — FAILING")
    lines = "".join(f"<tr><td><code>{name}</code></td><td>{text(shows[name])}</td><td>{passed(found[name]['r'])}</td><td>{passed(found[name]['python'])}</td></tr>" for name in shows)
    return f"<table class='status'><thead><tr><th>Check</th><th>What it shows</th><th>Passed in R</th><th>Passed in Python</th></tr></thead><tbody>{lines}</tbody></table>"


def principles(rationale, link, heading) -> str:
    """The principles, read from the table under `heading` in a RATIONALE.md, so
    the page and the document cannot differ."""
    out = []
    for line in Path(rationale).read_text().split(f"\n{heading}\n")[1].split("\n---")[0].splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split(" | ")]
        if len(cells) == 4 and cells[0].isdigit():
            number, principle, why, follows = cells
            follows = follows.replace("](#", f"]({link}#").replace("](../", f"]({link.rsplit('/', 2)[0]}/")
            out.append(f'<details class="fold">\n<summary><strong>{number}. {html.escape(principle.strip("*"))}</strong></summary>\n\n'
                       f"- **Why.** {why}\n- **What follows.** {follows}\n\n</details>\n")
    return "\n".join(out)
