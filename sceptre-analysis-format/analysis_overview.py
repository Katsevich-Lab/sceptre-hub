"""The tree of an analysis folder, for this format's overview page. The rest of
the page (the table of check results, the list of principles, the pieces a tree
is drawn from) comes from the data file format's overview_tools.py. The page sets
ANALYSIS and table_spec before calling analysis_tree()."""
import json
import sys

import pyarrow.parquet as pq

sys.path.insert(0, "../perturb-seq-data-format")
import overview_tools as nt

# Set by the page: the example analysis it shows, and spec/tables.json, from
# which every description on the page is read.
ANALYSIS = table_spec = None


def table_leaf(path, what, depth) -> str:
    file = pq.ParquetFile(path)
    frame = file.read().to_pandas()
    body = nt.rows(frame)
    identity = (file.metadata.metadata or {}).get(table_spec["fits"]["identity_key"].encode())
    if identity is not None:
        body = "<p class='more'>How these fits were made, stored inside the file:</p>" + nt.as_json(json.loads(identity)) + body
    return nt.leaf(path.name, what, f"{nt.n(len(frame), 'row')} × {nt.n(frame.shape[1], 'column')}", body, depth)


def analysis_tree() -> str:
    says = table_spec["descriptions"]
    of_table = {t["name"]: t["description"] for t in table_spec["tables"]}
    children = [nt.leaf("analysis.json", says["analysis.json"], "", nt.as_json(json.loads((ANALYSIS / "analysis.json").read_text())), 1)]
    for stage in table_spec["stages"]:
        folder = ANALYSIS / stage
        if not folder.is_dir():
            continue
        inside = [table_leaf(folder / f"{name.split('/')[1]}.parquet", what, 2)
                  for name, what in of_table.items() if name.startswith(f"{stage}/") and (folder / f"{name.split('/')[1]}.parquet").exists()]
        inside.append(nt.leaf("record.json", says["record.json"], "", nt.as_json(json.loads((folder / "record.json").read_text())), 2))
        children.append(nt.branch(f"{stage}/", says["stages"][stage], nt.n(len(inside), "file"), inside, 1))
    fits = ANALYSIS / table_spec["fits"]["folder"]
    if fits.is_dir():
        parts = [table_leaf(p, says["part"], 2) for p in sorted(fits.glob("*.parquet"))]
        children.append(nt.branch(f"{table_spec['fits']['folder']}/", says["fits"], nt.n(len(parts), "file"), parts, 1))
    return nt.as_tree(nt.branch(f"{ANALYSIS.name}/", says["analysis"], nt.kilobytes(ANALYSIS), children, 0))


