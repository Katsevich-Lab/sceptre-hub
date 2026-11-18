"""Encoding matrix, Python side. `write`: each kind of value, written with the
standard Python libraries (not the prototype's writer), one small data file and
one Parquet file per kind. `read`: the files R wrote with its standard
libraries, read here and compared with what was put in.

Kinds the format allows must come back correct. Kinds it forbids are recorded as
observed and compared with verify/expected_encodings.json: they show why the
restriction exists, and a change there means an upstream library has changed.

    python verify/encodings.py write
    python verify/encodings.py read
"""
import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse

sys.path.insert(0, "verify")
from common import OUT, Results

warnings.simplefilter("ignore")
FOLDER = OUT / "encodings"
CELLS = [f"cell_{i}" for i in range(1, 7)]
NAN, INF = float("nan"), float("inf")
# kind -> (allowed by the format, the values put in). Cell-table kinds have six
# values, feature-table kinds two.
KINDS = {
    "double": (True, [0.5, NAN, 1.5, INF, -INF, 3.5]),
    "integer_32": (True, [1, 2, 3, 4, 5, 6]),
    "factor_levels_not_alphabetical": (True, ["b", "a", None, "a", "b", "a"]),
    "text_in_feature_table": (True, ["x", "y"]),
    "boolean": (False, [True, False, True, True, False, False]),
    "integer_64": (False, [1, 2, 3, 4, 5, 2**40]),
    "integer_with_missing": (False, [1, None, 3, 4, None, 6]),
    "boolean_with_missing": (False, [True, None, True, True, None, False]),
    "text_with_missing_in_feature_table": (False, ["x", None]),
    "text_covariate": (False, ["L2", "L1", "L2", "L1", "L2", "L1"]),
}


def column(kind, values):
    if kind == "double":
        return np.array(values, dtype=np.float64)
    if kind == "integer_32":
        return np.array(values, dtype=np.int32)
    if kind == "factor_levels_not_alphabetical":
        return pd.Categorical(values, categories=["b", "a"])
    if kind == "boolean":
        return np.array(values, dtype=bool)
    if kind == "integer_64":
        return np.array(values, dtype=np.int64)
    if kind == "integer_with_missing":
        return pd.array(values, dtype="Int32")
    if kind == "boolean_with_missing":
        return pd.array(values, dtype="boolean")
    return np.array(values, dtype=object)


def data_file(kind, values, path, default_strings=False):
    counts = sparse.csc_matrix(np.array([[3, 0], [1, 0], [0, 2], [5, 1], [0, 0], [2, 4]]))
    obs, var = pd.DataFrame(index=CELLS), pd.DataFrame(index=["gene_1", "gene_2"])
    (var if len(values) == 2 else obs)["value"] = column(kind, values)
    mdata = md.MuData({
        "rna": ad.AnnData(X=counts, obs=obs, var=var),
        "grna": ad.AnnData(X=counts.copy(), obs=pd.DataFrame(index=CELLS), var=pd.DataFrame(index=["grna_1", "grna_2"])),
    })
    with ad.settings.override(allow_write_nullable_strings=default_strings or None if default_strings else False):
        mdata.write(path)


def same(read, put):
    """Whether what was read equals what was put in, element by element."""
    if len(read) != len(put):
        return False
    for r, p in zip(read, put):
        r_missing = r is None or r is pd.NA or (isinstance(r, float) and np.isnan(r))
        p_missing = p is None or (isinstance(p, float) and np.isnan(p))
        if r_missing or p_missing:
            if not (r_missing and p_missing):
                return False
        elif isinstance(p, bool) != isinstance(r, (bool, np.bool_)) or r != p:
            return False
    return True


def describe(series, put):
    """ "correct": the same values in the kind put in. "same values, read as X": the
    values survive and the kind does not."""
    categorical = isinstance(series.dtype, pd.CategoricalDtype)
    values = [None if pd.isna(v) else v for v in series] if categorical else list(series.astype(object).where(series.notna(), None))
    kind = str(series.dtype) + (f", levels {list(series.cat.categories)}" if categorical else "")
    if same(values, put):
        return f"correct (read as {kind})"
    as_text = lambda xs: [None if x is None or (isinstance(x, float) and np.isnan(x)) else str(x) for x in xs]
    loose = [None if v is None else (v.item() if isinstance(v, np.generic) else v) for v in values]
    if as_text(loose) == as_text(put) or [None if v is None else float(v) for v in loose if not isinstance(v, str)] == [None if v is None else float(v) for v in put if not isinstance(v, str)] and all(not isinstance(v, str) for v in loose):
        return f"same values, read as {kind}"
    return f"values differ (read as {kind})"


results = Results("encodings")
expected = json.loads(Path("verify/expected_encodings.json").read_text()) if Path("verify/expected_encodings.json").exists() else {}


def record(name, allowed, outcome):
    if allowed:
        results.add(name, outcome.startswith("correct"), outcome)
    elif name in expected:
        results.add(name, True if expected[name] == outcome else f"CHANGED: {outcome}", f"forbidden kind; as expected: {expected[name]}")
    else:
        results.add(name, outcome, "forbidden kind; no expectation recorded")
    observed[name] = outcome


observed = {}
if sys.argv[1] == "write":
    (FOLDER / "python").mkdir(parents=True, exist_ok=True)
    for kind, (allowed, values) in KINDS.items():
        data_file(kind, values, FOLDER / "python" / f"{kind}.h5mu")
    data_file("double", KINDS["double"][1], FOLDER / "python" / "default_string_settings.h5mu", default_strings=True)
    table = pa.table({
        "text": pa.array(["a", "b", "c"]), "text_dictionary": pa.array(["a", "b", "a"]).dictionary_encode(),
        "boolean_with_missing": pa.array([True, None, False]), "integer_32": pa.array([1, 2, 3], type=pa.int32()),
        "integer_64": pa.array([1, 2, 2**40], type=pa.int64()), "double": pa.array([0.5, None, INF]),
        "list_of_doubles": pa.array([[1.5, 2.5], [], [3.5]], type=pa.list_(pa.float64())),
    })
    pq.write_table(table, FOLDER / "python" / "table.parquet")
    print("Python wrote", len(list((FOLDER / "python").iterdir())), "files")
    sys.exit(0)

# read what R wrote
for kind, (allowed, values) in KINDS.items():
    if kind == "integer_64":
        continue   # R has no 64-bit integer to write
    path = FOLDER / "r" / f"{kind}.h5mu"
    try:
        mdata = md.read_h5mu(path)
        frame = mdata["rna"].var if len(values) == 2 else mdata["rna"].obs
        outcome = describe(frame["value"], values)
    except Exception as error:
        outcome = f"cannot open ({type(error).__name__})"
    record(f"data file, {kind}: R writes, Python reads", allowed, outcome)

put = {"text": ["a", "b", "c"], "text_dictionary": ["a", "b", "a"], "boolean_with_missing": [True, None, False], "integer_32": [1, 2, 3],
       "double": [0.5, None, INF], "list_of_doubles": [[1.5, 2.5], [], [3.5]]}
readers = {
    "pyarrow": lambda p: pq.read_table(p).to_pandas(),
    "polars": lambda p: __import__("polars").read_parquet(p).to_pandas(),
    "duckdb": lambda p: __import__("duckdb").sql(f"select * from read_parquet('{p}')").df(),
}
for reader, read in readers.items():
    try:
        frame = read(FOLDER / "r" / "table.parquet")
    except Exception as error:
        results.add(f"table: R (nanoparquet) writes, {reader} reads", False, f"{type(error).__name__}: {error}")
        continue
    for name, values in put.items():
        got = frame[name]
        if name == "list_of_doubles":
            ok = [list(map(float, v)) for v in got] == values
        else:
            ok = same([None if pd.isna(v) else (v.item() if isinstance(v, np.generic) else v) for v in got.astype(object)], values)
        results.add(f"table, {name}: R (nanoparquet) writes, {reader} reads", bool(ok), f"read as {got.dtype}")
if not expected:
    Path(OUT / "encodings" / "observed-python.json").write_text(json.dumps(observed, indent=1) + "\n")
results.save()
