"""Encoding matrix for the tables, Python side. `write`: one Parquet file holding
each kind of column, written with pyarrow (not the prototype's writer). `read`:
the file R wrote with nanoparquet, read here with pyarrow, polars and duckdb and
compared with what was put in. verify/encodings.R is the R side.

Kinds the format allows must come back correct. A kind it forbids is recorded as
observed and compared with verify/expected_encodings.json: it shows why the
restriction exists, and a change there means an upstream library has changed.

    python verify/encodings.py write
    python verify/encodings.py read
"""
import sys
import warnings

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, "../perturb-seq-data-format/verify")
from common import OUT, Results

warnings.simplefilter("ignore")
FOLDER = OUT / "encodings"
INF = float("inf")

if sys.argv[1] == "write":
    (FOLDER / "python").mkdir(parents=True, exist_ok=True)
    table = pa.table({
        "text": pa.array(["a", "b", "c"]), "text_dictionary": pa.array(["a", "b", "a"]).dictionary_encode(),
        "boolean_with_missing": pa.array([True, None, False]), "integer_32": pa.array([1, 2, 3], type=pa.int32()),
        "integer_64": pa.array([1, 2, 2**40], type=pa.int64()), "double": pa.array([0.5, None, INF]),
        "list_of_doubles": pa.array([[1.5, 2.5], [], [3.5]], type=pa.list_(pa.float64())),
    })
    pq.write_table(table, FOLDER / "python" / "table.parquet")
    print("Python wrote", FOLDER / "python" / "table.parquet")
    sys.exit(0)


def same(read, put):
    """Whether what was read equals what was put in, element by element."""
    if len(read) != len(put):
        return False
    for r, p in zip(read, put):
        r_missing = r is None or r is pd.NA or (isinstance(r, float) and np.isnan(r))
        if r_missing or p is None:
            if not (r_missing and p is None):
                return False
        elif isinstance(p, bool) != isinstance(r, (bool, np.bool_)) or r != p:
            return False
    return True


# read what R wrote
results = Results("encodings")
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
results.save()
