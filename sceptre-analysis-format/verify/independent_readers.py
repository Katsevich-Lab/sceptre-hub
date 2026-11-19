"""Readers we did not write. The tables of the example analyses are read with
polars and duckdb, and what they return is compared with pyarrow.

    python verify/independent_readers.py
"""
import json
import sys
import warnings
from pathlib import Path

import duckdb
import polars as pl
import pyarrow.parquet as pq

sys.path.insert(0, "../perturb-seq-data-format/verify")
from common import Results

warnings.simplefilter("ignore")
results = Results("independent-readers")
for folder in ("example/highmoi/analyses/left_sided", "example/highmoi/analyses/singleton", "example/lowmoi/analyses/default",
               "conformance/valid/analyses/left_sided_written_by_python"):
    for path in sorted(Path(folder).rglob("*.parquet")):
        if "fits" in path.parts:
            continue
        name = f"{path.relative_to(Path(folder).parent)}"
        arrow = pq.read_table(path)
        expected = arrow.to_pylist()
        polars = pl.read_parquet(path)
        results.add(f"polars reads {name} to the same rows as pyarrow", polars.columns == arrow.column_names and polars.to_dicts() == expected)
        duck = duckdb.sql(f"select * from read_parquet('{path}')").to_arrow_table()
        results.add(f"duckdb reads {name} to the same rows as pyarrow", duck.column_names == arrow.column_names and duck.to_pylist() == expected)
    # the fit cache is a folder of part files: both read it as one table
    files = sorted(Path(folder).glob("fits/response/*.parquet"))
    if files:
        parts = f"{folder}/fits/response/*.parquet"
        rows = sum(pq.ParquetFile(p).metadata.num_rows for p in files)
        results.add(f"polars reads the part files of the fit cache of {folder} as one table", pl.read_parquet(parts).height == rows)
        results.add(f"duckdb reads the part files of the fit cache of {folder} as one table",
                    duckdb.sql(f"select count(*) from read_parquet('{parts}')").fetchall()[0][0] == rows)
        # the folder holds part files and nothing else (T6.1), so readers open the folder itself
        whole = f"{folder}/fits/response"
        def opens(read):
            try:
                return read() == rows
            except Exception as error:
                return f"FAILS ({type(error).__name__})"
        import pandas as pd
        import pyarrow.dataset as ds
        for reader, read in {
            "pyarrow.dataset": lambda: ds.dataset(whole, format="parquet").count_rows(),
            "pyarrow.parquet.read_table": lambda: pq.read_table(whole).num_rows,
            "pandas.read_parquet": lambda: len(pd.read_parquet(whole)),
            "polars.scan_parquet": lambda: pl.scan_parquet(whole).collect().height,
            "duckdb read_parquet('folder/*')": lambda: duckdb.sql(f"select count(*) from read_parquet('{whole}/*')").fetchall()[0][0],
        }.items():
            outcome = opens(read)
            results.add(f"{reader} opens the fit cache folder of {folder} itself", outcome is True, "" if outcome is True else outcome)
        # the identity of the fits is in each part's own metadata, which pyarrow and duckdb show
        key = "sceptre_fit_identity"
        seen = pq.ParquetFile(files[0]).metadata.metadata.get(key.encode())
        results.add(f"pyarrow shows the identity inside a part of the fit cache of {folder}", seen is not None and "cells_in_use_sha256" in json.loads(seen))
        shown = duckdb.sql(f"select decode(value) from parquet_kv_metadata('{files[0]}') where decode(key) = '{key}'").fetchall()
        results.add(f"duckdb shows the identity inside a part of the fit cache of {folder}", len(shown) == 1 and "cells_in_use_sha256" in json.loads(shown[0][0]))
results.save()
