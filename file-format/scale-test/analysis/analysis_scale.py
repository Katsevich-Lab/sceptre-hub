"""Made-up tables of an analysis folder at real size: results, resampled
statistics and the fit cache. No statistics are computed; only the shapes,
types and sizes are those of a real analysis.

    python analysis_scale.py run <work directory> <results.json>

Each measurement runs in a process of its own, so that its peak memory is its
own. Before a read, the operating system is asked to drop the file from its
cache; on a network file system that request may not be honoured."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

HERE = Path(__file__).resolve()
N_RESAMPLES = 499
N_COEFFICIENTS = 20
N_FIT_RESPONSES = 20_000
IDENTITY = {
    "screen_id": "made-up", "formula": "~ log(response_n_umis) + batch",
    "cells_in_use_sha256": "0" * 64, "response_fit_method": "glm_nb", "control_group": "complement",
}


def names(prefix, width, index, n_levels):
    """Identifiers as text, in sorted order of `index`."""
    levels = pa.array([f"{prefix}{i:0{width}d}" for i in range(n_levels)])
    return pa.DictionaryArray.from_arrays(pa.array(index.astype(np.int32)), levels).cast(pa.string())


def keys(start, n, n_responses, n_targets):
    index = np.arange(start, start + n)
    return [names("target_", 6, index // n_responses, n_targets), names("ENSG", 11, index % n_responses, n_responses)]


# ------------------------------------------------------------------ writers


def write_results(path, n, n_responses):
    rng = np.random.default_rng(1)
    schema = pa.schema([("grna_target", pa.string()), ("response_id", pa.string()), ("pass_qc", pa.bool_()),
                        ("p_value", pa.float64()), ("fold_change", pa.float64()), ("se_fold_change", pa.float64())])
    with pq.ParquetWriter(path, schema, compression="zstd") as writer:
        for start in range(0, n, 1_000_000):
            m = min(1_000_000, n - start)
            columns = keys(start, m, n_responses, n // n_responses) + [
                pa.array(rng.random(m) < 0.9), pa.array(rng.random(m)),
                pa.array(np.exp(rng.normal(0, 0.2, m))), pa.array(rng.random(m) * 0.1)]
            writer.write_table(pa.Table.from_arrays(columns, schema=schema), row_group_size=m)
    return {"rows": n}


def write_resamples(path, n, n_responses, layout, row_group, compression):
    """layout: `list` (the format), `list32` (32-bit values), `wide` (one column
    per resample, as packages hold them in memory)."""
    rng = np.random.default_rng(2)
    value_type = pa.float32() if layout == "list32" else pa.float64()
    key_fields = [("grna_target", pa.string()), ("response_id", pa.string())]
    if layout == "wide":
        schema = pa.schema(key_fields + [(f"z_null_{j + 1}", pa.float64()) for j in range(N_RESAMPLES)])
    else:
        schema = pa.schema(key_fields + [("z_null", pa.list_(value_type))])
    with pq.ParquetWriter(path, schema, compression=compression) as writer:
        for start in range(0, n, row_group):
            m = min(row_group, n - start)
            z = rng.standard_normal((m, N_RESAMPLES))
            if layout == "wide":
                columns = [pa.array(z[:, j]) for j in range(N_RESAMPLES)]
            else:
                offsets = pa.array(np.arange(0, (m + 1) * N_RESAMPLES, N_RESAMPLES, dtype=np.int64))
                values = pa.array(z.ravel().astype(value_type.to_pandas_dtype()))
                columns = [pa.LargeListArray.from_arrays(offsets, values).cast(pa.list_(value_type))
                           if m * N_RESAMPLES > 2**31 - 1 else
                           pa.ListArray.from_arrays(offsets.cast(pa.int32()), values)]
            writer.write_table(pa.Table.from_arrays(keys(start, m, n_responses, n // n_responses) + columns,
                                                    schema=schema), row_group_size=m)
    return {"rows": n}


def write_fits(folder, n_parts):
    rng = np.random.default_rng(3)
    folder = Path(folder)
    folder.mkdir(parents=True)
    schema = pa.schema([("response_id", pa.string()), ("coefficients", pa.list_(pa.float64())), ("theta", pa.float64())],
                       metadata={"sceptre_fit_identity": json.dumps(IDENTITY, sort_keys=True)})
    per_part = N_FIT_RESPONSES // n_parts
    for p in range(n_parts):
        index = np.arange(p * per_part, (p + 1) * per_part)
        offsets = pa.array(np.arange(0, (per_part + 1) * N_COEFFICIENTS, N_COEFFICIENTS, dtype=np.int32))
        table = pa.Table.from_arrays([
            names("ENSG", 11, index, N_FIT_RESPONSES),
            pa.ListArray.from_arrays(offsets, pa.array(rng.standard_normal(per_part * N_COEFFICIENTS))),
            pa.array(rng.random(per_part) * 10)], schema=schema)
        pq.write_table(table, folder / f"part-{p:012d}.parquet")
    return {"parts": n_parts, "rows": per_part * n_parts}


# ------------------------------------------------------------------ readers


def read_all(path):
    table = pq.read_table(path)
    return {"rows": table.num_rows, "gb_in_memory": round(table.nbytes / 1e9, 2)}


def read_keys(path):
    return {"rows": pq.read_table(path, columns=["grna_target", "response_id"]).num_rows}


def read_where(path, column, value, value_2=None):
    filters = [(column, "=", value)] + ([("response_id", "=", value_2)] if value_2 else [])
    return {"rows": pq.read_table(path, filters=filters).num_rows}


def fits_identities(folder):
    files = sorted(Path(folder).glob("*.parquet"))
    found = {pq.ParquetFile(f).metadata.metadata[b"sceptre_fit_identity"] for f in files}
    return {"parts": len(files), "identities": len(found)}


def fits_read_each(folder):
    """As the prototype reads: the identity of each part, then each part."""
    files = sorted(Path(folder).glob("*.parquet"))
    for f in files:
        pq.ParquetFile(f).metadata.metadata[b"sceptre_fit_identity"]
    return {"rows": pa.concat_tables([pq.read_table(f) for f in files]).num_rows}


def fits_read_dataset(folder):
    return {"rows": ds.dataset(folder, format="parquet").to_table().num_rows}


OPERATIONS = {f.__name__: f for f in (write_results, write_resamples, write_fits, read_all, read_keys, read_where,
                                      fits_identities, fits_read_each, fits_read_dataset)}


# ------------------------------------------------------------------ driver


def size_of(path):
    path = Path(path)
    files = list(path.glob("*")) if path.is_dir() else [path]
    return sum(f.stat().st_size for f in files)


def drop_from_cache(path):
    path = Path(path)
    for f in (path.glob("*") if path.is_dir() else [path]):
        fd = os.open(f, os.O_RDONLY)
        os.fsync(fd)
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        os.close(fd)


def measure(command):
    """Run one measurement in a process of its own: its answer, seconds and peak memory."""
    done = subprocess.run(["/usr/bin/time", "-f", "%M"] + command, capture_output=True, text=True)
    if done.returncode != 0:
        return {"failed": done.stderr.strip().splitlines()[-3:]}
    out = json.loads(done.stdout.strip().splitlines()[-1])
    out["peak_gb"] = round(int(done.stderr.strip().splitlines()[-1]) / 1e6, 2)
    return out


def python(operation, *arguments):
    return measure([sys.executable, str(HERE), "one", operation] + [str(a) for a in arguments])


def r(operation, *arguments):
    return measure(["Rscript", str(HERE.with_suffix(".R")), operation] + [str(a) for a in arguments])


def run(work, out_file):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    results = {"libraries": {"python": sys.version.split()[0], "pyarrow": pa.__version__, "numpy": np.__version__},
               "resamples_per_pair": N_RESAMPLES, "results": [], "resamples": [], "fits": []}

    def save():
        Path(out_file).write_text(json.dumps(results, indent=1) + "\n")

    def reads(path, entry, whole, in_r):
        for name, arguments in (("read_keys", ()), ("read_where", ("grna_target", "target_000007")),
                                ("read_where", ("grna_target", "target_000007", "ENSG00000000042")),
                                ("read_where", ("response_id", "ENSG00000000042"))) + ((("read_all", ()),) if whole else ()):
            drop_from_cache(path)
            label = name if name != "read_where" else (
                "one_response" if arguments[0] == "response_id" else "one_pair" if len(arguments) == 3 else "one_target")
            entry[label] = python(name, path, *arguments)
            print(path.name, label, entry[label], flush=True)
        if in_r:
            for reader in ("nanoparquet", "arrow"):
                drop_from_cache(path)
                entry[f"r_{reader}_read_all"] = r("read_all", reader, path)
                print(path.name, "R", reader, entry[f"r_{reader}_read_all"], flush=True)

    for n in (20_000_000, 100_000_000):
        path = work / f"results-{n}.parquet"
        entry = {"rows": n, "write": python("write_results", path, n, 20_000)}
        entry["gb_on_disk"] = round(size_of(path) / 1e9, 2)
        reads(path, entry, whole=True, in_r=True)
        results["results"].append(entry)
        path.unlink()
        save()

    cases = [(1_000_000, "list", 1_000_000, "snappy"), (1_000_000, "list", 1_000_000, "zstd"),
             (1_000_000, "list", 10_000, "zstd"), (1_000_000, "list32", 10_000, "zstd"),
             (1_000_000, "wide", 10_000, "zstd"),
             (10_000_000, "list", 1_000_000, "zstd"), (10_000_000, "list", 10_000, "zstd")]
    for n, layout, row_group, compression in cases:
        path = work / f"resamples-{n}-{layout}-{row_group}-{compression}.parquet"
        entry = {"rows": n, "layout": layout, "rows_per_row_group": row_group, "compression": compression,
                 "write": python("write_resamples", path, n, 10_000, layout, row_group, compression)}
        entry["gb_on_disk"] = round(size_of(path) / 1e9, 2)
        reads(path, entry, whole=True, in_r=(n == 1_000_000 and layout == "list"))
        results["resamples"].append(entry)
        path.unlink()
        save()

    for n_parts in (1, 100, 1_000, 10_000, 20_000):
        folder = work / f"fits-{n_parts}"
        entry = {"parts": n_parts, "rows_per_part": N_FIT_RESPONSES // n_parts, "write": python("write_fits", folder, n_parts)}
        entry["mb_on_disk"] = round(size_of(folder) / 1e6, 2)
        for name in ("fits_identities", "fits_read_each", "fits_read_dataset"):
            drop_from_cache(folder)
            entry[name] = python(name, folder)
        drop_from_cache(folder)
        entry["r_nanoparquet_read_each"] = r("fits_read_each", "nanoparquet", folder)
        entry["warm_fits_read_each"] = python("fits_read_each", folder)
        print("fits", entry, flush=True)
        results["fits"].append(entry)
        subprocess.run(["rm", "-rf", str(folder)])
        save()


def run_r(work, out_file):
    """R's arrow on the same tables: whole reads and reads of chosen rows."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    results = []
    cases = [("results", 100_000_000, None, None), ("resamples", 1_000_000, 1_000_000, True),
             ("resamples", 1_000_000, 10_000, True), ("resamples", 10_000_000, 10_000, False)]
    for table, n, row_group, whole in cases:
        path = work / f"{table}-{n}-{row_group}.parquet"
        if table == "results":
            python("write_results", path, n, 20_000)
        else:
            python("write_resamples", path, n, 10_000, "list", row_group, "zstd")
        entry = {"table": table, "rows": n, "rows_per_row_group": row_group or 1_000_000}
        for label, arguments in (("one_target", ("grna_target", "target_000007")),
                                 ("one_pair", ("grna_target", "target_000007", "ENSG00000000042")),
                                 ("one_response", ("response_id", "ENSG00000000042"))):
            drop_from_cache(path)
            entry[label] = r("read_where", "arrow", path, *arguments)
        for reader in ("arrow_table", "arrow", "nanoparquet") if whole is not False else ():
            drop_from_cache(path)
            entry[f"read_all_{reader}"] = r("read_all", reader, path)
        print(entry, flush=True)
        results.append(entry)
        path.unlink()
        Path(out_file).write_text(json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "run_r":
        run_r(sys.argv[2], sys.argv[3])
    else:
        arguments = [int(a) if a.isdigit() else a for a in sys.argv[3:]]
        start = time.perf_counter()
        answer = OPERATIONS[sys.argv[2]](*arguments)
        answer["seconds"] = round(time.perf_counter() - start, 2)
        print(json.dumps(answer))
