"""Files of the sceptre format <-> pysceptre's Analysis.

This is the package-dependent half of the prototype. It needs the sceptre3 build
of pysceptre, an unreviewed branch build, and reaches into its internals. It is
here to test that the files hold everything an analysis needs, and is meant to
move into the package when the package adopts the format. Nothing in format/ or
verify/ uses it.

    write_sceptre_screen(), write_sceptre_analysis(), read_sceptre_object()

pysceptre calls the response modality ``response`` whatever was measured, and the
gRNA table's target coordinates ``chr``, ``start``, ``end``. The functions here
translate between those names and the file's.

Everything of format/sceptre_format.py is available from this module too, so a
notebook can import this one module.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "format"))
import sceptre_format as _format

globals().update({name: value for name, value in vars(_format).items() if not name.startswith("__")})

def design_matrix(formula: str, covariates: pd.DataFrame):
    """pysceptre's `model_matrix`, for covariate names with dots as well.

    pysceptre's formula parser takes names of letters, digits and underscores; the
    format also allows dots (S5.4). A dotted name is given a dot-free alias for the
    call and put back in the names of the terms."""
    from pysceptre.pipeline.formula import model_matrix

    dotted = [c for c in covariates.columns if "." in str(c)]
    alias = {}
    for name in dotted:
        candidate = str(name).replace(".", "_DOT_")
        while candidate in covariates.columns or candidate in alias.values():
            candidate += "_"
        alias[name] = candidate

    def swap(text, mapping):
        for old, new in mapping.items():
            text = re.sub(rf"(?<![A-Za-z0-9_.]){re.escape(old)}(?![A-Za-z0-9_.])", new, text)
        return text

    from types import SimpleNamespace

    built = model_matrix(swap(formula, alias), covariates.rename(columns=alias))
    back = {v: k for k, v in alias.items()}
    return SimpleNamespace(
        matrix=built.matrix,
        terms=[swap(t, back) for t in built.terms],
        assign=built.assign,
        contrasts={back.get(k, k): v for k, v in built.contrasts.items()},
    )




# The user-facing name "spa" is the default saddlepoint route of each mechanism.
SPA_DEFAULT = {"crt": "crt_spa_fast", "permutations": "rpt_spa_fast"}
RESULT_STAGES = ("calibration", "power", "discovery")


def _r_value(value):
    """A plain JSON value as the R value pysceptre's state holds."""
    from pysceptre.analysis.state import r_chr, r_dbl, r_int, r_lgl, r_list

    if isinstance(value, dict):
        return r_list({k: _r_value(v) for k, v in value.items()})
    values = value if isinstance(value, list) else [value]
    if all(isinstance(v, bool) for v in values):
        return r_lgl(values)
    if all(isinstance(v, int) for v in values):
        return r_int(values)
    if all(isinstance(v, (int, float)) for v in values):
        return r_dbl(values)
    return r_chr(values)


def _plain(value):
    """An R value of pysceptre's state as plain JSON values."""
    from pysceptre.io.r_json import to_python

    def convert(v):
        if isinstance(v, dict):
            return {k: convert(x) for k, x in v.items()}
        if isinstance(v, np.ndarray):
            v = v.tolist()
        if isinstance(v, (list, tuple)):
            v = [convert(x) for x in v]
            return v[0] if len(v) == 1 else v
        return v.item() if isinstance(v, np.generic) else v

    return convert(to_python(value))


def _swap_prefix(text, old: str, new: str):
    """`old_` -> `new_` at the start of each name in a formula, term or column name."""
    if isinstance(text, (list, tuple, pd.Index)):
        return [_swap_prefix(t, old, new) for t in text]
    return re.sub(rf"(?<![A-Za-z0-9_.]){re.escape(old)}_", f"{new}_", str(text))


def to_file_names(text, modality: str):
    """pysceptre prefixes the response modality's covariates with `response_`; the
    file uses the modality's name."""
    return _swap_prefix(text, "response", modality)


def to_package_names(text, modality: str):
    return _swap_prefix(text, modality, "response")


def _unit_group(table: pd.DataFrame, strategy: str) -> np.ndarray:
    return table["grna_target" if strategy == "union" else "grna_id"].to_numpy(dtype=object)


def _unit_columns(grna_group, grna_table: pd.DataFrame, strategy: str) -> dict:
    grna_group = np.asarray(grna_group, dtype=object)
    if strategy == "union":
        return {"grna_target": grna_group}
    target_of = dict(zip(grna_table["grna_id"], grna_table["grna_target"], strict=True))
    return {"grna_target": np.array([target_of[g] for g in grna_group], dtype=object), "grna_id": grna_group}


def _pysceptre_screen(mdata):
    """The MuData in the layout pysceptre's Analysis takes: every covariate in the
    response modality's obs, and the top gRNA of each cell in the gRNA modality's."""
    import anndata as ad
    import mudata as md

    main = response_of(mdata)
    summary = top_grna(mdata[GRNA].X, mdata[GRNA].var_names)
    summary.index = mdata[GRNA].obs_names
    covariates = covariate_table(mdata)
    covariates.columns = to_package_names(covariates.columns, main)
    grna_table = mdata[GRNA].var.rename(columns={v: k for k, v in GRNA_COORDINATES.items()})
    out = md.MuData(
        {
            "response": ad.AnnData(X=mdata[main].X, obs=covariates, var=mdata[main].var),
            "grna": ad.AnnData(X=mdata[GRNA].X, obs=summary, var=grna_table),
        }
    )
    out.uns = {"low_moi": np.array([mdata.uns["sceptre"]["moi"] == "low"])}
    out.uns["grna_target_data_frame"] = grna_table.reset_index(names="grna_id")
    return out


def _fit_identity(analysis, mdata) -> dict:
    """What makes a saved fit valid (T6)."""
    parameters = analysis.parameters
    cells = np.asarray(mdata[GRNA].obs_names, dtype=object)
    return {
        "screen_id": mdata.uns["sceptre"]["id"],
        "formula": to_file_names(parameters.formula, response_of(mdata)),
        "cells_in_use_sha256": sha256_ids(cells[analysis.qc.cells_in_use]),
        "response_fit_method": parameters.response_fit_method,
        "control_group": parameters.control_group,
    }


def _fit_mismatch(saved: dict, current: dict):
    """Why cached fits cannot be reused, or None. A method name means the same
    estimator in every language that lists it; a method this language does not
    share is never reused."""
    for name in current:
        if _identity_value(saved, name) != _identity_value(current, name):
            return f"their `{name}` differs"
    method = current["response_fit_method"]
    if "Python" not in FIT_METHODS.get(method, ()):
        return f"the method '{method}' is not one Python shares"
    return None


def read_sceptre_object(screen_file, directory=None, *, response_modality=None, verbose=True):
    """A data file, and optionally one analysis of it, as a pysceptre Analysis.
    Without `directory` the analysis is as after ``import_data()``. The response
    modality is the one the analysis names; without an analysis, `response_modality`
    or the file's only modality besides ``grna``."""
    from pysceptre import Analysis
    from pysceptre.analysis.records import (
        AnalysisParameters,
        CalibrationSettings,
        CellQC,
        DesignMatrix,
        GuideAssignment,
    )
    if directory is not None:
        response_modality = analysis_response_modality(directory)
    mdata = read_screen(screen_file, response_modality)
    main = response_of(mdata)
    cells = pd.Index(mdata[main].obs_names)
    response_ids = pd.Index(mdata[main].var_names)
    covariates = covariate_table(mdata)
    covariates.columns = to_package_names(covariates.columns, main)
    low_moi = mdata.uns["sceptre"]["moi"] == "low"
    analysis = Analysis.from_screen(_pysceptre_screen(mdata))
    analysis.screen_id = mdata.uns["sceptre"]["id"]
    analysis.screen_metadata = {k: mdata.uns["sceptre"][k] for k in OPTIONAL_METADATA if k in mdata.uns["sceptre"]}
    analysis.format_provenance = {}
    if directory is None:
        return analysis
    saved = read_analysis(directory)
    check_link(saved["screen"], mdata, directory)
    done, records, tables = list(saved["stages"]), saved["stages"], saved["tables"]
    analysis.format_provenance = {s: r["provenance"] for s, r in records.items()}

    if "parameters" in done:
        settings = records["parameters"]["settings"]
        strategy = settings["grna_integration_strategy"]
        grna_table = analysis.grna_target_data_frame.copy()
        grna_table["grna_group"] = grna_groups(grna_table, strategy)
        analysis.set_grna_target_data_frame(grna_table)
        approximation = settings["resampling_approximation"]
        if approximation == "spa":
            approximation = SPA_DEFAULT[settings["resampling_mechanism"]]
        # in the file, a covariate carries its modality's name
        formula = to_package_names(settings["formula"], main)
        built = design_matrix(formula, covariates)
        recorded = records["parameters"]["design"]
        if to_file_names(list(built.terms), main) != list(recorded["terms"]) or design_sha256(built.matrix) != recorded["sha256"]:
            raise ValueError("The design matrix rebuilt from the formula is not the one the analysis was run with.")
        analysis.set_parameters(
            AnalysisParameters(
                formula=formula,
                side=settings["side"],
                resampling_approximation=approximation,
                control_group=settings["control_group"],
                resampling_mechanism=settings["resampling_mechanism"],
                B1=int(settings["B1"]),
                B2=int(settings["B2"]),
                B3=int(settings["B3"]),
                multiple_testing_method=settings["multiple_testing_method"],
                multiple_testing_alpha=float(settings["multiple_testing_alpha"]),
                grna_integration_strategy=strategy,
                # counted as the user gave them, before any expansion by gRNA
                n_discovery_pairs=len(tables["parameters/discovery_pairs"]),
                n_positive_control_pairs=len(tables["parameters/positive_control_pairs"]),
                response_fit_method=settings["response_fit_method"],
                grna_fit_method=settings["grna_fit_method"],
                seed=settings["seed"],
            ),
            DesignMatrix(
                matrix=built.matrix,
                terms=list(built.terms),
                assign=list(built.assign),
                contrasts={k: [v] for k, v in built.contrasts.items()} or None,
            ),
            expand_pairs(tables["parameters/discovery_pairs"], grna_table),
            expand_pairs(tables["parameters/positive_control_pairs"], grna_table),
        )

    if "assignment" in done:
        settings = records["assignment"]["settings"]
        initial = cells_by_grna(tables["assignment/assignments"], analysis.grna_ids, cells)
        raw_groups, raw_nt = group_cells(initial, grna_table)
        flagged = None
        if low_moi:
            listed = tables.get("assignment/flagged_cells")
            listed = cells.get_indexer(listed["cell_id"]) if listed is not None else np.empty(0, dtype=np.int64)
            flagged = np.union1d(not_one_grna(initial, len(cells)), listed).astype(np.int64)
        analysis.set_assignment(
            GuideAssignment(
                method=settings["method"],
                hyperparameters=_r_value(settings["hyperparameters"]),
                initial=initial,
                group_cells=raw_groups,
                nt_cells=raw_nt,
                mean_cells_per_grna=float(np.mean([len(v) for v in initial.values()])),
                # the "maximum" method does not count gRNAs per cell
                grnas_per_cell=None if settings["method"] == "maximum" else grnas_per_cell(initial, len(cells)),
                cells_w_zero_or_twoplus_grnas=flagged if flagged is not None and flagged.size else None,
            )
        )

    if "qc" in done:
        settings = records["qc"]["settings"]
        in_use = np.flatnonzero(tables["qc/cells"].set_index("cell_id").loc[cells, "in_use"].to_numpy(dtype=bool)).astype(np.int64)
        nt_control = records["parameters"]["settings"]["control_group"] == "nt_cells"
        groups, nt, all_nt = qc_assignments(raw_groups, raw_nt, in_use, nt_control)
        info = {}
        for pair_type in ("discovery", "positive_control"):
            table = tables[f"qc/{pair_type}_pair_qc"]
            frame = pd.DataFrame(
                {
                    "response_id": table["response_id"].to_numpy(dtype=object),
                    "grna_group": _unit_group(table, strategy),
                    "pass_qc": table["pass_qc"].to_numpy(dtype=bool),
                }
            )
            # the order a package tests in: by response, in the data file's order
            order = np.argsort(response_ids.get_indexer(frame["response_id"]), kind="stable")
            info[pair_type] = frame.iloc[order].reset_index(drop=True)
        # per response, the cells of the control pool that express it
        pool = in_use[all_nt] if nt_control else in_use
        n_nonzero = np.asarray((mdata[main].X.tocsr()[pool] > 0).sum(axis=0)).reshape(-1)
        analysis.set_qc(
            CellQC(
                response_n_umis_range=tuple(settings["response_n_umis_range"]),
                response_n_nonzero_range=tuple(settings["response_n_nonzero_range"]),
                p_mito_threshold=float(settings["p_mito_threshold"]),
                cell_removal_metrics={f"n_cells_rm_{k}": int(v) for k, v in records["qc"]["summaries"]["cells_removed"].items()},
                cells_in_use=in_use,
                group_cells=groups,
                nt_cells=nt,
                n_nonzero_tot=n_nonzero.astype(np.int64),
                n_ok_discovery_pairs=n_pairs_passing(tables["qc/discovery_pair_qc"], strategy) if len(info["discovery"]) else None,
                n_ok_positive_control_pairs=n_pairs_passing(tables["qc/positive_control_pair_qc"], strategy) if len(info["positive_control"]) else None,
                all_nt=all_nt,
                pairwise_qc=bool(settings["pairwise_qc"]),
            ),
            info["discovery"],
            info["positive_control"],
            B3=int(settings["B3"]),
        )

    def result(stage):
        table = add_log_2_fold_change(tables[f"{stage}/results"])
        table = add_resamples(table, tables.get(f"{stage}/resamples"))
        if stage != "power":
            table = add_significant(table, records[stage].get("significance_threshold"))
            table["significant"] = table["significant"].astype(object).where(table["significant"].notna(), None)
        # as sceptre returns results: by p-value, then response
        table = table.sort_values(["p_value", "response_id", *unit_key(table)], kind="stable", na_position="last")
        return table.reset_index(drop=True)

    if "calibration" in done:
        settings = records["calibration"]["settings"]
        calibration = tables["calibration/results"]
        negative = pd.DataFrame(
            {
                "response_id": calibration["response_id"].to_numpy(dtype=object),
                "grna_group": _unit_group(calibration, strategy),
                "pass_qc": calibration["pass_qc"].to_numpy(dtype=bool),
            }
        )
        analysis.set_calibration(
            CalibrationSettings(int(settings["calibration_group_size"]), int(settings["n_calibration_pairs"])),
            result("calibration"),
            negative,
            B3=int(settings["B3"]),
        )
    for stage in ("power", "discovery"):
        if stage in done:
            analysis.set_result(STAGE_FUNCTIONS[stage], result(stage))

    # cached fits are reused only if they were made as this analysis would make them (T6)
    identity = saved.get("fits_identity")
    if identity is not None and "qc" in done and FITS in tables:
        reason = _fit_mismatch(identity, _fit_identity(analysis, mdata))
        if reason:
            if verbose:
                print(f"The fits saved in {FITS} are not reused: {reason}.")
        else:
            table = tables[FITS]
            terms = list(analysis.design.terms)
            long = pd.DataFrame(
                {
                    "response_id": np.repeat(table["response_id"].to_numpy(dtype=object), len(terms) + 1),
                    "term": np.tile(np.array(terms + [None], dtype=object), len(table)),
                    "value": np.concatenate([np.append(np.asarray(c, dtype=float), t) for c, t in zip(table["coefficients"], table["theta"], strict=True)]),
                }
            )
            analysis.add_fits(long)
    return analysis


def analysis_screen(analysis, id=None, response_modality=DEFAULT_RESPONSE_MODALITY):
    """The data of a pysceptre Analysis as the MuData `write_screen` takes."""
    import anndata as ad
    import mudata as md

    screen = analysis.screen
    grna_table = analysis.grna_target_data_frame.drop(columns="grna_group", errors="ignore")
    if grna_table["grna_id"].duplicated().any():
        raise ValueError("A gRNA listed under two targets cannot be written in this format.")
    grna_table = grna_table.set_index("grna_id").loc[list(screen["grna"].var_names)]
    grna_table.index.name = None
    grna_table = grna_table.rename(columns=GRNA_COORDINATES)
    parts = split_covariates(screen["response"].obs, ("response", "grna"))
    out = md.MuData(
        {
            response_modality: ad.AnnData(X=screen["response"].X, obs=parts["response"], var=screen["response"].var),
            GRNA: ad.AnnData(X=screen["grna"].X, obs=parts["grna"], var=grna_table),
        }
    )
    for name in parts["shared"].columns:
        out.obs[name] = parts["shared"][name].array
    out.uns["sceptre"] = {
        "moi": analysis.moi,
        "id": id or getattr(analysis, "screen_id", None),
        "response_modality": response_modality,
        **getattr(analysis, "screen_metadata", {}),
    }
    return out


def write_sceptre_screen(analysis, file, id=None, response_modality=DEFAULT_RESPONSE_MODALITY) -> None:
    write_screen(analysis_screen(analysis, id, response_modality), file)


def stages_called(analysis) -> list[str]:
    called = analysis.functions_called
    return [s for s in STAGES if called.get(STAGE_FUNCTIONS[s])]


def write_sceptre_analysis(analysis, screen_file, directory, provenance=None, response_modality=None) -> Path:
    """Write the analysis of a pysceptre Analysis into a new folder, as the stages
    would have written it one after another."""
    from importlib.metadata import version

    mdata = read_screen(screen_file, response_modality)
    main = response_of(mdata)
    cells = pd.Index(mdata[main].obs_names)
    if list(analysis.screen["response"].obs_names) != list(cells):
        raise ValueError("The analysis and the data file have different cells.")
    done = stages_called(analysis)
    provenance = provenance or getattr(analysis, "format_provenance", {}) or {}
    package = {"package": "pysceptre", "version": version("pysceptre"), "language": "Python"}
    init_analysis(directory, screen_link(mdata, screen_file, directory), main)

    def write(stage, record, tables=None):
        tables = {k: v for k, v in (tables or {}).items() if v is not None}
        write_stage(directory, stage, record, tables, provenance.get(stage, package))

    if "parameters" in done:
        p = analysis.parameters
        strategy = p.grna_integration_strategy
        grna_table = analysis.grna_target_data_frame
        approximation = p.resampling_approximation
        if approximation == SPA_DEFAULT[p.resampling_mechanism]:
            approximation = "spa"
        if not formula_in_subset(p.formula):
            raise ValueError(
                f"The formula {p.formula} cannot be written: the format takes additive terms only, "
                "each a covariate name, log(x) or log(x + 1)."
            )
        design = analysis.design
        if not np.array_equal(design_matrix(p.formula, analysis.cell_covariates).matrix, design.matrix):
            raise ValueError("The analysis's design matrix is not the one its formula and covariates give.")
        pairs = {
            t: analysis.tables[f"{t}_pairs"][["grna_target", "response_id"]].drop_duplicates().reset_index(drop=True)
            for t in ("discovery", "positive_control")
        }
        write(
            "parameters",
            {
                "settings": {
                    # in the file, a covariate carries its modality's name
                    "formula": to_file_names(p.formula, main),
                    "side": p.side,
                    "resampling_mechanism": p.resampling_mechanism,
                    "resampling_approximation": approximation,
                    "control_group": p.control_group,
                    "grna_integration_strategy": strategy,
                    "B1": int(p.B1),
                    "B2": int(p.B2),
                    "B3": int(p.B3),
                    "multiple_testing_method": p.multiple_testing_method,
                    "multiple_testing_alpha": float(p.multiple_testing_alpha),
                    "seed": p.seed,
                    "response_fit_method": p.response_fit_method,
                    "grna_fit_method": p.grna_fit_method,
                },
                "design": {"terms": to_file_names(list(design.terms), main), "sha256": design_sha256(design.matrix)},
            },
            {"parameters/discovery_pairs": pairs["discovery"], "parameters/positive_control_pairs": pairs["positive_control"]},
        )

    if "assignment" in done:
        a = analysis.assignment
        initial = a.initial
        held = a.cells_w_zero_or_twoplus_grnas
        flagged = np.empty(0, dtype=np.int64)
        if analysis.low_moi and held is not None:
            # cells the method itself flags, beyond what the assignment table shows
            flagged = np.setdiff1d(held, not_one_grna(initial, len(cells)))
        write(
            "assignment",
            {"settings": {"method": a.method, "hyperparameters": _plain(a.hyperparameters) or {}}},
            {
                "assignment/assignments": pd.DataFrame(
                    {
                        "grna_id": np.repeat(np.array(list(initial), dtype=object), [len(v) for v in initial.values()]),
                        "cell_id": cells[np.concatenate(list(initial.values())).astype(np.int64)],
                    }
                ),
                "assignment/flagged_cells": pd.DataFrame({"cell_id": cells[flagged]}) if flagged.size else None,
            },
        )

    if "qc" in done:
        qc = analysis.qc
        in_use = np.zeros(len(cells), dtype=bool)
        in_use[qc.cells_in_use] = True
        pair_qc = {}
        for pair_type in ("discovery", "positive_control"):
            info = analysis.tables[f"{pair_type}_pairs_with_info"]
            pair_qc[pair_type] = pd.DataFrame(
                {
                    **_unit_columns(info["grna_group"], grna_table, strategy),
                    "response_id": info["response_id"].to_numpy(dtype=object),
                    "pass_qc": info["pass_qc"].to_numpy(dtype=bool),
                }
            )
        write(
            "qc",
            {
                "settings": {
                    "response_n_umis_range": list(qc.response_n_umis_range),
                    "response_n_nonzero_range": list(qc.response_n_nonzero_range),
                    "p_mito_threshold": qc.p_mito_threshold,
                    "pairwise_qc": qc.pairwise_qc,
                    "B3": int(analysis.parameters.B3),
                },
                "summaries": {
                    "cells_removed": {k.removeprefix("n_cells_rm_"): int(v) for k, v in qc.cell_removal_metrics.items()}
                },
            },
            {
                "qc/cells": pd.DataFrame({"cell_id": cells, "in_use": in_use}),
                "qc/discovery_pair_qc": pair_qc["discovery"],
                "qc/positive_control_pair_qc": pair_qc["positive_control"],
            },
        )

    fits = analysis.fits

    def fit_table():
        terms = list(analysis.design.terms)
        rows = []
        for response_id, fit in analysis.response_fits().items():
            if list(fit["fitted_coefs"]) != terms:
                raise ValueError("A saved fit's terms are not the design's.")
            rows.append({"response_id": response_id, "coefficients": list(fit["fitted_coefs"].values()), "theta": fit["theta"]})
        return pd.DataFrame(rows, columns=["response_id", "coefficients", "theta"])

    def result_table(stage):
        result = analysis.results[stage].drop(columns="significant", errors="ignore")
        # stored only where it does not follow from `fold_change` (T5)
        if log_2_follows(result):
            result = result.drop(columns="log_2_fold_change")
        return result.drop(columns=resample_columns(result))

    def resamples(stage):
        return resamples_table(analysis.results[stage])

    def threshold(stage):
        return significance_threshold(analysis.results[stage])

    if "calibration" in done:
        c = analysis.calibration
        groups = list(dict.fromkeys(analysis.tables["negative_control_pairs"]["grna_group"]))
        # today a control group's members exist only in its name (T4)
        members = [g.split("&") for g in groups]
        write(
            "calibration",
            {
                "settings": {
                    "calibration_group_size": int(c.calibration_group_size),
                    "n_calibration_pairs": int(c.n_calibration_pairs),
                    "B3": int(analysis.parameters.B3),
                },
                "significance_threshold": threshold("calibration"),
            },
            {
                "calibration/negative_control_groups": pd.DataFrame(
                    {
                        "grna_group": np.repeat(np.array(groups, dtype=object), [len(m) for m in members]),
                        "grna_id": np.array([g for m in members for g in m], dtype=object),
                    }
                ),
                "calibration/results": result_table("calibration"),
                "calibration/resamples": resamples("calibration"),
            },
        )
    if "power" in done:
        write("power", {"settings": {}}, {"power/results": result_table("power"), "power/resamples": resamples("power")})
    if "discovery" in done:
        write(
            "discovery",
            {"settings": {}, "significance_threshold": threshold("discovery")},
            {"discovery/results": result_table("discovery"), "discovery/resamples": resamples("discovery")},
        )
    # the fits the analysis holds go to the cache, which no stage owns (T6)
    if fits is not None and "qc" in done:
        write_fits(directory, fit_table(), _fit_identity(analysis, mdata))
    return Path(directory)


def compare_analyses(a, b) -> pd.Series:
    """Compare two pysceptre Analysis objects part by part. Tables are compared with
    their rows sorted and per-cell values as sets, since the files keep neither order."""

    def equal(x, y):
        try:
            np.testing.assert_equal(x, y)
            return True
        except AssertionError:
            return False

    cell_slots = set(a._state.cells)
    verdict = {}
    settings_a, settings_b = a.settings, b.settings
    for name in sorted((set(settings_a) | set(settings_b)) - cell_slots):
        verdict[f"slot: {name}"] = equal(settings_a.get(name), settings_b.get(name))
    for name in sorted(cell_slots):
        verdict[f"cells: {name}"] = np.array_equal(np.sort(a._state.cells[name]), np.sort(b._state.cells[name]))
    if a.design is not None or b.design is not None:
        verdict["design matrix"] = a.design is not None and b.design is not None and np.array_equal(a.design.matrix, b.design.matrix)
    for kind, x, y in (("table", a.tables, b.tables), ("result", a.results, b.results)):
        for name in sorted(set(x) | set(y)):
            verdict[f"{kind}: {name}"] = name in x and name in y and same_table(x[name], y[name])
    verdict["fits"] = a.response_fits() == b.response_fits()
    for modality in ("response", "grna"):
        verdict[f"counts: {modality}"] = (a.screen[modality].X != b.screen[modality].X).nnz == 0
    verdict["covariates"] = a.cell_covariates.equals(b.cell_covariates)
    return pd.Series(verdict)


def comparison_summary(verdict: pd.Series) -> str:
    different = list(verdict.index[~verdict.astype(bool)])
    return f"{int(verdict.sum())} of {len(verdict)} parts equal" + (f"; different: {different}" if different else "")


