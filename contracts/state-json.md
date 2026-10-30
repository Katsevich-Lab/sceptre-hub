# Analysis files: state JSON, format stamp and execution plans

R sceptre and pysceptre save an analysis to the same directory layout, and each
reads what the other writes. This contract fixes the parts of that layout that
are not a standard format on their own: the analysis state in `state.h5ad`,
the file format stamp, and `execution.json`. It describes **file format 2**.

The implementations are `R/state_json.R` and `R/execution_plans.R` in sceptre,
and `src/pysceptre/io/state_json.py` in pysceptre. Both test the same pinned
state text (sceptre `tests/testthat/test-analysis-objects.R`, "the state is
plain JSON, one key per value"; pysceptre `tests/validation/test_state_json.py`),
so for the same state they write the same bytes.

## The directory

| File | Holds | Format stamp |
| --- | --- | --- |
| `screen.h5mu` | response and gRNA counts, cell covariates, MOI (MuData; modalities `response` and `grna`) | `uns["sceptre_format"]` |
| `grna_targets.parquet` | the gRNA table | Parquet schema metadata `sceptre_format` |
| `state.h5ad` | per-cell stage values in `obs`, the design matrix in `obsm`, every other value in `uns["state"]` (below) | `uns["sceptre_format"]` |
| `tables/NAME.parquet` | pair tables and other tables of the analysis | schema metadata `sceptre_format` |
| `results/NAME.parquet` | calibration, power and discovery results | schema metadata `sceptre_format` |
| `fits.parquet` | response model fits, one row per response and term | schema metadata `sceptre_format` |
| `execution.json` | how each result was computed (below) | top-level key `sceptre_format` |
| `checksums.tsv` | size and MD5 of every other file | none (a manifest) |

**Shared screens.** An analysis whose response counts are read from disk
(backed) is saved without copying the counts:

- **In place.** Saved to the directory it was read from, or that holds its
  `screen.h5mu`, it replaces that directory's analysis files (`state.h5ad`,
  `tables/`, `results/`, `fits.parquet`, `execution.json`) and leaves
  `screen.h5mu` and `grna_targets.parquet` untouched. The new files are written
  to a hidden staging folder in the directory and renamed into place, so a
  failed save leaves the previous analysis readable, and no table or result of
  the previous save survives. `checksums.tsv` is written last, under another
  name, then renamed.
- **By reference.** Saved to another directory, it writes the analysis files
  and no screen: `uns["state"]` key `screen_reference` (below) says where the
  screen is.
- In both cases the analysis's gRNA table (which `set_analysis_parameters()`
  extends with `grna_group`) is `tables/grna_target_data_frame.parquet`, and
  readers use it in place of the screen's `grna_targets.parquet`.
- A self-contained copy (`self_contained` in both packages) writes the screen
  into the directory. An analysis held in memory is always written with its
  screen.

Readers resolve `screen_reference` to the directory given relative to the
analysis directory, else to the absolute one; the files there must have the
recorded sizes (and, when the reader verifies, MD5s). Otherwise the read stops
with an error that says where the screen was and how to put it back.

## The format stamp

Every file but `checksums.tsv` carries the same stamp, a JSON text:

```json
{"format":2,"readers":{"pysceptre":"0.6.0","sceptre":"0.99.1"}}
```

`format` is the file format. `readers` gives, for each package, the first
release that reads it.

- A file without a stamp is **format 1**: written before the stamp existed.
  Its `uns["state"]` is R's `jsonlite::serializeJSON()` text. Both packages
  still read format 1 (sceptre with `jsonlite::unserializeJSON()`, pysceptre
  with `pysceptre.io.r_json`), and write format 2.
- A reader refuses a file whose format is newer than it reads, and says which
  release reads it, from the stamp's `readers`:

  > state.h5ad was written in sceptre file format 3; this version of
  > pysceptre (0.6.0) reads formats up to 2. Upgrade to pysceptre >= 0.7.0 to
  > read it.

- Within a format it reads, a reader stays strict: an unknown state key is an
  error, never dropped silently.

**Changing the format.** Adding, removing or retyping a state key, or changing
how a file is laid out, is a new format. Bump `format` in both packages
together, update both readers, keep reading the older formats, and set
`readers` to the first releases that read the new format. Changes to
`execution.json` plans that readers ignore (new keys inside a plan) are not.

## `uns["state"]`: plain JSON, one key per value

`state.h5ad` is a standard h5ad (AnnData) file without features; sceptre writes
it with `rhdf5`, following the AnnData on-disk format, and pysceptre with
`anndata`.

`uns["state"]` is one JSON object (a string in `state.h5ad`), with one key per
value of the analysis state, in the order of the table below. It holds every
value that has no other home: per-cell values that are a cell mask or one value
per cell are columns of `obs`, the design matrix is `obsm["covariate_matrix"]`,
and tables are Parquet files.

### Encoding

Values keep the form of R values:

| R value | JSON |
| --- | --- |
| `NULL` | `null` |
| atomic vector, unnamed | an array of its elements, of any length (a length-one vector is a one-element array; an empty vector is `[]`) |
| atomic vector, named | an object, one member per element, in order |
| missing element (`NA`) | `null` |
| logical element | `true`, `false` |
| integer element | a JSON integer |
| double element | a number with 17 significant digits (C `%.17g`), always containing `.` or an exponent (`2.0`, not `2`); `NaN`, `Inf` and `-Inf` are the strings `"NaN"`, `"Inf"`, `"-Inf"` |
| character element | a JSON string, escaped as Python's `json.dumps(ensure_ascii=False)` escapes it (`\"`, `\\`, `\b`, `\f`, `\n`, `\r`, `\t`, other control characters as `\u00XX` in lowercase hex; non-ASCII characters as UTF-8) |
| list, unnamed | an array of its encoded elements |
| list, named | an object of its encoded elements, in order |
| formula | its text (R's `deparse()`), a bare JSON string; `""` is the empty formula |

Writers emit no whitespace. Names of named values are unique and non-empty.
Other R attributes (dimensions, classes other than formula, factor levels)
cannot be stored, and writing them is an error.

### Types

Each key has a type, which a reader uses where the JSON alone is ambiguous:
the type of an empty or all-missing vector, and of the vectors inside the gRNA
assignment lists.

- **logical, integer, double, character:** an atomic vector of that type.
  Writers write that type; readers refuse another.
- **numeric:** integer or double, as the elements are written.
- **formula:** a formula.
- **integer_list:** a list (named or not) of integer vectors.
- **assignments:** a named list whose members `grna_group_idxs` and
  `indiv_nt_grna_idxs` are integer_lists and whose member `all_nt_idxs` is an
  integer vector; other members are free-form.
- **any:** free-form. An array or object whose elements are all scalars or
  `null` is an atomic vector: double if any element is a non-integer number,
  else character if any is a string, else integer if any is an integer, else
  logical. Any other array or object is a list; a bare string is a formula.
  Because the JSON alone decides, writers check that the value reads back
  exactly and refuse it otherwise (for example an integer vector of only
  missing values, which would read back as logical).

### Keys

In this order. "Empty" is the value before the step that sets it; per-cell
keys are present here only when they cannot be an `obs` column (empty, or an
index list that is not strictly increasing).

| Key | Type | Value |
| --- | --- | --- |
| `covariate_matrix_attributes` | any | the design matrix's column names (`terms`), `assign` and `contrasts` (or `null`); `null` before `set_analysis_parameters()` builds the design matrix |
| `covariate_names` | character | names of the cell covariates |
| `formula_object` | formula | the model formula |
| `side_code` | integer | test side: -1 left, 0 both, 1 right |
| `resampling_approximation` | character | `"spa"`, `"skew_normal"` or `"no_approximation"` |
| `control_group_complement` | logical | whether the control group is the complement set |
| `run_permutations` | logical | permutations (true) or the conditional randomization test |
| `B1`, `B2`, `B3` | integer | resampling budgets |
| `grna_integration_strategy` | character | `"union"`, `"singleton"` or `"bonferroni"` |
| `grna_assignment_method` | character | the gRNA assignment method |
| `grna_assignment_hyperparameters` | any | the method's hyperparameters (named list; a mixture's includes `formula_object`) |
| `multiple_testing_alpha` | numeric | the level of the multiple-testing correction and of pairwise QC |
| `multiple_testing_method` | character | the `p.adjust()` method |
| `cell_removal_metrics` | integer | named: cells removed by each cellwise QC rule |
| `cellwise_qc_thresholds` | any | named list of the cellwise QC thresholds |
| `mitochondrial_gene` | logical | unused; kept empty |
| `n_nonzero_tot_vector` | integer | named by response: cells expressing it among the cells of the test |
| `initial_grna_assignment_list` | integer_list | named by gRNA: the cells assigned each gRNA (1-based) |
| `grna_assignments_raw` | assignments | cells of each gRNA group and non-targeting gRNA, before cellwise QC |
| `grna_assignments` | assignments | the same after cellwise QC, and `all_nt_idxs` |
| `grnas_per_cell` | integer | per cell (usually an `obs` column) |
| `cells_w_zero_or_twoplus_grnas` | integer | 1-based cell indices (usually an `obs` mask) |
| `cells_in_use` | integer | 1-based cell indices (usually an `obs` mask) |
| `n_discovery_pairs`, `n_ok_discovery_pairs` | integer | numbers of discovery pairs, and of those passing pairwise QC |
| `n_positive_control_pairs`, `n_ok_positive_control_pairs` | integer | the same for positive-control pairs |
| `calibration_group_size` | integer | non-targeting gRNAs per undercover group |
| `n_calibration_pairs` | integer | negative-control pairs of the calibration check |
| `mean_cells_per_grna` | numeric | mean cells per gRNA |
| `last_function_called` | character | the last step run |
| `functs_called` | logical | named by step: whether it has run |
| `nuclear`, `nf_pipeline` | logical | flags of the old on-disk pipeline |
| `integer_id` | integer | kept from the old on-disk pipeline |
| `elements_to_analyze` | character | kept from the old on-disk pipeline |
| `response_fit_method` | character | the response model fitter; empty in files that predate it |
| `grna_fit_method` | character | the gRNA model fitter; empty in files that predate it |
| `seed` | integer | the resampling seed; empty means the default (4) |
| `pairwise_qc` | logical | whether pairwise QC (testability) is on; empty means on |
| `screen_reference` | any | `null` when the directory holds its screen; otherwise the shared screen: `directory` (absolute path), `relative_directory` (relative to the analysis directory) and `files`, an object with, for `screen.h5mu` and `grna_targets.parquet`, `size` (bytes, a double) and `md5` |

### Example

```json
{"covariate_matrix_attributes":null,"covariate_names":["batch","grna_n_nonzero"],"formula_object":"","side_code":[],"B1":[],"multiple_testing_alpha":[],"functs_called":{"import_data":true,"set_analysis_parameters":false},"elements_to_analyze":[null],"seed":[]}
```

(abridged: a real state holds every key of the table.)

## `execution.json`: how each result was computed

`execution.json` is run provenance. It never changes a result, and results are
identical whatever it says. It is a JSON object with the key `sceptre_format`
(the stamp, as a JSON object) and one key per result in `results/`
(`calibration`, `power`, `discovery`) whose value is the **plan** of the run
that computed it. A directory whose results were all computed before plans were
recorded has no `execution.json`.

A plan is a JSON object with these keys, each a JSON object:

| Key | Holds |
| --- | --- |
| `software` | `package` (`"sceptre"` or `"pysceptre"`) and `version` |
| `options` | the execution options the run was given, in that package's own names (R: `parallel`, `n_processors`; pysceptre: `ExecutionOptions`) |
| `engine` | `score_backend`, `storage` (`"in memory"` or `"backed"`), `output_amount` where known, and the package's fitter settings |
| `host` | the layout, workers and chunk or tile sizes, with reasons where the package records them |
| `resampling_rng` | `rng_protocol`, `resampling_seed`, `description` (pysceptre adds the protocol's hashes) |

Packages may add keys (sceptre adds `timings`; pysceptre adds `streaming`,
`device` and `scores` for the runs that use them). Readers ignore keys they do
not know inside a plan, keep the plans of results they read, and replace a
result's plan when they compute that result again. `execution.json` is read
only when `checksums.tsv` lists it.
