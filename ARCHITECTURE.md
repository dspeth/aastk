# AASTK Architecture

This document captures the design decisions and conventions that contributors should follow.
For a user-facing overview of the tools and their scientific purpose, see [README.md](README.md).

---

## Table of contents

1. [Repository layout](#repository-layout)
2. [Module responsibilities](#module-responsibilities)
3. [Core design patterns](#core-design-patterns)
4. [CLI conventions](#cli-conventions)
5. [Adding a new tool](#adding-a-new-tool)
6. [Versioning and releases](#versioning-and-releases)

---

## Repository layout

```
aastk/
├── aastk/
│   ├── __init__.py       # empty; marks the package
│   ├── __main__.py       # entry point: argument dispatch, print_help()
│   ├── cli.py            # argparse definitions for every subcommand
│   ├── log.py            # logger_setup() – dual console/file logging
│   ├── util.py           # shared helpers used by all modules
│   ├── version.py        # __version__, __author__, __copyright__
│   ├── database.py       # SQLite schema, database construction, meta, export_fasta
│   ├── pasr.py           # PASR workflow and subcommand functions
│   ├── casm.py           # CASM workflow and subcommand functions
│   └── cugo.py           # CUGO workflow and subcommand functions
├── pyproject.toml        # build config, dependencies, entry point
├── meta.yaml             # conda recipe (bioconda)
└── .github/workflows/
    └── release.yml       # automated PyPI publish on GitHub release
```

---

## Module responsibilities

| Module | Owns |
|---|---|
| `__main__.py` | `main()` dispatcher, `print_help()` |
| `cli.py` | All `argparse` definitions; no business logic |
| `log.py` | `logger_setup()` only |
| `util.py` | Functions shared by ≥ 2 tool modules |
| `database.py` | SQLite schema and population; `meta()` for metadata retrieval; `export_fasta()`. The database stores protein sequences (zlib-compressed), genomic context (CUGO assignments, strand, position), functional annotations (COG, KEGG, Pfam), transmembrane helix counts, and genome-level metadata (taxonomy, environment, culture collection). |
| `pasr.py` | Protein Alignment Score Ratio: identifies homologous proteins in a large query dataset by aligning against a small seed set with DIAMOND and computing a score ratio relative to each sequence's theoretical maximum BLOSUM score. |
| `casm.py` | Clustering Alignment Score Matrix: characterizes the structure of a protein superfamily by building an N×n alignment score matrix, reducing it to 2D with t-SNE (via openTSNE), and calling clusters with DBSCAN. |
| `cugo.py` | Colocated Unidirectional Gene Organization: retrieves the genomic neighborhood of each query protein from the database (bounded by strand changes or contig ends) and visualizes consensus context across annotation, sequence length, and transmembrane helix dimensions. |

The separation between `cli.py` and the tool modules is intentional: `cli.py` contains
only argument definitions, and tool modules contain only logic. This makes it possible to
call any subcommand function programmatically without touching argparse.

---

## Core design patterns

### Dual-mode functions

Each tool is implemented as a set of individually callable subcommand functions plus a
top-level orchestrator that chains them. For example, `pasr.py` exposes `build()`,
`search()`, `get_hit_seqs()`, `max_score()`, `bsr()`, and `pasr_plot()` as standalone
functions, and `pasr()` calls them in sequence.

This means a user can run the full workflow with `aastk pasr ...` or run individual steps
separately—useful for debugging and for resuming a partially complete run.

### Output path management

All functions that write files call `ensure_path()` from `util.py` before writing.
`ensure_path()` creates parent directories as needed, raises `FileExistsError` if a file
would be overwritten without `force=True`, and returns the resolved path string. Do not
replicate this logic inline.

### External tool invocation

DIAMOND and SeqKit are invoked via `subprocess`. Always call
`check_dependency_availability(command)` from `util.py` before the subprocess call; it
raises `FileNotFoundError` with a clear message if the binary is missing. Subprocess
`stderr` is captured and written to the subprocess log file at level 99
(see [Logging](#logging) below).

### Parallel processing

CPU-bound batch work (e.g., batch database fetches, sequence retrieval) uses
`concurrent.futures.ThreadPoolExecutor`. Pass `threads` through from the CLI argument
all the way to the executor. Long-running operations expose progress via `tqdm`.

### Database access

The SQLite database is opened with `sqlite3.connect()`. Sequences are stored zlib-compressed
(`compress_sequence` / `decompress_sequence` in `util.py`). Queries that touch many rows
use batch fetching (typically 500–900 IDs per query) to bound memory usage. Never read
the entire `protein_data` table into memory at once—use `stream_all_proteins()` from
`database.py` when iterating over all sequences.

### Logging

`logger_setup()` in `log.py` configures a root logger with two handlers:

- **Console handler** — `INFO` level, written to `stdout`. Use for user-visible progress
  messages (`logger.info(...)`).
- **File handler** — level 99 (custom `SUBPROCESS` pseudo-level), written to a
  timestamped `aastk-subprocesses-<timestamp>.log` next to the output directory. Use for
  subprocess stderr (`logger.log(99, ...)`).

Pass `--silent` to suppress console output entirely (sets root level to `CRITICAL + 1`).
Do not use `print()` for output that belongs in the log.

### Visualisation

All plots use `matplotlib`. Functions that produce plots accept an `svg: bool` parameter
and write both `.png` and `.svg` when it is `True`, otherwise only `.png`. Colour mapping
for categorical data (clusters, annotations) is determined at plot time from the data; do
not hard-code colour lists.

---

## CLI conventions

### Argument helpers

Every reusable argument is defined as a private module-level function in `cli.py` with the
naming convention `__<argument_name>(group, required=False)`. The function calls
`group.add_argument(...)` and returns nothing. This keeps subparser definitions concise
and ensures that shared arguments (e.g. `--output`, `--threads`, `--force`) have
consistent flag letters and help text across all subcommands.

When adding a new argument that will appear in more than one subcommand, add a helper
function. When an argument is specific to a single subcommand, it is acceptable to define
it inline in the subparser block.

### Subparser structure

Use the three context-manager helpers defined at the top of `cli.py`:

```python
with subparser(sub_parsers, 'my_cmd', 'Short description') as parser:
    with arg_group(parser, 'Required arguments') as grp:
        __fasta(grp, required=True)
    with arg_group(parser, 'Optional') as grp:
        __output(grp)
        __force(grp)
```

Use `mutex_group(parser, required=True/False)` when two arguments are mutually exclusive:

```python
with mutex_group(parser, required=True) as grp:
    __fasta(grp)
    __id_list(grp)
```

All subparsers should include `__force(grp)` in the optional group whenever the subcommand
writes files, and `__output(grp)` unless the output location is fixed by another argument.

---

## Adding a new tool

Follow these steps to integrate a new top-level tool (e.g. `mytool`).

### 1. Create the module

Create `aastk/mytool.py`. Structure it like the existing tool modules:

- One function per subcommand step (e.g. `step_a()`, `step_b()`).
- One top-level orchestrator (e.g. `mytool()`) that calls the steps in sequence.
- All functions accept explicit keyword arguments—no `argparse.Namespace` objects inside
  the module. Keep argparse in `cli.py` and `__main__.py`.
- Use `ensure_path()` for all file output.
- Use `logger.info()` for progress messages; `logger.error()` for failures.
- Raise `RuntimeError` for unrecoverable errors (caught by the top-level handler in
  `__main__.py`).

### 2. Add argument helpers to `cli.py`

Add a `__myarg(group, required=False)` function for every new argument that is not already
covered by an existing helper. Place new helpers in alphabetical order among the existing
ones.

### 3. Register subparsers in `cli.py`

Add subparser blocks inside `get_main_parser()`, grouped with a `### PARSER FOR MYTOOL ###`
comment:

```python
### PARSER FOR MYTOOL FUNCTIONALITIES AND WORKFLOW ###
with subparser(sub_parsers, 'step_a', 'Description of step A') as parser:
    with arg_group(parser, 'Required arguments') as grp:
        __fasta(grp, required=True)
    with arg_group(parser, 'Optional') as grp:
        __output(grp)
        __threads(grp)
        __force(grp)

with subparser(sub_parsers, 'mytool', 'MYTOOL: short description') as parser:
    ...
```

### 4. Register dispatch cases in `__main__.py`

Add the import at the top of `__main__.py`:

```python
from aastk.mytool import *
```

Add `elif` branches in `main()` for each new subcommand, grouped under a matching
`### PARSER FOR MYTOOL ###` comment:

```python
### PARSER FOR MYTOOL FUNCTIONALITIES AND WORKFLOW ###
elif args.subparser_name == 'step_a':
    step_a(
        fasta=args.fasta,
        output_dir=args.output,
        threads=args.threads,
        force=args.force
    )

elif args.subparser_name == 'mytool':
    mytool(
        fasta=args.fasta,
        output_dir=args.output,
        threads=args.threads,
        force=args.force
    )
```

Call functions with explicit keyword arguments, never by positional order.

### 5. Update `print_help()` in `__main__.py`

Add the new tool and its helper subcommands to the help text under the appropriate
section (Main Tools, Helper tools, or Subcommands).

---

## Versioning and releases

The single source of truth for the version is `aastk/version.py`. Update `__version__`
there before tagging a release. The version is pulled into `pyproject.toml` dynamically
via `[tool.setuptools.dynamic]` and into the `meta.yaml` conda recipe manually.

Releases are published automatically by `.github/workflows/release.yml` when a GitHub
Release is created. The workflow publishes to TestPyPI first (optional step), then to
PyPI, using OIDC trusted publishing (no stored API tokens required).

The conda package on bioconda is maintained separately via the bioconda-recipes repository
and uses `meta.yaml` as the recipe source.
