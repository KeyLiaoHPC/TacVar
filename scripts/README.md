# TacVar visualization scripts

Python helpers for inspecting TacVar measurement CSVs from NPB-MPI (and related suites that share the same CSV layout).

## Dependencies

- Python 3
- `numpy`, `pandas`, `matplotlib`

Open notebooks from this directory so `import tacvar_visual` resolves:

```bash
cd scripts
jupyter notebook visualize_npb-mpi.ipynb
```

## CSV layout

Measurement runs write:

```text
<data_root>/timer_info.csv
<data_root>/<kernel_name>.<kernel_class>/<short_host>_rRRRR_tTTTT_pPID.csv
```

See root `README.md` §3.5 for the full column schema. `timer_info.csv` is the run catalog (`suite,benchmark,class,test_tag,timer,region_id,nloc,name`); loaders filter rows to the requested kernel. Heatmaps and distribution plots use `region_id`, `loc_id`, `rank`, and the chosen figure-of-merit column (default `elapsed_ns`).

## API: `tacvar_visual.draw_heatmap_npb_mpi`

Draw one discrete green–white–red heatmap per `(region_id, loc_id)` pair that exists in the data. `data_root` may be a single `data_*` path or a list of paths; each root is plotted independently (figures concatenated in list order). With two or more roots, the title includes the data-directory basename.

- **X axis**: sample index in CSV row order (after filtering to that pair). Tick labels keep the original sample indices when `xrange` selects a subset.
- **Y axis**: MPI `rank` (ascending; rank 0 at the top).
- **Color**: discrete bin of `value / ref` — more negative vs ref → deeper green; bin containing 0% → white; more positive → deeper red.
- **Figure height**: 1024 pixels (`dpi=100`, height 10.24 in). Academic serif style (Cambria with serif fallbacks).

Before loading rank CSVs the function prints `timer_info.csv`. Before each figure it prints the rank count and data-point count for that `(region_id, loc_id)`.

Title pattern:

```text
{fom} heatmap of {Kernel.Class} on {host} (region={region_name}, loc_id={id})
```

Remaining arguments (`hosts`, `region_ids`, `loc_ids`, `fom`, …) are **keyword-only**. After changing the module, re-run `importlib.reload` or restart the Jupyter kernel.

```python
import importlib
import tacvar_visual as tvvis
importlib.reload(tvvis)

figs = tvvis.draw_heatmap_npb_mpi(
    data_root="path/to/data_YYYYMMDDTHHmmss",
    kernel_name="is",
    kernel_class="S",
    hosts=["all"],
    region_ids=range(0),    # empty → all regions
    loc_ids=range(0),       # empty → all loc_ids
    fom="elapsed_ns",
    xrange=-1,
    ref_section_base=True,
    ref_key="median",
    ref_base="all",
    bars=[-15, -5, -1, 1, 5, 15],
)
```

| Argument | Default | Meaning |
|----------|---------|---------|
| `data_root` | (required) | Path to a `data_*` directory, or a list of such paths (one figure set per root) |
| `kernel_name` | (required) | Benchmark name, e.g. `is`, `cg` |
| `kernel_class` | (required) | NPB class letter, e.g. `S` |
| `hosts` | `["all"]` | `["all"]` / `[""]` / empty → all hosts; otherwise match CSV `short_host` prefixes |
| `region_ids` | `range(0)` | Region ids; empty iterable → all regions found |
| `loc_ids` | `range(0)` | Loc ids; empty iterable → all loc_ids found; only pairs present in data are plotted |
| `fom` | `"elapsed_ns"` | Numeric CSV column (`elapsed_ns` or a counter `*_delta`) |
| `xrange` | `-1` | Columns to draw (see below) |
| `ref_section_base` | `True` | If True, compute refs on the `xrange` window; if False, on all samples (display still uses `xrange`) |
| `ref_key` | `"median"` | Reference statistic for `value / ref` (see below) |
| `ref_base` | `"all"` | How that sample set is aggregated: `all` / `x` / `y` (see below) |
| `bars` | `[-15, -5, -1, 1, 5, 15]` | Percent-deviation cut points (see bins below) |

Returns a `list` of `matplotlib.figure.Figure`.

### `xrange`

Selects sample columns (X) to draw. Dispatch is by type: a `list` is fancy indexing, a `tuple` is a slice.

- `-1` — all columns (`0 .. n-1`).
- `list` of integers — those indices, in the given order. Python-style negatives are allowed. Empty list or out-of-range index raises `ValueError`.
- `tuple (a, b)` of length 2:
  - both `int` and `a <= b`: half-open slice `[a, b)` (Python clipping/negatives, like `arr[a:b]`).
  - both `float` with `0 <= a <= b <= 1`: `[int(n*a), int(n*b))`. `b=1.0` is allowed so `(0.0, 1.0)` is the full range.
  - mixed types, `a > b`, or an empty selection raise `ValueError`.

Integer `(0, 1)` is only column `0`. Use `(0.0, 1.0)` for a fractional full range.

### `ref_key`

Statistic applied to the sample set chosen by `ref_section_base` and `ref_base`:

- `max`, `min`, `median`, `average` (`mean` is accepted as an alias of `average`)
- `pxx` — percentile with `xx` in `[0, 100]`, e.g. `p50`, `p90`, `p99` (`numpy.percentile`, linear interpolation)

Raises `ValueError` if the key is unknown, there are no samples, or `ref == 0`.

### `ref_section_base` and `ref_base`

`ref_section_base` chooses which columns feed the statistic. Display always crops to `xrange`. When `xrange=-1` the two settings match.

- `True` — refs from the `xrange` window only.
- `False` — refs from the full rank-by-sample matrix.

`ref_base` then aggregates that sample set for each drawn cell `(rank i, sample j)`:

- `all` — one ref over all finite values in the chosen sample set.
- `x` — ref from column `j` (all ranks of this seq). Unaffected by `ref_section_base`, because a column is either drawn or not.
- `y` — ref from row `i`. With `ref_section_base=True`, only xrange seqs of this rank; with `False`, all seqs of this rank.

All-NaN row/column cells stay NaN (gray). Unknown `ref_base` raises `ValueError`.

### `bars`

`bars` are **percent deviation** cut points. They must be strictly increasing. Color is discrete: number of bins is `len(bars) + 1` after wrapping with `±inf`.

1. Convert each bar `p` to a ratio edge `1 + p/100`.
2. Add the two infinite ends.

Example for `bars=[a, b, c]` (after percent→ratio conversion):

1. `(-inf, a]`
2. `(a, b]`
3. `(b, c]`
4. `(c, +inf)`

The bin that contains ratio `1.0` (zero deviation vs ref) is **always white**. Left of white: deep green → light green; right of white: light red → deep red.

Default `bars=[-15, -5, -1, 1, 5, 15]` therefore uses **7 bins**, with white = `(-1%, 1%]`.

## API: pooled distribution plots

`draw_histogram_npb_mpi`, `draw_pdf_npb_mpi`, and `draw_cdf_npb_mpi` share the heatmap data path (`hosts`, `region_ids`, `loc_ids`, `fom`, `xrange`, `ref_section_base`, `ref_key`, `ref_base`, `bars`) and add `ranks`. `data_root` may be a path or a list of paths. Histogram always emits one figure per `(region_id, loc_id)` per root. PDF and CDF do the same unless `overlap=True`. Selected ranks are **pooled**; there are no per-rank series.

- **Histogram X**: percent-vs-ref bins from `bars` (same quantity as the heatmap colorbar).
- **PDF/CDF X**: bottom axis is the FOM column (`elapsed_ns`, …); top axis is percent-vs-ref `bars`, mapped with a single pooled ref (`ref_key` on the `ref_section_base` window). GWR bands and dashed guides are drawn in FOM coordinates. `overlap=True` overlays listed roots on one axes per `(region_id, loc_id)` and **ignores** `bars` (no twin axis, no GWR).
- **Figure height**: 1024 pixels (`dpi=100`, height 10.24 in). Academic serif style.
- **Prints**: `timer_info.csv` once per root, then rank count and data-point count per `(region_id, loc_id)` after `ranks` and `xrange`.

Title pattern (`kind` is `histogram`, `pdf`, or `cdf`):

```text
{fom} {kind} of {Kernel.Class} on {host} (region={region_name}, loc_id={id})
```

Remaining arguments are **keyword-only**. After changing the module, re-run `importlib.reload` or restart the Jupyter kernel.

```python
import importlib
import tacvar_visual as tvvis
importlib.reload(tvvis)

common = dict(
    data_root="path/to/data_YYYYMMDDTHHmmss",
    kernel_name="is",
    kernel_class="S",
    hosts=["all"],
    region_ids=range(0),
    loc_ids=range(0),
    ranks=range(0),          # empty → all ranks
    fom="elapsed_ns",
    xrange=-1,
    ref_section_base=True,
    ref_key="median",
    ref_base="all",
    bars=[-15, -5, -1, 1, 5, 15],
)
hfigs = tvvis.draw_histogram_npb_mpi(**common)
pfigs = tvvis.draw_pdf_npb_mpi(**common, overlap=False)
cfigs = tvvis.draw_cdf_npb_mpi(**common, overlap=False)
```

| Argument | Default | Meaning |
|----------|---------|---------|
| `data_root` | (required) | Path to a `data_*` directory, or a list of such paths |
| `kernel_name` | (required) | Benchmark name, e.g. `is`, `cg` |
| `kernel_class` | (required) | NPB class letter, e.g. `S` |
| `hosts` | `["all"]` | `["all"]` / `[""]` / empty → all hosts; otherwise match CSV `short_host` prefixes |
| `region_ids` | `range(0)` | Region ids; empty iterable → all regions found |
| `loc_ids` | `range(0)` | Loc ids; empty iterable → all loc_ids found; only pairs present in data are plotted |
| `ranks` | `range(0)` | MPI ranks to pool; empty iterable → all ranks found |
| `fom` | `"elapsed_ns"` | Numeric CSV column (`elapsed_ns` or a counter `*_delta`) |
| `xrange` | `-1` | Sample columns to include (same rules as the heatmap) |
| `ref_section_base` | `True` | If True, compute refs on the `xrange` window; if False, on all samples |
| `ref_key` | `"median"` | Reference statistic for `value / ref` (same as heatmap) |
| `ref_base` | `"all"` | How that sample set is aggregated: `all` / `x` / `y` (per cell, then pooled). PDF/CDF twin axis always uses one pooled ref |
| `bars` | `[-15, -5, -1, 1, 5, 15]` | Percent-deviation cut points. PDF/CDF: top-axis ticks. Ignored when `overlap=True` |
| `overlap` | `False` | PDF/CDF only. If True, overlay listed roots on one axes per `(region_id, loc_id)` and ignore `bars` |

Returns a `list` of `matplotlib.figure.Figure`. `xrange`, `ref_key`, `ref_section_base`, `ref_base`, and `bars` follow the heatmap sections above. `ref_base='y'` still forms per-rank refs for each cell, then those ratios are flattened (histogram only). PDF and CDF plot raw FOM values; the percent top axis uses one pooled `ref_key` statistic.

### `ranks`

Empty iterable (`range(0)`, `[]`) means every rank present in the loaded CSVs. A non-empty list is intersected with available ranks (same `select_ids` behavior as `region_ids`). After that filter, histogram bins all finite `value/ref` cells in the `xrange` window; PDF and CDF use the finite FOM values in that window.

### Histogram (`draw_histogram_npb_mpi`)

Discrete counts in the same `bars` bins as the heatmap (`len(bars)+1` after `±inf`). X tick labels are the ratio-edge labels (`-inf`, `-15%`, …, `+inf`). Bar faces use the green–white–red colormap; the bin containing 0% is white. Y = count. No `overlap` argument; a list of `data_root`s yields one figure per root.

### PDF (`draw_pdf_npb_mpi`)

One density histogram (`numpy.histogram(..., density=True)` + stairs) of the FOM column over a finite x-range (data min/max, padded to include the finite `bars` mapped through the pooled ref). Bottom x = FOM; top x = percent-vs-ref `bars`. Light GWR bands and dashed guides are drawn in FOM coordinates. Y = probability density.

With `overlap=True`, all listed roots share one axes per `(region_id, loc_id)` (tab10 colors, light-alpha fill, legend = data-directory basename). `bars` is ignored: no twin axis and no GWR.

### CDF (`draw_cdf_npb_mpi`)

One empirical CDF of the FOM column (`sort` + `arange(1, n+1)/n`, step). Same dual x-axis and `bars` guides as the PDF when `overlap=False`. Y in `[0, 1]`. A dotted horizontal line marks cumulative probability 0.5.

`overlap=True` overlays roots as colored step lines on one axes per `(region_id, loc_id)` and ignores `bars`. The 0.5 guideline remains.
