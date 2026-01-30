# agroclimExtremes

**R workflow** associated to the manuscript "Global hotspots of co-occurring extreme droughts in agriculture". It contains the R scripts ordered ordered by step to compute extreme agricultural droughts and their spatial dependence over agricultural production areas worldwide at 0.25° resolution.

> **Status:** research code (intended for reproducible analysis accompanying a manuscript / report).  
> **License:** GPL-3.0

---

## What this repository does

This repository brings together scripts to:

- **Compute agroclimatic indices** (e.g., drought-related indices) from climate inputs.
- **Identify extreme events** based on user-defined thresholds and time windows.
- **Summarize and export gridded outputs** for downstream analysis and visualization.
- **Explore results interactively** via Leaflet and Shiny maps.

---

## Repository structure

- `R/`  
  Main analysis pipeline. Scripts are organized **by steps** (run in sequence). Each step reads the outputs of the previous step and writes intermediate results to disk.

- `data/`  
  Example data and/or lightweight inputs used by the scripts. (Large raw datasets are typically not stored in the repository.)

- `interactive/`  
  Resources used by the interactive components (e.g., Shiny app assets, helper files).

- Top-level scripts:
  - `map_leaflet.R` — lightweight Leaflet map (standalone).
  - `map_shiny.R` — Shiny app for interactive exploration.
  - `spatial_dependence_simulation_example.R` — example / sandbox for spatial-dependence simulations.

---

## How to run the pipeline (R scripts by steps)

1. **Clone the repository**
   ```bash
   git clone https://github.com/haachicanoy/agroclimExtremes.git
   cd agroclimExtremes
   ```

2. **Open an R session in the project root** (RStudio recommended).

3. **Run the scripts in `R/` in step order**  
   The scripts in `R/` are numbered / grouped by steps. Run them sequentially (Step 0 → Step N).  
   Each script should document:
   - required inputs
   - parameters to set (paths, years, thresholds, etc.)
   - outputs written to disk

> Tip: if you’re setting this up for publication, keep your configuration (paths, thresholds, time window, scenario choices) in a single “settings” section or config file to make reruns easy.

---

## Running the interactive maps

### Leaflet (standalone)
In R:
```r
source("map_leaflet.R")
```

### Shiny app
In R:
```r
source("map_shiny.R")
```
If the script defines a `shinyApp()` object, it should launch automatically; otherwise, follow the instructions inside the script.

---

## Dependencies

This repository is written in **R**. Package requirements depend on which step(s) you run, but commonly include packages for:

- data wrangling (e.g., `dplyr`, `data.table`)
- spatial data (e.g., `sf`, `terra`/`raster`)
- interactive mapping (e.g., `leaflet`, `shiny`)
- plotting (e.g., `ggplot2`)

If you want a quick “install everything” approach:
```r
pkgs <- c("sf","terra","raster","dplyr","data.table","ggplot2","leaflet","shiny")
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

---

## Reproducibility (recommended for journals / Code Ocean)

If you are packaging this repository for **Code Ocean** (or any containerized environment):

- Pin the **R version** and the **package versions** used in the analysis.
- Prefer a lockfile workflow (e.g., `renv`) so the capsule can restore dependencies deterministically.
- Keep large climate datasets outside the repo and document how they are obtained (download links, preprocessing, checksums).

A minimal, publication-friendly capsule structure is:

- `R/` (analysis steps)
- `data/` (small example inputs; *or* a script that downloads inputs)
- `results/` (generated outputs)
- `run_all.R` (optional driver script to run steps end-to-end)

---

## Citation

If you use this code in academic work, please cite the associated manuscript / software record.

**Suggested citation (update as needed):**
> Achicanoy, H.; Ramirez Villegas, J.; Meuwissen, M.; Dalhaus, T. *Interactive platform for mapping extreme drought in agricultural systems*.

---

## Contributing

Issues and pull requests are welcome. For substantial changes, please open an issue first to discuss scope and approach.

---

## Acknowledgements

This repository may build upon concepts and methods from existing open-source tools and published literature. Please see the manuscript and code comments for method-specific attributions.