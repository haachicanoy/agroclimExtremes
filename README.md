# agroclimExtremes

**R workflow** associated to the manuscript "_Global hotspots of co-occurring extreme droughts in agriculture_". It contains the R scripts ordered ordered by step to compute extreme agricultural droughts and their spatial dependence over agricultural production areas worldwide at 0.25° resolution.

> **Status:** research code (intended for reproducible analysis accompanying a manuscript / report).  
> **License:** GPL-3.0

---

## What this repository does

This repository brings together scripts to:

- **Compute the SPEI (agricultural drought index)** from daily weather inputs.
- **Identify extreme events during growing seasons**.
- **Summarize and export gridded outputs** for downstream analysis and visualization.
- [**Explore results interactively**](https://haachicanoy.shinyapps.io/extreme_drought_clusters/) via Leaflet and Shiny maps.

---

## Repository structure

- `R/`  
  Main analysis pipeline. Scripts are organized **by steps** (run in sequence). Each step reads the outputs of the previous step and writes intermediate results to disk.

- `data/`  
  Example data and/or lightweight inputs used by the scripts. (Large raw datasets are not stored in the repository.)

- `interactive/`  
  Resources used by the interactive components (e.g., Shiny app assets, helper files).

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
   - parameters to set (paths, years, etc.)
   - outputs written to disk

> Tip: if you’re setting this up for publication, keep your configuration (paths, thresholds, time window, scenario choices) in a single “settings” section or config file to make reruns easy.

---

## Dependencies

This repository is written in **R**. Package requirements depend on which step(s) you run, but commonly include packages for:

- data wrangling (e.g., `tidyverse`)
- spatial data (e.g., `sf`, `terra`)
- interactive mapping (e.g., `leaflet`, `shiny`)
- plotting (e.g., `ggplot2`)

If you want a quick “install everything” approach:
```r
if(!require(pacman)) {install.packages('pacman'); library(pacman)} else {library(pacman)}
pacman::p_load(tidyverse,sf,terra,leaflet,shiny)
```

## Citation

If you use this code in academic work, please cite the associated manuscript / software record.

**Suggested citation (update as needed):**
> Achicanoy, H.; Ramirez Villegas, J.; Meuwissen, M.; Mehrabi, Z.; Dalhaus, T. *Global hotspots of co-occurring extreme droughts in agriculture*.

---

## Contributing

Issues and pull requests are welcome. For substantial changes, please open an issue first to discuss scope and approach.

---

## Acknowledgements

This repository may build upon concepts and methods from existing open-source tools and published literature. Please see the manuscript and code comments for method-specific attributions.