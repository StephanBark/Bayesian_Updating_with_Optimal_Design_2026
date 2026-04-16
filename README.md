# A Bayesian Updating Framework for Long-term Multi-Environment Trial Data in Plant Breeding - Electronic Appendix (2026)

This repository contains the R / Stan / Gurobi workflow used.
Its purpose is reproducabilty of our paper content as well as a generalized pipeline 
to play around with Bayesian LMM model specification to generate in depth HTML reports for inspection.

Author and maintainer of this repository is Stephan Bark. Feel free to contact me!

It supports three working styles:

1) Basic git pulls and local work in the repository via R / RStudio / VS Code ect.
2) As a VS Code Dev Container for editing and running the R workflow in VS Code.
3) As a separate Docker Compose service for RStudio Server in the browser.

The analysis-ready pre-processed MET data file `0 Data/yield_winter_medium.csv` is included in the
repository, so readers can reproduce the main workflow directly after cloning.
However, note that genotype names are anonymized!

## Project Structure

```
### Infrastructure containing docker engine ------------------------------
├── Dockerfile                 -> Shared container build with dedicated VS Code and RStudio targets
├── .dockerignore              -> Ignored components for container builds
├── .gitignore                 -> Ignore components from beeing pushed to git
├── docker-compose.yml         -> Docker compose option
├── .vscode/                   -> VS code setting.json
├── .devcontainer/             -> Docker devcontainer option
└── cmunrm.ttf                 -> Figture writting style

### Content --------------------------------------------------------------
├── read_in_data.R             -> Main entry point after working style is defined
├── 0 Data/                    -> Analysis-ready input data
├── 1 Functions/               -> User R functions to work with
├── 2 Descriptive Analysis/    -> Just a full Rahman et al. 2023 data describtive plot: NOTE no full data access through here!
├── 3 Frequentist Framework/   -> asreml-based mixed-model workflow
└── 4 Bayesian Framework/      -> rstan-based sequential Bayesian updating workflow
```

### Entry Point

`read_in_data.R` is the normal starting point. It loads packages, sources the
reusable functions in `1 Functions/`, and imports `yield_winter_medium`.

### Frequentist Pipeline

The frequentist workflow is built around the scripts in `3 Frequentist Framework/`
and two reusable helper functions from `1 Functions/`.

1. Run `read_in_data.R`.
2. Fit the publication mixed model in
   `3 Frequentist Framework/Frequentist LMM and Optimal Design Analysis.R`.
3. Use `my_residuals()` from `1 Functions/My Residuals Function.R` to extract
   marginal and conditional residuals for assumption checks.
4. Inspect the residual plots.
5. Use `grid_design()` from `1 Functions/Grid of Optimal Design Function.R` to
   translate fitted variance components into approximate and exact optimal
   design allocations.


### Bayesian Pipeline

The Bayesian workflow is organised around sequential Bayesian updating across cycles 
in `4 Bayesian Framework/` and uses three reusable functions from `1 Functions/`.

1. Run `read_in_data.R`.
2. Fit the sequential Bayesian mixed model with
   `fit_general_multiple_cycles()` from
   `1 Functions/Fit general multiple Cycles Function.R`.
3. Render an HTML summary with
   `render_fit_general_multiple_cycles_report()` from
   `1 Functions/Render Bayesian Cycle Report Function.R`.
4. Reproduce exact paper Bayesian histogram figures from our paper 
   following the scribt.
5. Propagate posterior uncertainty into design recommendations with
   `grid_bayes_design()` from
   `1 Functions/Grid of Optimal Bayes Design Function.R`.


## Software Requirements

### Required

- R >= 4.4
- C++17 toolchain
  - Windows: Rtools
  - macOS: Xcode Command Line Tools
  - Linux: `build-essential`

### Optional License

`Gurobi` -> Required for the `od_MISOCP()`-based optimal design workflows used in this repository -> Required for `grid_design()` and `grid_bayes_design()` R-functions.
`asreml` -> Frequentist mixed-model fitting in the publication scripts -> Required for the Frequentist Pipeline.

## Option 1 - Local Setup

Use this route if you want to work directly in RStudio or VS Code without
Docker.

1. Clone the repository.

   ```bash
   git clone https://github.com/<your-username>/Bayesian_Updating_with_Optimal_Design_2026.git
   cd Bayesian_Updating_with_Optimal_Design_2026
   ```

2. Install the core R packages.

   ```r
   install.packages(c(
     "data.table", "readxl", "writexl", "reshape2", "tidyr", "dplyr",
     "lubridate", "purrr", "stringr", "fitdistrplus", "matrixcalc",
     "viridis", "showtext", "ggplot2", "gridExtra", "patchwork", "GGally",
     "MCMCpack", "coda", "ggmcmc", "SpATS", "OptimalDesign", "MASS",
     "knitr", "htmltools", "base64enc", "rstan"
   ))

   # Optional, if you hold a licence:
   # install.packages("asreml", repos = "https://vsni.r-universe.dev")
   ```

3. If you want Gurobi support, install Gurobi locally and then install the R
   interface shipped with your Gurobi installation.

   ```r
   install.packages(
     file.path(Sys.getenv("GUROBI_HOME"), "R", "gurobi_13.0-1.tar.gz"),
     repos = NULL,
     type = "source"
   )
   ```

4. Run `read_in_data.R`.

5. Execute the scripts in `3 Frequentist Framework/` or `4 Bayesian Framework/`.

## Option 2 - VS Code Dev Container

Use this route if you want the fully prepared environment inside VS Code.

The Dev Container is separated from the RStudio Server setup.
It does not start RStudio Server, does not forward port 8787, and is not meant to
manage Docker from inside the container.

1. Install Docker Desktop.
2. Install the Dev Containers extension in VS Code.
3. Open the repository in VS Code.
4. Choose **Dev Containers: Reopen in Container**.
5. Wait for the first image build, then work inside the container as usual.

Notes:

- If the Docker extension previously showed `Failed to connect. Is Docker installed?`,
   that was likely because the extension was running in the remote Dev Container,
   where `docker` is not available.
- Use the Dev Container for code editing, package installation, and running the
   R scripts in VS Code.
- Use Docker Compose separately when you want an RStudio Server session.

## Option 3 - Docker Compose + RStudio Server

Use this route if you want a browser-based RStudio session without opening VS
Code inside the container.

This workflow is kept separate from the Dev Container on purpose.

The base compose workflow works without Gurobi.

```bash
docker compose up --build -d
```

Then open `http://localhost:8787` in your browser.

To stop the service:

```bash
docker compose down
```

### Optional Gurobi in Docker

If you want Gurobi inside the Docker workflows:

1. Place the Linux installer archive, for example
   `gurobi13.0.1_linux64.tar.gz`, in the project root.
2. Place your `gurobi.lic` licence file in the project root.
3. Uncomment the Gurobi licence mount in `docker-compose.yml`.
4. Rebuild the image.

If you skip these steps, the image still builds, but the `od_MISOCP()`-based
optimal-design scripts in this repository will not run.

### Optional asreml in Docker

`asreml` remains commented out in the Docker image. Please modify with your licence.
