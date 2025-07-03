# SpaCo Docker RStudio Environment

This project provides a Docker-based RStudio environment for running SpaCo scripts and generating paper figures.

## Build the Docker Image

Make sure you're in the same directory as the Dockerfile. Then run:

```
docker build -t <your-image-name> .
```

Replace `<your-image-name>` with a name of your choice.

## Run the Docker Container

From the same directory, run:

```
docker run --rm -ti \
  -v "$(pwd):/home/rstudio/data_dir" \
  -e ROOT=true \
  -e PASSWORD=123 \
  -p 8787:8787 \
  <your-image-name>
```

- This mounts your current directory into the container at `/home/rstudio/data_dir`.
- RStudio Server will be available at `http://localhost:8787`.
- Login with:
  - Username: `rstudio`
  - Password: `123`

**Important:** The `R_scripts` directory must be in the current working directory when you run the container. If it's not, file paths inside the container will not resolve correctly due to how the volume is mounted.

## Directory Guide

- `R_scripts/` — Main analysis scripts. Must be present in your working directory when running the container.
- `fig2_scripts/` — R scripts used to generate the paper figures.
- `SPACO_paper_data/` — Contains saved R objects for benchmarking and original Seurat objects.
- `SpaCo_Kia_R/` — The updated SpaCo package source that works in this Docker environment.
