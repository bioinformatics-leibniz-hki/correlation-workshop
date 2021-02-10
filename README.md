# Correlation workshop

This is the repository for the workshop at [MiCom 2021](https://www.micom.uni-jena.de/) about taxonomic correlation analyses in metagenomics using R.

## Contents

- [Raw data directory](data) with example data and meta data about species and pathway abundances in helthy and cancer subjects based on [curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/)
- [Installation script](src/install.R) for R packages including [BAnOCC](https://github.com/biobakery/banocc), [SpiecEasi](https://github.com/zdk123/SpiecEasi) and [drake](https://github.com/ropensci/drake)
- [Dockerfile](Dockerfile) to create a container with all tools needed for the analysis. This is based on [rocker/verse](https://hub.docker.com/r/rocker/verse) containing RStudio Server.
- [Makefile](Makefile) to build and run the docker container
- [Run script](src/analyze.R) to run the drake workflow and build the HTML report
- [Setup script](src/setup.R) to create the R environment with loaded packages and default options
- [Drake plans](src/plans) defining an example analysis workflow
- [Library](src/lib) with wrapper functions to call various correlation methods with a uniform interface
- [BibTeX](src/literature.bib) containing literature used

## Get started

Clone this repository with these shell commands.
GNU Make can be used to build and run the docker container:

```{bash}
git clone https://github.com/bioinformatics-leibniz-hki/correlation-workshop
cd correlation-workshop
make
```

This will deploy a webserver which can be accessed in your browser.
Go to the URL displayed and log in to RStudio Server using username `rstudio` and password `password`.

## Websites

- [Tutorial](https://bioinformatics-leibniz-hki.github.io/correlation-workshop)
- [Example report](https://bioinformatics-leibniz-hki.github.io/correlation-workshop/src/report.html)
