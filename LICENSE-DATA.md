# Licence for derived results and tables

The **analysis code** in this repository is licensed under the MIT License — see `LICENSE`.

The **derived results and tables** in this repository — i.e. everything under `results/` and `tables/` and the result values quoted in `README.md` and `RUNNING.md` — are licensed under the **Creative Commons Attribution 4.0 International licence (CC BY 4.0)**:

- Licence deed: https://creativecommons.org/licenses/by/4.0/
- Legal code: https://creativecommons.org/licenses/by/4.0/legalcode

You are free to share and adapt the derived results and tables for any purpose, provided that you give appropriate credit (cite the manuscript and this archive), link to the licence, and indicate if changes were made.

## Third-party data

This repository does **not** redistribute any third-party data. The analyses use data obtained from:

- **TCGA / GDC** (https://portal.gdc.cancer.gov) — TCGA data are subject to the NIH Genomic Data Sharing Policy and the TCGA data-use terms; the data used here are open-access (level 3) tumour RNA-seq, clinical, copy-number, methylation and mutation data.
- **UCSC Xena** (https://xenabrowser.net) — TCGA-COADREAD processed matrices derived from TCGA; use subject to the Xena terms and the underlying TCGA terms.
- **cBioPortal** (https://www.cbioportal.org) — PanCancer Atlas clinical and mutation data; use subject to the cBioPortal terms and the underlying TCGA publication terms.
- **NCBI GEO** (https://www.ncbi.nlm.nih.gov/geo) — series GSE39582 and GSE17536; use subject to the GEO terms and the original series authors' terms.
- **PanCanAtlas / Thorsson et al. 2018 and Liu et al. 2018** (https://gdc.cancer.gov/about-data/publications/panimmune) — immune-subtype and TCGA-CDR data; use subject to the PanCanAtlas publication terms.
- **MSigDB** via the R package `msigdbr` (https://www.gsea-msigdb.org) — gene-set collections; use subject to the MSigDB licence.
- **Charoentong et al. 2017** (https://github.com/msfuji/immunoduct) — immune cell-type gene sets; cite the original publication.
- **MCP-counter / Becht et al. 2016** (https://github.com/ebecht/MCPcounter) — marker signatures; cite the original publication.

The file `results/data_manifest.csv` records the source URL of every input used. Users of this repository are responsible for complying with the terms of the data providers above.

## Citation

If you use the derived results or tables, please cite the manuscript and this archive; the machine-readable citation is in `CITATION.cff`.
