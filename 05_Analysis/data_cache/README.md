# data_cache

Input and intermediate files used by the pipeline. Only the small files marked **included** are in the repository; the others are public downloads or are rebuilt by the scripts. Sources, URLs, sizes and MD5 checksums are in `../data_manifest.csv`.

| File | Source | Size (MB) | Status |
|---|---|---|---|
| `ABSOLUTE_mastercalls.txt` | GDC PanCanAtlas supplement (ABSOLUTE purity/ploidy calls; Taylor 2018) | 0.90 | download |
| `cbio_attrs.json` | cBioPortal REST API, study coadread_tcga_pan_can_atlas_2018 | 0.02 | download |
| `cbio_mutations_selected.json` | cBioPortal REST API, study coadread_tcga_pan_can_atlas_2018 | 2.35 | download |
| `cbio_pancan_coadread.rds` | cBioPortal REST API, study coadread_tcga_pan_can_atlas_2018 | 0.06 | download |
| `cbio_pancan_patient.json` | cBioPortal REST API, study coadread_tcga_pan_can_atlas_2018 | 4.62 | download |
| `cbio_pancan_sample.json` | cBioPortal REST API, study coadread_tcga_pan_can_atlas_2018 | 3.47 | download |
| `Charoentong2017.gmt` | Charoentong 2017 Table S6 gene sets (GitHub mirror msfuji/immunoduct) | 0.01 | included |
| `COAD_clinicalMatrix.tsv` | UCSC Xena legacy TCGA hub (TCGA.COAD/READ.sampleMap) | 0.52 | download |
| `COAD_gistic_thresholded.tsv.gz` | UCSC Xena legacy TCGA hub (TCGA.COAD/READ.sampleMap) | 0.52 | download |
| `COAD_HM450.tsv.gz` | UCSC Xena legacy TCGA hub (TCGA.COAD/READ.sampleMap) | 308.84 | download |
| `cohort.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 39.84 | rebuilt by the scripts |
| `de_paired.rds` | GDC API (STAR counts of 32 matched normal/tumour pairs) | 1.31 | download |
| `de_tumor_normal.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 1.66 | rebuilt by the scripts |
| `gdc_normal_tumour_file_listing.csv` | GDC API (STAR counts of 32 matched normal/tumour pairs) | 0.02 | included |
| `gdc_paired_manifest.csv` | GDC API (STAR counts of 32 matched normal/tumour pairs) | 0.01 | included |
| `geo_cohorts.rds` | NCBI GEO (GSE39582, GSE17536; GPL570) | 130.94 | download |
| `gistic_thresholded_coadread.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 0.26 | rebuilt by the scripts |
| `GPL570_table.txt` | NCBI GEO (GSE39582, GSE17536; GPL570) | 79.50 | download |
| `GSE17536_series_matrix.txt.gz` | NCBI GEO (GSE39582, GSE17536; GPL570) | 52.24 | download |
| `GSE39582_series_matrix.txt.gz` | NCBI GEO (GSE39582, GSE17536; GPL570) | 172.16 | download |
| `hm450_probe_annotation_subset.csv` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 0.02 | included |
| `hm450_subset.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 0.45 | rebuilt by the scripts |
| `illuminaMethyl450_probeMap.tsv` | UCSC Xena HM450 probe map | 18.20 | download |
| `log2cpm_tumor_591.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 40.83 | rebuilt by the scripts |
| `MCPcounter_genes.txt` | MCP-counter marker signatures (Becht 2016) | 0.01 | included |
| `MCPcounter_probesets.txt` | MCP-counter marker signatures (Becht 2016) | 0.01 | included |
| `READ_clinicalMatrix.tsv` | UCSC Xena legacy TCGA hub (TCGA.COAD/READ.sampleMap) | 0.18 | download |
| `READ_gistic_thresholded.tsv.gz` | UCSC Xena legacy TCGA hub (TCGA.COAD/READ.sampleMap) | 0.23 | download |
| `READ_HM450.tsv.gz` | UCSC Xena legacy TCGA hub (TCGA.COAD/READ.sampleMap) | 99.25 | download |
| `Thorsson_2018_immune_landscape.xlsx` | PanCanAtlas immune landscape (Thorsson 2018) incl. TCGA-CDR (Liu 2018) via GDC/Xena | 2.95 | download |
| `Thorsson_download.zip` | PanCanAtlas immune landscape (Thorsson 2018) incl. TCGA-CDR (Liu 2018) via GDC/Xena | 2.95 | download |
| `Thorsson_immune_subtypes.tsv.gz` | PanCanAtlas immune landscape (Thorsson 2018) incl. TCGA-CDR (Liu 2018) via GDC/Xena | 2.95 | download |
| `Thorsson_Subtype_Immune_Model_Based.txt.gz` | PanCanAtlas immune landscape (Thorsson 2018) incl. TCGA-CDR (Liu 2018) via GDC/Xena | 0.05 | download |
| `vst_all_624.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 40.96 | rebuilt by the scripts |
| `vst_tumor_591.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 38.50 | rebuilt by the scripts |
| `xena_clinical_coadread.rds` | Derived by the scripts in 05_Analysis/scripts from the supplied TCGA count matrix (TCGA_CRC_combined.rds) and the sources above | 0.16 | rebuilt by the scripts |

The primary RNA-seq input, `TCGA_CRC_combined.rds` (GDC STAR - Counts gene-level counts and clinical data for TCGA-COAD and TCGA-READ, 625 samples, 210 MB), is read from the folder given by `GLP2R_ORIG`; it is not included because of its size.
