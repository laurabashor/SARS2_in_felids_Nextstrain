# Data and code to support: Rapid evolution of SARS-CoV-2 in domestic cats
### Laura Bashor, Roderick B Gagne, Angela Bosco-Lauth, Mark Stenglein, and Sue VandeWoude

Data and R code for processing, analyzing & visualizing results after running https://github.com/stenglein-lab/viral_variant_caller pipeline.

Three cohorts of cats (N=23) were inoculated with SARS-CoV-2 USA-WA1/2020 or infected via cat-to-cat contact transmission. Full viral genomes were recovered from RNA obtained from nasal washes 1-3 days post-infection, sequenced and analyzed for within-host viral variants. 

Raw sequence data is publicly available in NCBI SRA database under Bioprojects PRJNA704947 and 
PRJNA842572.

The unmodified excel spreadsheet variant table output by the pipeline listing all the within-host variants detected within the cats and the P1, P2, and P3 stocks used for inoculations is available in the folder named "viral_variant_caller pipeline output" with the file name "variant_summary_0.001_cutoff_7.13.21.xlsx". Processed variant tables containing only variants detected in both technical sequencing replicates above allele frequency cutoffs of 0.1% and 3% are available in the folder named "processed variant tables" with the file names "variant_summary_0.001_processed.xlsx" and "variant_summary_0.03_processed.xlsx." These files include the following columns:

| Column | Description |
| --- | --- |
| `reference_sequence` | The GenBank accession number of the viral reference sequence used to call variants |
| `position` | Nucleotide position of the variant within the SARS2 genome |
| `gene` | The open reading frame or gene location of the variant  |
| `codon` | The codon location of the variant |
| `indel` | TRUE/FALSE if the gene is a structural variant (indel) |
| `variant` | The variant name with standard amino acid mutation nomenclature  |
| `reference_base` | The base in the reference sequence |
| `variant_base` | The new base in the sample |
| `effect` | The predicted effect of the variant |
| `featureid` | The presence of this variant in SARS2 Variant of Concern (VOC) lineages (not comprehensive) |
| `[Dataset ID]` | Each remaining column is named for the sample (cat or viral stock inoculum passage) and contains the allele frequency at which the variant was detected within that sample |

########################

Custom Nextstrain builds for SARS-CoV-2 sequences recovered from domestic cats and other felids (all sequences obtained from GISAID.org):
- Domestic cats: https://nextstrain.org/community/laurabashor/SARS2felids/catphylogeny
- All felids: https://nextstrain.org/community/laurabashor/SARS2felids/felidphylogeny

