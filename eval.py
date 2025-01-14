#!/bin/env python
import polars as pl
def fitter(file_name):
    m5c_candidates = pl.read_csv(file_name + "/" + file_name + ".filtered.tsv",
                                separator="\t",
                                dtypes={
                                    "ref": pl.Utf8,
                                    "pos": pl.Int64,
                                    "strand": pl.Utf8
                                },)
    passed_candidates = m5c_candidates.filter(pl.col("passed") == True).filter(pl.col("pval") < 1e-6)
    passed_candidates.write_csv("passed_" + file_name + ".tsv", separator="\t", include_header=True)
    return passed_candidates

gse_data = pl.read_csv(
	"GSE225614_HeLa-WT_sites.tsv",
	separator="\t",
	dtypes={
		"chromosome": pl.Utf8,
		"position": pl.Int64,
		"strand": pl.Utf8,
		"gene_pos": pl.Int64,
	},
	null_values=["."],  # add '.' to the null_values list
)
gse_data = gse_data.rename({"chromosome": "ref", "position": "pos"})

for file_name in ["SRR23538290", "SRR23538291", "SRR23538292"]:
    passed_candidates = fitter(file_name)
    merged_data = passed_candidates.join(gse_data, on=["ref", "pos"], how="inner")
    print(merged_data.shape[0]/passed_candidates.shape[0])
    pearson_corr = merged_data.select(pl.corr("ratio", "ur"))
    print(pearson_corr)