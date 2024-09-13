.PHONY: data hgdp

data: hgdp_tgp_sgdp_high_cov_ancients_chr18_p.dated.trees
	Rscript --vanilla make-data.R

hgdp:
	curl -O https://zenodo.org/records/5512994/files/hgdp_tgp_sgdp_high_cov_ancients_chr18_p.dated.trees.tsz?download=1
	tsunzip hgdp_tgp_sgdp_high_cov_ancients_chr18_p.dated.trees.tsz
