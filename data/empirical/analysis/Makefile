.PHONY: mpr ancestry flux geoarg

mpr:
	Rscript --vanilla mpr.R
	Rscript --vanilla mpr-subsets.R

ancestry:
	Rscript --vanilla ancestry.R 100
	Rscript --vanilla ancestry-subsets.R 10

flux:
	Rscript --vanilla flux.R 100
	Rscript --vanilla flux-subsets.R 10

geoarg:
	Rscript --vanilla geoarg.R