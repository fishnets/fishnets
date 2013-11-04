# A shell script to concatenate all code files together
sed -s '$a \\n###################################################################################\n' \
	object.R \
	node.R \
	species-any.R \
	taxon-lookupper.R \
	taxonomic-imputer.R \
	glmer.R \
	svmer.R \
	distributions.R \
	fishnet.R \
	tests.R \
	> all.R
