anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
	-o concatenated-proteins.fa --hmm-source Campbell_et_al \
	--gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 --return-best-hit --get-aa-sequences --align-with muscle --concatenate

anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
	-o concatenated-proteins.fa --hmm-source Campbell_et_al \
	--gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 --return-best-hit --align-with muscle --concatenate
