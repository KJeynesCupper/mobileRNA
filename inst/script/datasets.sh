## Here is listed the scripts used to generate the data files

#### reduced_chr12_Eggplant.gff ####
wget https://solgenomics.net/ftp/genomes/Solanum_melongena_V4.1/Eggplant_V4.1_function_IPR_final.gff
grep '^12\s' Eggplant_V4.1_function_IPR_final.gff | head -n 8 > reduced_chr12_Eggplant.gff

#### reduced_chr12_Eggplant.fa ####
wget https://solgenomics.net/ftp/genomes/Solanum_melongena_V4.1/Eggplant_V4.1.fa
grep -A 49 '>12' Eggplant_V4.1.fa > reduced_chr12_Eggplant.fa # select chromosome 12, and 50 lines

#### reduced_chr2_Tomato.gff ####
wget https://solgenomics.net/ftp/tomato_genome/Heinz1706/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff
grep -E 'ID=gene:Solyc02g063520\.4|ID=gene:Solyc02g078570\.3|ID=gene:Solyc02g089300\.3' ITAG4.1_gene_models.gff > reduced_chr2_Tomato.gff


#### reduced_chr2_Tomato.fa ####
wget https://solgenomics.net/ftp/tomato_genome/Heinz1706/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa
sed -i '' 's/\.//g' S_lycopersicum_chromosomes.4.00.fa && grep -A 49 'SL40ch02' S_lycopersicum_chromosomes.4.00.fa > reduced_chr2_Tomato.fa
