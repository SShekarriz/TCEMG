#!/bin/bash
####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i g_and_g_TCEMG_${DATE}.log)
exec 2>&1
####################################################################################
#activate conda env:
source /home/shekas3/anaconda3/bin/activate anvio-6.2
anvi-gen-contigs-database

inter=$PWD/trim_interleaved
inter_pt=$PWD/trim_inter_PerPT
assembly=$PWD/trim_inter_PerPT_assembled
Samtools="/dataone/common/software/samtools-1.11/samtools"
binning=$PWD/trim_inter_PerPT_assem_binned
path1="/dataone/common/software/metabat2-V2021"


LongExractor=$PWD/bin/extract_LongContig_fromFasta.py
cover_dir=$PWD/trim_inter_PerPT_assem_coverage
ProkkaOut=$PWD/trim_inter_PerPT_assembled_prokka
bin_info_sorter=$PWD/bin/bin_info_sorter.R
bin_parser=$PWD/bin/bin_parser.R


AMJ154=(UCFMT1540 UCFMT1541 UCFMT1614)
BHI147=(UCFMT1608 UCFMT1609 UCFMT1707)
HMT141=(UCFMT1112 UCFMT1113 UCFMT1380)
JFM135=(UCFMT1044 UCFMT1045 UCFMT1148)
JLS151=(UCFMT1711 UCFMT1712 UCFMT1810)
JM150=(UCFMT1665 UCFMT1666 UCFMT1713)
JWC134=(UCFMT1042 UCFMT1043 UCFMT1116)
LRW144=(UCFMT1494 UCFMT1495 UCFMT1585)
MM148=(UCFMT1610 UCFMT1611 UCFMT1708)
TDW140=(UCFMT1114 UCFMT1115 UCFMT1379)

PATIENTS=(AMJ154 BHI147 HMT141 JFM135 JLS151 JM150 JWC134 LRW144 MM148 TDW140)
#PATIENTS=(HMT141)
#PATIENTS=(BHI147)
for PT in "${PATIENTS[@]}"; do
	FMTnum="${PT}[@]"
	echo "Patient: $PT"
        contig_file=$assembly/${PT}_inter_Cat3Sample/contigs.fasta
	echo "Ref assembly: $contig_file"

	echo "################################"
	echo "####### indexing contigs #######"
	echo "################################"

	if [ -f ${contig_file}.fai ]; then
   	echo "The file '${contig_file}.fai' exists."
	else
   	echo "The file '${contig_file}.fai' is not found."
   	echo "Indexing the contig file now"
   	bwa index $contig_file
   	samtools faidx $contig_file
	fi

	echo "################################"
        echo "## Running Prokka on contigs ###"
        echo "################################"

	ptProk=$ProkkaOut/${PT}_inter_Cat3Sample/${PT}.ffn
        if [ -f $ptProk ]; then
        echo "The file '$ptProk' exists."
        else
        echo "The file '$ptProk' is not found."
        echo "Running Prokka on the contig file now"
        prokka --cpus 15 --outdir $ProkkaOut/${PT}_inter_Cat3Sample --prefix $PT $contig_file
	fi

        echo "################################"
        echo "## Indexing Prokka ffn file  ###"
        echo "################################"
	
        if [ -f $ptProk.fai ]; then
        echo "The file '$ptProk.fai' exists."
        else
        echo "The file '$ptProk.fai' is not found."
        echo "Indexing  Prokka ffn file now"
        bwa index $ptProk
        samtools faidx $ptProk
        fi

	echo "################################"
        echo "#####selecting large contigs####"
        echo "################################"
	
	ptRef=${contig_file%.fasta}_1k.fa
	if [ -f $ptRef ]; then
        echo "The file '$ptRef' exists."
        else
        echo "The file '$ptRef' is not found."
        echo "Selecting contig >= 1kb"
        python2.7 $LongExractor $contig_file
	mv ${contig_file}_out $ptRef
	fi

	echo "################################"
        echo "#indexing contigs Large contigs#"
        echo "################################"

        if [ -f ${ptRef}.fai ]; then
        echo "The file '${ptRef}.fai' exists."
        else
        echo "The file '${ptRef}.fai' is not found."
        echo "Indexing the contig file now"
        bwa index $ptRef
        samtools faidx $ptRef
        fi
	
	echo "################################"
        echo "#####Building anvio contigdb####"
        echo "################################"
	
	ptAnvi=${ptRef%.fa}.db
        if [ -f $ptAnvi ]; then
        echo "The file '$ptAnvi' exists."
        else
        echo "The file '$ptAnvi' is not found."
        echo "Building anvio contig database"
	anvi-gen-contigs-database -f $ptRef -o $ptAnvi
	echo "Building HMMs models"
	anvi-run-hmms -c $ptAnvi      
	fi

	echo "################################"
        echo "####### Short-read mapping #####"
        echo "################################"

	mkdir -p $cover_dir
	cover=${cover_dir}/${PT}_info
	prokcover=${cover_dir}/${PT}_prok_ffn
	mkdir -p $cover $prokcover
	for sample in "${!FMTnum}"; do
		echo "Mapping $sample to $PT as follow:"
		sample_file=$inter/${sample}_inter.fastq
		echo "Sample:$sample_file"
		echo "pt:$ptRef"
	
		if [ -f $cover/${sample}.bam ]; then
        	echo "The file '$cover/${sample}.bam' exists."
        	else
        	echo "The file '$cover/${sample}.bam' is not found."
        	echo "Mapping short-reads to assembly"
        	bwa mem -t 15 -p $ptRef $sample_file | $Samtools sort -@ 15 | $Samtools view -@ 15 -bS -o $cover/${sample}.bam
		fi
		
		if [ -f $cover/${sample}.bam.bai ]; then
                echo "The file '$cover/${sample}.bam.bai' exists."
                else
                echo "The file '$cover/${sample}.bam.bai' is not found."
                $Samtools index $cover/${sample}.bam
		fi
		
		if [ -f $cover/${sample}.coverage.txt ]; then
                echo "The file '$cover/${sample}.coverage.txt' exists."
                else
                echo "The file '$cover/${sample}.coverage.txt' is not found."
                echo "Calculating the number of mapped reads"
		$Samtools coverage $cover/${sample}.bam > $cover/${sample}.coverage.txt	
                fi
		
		if [ -f $cover/${sample}.bam-ANVIO_PROFILE/PROFILE.db ]; then
                echo "The file '$cover/${sample}.bam-ANVIO_PROFILE/PROFILE.db' exists."
                else
                echo "The file '$cover/${sample}.bam-ANVIO_PROFILE' is not found."
                echo "Building anvio sample database"
                anvi-profile -i $cover/${sample}.bam -c $ptAnvi
		fi

		if [ -f $prokcover/${sample}.bam ]; then
                echo "The file '$prokcover/${sample}.bam' exists."
                else
                echo "The file '$prokcover/${sample}.bam' is not found."
                echo "Mapping short-reads to Prokka ffn file for gene coverage"
                bwa mem -t 15 -p $ptProk $sample_file | $Samtools sort -@ 15 | $Samtools view -@ 15 -F 4 -o $prokcover/${sample}.bam
                fi

		if [ -f $prokcover/${sample}.bam.bai ]; then
                echo "The file '$prokcover/${sample}.bam.bai' exists."
                else
                echo "The file '$prokcover/${sample}.bam.bai' is not found."
                $Samtools index $prokcover/${sample}.bam
                fi

                if [ -f $prokcover/${sample}_depth.txt ]; then
                echo "The file '$prokcover/${sample}_depth.txt' exists."
                else
                echo "The file '$prokcover/${sample}_depth.txt' is not found."
                echo "Calculating per base coverage for prokka ffn genes"
                $Samtools depth -aa $prokcover/${sample}.bam > $prokcover/${sample}_depth.txt
		fi

	done

	echo "########################################"
	echo "Binning the cotnigs using Metabat"
	echo "########################################"
        
	BIN=${binning}/${PT}_bins
        mkdir -p $BIN
	
	metabat_tbl=$binning/${PT}_depthForMetabat.txt
	if [ -f $metabat_tbl ]; then
        echo "The file '$metabat_tbl' exists."
        else
        echo "The file '$metabat_tbl' is not found."
        echo "Making a coverage table from bam file for metabat2"
        $path1/jgi_summarize_bam_contig_depths --outputDepth $metabat_tbl $cover/*.bam
        fi

	if [ -f $BIN/bin.1.fa ]; then
        echo "The file '$BIN/bin.1.fa' exists."
        else
        echo "The file '$BIN/bin.1.fa' is not found."
        echo "Binning the contigs"
        $path1/metabat2 -i $ptRef -a $metabat_tbl -o $BIN/bin
	fi
	
	binInfo=$BIN/contig_in_bins_anvi.txt
	if [ -f $binInfo ]; then
	echo "The file '$binInfo' exists."
	else
	echo "The file '$binInfo' is not found."
	grep -e "^>" $BIN/bin* > $BIN/contig_in_bins.txt
	Rscript $bin_info_sorter $BIN/contig_in_bins.txt $binInfo
	rm $BIN/contig_in_bins.txt	
	fi

	echo "################################"
        echo "#####Building anvio sampledb####"
        echo "################################"

        ptAnviProfile=$cover/SAMPLES-MERGED/PROFILE.db
        if [ -f $ptAnviProfile ]; then
        echo "The file '$ptAnviProfile' exists."
        else
        echo "The file '$ptAnviProfile' is not found."
	anvi-merge $cover/*ANVIO_PROFILE/PROFILE.db -o $cover/SAMPLES-MERGED -c $ptAnvi
        echo "Importing metabat2 bin info"
	anvi-import-collection $binInfo -p $ptAnviProfile -c $ptAnvi --contigs-mode -C metabat2
	anvi-script-add-default-collection -p $ptAnviProfile -c $ptAnvi -C default
	fi	
        
	echo "################################"
	echo "#####Exporting coverage table###"
        echo "################################"		
	
	ptAnviSUM_M=$cover/SAMPLES-MERGED-SUM-metabat2/bins_summary.txt
	if [ -f $ptAnviSUM_M ]; then
        echo "The file '$ptAnviSUM_M' exists."
        else
        echo "The file '$ptAnviSUM_M' is not found."
        anvi-summarize -p $ptAnviProfile -c $ptAnvi -C metabat2 -o $cover/SAMPLES-MERGED-SUM-metabat2 --init-gene-coverages
        fi

	ptAnviSUM_D=$cover/SAMPLES-MERGED-SUM-default/bins_summary.txt
        if [ -f $ptAnviSUM_D ]; then
        echo "The file '$ptAnviSUM_D' exists."
        else
        echo "The file '$ptAnviSUM_D' is not found."
        anvi-summarize -p $ptAnviProfile -c $ptAnvi -C default -o $cover/SAMPLES-MERGED-SUM-default --init-gene-coverages
        fi


	echo "#########################################"
	echo "check the quality of BINS using checkM"
	echo "##########################################"

        BINQ=${binning}/${PT}_bins_quality
        mkdir -p $BINQ
	echo "CHANGING THE CONDA ENVIRONMENT!!!!"
	source /home/shekas3/anaconda3/bin/activate py2
	
	if [ -f $BINQ/lineage_wf/lineage.ms ]; then
	echo "The file '$BINQ/lineage_wf/lineage.ms' exists."
        else
	echo "The file '$BINQ/lineage_wf/lineage.ms' is not found."
	checkm lineage_wf -x fa -t 10 $BIN $BINQ/lineage_wf
	fi

	if [ -f $BINQ/lineage_wf/${PT}_bin_qc.txt ]; then
        echo "The file '$BINQ/lineage_wf/${PT}_bin_qc.txt' exists."
        else
        echo "The file '$BINQ/lineage_wf/${PT}_bin_qc.txt' is not found."
        checkm qa $BINQ/lineage_wf/lineage.ms \
        $BINQ/lineage_wf -o 1 --tab_table > $BINQ/lineage_wf/${PT}_bin_qc.txt
        fi

	echo "######################################"
	echo "Running gtdbtk (taxonomy) for each bin"
	echo "######################################"

	BINT=${binning}/${PT}_bins_taxonomy
	mkdir -p $BINT
	echo "CHANGING THE CONDA ENVIRONMENT!!!!"
        source /home/shekas3/anaconda3/bin/activate gtdbtk
	
	gtdb_out="$BINT/classify/${PT}.bac120.summary.tsv"
	if [ -f $gtdb_out ]; then
	echo "The file '$gtdb_out' exists."
        else
        echo "The file '$gtdb_out' is not found."
	gtdbtk classify_wf --genome_dir $BIN --cpus 20 --out_dir $BINT --extension fa --prefix $PT
	fi

 	echo "######################################"
        echo "Parsing all bin info into one"
        echo "######################################"
	
	BINSum=${binning}/${PT}_bins_summary
        mkdir -p $BINSum
	
        checkm=$BINQ/lineage_wf/${PT}_bin_qc.txt
        pathbincov=$cover/SAMPLES-MERGED-SUM-metabat2/bins_across_samples
        mkdir -p $pathbincov/other_files


        if [ -f $pathbincov/other_files/bins_percent_recruitment.txt ]; then
        echo "The file 'bins_percent_recruitment.txt' exists."
        else
        echo "The file 'bins_percent_recruitment.txt' is not found."
        mv $pathbincov/bins_percent_recruitment.txt $pathbincov/other_files/
        mv $pathbincov/hmm*txt $pathbincov/other_files/
	fi

	summary="$BINSum/${PT}_bin_sum_cov.txt"
	if [ -f $summary ]; then
        echo "The file '$summary' exists."
        else
        echo "The file '$summary' is not found."
        Rscript $bin_parser $checkm $gtdb_out $pathbincov $summary
	fi


	echo "######################################"
        echo "Finding plasmids from the assemblies #"
        echo "######################################"

	echo "CHANGING THE CONDA ENVIRONMENT!!!!"
	source /home/shekas3/anaconda3/bin/activate SCAPP
	
	assm_dir=$assembly/${PT}_inter_Cat3Sample
	graph=$assm_dir/assembly_graph.fastg
	cor=corrected/${PT}_inter_Cat3Sample
	cor_R1=$assm_dir/${cor}_1.00.0_0.cor.fastq.gz
	cor_R2=$assm_dir/${cor}_2.00.0_0.cor.fastq.gz

	out=${assm_dir}_plasmid/assembly_graph.confident_cycs.fasta
	if [ -f $out ]; then
        echo "The file '$out' exists."
        else
        echo "The file '$out' is not found."
	scapp -g $graph -o ${assm_dir}_plasmid -k 55 -r1 $cor_R1 -r2 $cor_R2 -p 25
	fi

done

echo "################### DONE ####################################"
date -u
echo "#############################################################"
