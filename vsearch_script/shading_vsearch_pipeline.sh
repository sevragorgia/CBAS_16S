#!/bin/sh                                                                       

# This is an example of a pipeline using vsearch to process data in the         
# Mothur 16S rRNA MiSeq SOP tutorial dataset to perform initial paired-end      
# read merging, quality filtering, chimera removal and OTU clustering.          

#taken from 
#note this exports my path and should commented or modified if required
export PATH=$PATH:/home/sergio/Software/vsearch-2.5.0/bin

THREADS=8
REF=../gold.fasta
PERL=$(which perl)
VSEARCH=$(which vsearch)

echo
echo Checking FASTQ format version for one file

$VSEARCH --threads $THREADS --fastq_chars $(ls -1 *.fastq.gz | head -1)

# Process samples                                                               

for f in *_R1_*.fastq; do

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(cut -d_ -f1 <<< "$f")

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================

    echo $f
    echo $r
    
    $VSEARCH --threads $THREADS \
        --fastq_mergepairs $f \
        --reverse $r \
        --fastq_minovlen 45 \
        --fastq_maxdiffs 0 \
        --fastqout $s.merged.fastq \
        --fastq_eeout 

    # Commands to demultiplex and remove tags and primers                       
    # using e.g. cutadapt may be added here.                                    

    echo
    echo Calculate quality statistics

    $VSEARCH --threads $THREADS \
        --fastq_eestats $s.merged.fastq \
        --output $s.stats
    echo
    echo Quality filtering

    $VSEARCH --threads $THREADS \
        --fastq_filter $s.merged.fastq \
        --fastq_maxee 0.5 \
        --fastq_minlen 225 \
        --fastq_maxlen 275 \
        --fastq_maxns 0 \
        --fastaout $s.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --threads $THREADS \
        --derep_fulllength $s.filtered.fasta \
        --strand plus \
        --output $s.derep.fasta \
        --sizeout \
        --uc $s.derep.uc \
        --relabel $s. \
        --fasta_width 0

done

echo Sum of unique sequences in each sample: $(cat *.derep.fasta | grep -c "^>")

# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               

echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

rm -f all.derep.fasta all.nonchimeras.derep.fasta
cat *.derep.fasta > all.fasta

echo
echo Dereplicate across samples and remove singletons

$VSEARCH --threads $THREADS \
    --derep_fulllength all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.derep.uc \
    --output all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" all.derep.fasta)

echo
echo Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
    --cluster_size all.derep.fasta \
    --id 0.98 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.preclustered.uc \
    --centroids all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" all.preclustered.fasta)

echo
echo De novo chimera detection

$VSEARCH --threads $THREADS \
    --uchime_denovo all.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" all.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection

$VSEARCH --threads $THREADS \
    --uchime_ref all.denovo.nonchimeras.fasta \
    --db $REF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" all.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

$PERL ../map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

$PERL ../map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids all.otus.fasta \
    --otutabout all.otutab.txt

echo
echo Number of OTUs: $(grep -c "^>" all.otus.fasta)

echo
echo Done

date
