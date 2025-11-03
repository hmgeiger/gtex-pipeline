# TOPMed @ NYGC RNA-seq pipeline, HiGlobin modification, harmonization summary

### *Note: this is a draft of the next pipeline iteration. Software and reference versions may change.*

This document is an extension of the documentation provided [here](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) by the Broad Institute.

However, one important note is that the NYGC version of the pipeline was frozen at an earlier time point, and so the version of STAR is different (v2.6.1d for NYGC, versus v2.7.10a for the Broad).

NYGC is not currently hosting the compiled reference files for this pipeline, in cases where they differ from those of the Broad. But these may be generated from publicly available files, using the instructions below.

The below instructions will focus on the steps needed to create the input for the STAR index creation (genomeGenerate) step. This does not include collapsing the annotation, which only applies to the RNA-SeQC step.


#### Genome reference

For RNA-seq analyses, a reference FASTA excluding ALT, HLA, and Decoy contigs was generated (as of the writing of this document, none of the RNA-seq pipeline components were designed to properly handle these regions).

*Note: The reference produced after filtering out ALT, HLA, and Decoy contigs is identical to the one used by ENCODE ([FASTA file](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)). However, the ENCODE reference does not  include ERCC spike-in sequences.*

1. The Broad Institute's [GRCh38 reference](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) can be obtained from the Broad Institute's FTP server (ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/) or from [Google Cloud](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0).
2. ERCC spike-in reference annotations were downloaded from [ThermoFisher](https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip) (the archive contains two files, ERCC92.fa and ERCC92.gtf) and processed as detailed below. The patched ERCC references are also available [here](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq/references).<br>The ERCC FASTA was patched for compatibility with RNA-SeQC/GATK using
    ```bash
    sed 's/ERCC-/ERCC_/g' ERCC92.fa > ERCC92.patched.fa
    ```
3. ALT, HLA, and Decoy contigs were excluded from the reference genome FASTA using the following Python code:
    ```python
    with open('Homo_sapiens_assembly38.fasta', 'r') as fasta:
        contigs = fasta.read()
    contigs = contigs.split('>')
    contig_ids = [i.split(' ', 1)[0] for i in contigs]

    # exclude ALT, HLA and decoy contigs
    filtered_fasta = '>'.join([c for i,c in zip(contig_ids, contigs)
        if not (i[-4:]=='_alt' or i[:3]=='HLA' or i[-6:]=='_decoy')])

    with open('Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta', 'w') as fasta:
        fasta.write(filtered_fasta)
    ```

4. ERCC spike-in sequences were appended:
    ```bash
    cat Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta ERCC92.patched.fa \
        > Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
    ```
5. The FASTA index and dictionary (required for Picard/GATK) were generated:
    ```bash
    samtools faidx Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
    java -jar picard.jar \
        CreateSequenceDictionary \
        R=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
        O=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict
    ```

#### Reference annotation
The reference annotations were prepared as follows:
1. The comprehensive gene annotation was downloaded from [GENCODE](https://www.gencodegenes.org/human/release_39.html): [gencode.v39.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz).

2. Gene- and transcript-level attributes were added to the ERCC GTF with the following Python code:
    ```python
    with open('ERCC92.gtf') as exon_gtf, open('ERCC92.patched.gtf', 'w') as gene_gtf:
        for line in exon_gtf:
            f = line.strip().split('\t')
            f[0] = f[0].replace('-', '_')  # required for RNA-SeQC/GATK (no '-' in contig name)
            f[5] = '.'

            attr = f[8]
            if attr[-1] == ';':
                attr = attr[:-1]
            attr = dict([i.split(' ') for i in attr.replace('"', '').split('; ')])
            # add gene_name, gene_type
            attr['gene_name'] = attr['gene_id']
            attr['gene_type'] = 'ercc_control'
            attr['gene_status'] = 'KNOWN'
            attr['level'] = 2
            for k in ['id', 'type', 'name', 'status']:
                attr[f'transcript_{k}'] = attr[f'gene_{k}']

            attr_str = []
            for k in ['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
                'transcript_type', 'transcript_status', 'transcript_name']:
                attr_str.append(f'{k} "{attr[k]}";')
            attr_str.append(f"level {attr['level']};")
            f[8] = ' '.join(attr_str)

            # write gene, transcript, exon
            gene_gtf.write('\t'.join(f[:2] + ['gene'] + f[3:]) + '\n')
            gene_gtf.write('\t'.join(f[:2] + ['transcript'] + f[3:]) + '\n')
            f[8] = ' '.join(attr_str[:2])
            gene_gtf.write('\t'.join(f[:2] + ['exon'] + f[3:]) + '\n')
    ```

3. The ERCC annotation was appended to the reference GTF:
    ```bash
    cat gencode.v39.GRCh38.annotation.gtf ERCC92.patched.gtf \
        > gencode.v39.GRCh38.annotation.ERCC.gtf
    ```

#### STAR index

The STAR index was originally built to accommodate 2x101bp paired-end reads, and so the STAR index overhang was set accordingly (read length - 1 = 100).

The code below uses 10 threads, but adjust accordingly based on the tradeoff between computational resource availability and runtime.

```bash
NB_CPU=10
STAR --runMode genomeGenerate \
    --genomeDir STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v39_oh100 \
    --genomeFastaFiles Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta ERCC92.patched.fa \
    --sjdbGTFfile gencode.v39.GRCh38.annotation.ERCC.gtf \
    --sjdbOverhang 100 --runThreadN $NB_CPU
```

#### Installation of custom STAR

The HiGlobin modification to this pipeline requires using a custom fork of STAR, that has been modified to be able to run 2-pass mode as two separate steps, with identical results to running as one single command (--twoPassMode).

See here (insert link here to the relevant fork link on Github once I get it from Manisha).

### Pipeline parameters

This section contains the full list of inputs and parameters provided to STAR.

The following variables must be defined:
* `star_index`: path to the directory containing the STAR index (STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v39_oh100 if following the code above)
* `fastq1` and `fastq2`: paths to the two FASTQ files (pipeline assumes they are zipped and can be unzipped with zcat)
* `sample_id`: sample identifier; this will be prepended to output files
* `star_bam_file`: name/file path of the BAM generated by the STAR aligner, by default `${sample_id}.Aligned.sortedByCoord.out.bam`.

The code below uses 8 threads, but adjust accordingly based on the tradeoff between computational resource availability and runtime.

The argument --twopassMode Basic, when running with the custom STAR, does not actually run both passes.

Rather, it runs the equivalent of only the first pass using the same parameters that would normally apply when running in two pass mode.

Then, we need to separately run the second pass.

```bash
NB_CPU=8
STAR --runMode alignReads \
 --runThreadN $NB_CPU \
 --genomeDir ${star_index} \
 --twopassMode Basic \
 --outFilterIntronMotifs None \
 --alignSJDBoverhangMin 1 \
 --alignMatesGapMax 1000000 \
 --outFilterMismatchNoverLmax 0.1 \
 --outFilterMismatchNmax 999 \
 --outFilterMultimapNmax 20 \
 --outSAMtype None \
 --runMode alignReads \
 --readFilesCommand zcat \
 --alignSJoverhangMin 8 \
 --alignIntronMin 20 \
 --outFilterScoreMinOverLread 0.33 \
 --outFilterMatchNminOverLread 0.33 \
 --alignIntronMax 1000000 \
 --limitSjdbInsertNsj 1200000 \
 --genomeLoad NoSharedMemory \
 --alignSoftClipAtReferenceEnds Yes \
 --readFilesIn ${fastq1} ${fastq2} \
 --outFileNamePrefix ${sample_id}-firstpass.

SJ=${sample_id}-firstpass.SJ.out.tab
filtered_sj=${sample_id}-firstpass-filtered.SJ.out.tab

#Remove unannotated junctions (field 6=0) if they start or end within the range chr11:5225464-5229395.
grep -w -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10' $sj > $filtered_sj
grep -w chr11 $sj | awk '$6 == 1 || !(($2 >= 5225464 && $2 <= 5229395) || ($3 >= 5225464 && $3 <= 5229395))' >> $filtered_sj
grep -v -E -w 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11' $sj >> $filtered_sj

STAR \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFilterIntronMotifs None \
        --alignSJDBoverhangMin 1 \
        --alignMatesGapMax 1000000 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterMismatchNmax 999 \
        --outFilterMultimapNmax 20 \
        --outSAMunmapped Within \
        --chimMainSegmentMultNmax 1 \
        --outSAMtype BAM Unsorted \
        --outFilterType BySJout \
        --runMode alignReads \
        --readFilesCommand zcat \
        --chimOutType Junctions WithinBAM SoftClip \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS nM NM ch \
        --alignSJoverhangMin 8 \
        --outFileNamePrefix ${sample_id}. \
        --chimJunctionOverhangMin 15 \
        --chimSegmentMin 15 \
        --alignIntronMin 20 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --alignIntronMax 1000000 \
        --limitSjdbInsertNsj 1200000 \
        --genomeLoad NoSharedMemory \
        --alignSoftClipAtReferenceEnds Yes \
        --runThreadN $NB_CPU \
        --genomeDir ${star_index} \
        --readFilesIn ${fastq1} ${fastq2} \
        --sjdbFileChrStartEnd $filtered_sj 
```
