
rule igetPacbio:
    output:
        expand('output/reads/dna/pacbio/{smrtCell}.subreads.bam', smrtCell = ['a', 'b', 'c', 'd', 'e']),
    shell:
        'bash scripts/igetPacbio.sh'
# local
# snakemake --printshellcmds --keep-going --reason output/reads/dna/pacbio/{a,b,c,d,e}.subreads.bam -n

rule igetIllumina:
    output:
        expand('output/reads/dna/illumina/{smrtCell}.subreads.bam', smrtCell = ['a', 'b', 'c', 'd', 'e']),
    shell:
        'bash scripts/igetIllumina.sh'
# local
# snakemake --printshellcmds --keep-going --reason output/reads/dna/illumina/{a,b,c,d,e}.subreads.bam -n

# Convert the BAM files containing the raw PacBio data into FASTQ files
rule bam2fastq:
    input:
        bam = 'output/reads/dna/pacbio/{smrtCell}.subreads.bam',
    output:
        fastq = 'output/reads/dna/pacbio/{smrtCell}.subreads.fq',
    wildcard_constraints:
        smrtCell = '[a-e]',
    benchmark:
        'output/benchmark/bam2fastq/{smrtCell}.tsv'
    log:
        'output/reads/dna/pacbio/{smrtCell}_bam2fastq.log'
    conda:
        'envs/conda/bamtools2_5_1.yaml'
    shell:
        'bamtools convert -format fastq -in {input.bam} -out {output.fastq} > {log} 2>&1'
# snakemake --use-conda --printshellcmds --reason output/reads/dna/pacbio/{a,b,c,d,e}.subreads.fq -n
# LSF
# snakemake --jobs 2000 --use-conda --printshellcmds --keep-going --reason output/reads/dna/pacbio/{a,b,c,d,e}.subreads.fq --cluster-config lsf.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -n

# Convert the PacBio FASTQ files to FASTA format
# PacBio does not provide base-level quality values so working
# with FASTQ files is unnecessary and takes more disk space
rule fastq2fastaPacbio:
    input:
        fastq = 'output/reads/dna/pacbio/{smrtCell}.subreads.fq',
    output:
        fasta = 'output/reads/dna/pacbio/{smrtCell}.subreads.fa',
    wildcard_constraints:
        smrtCell = '[a-e]'
    benchmark:
        'output/benchmark/fastq2fastaPacbio/{smrtCell}.tsv'
    log:
        'output/reads/dna/pacbio/{smrtCell}_fastq2fastaPacbio.stderr'
    conda:
        'envs/conda/seqtk1_3.yaml'
    shell:
        'seqtk seq -a {input.fastq} > {output.fasta} 2> {log}'
# LSF
# snakemake --jobs 2000 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason output/reads/dna/pacbio/{a,b,c,d,e}.subreads.fa --cluster-config lsf.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -n

# Merge PacBio FASTA files together
rule mergePacbioFastas:
    input:
        expand('output/reads/dna/pacbio/{smrtCell}.subreads.fa', smrtCell = ['a', 'b', 'c', 'd', 'e']),
    output:
        'output/reads/dna/pacbio/all.subreads.fa',
    benchmark:
        'output/benchmark/mergePacbioFastas/benchmark.tsv'
    log:
        'output/reads/dna/pacbio/mergePacbioFastas.stderr'
    shell:
        'cat {input} > {output} 2> {log}'
# local
# snakemake --cores 7 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason output/reads/dna/pacbio/all.subreads.fa -n

# Measure characteristics of raw PacBio reads
rule seqtkComp:
    input:
        reads = 'output/reads/dna/pacbio/all.subreads.fa',
    output:
        seqtkComp = 'output/reads/dna/pacbio/all.subreads.seqtkComp.txt'
    benchmark:
        'output/benchmark/seqtkComp/benchmark.tsv'
    log:
        'output/reads/dna/pacbio/seqtkComp.stderr'
    conda:
        'envs/conda/seqtk1_3.yaml'
    shell:
        'seqtk comp {input.reads} > {output.seqtkComp} 2> {log}'
# LSF
# snakemake --jobs 2000 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason --rerun-incomplete output/reads/dna/pacbio/all.subreads.seqtkComp.txt --cluster-config lsf.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -n

# https://gatkforums.broadinstitute.org/gatk/discussion/6484/how-to-generate-an-unmapped-bam-from-fastq-or-aligned-bam
rule fastq2bam:
    input:
        fastq1 = '../other/input/data/data/projects/caenorhabditis/raw/11400_Blaxter_Mark/Meloidogyne_Slovenia_samples/V13/{fileName}_1.fastq.gz',
        fastq2 = '../other/input/data/data/projects/caenorhabditis/raw/11400_Blaxter_Mark/Meloidogyne_Slovenia_samples/V13/{fileName}_2.fastq.gz',
    output:
        bam = 'output/reads/dna/illumina/V13/{fileName}.bam',
    benchmark:
        'output/benchmark/fastq2bam/{fileName}.tsv'
    log:
        'output/reads/dna/illumina/V13/{fileName}_fastq2bam.log'
    conda:
        'envs/conda/picard2_21_3.yaml'
    shell:
        'picard FastqToSam FASTQ={input.fastq1} FASTQ2={input.fastq2} OUTPUT={output.bam} READ_GROUP_NAME=bogus SAMPLE_NAME=V13 LIBRARY_NAME=TruSeqDNA_PCRfree PLATFORM_UNIT=bogus PLATFORM=illuminaHiSeqX SEQUENCING_CENTER=EDINBURGHGENOMICS RUN_DATE=2018-11-22T00:00:00-0400'
# local: limit CPU to two cores to not use up the laptop's memory
# snakemake --cores 2 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason --rerun-incomplete output/reads/dna/illumina/V13/180907_E00306_0370_AHMNKHCCXY_4_11400BM0014L01.bam output/reads/dna/illumina/V13/180907_E00306_0370_AHMNKHCCXY_5_11400BM0014L01.bam output/reads/dna/illumina/V13/180907_E00306_0370_AHMNKHCCXY_6_11400BM0014L01.bam output/reads/dna/illumina/V13/181005_E00306_0379_AHMW7NCCXY_7_11400BM0014L01.bam output/reads/dna/illumina/V13/181019_E00306_0386_AHT5G3CCXY_7_11400BM0014L01.bam output/reads/dna/illumina/V13/181005_E00306_0379_AHMW7NCCXY_8_11400BM0014L01.bam -n
# LSF
# snakemake --jobs 2000 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason --rerun-incomplete output/reads/dna/illumina/V13/180907_E00306_0370_AHMNKHCCXY_4_11400BM0014L01.bam output/reads/dna/illumina/V13/180907_E00306_0370_AHMNKHCCXY_5_11400BM0014L01.bam output/reads/dna/illumina/V13/180907_E00306_0370_AHMNKHCCXY_6_11400BM0014L01.bam output/reads/dna/illumina/V13/181005_E00306_0379_AHMW7NCCXY_7_11400BM0014L01.bam output/reads/dna/illumina/V13/181019_E00306_0386_AHT5G3CCXY_7_11400BM0014L01.bam output/reads/dna/illumina/V13/181005_E00306_0379_AHMW7NCCXY_8_11400BM0014L01.bam --cluster-config lsf.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -n

# cp ../other/input/data/V13.polished.fasta ../other/input/data/V13.polished2.fasta
# gzip ../other/input/data/V13.polished.fasta
# mv ../other/input/data/V13.polished2.fasta ../other/input/data/V13.polished.fasta
