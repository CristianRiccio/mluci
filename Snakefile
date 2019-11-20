
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
# snakemake --jobs 2000 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason output/reads/dna/pacbio/all.subreads.seqtkComp.txt --cluster-config lsf.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -n
