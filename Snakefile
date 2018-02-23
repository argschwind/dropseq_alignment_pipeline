configfile: "config.yml"

# import modules and define variables ------------------------------------------

import glob

dropseq_root = "Drop-seq_tools-1.13"
picard_jar = dropseq_root + "/3rdParty/picard/picard.jar"

# define input function --------------------------------------------------------

# infer input fastq files from dirs in config file and sample wildcard
def get_fastq_files(wildcards):
    indir = config["samples"][wildcards.sample]
    f1 = glob.glob(indir + "/*" + wildcards.sample  + "_1_sequence.txt.gz")
    f2 = glob.glob(indir + "/*" + wildcards.sample  + "_2_sequence.txt.gz")
    return {"fastq1" : f1, "fastq2" : f2}

# pipeline rules ---------------------------------------------------------------

# run entire pipeline for all samples specified in config file
rule all:
  input: expand("{sample}/report.html", sample = config["samples"])

# convert fastq input files into one unmapped bam file
rule fastq_to_bam:
  input:
    unpack(get_fastq_files)
  output:
    temp("{sample}/unmapped.bam")
  log:
    "{sample}/logs/fastq_to_bam.log"
  shell:
    "java -jar {picard_jar} FastqToSam "
    "FASTQ={input.fastq1} "
    "FASTQ2={input.fastq2} "
    "OUTPUT={output} "
    "SAMPLE_NAME={wildcards.sample} "
    "2> {log}"

# tag genome reads with CELL barcodes
rule tag_cell_barcodes:
  input:
    "{sample}/unmapped.bam"
  output:
    bam = temp("{sample}/cell_tagged_unmapped.bam"),
    summary = "{sample}/cell_tags_summary.txt"
  log:
    "{sample}/logs/tag_cell_barcodes.log"
  params:
    base_qual = config["tag_cell_barcodes"]["base_quality"],
    bases_below_qual = config["tag_cell_barcodes"]["num_bases_below_quality"]
  shell:
    "{dropseq_root}/TagBamWithReadSequenceExtended "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "SUMMARY={output.summary} "
    "BASE_QUALITY={params.base_qual} "
    "NUM_BASES_BELOW_QUALITY={params.bases_below_qual} "
    "BASE_RANGE=1-12 "
    "BARCODED_READ=1 "
    "DISCARD_READ=false "
    "TAG_NAME=XC "
    "2> {log}"

# tag genome reads with MOLECULE barcodes
rule tag_molecule_barcodes:
  input:
    "{sample}/cell_tagged_unmapped.bam"
  output:
    bam = temp("{sample}/mol_tagged_unmapped.bam"),
    summary = "{sample}/mol_tags_summary.txt"
  log:
    "{sample}/logs/tag_molecule_barcodes.log"
  params:
    base_qual = config["tag_cell_barcodes"]["base_quality"],
    bases_below_qual = config["tag_cell_barcodes"]["num_bases_below_quality"]
  shell:
    "{dropseq_root}/TagBamWithReadSequenceExtended "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "SUMMARY={output.summary} "
    "BASE_QUALITY={params.base_qual} "
    "NUM_BASES_BELOW_QUALITY={params.bases_below_qual} "
    "BASE_RANGE=13-20 "
    "BARCODED_READ=1 "
    "DISCARD_READ=true "
    "TAG_NAME=XM "
    "2> {log}"

# filter reads marked as 'rejected' by TagBamWithReadSequenceExtended
rule filter_bam:
  input:
    "{sample}/mol_tagged_unmapped.bam"
  output:
    temp("{sample}/filt_unmapped.bam")
  log:
    "{sample}/logs/filter_bam.log"
  shell:
    "{dropseq_root}/FilterBAM "
    "INPUT={input} "
    "OUTPUT={output} "
    "TAG_REJECT=XQ "
    "2> {log}"

# trim SMART adapter sequences from 5'
rule trim_starting_sequence:
  input:
    "{sample}/filt_unmapped.bam"
  output:
    bam = temp("{sample}/adapter_trimmed_unmapped.bam"),
    summary = "{sample}/adapter_trimming_report.txt"
  log:
    "{sample}/logs/trim_starting_sequence.log"
  params:
    adapter_sequence = config["trim_starting_sequence"]["adapter_sequence"],
    mismatches = config["trim_starting_sequence"]["mismatches"],
    num_bases = config["trim_starting_sequence"]["num_bases"]
  shell:
    "{dropseq_root}/TrimStartingSequence "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "OUTPUT_SUMMARY={output.summary} "
    "SEQUENCE={params.adapter_sequence} "
    "MISMATCHES={params.mismatches} "
    "NUM_BASES={params.num_bases} "
    "2> {log}"

# trim polyA sequences from 3'
rule trim_polyA:
  input:
    "{sample}/adapter_trimmed_unmapped.bam"
  output:
    bam = temp("{sample}/polyA_trimmed_unmapped.bam"),
    summary = "{sample}/polyA_trimming_report.txt"
  log:
    "{sample}/logs/trim_polyA.log"
  params:
    mismatches = config["trim_polyA"]["mismatches"],
    num_bases = config["trim_polyA"]["num_bases"]
  shell:
    "{dropseq_root}/PolyATrimmer "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "OUTPUT_SUMMARY={output.summary} "
    "MISMATCHES={params.mismatches} "
    "NUM_BASES={params.num_bases} "
    "2> {log}"

# convert to fastq for STAR read aligner
rule sam_to_fastq:
  input:
    "{sample}/polyA_trimmed_unmapped.bam"
  output:
    temp("{sample}/polyA_trimmed_unmapped.fastq")
  log:
    "{sample}/logs/sam_to_fastq.log"
  shell:
    "java -jar {picard_jar} SamToFastq "
    "INPUT={input} "
    "FASTQ={output} "
    "2> {log}"

# align reads using STAR
rule star_align:
  input:
    "{sample}/polyA_trimmed_unmapped.fastq"
  output:
    temp("{sample}/star.Aligned.out.bam"),
    "{sample}/star.Log.final.out"
  params:
    genomedir = config["star_align"]["genomedir"],
    outprefix = "{sample}/star."
  threads: config["star_align"]["threads"]
  shell:
    "STAR --runThreadN {threads} "
    "--genomeDir {params.genomedir} "
    "--readFilesIn {input} "
    "--outFileNamePrefix {params.outprefix} "
    "--outSAMtype BAM Unsorted ; "
    # move STAR "progress" logs into log directory
    "mv {wildcards.sample}/star.Log.progress.out {wildcards.sample}/logs ; "
    "mv {wildcards.sample}/star.Log.out {wildcards.sample}/logs"


# sort aligned reads
rule sort_aligned:
  input:
    "{sample}/star.Aligned.out.bam"
  output:
    temp("{sample}/star.Aligned.sorted.bam")
  log:
    "{sample}/logs/sort_aligned.log"
  shell:
    "java -jar {picard_jar} SortSam "
    "INPUT={input} "
    "OUTPUT={output} "
    "SORT_ORDER=queryname "
    "2> {log}"

# merge aligned and unaligned reads to add tags to aligned reads
rule merge_bam:
  input:
    aligned = "{sample}/star.Aligned.sorted.bam",
    unaligned = "{sample}/polyA_trimmed_unmapped.bam"
  output:
    temp("{sample}/merged_aligned.bam")
  log:
    "{sample}/logs/merge_bam.log"
  params:
    reference = config["merge_bam"]["reference"]
  shell:
    "java -jar {picard_jar} MergeBamAlignment "
    "ALIGNED_BAM={input.aligned} "
    "UNMAPPED_BAM={input.unaligned} "
    "REFERENCE_SEQUENCE={params.reference} "
    "OUTPUT={output} "
    "INCLUDE_SECONDARY_ALIGNMENTS=false "
    "PAIRED_RUN=false "
    "2> {log}"

# tag reads with gene exons
rule tag_with_gene_exon:
  input:
    "{sample}/merged_aligned.bam"
  output:
    temp("{sample}/gene_tagged_aligned.bam")
  log:
    "{sample}/logs/tag_with_gene_exon.log"
  params:
    annot = config["tag_with_gene_exon"]["annotations_file"]
  shell:
    "{dropseq_root}/TagReadWithGeneExon "
    "INPUT={input} "
    "OUTPUT={output} "
    "ANNOTATIONS_FILE={params.annot} "
    "TAG=GE "
    "2> {log}"

# filter UMI error
rule bead_synthesis_error:
  input:
    "{sample}/gene_tagged_aligned.bam"
  output:
    bam = "{sample}/filt_gene_tagged_aligned.bam",
    summary = "{sample}/bead_synthesis_error_summary.txt",
    detail = "{sample}/bead_synthesis_error_detail.txt"
  log:
    "{sample}/logs/filter_synthesis_error.log"
  params:
    num_bcs = lambda wildcards:
      config["expect_cell_numbers"][wildcards.sample] * 2,
    min_umis_per_cell = config["bead_synthesis_error"]["min_umis_per_cell"],
    max_num_errors = config["bead_synthesis_error"]["max_num_errors"],
    read_mq = config["bead_synthesis_error"]["read_mq"],
    primer_sequence = config["bead_synthesis_error"]["primer_sequence"],
    edit_distance = config["bead_synthesis_error"]["edit_distance"]
  shell:
    "{dropseq_root}/DetectBeadSynthesisErrors "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "NUM_BARCODES={params.num_bcs} "
    "SUMMARY={output.summary} "
    "OUTPUT_STATS={output.detail} "
    "MIN_UMIS_PER_CELL={params.min_umis_per_cell} "
    "MAX_NUM_ERRORS={params.max_num_errors} "
    "READ_MQ={params.read_mq} "
    "PRIMER_SEQUENCE={params.primer_sequence} "
    "EDIT_DISTANCE={params.edit_distance} "
    "CREATE_INDEX=true "
    "2> {log}"

# calculate reads per cell barcode
rule reads_per_cell:
  input:
    "{sample}/filt_gene_tagged_aligned.bam"
  output:
    "{sample}/reads_per_cell_barcode.txt"
  log:
    "{sample}/logs/reads_per_cell.log"
  shell:
    "{dropseq_root}/BAMTagHistogram "
    "INPUT={input} "
    "OUTPUT={output} "
    "TAG=XC "
    "2> {log}"

# compile R Markdown report in html format
rule report:
  input:
    cell_bcs = "{sample}/cell_tags_summary.txt",
    mol_bcs = "{sample}/mol_tags_summary.txt",
    star_smry = "{sample}/star.Log.final.out",
    adapt_trim = "{sample}/adapter_trimming_report.txt",
    polyA_trim = "{sample}/polyA_trimming_report.txt",
    synthesis_error = "{sample}/bead_synthesis_error_summary.txt",
    reads_per_cell = "{sample}/reads_per_cell_barcode.txt"
  output:
    "{sample}/report.html"
  params:
    expect_cells = lambda wildcards:
      config["expect_cell_numbers"][wildcards.sample]
  script:
    "scripts/report.Rmd"
