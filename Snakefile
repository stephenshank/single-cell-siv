import pandas as pd
from Bio import SeqIO


with open('data/cells.txt') as f:
  barcodes = [line.strip() for line in f.readlines()]

rule fastq_to_fasta:
  input:
    'data/{dataset}.fastq'
  output:
    'data/{dataset}.fasta'
  run:
    SeqIO.write(
      SeqIO.parse(input[0], 'fastq'),
      output[0],
      'fasta'
    )

rule blast:
  input:
    rules.fastq_to_fasta.output[0]
  output:
    "data/{dataset}-blast.csv"
  shell:
    "blastn -db /Volumes/Data/macaca-mulatta/blast/db -outfmt 10 -query {input} -out {output}"

rule barcode_fasta:
  input:
    rules.fastq_to_fasta.output[0]
  output:
    'data/{barcode}/{dataset}.fasta'
  shell:
    'seqkit grep --by-seq --max-mismatch 1 --pattern "{wildcards.barcode}" {input} > {output}'

rule barcode_csv:
  input:
    all_records=rules.fastq_to_fasta.output[0],
    barcoded_records=rules.barcode_fasta.output[0]
  output:
    'data/{barcode}/{dataset}.csv'
  run:
    headers = [record.id for record in SeqIO.parse(input.all_records, 'fasta')]
    barcode_hash = SeqIO.to_dict(SeqIO.parse(input.barcoded_records, 'fasta'))
    contains_barcode = [header in barcode_hash for header in headers]
    pd.DataFrame({wildcards.barcode: contains_barcode}).to_csv(output[0])

rule all_barcodes:
  input:
    all_barcodes=expand('data/{barcode}/{{dataset}}.csv', barcode=barcodes),
    all_records=rules.fastq_to_fasta.output[0],
  output:
    'data/{dataset}-barcodes.csv'
  run:
    headrs = [record.id for record in SeqIO.parse(input.all_records, 'fasta')]
    pd.concat([
      pd.read_csv('data/%s/%s.csv' % (barcode, wildcards.dataset))[barcode]
      for barcode in barcodes
    ], axis=1).to_csv(output[0])
