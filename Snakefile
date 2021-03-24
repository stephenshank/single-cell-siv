import os

import pandas as pd
from Bio import SeqIO
from Bio import Entrez


wildcard_constraints:
  reads="[^/]+",
  barcode="[^/]+",

macaque_accessions = [
  "NC_027914.1",
  "NC_041754.1",
  "NC_041755.1",
  "NC_041756.1",
  "NC_041757.1",
  "NC_041758.1",
  "NC_041759.1",
  "NC_041760.1",
  "NC_041761.1",
  "NC_041762.1",
  "NC_041763.1",
  "NC_041764.1",
  "NC_041765.1",
  "NC_041766.1",
  "NC_041767.1",
  "NC_041768.1",
  "NC_041769.1",
  "NC_041770.1",
  "NC_041771.1",
  "NC_041772.1",
  "NC_041773.1",
  "NC_041774.1"
]

with open('data/input/cells.txt') as f:
  barcodes = [line.strip() for line in f.readlines()]

rule unzip_reads:
  input:
    "data/input/{reads}.fastq.gz"
  output:
    "data/{reads}.fastq"
  shell:
    "gunzip -c {input} > {output}"

rule fastq_to_fasta:
  input:
    rules.unzip_reads.output[0]
  output:
    'data/{reads}.fasta'
  run:
    SeqIO.write(
      SeqIO.parse(input[0], 'fastq'),
      output[0],
      'fasta'
    )

rule ncbi_download:
  output:
    ['data/macaque/%s.fasta' % accession for accession in macaque_accessions]
  run:
    Entrez.email = os.environ.get('email') or 'sshank@temple.edu'
    for accession in macaque_accessions:
        record = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        filename = 'data/macaque/{}.fasta'.format(accession)
        print('Writing:{}'.format(filename))
        with open(filename, "w") as f:
            f.write(record.read())

rule full_macaque:
  input:
    rules.ncbi_download.output
  output:
    'data/macaque/genome.fasta'
  shell:
    'cat data/macaque/NC_*.fasta > {output}'

rule macaque_blast_db:
  input:
    rules.full_macaque.output[0]
  output:
    "data/macaque/blast.nhr",
    "data/macaque/blast.nin",
    "data/macaque/blast.nsq"
  shell:
    "makeblastdb -in {input} -dbtype nucl -out data/macaque/blast"

rule make_header:
  output:
    "data/header.csv"
  shell:
    "echo qaccver,saccver,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore > {output}"

rule macaque_blast:
  input:
    fasta=rules.fastq_to_fasta.output[0],
    macaque_db=rules.macaque_blast_db.output
  output:
    no_header=temp("data/{reads}-blast-no_header.csv"),
    header="data/{reads}-blast.csv",
  shell:
    """
      blastn -db data/macaque/blast -outfmt 10 -query {input.fasta} -out {output.no_header}
      cat {input} {output.no_header} > {output.header}
    """
