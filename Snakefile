import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

constant_region = 'AGTACGTACGAGTC'


rule all:
  input:
    "data/barcodes/reads-blast.csv",
    "data/barcodes/mates-blast.csv",
    "data/macaque/reads-blast.csv",
    "data/macaque/mates-blast.csv"

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
    macaque_db=rules.macaque_blast_db.output,
    header=rules.make_header.output[0]
  output:
    no_header=temp("data/macaque/{reads}-blast-no_header.csv"),
    header="data/macaque/{reads}-blast.csv"
  shell:
    """
      blastn -db data/macaque/blast -outfmt 10 -query {input.fasta} -word_size 64 -evalue 1000 -out {output.no_header}
      cat {input.header} {output.no_header} > {output.header}
    """

rule cells_fasta:
  input:
    "data/input/cells.txt"
  output:
    "data/cells.fasta"
  run:
    with open(input[0]) as cells_file:
      barcodes = [line.strip() for line in cells_file.readlines()]
    SeqIO.write(
      [
        SeqRecord(
          Seq(barcode[:9] + constant_region + barcode[-9:]),
          id='barcode_%d' % i,
          description=''
        )
        for i, barcode in enumerate(barcodes)
      ],
      output[0],
      'fasta'
    )

rule barcode_blast_db:
  input:
    rules.cells_fasta.output[0]
  output:
    "data/barcodes/blast.nhr",
    "data/barcodes/blast.nin",
    "data/barcodes/blast.nsq"
  shell:
    "makeblastdb -in {input} -dbtype nucl -out data/barcodes/blast"

rule extract_best_hits:
  input:
    rules.macaque_blast.output.header
  output:
    "data/macaque/{reads}-besthits.csv"
  run:
    all_hits = pd.read_csv(input[0])
    best_hit_indices = all_hits.groupby('qaccver').evalue.idxmin()
    all_hits.loc[best_hit_indices].to_csv(output[0])

rule barcode_blast:
  input:
    fasta=rules.fastq_to_fasta.output[0],
    blast_db=rules.barcode_blast_db.output,
    header=rules.make_header.output[0]
  output:
    no_header=temp("data/barcodes/{reads}-blast-no_header.csv"),
    header="data/barcodes/{reads}-blast.csv"
  shell:
    """
      blastn -db data/barcodes/blast -outfmt 10 -query {input.fasta} -word_size 16 -evalue 1000 -out {output.no_header}
      cat {input.header} {output.no_header} > {output.header}
    """
