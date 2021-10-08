# MetaDecoder

An algorithm for clustering metagenomic sequences.

## Description

### Five reference genomes

* **A**cidithiobacillus_**f**errivorans.fasta.gz
* **B**eijerinckia_**i**ndica.fasta.gz
* **C**yanobacterium_**a**poninum.fasta.gz
* **D**actylococcopsis_**s**alina.fasta.gz
* **E**ggerthella_**l**enta.fasta.gz

### Two simulated datasets

* Af3Bi6Ca4Ds9El2.*.fastq.gz with coverages **Af** 3X **Bi** 6X **Ca** 4X **Ds** 9X **El** 2X
* Af6Bi2Ca5Ds7El8.*.fastq.gz with coverages **Af** 6X **Bi** 2X **Ca** 5X **Ds** 7X **El** 8X

* Af3Bi6Ca4Ds9El2.sam.gz (map Af3Bi6Ca4Ds9El2.fastq to the assembly)
* Af6Bi2Ca5Ds7El8.sam.gz (map Af6Bi2Ca5Ds7El8.fastq to the assembly)

### The assembly file

* AfBiCaDsEl.fasta.gz

## Run MetaDecoder

* Download **AfBiCaDsEl.fasta.gz**, **Af3Bi6Ca4Ds9El2.sam.gz**, and **Af6Bi2Ca5Ds7El8.sam.gz** from **[Google Drive](https://drive.google.com/drive/folders/10a3e3W3ei6OV8d5aIpzhvrL1sNdjzwRg?usp=sharing)** and decompress them.

* Obtain the coverages of contigs.

```shell
metadecoder coverage -s Af3Bi6Ca4Ds9El2.sam Af6Bi2Ca5Ds7El8.sam -o AfBiCaDsEl.metadecoder.coverage
```

* Map single-copy marker genes to the assembly.

```shell
metadecoder seed --threads 50 -f AfBiCaDsEl.fasta -o AfBiCaDsEl.metadecoder.seed
```

* Run MetaDecoder algorithm to cluster contigs.

```shell
metadecoder cluster --min_sequence_length 1000 -f AfBiCaDsEl.fasta -c AfBiCaDsEl.metadecoder.coverage -s AfBiCaDsEl.metadecoder.seed -o AfBiCaDsEl.metadecoder
```
