# MetaDecoder

An algorithm for clustering metagenomic sequences.

Cite [MetaDecoder](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01237-8):

Liu, CC., Dong, SS., Chen, JB. et al. MetaDecoder: a novel method for clustering metagenomic contigs. Microbiome 10, 46 (2022). [https://doi.org/10.1186/s40168-022-01237-8](https://doi.org/10.1186/s40168-022-01237-8)

## Benchmarks

All datasets mentioned in text and some NEWLY ADDED datasets are available in **[Google Drive](https://drive.google.com/drive/folders/1_mybcewf3VE-7dte6oA-vDmlRx2ugzyD?usp=sharing)**.

Benchmarks for all datasets are available in **benchmarks** directory.

## Dependencies

* [python (version >= 3.8)](https://www.python.org/)
* [numpy](https://pypi.org/project/numpy/)
* [scipy](https://pypi.org/project/scipy/)
* [scikit-learn](https://pypi.org/project/scikit-learn/)
* [threadpoolctl](https://pypi.org/project/threadpoolctl/)
* [fraggenescan (version 1.31)](https://sourceforge.net/projects/fraggenescan/)
* [hmmer (version 3.2.1)](http://www.hmmer.org/)

## Installation

### Download and install MetaDecoder (Do not clone this repository)

```shell
# You may need to install pip3 before. #
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py

# You may need to install or upgrade setuptools and wheel using pip3 before. #
pip3 install --upgrade setuptools wheel

# Download and install MetaDecoder version 1.0.15 #
pip3 install https://github.com/liu-congcong/MetaDecoder/releases/download/v1.0.15/metadecoder-1.0.15-py3-none-any.whl
```

Make sure you have a good internet connection, as MetaDecoder will install the required Python dependencies, you can also install the dependencies manually:

* **numpy**
* **scipy**
* **scikit-learn**
* **threadpoolctl**

```shell
pip3 install numpy scipy scikit-learn threadpoolctl
```

MetaDecoder uses **FragGeneScan** and **Hmmer** for predicting protein coding genes and mapping single-copy marker genes to contigs, respectively.

MetaDecoder has included the compiled FragGeneScan (version 1.31) and Hmmer (version 3.2.1).

### The GPU version of MetaDecoder

MetaDecoder can be accelerated using the GPU on the basis of **CUDA** and **[CuPy](https://cupy.dev/)**.

To use the GPU version of Metadecoder, you need to have a compatible driver installed for your GPU (**CUDA**), and then install CuPy using pip3:

```shell
# You may need to install and upgrade setuptools and wheel using pip3 before. #
pip3 install --upgrade setuptools wheel
# Please note that XXX is the CUDA version. e.g. cupy-cuda101 means CuPy with CUDA 10.1. #
pip3 install cupy-cudaXXX
```

Please be careful not to install multiple CuPy packages at the same time.

MetaDecoder will automatically enable GPU if it is available.

And moreover, CuPy can use additional CUDA library (**cuTENSOR**) to accelerate tensor operations: **UNTESTED**

```shell
# Please note that XXX is the CUDA version. #
python3 -m cupyx.tools.install_library --cuda XXX --library cutensor
# Setting the environment variable to activate some CUDA features in CuPy. #
echo 'export CUPY_ACCELERATORS="cutensor"' >> ~/.bashrc
```

### Run MetaDecoder on MacOS (only for MacOS 11 now)

We provided the compiled fraggenescan and hmmsearch for MacOS users.

To run MetaDecoder on a Mac, please:

* Install MetaDecoder.

* Download and rename **fraggenescan (version 1.31 on MacOS)** and **hmmsearch (version 3.2.1 on MacOS)** to **fraggenescan** and **hmmsearch**, respectively.

```shell
mv "fraggenescan (version 1.31 on MacOS)" fraggenescan
mv "hmmsearch (version 3.2.1 on MacOS)" hmmsearch
```

* Move them to the metadecoder folder.

```shell
mv fraggenescan hmmsearch /Library/Frameworks/Python.framework/Versions/3.*/lib/python3.*/site-packages/metadecoder/
```

* Assign appropriate permissions to the following files such as 755.

```shell
chmod 755 /Library/Frameworks/Python.framework/Versions/3.*/lib/python3.*/site-packages/metadecoder/fraggenescan
chmod 755 /Library/Frameworks/Python.framework/Versions/3.*/lib/python3.*/site-packages/metadecoder/hmmsearch
```

## Usage

### Preparations

**Before running MetaDecoder, you may need to prepare some files by yourself.**

* A FASTA formatted assembly file: **ASSEMBLY.FASTA**

* Some SAM formatted read files with the same headers: **SAMPLE1.SAM**, **SAMPLE2.SAM** ...

### Run MetaDecoder

#### Obtain the coverages of contigs

Input: **SAMPLE1.SAM**, **SAMPLE2.SAM**, **...**

Output: **METADECODER.COVERAGE**

```shell
metadecoder coverage -s SAMPLE1.SAM SAMPLE2.SAM ... -o METADECODER.COVERAGE
```

#### Map single-copy marker genes to the assembly

Input: **ASSEMBLY.FASTA**

Output: **METADECODER.SEED**

```shell
metadecoder seed --threads 50 -f ASSEMBLY.FASTA -o METADECODER.SEED
```

#### Run MetaDecoder algorithm to cluster contigs

Input: **ASSEMBLY.FASTA**, **METADECODER.COVERAGE**, **METADECODER.SEED**

Output: **METADECODER.1.FASTA**, **METADECODER.2.FASTA**, ...

```shell
metadecoder cluster -f ASSEMBLY.FASTA -c METADECODER.COVERAGE -s METADECODER.SEED -o METADECODER
```

### A simple example to use MetaDecoder is available in MetaDecoder/example/

## Change logs

* 1.0.3 (20211008): Initial version.

* 1.0.4 (20211029): Optimize the calculation of distance of pairwise kmer frequency.

* 1.0.5 (20211105): Optimize the counting process of kmers.

* 1.0.6 (20211207): Added an option (--no_clusters) to output only sequence IDs and the corresponding cluster IDs instead of sequences.

* 1.0.7 (20220206): Minor bugs fixed.

* 1.0.8 (20220315): Minor bugs fixed.

* 1.0.9 (20220508): Fix a bug that causes abnormal coverage when the "--mapq" parameter is set to 0.

* 1.0.10 (20220514): Minor bugs fixed.

* 1.0.11 (20220518): Minor bugs fixed.

* 1.0.12 (20220708): Minor bugs fixed.

* 1.0.13 (20220712): Minor bugs fixed.

* 1.0.14 (20220914): Minor bugs fixed.

* 1.0.15 (20221103): Minor bugs fixed.

## References

* Mina Rho, Haixu Tang, and Yuzhen Ye. FragGeneScan: Predicting Genes in Short and Error-prone Reads. Nucl. Acids Res., 2010 doi: 10.1093/nar/gkq747.

* nhmmer: DNA Homology Search With Profile HMMs. T. J. Wheeler, S. R. Eddy. Bioinformatics, 29:2487-2489, 2013.
