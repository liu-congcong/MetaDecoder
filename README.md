# MetaDecoder

An algorithm for clustering metagenomic sequences.

All datasets mentioned in text and some NEWLY ADDED datasets are available in **[Google Drive](https://drive.google.com/drive/folders/1_mybcewf3VE-7dte6oA-vDmlRx2ugzyD?usp=sharing)**.

Benchmarks for all datasets are available in **benchmarks** directory.

## Dependencies

* [python (version 3.8.6 or >= 3.8)](https://www.python.org/)
* [numpy (version 1.18.5)](https://pypi.org/project/numpy/)
* [scipy (version 1.5.4)](https://pypi.org/project/scipy/)
* [scikit-learn (version 0.23.2)](https://pypi.org/project/scikit-learn/)
* [threadpoolctl](https://pypi.org/project/threadpoolctl/)
* [fraggenescan (version 1.31)](https://sourceforge.net/projects/fraggenescan/) | [prodigal (version 2.6.3)](https://github.com/hyattpd/Prodigal/)
* [hmmer (version 3.2.1)](http://www.hmmer.org/)

## Installation

### Download and install MetaDecoder (Do not clone this repository)

```shell
# You may need to install pip3 before. #
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py

# You may need to install or upgrade setuptools and wheel using pip3 before. #
pip3 install --upgrade setuptools wheel

# Download and install MetaDecoder (MD5: dda13f738dfedcdfc226232b089633c9) #
wget https://github.com/liu-congcong/MetaDecoder/releases/download/v1.0.3/metadecoder.1.0.3.tar.gz
pip3 install metadecoder*.tar.gz
```

Make sure you have a good internet connection, as MetaDecoder will install the required Python dependencies, you can also install the dependencies manually:

* **numpy (version 1.18.5)**
* **scipy (version 1.5.4)**
* **scikit-learn (version 0.23.2)**
* **threadpoolctl**

```shell
pip3 install numpy==1.18.5 scipy==1.5.4 scikit-learn==0.23.2 threadpoolctl
```

MetaDecoder uses **FragGeneScan** or **Prodigal** and **Hmmer** for predicting protein coding genes and mapping single-copy marker genes to contigs, respectively.

MetaDecoder has included the compiled FragGeneScan (version 1.31) and Hmmer (version 3.2.1).

If any errors occurs during the process, you may need to compile the following programs manually:

* **FragGeneScan (version 1.31)** or **Prodigal (version 2.6.3)**

* **Hmmer (version 3.2.1)**

If you encounter any installation problems, please:

* Download and decompress **metadecoder.*.tar.gz**.

* Move **metadecoder** from **bin** folder to **python3.X/bin/**.

```shell
mv metadecoder-1.0/bin/metadecoder python3.X/bin/
```

* Move the **metadecoder** folder to **python3.X/lib/python3.X/site-packages/**.

```shell
mv metadecoder-1.0/metadecoder python3.X/lib/python3.X/site-packages/
```

* Assign appropriate permissions to the following files such as 755.

```shell
chmod 755 python3.X/bin/metadecoder
chmod 755 python3.X/lib/python3.X/site-packages/metadecoder/fraggenescan
chmod 755 python3.X/lib/python3.X/site-packages/metadecoder/hmmsearch
```

To verify your installation, please run **python3 -c "import metadecoder"**.

If there is no error, you have successfully installed MetaDecoder.

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

* Some **UNSORTED** (the raw output by BWA/Bowtie2/...) SAM formatted aligned reads files: **SAMPLE1.SAM**, **SAMPLE2.SAM** ...

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

* To do: Add GPU version KMeans++ in DPGMM.

## References

* Mina Rho, Haixu Tang, and Yuzhen Ye. FragGeneScan: Predicting Genes in Short and Error-prone Reads. Nucl. Acids Res., 2010 doi: 10.1093/nar/gkq747.

* Hyatt, Doug & Chen, G.L. & Locascio, Phil & Land, Miriam & Larimer, F.W. & Hauser, Loren. (2010). Prodigal prokaryotic dynamic programming genefinding algorithm. BMC Bioinformatics. 11.

* nhmmer: DNA Homology Search With Profile HMMs. T. J. Wheeler, S. R. Eddy. Bioinformatics, 29:2487-2489, 2013.
