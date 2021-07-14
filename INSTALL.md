#Installation instructions for DAS Tool version 1.1

# 1. Dependencies

- R (>= 3.2.3): https://www.r-project.org
- R-packages: data.table (>= 1.9.6), doMC (>= 1.3.4)
- ruby (>= v2.3.1): https://www.ruby-lang.org
- Pullseq (>= 1.0.2): https://github.com/bcthomas/pullseq
- Prodigal (>= 2.6.3): https://github.com/hyattpd/Prodigal
- coreutils (only macOS/ OS X): https://www.gnu.org/software/coreutils
- One of the following search engines:
	- USEARCH (>= 8.1): http://www.drive5.com/usearch/download.html
	- DIAMOND (>= 0.9.14): https://ab.inf.uni-tuebingen.de/software/diamond
	- BLAST+ (>= 2.5.0): https://blast.ncbi.nlm.nih.gov/Blast.cgi



# 2. Installation

## Github
```
# Download and extract DASTool.zip archive:
unzip DAS_Tool-1.x.x.zip
cd ./DAS_Tool-1.x.x

# Install R-packages:
R CMD INSTALL ./package/DASTool_1.x.x.tar.gz

# Unzip SCG database:
unzip ./db.zip

# Call DAS Tool:
./DAS_Tool -h
```

## Bioconda

 Bioconda repository: https://bioconda.github.io/recipes/das_tool/README.html. Thanks @[keuv-grvl]("https://github.com/keuv-grvl") and @[silask]("https://github.com/SilasK")!.

Add bioconda channel:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Install DAS Tool using conda:
```
conda install -c bioconda das_tool
```

## Homebrew

Homebrew-bio repository: https://github.com/brewsci/homebrew-bio. Thanks @[gaberoo]("https://github.com/gaberoo")!

Install DAS Tool using Homebrew:
```
brew install brewsci/bio/das_tool
```


For detailed instructions please read the documentation.
