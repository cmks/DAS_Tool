#Installation instructions for DAS Tool version 1.1

# 1. Dependencies

- R (>= 3.2.3): https://www.r-project.org
- R-packages: data.table (>= 1.9.6), doMC (>= 1.3.4), ggplot2 (>= 2.1.0)
- ruby (>= v2.3.1): https://www.ruby-lang.org
- Pullseq (>= 1.0.2): https://github.com/bcthomas/pullseq
- Prodigal (>= 2.6.3): https://github.com/hyattpd/Prodigal
- coreutils (only macOS/ OS X): https://www.gnu.org/software/coreutils
- One of the following search engines:
	- USEARCH (>= 8.1): http://www.drive5.com/usearch/download.html
	- BLAST+ (>= 2.5.0): https://blast.ncbi.nlm.nih.gov/Blast.cgi


# 2. Installation

``` 
# Download and extract DASTool.zip archive:
unzip DAS_Tool.v1.1.zip
cd ./DAS_Tool.v1.1

# Install R-packages:
R CMD INSTALL ./package/DASTool_1.1.0.tar.gz

# Download SCG database:
wget http://banfieldlab.berkeley.edu/~csieber/db.zip
unzip db.zip

# Call DAS Tool:
./DAS_Tool -h
``` 
For detailed instructions please read the documentation.
