FROM ipython/scipystack

MAINTAINER Zichen Wang <zichen.wang@mssm.edu>

# Copy the application folder inside the container
ADD . /notebook

# Set the default directory where CMD will execute
WORKDIR /notebook
# Set environment variable
ENV HOME /notebook

# Install additional python packages
RUN pip2 install -r requirements.txt 

# Install wget and unzip
RUN apt-get update -qq && apt-get install -y \
	wget \
	unzip \
	# Install java for fastQC
	default-jre \
	# Install R
	r-base \
	r-base-dev \
&& rm -rf /var/lib/apt/lists/*

# Install R pacakges
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("edgeR");'


# Install SRA-tookit, fastQC, kallisto
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.6.2/sratoolkit.2.6.2-ubuntu64.tar.gz
RUN tar zxvf sratoolkit.2.6.2-ubuntu64.tar.gz
RUN ln -s /notebook/sratoolkit.2.6.2-ubuntu64/bin/* /usr/local/bin/

RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
RUN unzip fastqc_v0.11.5.zip
RUN chmod 755 FastQC/fastqc && ln -s /notebook/FastQC/fastqc /usr/local/bin/

RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
RUN tar zxvf kallisto_linux-v0.44.0.tar.gz
RUN ln -s /notebook/kallisto_linux-v0.44.0/kallisto /usr/local/bin/


# Clean-ups
RUN rm *.gz && rm *.zip

# Expose port
EXPOSE 8888

ADD notebook.sh /

# Start notebook server
CMD ["/notebook.sh"]
