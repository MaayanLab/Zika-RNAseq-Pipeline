FROM ipython/scipystack

# Copy the application folder inside the container
ADD . /notebook

# Set the default directory where CMD will execute
WORKDIR /notebook

# Install additional python packages
RUN pip2 install -r requirements.txt 

# Install wget
RUN apt-get install -y wget

# 
RUN apt-get update -qq

RUN apt-get install -y unzip

# Install java for fastQC
RUN apt-get install -y default-jre


# Install R
RUN \
  apt-get install -y r-base r-base-dev && \
  rm -rf /var/lib/apt/lists/*

RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("edgeR");'



# Install SRA-tookit, fastQC STAR, featureCounts
RUN cd /home
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.6.2/sratoolkit.2.6.2-ubuntu64.tar.gz
RUN tar zxvf sratoolkit.2.6.2-ubuntu64.tar.gz
RUN ln -s /home/sratoolkit.2.6.2-ubuntu64/bin/* /usr/local/bin/

RUN cd /home
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
RUN unzip fastqc_v0.11.5.zip
RUN chmod 755 FastQC/fastqc && ln -s /home/FastQC/fastqc /usr/local/bin/


RUN cd /home
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.1c.tar.gz
RUN tar zxvf STAR_2.4.1c.tar.gz
RUN cd STAR-STAR_2.4.1c/source && make STAR
RUN ln -s /home/STAR-STAR_2.4.1c/bin/Linux_x86_64/STAR /usr/local/bin/

RUN tar zxvf subread-1.4.6-p2-Linux-x86_64.tar.gz

# Expose port
EXPOSE 8888

ADD notebook.sh /

# Start notebook server
CMD ["/notebook.sh"]
