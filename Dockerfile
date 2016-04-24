FROM ipython/scipystack

RUN pip install -r requirements.txt 

# Update the repositories
RUN \
  apt-get update -qq

# Install R
RUN \
  apt-get install -y r-base r-base-dev && \
  rm -rf /var/lib/apt/lists/*

RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("edgeR");'

# Copy the application folder inside the container
ADD . /notebook

# Set the default directory where CMD will execute
WORKDIR /notebook

# Install SRA-tookit, STAR, featureCounts
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.6.2/sratoolkit.2.6.2-ubuntu64.tar.gz
RUN tar zxvf sratoolkit.2.6.2-ubuntu64.tar.gz

RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.1c.tar.gz
RUN tar zxvf STAR_2.4.1c.tar.gz
RUN cd STAR_2.4.1c && make STAR

RUN wget http://downloads.sourceforge.net/project/subread/subread-1.5.0-p2/subread-1.5.0-p2-Linux-x86_64.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.5.0-p2%2F&ts=1461463153&use_mirror=pilotfiber
RUN tar zxvf subread-1.5.0-p2-Linux-x86_64.tar
