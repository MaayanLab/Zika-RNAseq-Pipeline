# An Open RNA-Seq Data Analysis Pipeline with an Example of Reprocessing Data from a Recent Zika Virus Studyfrom NCBI's GEO

_Zichen Wang_ and _Avi Ma'ayan_

_BD2K-LINCS Data Coordination and Integration Center (DCIC)_
_Icahn School of Medicine at Mount Sinai, New York, NY 10029 USA_

[![DOI](https://zenodo.org/badge/22891/MaayanLab/Zika-RNAseq-Pipeline.svg)](https://zenodo.org/badge/latestdoi/22891/MaayanLab/Zika-RNAseq-Pipeline)

## Guide to downloading and running this Docker image

### Step 1: Install Docker

To download and run this Docker image, you first need to set up Docker on your machine. The easiest way to start with Docker is to install the [Docker Toolbox](https://www.docker.com/products/docker-toolbox) by simply downloading and clicking the installer which is available for both Mac OSX and Windows. For Linux users, follow the instructions [here](https://docs.docker.com/linux/step_one/). 

### Step 2: Download and run the Docker image

#### Option 1: Through a Command Line Interface (CLI)

The image can be downloaded and executed through the CLI of Docker's `Docker Quickstart Terminal` in the [Docker Toolbox](https://www.docker.com/products/docker-toolbox) with the following commands:

1. Pull(download) the Docker image:
	```
	$ docker pull maayanlab/zika
	```
2. Run the Docker image. The Docker container requires host to mount two directories as data volumes: the reference genome directory (`/notebook/genomes`) and the data directory (`/notebook/data`). This can be done by specifying the -v tag when running Docker:
	`-v /host/path/to/genomes:/notebook/genomes -v /host/path/to/data:/notebook/data`
	```
	$ docker run -d -p 80:8888 -e "PASSWORD=YourPassword" -e "USE_HTTP=1" -v /host/path/to/genomes:/notebook/genomes -v /host/path/to/data:/notebook/data maayanlab/zika
	```
3. Get the IP of your Docker machine:   
	```
	$ docker-machine ip
	```
4. Open a browser and go to http://your.docker-machine.ip and enter the password you set to run the RNA-Seq pipeline. 

More detailed instructions on how to open the Docker Quickstart Terminal are available for [Mac OSX](https://docs.docker.com/mac/step_one/) and [Windows](https://docs.docker.com/windows/step_one/).


#### Option 2: Through the Graphical User Interface (GUI)

**Note:** Kitematic currently does not support mounting host directories as data volumes of the Docker container. Therefore it is suggested to use CLI to run the Docker image if you need to analyze new data with this pipeline.

1. Open [Kitematic](https://www.docker.com/products/docker-kitematic)
2. Search for maayanlab/zika in the search box and download the Docker image
3. Set the following variables under the 'Setting' tab:
	+ Environment Varables: 
		+ PASSWORD: your password
		+ USE_HTTP: 1
	+ Docker port: 8888
4. Click 'START' then click the maximize button on the top right corner of 'WEB PREVIEW'

## Deploy this Docker image onto your cloud

You can download and deploy this Docker image with your cloud provider such as [Amazon Web Services](https://www.docker.com/aws), [HP Enterprise](https://www.docker.com/aws), [IBM](https://www.docker.com/IBM), [Microsoft Azure Cloud](https://www.docker.com/microsoft) or others.

