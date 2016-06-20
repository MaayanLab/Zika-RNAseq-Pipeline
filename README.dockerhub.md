# RNA-Seq data analysis pipeline: an example of processing a Zika virus study from NCBI's GEO

_Zichen Wang_ and _Avi Ma'ayan_

## Guide to downloading and running this Docker image

### Step 1: Install Docker

To download and run this Docker image, you first need to set up Docker on your machine. The easiest way to start with Docker is to install [Docker Toolbox](https://www.docker.com/products/docker-toolbox) by simply downloading and click the installer which is available for both Mac OSX and Windows. For Linux users, follow the instruction [here](https://docs.docker.com/linux/step_one/) to set up Docker on your machine. 

### Step 2: Download and run the Docker image

#### Option 1: Through Command Line Interface (CLI)

The image can be downloaded and executed through the CLI of Docker's `Docker Quickstart Terminal` in the [Docker Toolbox](https://www.docker.com/products/docker-toolbox) with the following commands:

1. Pull(download) the Docker image:
	```
	$ docker pull maayanlab/zika
	```
2. Run the Docker image
	```
	$ docker run -d -p 80:8888 -e "PASSWORD=YourPassword" -e "USE_HTTP=1" maayanlab/zika
	```
3. Get the IP of your Docker machine:   
	```
	$ docker-machine ip
	```
4. Open a browser and go to http://your.docker-machine.ip and enter the password you set to run the RNA-seq pipeline. 

More detailed instructions on how to open the Docker Quickstart Terminal is available for [Mac OSX](https://docs.docker.com/mac/step_one/) and [Windows](https://docs.docker.com/windows/step_one/).


#### Option 2: Through Graphical User Interface (GUI)

1. Open [Kitematic](https://www.docker.com/products/docker-kitematic)
2. Search for maayanlab/zika in the search box and download the Docker image
3. Set the following variables under the 'Setting' tab:
	+ Environment Varables: 
		+ PASSWORD: your password
		+ USE_HTTP: 1
	+ Docker port: 8888
4. Click 'START' then click the maximize button on the top right corner of 'WEB PREVIEW'

## Deploy this Docker image onto your cloud

You can downloand and deploy this Docker with your cloud provider such as [Amazon Web Services](https://www.docker.com/aws), [HP Enterprise](https://www.docker.com/aws), [IBM](https://www.docker.com/IBM), [Microsoft Azure Cloud](https://www.docker.com/microsoft) and etc.

