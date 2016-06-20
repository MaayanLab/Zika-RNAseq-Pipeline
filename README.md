# RNA-Seq data analysis pipeline: an example of processing a Zika virus study from NCBI's GEO

_Zichen Wang_ and _Avi Ma'ayan_

## Browse the notebook

To view the IPython notebook with richer display, click [here](http://nbviewer.jupyter.org/github/maayanlab/Zika-RNAseq-Pipeline/blob/master/Zika.ipynb)

## Run the RNA-seq pipeline interactively

We have created a Docker image ([maayanlab/zika](https://hub.docker.com/r/maayanlab/zika/)) packaging all the dependencies (command line tools, R and Python packages) for the pipeline, which is publically available on [Dockerhub](https://hub.docker.com/). There are several options to run the Docker image:

1. From our server
	We have deployed the Docker image on our [Mesos](http://mesos.apache.org/) cluster available at: http://isabella.1425mad.mssm.edu:31516/.

2. On your local machine
	1. Through command line
		1. Install Docker Toolbox following the instructions [here](https://www.docker.com/products/docker-toolbox) 
		2. Pull our Docker image from Dockerhub   
			```
			$ docker pull maayanlab/zika
			```
		3. Run the Docker image   
			```
			$ docker run -d -p 80:8888 -e "PASSWORD=YourPassword" -e "USE_HTTP=1" maayanlab/zika
			```
		4. Get the IP of your Docker machine:   
			```
			$ docker-machine ip
			```
			use `boot2docker` if you are using an earlier version of Docker initiated by `boot2docker`
			```
			$ boot2docker ip
			```
		5. Open a browser and go to http://your.docker.ip

	2. Through Graphical User Interface (GUI)
		1. Install Docker Toolbox following the instructions [here](https://www.docker.com/products/docker-toolbox) 
		2. Open [Kitematic](https://www.docker.com/products/docker-kitematic)
		3. Search for maayanlab/zika in the search box and download the Docker image
		4. Set the following variables under the 'Setting' tab:
			+ Environment Varables: 
				+ PASSWORD: your password
				+ USE_HTTP: 1
			+ Docker port: 8888
		5. Click 'START' then click the maximize button on the top right corner of 'WEB PREVIEW'

3. Deploy with your cloud provider such as [Amazon Web Services](https://www.docker.com/aws), [HP Enterprise](https://www.docker.com/aws), [IBM](https://www.docker.com/IBM), [Microsoft Azure Cloud](https://www.docker.com/microsoft) and etc.

