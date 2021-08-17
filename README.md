<p align="center">
  <img src="viewer/static/img/decopath_logo.png" alt="DecoPath">
</p>

A web application for visualizing and exploring the results of pathway enrichment analysis.

## Table of Contents

* [Installation](#installation)
* [Issues](#issues)
* [Citation](#citation)
* [Disclaimer](#disclaimer)

## Installation

### Classical Installation

You can host DecoPath locally using Docker by following the steps below.

#### Prerequisites
In order to host DecoPath, make sure you have Docker and Docker Compose installed on your system (see [how to 
install Docker](https://docs.docker.com/engine/install/) and/or [how to install Docker Compose](https://docs.docker.com/compose/install/)).

#### Setting up DecoPath
1. Clone and download DecoPath from GitHub.
```bash
$ git clone https://github.com/decopath/decopath.git
$ cd decopath
```
2. Generate a SECRET_KEY for password encryption using the code below and add the generated text to the 
   [settings file](./DecoPath/settings.py) on line 23.
```bash
$ python3.7 generate_secret_key.py
```
3. If you are not using a domain name, remove the string *domain.com* from line 29 of the [settings](./DecoPath/settings.py)
   file and line 45 of the [docker-compose](./docker-compose.yaml) file.
   
   Additionally, you can choose to update the admin email address in the [docker-compose](./docker-compose.yaml) file on
   line 14. Note the default admin password: admin.
   

4. Build the DecoPath image using Docker.
```bash
$ docker image build -t decopath:latest .
```
5. Start the docker image using docker-compose.
```bash
$ docker-compose up -d
```
6. If you are not using a domain name, navigate to the local host (e.g., localhost:8000) from your browser and login 
with the credentials:
   - **username:** user@domain.com
   - **password:** admin

   Alternatively, if you changed the admin email address, login with your chosen username and the default password (admin).

At this point, you should be all set to host DecoPath locally.

**Optional:**
- You can obtain a domain name in case you are installing DecoPath for your entire institution
- If you wish to host DecoPath using your own domain, you might need to change the ports in
the [docker-compose](./docker-compose.yaml) file to reflect your port-forwarding configuration (See [Networking in 
Docker](https://docs.docker.com/compose/networking/)).
- You must then change "domain.com" with your domain name in the [settings](./DecoPath/settings.py) file and the 
  [docker-compose](./docker-compose.yaml) file. 

#### Deployment using a Web Server
DecoPath, out of the box, is installed in development mode. In order to deploy it, you can change the "DEBUG" variable 
in the [settings file](./DecoPath/settings.py) from "True" to "False" and serve static files along with the website 
using a web server of your choosing, such as
[nginx](https://www.cleanpython.com/how-to-use-nginx-with-your-django-docker-image/).

### For developers

DecoPath is provided as a Python 3.7 package that you can customize and edit the functionality of. You can clone the 
repository from [GitHub](https://github.com/decopath/decopath) and set up the website following the commands below:

```bash
$ git clone https://github.com/decopath/decopath.git
$ cd decopath
$ python3.7 -m pip install -r requirements.txt
$ python3.7 manage.py makemigrations
$ python3.7 manage.py migrate
$ python3.7 manage.py load_db
$ python3.7 manage.py init_user --superuser "user@domain.com"
$ python3.7 manage.py runserver
```
**Login credentials:**
   - **username:** user@domain.com
   - **password:** admin

**Note**: You can use replace "user@domain.com" with your own administrator Email ID (default admin password: admin).

In a separate terminal, RabbitMQ & Celery can be setup using:

```bash
$ docker run -d --name rabbitmq -e RABBITMQ_DEFAULT_USER=user -e RABBITMQ_DEFAULT_PASS=password -e RABBITMQ_DEFAULT_VHOST=vhost -p 8080:15672 -p 5672:5672 rabbitmq:management
$ pip install celery
$ celery -A DecoPath worker -l info
```

**Note**: To obtain fold changes, the DESeq2 package must first be installed. To do so, start R and enter:

```bash
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

#### GSEA and large gene sets

To speed up a GSEA run with large gene sets, you can clone this
[branch](https://github.com/DecoPath/DecoPath/tree/increase_processes) to enable multi-processing. Then, alter the
number of processes [here](https://github.com/DecoPath/DecoPath/blob/increase_processes/viewer/src/gsea.py).

## Issues

If you have difficulties using DecoPath, please open an issue at our bug-tracker
on [GitHub](https://github.com/DecoPath/DecoPath/issues/new).

## Citation

If you have found DecoPath useful in your work, please consider citing:

Mubeen, S., Bharadhwaj, V.S., Gadiya, Y., Hofmann-Apitius, M., Kodamullil, A.T., & Domingo-Fernández, D. (2021). 
DecoPath: A web application for decoding pathway enrichment analysis. bioRxiv. 
https://doi.org/10.1101/2021.05.22.445243

## Disclaimer

DecoPath is a scientific web application that has been developed in an academic capacity and thus comes with no
warranty or guarantee of maintenance, support, or back-up of data.
