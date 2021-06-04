<p align="center">
  <img src="viewer/static/img/decopath_logo.png">
</p>

A web application for visualizing and exploring the results of pathway enrichment analysis.

## Table of Contents

* [General Info](#general-info)
* [Installation](#installation)
* [Documentation](#documentation)
* [Usage](#usage)
* [Issues](#issues)
* [Citation](#citation)
* [Disclaimer](#disclaimer)

## Installation

### For a classical installation

You can host DecoPath using Docker. To do so, run the following commands,

```bash
$ docker image build -t decopath:latest .
$ docker-compose up
```

**Note**: If you wish to host DecoPath locally, you might need to change the ports in
the [docker-compose](./docker-compose.yaml)
file to reflect your port-forwarding configuration.

### For developers

You can clone the repository from [GitHub](https://github.com/decopath/decopath).

```bash
$ git clone https://github.com/decopath/decopath.git
$ cd decopath
$ python -m pip install -r requirements.txt
$ python manage.py makemigrations
$ python manage.py migrate
$ python manage.py load_db
$ python manage.py runserver
```

And RabbitMQ & Celery can be setup using,

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

## GSEA and large gene sets

To speed up a GSEA run with large gene sets, you can clone this
[branch](https://github.com/DecoPath/DecoPath/tree/increase_processes) to enable multi-processing. Then, alter the
number of processes [here](https://github.com/DecoPath/DecoPath/blob/increase_processes/viewer/src/gsea.py).

## Issues

If you have difficulties using DecoPath, please open an issue at our bug-tracker
on [GitHub](https://github.com/DecoPath/DecoPath/issues/new).

## Citation

If you have found DecoPath useful in your work, please consider citing:

Mubeen, S., Bharadhwaj, V. S., Gadiya, Y., Hofmann-Apitius, M., & Domingo-Fernández, D. (2021). DecoPath: A web
application for decoding pathway enrichment analysis. bioRxiv. https://doi.org/10.1101/2021.05.22.445243



## Disclaimer

DecoPath is a scientific web application that has been developed in an academic capacity, and thus comes with no
warranty or guarantee of maintenance, support, or back-up of data.
