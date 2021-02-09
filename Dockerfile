# matching rstan -> banocc
FROM rocker/verse:4.0.2 
WORKDIR /build

ENV DEBIAN_FRONTEND noninteractive
ENV TZ Europe/Berlin

# install debian packages
RUN apt-get update && \
	apt-get install -y git parallel && \
	rm -rf /var/lib/apt/lists/*

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.2-Linux-x86_64.sh && \
	bash Miniconda3-py38_4.8.2-Linux-x86_64.sh -b -p /miniconda3 && \
 	/miniconda3/bin/conda init && \
	rm Miniconda*.sh && \
	chmod a+rwx -R /miniconda3
ENV PATH="/miniconda3/bin:${PATH}"

# install conda packages
RUN conda init bash &&  \
	conda create -y -n fastspar -c bioconda -c conda-forge fastspar

# install R packages
COPY src/install.R .
RUN Rscript install.R

# setup
USER rstudio
RUN conda init bash
EXPOSE 8787
SHELL ["/bin/bash"]
CMD ["/init"]
