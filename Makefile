#!/usr/bin/make -f

#
# Makefile providing commands to build docker images
# and run teh analysis in the container
#

SHELL := /bin/bash

USER_ID := $(shell id -u)
GROUP_ID := $(shell id -g)
PROJECT_ROOT := $(shell pwd)

DOCKER_NAME := correlation
# ongoing numbering of deploy instances
DOCKER_DEPLOY_NR := $(shell docker container ls -a -f name=^/$(DOCKER_NAME)_deploy$ | wc -l)
DOCKER_DEPLOY_PORT := $(shell comm -23 <(seq 49152 65535 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | shuf | head -n 1)
DOCKER_DEPLOY_NAME := $(DOCKER_NAME)_deploy_$(DOCKER_DEPLOY_NR)
DOCKER_TAG := latest
DOCKER_IMAGE := $(DOCKER_NAME):$(DOCKER_TAG)

define DOCKER_RUN_ARGS
	--volume ${PWD}:/analysis \
	--volume $(PROJECT_ROOT)/src/Rprofile.R:/usr/local/lib/R/etc/Rprofile.site \
	--hostname $(DOCKER_NAME)
endef

# Start docker container for interactive analysis
.PHONY: deploy
deploy: build
	docker run -d --rm \
		--user root:root \
		-e PASSWORD=password \
		-p $(DOCKER_DEPLOY_PORT):8787 \
		--name $(DOCKER_DEPLOY_NAME) \
		$(DOCKER_RUN_ARGS) $(DOCKER_IMAGE) \
		bash -c " \
			groupadd -g $(GROUP_ID) GROUP || true && \
			usermod -u $(USER_ID) -g $(GROUP_ID) -d /analysis rstudio && \
			env >> /etc/R/Renviron.site && \
			/init \
			" && \
	echo SUCCESS: Interactive analysis container $(DOCKER_DEPLOY_NAME) started. && \
	echo Access: http://$(shell hostname):$(DOCKER_DEPLOY_PORT)

# Run docker container. This will also run the analysis by default
.PHONY: run
run: build
	docker run --rm \
		--user $(USER_ID):$(GROUP_ID) \
		$(DOCKER_RUN_ARGS) $(DOCKER_IMAGE) \
		make analyse

# Status of DVC, git and make
.PHONY: status
status:
	@echo "\n\n----------- git status -----------"
	git status
	@echo "\n\n---------- make status -----------"
	make -n

.PHONY: build
build:
	mkdir -p log && \
	docker build $(PROJECT_ROOT) --tag $(DOCKER_IMAGE) \
		2>&1 | tee log/docker_build.log

# Main entry point for the analysis
.PHONY: analyse
analyse:
	@echo "Start anlysis at ${HOSTNAME}"
	ls results | xargs -i rm -rf results/{}
	rm -rf log/*
	set -o pipefail && \
	Rscript src/analyze.R 2>&1 | tee log/analyse.log
