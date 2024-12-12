SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools 
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y --override-channels -c conda-forge -c bioconda samtools bcftools minimap2 NanoPlot svim-asm && touch .tools

clean:
	rm -rf mamba/ .mamba .tools
