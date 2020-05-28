.PHONY: build

JULIA ?= julia
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

build:
	cd $(realpath $(mkfile_path)..) && \
	$(JULIA) --startup-file=no --project=docs/ \
		-e 'using Pkg; Pkg.develop(PackageSpec(path=pwd()));' && \
	$(JULIA) --startup-file=no --project=docs/ \
		docs/make.jl
