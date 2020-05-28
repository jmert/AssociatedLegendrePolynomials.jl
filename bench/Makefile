.PHONY: benchmark

JULIA ?= julia
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

benchmark:
	cd $(realpath $(mkfile_path)..) && \
	$(JULIA) --startup-file=no --project=bench/ \
		-e 'using Pkg; Pkg.develop(PackageSpec(path=pwd()));' && \
	$(JULIA) --startup-file=no --project=bench/ \
		bench/benchmark.jl -- --save
