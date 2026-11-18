# The C programs.  They need gcc with OpenMP, and GMP.
#   src/     the search: k3pts (points of the surface), k3inc (scores of the lines), linesearch (restriction of the sextic
#            to a line), g2rp (rational points of a curve of genus 2), cubic (the four-point closure)
#   models/  the construction of the models: conicmodp (conics modulo p), sqfilter (square test), frames (frame search of
#            the reduction), canon (canonical form of a configuration of lines)
# The four programs of src/ that work on 512-bit words use the AVX-512 instructions of the processor when it has them
# (F, BW, VL, VBMI2, IFMA) and plain C otherwise (src/simd.h); the results are the same.  "make MARCH=x86-64-v3" builds
# the plain C versions on a machine that has AVX-512.
CC = gcc
MARCH ?= native
CFLAGS = -O3 -march=$(MARCH) -fopenmp -Wall -Wno-psabi -Wno-misleading-indentation -Wno-maybe-uninitialized
SEARCH = src/k3pts src/k3inc src/linesearch src/g2rp src/cubic
MODELS = models/conicmodp models/sqfilter models/frames models/canon

all: $(SEARCH) $(MODELS)

src/k3inc: src/k3inc.c src/simd.h
	$(CC) $(CFLAGS) -o $@ $< -lm
src/%: src/%.c src/simd.h
	$(CC) $(CFLAGS) -o $@ $< -lgmp -lm
models/conicmodp: models/conicmodp.c
	$(CC) -O3 -march=$(MARCH) -fopenmp -o $@ $<
models/sqfilter: models/sqfilter.c
	$(CC) -O2 -o $@ $< -lgmp
models/frames: models/frames.c
	$(CC) -O3 -march=$(MARCH) -o $@ $< -lgmp -lm
models/canon: models/canon.c
	$(CC) -O3 -march=$(MARCH) -o $@ $<

# a short test: the search on three models, compared with the stored results (see check.sh)
check: all
	./check.sh

clean:
	rm -f $(SEARCH) $(MODELS)

.PHONY: all check clean
