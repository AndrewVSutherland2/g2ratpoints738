# The C programs.  They need gcc with OpenMP, GMP, and a processor with AVX-512 (including VBMI2 and IFMA).
#   src/     the search: k3pts (points of the surface), k3inc (scores of the lines), linesearch (restriction of the sextic
#            to a line), g2rp (rational points of a curve of genus 2), cubic (the four-point closure)
#   models/  the construction of the models: conicmodp (conics modulo p), sqfilter (square test), frames (frame search of
#            the reduction), canon (canonical form of a configuration of lines)
CC = gcc
CFLAGS = -O3 -march=native -fopenmp -Wall -Wno-misleading-indentation -Wno-maybe-uninitialized
SEARCH = src/k3pts src/k3inc src/linesearch src/g2rp src/cubic
MODELS = models/conicmodp models/sqfilter models/frames models/canon

all: $(SEARCH) $(MODELS)

src/k3inc: src/k3inc.c
	$(CC) $(CFLAGS) -o $@ $< -lm
src/%: src/%.c
	$(CC) $(CFLAGS) -o $@ $< -lgmp -lm
models/conicmodp: models/conicmodp.c
	$(CC) -O3 -march=native -fopenmp -o $@ $<
models/sqfilter: models/sqfilter.c
	$(CC) -O2 -o $@ $< -lgmp
models/frames: models/frames.c
	$(CC) -O3 -march=native -o $@ $< -lgmp -lm
models/canon: models/canon.c
	$(CC) -O3 -march=native -o $@ $<

clean:
	rm -f $(SEARCH) $(MODELS)
