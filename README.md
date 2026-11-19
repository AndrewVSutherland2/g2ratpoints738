# A genus 2 curve with 738 rational points: the computation

This repository holds the programs, the input and the output of the computation described in the note
*A genus 2 curve with 738 rational points* (A. V. Sutherland, November 2026), which is at
<https://math.mit.edu/~drew/g2ratpoints738.html> and, as a paper, in `curve738.pdf` here.  The curve

    C_738:  y^2 = 743425 x^6 - 16144220 x^5 + 60026500 x^4 + 172162750 x^3 + 40133975 x^2 - 46160230 x + 6838225

is a line section of a double-plane model w^2 = S(z0, z1, z2) of the K3 surface X_163 of Picard number 20 and discriminant
-163.  The computation starts from one sextic, Elkies' model of that surface, and nothing else:

1. **Models.**  From Elkies' sextic it constructs a sextic, in minimal coordinates, for every model of the surface with one,
   two or three nodes that is the sum of two conics of a one-node model: 7 + 186 + 2188 = 2381 models.
2. **Search.**  On each model it lists the rational points of the surface of height at most 4000, counts them on every line
   of height at most 400 (the *score* of the line), counts rational points on the curves cut out by the lines with a
   score of at least 40, and closes the sets of points of the best curves under the relation given by a cubic through
   four points.

It finds 1061 curves with at least 440 points of height at most 10^5, sixteen of them with at least 600 points after
closure: C_738, a second new curve with 654 points, and the curve with 642 points of Elkies and Stoll among them.

## Running it

    make                                 # the C programs
    search/search_model.sh t0453         # a first check: the search on the model of C_738 (output in out_t0453/)
    ./run_all.sh WORK                    # everything: the models, the search, the tables; compares the result with data/
    ./run_all.sh WORK data/models.txt    # the search only, on the stored models (no Magma needed)

`search/search_model.sh t0453` takes about 15 seconds; `out_t0453/summary.txt` should say
`curves closed: 1; best 738 points (662 of height <= 10^5) on the line (66, 31, -52)`.  `make check` is a longer test
(a few minutes): it searches three models, one of each kind, and compares the numbers and the curves with the stored ones.

`run_all.sh` writes only into the directory `WORK` (about 1 GB).  Nothing in the computation is random: the output is
determined by the input, the settings and the versions of the programs named below, and the last line that `run_all.sh`
prints says whether the models, the numbers of the search and the curves are those of `data/`:

    the models, the numbers of the search and the curves closed are those of .../data

On a machine with a 16-core AMD Ryzen AI Max+ 395 processor (32 threads) the whole run takes 62 minutes, 7 of them for
the models; `data/search_summary.txt` has the times of the run stored here, part by part.  The search runs
`nproc`/2 models at a time with 2 threads each (`JOBS` and `P` in the environment change that), so the time is roughly
inversely proportional to the number of cores.  The search of one model needs about 0.6 GB of memory.

### What it needs

- A 64-bit Linux system with bash, GNU make, awk, xargs, bc, flock and the GNU core utilities.
- gcc with OpenMP (13.3 was used) and the GMP library with its header (`libgmp-dev`).
- PARI/GP 2.17 or later (2.18.1 was used; 2.17.2 gives the same results, 2.15 does not: see below).
- Python 3; numpy for the construction of the models.
- Magma for the construction of the models (V2.29-10 was used): its `MinRedTernaryForm` computes the minimal
  coordinates, and it finds the singular points of each sextic.  The search does not use Magma.  The certificates
  (below) are Magma scripts, apart from those for the 738 points.

### What the result depends on

- **The processor.**  Four of the five programs of the search work on 512-bit words.  They use the AVX-512 instructions
  of the processor when it has them (the subsets F, BW, VL, VBMI2 and IFMA: AMD Zen 4 and later, Intel Ice Lake and
  later server and some client processors), and plain C otherwise (`src/simd.h`).  The output is the same: the whole
  search was repeated with the plain C programs (`make MARCH=x86-64-v3`), and the output files of each of the 2381
  models are those of the stored run, apart from the running times that they record.  On the machine above the plain C
  programs take about twice as long (26 seconds instead of 12.5 on the model of C_738).
- **PARI/GP.**  The search counts the points of bounded height on a reduced equation of each curve, computed by
  `hyperellred`.  Since PARI 2.17 that function reduces a binary quadratic form exactly; before, it used a
  floating-point LLL.  On the 6,605,562 curves of the stored run, 2.17.2 returns the same equation as 2.18.1 for every
  curve.  2.15.4 (the version packaged in Ubuntu 24.04) returns another equation for 2,288,212 of them; for all but
  39,671 (0.6%) the two equations differ by x -> -x or x -> 1/x, which does not change the counts, and for those 39,671
  some of the counts differ.  The construction of the models depends on the version of PARI/GP as well.  With 2.17.2 it
  gives the same 2381 sextics as with 2.18.1.  With 2.15.4 everything runs, but 78 of the 2381 models come out in other
  coordinates, so the results cannot be compared with `data/` row by row; `run_all.sh` prints a warning.
- **Magma.**  The sextic of a model depends on the coordinates that `MinRedTernaryForm` returns, so another version of
  Magma could give other coordinates for some models, and then other lines of height at most 400.  This has not been
  tried: only V2.29-10 was available.  `run_all.sh WORK data/models.txt` searches the stored models and does not depend
  on it.

### Parts of the run by themselves

    models/build_models.sh WORK          # the models only: WORK/models.txt (7 minutes)
    search/search_model.sh NAME          # the search on one model of data/models.txt
    search/line.sh t0453 66 31 -52       # one line of one model: its curve, the points of height <= 10^5, the closure
    search/sample.sh WORK2               # 87 models again with lower cuts: what the selection loses (5 processor hours)

## What is where

| | |
|---|---|
| `data/elkies163.txt` | the input: the 28 coefficients of Elkies' sextic |
| `inose/` | the surface: explicit birational maps, in both directions, between Elkies' sextic and Inose's equation y^2 = x^3 - 3a^2 x - (t + 1/t + 2a^3 + 2), a = 53360, with proofs (`inose/README.md`) |
| `models/` | the construction of the models (`models/README.md`) |
| `search/` | the search: `search_model.sh` (one model), `search_all.sh` (all), `reduce.gp` (reduced equations of the curves), `collect.py` (tables), `time_by_part.py`, `line.sh` (one line), `model_lines.py` (all lines kept on one model), `sample.sh` and `sample_table.py` (the sample with lower cuts) |
| `src/` | the five C programs of the search: `k3pts` (points of the surface), `k3inc` (scores), `linesearch` (restriction of the sextic to a line), `g2rp` (rational points of a curve of genus 2), `cubic` (the closure) |
| `data/models.txt` | output of step 1: the 2381 sextics with their singular points and the lines on which they are squares |
| `data/elkies163_lines.txt`, `data/elkies163_conics.txt`, `data/minimisation.txt` | also from step 1: the 52 lines and the 1241 conics on which Elkies' sextic is a square, and what the minimisation changed for each sextic |
| `data/search_models.txt` | output of step 2, one row per model: the numbers at each stage and the times |
| `data/search_curves.txt` | the 1061 curves closed: points, model, line, sextic, score, Galois group, isomorphism class |
| `data/S3_search.txt` | the scores and counts of all the lines kept on the model of C_738 |
| `data/search_sample.txt` | a sample of 87 models searched again with lower cuts (score at least 30, thresholds 200 and 300), to measure what the selection loses (`search/sample.sh`) |
| `data/search_summary.txt` | machine, programs, settings and times of the stored run |
| `curve738.pdf`, `curve738.tex`, `figs/` | the note as a paper, and its source (which also reads `k3model.m`) |

## Certificates

Each script says in its first lines what it proves and how to run it; the `.log` files are their outputs.

| statement | files |
|---|---|
| C_738 has at least 738 rational points | `closure.js` regenerates them from five x-values (it is the program of the HTML version of the note; `test_closure.py` runs it in headless Firefox); `data/witnesses_738.json` holds one relation for each of the other 364 x-values, verified without any search by `check_witnesses.py`; the list is `data/points_738.txt` |
| Elkies' sextic is a model of the surface defined by Inose's equation | `inose/verify.m`, `inose/verify_inverse.m` (Magma) |
| each of the 2381 sextics is a model of X_163 with exactly the nodes and lines listed, and they are 2381 different models | `models.m` (Magma, about ten minutes) |
| the 2381 models are classes of the lattice census quoted in the note: all those with one or two nodes, and 2188 of the 2192 with three | `census.m` (Magma, after `magma -b gramout:=grams.txt models.m`), `data/census.txt` |
| the model of C_738: three nodes, 39 lines, 810 conics, point counts | `k3model.m` (Magma, seconds), `data/S3_coefficients.txt`, `data/S3_lines.txt`, `data/S3_conics_810.txt`; `data/model_map.txt` is the map from Elkies' model |
| C_738 is a line section of exactly one of the 2381 sextics, by exactly one line, over an algebraic closure of Q | `linesection.m` (Magma, a few hours in all; it can be run in pieces) |
| each curve of `data/search_curves.txt` is the stated line section; Galois groups; isomorphism classes | `curves.m` (Magma, five minutes) |
| the Jacobian of C_738: no torsion, rank at least 21, a basis of the group generated by the known points, saturation at 2, ..., 13, regulator | `mw.m` (Magma, 20 seconds), `data/mw_basis.txt`, `data/mw_coords.txt` |
| rank at most 21, assuming GRH | `selmer.m` (Magma, about four hours) |
| the masses of the two genera quoted in the note | `mass.m` (Magma) |
| the formulas displayed in the note: the equations are those of `data/`; the map from Elkies' model to the model of C_738 satisfies the identity stated; the line of C_738; the discriminant, the reductions and the numbers of points by height quoted for it; the example of the closure | `verify.py` (Python and PARI/GP, seconds), which takes the formulas from the source of the note, `curve738.tex` |
| the figures that the note quotes for the computation are those of the files in `data/` and of the logs above | `check_census.py` (Python, seconds).  It runs 29 checks here; 14 more compare the note with the files of the lattice census, which are not here, and `check_census.log` is the output of all 43 |

Not here: the random walks of the lattice census that the note quotes for the number of models of each kind (the
computation above does not use the census; `data/census.txt` holds the classes it found with at most three nodes), and the
enumeration of points of bounded canonical height that the note gives as evidence that the list of 738 points is complete.

## Licence

The programs and the data are under the MIT licence (`LICENSE`).  The note itself (`curve738.pdf`, `curve738.tex` and the
figures in `figs/`) is not covered by it: copyright 2026 Andrew V. Sutherland.
