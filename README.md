# A genus 2 curve with 738 rational points: the computation

This repository holds the programs, the input and the output of the computation described in the note
*A genus 2 curve with 738 rational points* (A. V. Sutherland, November 2026).  The curve

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

    make                 # the C programs (gcc with OpenMP, GMP; a processor with AVX-512)
    ./run_all.sh WORK    # models, search, tables; compares the result with data/

`run_all.sh` writes only into the directory `WORK` (about 1 GB).  On a machine with a 16-core AMD Ryzen AI Max+ 395
processor (32 threads) the whole run takes 62 minutes, 7 of them for the models; `data/search_summary.txt`
has the times of the run stored here, part by part.  Nothing in the computation is random: the output is determined by the input and the settings, and the last
line of `run_all.sh` says whether the models, the numbers of the search and the curves are those of `data/`.

It needs, besides gcc and GMP: PARI/GP (2.18.1 was used), Magma (V2.29-10 was used; its `MinRedTernaryForm` computes the
minimal coordinates, and it finds the singular points of each sextic), Python 3 with numpy, and the usual shell tools
(bash, awk, xargs, bc).  The sextic of a model depends on the coordinates that `MinRedTernaryForm` returns, so another
version of Magma could give other coordinates for some models, and then other lines of height at most 400.

Parts of the run by themselves:

    models/build_models.sh WORK          # the models only: WORK/models.txt (7 minutes)
    search/search_model.sh t0453         # the search on one model of data/models.txt: the model of C_738 (20 seconds)
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

Not here: the random walks of the lattice census that the note quotes for the number of models of each kind (the
computation above does not use the census; `data/census.txt` holds the classes it found with at most three nodes), and the
enumeration of points of bounded canonical height that the note gives as evidence that the list of 738 points is complete.
