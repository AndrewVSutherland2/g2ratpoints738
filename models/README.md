# The construction of the models

`build_models.sh WORK` goes from Elkies' sextic (`../data/elkies163.txt`) to the list `WORK/models.txt` of the models of
X_163 with one, two or three nodes that are sums of two conics of a one-node model: 7 + 186 + 2188 sextics, each with its
singular points and the lines on which it is a square.  The programs of this directory are copied to `WORK/root/models`
and run there; nothing outside `WORK` is written.  It takes about seven minutes on 16 cores.

## The steps

1. **The lines of Elkies' sextic** (`sqlines.sh`, `sqfilter.c`, `../src/linesearch`): all lines of height at most 150 on
   which the sextic is a square.  There are 52: 43 tritangents and 9 lines through the node.
2. **A one-node model as a base** (`base.sh`).
   * Its conics: all conics modulo 41 on which the sextic is a constant times a square (`conicmodp.c`), lifted to
     characteristic 0 (`lift_conics.py`), verified exactly and parametrised (`verify_conics.gp`, `conic_table.py`).
   * Its lattice (`base.gp`): the lifts of the lines and the exceptional curve of the node span a lattice of rank 20 and
     discriminant -163, the Neron-Severi lattice.  The class of each conic is determined by its intersection numbers with
     the lifts of the lines.  All classes v with v^2 = -2 of degree 0, 1 and 2 are enumerated: those of degree 1 are
     exactly the lifts of the lines, and the irreducible ones of degree 2 are exactly the lifts of the conics, so both
     tables are complete.  (If the conics found modulo 41 did not fill the count, a second prime would be added.)
   * Its conic pairs (`pairs.py`): the classes D' = C_a + C_b for lifts of conics with C_a.C_b = 3, so that D'^2 = 2,
     and the number of roots of the lattice orthogonal to D' (they are among the classes of degree at most 2).  The
     pairs with at most three roots, pairwise orthogonal, are kept: models with one, two or three nodes.
   * Their classes (`classify.sh`, `config.gp`, `canon.c`): for each D', the contracted curves and the classes of the
     lifts of the lines of the new model are found by enumeration in the lattice, and the matrix of their intersection
     numbers is brought to a canonical form.  These classes generate the Neron-Severi lattice (the index is checked to be
     1 every time), so two polarizations have the same canonical form exactly when an isometry of the lattice takes one
     to the other, that is, when the two models are isomorphic.
3. **The other one-node models.**  Six of the classes found on Elkies' model have one node and are not its own class.
   Each is built as a sextic (step 5), moved to coordinates with its node at (0:0:1) (`node_frame.gp`), and step 2 is
   repeated on it.  No further one-node class appears on any of the seven.
4. **One pair for each class** (`choose_pairs.py`): on the first of the seven bases on which the class occurs, the pair
   with the fewest conics through the node, then of smallest height.
5. **The sextic of a pair** (`pipeline.sh`): the three sections of D' (`lifts.gp`, `sections.gp`, `interlib.gp`), the new
   branch sextic and the images of the known curves (`model.gp`, `modellib.gp`), small coordinates from a frame of three
   of the new lines (`reduce.gp`, `reducelib.gp`, `frames.c`), minimisation at every prime and reduction with Magma's
   `MinRedTernaryForm` (`minred.m`, `minred.gp`), the singular points (`sing.m`), and a final check with output
   (`finalize.gp`): the images of the known curves are lines on which the new sextic is a square and singular points of
   it, and the point counts of the new surface modulo the primes between 47 and 100 agree with the CM newform.
6. **The list** (`write_models.py`).  The number of nodes and the number of lines of every sextic are compared with the
   numbers the lattice predicts for its class.

Elkies' own sextic enters the list in minimal coordinates as well (`minred_file.sh`, `minred_file.gp`).

`../models.m` verifies the list afterwards from the sextics alone, without any of this: each is a model of X_163 with
exactly the singular points and lines listed, and the 2381 models are pairwise non-isomorphic.

## Other files

`run.sh` runs a GP script and fails if GP reports an error; `cplib.gp` holds common definitions; `base_class.sh` computes
the canonical form of the polarization of a base itself.
