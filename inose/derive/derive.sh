#!/bin/bash
# derive.sh : recompute map_data.m and inverse_data.m in a work directory (about two minutes on 16 cores).
#   usage: derive.sh WORKDIR        (WORKDIR must not exist)
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd); W=$1; mkdir "$W"; cd "$W"
cp "$D"/elkdata.m "$D"/pencil.m "$D"/t_q.m "$D"/u_q.m "$D"/fit.m "$D"/mkdata.m "$D"/invlib.m "$D"/inv_prime.m "$D"/crt.py "$D"/lll.m "$D"/mkinv.m .
magma -b t_q.m < /dev/null > t_q.log 2>&1 &            # the pencil of the Inose fibration: forms of degree 32
magma -b u_q.m < /dev/null > u_q.log 2>&1 &            # the second pencil: forms of degree 24
wait
magma -b fit.m < /dev/null > fit.log 2>&1              # the two functions; norm identities
magma -b mkdata.m < /dev/null > mkdata.log 2>&1        # map_data.m
seq 1 32 | xargs -P 16 -I{} sh -c 'magma -b I:={} inv_prime.m < /dev/null > log_p{}.txt 2>&1'   # the inverse modulo 32 primes
python3 crt.py 30 > crt.log                            # rational reconstruction from 30 primes, checked on 2 more
magma -b lll.m < /dev/null > lll.log 2>&1              # small coefficients
magma -b mkinv.m < /dev/null > mkinv.log 2>&1          # inverse_data.m
echo "derive.sh: done; compare map_data.m and inverse_data.m with the stored files"
