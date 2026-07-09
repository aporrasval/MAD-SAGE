#!/bin/bash
# Build darksage_coherent: same as darksage but with coherentonlylcode/model_mad.c (coherent=1 always)
set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"
GSLINC="-I/usr/local/Cellar/gsl/2.6/include"
GSLLIB="-L/usr/local/Cellar/gsl/2.6/lib -lgsl -lgslcblas -Xlinker -rpath -Xlinker /usr/local/Cellar/gsl/2.6/lib"
CFLAGS="-g -O2 -Wall -Wextra -Wshadow -I${ROOT}/code ${GSLINC}"

SRCS=(
  code/main.c
  code/core_read_parameter_file.c
  code/core_init.c
  code/core_io_tree.c
  code/core_cool_func.c
  code/core_build_model.c
  code/core_save.c
  code/core_mymalloc.c
  code/core_allvars.c
  code/model_infall.c
  code/model_cooling.c
  code/model_starformation_and_feedback.c
  code/model_disk_instability.c
  code/model_reincorporation.c
  code/model_mergers.c
  code/model_misc.c
)

OBJS=()
cd "$ROOT"

echo "Compiling shared object files..."
for f in "${SRCS[@]}"; do
  obj="${f%.c}_co.o"
  clang $CFLAGS -c "$f" -o "$obj"
  OBJS+=("$obj")
done

echo "Compiling coherent-only model_mad.c..."
clang $CFLAGS -c coherentonlylcode/model_mad.c -o coherentonlylcode/model_mad_co.o
OBJS+=("coherentonlylcode/model_mad_co.o")

echo "Linking darksage_coherent..."
clang "${OBJS[@]}" -lm $GSLLIB -o darksage_coherent

echo "Cleaning object files..."
rm -f "${OBJS[@]}"

echo "Done: darksage_coherent"
