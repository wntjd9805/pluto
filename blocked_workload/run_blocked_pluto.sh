#!/bin/bash
set -euo pipefail

# Override with env if needed: MARCH=sapphirerapids ./run_blocked_pluto.sh
: "${MARCH:=znver4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="${BENCH_ROOT:-$SCRIPT_DIR}"
POLYBENCH_DIR="${POLYBENCH_DIR:-$BENCH_ROOT/utilities}"

OMP_THREADS="${OMP_THREADS:-48}"
OMP_BIND="${OMP_BIND:-CLOSE}"
RUNS_PER_ROUND="${RUNS_PER_ROUND:-11}"
ROUNDS="${ROUNDS:-1}"

LOG_ROOT="${LOG_ROOT:-$BENCH_ROOT/pluto_experiment_logs/$(date +%Y%m%d_%H%M%S)}"
SUMMARY_CSV="$LOG_ROOT/summary.csv"

declare -a FILES=(
  # "batched_gemm/small/batched_gemm_small_pluto.c"
  # "batched_gemm/standard/batched_gemm_standard_pluto.c"
  # "batched_gemm/large/batched_gemm_large_pluto.c"
  # "conv_im2col/small/conv_im2col_small_pluto.c"
  # "conv_im2col/standard/conv_im2col_standard_pluto.c"
  # "conv_im2col/large/conv_im2col_large_pluto.c"
  # "geqrf/small/geqrf_small_pluto.c"
  # "geqrf/standard/geqrf_standard_pluto.c"
  # "syrk/small/syrk_small_pluto.c"
  # "syrk/standard/syrk_standard_pluto.c"
  # "syrk/large/syrk_large_pluto.c"
  # "trmm_block/small/trmm_blocked_pluto.c"
  # "trmm_block/standard/trmm_blocked_pluto.c"
  # "trmm_block/large/trmm_blocked_pluto.c"
  # "geqrf/large/geqrf_large_pluto.c"
  "ruiz_equilibration/small/ruiz_equilibration_small_pluto.c"
  "ruiz_equilibration/standard/ruiz_equilibration_standard_pluto.c"
  "ruiz_equilibration/large/ruiz_equilibration_large_pluto.c"
  "Multigrid/small/MG_small_pluto.c"
  "Multigrid/standard/MG_standard_pluto.c"
  "Multigrid/large/MG_large_pluto.c"
)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "FILES array is empty."
  exit 1
fi
if [[ ! -d "$BENCH_ROOT" ]]; then
  echo "BENCH_ROOT not found: $BENCH_ROOT"
  exit 1
fi
if [[ ! -f "$POLYBENCH_DIR/polybench.c" ]]; then
  echo "polybench.c not found: $POLYBENCH_DIR/polybench.c"
  exit 1
fi

mkdir -p "$LOG_ROOT"
echo "benchmark,size,src,round1_median,round2_median,round3_median,final_median,status" > "$SUMMARY_CSV"

echo "LOG_ROOT=$LOG_ROOT"
echo "BENCH_ROOT=$BENCH_ROOT"
echo "POLYBENCH_DIR=$POLYBENCH_DIR"
echo "MARCH=$MARCH"
echo

median_of_numbers() {
  if (($# == 0)); then
    echo "NA"
    return
  fi

  local sorted=()
  mapfile -t sorted < <(printf '%s\n' "$@" | sort -n)
  local n=${#sorted[@]}

  if ((n % 2 == 1)); then
    echo "${sorted[$((n / 2))]}"
  else
    local a="${sorted[$((n / 2 - 1))]}"
    local b="${sorted[$((n / 2))]}"
    echo "scale=6; ($a + $b) / 2" | bc
  fi
}

extract_time() {
  echo "$1" | grep -E '^[0-9]+\.[0-9]+' | head -n 1 || true
}

MEASURE_R1="NA"
MEASURE_R2="NA"
MEASURE_R3="NA"
MEASURE_FINAL="NA"

measure_exec() {
  local exec_name=$1
  local log_file=$2
  local round_medians=()

  : > "$log_file"

  for round in $(seq 1 "$ROUNDS"); do
    local times=()
    echo "===== ROUND $round =====" >> "$log_file"

    for i in $(seq 1 "$RUNS_PER_ROUND"); do
      echo "[$exec_name] Run #$i" >> "$log_file"
      output=$(OMP_NUM_THREADS=$OMP_THREADS OMP_PROC_BIND=$OMP_BIND ./"$exec_name" 2>&1 || true)
      echo "$output" >> "$log_file"
      echo >> "$log_file"

      t=$(extract_time "$output")
      if [[ -n "$t" ]]; then
        times+=("$t")
      fi
    done

    if [[ ${#times[@]} -eq 0 ]]; then
      round_medians+=("NA")
      echo "ROUND_${round}_MEDIAN=NA" >> "$log_file"
    else
      m=$(median_of_numbers "${times[@]}")
      round_medians+=("$m")
      echo "ROUND_${round}_MEDIAN=$m" >> "$log_file"
    fi
    echo >> "$log_file"
  done

  MEASURE_R1="${round_medians[0]:-NA}"
  MEASURE_R2="${round_medians[1]:-NA}"
  MEASURE_R3="${round_medians[2]:-NA}"

  local valid=()
  [[ "$MEASURE_R1" != "NA" ]] && valid+=("$MEASURE_R1")
  [[ "$MEASURE_R2" != "NA" ]] && valid+=("$MEASURE_R2")
  [[ "$MEASURE_R3" != "NA" ]] && valid+=("$MEASURE_R3")
  MEASURE_FINAL=$(median_of_numbers "${valid[@]}")
}

for FILE in "${FILES[@]}"; do
  REL_DIR=$(dirname "$FILE")
  BASE=$(basename "$FILE")
  SRC=${BASE%.c}
  SIZE=$(basename "$REL_DIR")
  BENCH=$(basename "$(dirname "$REL_DIR")")
  WORK_DIR="$BENCH_ROOT/$REL_DIR"
  TAG="$BENCH/$SIZE"

  echo "============================================================"
  echo "[$TAG] src=$SRC"

  if [[ ! -d "$WORK_DIR" ]]; then
    echo "[$TAG] skip: directory not found ($WORK_DIR)"
    echo "$BENCH,$SIZE,$SRC,NA,NA,NA,NA,missing_dir" >> "$SUMMARY_CSV"
    continue
  fi
  if [[ ! -f "$WORK_DIR/$BASE" ]]; then
    echo "[$TAG] skip: source not found ($WORK_DIR/$BASE)"
    echo "$BENCH,$SIZE,$SRC,NA,NA,NA,NA,missing_src" >> "$SUMMARY_CSV"
    continue
  fi

  BUILD_LOG_DIR="$LOG_ROOT/${BENCH}__${SIZE}"
  mkdir -p "$BUILD_LOG_DIR"

  BUILD_LOG="$BUILD_LOG_DIR/build_make.log"
  if ! (cd "$WORK_DIR" && make "$SRC.par.c") > "$BUILD_LOG" 2>&1; then
    echo "[$TAG] FAIL: make $SRC.par.c"
    echo "$BENCH,$SIZE,$SRC,NA,NA,NA,NA,make_failed" >> "$SUMMARY_CSV"
    continue
  fi

  if [[ ! -f "$WORK_DIR/$SRC.par.c" ]]; then
    echo "[$TAG] FAIL: $SRC.par.c not found"
    echo "$BENCH,$SIZE,$SRC,NA,NA,NA,NA,par_missing" >> "$SUMMARY_CSV"
    continue
  fi

  EXE="omp_exec_${BENCH}_${SIZE}"
  COMPILE_LOG="$BUILD_LOG_DIR/compile.log"
  SIZE_OPT=""
  if [[ "$BENCH" == "trmm_block" ]]; then
    case "$SIZE" in
      small) SIZE_OPT="-DSMALL_DATASET" ;;
      standard) SIZE_OPT="-DSTANDARD_DATASET" ;;
      large) SIZE_OPT="-DLARGE_DATASET" ;;
      *)
        echo "[$TAG] FAIL: unknown TRMM size '$SIZE'"
        echo "$BENCH,$SIZE,$SRC,NA,NA,NA,NA,unknown_size" >> "$SUMMARY_CSV"
        continue
        ;;
    esac
  fi

  if ! (
    cd "$WORK_DIR"
    echo "Compiling $EXE ${SIZE_OPT:+($SIZE_OPT)}"
    clang -O3 -march="$MARCH" -fopenmp -DPOLYBENCH_TIME \
      ${SIZE_OPT:+$SIZE_OPT} \
      -I"$POLYBENCH_DIR" \
      -I"$BENCH_ROOT" \
      "$SRC.par.c" \
      "$POLYBENCH_DIR/polybench.c" \
      -o "$EXE" -lstdc++ -lm
  ) > "$COMPILE_LOG" 2>&1; then
    echo "[$TAG] FAIL: compile"
    echo "$BENCH,$SIZE,$SRC,NA,NA,NA,NA,compile_failed" >> "$SUMMARY_CSV"
    continue
  fi

  RUN_LOG="$BUILD_LOG_DIR/run.log"
  pushd "$WORK_DIR" >/dev/null
  measure_exec "$EXE" "$RUN_LOG"
  popd >/dev/null

  echo "$BENCH,$SIZE,$SRC,$MEASURE_R1,$MEASURE_R2,$MEASURE_R3,$MEASURE_FINAL,ok" >> "$SUMMARY_CSV"
  echo "[$TAG] medians: $MEASURE_R1, $MEASURE_R2, $MEASURE_R3 -> final=$MEASURE_FINAL"
  echo
done

echo "Done."
echo "Summary: $SUMMARY_CSV"
echo "Logs: $LOG_ROOT"
