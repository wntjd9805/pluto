#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

LLVM_PREFIX="${LLVM_PREFIX:-/usr/local}"
CLANG_PREFIX="${CLANG_PREFIX:-}"
HYPERF_HOME="${HYPERF_HOME:-/root/hitune}"
FILECHECK_BIN="${FILECHECK_BIN:-}"
INSTALL_PREFIX="${INSTALL_PREFIX:-/usr/local}"
RUN_TESTS=1
DO_INSTALL=1
PERSIST_ENV=1
INSTALL_DEPS=1

usage() {
  cat <<'EOF'
Usage: ./scripts/install.sh [options]

Options:
  --llvm-prefix <path>   LLVM prefix for llvm-config/libs/bin (default: /usr/local)
  --clang-prefix <path>  prefix passed to --with-clang-prefix (default: /usr/local)
  --filecheck-bin <path> PATH entry that contains FileCheck binary
  --prefix <path>        install prefix for make install (default: /usr/local)
  --skip-deps            skip apt dependency installation
  --skip-tests           skip make check-pluto
  --no-install           skip make install
  --no-persist-env       do not append LLVM exports to ~/.bashrc
  -h, --help             show this help

Env vars:
  LLVM_PREFIX            same as --llvm-prefix
  CLANG_PREFIX           same as --clang-prefix
  HYPERF_HOME            used for default FileCheck path (default: /root/hitune)
  FILECHECK_BIN          same as --filecheck-bin
  INSTALL_PREFIX         same as --prefix
  INSTALL_DEPS           set 0 to skip apt dependency installation
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --llvm-prefix)
      LLVM_PREFIX="${2:-}"
      shift 2
      ;;
    --clang-prefix)
      CLANG_PREFIX="${2:-}"
      shift 2
      ;;
    --filecheck-bin)
      FILECHECK_BIN="${2:-}"
      shift 2
      ;;
    --prefix)
      INSTALL_PREFIX="${2:-}"
      shift 2
      ;;
    --skip-deps)
      INSTALL_DEPS=0
      shift
      ;;
    --skip-tests)
      RUN_TESTS=0
      shift
      ;;
    --no-install)
      DO_INSTALL=0
      shift
      ;;
    --no-persist-env)
      PERSIST_ENV=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

install_deps_if_needed() {
  if [[ "${INSTALL_DEPS}" -eq 0 ]]; then
    echo "[pluto-install] skipping dependency install (--skip-deps)"
    return
  fi

  if ! command -v apt-get >/dev/null 2>&1; then
    echo "[pluto-install] apt-get not found, skipping dependency install"
    return
  fi

  local deps=(libtool libtool-bin)
  local missing=()
  local pkg
  for pkg in "${deps[@]}"; do
    if command -v dpkg >/dev/null 2>&1; then
      dpkg -s "${pkg}" >/dev/null 2>&1 || missing+=("${pkg}")
    else
      missing+=("${pkg}")
    fi
  done

  if [[ "${#missing[@]}" -eq 0 ]]; then
    echo "[pluto-install] apt deps already installed: ${deps[*]}"
    return
  fi

  echo "[pluto-install] installing apt deps: ${missing[*]}"
  if [[ "${EUID}" -eq 0 ]]; then
    apt-get install -y "${missing[@]}"
  elif command -v sudo >/dev/null 2>&1; then
    sudo apt-get install -y "${missing[@]}"
  else
    echo "Need root/sudo to install packages: ${missing[*]}" >&2
    exit 1
  fi
}

install_deps_if_needed

if [[ -z "${FILECHECK_BIN}" && -d "${HYPERF_HOME}/llvm-project/build/bin" ]]; then
  FILECHECK_BIN="${HYPERF_HOME}/llvm-project/build/bin"
fi

if [[ -z "${CLANG_PREFIX}" ]]; then
  CLANG_PREFIX="/usr/local"
fi

if [[ ! -d "${CLANG_PREFIX}" ]]; then
  echo "Invalid clang prefix: ${CLANG_PREFIX} (directory does not exist)" >&2
  exit 1
fi

if [[ ! -x "${LLVM_PREFIX}/bin/llvm-config" ]]; then
  echo "Invalid LLVM prefix: ${LLVM_PREFIX} (missing ${LLVM_PREFIX}/bin/llvm-config)" >&2
  exit 1
fi

if [[ -n "${FILECHECK_BIN}" ]]; then
  export PATH="${FILECHECK_BIN}:${PATH}"
fi

export HYPERF_HOME
export LLVM_PREFIX
export PATH="${LLVM_PREFIX}/bin:${PATH}"
if [[ -d "${LLVM_PREFIX}/lib64" ]]; then
  export LD_LIBRARY_PATH="${LLVM_PREFIX}/lib64:${LD_LIBRARY_PATH:-}"
fi
if [[ -d "${LLVM_PREFIX}/lib" ]]; then
  export LD_LIBRARY_PATH="${LLVM_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
fi

if ! command -v FileCheck >/dev/null 2>&1; then
  cat <<'EOF' >&2
FileCheck not found in PATH.
Use one of:
  1) export PATH=$HYPERF_HOME/llvm-project/build/bin:$PATH
  2) ./scripts/install.sh --filecheck-bin /path/to/llvm-project/build/bin
  3) sudo ln -s /path/to/llvm-project/build/bin/FileCheck /usr/local/bin/FileCheck
EOF
  exit 1
fi

echo "[pluto-install] CLANG_PREFIX=${CLANG_PREFIX}"
echo "[pluto-install] HYPERF_HOME=${HYPERF_HOME}"
echo "[pluto-install] FILECHECK_BIN=${FILECHECK_BIN:-<none>}"
echo "[pluto-install] LLVM_PREFIX=${LLVM_PREFIX}"
echo "[pluto-install] INSTALL_PREFIX=${INSTALL_PREFIX}"

if [[ ! -x ./configure ]]; then
  echo "[pluto-install] configure not found, running ./autogen.sh"
  ./autogen.sh
fi

if [[ -f .gitmodules ]]; then
  echo "[pluto-install] syncing submodules"
  git submodule update --init --recursive
fi

echo "[pluto-install] running configure"
configure_args=(
  "--prefix=${INSTALL_PREFIX}"
  "--with-clang-prefix=${CLANG_PREFIX}"
)
./configure "${configure_args[@]}"

echo "[pluto-install] normalizing autotools timestamps"
find . -type f \( -name 'Makefile.in' -o -name 'configure' -o -name 'aclocal.m4' \) -exec touch {} +

echo "[pluto-install] building"
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"

if [[ "${RUN_TESTS}" -eq 1 ]]; then
  echo "[pluto-install] running tests (check-pluto)"
  make check-pluto
fi

if [[ "${DO_INSTALL}" -eq 1 ]]; then
  echo "[pluto-install] installing"
  if ! make install; then
    cat <<EOF >&2
make install failed (likely permission issue for prefix: ${INSTALL_PREFIX}).
Try one of:
  1) sudo make install
  2) ./scripts/install.sh --prefix \$HOME/.local
EOF
    exit 1
  fi
fi

if [[ "${PERSIST_ENV}" -eq 1 ]]; then
  BASHRC="${HOME}/.bashrc"
  MARK_BEGIN="# >>> pluto llvm >>>"
  MARK_END="# <<< pluto llvm <<<"
  TMP_FILE="$(mktemp)"

  if [[ -f "${BASHRC}" ]]; then
    awk -v b="$MARK_BEGIN" -v e="$MARK_END" '
      $0==b {skip=1; next}
      $0==e {skip=0; next}
      skip!=1 {print}
    ' "${BASHRC}" > "${TMP_FILE}"
  fi

  {
    [[ -f "${TMP_FILE}" ]] && cat "${TMP_FILE}"
    echo "${MARK_BEGIN}"
    echo "export HYPERF_HOME=\"${HYPERF_HOME}\""
    if [[ -n "${FILECHECK_BIN}" ]]; then
      echo "export PATH=\"${FILECHECK_BIN}:\${PATH}\""
    fi
    echo "export LLVM_PREFIX=\"${LLVM_PREFIX}\""
    echo 'export PATH="${LLVM_PREFIX}/bin:${PATH}"'
    echo 'if [ -d "${LLVM_PREFIX}/lib64" ]; then export LD_LIBRARY_PATH="${LLVM_PREFIX}/lib64:${LD_LIBRARY_PATH:-}"; fi'
    echo 'if [ -d "${LLVM_PREFIX}/lib" ]; then export LD_LIBRARY_PATH="${LLVM_PREFIX}/lib:${LD_LIBRARY_PATH:-}"; fi'
    echo "${MARK_END}"
  } > "${BASHRC}"
  rm -f "${TMP_FILE}"

  echo "[pluto-install] wrote LLVM exports to ${BASHRC}"
fi

echo "[pluto-install] done"
