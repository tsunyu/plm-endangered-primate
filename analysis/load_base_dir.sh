# Resolve PLM_BASE_DIR / BASE_DIR for shell pipelines.
# Source from a script in analysis/<phase>/ :
#   source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"

_ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${PLM_BASE_DIR:-}" && -f "${_ANALYSIS_DIR}/base_dir.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "${_ANALYSIS_DIR}/base_dir.env"
  set +a
fi

if [[ -z "${PLM_BASE_DIR:-}" ]]; then
  echo "PLM_BASE_DIR is not set." >&2
  echo "Run:  bash configure_base_dir.sh /path/to/analysis_root" >&2
  echo "  or: export PLM_BASE_DIR=/path/to/analysis_root" >&2
  echo "The analysis root is the directory that contains data/ and output/." >&2
  exit 1
fi

BASE_DIR="${PLM_BASE_DIR}"
export PLM_BASE_DIR BASE_DIR
export PLM_ANALYSIS_DIR="${_ANALYSIS_DIR}"
unset _ANALYSIS_DIR
