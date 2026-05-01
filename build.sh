#!/usr/bin/env bash
set -euo pipefail

IMAGE=luck-latex
BASE=devcontainer-latex:latest
BASE_DIR=../devcontainer-latex

cd "$(dirname -- "${BASH_SOURCE[0]}")"

if ! docker image inspect "$BASE" >/dev/null 2>&1; then
  echo "base image $BASE not found locally; building it from $BASE_DIR..."
  ( cd "$BASE_DIR" && ./build.sh )
fi

docker build -q -t "$IMAGE" -f .devcontainer/Dockerfile .devcontainer/

docker run --rm -v "$(pwd):/workspace" -w /workspace "$IMAGE" \
  make -B "${@:-all}"
