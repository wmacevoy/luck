#!/usr/bin/env bash
# Scan a LaTeX .log for warnings.
# Fatal patterns => print, delete the .pdf, exit 1 (build fails).
# Soft patterns  => print as a summary so they don't get lost in build noise.

set -u

log=${1:?usage: check-log.sh <log> <pdf>}
pdf=${2:?usage: check-log.sh <log> <pdf>}

if [ ! -f "$log" ]; then
  echo "check-log: $log not found" >&2
  exit 1
fi

# Silent content loss or broken output. A match in any of these = build fails.
fatal_re='Missing character|LaTeX Warning: (Reference|Citation) .* undefined|There were undefined (references|citations)|LaTeX Font Warning: Font shape .* undefined|Some font shapes were not available'

# Other things worth surfacing but not failing on (typography, deprecations).
soft_re='^(Overfull|Underfull)|LaTeX( Font)? Warning:|Package [A-Za-z]+ Warning:'

# Known-noisy lines we deliberately ignore.
soft_ignore='inputenc package ignored with utf8 based engines|hooks Warning: Generic hook .* is deprecated|csquotes Warning: No multilingual support'

fatal=$(grep -E "$fatal_re" "$log" 2>/dev/null | sort -u || true)
if [ -n "$fatal" ]; then
  echo
  echo "ERROR: fatal warnings in $log (build would silently lose content):"
  echo "$fatal" | sed 's/^/  /'
  rm -f "$pdf"
  exit 1
fi

soft=$(grep -E "$soft_re" "$log" 2>/dev/null | grep -vE "$soft_ignore" | sort -u || true)
if [ -n "$soft" ]; then
  echo
  echo "warnings in $log (non-fatal):"
  echo "$soft" | sed 's/^/  /'
fi
