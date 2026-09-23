#!/bin/bash
# Double-click: installs/updates the environment if needed, then runs the offline self-test on the bundled
# LMNA example (every analysis and assistant tool, no internet). Nothing is started.
cd "$(dirname "$0")" || exit 1
APP="$(pwd)/Structure Bench.app/Contents/Resources/app"
[ -d "$APP" ] || APP="$(cd .. && pwd)"
MODE=terminal exec /bin/bash "$APP/launch.sh" --selftest
