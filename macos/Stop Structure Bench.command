#!/bin/bash
# Stops the Structure Bench server if it is running (same as "Quit Structure Bench" in the sidebar).
PORT="${PORT:-8767}"
if curl -s -m 2 "http://localhost:$PORT/api/settings" >/dev/null 2>&1; then
  curl -s -m 2 -X POST "http://localhost:$PORT/api/quit" >/dev/null 2>&1; sleep 1
  lsof -ti tcp:"$PORT" 2>/dev/null | xargs kill 2>/dev/null
  echo "Structure Bench stopped."
else
  echo "Structure Bench was not running."
fi
sleep 1
