#!/bin/bash

set -e
python3 --version || { echo "Python 3 is required"; exit 1; }
python3 -c "import Bio" 2>/dev/null || pip3 install biopython
