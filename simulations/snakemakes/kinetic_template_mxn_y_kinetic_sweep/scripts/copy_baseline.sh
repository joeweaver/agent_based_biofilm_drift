#!/usr/bin/env bash

set -euo pipefail

source=$(dirname "$1")
dest=$(dirname "$2")

echo $source
echo $dest

rsync -a --exclude '*.vt*' --exclude 'Inputscript.lmp.bak' --exclude '*.log' --exclude 'log.lammps' --exclude 'done.tkn' $source/ $dest

mv $dest/Inputscript.lmp $dest/Inputscript.orig
