#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

echo "-----------"
echo "0. Clean temp files"
rm -rf .snakemake/shadow
rm -rf tmp
mkdir -p tmp
echo "[DONE]"

echo "-----------"
echo "1. pipeline intermediate outputs not affected by --touch"
echo "[SKIP]"

echo "-----------"
echo "2. pipeline logs are not handled by --touch"
find logs -type f -exec touch {} +
echo "[DONE]"

echo "-----------"
echo "3. touch snakemake pipeline results:"
snakemake --use-conda --cores 1 -q --touch --forceall -n
snakemake --use-conda --cores 1 --touch --forceall
echo "[DONE]"

echo "-----------"
echo "4. rules outputs marked as directories only"
# Snakemake updates `.snakemake_timestamp` file timestamp, find other files in such directories and touch them:
find results -name ".snakemake_timestamp" -exec dirname {} + | xargs -I FOLDER find FOLDER -type f \! -name ".snakemake_timestamp" -exec touch -h {} +
echo "[DONE]"

echo "-----------"
echo "5. .snakemake folder files:"
find .snakemake \( -type f -or -type l \) -exec touch -h {} +
echo "[DONE]"

echo "-----------"
echo "6. pipeline input (because it is copied to 'reads' folder)"
find reads -type f -exec touch -h {} +
echo "[DONE]"

echo "-----------"
echo "7. pipeline sources:"
# E.g. .git, workflow, config, images, envs, ..
SRC_FOLDERS=".git
config
images
workflow"

for F in $SRC_FOLDERS; do
  echo "  touch all in: $F/"
  find "$F" -type f -exec touch -h {} +
done

# E.g. .gitignore, environment.yaml, LICENSE.md, README.md, Snakefile, ...
SRC_FILES="./.gitignore
./environment.yaml
./LICENSE.md
./LICENSE.md
./README.md"
for F in $SRC_FILES; do
  echo "  touch: $F"
  touch -h "$F"
done
echo "[DONE]"

echo "-----------"
echo "8. pipeline extra files:"
touch ./lsf_docker_env_file.env
echo "[DONE]"

echo "-----------"
echo "9. CHECK files >25 days old: (first 10 files)"
#find . \! -type d -mtime +25 | head

# or ignoring links (RIS doesn't check symlinks modification data):
echo "* Files number:"
find . \! \( -type d -or -type l \) -mtime +25 | wc -l
echo "* First 10:"
find . \! \( -type d -or -type l \) -mtime +25 | head

echo "-----------"
echo "10. CHECK files >2 days old: (first 10 files)"
#find . \! -type d -mtime +2 | head

# or ignoring links:
echo "* Files number:"
find . \! \( -type d -or -type l \) -mtime +2 | wc -l
echo "* First 10:"
find . \! \( -type d -or -type l \) -mtime +2 | head

echo "====== DONE ==========="
