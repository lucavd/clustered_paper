#!/usr/bin/env bash
# public_sync.sh -- public copies of this private repository, updated in parallel with it.
#   paper   -> https://github.com/lucavd/clustered_paper : allow-list of paths of branch morris-v4 (PATHS_PAPER)
#   package -> https://github.com/lucavd/fewclust        : directory fewclust/ of branch package (~/MMM-package)
# The public histories start empty and get one commit per sync ("Sync from private repository, commit <hash>");
# the private history (which contains the individual data of the neonatal cohort) is never pushed.
# Usage: bash tools/public_sync.sh paper|package [--push]
#   without --push: builds the staging copy in ~/public_staging/<repo>, runs the checks, writes the manifest;
#   with --push: the same, then commits and pushes to the public repository (branch main).
# Checks (any failure stops before the commit): no file of the cohort data (data_raw/, LURID, lurid_patient, .sav,
#   .xlsx), no row header of the cohort variables, no Telegram token or other credential pattern, no file over 50 MB.
set -euo pipefail
WHAT=${1:?paper or package}; PUSH=${2:-}
ST=/home/user/public_staging; mkdir -p $ST
PATHS_PAPER="R py tools tables manuscript reports/morris docs .rundir MMM_sim_manual rev2 baseline renv/activate.R renv/settings.json renv/.gitignore renv.lock .Rprofile LICENSE .gitattributes MMM.Rproj"
case $WHAT in
  paper)   SRC=/home/user/MMM; REF=morris-v4; REPO=clustered_paper; SUB="" ;;
  package) SRC=/home/user/MMM-package; REF=package; REPO=fewclust; SUB="fewclust" ;;
  *) echo "unknown target $WHAT"; exit 2 ;;
esac
HASH=$(git -C $SRC rev-parse --short $REF)
D=$ST/$REPO
[ -d $D/.git ] || git clone -q https://github.com/lucavd/$REPO $D
git -C $D fetch -q origin 2>/dev/null || true
if git -C $D rev-parse -q --verify origin/main >/dev/null; then git -C $D checkout -q -B main origin/main; else git -C $D checkout -q --orphan main 2>/dev/null || true; fi
find $D -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
if [ "$WHAT" = paper ]; then
  git -C $SRC archive --format=tar $REF -- $PATHS_PAPER | tar -x -C $D
  cp $SRC/tools/public/README_paper.md $D/README.md
  # Telegram notifications of the original runs: credentials replaced by environment variables
  grep -rlE 'Bot\(token *= *"' $D --include='*.R' | while read f; do
    sed -E -i 's/Bot\(token *= *"[^"]*"\)/Bot(token = Sys.getenv("TELEGRAM_BOT_TOKEN"))/; s/chat_id *<- *"[0-9-]+"/chat_id <- Sys.getenv("TELEGRAM_CHAT_ID")/' "$f"; done
else
  git -C $SRC archive --format=tar $REF $SUB | tar -x -C $D --strip-components=1
fi
# ---- checks ----
FAIL=0
BAD=$(cd $D && find . -path ./.git -prune -o -type f \( -ipath '*data_raw/*' -o -iname '*lurid*' -o -iname '*.sav' -o -iname '*.xlsx' \) -print)
[ -z "$BAD" ] || { echo "FAIL cohort files:"; echo "$BAD"; FAIL=1; }
HDR=$(cd $D && grep -rlIE '(^|[,;\t ])(NIP|GAweeks)[,;\t ].*(BWg|CRIBII|SNAPPE)|(^|,)id,age,age_weeks,disease' --exclude-dir=.git . || true)
[ -z "$HDR" ] || { echo "FAIL cohort headers:"; echo "$HDR"; FAIL=1; }
SEC=$(cd $D && grep -rlIE '[0-9]{8,10}:[A-Za-z0-9_-]{30,}|ghp_[A-Za-z0-9]{30,}|github_pa[t]_|sk-[A-Za-z0-9]{30,}|BEGIN [A-Z ]*PRIVATE KEY|AKIA[0-9A-Z]{16}|chat_id *<- *"[0-9]{5,}"' --exclude-dir=.git . || true)
[ -z "$SEC" ] || { echo "FAIL credential patterns:"; echo "$SEC"; FAIL=1; }
BIG=$(cd $D && find . -path ./.git -prune -o -type f -size +50M -print)
[ -z "$BIG" ] || { echo "FAIL files over 50 MB:"; echo "$BIG"; FAIL=1; }
M=$ST/manifest_$REPO.txt
{ echo "source: $SRC $REF $HASH"; echo "files: $(cd $D && find . -path ./.git -prune -o -type f -print | wc -l)  size: $(du -sh --exclude=.git $D | cut -f1)"
  (cd $D && find . -path ./.git -prune -o -type f -printf '%s %h\n' | awk '{split($2,a,"/"); t=(a[2]==""?"(root)":a[2]); n[t]++; s[t]+=$1} END{for(k in n) printf "  %6d files %9.1f MB  %s\n", n[k], s[k]/1e6, k}' | sort -k3 -nr)
  if [ "$WHAT" = paper ]; then echo "excluded top-level paths of $REF:"; git -C $SRC ls-tree --name-only $REF | while read p; do keep=0; for q in $PATHS_PAPER; do [ "${q%%/*}" = "$p" ] && keep=1; done; [ $keep = 1 ] || echo "  $p"; done
    echo "redacted (Telegram credentials -> environment variables):"; (cd $D && grep -rlE 'Sys.getenv\("TELEGRAM_BOT_TOKEN"\)' --include='*.R' . | sed 's/^/  /'); fi
  echo "checks: $([ $FAIL = 0 ] && echo pass || echo FAIL)"; } > $M
cat $M
[ $FAIL = 0 ] || exit 1
if [ "$PUSH" = --push ]; then
  cd $D; git add -A
  if git diff --cached --quiet; then echo "nothing to commit"; else
    git -c user.name="Luca Vedovelli" -c user.email="$(git -C $SRC config user.email)" commit -q -m "Sync from private repository, commit $HASH (branch $REF)"
    git push -q origin main; echo "pushed $(git log --oneline -n 1)"; fi
fi
