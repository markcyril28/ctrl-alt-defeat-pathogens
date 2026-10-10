#!/usr/bin/env bash
# clone_repo_with_all_branches.sh — clone a repository and end up with every one of its branches
# Run from: Ubuntu/WSL terminal  or  Mac Terminal, in the folder where you want the copy to appear
# Usage:    bash clone_repo_with_all_branches.sh <repository-url>              # clone into a folder named after the repo
#           bash clone_repo_with_all_branches.sh <repository-url> my_folder    # clone into my_folder instead
#           bash clone_repo_with_all_branches.sh <repository-url> --no-local   # leave the branches as origin/<name> only
#           bash clone_repo_with_all_branches.sh <repository-url> --mirror     # bare copy of every ref, for a backup
#           bash clone_repo_with_all_branches.sh --fix                         # repair the clone you are standing in
#
# Why this exists: "git clone" already downloads every branch, but it checks out one of them and leaves the
# rest as remote-tracking names (origin/metagenomics), which plain "git branch" does not list. So a fresh
# clone looks like it only has main, and people clone again once per branch. This script does the one missing
# step — it gives each remote branch a local branch of the same name, tracking the remote one — so that
# "git branch" lists them all and "git switch metagenomics" works straight away.
#
# What a normal run does, in order:
#   1. refuses to start if the destination folder already exists with something in it
#   2. clones with --no-single-branch and --tags, so the history of every branch and every tag comes down
#   3. fetches once more with --prune, so the remote-tracking branches are exactly what the remote has
#   4. creates a local branch for each remote branch, except the one git already checked out
#   5. prints every branch with the remote it tracks, and the command to move between them
#
# --mirror is the other kind of complete copy: a bare repository (no working tree, no files to edit) holding
# every ref the remote has. That is the one to use for a backup, or to push a repository to a new host:
#   git clone <the mirror folder> working_copy      # get a normal working copy back out of it
#   git push --mirror <new-remote-url>              # move the whole repository somewhere else
#
# --fix is for a clone you already have that only brought one branch down, which happens when it was made
# with --single-branch, or with --depth (shallow clones imply a single branch). It rewrites that clone's
# fetch setting to cover all branches, unshallows it if needed, fetches, and then creates the local branches.
# Point it at nothing — it works on the repository the terminal is currently inside.
#
# Nothing here changes any commit: it clones, fetches, and adds branch names. It never deletes a branch.

set -euo pipefail

# ── Settings ────────────────────────────────────────────────────────
MAKE_LOCAL=1      # 1 = give every remote branch a local branch of the same name   --no-local turns this off
MIRROR=0          # 1 = bare mirror of every ref instead of a working copy         --mirror
FIX=0             # 1 = work on the existing clone here, do not clone anything     --fix
URL=""            # the repository to clone
DEST=""           # the folder to clone into (default: named after the repository)

# ── Helpers ─────────────────────────────────────────────────────────
GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  BLUE='\033[0;34m'  NC='\033[0m'

info()  { printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn()  { printf "${YELLOW}[SKIP]${NC}  %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC}  %s\n" "$*"; }
step()  { printf "\n${BLUE}==>${NC} %s\n" "$*"; }

# ── Read the options ────────────────────────────────────────────────
for arg in "$@"; do
    case "$arg" in
        --no-local)   MAKE_LOCAL=0 ;;
        --mirror)     MIRROR=1 ;;
        --fix)        FIX=1 ;;
        -h|--help)    sed -n '2,8p' "$0"; exit 0 ;;
        -*)           fail "unknown option: $arg  (try --help)"; exit 2 ;;
        *)            if [[ -z "$URL" ]]; then URL="$arg"; else DEST="$arg"; fi ;;
    esac
done

if (( ! FIX )) && [[ -z "$URL" ]]; then
    fail "give the repository to clone, for example:"
    echo "    bash $(basename "$0") https://github.com/markcyril28/Ologist_Workshop.git"
    echo "    bash $(basename "$0") --fix        # instead, repair the clone you are standing in"
    exit 2
fi

# ── Step 1: get a repository to work in ─────────────────────────────
# Either the clone this script makes, or — with --fix — the one the terminal is already inside.
if (( FIX )); then
    step "Repairing the clone you are standing in"
    REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo "")"
    if [[ -z "$REPO_ROOT" ]]; then
        fail "this folder is not inside a git repository — run --fix from inside the clone you want repaired"
        exit 1
    fi
    cd "$REPO_ROOT"
    info "repository: $REPO_ROOT"
else
    # Work out the folder name from the url when none was given: the last part, without the .git ending.
    if [[ -z "$DEST" ]]; then
        DEST="$(basename "$URL" .git)"
        (( MIRROR )) && DEST="$DEST.git"
    fi

    if [[ -e "$DEST" ]] && [[ -n "$(ls -A "$DEST" 2>/dev/null)" ]]; then
        fail "$DEST already exists and is not empty — move it aside, or give another folder name"
        exit 1
    fi

    if (( MIRROR )); then
        step "Mirroring $URL into $DEST"
        # --mirror brings down every ref there is: branches, tags, notes. No working tree, so no files to edit.
        git clone --mirror "$URL" "$DEST"
        cd "$DEST"
        info "bare mirror written to $(pwd)"
        step "Summary"
        echo "  branches in the mirror:"
        git for-each-ref --format='    %(refname:lstrip=2)' refs/heads/
        echo "  tags: $(git for-each-ref refs/tags/ | wc -l | tr -d ' ')"
        echo
        echo "A mirror has no files to open. To work in it, or to move it elsewhere:"
        echo "    git clone $(pwd) working_copy          get a normal working copy out of it"
        echo "    git --git-dir=$(pwd) log --oneline     read its history in place"
        echo "    cd $(pwd) && git push --mirror <url>   push the whole repository to a new host"
        exit 0
    fi

    step "Cloning $URL into $DEST"
    # --no-single-branch is the default already, but spelling it out says what the script is for, and keeps
    # the clone correct if git's own defaults are changed by a ~/.gitconfig somewhere.
    git clone --no-single-branch --tags "$URL" "$DEST"
    cd "$DEST"
    REPO_ROOT="$(pwd)"
    info "repository: $REPO_ROOT"
fi

# The remote is usually called origin, but not always — this repository's own remote is called main.
REMOTE="$(git remote | head -1)"
if [[ -z "$REMOTE" ]]; then
    fail "this repository has no remote, so there are no other branches to fetch"
    exit 1
fi
info "remote: $REMOTE"

# ── Step 2: make sure the clone covers every branch, not one ────────
step "Fetching every branch and tag from $REMOTE"

# A --single-branch clone has a fetch setting that names one branch; this widens it to all of them.
# Harmless to run on a clone that was already complete.
git remote set-branches "$REMOTE" '*'

# A shallow clone has only the tip commits, and cannot hand out full branch histories until it is unshallowed.
if [[ "$(git rev-parse --is-shallow-repository)" == "true" ]]; then
    warn "this clone is shallow — fetching the full history, which takes longer"
    git fetch --unshallow "$REMOTE"
fi

git fetch --all --tags --prune

# ── Step 3: give every remote branch a local branch ─────────────────
# This is the step a plain clone leaves out. Branch names never contain spaces, so a simple for-loop is safe.
CREATED=()     # local branches this run made
ALREADY=()     # local branches that were already here
LEFT_OUT=()    # refs that are not branches of this repository

if (( MAKE_LOCAL )); then
    step "Creating a local branch for each branch on $REMOTE"
    CURRENT="$(git symbolic-ref --short HEAD 2>/dev/null || echo "")"

    for ref in $(git for-each-ref --sort=refname --format='%(refname)' "refs/remotes/$REMOTE/"); do
        branch="${ref#"refs/remotes/$REMOTE/"}"

        # origin/HEAD is a pointer to the default branch, not a branch, and GitHub's pr/<number>/head refs
        # are copies of other people's branches — neither should become a local branch here.
        if [[ "$branch" == "HEAD" ]] || [[ "$branch" == pr/* ]]; then
            LEFT_OUT+=("$REMOTE/$branch")
            continue
        fi

        if git show-ref --verify --quiet "refs/heads/$branch"; then
            if [[ "$branch" == "$CURRENT" ]]; then
                info "$branch is the branch git checked out"
            else
                warn "$branch is already a local branch"
            fi
            ALREADY+=("$branch")
            continue
        fi

        git branch --quiet --track "$branch" "$ref"
        info "created $branch, tracking $REMOTE/$branch"
        CREATED+=("$branch")
    done
else
    step "Leaving the branches as $REMOTE/<name> only (--no-local)"
    git for-each-ref --sort=refname --format='  %(refname:lstrip=2)' "refs/remotes/$REMOTE/"
fi

# ── Summary ─────────────────────────────────────────────────────────
step "Summary"
echo "  repository:        $REPO_ROOT"
echo "  remote:            $REMOTE"
echo "  branches created:  ${CREATED[*]:-none}"
echo "  already here:      ${ALREADY[*]:-none}"
echo "  not branches:      ${LEFT_OUT[*]:-none}"
echo "  tags:              $(git for-each-ref refs/tags/ | wc -l | tr -d ' ')"

if (( MAKE_LOCAL )); then
    echo
    echo "Every branch of the repository is now here:"
    git branch -vv | sed 's/^/    /'
fi

echo
echo "To move between them:"
echo "    git branch -a                     list local and remote branches"
echo "    git switch <branch>               move to one; your files change to that branch's version"
echo "    git log --oneline --graph --all    see how the branches relate"
echo "    git fetch --all --prune           later on, bring down what changed on the remote"
