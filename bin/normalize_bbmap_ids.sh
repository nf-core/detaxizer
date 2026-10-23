#!/usr/bin/env bash
#
# Normalize read ids for bbmap filterbyname.sh.
#
# filterbyname.sh matches the supplied names only against the first mate, so
# second-mate identifiers have to be mapped onto their first-mate counterpart.
# Supported forms (mate 2 -> mate 1):
#   CASAVA:        read/2
#   CASAVA 1.8+:   <coords> 2:N:<control>:<index>   (control number is any integer)
#   SRA/ENA:       read.2
# A first-mate id can itself end in .2, so that form keeps its original id next
# to the normalised one. Duplicates are removed and whitespace within a line is
# preserved.
#
# Usage: normalize_bbmap_ids.sh -o <output.txt> <ids.txt> [<ids2.txt> ...]
set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") -o <output.txt> <ids.txt> [<ids2.txt> ...]" >&2
    exit 1
}

output=""
inputs=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        -o)
            [ "$#" -ge 2 ] || usage
            output="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Error: unknown option '$1'" >&2
            usage
            ;;
        *)
            inputs+=("$1")
            shift
            ;;
    esac
done

if [ -z "$output" ] || [ "${#inputs[@]}" -eq 0 ]; then
    usage
fi

for f in "${inputs[@]}"; do
    if [ ! -f "$f" ]; then
        echo "Error: input file '$f' does not exist" >&2
        exit 1
    fi
done

awk '
{
    orig = $0
    sub(/\r$/, "", orig)
    if (orig == "") next

    line = orig
    ambiguous = 0

    # First token: CASAVA "/2", SRA/ENA ".2", or CASAVA 1.8+ mate placed first.
    if (match(line, /^[^ \t]+/)) {
        first = substr(line, 1, RLENGTH)
        rest  = substr(line, RLENGTH + 1)
        nf = first
        if (nf ~ /\/2$/) {
            sub(/\/2$/, "/1", nf)
        } else if (nf ~ /\.2$/) {
            sub(/\.2$/, ".1", nf)
            ambiguous = 1
        } else if (nf ~ /^2:[YN]:[0-9]+/) {
            sub(/^2:/, "1:", nf)
        }
        if (nf != first) line = nf rest
    }

    # Second token: CASAVA 1.8+ mate field when coordinates come first.
    if (match(line, /^[^ \t]+[ \t]+[^ \t]+/)) {
        pair = substr(line, 1, RLENGTH)
        rest = substr(line, RLENGTH + 1)
        if (match(pair, /[ \t]+/)) {
            t1 = substr(pair, 1, RSTART - 1)
            ws = substr(pair, RSTART, RLENGTH)
            t2 = substr(pair, RSTART + RLENGTH)
            nt2 = t2
            if (nt2 ~ /^2:[YN]:[0-9]+/) sub(/^2:/, "1:", nt2)
            if (nt2 != t2) line = t1 ws nt2 rest
        }
    }

    if (ambiguous && !seen[orig]++) print orig
    if (!seen[line]++) print line
}' "${inputs[@]}" > "$output"
