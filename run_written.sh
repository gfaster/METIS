#!/usr/bin/env bash

# small utility to run the "to partition" command from a GraphBuilder

set -euo pipefail

USAGE="usage: $0 <GRAPH> [rust|env|diff|norm(default)|norm-gdb|rust-gdb]"
FILE="${1:?"$USAGE"}"

# gets the program name
CMDLINE=$(awk '/^%/ { $1 = ""; bin = $0 }; /^[^%]/ { next }; END { print bin }' "$FILE" | cut -d' ' -f 2-)

# echo "$CMDLINE"

NAME="$(echo "$CMDLINE" | cut -d' ' -f 1)"

ARGS="$(echo "$CMDLINE" | cut -d' ' -f 2-)"

case "$NAME" in
	'gpmetis' | 'mpmetis' | 'ndmetis' | 'm2gmetis')
		;;
	*)
		echo "Last comment is not a metis program (\`$NAME'), was it created with GraphBuilder?"
		exit 1
esac

OP="${2:-norm}"
case "$OP" in
	'rust')
		cargo b
		BIN_BASE="target/debug/"
		;;
	'rust-gdb')
		cargo b
		BIN_BASE="rust-gdb --args target/debug/"
		;;
	'diff')
		cargo b
		: ${METIS_NORM:?"METIS_NORM is not set, are you in nix-shell?"}
		;;
	'env')
		BIN_BASE=""
		;;
	'norm')
		: ${METIS_NORM:?"METIS_NORM is not set, are you in nix-shell?"}
		BIN_BASE="$METIS_NORM/bin/"
		;;
	'norm-gdb')
		: ${METIS_NORM:?"METIS_NORM is not set, are you in nix-shell?"}
		BIN_BASE="norm-gdb --args $METIS_NORM/bin/"
		;;
	*)
	echo "$USAGE"
	exit 1
	;;
esac

if [ "$OP" = 'diff' ]; then
	diff -d --color=always -s -b -a -y -W130 \
		--label "norm" <(eval $(printf "$METIS_NORM/bin/%s %s" "$NAME" "$ARGS")) \
		--label "rust" <(eval $(printf "target/debug/%s %s" "$NAME" "$ARGS"))
	exit 0
fi

# same as one defined in shell.nix
function norm-gdb() {
      gdb -d "$METIS_NORM_SRC/programs" \
	      -d "$METIS_NORM_SRC/libmetis" \
	      -d "$METIS_NORM_SRC/include" \
	      -x gdb/options_pp.scm "$@"
}



CMD=$(printf "%s%s %s" "$BIN_BASE" "$NAME" "$ARGS")

printf "%s\n" "$CMD"

eval $CMD
