#!/usr/bin/env bash

set -euo pipefail

set -x

FILE="${1:?'requires file'}"

if [ ! -f "$FILE" ]; then
	echo "$FILE does not exist!"
	exit 1
fi

if [ ! -z "$(git diff "$FILE")" ]; then
	echo "$FILE is dirty!"
	exit 1
fi

trap 'git restore "$FILE"' EXIT

function s() {
	: ${1:?'requires sed expression'}
	sed -iE -e "$1" "$FILE"
}

function v() {
	: ${1:?'requires sed expression'}
	nvim -Es "$FILE" -c "$1" -c ":w"
}



s 's/#include.*//'
s 's/where/where_/g'
s 's/->/./g'
s 's/\/\*!/\/*/'

# Functions on one line
v ':g/\v^\w+ \**\w+\([^)]+$/norm vibVJ'

# Function argument syntax
v ':g/\v^\w+/:s/\v(\w+_t) (\**)(\w+)([,)])/\3: \2\1\4/g'
v ':g/\v^\w+/:s/\V*/*mut /g'


# slice indexing
v ':%s/\v\[([^=\]]*[-+*/ ][^=\]]*)\]/[(\1) as usize]/g'
v ':%s/\v\[(\w+)\]/[\1 as usize]/g'
v ':%s/\]\]/] as usize]/' || true

# the rest of the function signature
v ':%s/\v^(\w+) (%(\*mut )*)([^{]*)/#[metis_func]\rpub extern "C" fn \3 -> \2\1/'

# bracketless if/for
v ':%s/\v(^ *%(for|if) ?\([^{]*)\n(.*;)/\1 {\r\2\r}'
v ':%s/\v(^ *else *\n([^{]*)\n(.*;)/\1 {\r\2\r}' || true


# switch -> match
s 's/switch/match/'
v ':%s/\vcase (\w+):/\1 =>' || true
v ':%s/default:/_ =>/' || true

# comment out variable declarations
v ':g/\v^  (int|(\w*_t)) /norm gcc'
s 's/^ *nvtxs *= graph\.nvtxs/let & as usize/'
s 's/^ *ncon *= graph\.ncon/let & as usize/'

# for loops
# extract extra initial variables
v ':%s/for (\s*\([a-zA-Z]\+\s*=\s*.\w*\)\s*,\s*/\1;\rfor (/' || true
v ':%s/for (\s*\([a-zA-Z]\+\s*=\s*.\w*\s*\),\s*/\1;\rfor (/' || true
# actual loops
v ':%s/\vfor \((\s*[a-zA-Z]+)\s*\=\s*(.{-1,});;\s*\1\+\+\)/for \1 in (\2)..' || true
v ':%s/\vfor \((\s*[a-zA-Z]+)\s*\=\s*(.{-1,});\s*\1\s*\<\=\s*(.{-1,});\s*\1\+\+\)/for \1 in (\2)..=(\3)' || true
v ':%s/\vfor \((\s*[a-zA-Z]+)\s*\=\s*(.{-1,});\s*\1\s*\<\s*(.{-1,});\s*\1\+\+\)/for \1 in (\2)..(\3)' || true
v ':%s/\vfor \((\s*[a-zA-Z]+)\s*\=\s*(.{-1,});\s*\1\s*\>\=\s*(.{-1,});\s*\1\-\-\)/for \1 in ((\3)..=(\2)).rev()' || true
v ':%s/\vfor \((\s*[a-zA-Z]+)\s*\=\s*(.{-1,});\s*\1\s*\>\s*(.{-1,});\s*\1\-\-\)/for \1 in ((\3)..(\2)).rev()' || true

# make sure for loops are gone
! grep -qE '\s*for.*;' "$FILE"

# post/pre-fix increment/decrement
s 's/++;\s*$/ += 1;/'

v ':g/++\w\+/norm f+xxyiw0pwi+=1;' || true
v ':g/\w\+++/norm f+xxhyiw$pA+=1;' || true

v ':g/--\w\+/norm f-xxyiw0pwi-=1;' || true
v ':g/\w\+--/norm f-xxhyiw$pA-=1;' || true

# printing
v ':%s/\v\%(.{-0,5})"%(PRIDX|PRREAL)"/{:\1}/g' || true
s 's/%s/{}/'
s 's/\<printf\>/print!/'

# temp allocations
v ':%s/\viset\(([^,]+), ([^,]+), iwspacemalloc\(ctrl, .+\)\);/vec![\2; \1 as usize];/g' || true
v ':%s/\viwspacemalloc\(ctrl, (.*)\);/vec![0; \1 as usize];/g' || true
v ':%s/\vrwspacemalloc\(ctrl, (.*)\);/vec![0.0; \1 as usize];/g' || true
v ':g/WCORE\%\(PUSH\|POP\)/norm gcc' || true

# misc
s 's/-> void//' || true
s 's/NULL/std::ptr::null_mut()/g' || true
s 's/ASSERTP/debug_assert!/'
s 's/ASSERT/debug_assert!/'
s 's/IFSET/ifset!/'


git diff "$FILE"

echo 'does this look good?'
select yn in Yes No; do
	case $yn in
		Yes ) 
			trap - EXIT
			exit 0
			;;
		No ) 
			exit 1
			;;
	esac
done
exit 1
