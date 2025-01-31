# METIS C to Rust translation steps

**NOTE**: I'm using this file to transcribe VIM macros, so there are a few non-printable characters

I recommend recording a macro that yanks and executes a command from this
document.

## 0. Setup

Before anything else, `git mv` the file to port to `src/ported` and then copy
the file back into `src` with an `.rs` ending. Then, add the file in
`src/lib.rs`.

## 1. Replace most necessary tokens

First, run these commands
```
:%s/#include.*//
:%s/where/where_/g
:%s/->/./g
:%s/\/\*!/\/*
```

Properly accounting for `->` replacement is done in a later step.

## 2. Convert function declarations
first, run this command to make all function declarations one line
```
:g/\v^\w+ \**\w+\([^)]+$/norm vibVJ
```

Then, to make the argument format correct
```
:g/\v^\w+/:s/\v(\w+_t) (\**)(\w+)([,)])/\3: \2\1\4/g
:g/\v^\w+/:s/\V*/*mut /g
```

Then finally, to declare the function and set the return type:
```
:%s/\v^(\w+) (%(\*mut )*)([^{]*)/#[metis_func]\rpub extern "C" fn \3 -> \2\1/
```


## 3. Fix bracketing and blocks
replace most bracketless if and for statements
```
g/\v^ *(for|if) ?\([^{]*\n.*;/norm jVSB
g/\v^ *else *\n[^{]*;/norm jVSB
```

Replacing `switch` with `match`
```
:s/switch/match
:%s/\vcase (\w+):/\1 => {
:%s/default:/_ =>
```
Generally don't need brackets on `default` case, need to manually replace
`break` statements with `},`. We can't use regex since there might be a loop
with a `break` inside of it. Also be weary of fallthrough cases adding an extra
open bracket with the above substitution.

Most ternary statements can be fixed using this:
```
%s/\v\(([^?:]*)\?([^?:]*):([^?:]*)\);$/(if \1 { \2 } else { \3 });
```

## 4. Fix variable declarations
In each function:
- Comment out original variable decls.
- make `graph` and `ctrl` to mutable references
```rust
let graph = graph.as_mut().unwrap();
let ctrl = ctrl.as_mut().unwrap();
```
- for blocks of graph pointer references, use the `get_graph_slices` and
  `get_graph_slices_mut` macros.

Comment out C variable declarations:
```
:g/\v^  (int|(\w*_t)) /norm gcc
```

Convert:
```c
idx_t i, ii, j, k, nvtxs, pwgts[2], zeromaxpwgt, from, me, 
      bestcut=0, icut, mincut, inbfs;
idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_;
idx_t *perm, *bestwhere;

// the arrows will actually already have been replaced with '.' in step 1
nvtxs  = graph->nvtxs;
xadj   = graph->xadj;
vwgt   = graph->vwgt;
adjncy = graph->adjncy;
adjwgt = graph->adjwgt;
```

To:

```rust
// idx_t i, ii, j, k, nvtxs, pwgts[2], zeromaxpwgt, from, me, 
//       bestcut=0, icut, mincut, inbfs;
// idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_;
// idx_t *perm, *bestwhere;

let nvtxs  = graph.nvtxs;
get_graph_slices!(graph => xadj vwgt adjncy adjwgt);
```
Side note: it will probably save you some hassle later to make the local
binding of `nvtxs` as a `usize` since it's used in slice indexing a ton.
```rust
let nvtxs = graph.nvtxs as usize;
```

There's more we'll do later, but we need to fix all the syntax error first.

## 5. Convert for loops
replacing `iter` with relevant function decls:
```
:%s/\vfor \(iter\=(.{-1,}); iter\<(.{-1,}); iter\+\+\)/for iter in (\1)..(\2)
```

Also need to replace infinite loops:
```
:%s/for (;;)/loop/
```

Now that loop iterators are gone (double check this), we can get rid of unary
increment and decrement
```
:g/++\w\+/norm f+xxyiw0pwi+=1;
:g/\w\+++/norm f+xxhyiw$pA+=1;

:g/--\w\+/norm f-xxyiw0pwi-=1;
:g/\w\+--/norm f-xxhyiw$pA-=1;
```

Print statment formatting can be mostly fixed with this:
```
:%s/\v\%(.{-0,5})"%(PR%(IDX|REAL))"/{:\1}/g
```

With this, all syntax errors should be fixed. After you run `cargo fmt`, this
is the end of the first pass.

## 6. Handle local heap allocations

METIS handles local heap allocations (that are freed at end of scope) as follows:
```C
WCOREPUSH;
var = iwspacemalloc(ctrl, cnt);
// ...
WCOREPOP;
```
All workspace allocations will be freed at end of scope. `gk_malloc`
allocations will be freed on `gk_malloc_cleanup`, but they can be freed
earlier. It'll be fine to replace all these calls at the end.



```
:%s/\viset\(([^,]+), (-?\d), iwspacemalloc\(ctrl, .+\)\);/vec![\2; \1 as usize];/g

:%s/\viwspacemalloc\(ctrl, (.*)\);/vec![0; \1 as usize];/g
:%s/\vrwspacemalloc\(ctrl, (.*)\);/vec![0.0; \1 as usize];/g
```
Note there is no use of `rset` in the codebase.

We then go over all those new `Vec`s and verify the type is correct and then
add a `let mut` at the beginning.

## 7. Remaining variable declarations.

We want to declare variables in the narrowest scope possible. Declare variables
at first assignment. For slices, try to do it at the function scope.

Use the `mkslice` and `mkslice_mut` macros for slices. Generally, the length of
the slice is revealed by checking a loop that it's accessed in. Members of
`graph_t` should probably use its macro instead (detailed in step 4).

Example:

```rust
// This should have been commented out in step 4
// idx_t *arr;

// the syntax was converted in step 5
for i in 0..nvtxs {
    arr[i] = i;
}
```

Likely converts to:

```rust
mkslice_mut!(arr, nvtxs);
// expands to approximately:
// let arr = std::slice::from_raw_parts_mut(arr, {nvtxs} as usize);
```

`mkslice` will declare the variable passed to it. If there is an array as a
member of a struct pointer (or reference), you can use the following syntax:

```rust
// you should propably use get_graph_slices!() instead, but this is an example.
mkslice!(graph->xadj, nvtxs);
// expands to approximately:
// let xadj = std::slice::from_raw_parts((*graph).xadj, {nvtxs} as usize);
```

becuase slice indexing is like, super annoying, try this:
```
:%s/\v\[([^=\]]*[-+*/ ][^=\]]*)\]/[(\1) as usize]/g
:%s/\v\[(\w+)\]/[\1 as usize]/g
:%s/\]\]/] as usize]/
:%s/\#\[metis_func as usize\]/#[metis_func]
```

To change the name (also works without field access):
```rust
mkslice!(cwhere: cgraph->where_, nvtxs);
// expands to approximately:
// let cwhere = std::slice::from_raw_parts((*cgraph).where_, {nvtxs} as usize);
```

Finally, at the top of the file, replace the include with:
```rust
use crate::*;
```

## 8. Replacing function calls
first, we want to replace function calls to existing modules with their qualified syntax:

for example: 

```rust
ProjectKWayPartition( /* ... */ );
```

To

```rust
kwayrefine::ProjectKWayPartition( /* ... */ );
```

Next, execute these replacements. The first group is nearly universal, but you
can pick and choose for the second.

```
:%s/IFSET/ifset!/I
:g/ifset.*gk_\%\(start\|stop\)/norm f(lvibgc
:%s/^ *\(goto \w*\);/todo!("\1");
:%s/ASSERTP\=/assert!/I
:%s/sprintf/write!/g
:%s/\v([^\w])printf/\1println!/g
:%s/\\n"/"/g
:g/\vgk_%(start|stop)cputimer/normal gcc
:g/WCORE\%\(PUSH\|POP\)/norm gcc
:%s/\vgk_errexit\(\w+, =/panic!(/
:%s/\vsizeof\((\w+)\)/std::mem::size_of::<\1>()/g
:%s/-> void/
:%s/^ *nvtxs *= graph\.nvtxs/let & as usize
:%s/^ *ncon *= graph\.ncon/let & as usize
:%s/NULL/std::ptr::null_mut()/g

:%s/\vBND\w+/&!/I
:%s/INC_DEC/\L&!/I
:%s/%d/{}/g
:%s/\vSWAP\((.{-}),(.{-}),.{-}\)/std::mem::swap(\&mut \1, \&mut \2)/I
```

At this point, search diagnostics for syntax errors, and fix them until there is
no more. Then, you can run `cargo fmt` to make it nicer to work on the remaining

After all functions are replaced, we need to declare the remaining in
`bindings.rs`. Here's the best way I found to do it:

over each identifer that was just pasted into the bindings block, execute this
```
:exe ":norm *N:e src/proto.h\<Enter>nvibVy:e#\<Enter>Vp"
```

then fix the format by highlighting the new functions and executing this
```
:'<,'>:g/\v\([^)]+$/join
2@:
:'<,'>:s/\v([,(] ?\w+ ?\**)([,)])/\1 _\2/g
:'<,'>:s/\v(\w+) (\**)(\w*)([,)])/\3: \2\1\4/g
:'<,'>:s/\V*/*mut /g
:'<,'>:s/\v^(\w+) (%(\*mut )*)(.*);/pub fn \3 -> \2\1;/
:'<,'>:s/void/std::ffi::c_void/g
```

### Some additional out of place functions:

`iargmax` is at `util::iargmax` and takes a slice and stride

A conservative replace would be this:
```
:%s/\viargmax\((.{-}),(.{-}),(.{-})\)/util::iargmax(std::slice::from_raw(\2, \1), \3)/
```
But since we've generally replaced pointers with slices, this will likely work:
```
:%s/\viargmax\(.{-},(.{-}),(.{-})\)/util::iargmax(\1, \2)/
```
---

The replacement for `MAKECSR` and `SHIFTCSR`, `util::make_csr` and
`util::shift_csr` is similar, but it takes a redundant variable that's used for
an iterator.

This replace is probably fine:

```
:%s/\vMAKECSR\(\w{-},(.{-}),(.{-})\);/util::make_csr(\1, &mut \2);/I
:%s/\vSHIFTCSR\(\w{-},(.{-}),(.{-})\);/util::shift_csr(\1, &mut \2);/I
```

The first argument in `util::make_csr` and `util::shift_csr` is possibly
redundant - see implementation for details.

### Some more useful replacements

Replacing block comments over function, the 2nd one is finicky.
```
:g/\v\/\*{10,}/norm dd
:g/^\/\*/norm WgbcgccI/
```

Replacing `iset` calls that aren't attached to any given `malloc` call:
```
:s/\viset\(.*,(.*),(.*)\);/\/\/ &\r\2.fill(\1);/
```

Replacing `icopy` calls (detailed in appendix)
```
:s/\vicopy\(.*,(.*),(.*)\);/\/\/ &\r\2.copy_from_slice(\&\1);/
```

## 9. Everything else

Still a lot to do get to it until it compiles. Every time errors get to zero
and jump back up is a pass.

Once we've committed a version that compiles, we can run `cargo fix` to help
with some simple style things like `if` statement parentheses.

## 10. Once we compile again

In order to test our changes, we utilize dynamic dispatch on a per function
basis.

There is no changes needed on the Rust side, but on the C version we need to
declare each function as an `ifunc` using the `IFUNC` macro defined in
`ifunc.h` (which should be included manually), and for the function definition
only to be prefixed with `c__libmetis__`.

See `src/ported/coarsen.c` for examples.

Macro to help:
```
0f(€ý5avabJyyPiIFUNC(f(€ý5bi,f(€ý5i, bbi A);j0f(€ý5bic__libmetis__q€kb
```

To use the old C functions, set `METIS_OVERRIDE_SYMS` to a comma separated list
of symbols with an optional colon to specify the language. Later occurrences
overwrite earlier ones. 

There is some magic in parsing symbols, and all of the following are equivalent:

```
METIS_OVERRIDE_SYMS="SetupCoarseGraph"
METIS_OVERRIDE_SYMS="SetupCoarseGraph:c"
METIS_OVERRIDE_SYMS="c__SetupCoarseGraph"
METIS_OVERRIDE_SYMS="libmetis__SetupCoarseGraph"
METIS_OVERRIDE_SYMS="libmetis__SetupCoarseGraph:c"
METIS_OVERRIDE_SYMS="c__libmetis__SetupCoarseGraph"
METIS_OVERRIDE_SYMS="c__libmetis__SetupCoarseGraph"
```

Which would be later reset to the Rust version if any of the following
equivalent settings are found:

```
METIS_OVERRIDE_SYMS="SetupCoarseGraph:rs"
METIS_OVERRIDE_SYMS="rs__SetupCoarseGraph"
METIS_OVERRIDE_SYMS="libmetis__SetupCoarseGraph:rs"
METIS_OVERRIDE_SYMS="rs__libmetis__SetupCoarseGraph"
```

If there's a `*` in the symbol, it will be treated as a glob. Exact matches
always have priority over globs, and the last glob matched is the version used.
Note that globs match against the base symbol, which excludes `c__libmetis__`.

```
# only use C versions of functions
METIS_OVERRIDE_SYMS="*"

# only use C versions of functions, except SetupCoarseGraph
METIS_OVERRIDE_SYMS="*,SetupCoarseGraph:rs"
```


## Things to look out for

- make sure that the `gk_malloc` calls have null-terminated strings

- we use a ton of casting to `usize`, but remember that `idx_t` is signed and
  often negative

- sometimes a loop iterator is used outside the loop, typically to do something if
  it didn't break. A for loop cannot be used here because the iterator will not
  ever equal the end bounds, even if it reaches the end. Use a while loop with
  increment at the end of the block.

- pwgts is sometimes a different length than what `get_graph_slices` says
