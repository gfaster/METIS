# METIS C to Rust translation steps

I recommend recording a macro that yanks and executes a command from this
document.

## 1. Replace most necessary tokens
At the top of the file, add the following line:
```rust
use crate::*;
```

Then run these commands
```
:%s/where/where_/g
:%s/->/./g
:%s/\*\/!/*\/
```
Properly accounting for `->` replacement is done in a later step.

## 2. Convert function declarations
first, run this command 4 times
```
:g/\v^\w+ ?\([^)]+$/join
```

Then, at each function declaration to transform arguments and pointers:
```
:s/\v(\w+_t) (\**)(\w+)([,)])/\3: \2\1\4/g
:s/\V*/*mut /g
```

Then finally, to declare the function and set the return type:
```
:%s/\v^(\w+) (%(\*mut )*)([^{]*)/#[metis_func]\rpub extern "C" fn \3 -> \2\1/
```


## 3. Fix bracketing and blocks
Helps to record macro for bracket-less `if` and `for` statements:
```
ys_B
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

With this, all syntax errors should be fixed. After you run `cargo fmt`, this
is the end of the first pass.

## 6. Handle local heap allocations

METIS handles local heap allocations (that are freed at end of scope) as follows:
```C
WSPACEPUSH;
var = iwspacemalloc(ctrl, cnt);
// ...
WSPACEPOP;
```
All workspace allocations will be freed at end of scope. I'm not quite sure if
this is the case for `imalloc` (wrapper for `gk_malloc`), but I don't quite
think that's the case. Regardless, we can typically replace these with
`Vecs<_>`.

Until I figure out for sure, **Don't use local `Vec`s for allocations made with
`gk_malloc`**

```
:%s/\viwspacemalloc\(ctrl, (.*)\);/vec![0; \1 as usize]/g
:%s/\vrwspacemalloc\(ctrl, (.*)\);/vec![0.0; \1 as usize]/g

:%s/\iset\((\w+), iwspacemalloc\(ctrl, (\w+)\)\);/vec![\1; \2 as usize]/g
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

To change the name (also works without field access):
```rust
mkslice!(cwhere: cgraph->where_, nvtxs);
// expands to approximately:
// let cwhere = std::slice::from_raw_parts((*cgraph).where_, {nvtxs} as usize);
```

## 8. Replacing function calls
first, we want to replace function calls to existing modules with their qualified syntax:
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
:%s/ASSERTP\=/assert!/I
:%s/\v\%(.{-0,5})"%(PR%(IDX|REAL))"/{:\1}/g
:%s/\v([^\w])printf/\1println!/g
:%s/\\n"/"/g
:g/\vgk_%(start|stop)cputimer/normal gcc
:%s/\vgk_errexit\(\w+, =/panic!(/
:%s/\vsizeof\((\w+)\)/std::mem::size_of::<\1>()/g

:%s/\vBND\w+/&!/I
:%s/INC_DEC/\L&!/I
:%s/%d/{}/g
:%s/\vSWAP\((.{-}),(.{-}),.{-}\)/std::mem::swap(\&mut \1, \&mut \2)/I
```

After all functions are replaced, we need to declare the remaining in
`bindings.rs`. Here's the best way I found to do it:

1. Copy the identifier for every not found function to the `#[metis_decl]`
   block in `bindings.rs`, one function per line.

2. Record a macro that:
    1. Yanks the identifier
    2. Opens `proto.h`
    3. Finds the yanked identifier
    4. Yanks the lines used by the arguments (`vibVy`)
    5. Returns to `bindings.rs` and pastes the yanked C declaration

3. Merge multiline declarations to a single line and fix it.

The last step can be accomplished with modified commands from step 2. Start
by highlighting the new declarations and then execute:
```
:'<,'>:g/\v\([^)]+$/join
:'<,'>:g/\v\([^)]+$/join
:'<,'>:g/\v\([^)]+$/join
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
:%s/\vMAKECSR\(.{-},(.{-}),(.{-})\)/util::make_csr(\1, \2)/I
:%s/\SHIFT\(.{-},(.{-}),(.{-})\)/util::shift_csr(\1, \2)/I
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
:s/\vicopy\(.*,(.*),(.*)\);/\/\/ &\r\2.copy_from_slice(&\1);/
```

## 9. Everything else

Still a lot to do get to it until it compiles. Every time errors get to zero
and jump back up is a pass.

Once we've committed a version that compiles, we can run `cargo fix` to help
with some simple style things like `if` statement parentheses.

## Things to look out for

- make sure that the `gk_malloc` calls have null-terminated strings

- we use a ton of casting to `usize`, but remember that `idx_t` is signed and
  often negative
