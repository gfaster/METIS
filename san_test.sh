#!/bin/sh

ulimit -c unlimited
export ASAN_OPTIONS=abort_on_error=1:disable_coredump=0:unmap_shadow_on_exit=1:detect_leaks=0
export RUSTFLAGS=-Zsanitizer=address
cargo +nightly t --lib -Zbuild-std --target x86_64-unknown-linux-gnu -j1 -q
