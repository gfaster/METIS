# Configuration options.
i64        = not-set
r64        = not-set
gdb        = not-set
assert     = not-set
assert2    = not-set
debug      = not-set
gprof      = not-set
valgrind   = not-set
openmp     = not-set
shared     = not-set
cc         = not-set
prefix     = ./install
gklib_path = ./GKlib/build/install


# Basically proxies everything to the builddir cmake.

BUILDDIR = build

IDXWIDTH  = "\#define IDXTYPEWIDTH 32"
REALWIDTH = "\#define REALTYPEWIDTH 32"

# Process configuration options.
CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=1
ifneq ($(gklib_path), not-set)
    CONFIG_FLAGS += -DGKLIB_PATH=$(abspath $(gklib_path))
endif
ifneq ($(prefix), not-set)
    CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX=$(prefix)
endif
ifneq ($(i64), not-set)
    IDXWIDTH  = "\#define IDXTYPEWIDTH 64"
endif
ifneq ($(r64), not-set)
    REALWIDTH = "\#define REALTYPEWIDTH 64"
endif
ifneq ($(gdb), not-set)
    CONFIG_FLAGS += -DGDB=$(gdb)
endif
ifneq ($(assert), not-set)
    CONFIG_FLAGS += -DASSERT=$(assert)
endif
ifneq ($(assert2), not-set)
    CONFIG_FLAGS += -DASSERT2=$(assert2)
endif
ifneq ($(debug), not-set)
    CONFIG_FLAGS += -DDEBUG=$(debug)
endif
ifneq ($(gprof), not-set)
    CONFIG_FLAGS += -DGPROF=$(gprof)
endif
ifneq ($(valgrind), not-set)
    CONFIG_FLAGS += -DVALGRIND=$(valgrind)
endif
ifneq ($(openmp), not-set)
    CONFIG_FLAGS += -DOPENMP=$(openmp)
endif
ifneq ($(shared), not-set)
    CONFIG_FLAGS += -DSHARED=1
endif
ifneq ($(cc), not-set)
    CONFIG_FLAGS += -DCMAKE_C_COMPILER=$(cc)
endif

RUSTSRC = $(wildcard libmetis/*.rs)
RUSTOBJ = $(addprefix $(BUILDDIR)/,$(notdir $(RUSTSRC:.rs=.o)))
CSRC = $(wildcard libmetis/*.c)
COBJ = $(addprefix $(BUILDDIR)/,$(notdir $(CSRC:.c=.o)))
ALLOBJ = $(COBJ) $(RUSTOBJ)

VERNUM=5.1.0
PKGNAME=metis-$(VERNUM)

define run-config
mkdir -p $(BUILDDIR)
mkdir -p $(BUILDDIR)/xinclude
echo $(IDXWIDTH) > $(BUILDDIR)/xinclude/metis.h
echo $(REALWIDTH) >> $(BUILDDIR)/xinclude/metis.h
cat include/metis.h >> $(BUILDDIR)/xinclude/metis.h
cp include/CMakeLists.txt $(BUILDDIR)/xinclude
cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

lib: $(BUILDDIR)/libmetis.a
	@echo "complete"

$(BUILDDIR)/libmacros.so: rust/macros/target/debug/libmacros.so
	cp $< $@

rust/macros/target/debug/libmacros.so: $(wildcard rust/macros/src/*) rust/macros/Cargo.toml rust/macros/Cargo.lock
	cd rust/macros; cargo build

$(BUILDDIR)/libbindings.rlib: rust/bindings/bindings.rs $(BUILDDIR)/libmacros.so
	rustc -g --crate-type lib \
	    --crate-name bindings \
	    --extern macros=$(BUILDDIR)/libmacros.so \
	    --out-dir $(BUILDDIR) \
	    $<

$(RUSTOBJ):  $(BUILDDIR)/%.o: libmetis/%.rs $(BUILDDIR)/libmacros.so $(BUILDDIR)/libbindings.rlib
	rustc -g --emit obj $< \
	    --crate-type lib \
	    --extern macros=$(BUILDDIR)/libmacros.so \
	    --extern bindings=$(BUILDDIR)/libbindings.rlib \
	    --out-dir $(BUILDDIR)

$(COBJ): $(BUILDDIR)/%.o: libmetis/%.c $(BUILDDIR)/include/metis.h
	@echo "$^"
	gcc -g -c -o $@ -I$(gklib_path)/include -I$(BUILDDIR)/include $<

$(BUILDDIR)/libmetis.a: $(ALLOBJ)
	ar -r $@ $?

$(BUILDDIR)/include/metis.h: include/metis.h
	mkdir -p $(BUILDDIR)/include
	echo $(IDXWIDTH) > $(BUILDDIR)/include/metis.h
	echo $(REALWIDTH) >> $(BUILDDIR)/include/metis.h
	cat include/metis.h >> $(BUILDDIR)/include/metis.h

clean:
	rm -r $(BUILDDIR)


.PHONY: config distclean all clean install uninstall remake dist rustbuild
