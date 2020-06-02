### Build-time user configurations ###
DEBUG ?= 1
TQLI ?= 1
LANCZOS ?= 0
CC ?= mpicc
CFLAGS ?= -O0
BUILDDIR ?= $(CURDIR)/build
DESTDIR ?=
MPI ?= 1

### Dependencies ###
GS_DIR = $(GSLIBPATH)

##### Don't touch what follows #####
### Some pre-processing ###
ifneq ($(strip $(DESTDIR)),)
  PREFIX=$(abspath $(DESTDIR))
else
  PREFIX=$(abspath ./build)
endif

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(TQLI),0)
  PP += -DGENMAP_TQLI
endif

ifneq ($(MPI),0)
  CFLAGS += -DEXA_MPI
endif

ifneq ($(LANCZOS),0)
	PP += -DGENMAP_LANCZOS
endif

ifeq ($(GS_DIR),)
	$(error Specify GSLIBPATH=<path to gslib>/build)
endif

### Meta info ###
SRCROOT=$(CURDIR)
SRCDIR=$(SRCROOT)/src
EXAMPLEDIR=$(SRCROOT)/examples
TESTDIR=$(SRCROOT)/tests

EXA_DIR=$(PREFIX)
EXASORT_DIR=$(PREFIX)

LIB=$(BUILDDIR)/libparRSB.so
PP=

### Flags ###
INCFLAGS=-I$(SRCDIR) -I$(EXAMPLEDIR) -I$(SRCDIR)/parCon \
  -I$(GS_DIR)/include -I$(EXA_DIR)/include -I$(EXASORT_DIR)/include
CFLAGS += -fPIC

LDFLAGS=-Wl,-rpath,$(GS_DIR) -L$(GS_DIR)/lib -lgs -lm
EXALDFLAGS=-Wl,-rpath,$(EXA_DIR)/lib -L$(EXA_DIR)/lib -lexa
EXASORTLDFLAGS=-Wl,-rpath,$(EXASORT_DIR)/lib\
  -L$(EXASORT_DIR)/lib -lexaSort

PARRSBSRC=$(wildcard $(SRCDIR)/*.c)
SRCOBJS =$(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/src/%.o,$(PARRSBSRC))

PARCONSRC=$(wildcard $(SRCDIR)/parCon/*.c)
SRCOBJS +=$(patsubst $(SRCDIR)/parCon/%.c,\
  $(BUILDDIR)/src/parCon/%.o,$(PARCONSRC))

EXAMPLESRC=$(EXAMPLEDIR)/parRSB.c $(EXAMPLEDIR)/parCon.c\
  $(EXAMPLEDIR)/parRCB.c
EXAMPLEOBJ=$(patsubst $(EXAMPLEDIR)/%.c,$(BUILDDIR)/examples/%,\
	$(EXAMPLESRC))
EXAMPLELDFLAGS= -Wl,-rpath,$(BUILDDIR) -L$(BUILDDIR) -lparRSB\
  $(EXASORTLDFLAGS) $(EXALDFLAGS) $(LDFLAGS)

TESTSRC=$(wildcard $(TESTDIR)/*.c)
TESTOBJ=$(patsubst $(TESTDIR)/%.c,$(BUILDDIR)/tests/%,$(TESTSRC))
TESTLDFLAGS= -Wl,-rpath,$(BUILDDIR) -L$(BUILDDIR) -lparRSB\
  $(EXASORTLDFLAGS) $(EXALDFLAGS) $(LDFLAGS)

.PHONY: default
default: build-deps lib install

.PHONY: all
all: build-deps lib examples tests install

.PHONY: build-deps
build-deps:
	@cd $(BUILDDIR) && MPI=$(MPI) GS_DIR=$(GS_DIR) PREFIX=$(PREFIX) \
		$(SRCROOT)/get_deps.sh

$(BUILDDIR)/src/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

$(BUILDDIR)/src/parCon/%.o: $(SRCDIR)/parCon/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: lib
lib: build-deps $(SRCOBJS)
	$(CC) -shared -o $(LIB) $(SRCOBJS) $(EXASORTLDFLAGS) $(EXALDFLAGS) \
		$(LDFLAGS)

.PHONY: install
install: lib
	@mkdir -p $(PREFIX)/lib 2>/dev/null
	@cp $(LIB) $(PREFIX)/lib/ 2>/dev/null
	@mkdir -p $(PREFIX)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SRCDIR)/parCon/*.h $(PREFIX)/include/ 2>/dev/null

.PHONY: examples
examples: lib $(EXAMPLEOBJ)

$(BUILDDIR)/examples/%: $(EXAMPLEDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(EXAMPLELDFLAGS)

.PHONY: tests
tests: lib $(TESTOBJ)

$(BUILDDIR)/tests/%: $(TESTDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(TESTLDFLAGS)

.PHONY: clean
clean:
	@rm -rf $(BUILDDIR)

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/src)
$(shell mkdir -p $(BUILDDIR)/src/parCon)
$(shell mkdir -p $(BUILDDIR)/examples)
$(shell mkdir -p $(BUILDDIR)/tests)
