### Build-time user configurations ###
DEBUG ?= 0
PAUL ?= 1
LANCZOS ?= 0
CC ?= mpicc
CFLAGS ?= -O2
BUILDDIR ?= $(CURDIR)/build
DESTDIR ?=

# Dependencies #
GS_DIR = $(GSLIBPATH)
EXA_DIR ?= $(PREFIX)
EXASORT_DIR ?= $(PREFIX)

### Meta info ###
SRCROOT=$(CURDIR)
SRCDIR=$(SRCROOT)/src
EXAMPLEDIR=$(SRCROOT)/examples
TESTDIR=$(SRCROOT)/tests

### Flags ###
INCFLAGS=-I$(SRCDIR) -I$(SRCDIR)/parCon -I$(GS_DIR)/include\
  -I$(EXA_DIR)/include -I$(EXASORT_DIR)/include
LDFLAGS=-Wl,-rpath,$(GS_DIR) -L$(GS_DIR)/lib -lgs -lm

EXALDFLAGS=-Wl,-rpath,$(EXA_DIR)/lib -L$(EXA_DIR)/lib -lexa
EXASORTLDFLAGS=-Wl,-rpath,$(EXASORT_DIR)/lib\
  -L$(EXASORT_DIR)/lib -lexaSort

PARRSBSRC=$(wildcard $(SRCDIR)/*.c)
SRCOBJS =$(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(PARRSBSRC))

PARCONSRC=$(wildcard $(SRCDIR)/parCon/*.c)
SRCOBJS +=$(patsubst $(SRCDIR)/parCon/%.c,\
  $(BUILDDIR)/parCon/%.o,$(PARCONSRC))

EXAMPLESRC=$(EXAMPLEDIR)/parRSB.c $(EXAMPLEDIR)/parCon.c\
  $(EXAMPLEDIR)/parRCB.c
EXAMPLEOBJ=$(patsubst $(EXAMPLEDIR)/%.c,$(BUILDDIR)/examples/%,\
	$(EXAMPLESRC))
EXAMPLELDFLAGS= -Wl,-rpath,$(BUILDDIR) -L$(BUILDDIR) -lparRSB\
  $(EXASORTLDFLAGS) $(EXALDFLAGS) $(LDFLAGS)

LIB=$(BUILDDIR)/libparRSB.so

PP=

ifneq (,$(strip $(DESTDIR)))
  PREFIX=$(realpath $(DESTDIR))
else
  PREFIX=$(realpath ./build)
endif

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(PAUL),0)
  PP += -DGENMAP_PAUL
endif

ifneq ($(LANCZOS),0)
	PP += -DGENMAP_LANCZOS
endif

ifeq ($(GS_DIR),)
	$(error Specify GSLIBPATH=<path to gslib>/build)
endif

CFLAGS += -fPIC

.PHONY: default
default: build-deps lib install

.PHONY: all
all: build-deps lib examples tests install

.PHONY: build-deps
build-deps:
	@cp 3rd_party/exa.install $(BUILDDIR)/
	@cd $(BUILDDIR) && GS_DIR=$(GS_DIR)\
		PREFIX=$(PREFIX) ./exa.install

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

$(BUILDDIR)/parCon/%.o: $(SRCDIR)/parCon/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: lib
lib: build-deps $(SRCOBJS)
	$(CC) -shared -o $(LIB) $(SRCOBJS) $(EXASORTLDFLAGS) \
    $(EXALDFLAGS) $(LDFLAGS)

.PHONY: install
install: lib
	@mkdir -p $(PREFIX)/lib 2>/dev/null
	@cp $(LIB) $(PREFIX)/lib/ 2>/dev/null
	@mkdir -p $(PREFIX)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SRCDIR)/parCon/*.h \
		$(PREFIX)/include/ 2>/dev/null

.PHONY: examples
examples: lib $(EXAMPLEOBJ)

$(BUILDDIR)/examples/%: $(EXAMPLEDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(EXAMPLELDFLAGS)

.PHONY: tests
tests: examples
	@cp -rf $(TESTDIR) $(BUILDDIR)/
	@cd $(BUILDDIR)/tests && BUILDDIR=$(BUILDDIR) ./test.sh --run

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

$(shell mkdir -p $(BUILDDIR))
$(shell mkdir -p $(BUILDDIR)/parCon)
$(shell mkdir -p $(BUILDDIR)/examples)
