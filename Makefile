### Build-time user configurations
DEBUG ?= 0
PAUL ?= 1
CC ?= mpicc
CFLAGS ?= -O2

### Dependencies
GSDIR = $(GSLIBPATH)

### Meta info
SRCROOT=${CURDIR}
SRCDIR=$(SRCROOT)/src
EXAMPLEDIR=$(SRCROOT)/examples
TESTDIR=$(SRCROOT)/tests

BUILDDIR ?=$(SRCROOT)/build
EXADIR ?= $(BUILDDIR)/exaCore/build
EXASORTDIR ?= $(BUILDDIR)/exaSort/build

### Flags
INCFLAGS=-I$(SRCDIR) -I$(SRCDIR)/parCon -I$(GSDIR)/include\
  -I$(EXADIR)/include -I$(EXASORTDIR)/include

LDFLAGS= -L$(GSDIR)/lib -lgs -lm
EXALDFLAGS= -L$(EXADIR)/lib -lexa
EXASORTLDFLAGS= -L$(EXASORTDIR)/lib -lexaSort

PARRSBSRC=$(wildcard $(SRCDIR)/*.c)
SRCOBJS =$(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(PARRSBSRC))

PARCONSRC=$(wildcard $(SRCDIR)/parCon/*.c)
SRCOBJS +=$(patsubst $(SRCDIR)/parCon/%.c,\
  $(BUILDDIR)/parCon/%.o,$(PARCONSRC))

EXAMPLESRC=$(EXAMPLEDIR)/parRSB.c $(EXAMPLEDIR)/parCon.c
EXAMPLEOBJ=$(patsubst $(EXAMPLEDIR)/%.c,$(BUILDDIR)/examples/%,\
	$(EXAMPLESRC))
EXAMPLELDFLAGS=-L$(BUILDDIR) -lparRSB $(EXASORTLDFLAGS)\
	$(EXALDFLAGS) $(LDFLAGS)

LIB=$(BUILDDIR)/libparRSB.so

PP=

ifneq (,$(strip $(DESTDIR)))
  INSTALL_ROOT = $(DESTDIR)
else
  INSTALL_ROOT = build
endif

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(PAUL),0)
  PP += -DGENMAP_PAUL
endif

ifeq ($(GSDIR),)
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
	@cd $(BUILDDIR) && GSDIR=$(GSDIR) ./exa.install

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
	@mkdir -p $(INSTALL_ROOT)/lib 2>/dev/null
	@cp $(LIB) $(INSTALL_ROOT)/lib/ 2>/dev/null
	@mkdir -p $(INSTALL_ROOT)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SRCDIR)/parCon/*.h \
		$(INSTALL_ROOT)/include/ 2>/dev/null

.PHONY: examples
examples: lib $(EXAMPLEOBJ)

$(BUILDDIR)/examples/%: $(EXAMPLEDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(EXAMPLELDFLAGS)

.PHONY: tests
tests: examples
	@cp -rf $(TESTDIR) $(BUILDDIR)/
	@cd $(BUILDDIR)/tests && ./test.sh --get-deps
	@cd $(BUILDDIR)/tests && EXADIR=$(EXADIR)\
		EXASORTDIR=$(EXASORTDIR) BUILDDIR=$(BUILDDIR)\
		./test.sh --run

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
