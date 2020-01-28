### Build-time user configurations
DEBUG ?= 0
PAUL ?= 1
CC ?= mpicc
CFLAGS ?= -O2

### Dependencies
GSLIBDIR ?= $(GSLIBPATH)
EXADIR ?= $(BUILDDIR)/3rd_party/exa/exaCore/build
EXASORTDIR ?= $(BUILDDIR)/3rd_party/exa/exaSort/build

### Meta info
SRCDIR=src
BUILDDIR=build
EXAMPLEDIR=examples

### Flags
INCFLAGS=-I$(SRCDIR) -I$(SRCDIR)/parCon -I$(GSLIBDIR)/include\
  -I$(EXADIR)/include -I$(EXASORTDIR)/include

LDFLAGS= -L$(GSLIBDIR)/lib -lgs -lm
EXALDFLAGS= -L$(EXADIR)/lib -lexa
EXASORTLDFLAGS= -L$(EXASORTDIR)/lib -lexaSort

PARRSBSRC=$(wildcard $(SRCDIR)/*.c)
SRCOBJS =$(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(PARRSBSRC))

PARCONSRC=$(wildcard $(SRCDIR)/parCon/*.c)
SRCOBJS +=$(patsubst $(SRCDIR)/parCon/%.c,$(BUILDDIR)/parCon/%.o,\
  $(PARCONSRC))

EXAMPLESRC=$(EXAMPLEDIR)/parRSB.c $(EXAMPLEDIR)/parCon.c
EXAMPLEOBJ=$(patsubst $(EXAMPLEDIR)/%.c,$(BUILDDIR)/examples/%,\
	$(EXAMPLESRC))
EXAMPLELDFLAGS=-L$(BUILDDIR)/lib -lparRSB $(EXASORTLDFLAGS)\
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

ifeq ($(GSLIBPATH),)
	$(error Specify GSLIBPATH=<path to gslib>/build)
endif

CFLAGS += -fPIC

.PHONY: default
default: deps lib install

.PHONY: all
all: deps lib examples install

.PHONY: deps
deps:
	@cp -rf 3rd_party $(BUILDDIR)/
	@cd $(BUILDDIR)/3rd_party/exa && GSDIR=$(GSLIBDIR) ./install

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

$(BUILDDIR)/parCon/%.o: $(SRCDIR)/parCon/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: lib
lib: deps $(SRCOBJS)
	$(CC) -shared -o $(LIB) $(SRCOBJS) $(EXASORTLDFLAGS) \
    $(EXALDFLAGS) $(LDFLAGS)

.PHONY: install
install: lib
	@mkdir -p $(INSTALL_ROOT)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALL_ROOT)/lib 2>/dev/null
	@mkdir -p $(INSTALL_ROOT)/include 2>/dev/null
	@cp $(SRCDIR)/parRSB.h $(INSTALL_ROOT)/include 2>/dev/null

.PHONY: examples
examples: lib install $(EXAMPLEOBJ)

$(BUILDDIR)/examples/%: $(EXAMPLEDIR)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(EXAMPLELDFLAGS)

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
