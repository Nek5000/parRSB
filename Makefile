## General build parameters ##
CC ?= mpicc
CFLAGS ?= -g
DEBUG ?= 0
MPI ?= 1
UNDERSCORE ?= 1
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASDIR ?=
BLASFLAGS ?= -lblas -llapack
OCCA ?= 0
OCCADIR ?=

########################## Don't touch what follows ###########################
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib build>)
endif

ifneq ($(OCCA), 0)
  ifeq ($(OCCADIR),)
    $(error Specify OCCADIR=<path to occa> when using OCCA=1)
  endif
endif

MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT_ ?= $(patsubst %/,%,$(dir $(MKFILEPATH)))
SRCROOT = $(realpath $(SRCROOT_))

SRCDIR = $(SRCROOT)/src
BUILDDIR = $(SRCROOT)/build
EXAMPLEDIR = $(SRCROOT)/examples
INSTALLDIR=$(SRCROOT)

INCFLAGS = -I$(SRCDIR) -I$(GSLIBPATH)/include
LDFLAGS = -L$(BUILDDIR)/lib -lparRSB -L$(GSLIBPATH)/lib -lgs -lm

SRCS  = $(wildcard $(SRCDIR)/*.c)
EXAMPLES = $(wildcard $(EXAMPLEDIR)/*.c)
LIB = $(BUILDDIR)/lib/libparRSB.a

SRCOBJS = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))
EXAMPLEOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%,$(EXAMPLES))

PP =

ifneq ($(DEBUG),0)
  PP += -DGENMAP_DEBUG
else
  CFLAGS += -O2
endif

ifneq ($(MPI),0)
  PP += -DMPI
endif

ifneq ($(UNDERSCORE),0)
  PP += -DGENMAP_UNDERSCORE
endif

ifneq ($(SYNC_BY_REDUCTION),0)
  PP += -DGENMAP_SYNC_BY_REDUCTION
endif

ifneq ($(BLAS),0)
  PP += -DGENMAP_BLAS
  ifneq ($(BLASDIR),)
    LDFLAGS+= -L$(BLASDIR)
  endif
  LDFLAGS += $(BLASFLAGS)
endif

ifneq ($(OCCA),0)
  PP += -DGENMAP_OCCA
  LDFLAGS+= -L$(OCCADIR)/lib -locca
  INCFLAGS+= -I$(OCCADIR)/include
  SRCS += $(wildcard $(SRCDIR)/occa/*.c)
endif

ifneq (,$(strip $(DESTDIR)))
  INSTALLDIR = $(realpath $(DESTDIR))
endif

.PHONY: all install lib examples clean format

all: examples

install: lib
	@mkdir -p $(INSTALLDIR)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLDIR)/lib 2>/dev/null
	@mkdir -p $(INSTALLDIR)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(INSTALLDIR)/include 2>/dev/null
	@cp -r $(SRCDIR)/occa/*.okl $(INSTALLDIR)/ 2>/dev/null

lib: $(SRCOBJS)
	@mkdir -p $(BUILDDIR)/lib
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

$(BUILDDIR)/%.o: $(SRCROOT)/src/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

examples: install $(EXAMPLEOBJS)

$(BUILDDIR)/examples/%: $(SRCROOT)/examples/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(LDFLAGS)

clean:
	@rm -rf $(BUILDDIR)

format:
	find . -iname *.h -o -iname *.c -o -iname *.okl | xargs clang-format -i

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/gencon)
$(shell mkdir -p $(BUILDDIR)/occa)
$(shell mkdir -p $(BUILDDIR)/examples)
