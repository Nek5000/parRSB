## General build parameters ##
CC ?= mpicc
CFLAGS ?= -g -O2
DEBUG ?= 1
MPI ?= 1
UNDERSCORE ?= 1
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASDIR ?=
BLASFLAGS ?= -lblas -llapack
OCCA ?= 0
OCCADIR ?=

## Don't touch what follows ##
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib>/build)
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

SRCS  = $(wildcard $(SRCDIR)/genmap*.c)
SRCS += $(wildcard $(SRCDIR)/parrsb*.c)
SRCS += $(wildcard $(SRCDIR)/sort/*.c)
SRCS += $(wildcard $(SRCDIR)/precond/*.c)
SRCS += $(wildcard $(SRCDIR)/gencon/*.c)

EXAMPLES = $(wildcard $(EXAMPLEDIR)/*.c)

INCFLAGS = -I$(SRCDIR) -I$(SRCDIR)/sort -I$(SRCDIR)/precond -I$(SRCDIR)/gencon \
	-I$(GSLIBPATH)/include
LDFLAGS = -L$(BUILDDIR)/lib -lparRSB -L$(GSLIBPATH)/lib -lgs -lm

LIB = $(BUILDDIR)/lib/libparRSB.a

PP =

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
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
  SRCS += $(wildcard $(SRCDIR)/occa*.c)
endif

INSTALLDIR=
ifneq (,$(strip $(DESTDIR)))
  INSTALLDIR = $(realpath $(DESTDIR))
endif

SRCOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SRCS))
EXAMPLEOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%,$(EXAMPLES))

.PHONY: all
all: lib examples install

.PHONY: install
install: lib
ifneq ($(INSTALLDIR),)
	@mkdir -p $(INSTALLDIR)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLDIR)/lib 2>/dev/null
	@mkdir -p $(INSTALLDIR)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SRCDIR)/sort/*.h $(SRCDIR)/precond/*.h \
		$(INSTALLDIR)/include 2>/dev/null
	@cp -r okl $(INSTALLDIR)/ 2>/dev/null
endif

.PHONY: lib
lib: $(SRCOBJS)
	@mkdir -p $(BUILDDIR)/lib
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

$(BUILDDIR)/src/%.o: $(SRCROOT)/src/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: examples
examples: install $(EXAMPLEOBJS)

$(BUILDDIR)/examples/%: $(SRCROOT)/examples/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -rf $(BUILDDIR)

.PHONY: format
format:
	find . -iname *.h -o -iname *.c | xargs clang-format -i

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/src/sort)
$(shell mkdir -p $(BUILDDIR)/src/precond)
$(shell mkdir -p $(BUILDDIR)/src/gencon)
$(shell mkdir -p $(BUILDDIR)/examples)
