CC ?= mpicc
CFLAGS ?=
DEBUG ?= 0
MPI ?= 1
UNDERSCORE ?= 1
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASDIR ?=
BLASFLAGS ?= -lblas -llapack

########################## Don't touch what follows ###########################
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib build>)
endif

MKFILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT := $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRCDIR = $(SRCROOT)/src
EXAMPLEDIR = $(SRCROOT)/examples
BUILDROOT = $(SRCROOT)/build
INSTALLROOT = $(BUILDROOT)/install
ifneq ($(strip $(DESTDIR)),)
  INSTALLROOT = $(realpath $(DESTDIR))
endif

SRCS = $(wildcard $(SRCDIR)/*.c)
SRCOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDROOT)/%.o,$(SRCS))
EXAMPLES = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLEBINS = $(patsubst $(SRCROOT)/%.c,$(BUILDROOT)/%,$(EXAMPLES))

ifneq ($(DEBUG),0)
  PP += -DPARRSB_DEBUG -g
else
  CFLAGS += -O2
endif

ifneq ($(MPI),0)
  PP += -DMPI
endif

ifneq ($(UNDERSCORE),0)
  PP += -DPARRSB_UNDERSCORE
endif

ifneq ($(SYNC_BY_REDUCTION),0)
  PP += -DPARRSB_SYNC_BY_REDUCTION
endif

ifneq ($(BLAS),0)
  PP += -DPARRSB_BLAS
  ifneq ($(BLASDIR),)
    LDFLAGS+= -L$(BLASDIR)
  endif
  LDFLAGS += $(BLASFLAGS)
endif

LIB = $(BUILDROOT)/lib/libparRSB.a
INCFLAGS = -I$(SRCDIR) -I$(GSLIBPATH)/include
LDFLAGS = -L$(INSTALLROOT)/lib -lparRSB -L$(GSLIBPATH)/lib -lgs -lm
CCCMD = $(CC) $(CFLAGS) $(INCFLAGS) $(PP)

.PHONY: all install lib examples clean format

all: lib install

lib: $(SRCOBJS)
	@mkdir -p $(BUILDROOT)/lib
	@$(AR) cr $(LIB) $?
	@ranlib $(LIB)

install: lib
	@mkdir -p $(INSTALLROOT)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLROOT)/lib 2>/dev/null
	@mkdir -p $(INSTALLROOT)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(INSTALLROOT)/include 2>/dev/null

examples: lib install $(EXAMPLEBINS)

format:
	find . -iname *.h -o -iname *.c -o -iname *.okl | xargs clang-format -i

clean:
	@$(RM) -rf $(BUILDROOT)

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(BUILDROOT)/%.o: $(SRCROOT)/%.c
	$(CCCMD) -c $< -o $@

$(BUILDROOT)/%: $(SRCROOT)/%.c | lib install
	$(CCCMD) $< -o $@ $(LDFLAGS)

$(shell mkdir -p $(BUILDROOT)/examples)
$(shell mkdir -p $(BUILDROOT)/src)
