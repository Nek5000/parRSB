CC ?= mpicc
CFLAGS ?= -Wall -Wextra -Wpedantic -Wno-unused-function -Wno-unused-parameter -std=c99 -g
LDFLAGS ?=
DEBUG ?= 0
UNDERSCORE ?= 1
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASDIR ?=
BLASFLAGS ?= -lblas -llapack
GSLIBPATH ?=

########################## Don't touch what follows ###########################
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib build>)
endif

MKFILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT := $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRCDIR = $(SRCROOT)/src
EXAMPLEDIR = $(SRCROOT)/examples
BUILDROOT = $(SRCROOT)/build
ifneq (,$(strip $(DESTDIR)))
INSTALLROOT = $(DESTDIR)
else
INSTALLROOT = $(SRCROOT)/install
endif

SRCS = $(wildcard $(SRCDIR)/*.c)
SRCOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDROOT)/%.o,$(SRCS))
EXAMPLES = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLEBINS = $(patsubst $(SRCROOT)/%.c,$(BUILDROOT)/%,$(EXAMPLES))

LIB = $(BUILDROOT)/lib/libparRSB.a

ifneq ($(DEBUG),0)
  PP += -DPARRSB_DEBUG
  CFLAGS += -g
else
  CFLAGS += -O2
endif

PP += -DPARRSB_MPI

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

INCFLAGS = -I$(SRCDIR) -I$(GSLIBPATH)/include
CCCMD = $(CC) $(CFLAGS) $(INCFLAGS) $(PP)
LDFLAGS += -lm

.PHONY: all lib install examples format clean

all: lib install examples

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
	$(CCCMD) $< -o $@ -L$(INSTALLROOT)/lib -lparRSB -L$(GSLIBPATH)/lib -lgs \
		$(LDFLAGS)

$(shell mkdir -p $(BUILDROOT)/examples)
$(shell mkdir -p $(BUILDROOT)/src)
