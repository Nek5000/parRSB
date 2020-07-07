DEBUG ?= 0
PAUL ?= 1
CC ?= mpicc
CFLAGS ?= -O2

MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT ?= $(patsubst %/,%,$(dir $(MKFILEPATH)))

GSLIBDIR=$(GSLIBPATH)

SRCDIR  =$(SRCROOT)/src
SORTDIR =$(SRCROOT)/src/sort
BUILDDIR=$(SRCROOT)/build
TESTDIR =$(SRCROOT)/example

TARGET=parRSB
TESTS=$(TESTDIR)/example
LIB=src/lib$(TARGET).a

INCFLAGS=-I$(SRCDIR) -I$(SORTDIR) -I$(GSLIBDIR)/include
LDFLAGS:=-L$(BUILDDIR)/lib -l$(TARGET) -L $(GSLIBDIR)/lib -lgs -lm

ifneq (,$(strip $(DESTDIR)))
  INSTALL_ROOT=$(DESTDIR)
else
  INSTALL_ROOT=$(SRCROOT)/build
endif

SRCS=$(wildcard $(SRCDIR)/*.c)
SORTSRCS=$(wildcard $(SORTDIR)/*.c)

SRCOBJS =$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SRCS))
SRCOBJS+=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SORTSRCS))

PP=

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(PAUL),0)
  PP += -DGENMAP_PAUL
endif

.PHONY: default
default: check lib install

.PHONY: all
all: check lib tests install

.PHONY: install
install: lib
	@mkdir -p $(INSTALL_ROOT)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALL_ROOT)/lib 2>/dev/null
	@mkdir -p $(INSTALL_ROOT)/include 2>/dev/null
	@cp $(SRCDIR)/parRSB.h $(INSTALL_ROOT)/include 2>/dev/null


.PHONY: $(TARGET)
lib: $(SRCOBJS)
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

.PHONY: check
check: 
ifeq ($(GSLIBPATH),)
	$(error Specify GSLIBPATH=<path to gslib>/build)
endif

$(BUILDDIR)/%.o: $(SRCROOT)/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: tests
tests: $(TESTS)

$(TESTS): lib install
	$(CC) $(CFLAGS) -I$(GSLIBDIR)/include -I$(BUILDDIR)/include $@.c -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -f $(BUILDDIR) $(LIB) $(TESTS) $(TESTS).o

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/src/sort)
