DEBUG ?= 0
CC ?= mpicc
CXX ?= mpic++
CXXFLAGS ?= -O2
CFLAGS ?= -O2
GPU ?= 1
PAUL ?= 1
PP=

ifndef GSLIBPATH
  $(error "Required variable GSLIBPATH not set.")
endif
GSLIBDIR =  $(GSLIBPATH)

SRCROOT  =  .
SRCDIR   =  $(SRCROOT)/src
TESTDIR  =  $(SRCROOT)/example
BUILDDIR ?= $(SRCROOT)/build

TARGET=parRSB
LIB=$(SRCDIR)/lib$(TARGET).a

INCFLAGS=-I$(SRCDIR) -I$(GSLIBDIR)/include
TESTINCFLAGS=-I$(GSLIBDIR)/include -I$(BUILDDIR)/include
TESTLDFLAGS:=-L$(BUILDDIR)/lib -l$(TARGET) -L $(GSLIBDIR)/lib -lgs -lm $(LDFLAGS)

CSRC:= \
  $(SRCDIR)/genmap.c \
  $(SRCDIR)/genmap-vector.c \
  $(SRCDIR)/genmap-handle.c \
  $(SRCDIR)/genmap-comm.c \
  $(SRCDIR)/genmap-eigen.c \
  $(SRCDIR)/genmap-laplacian.c \
  $(SRCDIR)/genmap-lanczos.c \
  $(SRCDIR)/genmap-rsb.c \
  $(SRCDIR)/genmap-chelpers.c \
  $(SRCDIR)/parRSB.c 

TESTS:= \
  $(TESTDIR)/example

ifeq ($(GPU),1)
  ifndef OCCA_DIR
    $(error "Required ENV-variable OCCA_DIR is not set.")
  endif
  ifndef PARRSB_DIR
    $(error "Required variable PARRSB_DIR not set.")
  endif

  CFLAGS += -DPARRSB_GPU -DPARRSB_OKL_DIR="\"$(PARRSB_DIR)/okl/\""
  INCFLAGS += -I$(OCCA_DIR)/include
  CSRC += $(SRCDIR)/genmap-occa.c

  TESTINCFLAGS += -I$(OCCA_DIR)
  TESTLDFLAGS += -Wl,-rpath,$(OCCA_DIR)/lib -L$(OCCA_DIR)/lib -locca

  CXXSRC =
endif

COBJS:=$(CSRC:.c=.c.o)
OBJS:=$(COBJS)

ifeq ($(GPU),1)
  CXXOBJS:=$(CXXSRC:.cpp=.cpp.o)
  OBJS += $(CXXOBJS)
endif

ifneq (,$(strip $(DESTDIR)))
  INSTALL_ROOT = $(DESTDIR)
else
  INSTALL_ROOT = $(BUILDDIR)
endif

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

.PHONY: lib
lib: $(OBJS)
	@$(AR) cr $(LIB) $(OBJS)
	@ranlib $(LIB)

.PHONY: check
check: 
ifeq ($(GSLIBPATH),)
	$(error Specify GSLIBPATH=<path to gslib>/build)
endif

$(COBJS): %.c.o: %.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@ $(LDFLAGS)

$(CXXOBJS): %.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PP) $(INCFLAGS) -c $< -o $@ $(LDFLAGS)

.PHONY: tests
tests: $(TESTS)

$(TESTS): lib install
	$(CC) $(CFLAGS) $(TESTINCFLAGS) $@.c -o $@ $(TESTLDFLAGS)

.PHONY: clean
clean:
	@rm -f $(SRCOBJS) $(LIB) $(TESTS) $(TESTS).o

.PHONY: astyle
astyle:
	astyle --style=google --indent=spaces=2 --max-code-length=80 \
	    --keep-one-line-statements --keep-one-line-blocks --lineend=linux \
            --suffix=none --preserve-date --formatted --pad-oper \
	    --unpad-paren example/*.[ch] src/*.[ch]

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true
