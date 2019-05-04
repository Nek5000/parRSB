DEBUG ?= 0
CC ?= mpicc
CXX ?= mpic++
CXXFLAGS ?= -O2
CFLAGS ?= -O2
GPU ?= 0
PAUL ?= 1
PP=

## parRSB configs
SRCROOT  = .
SRCDIR   = $(SRCROOT)/src
TESTDIR  = $(SRCROOT)/example
BUILDDIR = $(SRCROOT)/build

## pre-build checks
ifndef GSLIBPATH
  $(error "Required variable GSLIBPATH not set.")
endif
GSLIBDIR =  $(GSLIBPATH)

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

## Build parRSB and targets
TARGET=parRSB
LIB=$(SRCDIR)/lib$(TARGET).a
INCFLAGS=-I$(SRCDIR) -I$(GSLIBDIR)/include

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
COBJ:=$(CSRC:.c=.c.o)

TESTSRC:= \
  $(TESTDIR)/example.c
TESTOBJ:=$(TESTSRC:.c=.c.o)
TESTEXE:=$(TESTSRC:.c=.out)

TESTINCFLAGS=-I$(GSLIBDIR)/include -I$(INSTALL_ROOT)/include
TESTLDFLAGS +=-L$(INSTALL_ROOT)/lib -l$(TARGET)

ifneq ($(GPU),0)
  ifndef OCCA_DIR
    $(error "Required ENV-variable OCCA_DIR is not set.")
  endif
  ifndef PARRSB_DIR
    $(error "Required variable PARRSB_DIR not set.")
  endif
  ifndef OGS_DIR
    $(error "Required variable OGS_DIR not set.")
  endif

  CXXSRC += $(SRCDIR)/parRSB-occa.cpp
  CXXOBJ:=$(CXXSRC:.cpp=.cpp.o)

  CXXFLAGS += -DPARRSB_GPU -DPARRSB_OKL_DIR="\"$(PARRSB_DIR)/okl/\""
  CFLAGS += -DPARRSB_GPU -DPARRSB_OKL_DIR="\"$(PARRSB_DIR)/okl/\""
  INCFLAGS += -I$(OCCA_DIR)/include -I$(OGS_DIR) -I$(OGS_DIR)/include

  TESTINCFLAGS += -I$(INSTALL_ROOT)/include
  ## Following order is important
  TESTLDFLAGS += -Wl,-rpath,$(OGS_DIR) -L$(OGS_DIR) -logs
  TESTLDFLAGS += -Wl,-rpath,$(OCCA_DIR)/lib -L$(OCCA_DIR)/lib -locca
endif

TESTLDFLAGS +=-L $(GSLIBDIR)/lib -lgs -lm

OBJ := $(COBJ) $(CXXOBJ)

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
lib: $(OBJ)
	@$(AR) cr $(LIB) $(OBJ)
	@ranlib $(LIB)

$(COBJ): %.c.o : %.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

$(CXXOBJ): %.cpp.o : %.cpp
	$(CXX) $(CXXFLAGS) $(PP) $(INCFLAGS) -c $< -o $@ 

.PHONY: tests
tests: lib install $(TESTEXE)

$(TESTEXE): $(TESTOBJ)
	$(CXX) $< -o $@ $(TESTLDFLAGS) $(LDFLAGS)

$(TESTOBJ): $(TESTSRC)
	$(CC) $(CFLAGS) $(TESTINCFLAGS) -c $< -o $@

.PHONY: check
check: 
ifeq ($(GSLIBPATH),)
	$(error "Specify GSLIBPATH=<path to gslib>/build")
endif
ifneq ($(GPU),0)
  ifndef OCCA_DIR
   	$(error "Required ENV-variable OCCA_DIR is not set.")
  endif
  ifndef PARRSB_DIR
   	$(error "Required variable PARRSB_DIR not set.")
  endif
  ifndef OGS_DIR
   	$(error "Required variable OGS_DIR not set.")
  endif
endif

.PHONY: clean
clean:
	@rm -f $(SRCOBJS) $(LIB) $(TESTS) $(TESTS).o

.PHONY: astyle
astyle:
	astyle --style=google --indent=spaces=2 --max-code-length=80 \
	    --keep-one-line-statements --keep-one-line-blocks --lineend=linux \
            --suffix=none --preserve-date --formatted --pad-oper \
	    --unpad-paren example/*.[ch] src/*.[ch] src/*.[ch]pp

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true
