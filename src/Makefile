# makefile for FEAST
# Author: Adam Pocock, apocock@cs.man.ac.uk
# Created: 29/06/2011
#
# Part of the FEAture Selection Toolbox, please reference
# "Conditional Likelihood Maximisation: A Unifying Framework for Information
# Theoretic Feature Selection"
# G. Brown, A. Pocock, M.-J. Zhao, M. Lujan
# Journal of Machine Learning Research (JMLR), 2012
#
# Please check www.cs.manchester.ac.uk/~gbrown/fstoolbox for updates.
#
# Copyright (c) 2010-2013, A. Pocock, G. Brown, The University of Manchester
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#   - Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#   - Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#   - Neither the name of The University of Manchester nor the names of its
#     contributors may be used to endorse or promote products derived from this
#     software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

PREFIX = /usr/local
CXXFLAGS = -O3 -fPIC -DCOMPILE_C
COMPILER = gcc
LINKER = gcc
ROOT_DIR = $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
MAKE_RUNABLE = 
OS = $(shell gcc -dumpmachine)
$(info $$OS is [${OS}])
ifneq (, $(findstring linux, $(OS)))
    #Linux
    MITOOLBOXPATH = $(ROOT_DIR)/MIToolbox/
    MAKE_RUNABLE = chmod 755 test
    $(info ************  LINUX ************)
else ifneq  (, $(findstring darwin, $(OS)))
    #Mac
    MITOOLBOXPATH = $(ROOT_DIR)/MIToolbox/
    MAKE_RUNABLE = chmod 755 test
    $(info ************  MAC ************)
else ifneq  (, $(findstring mingw, $(OS)))
    #mingw
    ROOT_DIR =.
    MITOOLBOXPATH = ./MIToolbox/
    RM =del
    $(info ************  MINGW ************)
else
    MITOOLBOXPATH = $(ROOT_DIR)\\MIToolbox\\
    $(info ************  WINDOWS ************)
endif
$(info Rootdir  $(ROOT_DIR))
objects = mRMR_D.o CMIM.o JMI.o DISR.o CondMI.o ICAP.o BetaGamma.o MIM.o FSToolboxR.o MutualInformation.o Entropy.o CalculateProbability.o ArrayOperations.o

LIB_A = libFSToolboxR.a
LIB_SO = libFSToolboxR.so
LIB_DYN = libFSToolboxR.dylib


all: shared

shared_test: $(LIB_SO) test.o

static: $(LIB_A)

shared: $(LIB_SO)

dynamic: $(LIB_DYN)

$(LIB_A): $(objects)
	$(RM) -f $(LIB_A)
	ar -csq $(LIB_A) $(objects)

$(LIB_SO): $(objects)
	$(RM) -f $(LIB_SO)
	$(LINKER) -shared -Wl,-soname,$(LIB_SO) -o $(LIB_SO) $(objects) -I$(MITOOLBOXPATH) -lm

$(LIB_DYN): $(objects)
	$(RM) -f $(LIB_DYN)
	$(COMPILER) $(CXXFLAGS) -dynamiclib -install_name $(LIB_DYN) -o $(LIB_DYN) $(objects) -I$(MITOOLBOXPATH) -lm

FSToolboxR.o : FSToolboxR.c
	$(COMPILER) $(CXXFLAGS) -c FSToolboxR.c -I$(MITOOLBOXPATH)

mRMR_D.o: mRMR_D.c  
	$(COMPILER) $(CXXFLAGS) -c mRMR_D.c -I$(MITOOLBOXPATH)

CMIM.o: CMIM.c 
	$(COMPILER) $(CXXFLAGS) -c CMIM.c -I$(MITOOLBOXPATH)

MIM.o: MIM.c 
	$(COMPILER) $(CXXFLAGS) -c MIM.c -I$(MITOOLBOXPATH)

JMI.o: JMI.c 
	$(COMPILER) $(CXXFLAGS) -c JMI.c -I$(MITOOLBOXPATH)

DISR.o: DISR.c 
	$(COMPILER) $(CXXFLAGS) -c DISR.c -I$(MITOOLBOXPATH)

CondMI.o: CondMI.c 
	$(COMPILER) $(CXXFLAGS) -c CondMI.c -I$(MITOOLBOXPATH)

ICAP.o: ICAP.c 
	$(COMPILER) $(CXXFLAGS) -c ICAP.c -I$(MITOOLBOXPATH)

BetaGamma.o: BetaGamma.c 
	$(COMPILER) $(CXXFLAGS) -c BetaGamma.c -I$(MITOOLBOXPATH)

MutualInformation.o: $(MITOOLBOXPATH)MutualInformation.c 
	$(COMPILER) $(CXXFLAGS) -c $(MITOOLBOXPATH)MutualInformation.c -I$(MITOOLBOXPATH)

Entropy.o: $(MITOOLBOXPATH)Entropy.c 
	$(COMPILER) $(CXXFLAGS) -c $(MITOOLBOXPATH)Entropy.c -I$(MITOOLBOXPATH)

CalculateProbability.o: $(MITOOLBOXPATH)CalculateProbability.c 
	$(COMPILER) $(CXXFLAGS) -c $(MITOOLBOXPATH)CalculateProbability.c -I$(MITOOLBOXPATH)

ArrayOperations.o: $(MITOOLBOXPATH)ArrayOperations.c 
	$(COMPILER) $(CXXFLAGS) -c $(MITOOLBOXPATH)ArrayOperations.c -I$(MITOOLBOXPATH)

test.o: test.c 
	$(COMPILER) $(CXXFLAGS) test.c -o test $(LIB_SO) -L$(ROOT_DIR) -Wl,-rpath=.
	$(MAKE_RUNABLE)
	
.PHONY : debug
debug:
	$(MAKE) shared_test "CXXFLAGS = -fPIC -DCOMPILE_C -g -ggdb -DDEBUG"

.PHONY : test
test:
	$(MAKE) shared_test "CXXFLAGS = -O3 -fPIC -DCOMPILE_C"
	
.PHONY : x86
x86:
	$(MAKE) shared "CXXFLAGS = -O3 -fPIC -DCOMPILE_C -m32"
	
.PHONY : x64
x64:
	$(MAKE) shared "CXXFLAGS = -O3 -fPIC -DCOMPILE_C -m64"
	
.PHONY : intel
intel:
	$(MAKE) shared "COMPILER = icc" "LINKER = icc" "CXXFLAGS = -O2 -fPIC -xHost"

.PHONY : clean
clean:
	$(RM) *.so *.o *.so.a test test.exe

.PHONY : install
install:
	$(MAKE)
	@echo "Installing FSToolboxR FSToolboxR.so to $(PREFIX)/lib"
	@cp -v libFSToolboxR.so $(PREFIX)/lib
	@echo "Installing FSToolboxR header files to $(PREFIX)/include/FSToolboxR"
	@mkdir -p $(PREFIX)/include/FSToolboxR
	@cp -v FSToolbox.h FSAlgorithms.h $(PREFIX)/include/FSToolboxR/
