#!/bin/bash
CC=g++

CFLAGS=-c -g -Wall -fPIC $(shell root-config --cflags)
LDFLAGS=$(shell root-config --glibs)
SOURCES=TestFSR.cpp FSR.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=TestFSR
LIBRARY=libFSR.so
BINDIR=./../bin
LIBDIR=./lib


all: $(SOURCES) $(EXECUTABLE) $(LIBRARY)

$(LIBRARY): $(OBJECTS) $(LIBDIR)
		$(CC) -shared $(OBJECTS) -o $(LIBDIR)/$@

$(EXECUTABLE): $(OBJECTS) $(BINDIR)
		$(CC) $(OBJECTS) $(LDFLAGS) -o $(BINDIR)/$@

$(BINDIR):
		mkdir -p $@

$(LIBDIR):
		mkdir -p $@

.cpp.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
		rm ./*~ ./*.o $(BINDIR)/$(EXECUTABLE) $(LIBDIR)/$(LIBRARY) 


