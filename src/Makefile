# Core verification is portable; optional backends and searches are explicit.
CC ?= cc
CXX ?= c++
CFLAGS ?= -O3 -Wall -Wextra -Wno-unused-parameter
CXXFLAGS ?= -O2 -std=c++14 -Wall
PKG_CONFIG ?= pkg-config
FLINT_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags flint 2>/dev/null)
FLINT_LIBS ?= $(shell $(PKG_CONFIG) --libs flint 2>/dev/null || echo -lflint)
FPLLL_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags fplll 2>/dev/null)
FPLLL_LIBS ?= $(shell $(PKG_CONFIG) --libs fplll 2>/dev/null || echo -lfplll)
PARI_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags libpari 2>/dev/null)
PARI_LIBS ?= $(shell $(PKG_CONFIG) --libs libpari 2>/dev/null || echo -lpari)
MPFI_CFLAGS ?=
MPFI_LIBS ?= -lmpfi
PLATT_INCLUDE ?= ../work/deps/platt
NUM_LIBS = $(FLINT_LIBS) -lmpfr -lgmp -lm
PLATTFLAGS = -O3 -msse2 -mfpmath=sse -frounding-math -I $(PLATT_INCLUDE) -I shim -DIALIB=PLATT

.PHONY: all search test test-verification clean
all: certify euclid
search: ne_sweep ne_chain

certify: certify.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FLINT_CFLAGS) $(FPLLL_CFLAGS) $< -o $@ $(LDFLAGS) $(FPLLL_LIBS) $(NUM_LIBS)

euclid: euclid.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FLINT_CFLAGS) $< -o $@ $(LDFLAGS) $(NUM_LIBS)

euclid_platt: euclid.cpp $(PLATT_INCLUDE)/int_double14.2.h shim/crlibm.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PLATTFLAGS) $(FLINT_CFLAGS) $< -o $@ $(LDFLAGS) $(NUM_LIBS)

euclid_mpfi: euclid.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FLINT_CFLAGS) $(MPFI_CFLAGS) -DIALIB=MPFI $< -o $@ $(LDFLAGS) $(MPFI_LIBS) $(NUM_LIBS)

ne_sweep: ne_sweep.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PARI_CFLAGS) $< -o $@ $(LDFLAGS) $(PARI_LIBS) -lm

ne_chain: ne_chain.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PARI_CFLAGS) $< -o $@ $(LDFLAGS) $(PARI_LIBS) -lm

test_rounding: test_rounding.cpp euclid.cpp $(PLATT_INCLUDE)/int_double14.2.h shim/crlibm.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PLATTFLAGS) $(FLINT_CFLAGS) $< -o $@ $(LDFLAGS) $(NUM_LIBS)

test:
	python3 -B test_verification.py
	python3 -B test_release.py

test-verification: test test_rounding
	./test_rounding

clean:
	rm -f certify euclid euclid_platt euclid_mpfi ne_sweep ne_chain test_rounding
