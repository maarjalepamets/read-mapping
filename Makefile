# Project: GMapper

VERSION = 4.0

CPP  = g++

# C Files

INDEXER_SOURCES = \
	indexer.c \
	utils.c

MAPPER_SOURCES = \
        mapper.c \
        utils.c \
        mappermethods.c \
        matcher.c

RELEASEFLAGS = -O3
DEBUGFLAGS = -O0 -g
LIBS = -lm
INCS = -I.
BINS  = indexer mapper
CXXFLAGS = $(INCS) $(DEBUGFLAGS) -D "VERSION=\"${VERSION}\"" -Wall
#CXXFLAGS = $(INCS) $(RELEASEFLAGS) -D "VERSION=\"${VERSION}\"" -Wall

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BINS) all-after

indexer: $(INDEXER_SOURCES)
	$(CPP) $(INDEXER_SOURCES) -o indexer $(LIBS) $(CXXFLAGS)

mapper: $(MAPPER_SOURCES)
	$(CPP) $(MAPPER_SOURCES) -o mapper $(LIBS) $(CXXFLAGS)

clean: clean-custom
	rm -f *.o $(BINS)

depend:
	$(CC) $(CFLAGS) -M *.c > .depend

include .depend
