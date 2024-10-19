CPPFLAGS=-I.
CFLAGS=-std=c2x -Wall -Wpedantic -Werror
LDLIBS=-lm

all : src/draft

src/%.o: surface/%.h

src/draft : src/draft.o libsurface.a

libsurface.a : src/vector.o src/surface.o
	$(AR) crsu $@ $^

clean:
	-rm libsurface.a src/vector.o src/surface.o
