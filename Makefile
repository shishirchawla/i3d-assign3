# $Id: $

# Linux
# CFLAGS = -ansi -Wall -pedantic -D_GNU_SOURCE -std=c99
# LDFLAGS = -lGL -lGLU -lglut -lm

# OS X
CFLAGS = 
LDFLAGS = -framework OpenGL -framework GLUT

PROG = ass3
SOURCES = i3d_ass3.c

all : $(PROG)

$(PROG) : $(SOURCES)
	gcc -o $@ $(SOURCES) $(CFLAGS) $(LDFLAGS)

clean :·
	rm -f $(PROG)

