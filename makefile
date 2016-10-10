CC=gcc

SOURCE= rt.c
OUTPUT=rt
OBJECT=rt.o
LDFLAGS= -lglut -lGL -lGLU -lm
CFLAGS = -I/usr/include/GL


$(OUTPUT): $(OBJECT)
	$(CC) $(OBJECT) -o $(OUTPUT) $(LDFLAGS) $(CFLAGS) 

$(OBJECT): $(SOURCE)
	$(CC) -c $(SOURCE) -o rt.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f rt
