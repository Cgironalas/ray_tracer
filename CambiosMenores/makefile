CC=gcc

SOURCE= rt.c
OUTPUT=rt
OBJECT=rt.o
LDFLAGS= -lm


$(OUTPUT): $(OBJECT)
	$(CC) $(OBJECT) -o $(OUTPUT) $(LDFLAGS) 

$(OBJECT): $(SOURCE)
	$(CC) -c $(SOURCE) -o rt.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f rt
