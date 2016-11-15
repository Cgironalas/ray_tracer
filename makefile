CC=gcc

SOURCE= rt.c
OUTPUT=rt
OBJECT=rt.o
LDFLAGS= -lpthread -lm


$(OUTPUT): $(OBJECT)
	$(CC) $(OBJECT) -o $(OUTPUT) $(LDFLAGS) 

$(OBJECT): $(SOURCE)
	$(CC) -c $(SOURCE) -o rt.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f rt
