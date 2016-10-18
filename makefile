CC=gcc

OUTPUT=rt
OBJECT=rt.o structures.o vectors.o  colorLighting.o TDA_LDE.o fileReader.o
LDFLAGS= -lglut -lGL -lGLU -lm
CFLAGS = -I/usr/include/GL



#Se compila cada archivo objeto de cada archivo fuente c y dependencias/headers
%.o: %.c 
	$(CC) -c $(LDFLAGS) $< 


$(OUTPUT): $(OBJECT)
	$(CC) $(OBJECT) -o $(OUTPUT) $(LDFLAGS) $(CFLAGS) 

clean:
	rm -f *.o
	rm -f rt
