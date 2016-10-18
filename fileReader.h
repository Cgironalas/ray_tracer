#ifndef  FILE_READER_H
#define	FILE_READER_H



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "TDA_LDE.h"
#include "structures.h"
#include "vectors.h"
#include "colorLighting.h"
#include "malloc.h"



//Esta variable sólo VIVE en este archivo.
static char* escenaFile;



long double *obtainPointFromString(char []);


long double *readValueFromLine(int , int*, char*, int*);


void createObjectFromData(long double*, int);

//Leer archivos con la escena
void getSceneObjects();



#endif /* FILE_READER_H */