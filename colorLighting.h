
#ifndef  COLORLIGHTING_H
#define	COLORLIGHTING_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "TDA_LDE.h"
#include "structures.h"


long double min(long double, long double);

long double max(long double, long double);





//////////////Ilumination models

//DONE
//Calcula el factor de atenuacion de una luz
long double getAtunuationFactor(struct Light,  long double);

//DONE
//Regresa el color base de un objeto al que se le aplica la intensidad de iluminacion difusa
struct Color colorXintensity(long double,  struct Color);




//////////////END Ilumination models




#endif	/* COLORLIGHTING_H */