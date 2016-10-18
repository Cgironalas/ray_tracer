
#ifndef  VECTORS_H
#define	VECTORS_H


#include "structures.h"

//////////////Vectors
//DONE
//Producto punto entre dos vectores
long double pointProduct(struct Vector*, struct Vector*);

//No recuerdo si esta es la formula correcta////////////////////
//Producto cruz entre dos vectores
struct Vector crossProduct(struct Vector, struct Vector);

//DONE
//Regresa la norma de un vector
long double getNorm(struct Vector*);

//DONE
//Regresa un vector normalizado
struct Vector normalize(struct Vector*);




#endif	/* VECTORS_H */