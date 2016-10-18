#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "structures.h"

//////////////Vectors
//DONE
//Producto punto entre dos vectores
long double pointProduct(struct Vector *a, struct Vector *	b){
	long double pp = 0;

	pp += (a->x * b->x);
	pp += (a->y * b->y);
	pp += (a->z * b->z);
	
	return pp;
}

//No recuerdo si esta es la formula correcta////////////////////
//Producto cruz entre dos vectores
struct Vector crossProduct(struct Vector a, struct Vector b){
	struct Vector newVector;

	newVector.x = (long double) (a.y * b.z) - (a.z * b.y);
	newVector.y = (long double) (a.z * b.x) - (a.x * b.z);
	newVector.z = (long double) (a.x * b.y) - (a.y * b.x);

	return newVector;
}

//DONE
//Regresa la norma de un vector
long double getNorm(struct Vector *vector){
	long double norm = sqrt(pow(vector->x ,2) + pow(vector->y ,2) + pow(vector->z ,2));
	return norm;
}

//DONE
//Regresa un vector normalizado
struct Vector normalize(struct Vector *vector){
	long double norm = getNorm(vector);
	struct Vector unitVector;

	unitVector.x = vector->x / norm;
	unitVector.y = vector->y / norm;
	unitVector.z = vector->z / norm;

	return unitVector;
}
//////////////END Vectors

