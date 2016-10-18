#ifndef  STRUCTURES_H
#define	STRUCTURES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "TDA_LDE.h"


//Image resolution
extern int Hres;
extern int Vres;

//Eye
extern long double Xe;
extern long double Ye;
extern long double Ze;

//Window
extern long double Xmax;
extern long double Ymax;
extern long double Xmin;
extern long double Ymin;

//Others
extern long double Ia;  //Iluminacion ambiente
extern long double e; //Epsilon 
extern time_t t;

extern struct Color { 
	long double r;
	long double g;
	long double b;
};

extern struct Vector{//Aqui podria ir un valor que sea de norma, por si es necesario usarlo varias veces
	long double x;
	long double y;
	long double z;
};

extern struct Light{
	long double Xp;
	long double Yp;
	long double Zp;
	long double c1; //Constante
	long double c2; // Lineal
	long double c3; // Cuadrática
	long double Ip;
};

extern struct Object {
	long double Xc;
	long double Yc;
	long double Zc;
	long double* otherData; //Radio para esferas y demás valores para otros objetos
	struct Intersection *(*intersectionFuncion)(struct Vector, struct Vector, struct Object);
	long double Kd; // ¿Qué tanta iluminación absorbe? Coef. Reflexíon DIfusa
	long double Ka; //¿Qué tanta iluminación ambiente absorbe? COef. Ilum. AMbiente.
	long double Kn; //FActor de reflexión especular
	struct Color color;
};

extern struct Intersection {
	struct Vector *(*normalVector)();
	long double Xi;
	long double Yi;
	long double Zi;
	long double distance;
	struct Object * object;
};

extern int numberObjects;
extern struct Color BACKGROUND; //Gris
extern struct Color Framebuffer[768][1366];



#endif	/* STRUCTURES_H */