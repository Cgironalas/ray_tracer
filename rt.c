#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "malloc.h"
#include <time.h>

//Image resolution
static int Hres = 1008;
static int Vres = 1008;

//Window
static long double Xmax = 400;
static long double Ymax = 400;
static long double Xmin = 0;
static long double Ymin = 0;

//Others
static long double Ia = 0.6;
static long double e = 0.05;
static time_t t;

struct Color { 
	long double r;
	long double g;
	long double b;
};

struct Vector{
	long double x;
	long double y;
	long double z;
};

struct Light{
	long double Xp;
	long double Yp;
	long double Zp;	
	long double c1;
	long double c2;
	long double c3;
	long double Ip;
};

struct Object {
	long double Xc;
	long double Yc;
	long double Zc;
	long double Kd;
	long double Ka;
	long double otherData;
	struct Color color;
	//struct Vector *(*normalVector)();
	//struct Intersection *(*intersectionFuncion)(struct Vector, struct Vector, struct Object);
};

struct Intersection {
	long double Xi;
	long double Yi;
	long double Zi;
	long double distance;
	struct Object object;
};

static struct Light Lights[1];
static struct Object Objects[2];

static struct Intersection tempIntersect;

static struct Vector eye = {200,200,-1500};
static struct Color Framebuffer[1008][1008];
static struct Color BACKGROUND = {0.3, 0.3, 0.3};

long double min(long double a, long double b){
    if(a < b) { return a; }
    else { return b; }
}
long double max(long double a, long double b){
	if(a > b) { return a; }
	else { return b; }
}



//////////////Vectors
//DONE
//Producto punto entre dos vectores
long double pointProduct(struct Vector a, struct Vector b){
	long double pp = 0;

	pp += (a.x * b.x);
	pp += (a.y * b.y);
	pp += (a.z * b.z);
	
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
long double getNorm(struct Vector vector){
	long double norm = sqrt(pow(vector.x ,2) + pow(vector.y ,2) + pow(vector.z ,2));
	return norm;
}

//DONE
//Regresa un vector normalizado
struct Vector normalize(struct Vector vector){
	long double norm = getNorm(vector);
	struct Vector unitVector;

	unitVector.x = vector.x / norm;
	unitVector.y = vector.y / norm;
	unitVector.z = vector.z / norm;

	return unitVector;
}
//////////////END Vectors



//////////////Files

//Leer archivos con la escena
void getSceneObjects(){
	//Pendiente
	//
}

//DONE
//Guarda el framebuffer en una imagen ppm
void saveFile(){
	int i, j;
	FILE *file;
	file = fopen("scene.ppm", "w");
	if(file == NULL){
		printf("Error creating/opening file!\n");
		exit(1);
	}

	//Formato necesario para un PPM
	fprintf(file, "%s\n", "P3");
	fprintf(file, "%i %i\n", Hres, Vres);
	fprintf(file, "%i\n", 255);


	//srand((unsigned) time(&t));
	//FOR parra recorrer el framebuffer escribiendo el color de cada pixel en el PPM
	//El formato es para que "se vea como matriz" en el PPM
	for (i = 0; i < Vres; i++){
		for (j = 0; j < Hres; j++){
			int R = (int) 255 * Framebuffer[i][j].r;//rand() % 255;//
			int G = (int) 255 * Framebuffer[i][j].g;
			int B = (int) 255 * Framebuffer[i][j].b;
			fprintf(file, "%i %i %i   ", R, G, B);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}
//////////////END Files



//////////////Ilumination models

//DONE
//Calcula el factor de atenuacion de una luz
long double getAtunuationFactor(struct Light light, long double distance){
	long double value = (long double) 1 / (light.c1 + (light.c2 * distance) + (light.c3 * pow(distance, 2)));
	return min(1.0, value);
}

//DONE
//Regresa el color base de un objeto al que se le aplica la intensidad de iluminacion difusa
struct Color colorXintensity(long double I, struct Color color){
	struct Color newColor;
	
	newColor.r = I * color.r;
	newColor.g = I * color.g;
	newColor.b = I * color.b;

	return newColor;
}
//////////////END Ilumination models



//////////////Ray Tracer Stuff

//Interseccion con esferas
//DONE
struct Intersection *sphereIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	long double t, t1, t2;

	long double Xdif = anchor.x - object.Xc;
	long double Ydif = anchor.y - object.Yc;
	long double Zdif = anchor.z - object.Zc;

	long double B = 2 * ((direction.x * Xdif) + (direction.y * Ydif) + (direction.z * Zdif));
	long double C = pow(Xdif, 2) + pow(Ydif, 2) + pow(Zdif, 2) - pow(object.otherData,2);

	long double discriminant = pow(B, 2) - (4 * C);

	if(discriminant >= 0){
		long double root = sqrt(discriminant);

		B *= -1;

		t1 = (B + root)/2;
		t2 = (B - root)/2;

		if(t1 > e){
			if(t2 > e){
				t = min(t1, t2);
			}else{
				t = t1;
			}
		}else{
			if(t2 > e){
				t = t2;
			}else{
				return NULL;
			}
		}
		
		tempIntersect.distance = t;
		tempIntersect.object = object;
		tempIntersect.Xi = anchor.x + (t * direction.x);
		tempIntersect.Yi = anchor.y + (t * direction.y);
		tempIntersect.Zi = anchor.z + (t * direction.z);

		return &tempIntersect;
	}else{
		return NULL;
	}
}

struct Vector sphereNormal(struct Object object, long double Xi, long double Yi, long double Zi) {
	struct Vector normal;
	
	normal.x = Xi - object.Xc;
	normal.y = Yi - object.Yc;
	normal.z = Zi - object.Zc;

	normal = normalize(normal);

	return normal;
}

//Revisar
struct Intersection getFirstIntersection(struct Vector anchor, struct Vector direction){
	int k;
	long double tmin;
	struct Intersection intersection;
	struct Intersection * tempIntersection;

	tmin = 100000;
	tempIntersection = NULL;
	intersection.distance = -1;

	int objectsAmount = sizeof(Objects)/sizeof(Objects[0]);
	for(k = 0; k < objectsAmount; k++){
		tempIntersection = sphereIntersection (anchor, direction, Objects[k]);
		if(tempIntersection != NULL && tempIntersection->distance > e && tempIntersection->distance < tmin){
			tmin = tempIntersection->distance;
			intersection.Xi = tempIntersection->Xi;
			intersection.Yi = tempIntersection->Yi;
			intersection.Zi = tempIntersection->Zi;
			intersection.object = tempIntersection->object;
			intersection.distance = tempIntersection->distance;
		}
		tempIntersection = NULL;
	}
	return intersection;
}

//Pendiente
//Funcion de que color del profe
struct Color getColor(struct Vector anchor, struct Vector direction){
	struct Color color;
	struct Intersection intersection;
	struct Intersection *tempIntersection;
	
	intersection = getFirstIntersection(anchor, direction);
	
	if(intersection.distance == -1){
		color = BACKGROUND;
	}else{
		int k;
		int lightsAmount = sizeof(Lights)/sizeof(Lights[0]);
		
		struct Vector L;
		struct Object Q = intersection.object;
		struct Vector intersectVector = {intersection.Xi, intersection.Yi, intersection.Zi};
		struct Vector N = sphereNormal(Q, intersection.Xi, intersection.Yi, intersection.Zi); 
		
		long double Fatt;
		long double I = 0.0;

		for(k = 0; k < lightsAmount; k++){
			struct Vector intresectionPoint;
			struct Intersection obstacle;
			struct Vector luz = {Lights[k].Xp - intersection.Xi, Lights[k].Yp - intersection.Yi, Lights[k].Zp - intersection.Zi};
			L = normalize(luz);
			obstacle = getFirstIntersection(intersectVector, L);
			if(obstacle.distance < e){
				long double pp = pointProduct(N, L);
				long double distanceToLight = getNorm(L);
				if(pp > 0.0){
					Fatt = getAtunuationFactor(Lights[k], distanceToLight);
					I = I + (pp * Q.Kd * Fatt * Lights[k].Ip);
				}
			}
		}
		I = I + Ia * Q.Ka;
		I = min(1.0, I);
		color = colorXintensity(I, Q.color);
	}
	return (color);
}

//DONE
int main(int argc, char *arcgv[]){
	int i, j;

	long double L;
	long double Xw, Yw;
	long double Xd, Yd, Zd;

	struct Color color;
	struct Vector direction;
	
	long double Xdif = Xmax - Xmin;
	long double Ydif = Ymax - Ymin;

	//Esfera ad hoc
	struct Object sphere1;
	sphere1.Xc = 200;
	sphere1.Yc = 200;
	sphere1.Zc = 650;
	sphere1.otherData = 200;
	

	sphere1.Kd = 0.8;
	sphere1.Ka = 0.4;

	sphere1.color.r = 1;
	sphere1.color.g = 0;	
	sphere1.color.b = 0;

	Objects[0] = sphere1;

	sphere1.Xc = 100;
	sphere1.Yc = 200;
	sphere1.Zc = 100;
	sphere1.otherData = 60;
	

	sphere1.Kd = 0.8;
	sphere1.Ka = 0.4;

	sphere1.color.r = 0;
	sphere1.color.g = 0;	
	sphere1.color.b = 1;

	Objects[1] = sphere1;

	struct Light light1;
	light1.Xp = -200;
	light1.Yp = 300;
	light1.Zp = -1000;
	light1.Ip = 1;
	light1.c1 = 1;
	light1.c2 = 0;
	light1.c3 = 0;
	Lights[0] = light1;

	Zd = -eye.z;
	for (i = 0; i < Vres; i++){
		Yw = (long double) ((i + (1/2)) * Ydif)/Vres + Ymin;
		Yd = Yw - eye.x;
		for (j = 0; j < Hres; j++){
			Xw = (long double) ((j + (1/2)) * Xdif)/Hres + Xmin;
			Xd = Xw - eye.y;

			L = sqrt(pow(Xd, 2) + pow(Yd, 2) + pow(Zd, 2));
		
			direction.x = Xd / L;
			direction.y = Yd / L;
			direction.z = Zd / L;

			color = getColor(eye, direction);

			Framebuffer[i][j] = color;
		}
	}
	saveFile();
}
//////////////END Ray Tracer Stuff