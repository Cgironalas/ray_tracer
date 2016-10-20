#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "malloc.h"
#include <time.h>

//Image resolution
static int Hres = 1366;
static int Vres = 768;

//Eye
static long double Xe = 0;
static long double Ye = 0;
static long double Ze = -200;

//Window
static long double Xmax = 100;
static long double Ymax = 100;
static long double Xmin = -100;
static long double Ymin = -100;

//Others
static long double Ia = 0.5;
static long double e = 0.005;
static time_t t;

static int cont = 0;

struct Color { 
	long double r;
	long double g;
	long double b;
};

struct Vector{//Aqui podria ir un valor que sea de norma, por si es necesario usarlo varias veces
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
	long double otherData;
	//struct Intersection *(*intersectionFuncion)(struct Vector, struct Vector, struct Object);
	long double Kd;
	long double Ka;
	struct Color color;
};


struct Intersection {
	//struct Vector *(*normalVector)();
	long double Xi;
	long double Yi;
	long double Zi;
	long double distance;
	struct Object object;
};

static struct Light *Lights;
static struct Object Objects[1];
static struct Color BACKGROUND = {0.3, 0.3, 0.3};
static struct Color Framebuffer[1366][768];
static struct Intersection tempIntersect;

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
}

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
//Pendiente
struct Intersection *sphereIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	long double t, t1, t2;

	long double Xdif = anchor.x - object.Xc;
	long double Ydif = anchor.y - object.Yc;
	long double Zdif = anchor.z - object.Zc;

	long double B = 2 * ((direction.x * Xdif) + (direction.y * Ydif) + (direction.z * Zdif));
	long double C = (long double)pow(Xdif, 2) + pow(Ydif, 2) + pow(Zdif, 2) - pow(object.otherData,2);

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
		printf("%LF\n", direction.y);

		tempIntersect.object = object;

		tempIntersect.Xi = anchor.x + t * direction.x;
		tempIntersect.Yi = anchor.y + t * direction.y;
		tempIntersect.Zi = anchor.z + t * direction.z;

		printf("%LF %LF %LF \n", tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi);
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
/*struct Intersection getFirstIntersection(struct Vector anchor, struct Vector direction){
	struct Intersection * intersection;
	struct Intersection * tempIntersection;
	long double tmin;

	intersection = NULL;
	tmin = 100000;

	int k;
	int objectsAmount = sizeof(Objects)/sizeof(Objects[0]);
	for(k = 0; k < objectsAmount; k++){
		tempIntersection = /*Objects[k]. sphereIntersection (anchor, direction, Objects[k]);
		if(tempIntersection && tempIntersection->distance > e && tempIntersection->distance < tmin){
			tmin = tempIntersection->distance;
			intersection = tempIntersection;
		}
	}
	return (intersection);
}*/

//Pendiente
//Funcion de que color del profe
struct Color getColor(struct Vector anchor, struct Vector direction){
	struct Color color;
	struct Intersection *tempIntersection;
	struct Intersection intersection;
	
	tempIntersection = sphereIntersection(anchor, direction, Objects[0]);
	
	if(!tempIntersection){
		color = BACKGROUND;
	}else{
		printf("heyoo");

		intersection.Xi = tempIntersection->Xi;
		intersection.Yi = tempIntersection->Yi;
		intersection.Zi = tempIntersection->Zi;
		intersection.distance = tempIntersection->distance;
		intersection.object = tempIntersection->object;

		int k;
		int lightsAmount = sizeof(Lights)/sizeof(Lights[0]);
		
		struct Object Q = intersection.object;
		struct Vector N = sphereNormal(intersection.object,intersection.Xi, intersection.Yi, intersection.Zi); 
		

		long double I = 0.0;
		
		long double Fatt;
		struct Vector L;
		for(k = 0; k < lightsAmount; k++){
			struct Vector intresectionPoint;
			//struct Intersection obstacle = NULL;
			//obstacle = getFirstIntersection({intersection.Xi, intersection.Yi, intersection.Zi},L);
			if(1 == 1){
				struct Vector luz = {Lights[k].Xp-intersection.Xi, Lights[k].Yp-intersection.Yi, Lights[k].Zp-intersection.Zi};
				L = normalize(luz);


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
	long double Zw = 0;

	long double Xd, Yd, Zd;

	struct Color color;
	struct Vector direction;
	struct Vector anchor = {Xe, Ye, Ze};
	
	long double Xdif = Xmax - Xmin;
	long double Ydif = Ymax - Ymin;


	//Esfera ad hoc
	struct Object esfera1;
	esfera1.Xc = 0;
	esfera1.Yc = 0;
	esfera1.Zc = 30;
	esfera1.otherData = 10;
	

	esfera1.Kd = 1;
	esfera1.Ka = 1;

	esfera1.color.r = 1;
	esfera1.color.g = 0;	
	esfera1.color.b = 0;

	Objects[0] = esfera1;

	for (i = 0; i < Hres; i++){
		for (j = 0; j < Vres; j++){
			Xw = (long double) ((i + (1/2)) * Xdif)/Hres + Xmin;
			Yw = (long double) ((j + (1/2)) * Ydif)/Vres + Ymin;

			Xd = Xw - Xe;
			Yd = Yw - Ye;
			Zd = Zw - Ze;


			L = sqrt(pow(Xd, 2) + pow(Yd, 2) + pow(Zd, 2));
		
			direction.x = Xd / L;
			direction.y = Yd / L;
			direction.z = Zd / L;

			color = getColor(anchor, direction);

			Framebuffer[i][j] = color;

		}
	}
	saveFile();
}


//////////////END Ray Tracer Stuff