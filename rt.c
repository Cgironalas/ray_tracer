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
static long double Ia = 0.002;

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
	long double Ip;
};

struct Object {
	long double Kd;
	long double Ka;
	struct Color color;
};

struct Intersection {
	long double Xi;
	long double Yi;
	long double Zi;
	long double distance;
	struct Object * object;
};

static struct Light *Lights;
static struct Object *Objects;
static struct Color BACKGROUND = {0.6, 0.6, 0.6};
static struct Color Framebuffer[1366][768];

long double min(long double a, long double b){
    if(a < b) { return a; }
    else { return b; }
}

void getSceneObjects(){

}

void saveFile(){
	int i, j;
	FILE *file;
	file = fopen("scene.ppm", "w");
	if(file == NULL){
		printf("Error creating/opening file!\n");
		exit(1);
	}

	fprintf(file, "%s\n", "P3");
	fprintf(file, "%i %i\n", Hres, Vres);
	fprintf(file, "%i\n", 255);
	for (i = 0; i < Hres; i++){
		for (j = 0; j < Vres; j++){
			int R = (int) Framebuffer[i][j].r / 255;
			int G = (int) Framebuffer[i][j].g / 255;
			int B = (int) Framebuffer[i][j].b / 255;
			fprintf(file, "%i %i %i   ", R, G, B);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

struct Intersection * getFirstIntersection(struct Vector a, struct Vector b){
	struct Intersection * intersection;
	long double tmin;

	intersection = NULL;
	tmin = 100000;

	int k;
	int nObjects = sizeof(Objects)/sizeof(Objects[0]);
	for(k = 0; k < nObjects; k++){
		//intersection = calcIntesection(Objects[k]->function, a);
		if(intersection){
			tmin = intersection->distance;
		}
	}
	return (intersection);
}

long double getAtunuationFactor(){
	return 0;
}

//DONE
long double pointProduct(struct Vector a, struct Vector b){
	long double pp = 0;
	pp += (a.x * b.x);
	pp += (a.y * b.y);
	pp += (a.z * b.z);
	return pp;
}

//DONE
struct Color colorXintensity(long double I, struct Color color){
	struct Color newColor;
	
	newColor.r = I * color.r;
	newColor.g = I * color.g;
	newColor.b = I * color.b;

	return newColor;
}

//DONE
struct Color getColor(struct Vector a, struct Vector b){
	struct Color color;
	struct Intersection * intersection;

	intersection = getFirstIntersection(a, b);

	if(!intersection){
		color = BACKGROUND;
	}else{
		int k;
		int nLights = sizeof(Lights)/sizeof(Lights[0]);
		
		struct Object * Q = intersection->object;
		struct Vector N; //= normal unitaria a Q en punto (Xi, Yi, Zi) AQUI
		long double I = 0.0;
		
		long double Fatt;
		struct Vector L;
		for(k = 0; k < nLights; k++){
			//L = getUnitaryVector(k); AQUI
			long double pp = pointProduct(N, L);
			if(pp > 0.0){
				Fatt = getAtunuationFactor();
				I = I + (pp * Q->Kd * Fatt * Lights[k].Ip);
			}
		}

		I = I + Ia * Q->Ka;
		I = min(1.0, I);
		color = colorXintensity(I, Q->color);
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
	
	for (i = 0; i < Hres; i++){
		for (j = 0; j < Vres; j++){
			Xw = (long double) ((i + (1/2)) * (Xmax - Xmin))/Hres + Xmin;
			Yw = (long double) ((i + (1/2)) * (Ymax - Ymin))/Vres + Ymin;
			
			L = sqrt(pow(Xw - Xe, 2) + pow(Yw - Ye, 2) + pow(Zw - Ze, 2));

			Xd = (long double) (Xw - Xe) / L;
			Yd = (long double) (Yw - Ye) / L;
			Zd = (long double) (Zw - Ze) / L;

			direction.x = Xd;
			direction.y = Yd;
			direction.z = Zd;

			color = getColor(anchor, direction);

			Framebuffer[i][j] = color;
		}
	}
	saveFile();
}