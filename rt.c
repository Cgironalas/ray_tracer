#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "TDA_LDE.h"
#include "structures.h"
#include "vectors.h"
#include "colorLighting.h"
#include "malloc.h"
#include "fileReader.h"



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




/////////////Ray Tracer Stuff

//Interseccion con esferas
//Pendiente
/*
struct Intersection * sphereIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	struct Intersection * intersection;
	long double t, t1, t2;

	long double Xdif = anchor.x - object.Xc;
	long double Ydif = anchor.y - object.Yc;
	long double Zdif = anchor.z - object.Zc;

	long double B = 2 * ((direction.x * Xdif) + (direction.y * Ydif) + (direction.z * Zdif));
	long double C = pow(Xdif, 2) + pow(Ydif, 2) + pow(Zdif, 2) + object.otherData[0];

	long double discriminant = pow(B, 2) - (4 * C);

	if(discriminant >= 0){
		long double root = sqrt(discriminant);
		B *= -1;
		t1 = (B + sqrt(discriminant))/2;
		t2 = (B - sqrt(discriminant))/2;

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

		intersection->distance = t;
		intersection->object = &object;
		intersection->Xi = anchor.x + t * direction.x;
		intersection->Yi = anchor.y + t * direction.y;
		intersection->Zi = anchor.z + t * direction.z;
	}else{
		return NULL;
	}
}

//Revisar
struct Intersection * getFirstIntersection(struct Vector anchor, struct Vector direction){
	struct Intersection * intersection;
	struct Intersection * tempIntersection;
	long double tmin;

	intersection = NULL;
	tmin = 100000; //Rango de visión máxima

	int k;
	int objectsAmount = sizeof(Objects)/sizeof(Objects[0]);
	for(k = 0; k < objectsAmount; k++){
		tempIntersection = Objects[k].intersectionFuncion (anchor, direction, Objects[k]);
		if(tempIntersection && tempIntersection->distance > e && tempIntersection->distance < tmin){
			tmin = tempIntersection->distance;
			intersection = tempIntersection;
		}
	}
	return (intersection);
}

//Pendiente
//Funcion de que color del profe
struct Color getColor(struct Vector anchor, struct Vector direction){
	struct Color color;
	struct Intersection * intersection;

	intersection = getFirstIntersection(anchor, direction);

	if(!intersection){
		color = BACKGROUND;
	}else{
		int k;
		int lightsAmount = sizeof(Lights)/sizeof(Lights[0]);
		
		struct Object *Q = intersection->object;
		struct Vector *N = intersection->normalVector(); //normal unitaria a Q en punto (Xi, Yi, Zi) AQUI
		
		long double I = 0.0;
		
		long double Fatt;
		struct Vector *L;
		for(k = 0; k < lightsAmount; k++){
			struct Vector intresectionPoint;
			struct Intersection * obstacle = NULL;
			//obstacle = getFirstIntersection({intersection.Xi, intersection.Yi, intersection.Zi},L);
			if(!obstacle){
				//L = getUnitaryVector(k); AQUI
				long double pp = pointProduct(N, L);
				long double distanceToLight = getNorm(L);
				if(pp > 0.0){
					Fatt = getAtunuationFactor(Lights[k], distanceToLight);
					I = I + (pp * Q->Kd * Fatt * Lights[k].Ip);
				}
			}
		}

		I = I + Ia * Q->Ka;
		I = min(1.0, I);
		color = colorXintensity(I, Q->color);
	}
	return (color);
}
*/


void initRayTracer(){

}



//DONE
int main(int argc, char *argv[]){
	deleteListObjects(Objects);
	deleteListLights(Lights);

	int i, j;

	long double L;
	long double Xw, Yw;
	long double Zw = 0;

	long double Xd, Yd, Zd;

	struct Color color;
	struct Vector direction;
	struct Vector anchor = {Xe, Ye, Ze};
	/*
	long double Xdif = Xmax - Xmin;
	long double Ydif = Ymax - Ymin;

	for (i = 0; i < Vres; i++){
		for (j = 0; j < Hres; j++){
			Xw = (long double) ((i + (1/2)) * Xdif)/Hres + Xmin;
			Yw = (long double) ((i + (1/2)) * Ydif)/Vres + Ymin;

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
	}*/
	//getSceneObjects();
	Ia = 0.002;
	
	long double valuesLight[7];
	long double valuesSphere[10];
	valuesLight[0] = 20.0; //Xp
	valuesLight[1] = 20.0; //Yp
	valuesLight[2] = 30.0; //Zp
	valuesLight[3] = 0.0; //c1
	valuesLight[4] = 0.3; //c2
	valuesLight[5] = 0.0; //c3
	valuesLight[6] = 0.7; //Ip

	valuesSphere[0] = 50.0; //Xc
	valuesSphere[1] = 50.0; //Yc
	valuesSphere[2] = 50.0; //Zc 
	valuesSphere[3] = 10.0; //radio  as OtherDAta
	valuesSphere[4] = 0.4; //Kd
	valuesSphere[5] = 0.4; //Ka
	valuesSphere[6] = 0.2; //Kn
	valuesSphere[7] = 0.6171; //R
	valuesSphere[8] = 1.0; //G
	valuesSphere[9] = 0.2929; //B  
	
	/*Verde HEX 9EFF4B*/


	createObjectFromData(valuesLight,1); //Luz
	createObjectFromData(valuesSphere,2); //Esfera
	/*
	saveFile();
	double long *value;
	char valueString[] = "20.0,20.0,30.0";
	value=obtainPointFromString(valueString);
	printf("%LF %LF %LF", value[0], value[1],value[2]);
	free(value);*/
	deleteListObjects(Objects);
	deleteListLights(Lights);
}
//////////////END Ray Tracer Stuff