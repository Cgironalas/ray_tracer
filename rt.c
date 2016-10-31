#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "malloc.h"
#include <time.h>

//Image resolution
static int Hres;
static int Vres;

//Window
static long double Xmax = 400;
static long double Ymax = 400;
static long double Xmin = -100;
static long double Ymin = -100;

//Others
static long double Ia;
static long double e;

static int debug = 0;
static int rec = 0;

struct Color { 
	long double r;
	long double g;
	long double b;
};

struct Point2D{
	long double u;
	long double v;
};

struct Point3D{
	long double x;
	long double y;
	long double z;
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
	long double Xc; //The center may also work as the anchor
	long double Yc;
	long double Zc;
	long double Kd;
	long double Ka;
	long double Kn;
	long double Ks;
	int pointAmount;
	long double other;  //Radius of the sphere and the cilinder
	struct Vector directionVector; //for the cilinder and cone
	long double D1; 
	long double D2; //Used to cut the cilinder and the cone
	long double K1;
	long double K2; //K1 and K2 are onlt for the cone.
	struct Color color;
	struct Point2D *points2D;
	struct Point3D *points3D;
	struct Vector (*normalVector)();
	struct Intersection (*intersectionFuncion)(struct Vector, struct Vector, struct Object);
};

struct Intersection {
	long double Xi;
	long double Yi;
	long double Zi;
	long double distance;
	struct Object object;
	long double null;
};

static struct Light *Lights;
static struct Object *Objects;
static int numberObjects = 0;
static int numberLights = 0;
static int lightIndex = 0;
static int objectIndex = 0;

static struct Vector V;

static struct Vector eye;
static struct Color **Framebuffer;
static struct Color background;

static char* escenaFile = "escena1.txt";

long double min(long double a, long double b){
    if(a < b) { return a; }
    else { return b; }
}

long double max(long double a, long double b){
	if(a > b) { return a; }
	else { return b; }
}

// Vectores: ======================================================

//Producto punto entre dos vectores
long double pointProduct(struct Vector a, struct Vector b){
	long double pp = 0;

	pp += (a.x * b.x);
	pp += (a.y * b.y);
	pp += (a.z * b.z);
	
	return pp;
}

//Producto cruz entre dos vectores
struct Vector crossProduct(struct Vector a, struct Vector b){
	struct Vector newVector;

	newVector.x = (a.y * b.z) - (a.z * b.y);
	newVector.y = (a.z * b.x) - (a.x * b.z);
	newVector.z = (a.x * b.y) - (a.y * b.x);

	return newVector;
}

//Regresa la norma de un vector
long double getNorm(struct Vector vector){
	long double norm = sqrt(pow(vector.x ,2) + pow(vector.y ,2) + pow(vector.z ,2));
	return norm;
}

//Regresa un vector normalizado
struct Vector normalize(struct Vector vector){
	long double norm = getNorm(vector);
	struct Vector unitVector;
	if(norm != 0){
		unitVector.x = vector.x / norm;
		unitVector.y = vector.y / norm;
		unitVector.z = vector.z / norm;	
	}else{
		unitVector.x = vector.x;
		unitVector.y = vector.y;
		unitVector.z = vector.z;
	}

	return unitVector;
}
// ==============================================================

// Guarda el framebuffer en una imagen ppm
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
	for (i = Vres-1; i >= 0; i--){
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
// ==============================================================

// Ilumination ==================================================

//Calcula el factor de atenuacion de una luz
long double getAttenuationFactor(struct Light light, long double distance){
	long double value = 1 / (light.c1 + (light.c2 * distance) + (light.c3 * pow(distance, 2)));
	return min(1.0, value);
}

//Regresa el color base de un objeto al que se le aplica la intensidad de iluminacion difusa
struct Color difusseColor(long double I, struct Color color){
	struct Color newColor;
	
	newColor.r = I * color.r;
	newColor.g = I * color.g;
	newColor.b = I * color.b;

	return newColor;
}

//Regresa el color difuso co reflexion especular
struct Color specularHighlight(long double E, struct Color color){
	struct Color newColor;

	newColor.r = color.r + (E * (1 - color.r));
	newColor.g = color.g + (E * (1 - color.g));
	newColor.b = color.b + (E * (1 - color.b));

	return newColor;
}
// ==============================================================

//RAY TRACER

//Esferas: =======================================================
struct Intersection sphereIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	long double t, t1, t2;

	long double Xdif = anchor.x - object.Xc;
	long double Ydif = anchor.y - object.Yc;
	long double Zdif = anchor.z - object.Zc;
	struct Intersection tempIntersect;
	tempIntersect.null = 0;
	
	long double B = 2 * ((direction.x * Xdif) + (direction.y * Ydif) + (direction.z * Zdif));
	long double C = pow(Xdif, 2) + pow(Ydif, 2) + pow(Zdif, 2) - pow(object.other,2);

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
				tempIntersect.null = 1;
			}
		}
		
		if (tempIntersect.null != 1) {
			tempIntersect.distance = t;
			tempIntersect.object = object;
			tempIntersect.Xi = anchor.x + (t * direction.x);
			tempIntersect.Yi = anchor.y + (t * direction.y);
			tempIntersect.Zi = anchor.z + (t * direction.z);
	
		}
		
		return tempIntersect;
	}else{
		tempIntersect.null = 1;
		return tempIntersect;
	}
}

struct Vector sphereNormal(struct Object object, struct Vector vector){
	struct Vector normal;
	
	normal.x = vector.x - object.Xc;
	normal.y = vector.y - object.Yc;
	normal.z = vector.z - object.Zc;

	return normal;
}
// ===============================================================

//Polígonos: =====================================================
int getSign(long double v){
	if(v >= 0){ return 1; }
	else{ return 0; }
}

struct Intersection polygonIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	
	long double numerator = -((anchor.x * object.Xc) + (anchor.y * object.Yc) + (anchor.z * object.Zc) + object.other);
	long double denominator = (direction.x * object.Xc) + (direction.y * object.Yc) + (direction.z * object.Zc);

	struct Intersection tempIntersect;
	tempIntersect.null = 0;

	if(denominator == 0){
		tempIntersect.null = 1;
		return tempIntersect;

	}else{
		long double t = numerator / denominator;
		tempIntersect.distance = t;
		tempIntersect.object = object;
		tempIntersect.Xi = anchor.x + (t * direction.x);
		tempIntersect.Yi = anchor.y + (t * direction.y);
		tempIntersect.Zi = anchor.z + (t * direction.z);

		long double maxA_B = max(fabs(object.Xc), fabs(object.Yc)); //maximo entre A y B
		long double maxA_B_C = max(maxA_B, fabs(object.Zc)); //maximo entre los tres
		long double u, v;
		
		if(maxA_B_C == fabs(object.Xc)){
			//A es maximo
			//Aplasto X
			u = tempIntersect.Zi;
			v = tempIntersect.Yi;
		}else if(maxA_B_C == fabs(object.Yc) ){
			//B es maximo
			//Aplasto Y
			u = tempIntersect.Xi;
			v = tempIntersect.Zi;
		}else if(maxA_B_C == fabs(object.Zc)){
			//C es maximo
			//Aplasto Z
			u = tempIntersect.Xi;
			v = tempIntersect.Yi;
		}

		//PENDING: hacer lo de revisar con intersecciones 2D
		int NC = 0;
		int NV = object.pointAmount;

		for(int i = 0; i < NV; i++){
			object.points2D[i].u = object.points2D[i].u - u;
			object.points2D[i].v = object.points2D[i].v - v;
		}

		//printf("%i\n", object.pointAmount);
		int SH = getSign(object.points2D[0].v);
		int NSH;
		int a = 0;
		int b = (a+1)%NV;
		for (a = 0; a < NV-1; a++){
			NSH = getSign(object.points2D[b].v);
			
			if(SH != NSH){
				if(object.points2D[a].u > 0 && object.points2D[b].u > 0){
					NC++;
				}else if(object.points2D[a].u > 0 || object.points2D[b].u > 0){
					long double N = (object.points2D[b].u - object.points2D[a].u);
					long double D = (object.points2D[b].v - object.points2D[a].v);
					if(D != 0){
						if(object.points2D[a].u - ((object.points2D[a].v * N)/D) > 0){
							NC++;
						}	
					}
				}
			}
			SH = NSH;
			b++;
		}

		if (NC%2 == 0){
			for(int i = 0; i < NV; i++){
				object.points2D[i].u = object.points2D[i].u + u;
				object.points2D[i].v = object.points2D[i].v + v;
			}
			tempIntersect.null = 1;
			return tempIntersect;
		}else{
			for(int i = 0; i < NV; i++){
				object.points2D[i].u = object.points2D[i].u + u;
				object.points2D[i].v = object.points2D[i].v + v;
			}
			return tempIntersect;
		}
	}
}

struct Vector polygonNormal(struct Object object, struct Vector vector){
	
	struct Point3D point0 = object.points3D[0];
	struct Point3D point1 = object.points3D[1];
	struct Point3D point2 = object.points3D[2];

	struct Vector vector1 = {point1.x - point0.x, point1.y - point0.y, point1.z - point0.z};
	struct Vector vector2 = {point2.x - point1.x, point2.y - point1.y, point2.z - point1.z};

	struct Vector normal = crossProduct(vector1, vector2);
	normal = normalize(normal);

	return normal;
}
// ===============================================================

//Cilindros: =====================================================
struct Intersection cilinderIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	struct Intersection tempIntersect;
	tempIntersect.null = 0;

	long double xo = object.Xc;
	long double yo = object.Yc;
	long double zo = object.Zc;

	long double xq = object.directionVector.x;
	long double yq = object.directionVector.y;
	long double zq = object.directionVector.z;

	long double xd = direction.x;
	long double yd = direction.y;
	long double zd = direction.z;

	long double xe = anchor.x;
	long double ye = anchor.y;
	long double ze = anchor.z;

	long double radius = object.other;

	long double coef1 = xd*xq*xq + yd*yq*xq + zd*zq*xq - xd;
	long double coef2 = xd*xq*yq + yd*yq*yq + zd*zq*yq - yd;
	long double coef3 = xd*xq*zq + yd*yq*zq + zd*zq*zq - zd;

	long double coef4 = xo + xe*xq*xq - xo*xq*xq + ye*yq*xq - yo*yq*xq + ze*zq*xq - zo*zq*xq - xe;
	long double coef5 = yo + xe*xq*yq - xo*xq*yq + ye*yq*yq - yo*yq*yq + ze*zq*yq - zo*zq*yq - ye;
	long double coef6 = zo + xe*xq*zq - xo*xq*zq + ye*yq*zq - yo*yq*zq + ze*zq*zq - zo*zq*zq - ze;

	long double A = pow (coef1,2) + pow (coef2,2) + pow (coef3,2);

	long double B =  2*(coef1*coef4 + coef2*coef5 + coef3*coef6);

	long double C = pow (coef4, 2) + pow (coef5 ,2) + pow (coef6 ,2) - pow(radius,2);

	long double discriminant = pow(B, 2) - (4 * A * C);

	long double t, t1, t2;

	if(discriminant >= 0){
		long double root = sqrt(discriminant);

		B *= -1;

		t1 = (B + root)/(2*A);
		t2 = (B - root)/(2*A);

		long double Xi;
		long double Yi;
		long double Zi;
		long double d1 = object.D1;
		long double d2 = object.D2;
		if(t1 > e){
			if(t2 > e){
				t = min(t1, t2);
				Xi = xe + (t * xd);
				Yi = ye + (t * yd);
				Zi = ze + (t * zd);

				if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){ 
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
				}else{
					t = max(t1,t2);
					Xi = xe + (t * xd);
					Yi = ye + (t * yd);
					Zi = ze + (t * zd);

					if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){ 
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
					}else{
						tempIntersect.null = 1;
						return tempIntersect;
					}
				}
			}else{
				t = t1;
				Xi = xe + (t * xd);
				Yi = ye + (t * yd);
				Zi = ze + (t * zd);

				if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
				}else{
					tempIntersect.null = 1;
					return tempIntersect;
				}
			}
		}else{
			if(t2 > e){
				t = t2;
				Xi = xe + (t * xd);
				Yi = ye + (t * yd);
				Zi = ze + (t * zd);
				if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){ //Pertenece al cilindro finito?
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
				}else{
					tempIntersect.null = 1;
					return tempIntersect;
				}
			}else{
				tempIntersect.null = 1;
				return tempIntersect;
			}
		}	
	}else{
		tempIntersect.null = 1;
		return tempIntersect;
	}
}

struct Vector cilinderNormal(struct Object object, struct Vector intersectionPoint){
	
	struct Vector normalCilinder;
	long double x = intersectionPoint.x;
	long double y = intersectionPoint.y; 
	long double z = intersectionPoint.z; 

	long double xo = object.Xc;
	long double yo = object.Yc;
	long double zo = object.Zc;

	long double xq = object.directionVector.x;
	long double yq = object.directionVector.y;
	long double zq = object.directionVector.z;

	long double parethesesWithXq = ((x-xo)*xq+(y-yo)*yq+(z-zo)*zq)*xq;
	long double parethesesWithYq = ((x-xo)*xq+(y-yo)*yq+(z-zo)*zq)*yq;
	long double parethesesWithZq = ((x-xo)*xq+(y-yo)*yq+(z-zo)*zq)*zq;

	normalCilinder.x = 
					2*(xo + parethesesWithXq - x)*(pow(xq,2)-1)+
					2*(yo + parethesesWithYq - y)*(yq*xq)+
					2*(zo + parethesesWithZq - z)*(zq*xq) -0;

	normalCilinder.y = 
					2*(xo + parethesesWithXq - x)*(xq*yq)+
					2*(yo + parethesesWithYq - y)*(pow(yq,2)-1)+
					2*(zo + parethesesWithZq - z)*(zq*yq) -0;

	normalCilinder.z = 
					2*(xo + parethesesWithXq - x)*(xq*zq)+
					2*(yo + parethesesWithYq - y)*(yq*zq)+
					2*(zo + parethesesWithZq - z)*(pow(zq,2)-1) -0;

	normalCilinder = normalize(normalCilinder);

	return normalCilinder;
}
// ===============================================================

//Conos: =========================================================
struct Intersection coneIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	struct Intersection tempIntersect;
	tempIntersect.null = 0;

	long double xo = object.Xc;
	long double yo = object.Yc;
	long double zo = object.Zc;

	long double xq = object.directionVector.x;
	long double yq = object.directionVector.y;
	long double zq = object.directionVector.z;

	long double xd = direction.x;
	long double yd = direction.y;
	long double zd = direction.z;

	long double k1 = object.K1;
	long double k2 = object.K2;

	long double xe = anchor.x;
	long double ye = anchor.y;
	long double ze = anchor.z;

	//Los mismos del cilindro
	long double coef1 = xd*xq*xq + yd*yq*xq + zd*zq*xq - xd;
	long double coef2 = xd*xq*yq + yd*yq*yq + zd*zq*yq - yd;
	long double coef3 = xd*xq*zq + yd*yq*zq + zd*zq*zq - zd;

	long double coef4 = xo + xe*xq*xq - xo*xq*xq + ye*yq*xq - yo*yq*xq + ze*zq*xq - zo*zq*xq - xe;
	long double coef5 = yo + xe*xq*yq - xo*xq*yq + ye*yq*yq - yo*yq*yq + ze*zq*yq - zo*zq*yq - ye;
	long double coef6 = zo + xe*xq*zq - xo*xq*zq + ye*yq*zq - yo*yq*zq + ze*zq*zq - zo*zq*zq - ze;

	//En lugar del radio usa una combinación de estos
	long double coefk = pow(k2/k1,2);
	long double coef7 = xd*xq + yd*yq + zd*zq;
	long double coef8 = xe*xq - xo*xq + ye*yq - yo*yq + ze*zq - zo*zq;

	long double A = pow(coef1,2) + pow(coef2,2) + pow(coef3,2) - (coefk * pow(coef7,2));
	long double B =  2*((coef1*coef4 + coef2*coef5 + coef3*coef6) - (coefk* coef7*coef8));
	long double C = pow (coef4, 2) + pow (coef5 ,2) + pow (coef6 ,2) - (coefk*pow(coef8,2));

	long double discriminant = pow(B, 2) - (4 * A * C);
	long double t, t1, t2;

	if(discriminant >= 0){
		long double root = sqrt(discriminant);

		B *= -1;

		t1 = (B + root)/(2*A);
		t2 = (B - root)/(2*A);

		long double Xi;
		long double Yi;
		long double Zi;
		long double d1 = object.D1;
		long double d2 = object.D2;
		if(t1 > e){
			if(t2 > e){
				t = min(t1, t2);
				Xi = xe + (t * xd);
				Yi = ye + (t * yd);
				Zi = ze + (t * zd);

				if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){ 
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
				}else{
					t = max(t1,t2);
					Xi = xe + (t * xd);
					Yi = ye + (t * yd);
					Zi = ze + (t * zd);

					if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){ 
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
					}else{
						tempIntersect.null = 1;
						return tempIntersect;
					}
				}
			}else{
				t = t1;
				Xi = xe + (t * xd);
				Yi = ye + (t * yd);
				Zi = ze + (t * zd);

				if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
				}else{
					tempIntersect.null = 1;
					return tempIntersect;
				}
			}
		}else{
			if(t2 > e){
				t = t2;
				Xi = xe + (t * xd);
				Yi = ye + (t * yd);
				Zi = ze + (t * zd);
				if ( d2 >= ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) && ((Xi-xo)*xq + (Yi-yo)*yq + (Zi-zo)*zq) >= d1){ //Pertenece al cilindro finito?
						tempIntersect.Xi = Xi;
						tempIntersect.Yi = Yi;
						tempIntersect.Zi = Zi;
						tempIntersect.distance = t;
						tempIntersect.object = object;
						return tempIntersect;
				}else{
					tempIntersect.null = 1;
					return tempIntersect;
				}
			}else{
				tempIntersect.null = 1;
				return tempIntersect;
			}
		}	
	}else{
		tempIntersect.null = 1;
		return tempIntersect;
	}
}

struct Vector coneNormal(struct Object object, struct Vector intersectionPoint){

	struct Vector normalCone;

	long double x = intersectionPoint.x;
	long double y = intersectionPoint.y; 
	long double z = intersectionPoint.z; 

	long double xo = object.Xc;
	long double yo = object.Yc;
	long double zo = object.Zc;

	long double xq = object.directionVector.x;
	long double yq = object.directionVector.y;
	long double zq = object.directionVector.z;

	long double k1 = object.K1;
	long double k2 = object.K2;

	long double parethesesWithXq = ((x-xo)*xq+(y-yo)*yq+(z-zo)*zq)*xq;
	long double parethesesWithYq = ((x-xo)*xq+(y-yo)*yq+(z-zo)*zq)*yq;
	long double parethesesWithZq = ((x-xo)*xq+(y-yo)*yq+(z-zo)*zq)*zq;


	long double factCoefsK = ( (yq*(y-yo) )+(zq*(z-zo))+(xq*(-xo+x)));

	long double lastFactorDerivedX = (((2*(pow(k2,2))*xq)*factCoefsK)/(pow(k1,2)));
	long double lastFactorDerivedY = (((2*(pow(k2,2))*yq)*factCoefsK)/(pow(k1,2)));
	long double lastFactorDerivedZ = (((2*(pow(k2,2))*zq)*factCoefsK)/(pow(k1,2)));

	normalCone.x = 2*(xo + parethesesWithXq - x)*(pow(xq,2)-1)+
					2*(yo + parethesesWithYq - y)*(xq*yq)+
					2*(zo + parethesesWithZq - z)*(xq*zq) -
					lastFactorDerivedX;

	normalCone.y = 2*(xo + parethesesWithXq - x)*(xq*yq)+
					2*(yo + parethesesWithYq - y)*(pow(yq,2)-1)+
					2*(zo + parethesesWithZq - z)*(yq*zq) -
					lastFactorDerivedY;

	normalCone.z = 2*(xo + parethesesWithXq - x)*(zq*xq)+
					2*(yo + parethesesWithYq - y)*(yq*zq)+
					2*(zo + parethesesWithZq - z)*(pow(zq,2)-1) 
					- lastFactorDerivedZ;

	normalCone = normalize(normalCone);

	return normalCone;
}
// ===============================================================

//Planos: ========================================================
long double whatsTheD(struct Object object){
	
	struct Point3D point = object.points3D[0];
	long double theD = -((object.Xc * point.x) + (object.Yc * point.y) + (object.Zc * point.z));
	
	return theD;
}

struct Object getABCD(struct Object object){
	struct Vector a;
	struct Vector normal = polygonNormal(object, a);

	object.Xc = normal.x; //A de la ecuacion del plano
	object.Yc = normal.y; //B de la ecuacion del plano
	object.Zc = normal.z; //C de la ecuacion del plano
	object.other = whatsTheD(object); //D de la ecuacion del plano

	long double L = getNorm(normal);
	object.Xc /= L;
	object.Yc /= L;
	object.Zc /= L;
	object.other /= L;	

	return object; //Probar esto, sino regresar un array con ABCD
}
//================================================================

struct Intersection getFirstIntersection(struct Vector anchor, struct Vector direction){
	int k;
	long double tmin;
	struct Intersection intersection;
	struct Intersection tempIntersection;

	tmin = 100000;
	tempIntersection.null = 1;
	intersection.distance = -1;

	int objectsAmount = numberObjects;
	for(k = 0; k < objectsAmount; k++){
		tempIntersection = Objects[k].intersectionFuncion(anchor, direction, Objects[k]);
			
		if(tempIntersection.null != 1 && tempIntersection.distance > e && tempIntersection.distance < tmin){
			
			tmin = tempIntersection.distance;
			intersection = tempIntersection;
		}
		tempIntersection.null = 1;
	}
	return intersection;
}

struct Color getColor(struct Vector anchor, struct Vector direction){
	struct Color color;
	struct Intersection intersection;
	struct Intersection *tempIntersection;
	
	intersection = getFirstIntersection(anchor, direction);

	if(intersection.distance == -1){
		color = background;
	}else{
	
		int k;
		int lightsAmount = numberLights;
		
		struct Object Q = intersection.object;

		struct Vector L;
		struct Vector intersectVector = {intersection.Xi, intersection.Yi, intersection.Zi};
		struct Vector N = normalize(Q.normalVector(Q, intersectVector));
		struct Vector R;

		if(pointProduct(N, direction) > 0){
			N.x *= -1;
			N.y *= -1;
			N.z *= -1;
		}

		long double Fatt;
		long double I = 0.0;
		long double E = 0.0;

		for(k = 0; k < numberLights; k++){
			struct Intersection obstacle;
			struct Vector light = {Lights[k].Xp - intersection.Xi, Lights[k].Yp - intersection.Yi, Lights[k].Zp - intersection.Zi};
			L = normalize(light);
			
			long double a = pointProduct(N, L);
	
			R.x = (2 * N.x * a) - L.x;
			R.y = (2 * N.y * a) - L.y;
			R.z = (2 * N.z * a) - L.z;
			R = normalize(R);

			//rec = 1;
			obstacle = getFirstIntersection(intersectVector, L);
			//rec = 0;

			long double distanceToLight = getNorm(light);

			if(obstacle.distance < e || (obstacle.distance > e && obstacle.distance > distanceToLight)) {
				long double pp = pointProduct(N, L);
				Fatt = getAttenuationFactor(Lights[k], distanceToLight);
				if(pp > 0.0){
					I = I + (pp * Q.Kd * Fatt * Lights[k].Ip);
				}
				long double pp2 = pointProduct(R, V);
				if(pp2 > 0.0){
					E = E + (pow(pp2, Q.Kn) * Q.Ks * Lights[k].Ip * Fatt);
				}
			}
		}
		I = I + Ia * Q.Ka;
		I = min(1.0, I);
		color = difusseColor(I, Q.color);

		E = min(1.0, E);
		color = specularHighlight(E, color);
	}
	return (color);
}

// Lectura de archivos ===========================================
void createObjectFromData(long double *data, int whichObjectCreate, int quantityData){
	/*whichObjectCreate indica qué objeto crear
		Si está en 0, altera valores de la escena
		SI está en 1, crea luces
		Si está en 2, crea esferas
		Si está en 3, crea polígonos, 
		SI está en 4, crea cilindros
		SI está en 5, crea conos.
	data es el arreglo de valores del objeto a crear. Se asume que estará completo.*/
	
	switch(whichObjectCreate){
		case 0: {
			if (debug == 1) {
				printf("Insertando datos de escena \n");
				printf("Iluminación ambiente: %LF \n", data[0]);
				printf("Plano de proyección (Xmin, Ymin) (Xmax, Ymax) : (%LF, %LF) (%LF, %LF) \n", data[1], data[2],data[3],data[4]);
				printf("Resolución:  %LFx%LF \n", data[5], data[6]);
				printf("Epsilon %LF \n", data[7]);
				printf("Ojo: (%LF, %LF, %LF) \n", data[8],data[9],data[10]);
				printf("Color background: (%LF, %LF, %LF) \n", data[11],data[12],data[13]);
			}

			Ia = data[0];
			Xmin = data[1];
			Ymin = data[2];
			Xmax = data[3];
			Ymax = data[4];
			Hres = data[5];
			Vres = data[6];
			e = data[7];

			eye.x = data[8];
			eye.y = data[9];
			eye.z = data[10];

			background.r = data[11];
			background.g = data[12];
			background.b = data[13];

			Framebuffer[Hres][Vres];
			Framebuffer = (struct Color **)malloc(Vres * sizeof(struct Color*));
			for (int i = 0; i<Vres; i++){
				Framebuffer[i] = (struct Color *)malloc(Hres*sizeof(struct Color));
			}
			
			//printf("Scene data read correctly.\n");
			return;
			}
		case 1: {
			if (debug == 1) {
			
				printf("Insertando Luz...\n");
	
				printf("Pos luz (%LF, %LF, %LF) \n", data[0],data[1],data[2]);
				printf("c1: %LF, c2: %LF, c3 %LF \n", data[3],data[4],data[5]);
				printf("Ip luz: %LF \n", data[6]);
			}

			struct Object polygon;
			struct Color colorPolygon;


			struct Light luz;
			luz.Xp=data[0];
			luz.Yp=data[1];
			luz.Zp=data[2];
			luz.c1=data[3];
			luz.c2=data[4];
			luz.c3=data[5];
			luz.Ip=data[6];

			Lights[lightIndex]=luz;
			lightIndex++;

			// printf("Light processed. \n \n");
			return;
			}
		case 2: {
			if (debug == 1) {
				printf("Insertando Esfera...");
				printf("Pos esfera (%LF, %LF, %LF) \n", data[0],data[1],data[2]);
				printf("Radio esfera: %LF \n", data[3]);
				printf("Esfera Kd: %LF \n", data[3]);
				printf("Esfera Ka: %LF \n", data[4]);
				printf("Esfera Kn: %LF \n", data[5]);
				printf("Esfera Ks: %LF \n", data[6]);
				printf("Color esfera (%LF, %LF, %LF) \n", data[8],data[9],data[10]);
				
			}
			
			struct Object polygon;
			struct Color colorPolygon;
			
			struct Object esfera;
			esfera.Xc=data[0];
			esfera.Yc=data[1];
			esfera.Zc=data[2];
			esfera.other=data[3]; //Radio
			esfera.Kd=data[4]; 
			esfera.Ka=data[5];
			esfera.Kn=data[6];
			esfera.Ks=data[7];
			esfera.normalVector = sphereNormal;
			esfera.intersectionFuncion = sphereIntersection;
			struct Color colorSphere;
			colorSphere.r = data[8];
			colorSphere.g = data[9];
			colorSphere.b = data[10];
			esfera.color=colorSphere;
			Objects[objectIndex]=esfera;
			objectIndex++;
			//printf("Sphere processed \n \n");
			return;
			}
		case 3: {
			if (debug == 1){
				printf("Insertando polígono...");
				printf("Color polígono (%LF, %LF, %LF) \n", data[0],data[1],data[2]);
				printf("Poligono Kd: %LF \n", data[3]);
				printf("Poligono Ka: %LF \n", data[4]);
				printf("Poligono Kn: %LF \n", data[5]);
				printf("Poligono Ks: %LF \n", data[6]);
			}

			//7 Elementos adicionales a los vertices
			//Creo un objeto temporal
			int vertexPolygonIndex = 0;
			struct Object temp;
			int numVertexesPolygon = (quantityData-7) / 3;
			temp.points3D = malloc(sizeof(struct Point3D)*3);
			for (int i =0; i+7 < quantityData;){
				if(vertexPolygonIndex==3){
					break;
				}
				struct Point3D vertex;
				//printf("Vertice (%LF, %LF, %Lf) \n", data[7+i],  data[7+i+1],  data[7+i+2] );
				vertex.x = data[7+i];
				i++;
				vertex.y = data[7+i];
				i++;
				vertex.z = data[7+i];
				i++;

				temp.points3D[vertexPolygonIndex]=vertex;
				vertexPolygonIndex++;
			}

			struct Object polygon;
			vertexPolygonIndex = 0;
			//printf("numVertexesPolygon %i \n", numVertexesPolygon);
			polygon = getABCD(temp);

			if (debug == 1) {
				printf("A del poligono %LF\n", polygon.Xc);
				printf("B del poligono %LF\n", polygon.Yc);
				printf("C del poligono %LF\n", polygon.Zc);
				printf("D del poligono %LF\n", polygon.other);
			}

			polygon.points3D = malloc(sizeof(struct Point3D)*numVertexesPolygon);
			polygon.points2D = malloc(sizeof(struct Point2D)*numVertexesPolygon);
			
			struct Color colorPolygon;
			colorPolygon.r = data[0];
			colorPolygon.g = data[1];
			colorPolygon.b = data[2];
			polygon.color =  colorPolygon;
			polygon.Kd = data[3];
			polygon.Ka = data[4];
			polygon.Kn = data[5];
			polygon.Ks = data[6];

			polygon.pointAmount = numVertexesPolygon;
			polygon.normalVector = polygonNormal;
			polygon.intersectionFuncion = polygonIntersection;
			
			long double u;
			long double v;
			
			//printf("abs A : %lf abs B: %lf abs C: %lf \n",fabs(polygon.Xc), fabs(polygon.Yc), fabs(polygon.Zc));
			
			long double maxA_B = max(fabs(polygon.Xc), fabs(polygon.Yc)); //maximo entre A y B
			long double maxA_B_C = max(maxA_B, fabs(polygon.Zc)); //maximo entre los tres
			
			for (int i =0; i+7 < quantityData;){
				struct Point3D vertex;
				
				vertex.x = data[7+i];
				i++;
				vertex.y = data[7+i];
				i++;
				vertex.z = data[7+i];
				i++;
				
				struct Point2D squashedVertex;

				if(maxA_B_C == fabs(polygon.Xc)){ //A es maximo
					u = vertex.z;
					v = vertex.y;
				}

				else if(maxA_B_C == fabs(polygon.Yc) ){ //B es maximo
					u = vertex.x;
					v = vertex.z;
				}

				else if(maxA_B_C == fabs(polygon.Zc)){ //C es maximo
					u = vertex.x;
					v = vertex.y;
				} 

				squashedVertex.u = u;
				squashedVertex.v = v;
				
				polygon.points3D[vertexPolygonIndex]=vertex;
				
				polygon.points2D[vertexPolygonIndex]=squashedVertex;
				vertexPolygonIndex++;
			}

			Objects[objectIndex] = polygon;
			objectIndex++;
			
			//printf("Polygon processed. \n \n");
			return;
		}
		case 4: { //Cilindros
			if (debug == 1) {
				printf("Insertando cilindro...");

				printf("Ancla: (%LF, %LF, %LF) \n", data[0], data[1],data[2]);
				printf("Vector: (%LF, %LF, %LF) \n", data[3], data[4],data[5]);
				printf("Cilindro Radio: %LF \n", data[6]);
				printf("Cilindro d1: %LF Cilindro d2: %LF \n", data[7],data[8]);
				printf("Cilindro Kd: %LF \n", data[9]);
				printf("Cilindro Ka: %LF \n", data[10]);
				printf("Cilindro Kn: %LF \n", data[11]);
				printf("Cilindro Ks: %LF \n", data[12]);
				printf("RGB Cilindro: (%LF, %LF, %LF) \n", data[13], data[14],data[15]);
	
			}
			
			struct Object cilinder;
			cilinder.Xc = data[0];
			cilinder.Yc = data[1];
			cilinder.Zc = data[2];

			struct Vector cilinderVector; 
			cilinderVector.x = data[3];
			cilinderVector.y = data[4];
			cilinderVector.z = data[5];
			cilinderVector = normalize(cilinderVector);
			//printf("Vector normalizado del: (%LF, %LF, %LF) \n", cilinderVector.x, cilinderVector.y,cilinderVector.z);
			cilinder.directionVector = cilinderVector;

			cilinder.other = data[6];
			cilinder.D1 = data[7];
			cilinder.D2 = data[8];
			cilinder.Kd = data[9];
			cilinder.Ka = data[10];
			cilinder.Kn = data[11];
			cilinder.Ks = data[12];
			cilinder.normalVector = cilinderNormal;
			cilinder.intersectionFuncion = cilinderIntersection;

			struct Color cilinderColor;
			cilinderColor.r = data[13];
			cilinderColor.g = data[14];
			cilinderColor.b = data[15];
			cilinder.color = cilinderColor;
			
			Objects[objectIndex] = cilinder;
			objectIndex++;
			//printf("Cylinder processed.\n \n");
			return;
		}
		case 5: { //Conos
			if (debug == 1) {
				printf("Insertando cono...");
				printf("Ancla: (%LF, %LF, %LF) \n", data[0], data[1],data[2]);
				printf("Vector: (%LF, %LF, %LF) \n", data[3], data[4],data[5]);
				printf("Cono k1: %LF COno k2: %LF \n", data[6],data[7]);
				printf("Cono d1: %LF COno d2: %LF \n", data[8],data[9]);
				printf("Cono Kd: %LF \n", data[10]);
				printf("Cono Ka: %LF \n", data[11]);
				printf("Cono Kn: %LF \n", data[12]);
				printf("Cono Ks: %LF \n", data[13]);
				printf("RGB Cono: (%LF, %LF, %LF) \n", data[14], data[15],data[16]);
				
			}
				
			struct Object cone;
			cone.Xc = data[0];
			cone.Yc = data[1];
			cone.Zc = data[2];
			
			struct Vector coneVector; 
			coneVector.x = data[3];
			coneVector.y = data[4];
			coneVector.z = data[5];
			coneVector = normalize(coneVector);
			
			cone.directionVector = coneVector;

			cone.K1 = data[6];
			cone.K2 = data[7];
			cone.D1 = data[8];
			cone.D2 = data[9];
			cone.Kd = data[10];
			cone.Ka = data[11];
			cone.Kn = data[12];
			cone.Ks = data[13];
			cone.intersectionFuncion = coneIntersection;
			cone.normalVector = coneNormal;

			struct Color coneColor;
			coneColor.r = data[14];
			coneColor.g = data[15];
			coneColor.b = data[16];
			cone.color = coneColor;
			
			Objects[objectIndex] = cone;
			objectIndex++;
			//printf("Cone processed.\n \n");
			return;
		}
	}
}

long double *obtainPointFromString(char stringPoint[]){
	/*Devuelve los tres valores long double de un punto tridimensional 
	a partir de la forma Xp,Yp,Zp. RECORDAR UTILIZAR free(valorDevuelto) tras
	 usarla.*/
	char *token;
	char *search = "=";
	long double numericValue;
	// Token will point to "Eye ".
	token = strtok(stringPoint, search);
	// Token will point to " 0.4, 0.5, 0.7".
	token = strtok(NULL, search);
	char *pch;
	long double *pointDimensions = malloc(sizeof(long double) * 3);
	int currentDimension=0;
	pch = strtok (token,",");
	while (pch != NULL)
	{
		sscanf(pch, "%LF", &pointDimensions[currentDimension]);
		pch = strtok (NULL, ",");
		currentDimension++;
	}
	return pointDimensions;
}

long double obtainSingleValueFromLine(char line[]){
	/*Devuelve el valor flotante leído de una línea, realizando el proceso de
	separación entre dicho valor y el valor a asignarse*/
	/*Ejemplo Kn = 0.3 */
	char *token;
	char *search = "=";
	long double numericValue;
	// Token will point to "Kn ".
	token = strtok(line, search);
	// Token will point to "0.3".
	token = strtok(NULL, search);
	sscanf(token, "%LF", &numericValue);
	return numericValue;
}

long double *readValueFromLine(int state, int *counterValueSegment, char* lineRead, int *numberValuesRead){
	/* RECIBE POR REFERENCIA EL VALOR DE counterValueSegment.
	Decide que hacer con cada linea de datos dependiendo de en cuál estado se encuentre
	el lector y en qué valor de dicho segmento se encuentra. RETORNA MÁXIMO 20 VALORES. 
	MODIFICAR PARA POLÍGONOS.  */ 
	//long double *values = malloc(sizeof(long double)*20);
	long double *values;
	//RECORDAR LIBERARLA TRAS USO.
	switch(state){
		case 0: //Escena
			if ((*counterValueSegment) >= 0 && (*counterValueSegment) <= 7){ //Lee Ilum. AMbiente
				values = malloc(sizeof(long double));
				values[0] = obtainSingleValueFromLine(lineRead);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
			}else if ((*counterValueSegment) >= 8 && (*counterValueSegment) <= 9){

				long double *point = obtainPointFromString(lineRead);
				values = malloc(sizeof(long double)*3);
				values[0] = point[0];
				values[1] = point[1];
				values[2] = point[2];
				*numberValuesRead = 3;
				free (point);
				if((*counterValueSegment) == 9){//ES Ip, último valor
					(*counterValueSegment) = 0;
				}else{
					(*counterValueSegment)++;
				}
				//printf("%LF %LF %LF \n", values[0], values[1],values[2]);
				return values;
			}
		case 1: //Luces
				if ((*counterValueSegment)== 0){
						long double *positionLight = obtainPointFromString(lineRead);
						values = malloc(sizeof(long double)*3);
						values[0] = positionLight[0];
						values[1] = positionLight[1];
						values[2] = positionLight[2];
						//memcpy(values, positionLight, 3);
						//printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
						(*counterValueSegment)++;
						*numberValuesRead = 3;
						free (positionLight);
						//printf("%LF %LF %LF \n", values[0], values[1],values[2]);
						return values;
				}else if((*counterValueSegment) >= 1 && (*counterValueSegment <= 4)){ 

					values = malloc(sizeof(long double));
					values[0] = obtainSingleValueFromLine(lineRead);
					if((*counterValueSegment) == 4){//ES Ip, último valor
						(*counterValueSegment) = 0;
					}else{
						(*counterValueSegment)++;
					}
					*numberValuesRead = 1;
					return values;
				}
		case 2: //ESferas
			if ((*counterValueSegment)== 0 || (*counterValueSegment) ==6){ //Nombre o posición de la esfera
				long double *positionSphere = obtainPointFromString(lineRead);
				values = malloc(sizeof(long double)*3);
				values[0] = positionSphere[0];
				values[1] = positionSphere[1];
				values[2] = positionSphere[2];

				free (positionSphere);
				*numberValuesRead = 3;
				if((*counterValueSegment)==6){
					(*counterValueSegment)=0;
				}else{
					(*counterValueSegment)++;
				}
				return values;
			}else if((*counterValueSegment) >= 1 && (*counterValueSegment )<= 5 ){ 
				values = malloc(sizeof(long double));
				values[0] = obtainSingleValueFromLine(lineRead);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
			}
		case 3: //Poligonos
			if ((*counterValueSegment)== 0){ //Nombre o RGB del poligono
				values = malloc(sizeof(long double)*3);	
				long double *rgbColors = obtainPointFromString(lineRead);
				values[0] = rgbColors[0];
				values[1] = rgbColors[1];
				values[2] = rgbColors[2];
				free (rgbColors);
				*numberValuesRead = 3;
				(*counterValueSegment)++;
				//printf("counterValueSegment tras asignar RGB %i \n",(*counterValueSegment));
				return values;
			}else if((*counterValueSegment) >= 1 && (*counterValueSegment) <= 4){ //Ks
				values = malloc(sizeof(long double));
				values[0] = obtainSingleValueFromLine(lineRead);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
			}else if((*counterValueSegment) == 5){ //Poligono Vertice
				if (strstr(lineRead,"END_Vertices")!=NULL){
					(*counterValueSegment)=0;
					return values;
				}else{
					values = malloc(sizeof(long double)*3);	
					long double *vertexPolygon = obtainPointFromString(lineRead);
					values[0] = vertexPolygon[0];
					values[1] = vertexPolygon[1];
					values[2] = vertexPolygon[2];
					free (vertexPolygon);
					*numberValuesRead = 3;
					//printf("counterValueSegment tras asignar RGB %i \n",(*counterValueSegment));
					return values;
				}
			}
		case 4: //Cilindros
			if(*counterValueSegment == 0 || *counterValueSegment == 1|| *counterValueSegment == 9){
				/*Lee el ancla, color del cilindro o el vector del cilindro. Tripleta*/
				long double *positionCilinder = obtainPointFromString(lineRead);
				values = malloc(sizeof(long double)*3);
				values[0] = positionCilinder[0];
				values[1] = positionCilinder[1];
				values[2] = positionCilinder[2];
				//memcpy(values, positionLight, 3);
				//printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
				if(*counterValueSegment == 9){
					(*counterValueSegment) = 0;
				}else{
					(*counterValueSegment)++;
				}
					
				*numberValuesRead = 3;
				free (positionCilinder);
				//printf("%LF %LF %LF \n", values[0], values[1],values[2]);
				return values;
			}else if(*counterValueSegment >= 2 && *counterValueSegment <= 8){
				values = malloc(sizeof(long double));
				values[0] = obtainSingleValueFromLine(lineRead);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
				/*Lee radio o d1 o d2 o Kd o Kd o Ka o Kn o Ks*/
			}
		case 5: {//Conos
			if(*counterValueSegment == 0 || *counterValueSegment == 1|| *counterValueSegment == 10){
				/*Lee el ancla, color del cono o el vector del cono. Tripleta*/
				long double *positionCone = obtainPointFromString(lineRead);
				values = malloc(sizeof(long double)*3);
				values[0] = positionCone[0];
				values[1] = positionCone[1];
				values[2] = positionCone[2];
				//memcpy(values, positionLight, 3);
				//printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
				if(*counterValueSegment == 10){
					(*counterValueSegment) = 0;
				}else{
					(*counterValueSegment)++;
				}
				*numberValuesRead = 3;
				free (positionCone);
				//printf("%LF %LF %LF \n", values[0], values[1],values[2]);
				return values;
			}else if(*counterValueSegment >= 2 && *counterValueSegment <= 9){
				values = malloc(sizeof(long double));
				values[0] = obtainSingleValueFromLine(lineRead);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
				/*Lee radio o k1 i k2 o d1 o d2 o Kd o Kd o Ka o Kn o Ks*/
			}
		}
	}
}
// ==============================================================

//Leer archivos con la escena ===================================
void getSceneObjects(){

	int i, j, c;
	int state = 0;
	/*  State indica qué segmento lee de la escena
		Si está en 0, lee variables de escena
		Si está en 1, lee luces
		Si está en 2, lee esferas
		Si está en 3, lee polígonos, 
		Si está en 4, lee cilindros
		Si está en 5, lee conos. */

	int counterValueSegment = 0; 
	// Contador que indica qué valor del objeto está leyendo del actual objeto de un segmento. 
	
	char temporalBuffer[100]; //Aquí se guardará lo leído cada línea
	long double *valuesRead;  //Se guarda los valores del objeto para finalmente crearlo.
	int indexValuesRead = 0;  //Pos de valuesRead 
	int currentTypeObjectReading = 1;
	
	/*  Si está en 1, crea luces
		Si está en 2, crea esferas
		Si está en 3, crea polígonos, 
		Si está en 4, crea cilindros
		Si está en 5, crea conos. */

	FILE* file; //archivo
	if (file = fopen(escenaFile, "r")){

		while (fgets(temporalBuffer, 100, file)!=NULL){ //Mientras el archivo siga teniendo algo

			if (temporalBuffer[0] == '\n'){
				continue;
			}
			if (strstr(temporalBuffer, "#")!=NULL){
				continue;
			}

			if (strstr(temporalBuffer, "Scene_Data")!=NULL){
			//Entra en state = 0. Escena
				state = 0;
				counterValueSegment = 0;
				indexValuesRead=0;
				valuesRead=NULL;
				valuesRead = malloc(sizeof(long double)*14);  
				currentTypeObjectReading=0;
				//printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer,"Light_Object")!=NULL){
			//Entra en state = 1. Luces
				state = 1;
				counterValueSegment = 0;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead = malloc(sizeof(long double)*7); //Las luces siempre seran 7 valores
				currentTypeObjectReading=1;
				//printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Sphere_Object")!=NULL){
			//Entra en state = 2. Esferas
				state = 2;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead=malloc(sizeof(long double)*11); //Esferas siempre serán 11 valores
				counterValueSegment = 0;
				currentTypeObjectReading=2;
				//printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Polygon_Object")!=NULL){
			//Entra en state = 3. Poligonos
				//printf("Leyó polígonos");
				state = 3;
				counterValueSegment = 0;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead=malloc(sizeof(long double)*200); //Polígonos max size 
				currentTypeObjectReading=3;
				//printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Cylinder_Object")!=NULL){
			//Entra en state = 4. Cilindros
				state = 4;
				counterValueSegment = 0;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead = malloc(sizeof(long double)*16);
				currentTypeObjectReading=4;
				//printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Cone_Object")!=NULL){
			//Entra en state = 5. Conos
				state = 5;
				counterValueSegment = 0;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead = malloc(sizeof(long double)*17);
				currentTypeObjectReading=5;
				//printf("%s",temporalBuffer);
				continue;
			}
			
			int numberValuesRead = 0;
			//printf("%s \n",temporalBuffer);
			//printf("state %i \n" ,state);
			//printf("counterValueSegment %i \n", counterValueSegment);
			long double *valuesReadTemp = readValueFromLine(state, &counterValueSegment, temporalBuffer, &numberValuesRead);
			if (valuesReadTemp == NULL){ //Se devolvió NULL
				//printf("readValueFromLine devolvió NULL \n");
				continue;
			}
			int i = 0;
			for (i = 0; i < numberValuesRead; i++){
				valuesRead[indexValuesRead+i] = valuesReadTemp[i];
			}
			indexValuesRead+=numberValuesRead;
			/*SI counterValueSegment vuelve como un 0, quiere decir que ya
			 se leyeron los datos de dicho objeto/luz, ergo, se crea.*/
			if (counterValueSegment == 0){
				createObjectFromData(valuesRead, currentTypeObjectReading, indexValuesRead);
				//printf("Sup \n");
				indexValuesRead=0;
			}
		}
		free(valuesRead);
	}
	fclose(file);	
	//Pendiente
}

void howManyObjectsLights(){
	/*Lee el número de objetos del archivo de texto para inicializar la 
	memoria de los arreglos globales de luces y objetos de manera
	dinámica*/
	char temporalBuffer[100]; //Aquí se guardará lo leído cada línea
	FILE* file; //archivo
	if (file = fopen(escenaFile, "r")){

		while (fgets(temporalBuffer, 100, file)!=NULL){//Mientras el archivo siga teniendo algo

			if (temporalBuffer[0] == '\n'){
				continue;
			}
			if (strstr(temporalBuffer, "#")!=NULL){
				continue;
			}
			if (strstr(temporalBuffer, "Light_Object")!=NULL){
				numberLights++;
				continue;
			}else if (strstr(temporalBuffer,"Sphere_Object")!=NULL){
				numberObjects++;
				continue;
			}else if (strstr(temporalBuffer, "Polygon_Object")!=NULL){
				numberObjects++;
				continue;
			}else if (strstr(temporalBuffer, "Cone_Object")!=NULL){
				numberObjects++;
				continue;
			}else if (strstr(temporalBuffer, "Cylinder_Object")!=NULL){
				numberObjects++;
				continue;
			}

		}
	}
	fclose(file);
	Objects = malloc(sizeof(struct Object)*numberObjects);
	Lights= malloc(sizeof(struct Light)*numberLights);
}
// ==============================================================

// OP main: =====================================================
int main(int argc, char *arcgv[]){
	howManyObjectsLights();
	printf("Lights: %i \n", numberLights);
	printf("Objects: %i \n", numberObjects);
	getSceneObjects();
	
	int i, j;

	long double L;
	long double Xw, Yw;
	long double Xd, Yd, Zd;

	struct Color color;
	
	struct Vector direction;
	
	long double Xdif = Xmax - Xmin;
	long double Ydif = Ymax - Ymin;

	Zd = -eye.z;

	printf("\nRay Tracing\n...\n...\n");
	for (i = 0; i < Vres; i++){
		Yw = (long double) ((i + (1/2)) * Ydif)/Vres + Ymin;
		Yd = Yw - eye.y;
		
		for (j = 0; j < Hres; j++){
			Xw = (long double) ((j + (1/2)) * Xdif)/Hres + Xmin;
			Xd = Xw - eye.x;

			L = sqrt(pow(Xd, 2) + pow(Yd, 2) + pow(Zd, 2));
		
			direction.x = Xd / L;
			direction.y = Yd / L;
			direction.z = Zd / L;

			direction = normalize(direction);

			V.x = -direction.x;
			V.y = -direction.y;
			V.z = -direction.z;

			V = normalize(V);

			color = getColor(eye, direction);

			Framebuffer[i][j] = color;
		}
	}
	saveFile();
	free(Objects);
	free(Lights);
	free(Framebuffer);
	printf("\nDONE.\n");
	

}
// ==============================================================
