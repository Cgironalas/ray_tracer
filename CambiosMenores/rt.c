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
static long double Ia = 0.2;
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
	long double Kn;
	long double Ks;
	long double other;
	struct Color color;
	struct Vector (*normalVector)();
	struct Intersection *(*intersectionFuncion)(struct Vector, struct Vector, struct Object);
};

struct Intersection {
	long double Xi;
	long double Yi;
	long double Zi;
	long double distance;
	struct Object object;
};

static struct Light *Lights;
static struct Object *Objects;
static int numberObjects = 0;
static int numberLights = 0;
static int lightIndex = 0;
static int objectIndex = 0;

static struct Vector V;
static struct Intersection tempIntersect;

static struct Vector eye = {200,200,-1500};
static struct Color Framebuffer[1008][1008];
static struct Color BACKGROUND = {0.3, 0.3, 0.3};

static char* escenaFile = "escena1.txt";


long double min(long double a, long double b){
    if(a < b) { return a; }
    else { return b; }
}
long double max(long double a, long double b){
	if(a > b) { return a; }
	else { return b; }
}








//DONE
//Producto punto entre dos vectores
long double pointProduct(struct Vector a, struct Vector b){
	long double pp = 0;

	pp += (a.x * b.x);
	pp += (a.y * b.y);
	pp += (a.z * b.z);
	
	return pp;
}

//DONE
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
struct Color difusseColor(long double I, struct Color color){
	struct Color newColor;
	
	newColor.r = I * color.r;
	newColor.g = I * color.g;
	newColor.b = I * color.b;

	return newColor;
}

struct Color specularHighlight(long double E, struct Color color){
	struct Color newColor;

	newColor.r = color.r + (E * (1 - color.r));
	newColor.g = color.g + (E * (1 - color.g));
	newColor.b = color.b + (E * (1 - color.b));

	return newColor;
}
//////////////END Ilumination models



//////////////Ray Tracer Stuff

//DONE
//Interseccion con esferas
struct Intersection *sphereIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	long double t, t1, t2;

	long double Xdif = anchor.x - object.Xc;
	long double Ydif = anchor.y - object.Yc;
	long double Zdif = anchor.z - object.Zc;

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

//DONE
//Sacar la normal normalizada de una esfera en un punto dado por un vector
struct Vector sphereNormal(struct Object object, struct Vector vector){
	struct Vector normal;
	
	normal.x = vector.x - object.Xc;
	normal.y = vector.y - object.Yc;
	normal.z = vector.z - object.Zc;

	normal = normalize(normal);

	return normal;
}

//DONE
//Funcion que pide la primer interseccion de un rayo dado.
struct Intersection getFirstIntersection(struct Vector anchor, struct Vector direction){
	int k;
	long double tmin;
	struct Intersection intersection;
	struct Intersection * tempIntersection;

	tmin = 100000;
	tempIntersection = NULL;
	intersection.distance = -1;

	int objectsAmount = sizeof(Objects)/sizeof(Objects[0]);
	for(k = 0; k < numberObjects; k++){
		tempIntersection = Objects[k].intersectionFuncion(anchor, direction, Objects[k]);
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
		
		struct Object Q = intersection.object;

		struct Vector L;
		struct Vector intersectVector = {intersection.Xi, intersection.Yi, intersection.Zi};
		struct Vector N = Q.normalVector(Q, intersectVector);
		struct Vector R;


		long double Fatt;
		long double I = 0.0;
		long double E = 0.0;

		for(k = 0; k < numberLights; k++){
			struct Intersection obstacle;
			struct Vector light = {Lights[k].Xp - intersection.Xi, Lights[k].Yp - intersection.Yi, Lights[k].Zp - intersection.Zi};
			L = normalize(light);
			
			long double a = pointProduct(N, L);
			//struct Vector another = {N.x - L.x, N.y - L.y, N.z - L.z};
			//R = crossProduct(R, another);
			R.x = (2 * N.x * a) - L.x;
			R.y = (2 * N.y * a) - L.y;
			R.z = (2 * N.z * a) - L.z;
			R = normalize(R);

			obstacle = getFirstIntersection(intersectVector, L);
			if(obstacle.distance < e){
				long double pp = pointProduct(N, L);
				long double distanceToLight = getNorm(L);
				Fatt = getAtunuationFactor(Lights[k], distanceToLight);
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
/**-------------LECTURA ARCHIVOS-------------------------**/

void createObjectFromData(long double *data, int whichObjectCreate){
	/*whichObjectCreate indica qué objeto crear
		SI está en 1, crea luces
		Si está en 2, crea esferas
		Si está en 3, crea polígonos, 
		SI está en 4, crea cilindros
		SI está en 5, crea conos.
	data es el arreglo de valores del objeto a crear. Se asume que estará completo.*/
	printf("%i \n",whichObjectCreate);
	switch(whichObjectCreate){
		case 1:
			{
				printf("Datos para crear luz %i ...\n",lightIndex);
				printf("Xp %LF \n", data[0]);
				printf("Yp %LF \n", data[1]);
				printf("Zp %LF \n", data[2]);
				printf("c1 %LF \n", data[3]);
				printf("c2 %LF \n", data[4]);
				printf("c3 %LF \n", data[5]);
				printf("Ip %LF \n", data[6]);
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
				printf("Luz insertado");
				return;
				/*Lo hago como si fuera una pila y ustedes no pueden detenerme.*/
				}
	case 2:
			{
				printf("Datos para crear esfera %i ...\n",objectIndex);
				printf("Xc %LF \n", data[0]);
				printf("Yc %LF \n", data[1]);
				printf("Zc %LF \n", data[2]);
				printf("R %LF \n", data[3]);
				printf("Kd %LF \n", data[4]);
				printf("Ka %LF \n", data[5]);
				printf("Kn %LF \n", data[6]);

				printf("Ks %LF \n", data[7]);
				printf("R %LF \n", data[8]);
				printf("G %LF \n", data[9]);
				printf("B %LF \n", data[10]);
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
				esfera .color=colorSphere;
				Objects[objectIndex]=esfera;
				objectIndex++;
				printf("Radio de esfera %i : %LF \n",objectIndex, Objects[objectIndex].other);
				/*Lo hago como si fuera una pila y ustedes no pueden detenerme.*/
			}
		case 3: 
			{}
		case 4:
			{}
		case 5:
			{}
	}
}

//DONE
long double *obtainPointFromString(char stringPoint[]){
	/*Devuelve los tres valores long double de un punto tridimensional 
	a partir de la forma Xp,Yp,Zp. RECORDAR UTILIZAR free(valorDevuelto) tras
	 usarla.*/
	char *pch;
	long double *pointDimensions = malloc(sizeof(long double) * 3);
	int currentDimension=0;
	pch = strtok (stringPoint,",");
	while (pch != NULL)
	{
		sscanf(pch, "%LF", &pointDimensions[currentDimension]);
		pch = strtok (NULL, ",");
		currentDimension++;
	}
	return pointDimensions;
}

//////////////Files
//TO-FIX. NO DEVUELVE VALUES CORRECTAMENTE.
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
			if (*counterValueSegment==0){ //Lee Ilum. AMbiente
				//printf("Lee iluminación ambiente \n");
 				sscanf(lineRead, "%LF", &Ia);
 				*numberValuesRead = 0;
 				return NULL;
			}
		case 1: //Luces
			switch(*counterValueSegment){
				case 0: //Nombre o posición de la luz
					if (strstr(lineRead,"Luz_")==NULL){ //No es el nombre
						long double *positionLight = obtainPointFromString(lineRead);
						values = malloc(sizeof(long double)*3);
						values[0] = positionLight[0];
						values[1] = positionLight[1];
						values[2] = positionLight[2];
						//memcpy(values, positionLight, 3);
						(*counterValueSegment)++;
						*numberValuesRead = 3;
						free (positionLight);
						printf("%LF %LF %LF \n", values[0], values[1],values[2]);
						return values;
					}
					*numberValuesRead = 0;
					return NULL;
				case 1: //Coeficiente atenunacion constante c1
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					(*counterValueSegment)++;
					*numberValuesRead = 1;
					return values;
				case 2: ////Coeficiente atenunacion lineal c2
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					*numberValuesRead = 1;
					(*counterValueSegment)++;
					return values;
				case 3: ////Coeficiente atenunacion cuadrático c3
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					*numberValuesRead = 1;
					(*counterValueSegment)++;
					return values;
				case 4: //Intensidad de la luz
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					*numberValuesRead = 1;
					(*counterValueSegment)=0;
					return values;
			}
		case 2: //ESferas
			switch(*counterValueSegment){
				case 0: //Nombre o posición de la esfera
					if (strstr(lineRead,"Esfera_")==NULL){ //No es el nombre
						long double *positionSphere = obtainPointFromString(lineRead);
						values = malloc(sizeof(long double)*3);
						values[0] = positionSphere[0];
						values[1] = positionSphere[1];
						values[2] = positionSphere[2];
						free (positionSphere);
						*numberValuesRead = 3;
						(*counterValueSegment)++;
						return values;
					}
					*numberValuesRead = 0;
					return NULL;
				case 1: //Radio de la esfera
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					(*counterValueSegment)++;
					*numberValuesRead = 1;
					return values;
				case 2: ////Kd de la esfera o coeficiente de reflexion difusa
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					(*counterValueSegment)++;
					*numberValuesRead = 1;
					return values;
				case 3: ////Coeficiente de iluminacion ambiente Ka
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					(*counterValueSegment)++;
					*numberValuesRead = 1;
					return values;
				case 4: //Factor de atenuación de reflexión especular Kn
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					(*counterValueSegment)++;
					*numberValuesRead = 1;
					return values;
				case 5: //Ks
					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					(*counterValueSegment)++;
					*numberValuesRead = 1;
					return values;
				case 6: //Tripleta RGB
					{
					values = malloc(sizeof(long double)*3);	
					long double *rgbColors = obtainPointFromString(lineRead);
					values[0] = rgbColors[0];
					values[1] = rgbColors[1];
					values[2] = rgbColors[2];
					free (rgbColors);
					*numberValuesRead = 3;
					(*counterValueSegment)=0;
					//printf("counterValueSegment tras asignar RGB %i \n",(*counterValueSegment));
					return values;
					}
			}

		/*case 3: //Poligonos
		case 4: //Cilindros
		case 5: //Conos
		*/
	}
	return values;
}








//Leer archivos con la escena
//TO-FIX
void getSceneObjects(){

	int i, j, c;
	int state = 0;
	/*State indica qué segmento lee de la escena
		Si está en 0, lee variables de escena
		SI está en 1, lee luces
		Si está en 2, lee esferas
		Si está en 3, lee polígonos, 
		SI está en 4, lee cilindros
		SI está en 5, lee conos.*/
	int counterValueSegment = 0; 
	//Contador que indica qué valor del objeto está leyendo del actual objeto de
	//un segmento. 
	char temporalBuffer[100]; //Aquí se guardará lo leído cada línea
	long double *valuesRead; //Se guarda los valores del objeto para finalmente crearlo.
	int indexValuesRead = 0; //Pos de valuesRead 
	int currentTypeObjectReading = 1;
	/*SI está en 1, crea luces
		Si está en 2, crea esferas
		Si está en 3, crea polígonos, 
		SI está en 4, crea cilindros
		SI está en 5, crea conos.
	*/
	FILE* file; //archivo
	if (file = fopen(escenaFile, "r")){

		while (fgets(temporalBuffer, 100, file)!=NULL){//Mientras el archivo siga teniendo algo

			if (strstr(temporalBuffer, "Escena")!=NULL){
			//Entra en state = 0. Escena
				state = 0;
				counterValueSegment = 0;
				indexValuesRead=0;
				valuesRead=NULL;
				printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer,"Luces")!=NULL){
			//Entra en state = 1. Luces
				state = 1;
				counterValueSegment = 0;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead = malloc(sizeof(long double)*7); //Las luces siempre seran 7 valores
				currentTypeObjectReading=1;
				printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Esferas ")!=NULL){
			//Entra en state = 2. Esferas
				state = 2;
				indexValuesRead=0;
				free(valuesRead);
				valuesRead=malloc(sizeof(long double)*11); //Esferas siempre serán 11 valores
				counterValueSegment = 0;
				currentTypeObjectReading=2;
				printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Poligonos")!=NULL){
			//Entra en state = 3. Poligonos
				state = 3;
				counterValueSegment = 0;
				indexValuesRead=0;
				valuesRead=NULL;
				currentTypeObjectReading=3;
				printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Cilindros")!=NULL){
			//Entra en state = 4. Cilindros
				state = 4;
				counterValueSegment = 0;
				indexValuesRead=0;
				valuesRead=NULL;
				currentTypeObjectReading=4;
				printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Conos")!=NULL){
			//Entra en state = 5. Conos
				state = 5;
				counterValueSegment = 0;
				indexValuesRead=0;
				valuesRead=NULL;
				currentTypeObjectReading=5;
				printf("%s",temporalBuffer);
				continue;
			}
			
			int numberValuesRead = 0;
			printf("%s \n",temporalBuffer);
			printf("state %i \n" ,state);
			printf("counterValueSegment %i \n", counterValueSegment);
			long double *valuesReadTemp = readValueFromLine(state, &counterValueSegment, temporalBuffer, &numberValuesRead);
			//printf("numberValuesRead %i \n",numberValuesRead);
			if (valuesReadTemp == NULL){ //Se devolvió NULL
				//printf("readValueFromLine devolvió NULL \n");
				continue;
			}
			
			if (state == 0){ //Se leyó iluminación ambiente y no se debe extraer nada.
				free(valuesReadTemp);
				continue;
			}
			int i = 0;
			//TO FIX- SEG-FAULT AL LEER DE valuesReadTemp.
			
			for (i = 0; i < numberValuesRead; i++){
				
				valuesRead[indexValuesRead+i] = valuesReadTemp[i];
				printf("Valor devuelto por readValueFromLine: %LF \n", valuesReadTemp[i]);

			}
			indexValuesRead+=numberValuesRead;
			free(valuesReadTemp);
			
			/*SI counterValueSegment vuelve como un 0, quiere decir que ya
			 se leyeron los datos de dicho objeto/luz, ergo, se crea.*/
		if (counterValueSegment == 0){
				
				createObjectFromData(valuesRead, currentTypeObjectReading);
				indexValuesRead=0;
			}
		}
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

			if (strstr(temporalBuffer, "Luz_")!=NULL){
				numberLights++;
				printf("%s", temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer,"Esfera_")!=NULL){
				numberObjects++;
				printf("%s", temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Poligono_")!=NULL){
				numberObjects++;
				printf("%s", temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Cono_")!=NULL){
				numberObjects++;
				printf("%s", temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Cilindro_")!=NULL){
				numberObjects++;
				printf("%s", temporalBuffer);
				continue;
			}

		}
	}
	fclose(file);
	Objects = malloc(sizeof(struct Object)*numberObjects);
	Lights= malloc(sizeof(struct Light)*numberLights);
}

/**------------------FIN LECTURA ARCHIVOS-------------------------**/



//DONE
int main(int argc, char *arcgv[]){
	//Esfera ad hoc


	howManyObjectsLights();
	printf("Luces %i \n", numberLights);
	printf("Objetos %i \n", numberObjects);
	getSceneObjects();

	printf("Ip de Luz 0 %LF \n", Objects[0].other);
	printf("RAdio de esfera 0 %LF \n", Objects[0].other);
	printf("Color R de esfera 0 %LF \n", Objects[0].color.r);
	/*struct Object sphere1;
	sphere1.Xc = 200;
	sphere1.Yc = 200;
	sphere1.Zc = 650;
	sphere1.other = 200;

	sphere1.color.r = 1;
	sphere1.color.g = 0;	
	sphere1.color.b = 0;
	
	sphere1.Kd = 0.8;
	sphere1.Ka = 0.4;
	sphere1.Kn = 0.8;
	sphere1.Ks = 0.3;
	sphere1.normalVector = sphereNormal;
	sphere1.intersectionFuncion = sphereIntersection;

	Objects[0] = sphere1;

	sphere1.Xc = 100;
	sphere1.Yc = 200;
	sphere1.Zc = 100;
	sphere1.other = 60;

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

	*/
	int i, j;

	long double L;
	long double Xw, Yw;
	long double Xd, Yd, Zd;

	struct Color color;
	
	struct Vector direction;
	
	long double Xdif = Xmax - Xmin;
	long double Ydif = Ymax - Ymin;

	
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

}
//////////////END Ray Tracer Stuff