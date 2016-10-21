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
	struct Vector *points;
	int numberVertexes;
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
/*Esta global llevara el contador de cuántos puntos/vertices contiene el poligono
actual del cual se está leyendo la información.*/

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



//////////////Ilumination

//DONE
//Calcula el factor de atenuacion de una luz
long double getAtunuationFactor(struct Light light, long double distance){
	long double value = 1 / (light.c1 + (light.c2 * distance) + (light.c3 * pow(distance, 2)));
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

//DONE
//Regresa el color difuso co reflexion especular
struct Color specularHighlight(long double E, struct Color color){
	struct Color newColor;

	newColor.r = color.r + (E * (1 - color.r));
	newColor.g = color.g + (E * (1 - color.g));
	newColor.b = color.b + (E * (1 - color.b));

	return newColor;
}
//////////////END Ilumination



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

	return normal;
}

//CHECK
struct Vector polygonNormal(struct Object object, struct Vector vector){
	struct Vector point0 = object.points[0];
	struct Vector point1 = object.points[1];
	struct Vector point2 = object.points[2];

	struct Vector vector1 = {point1.x - point0.x, point1.y - point0.y, point1.z - point0.z};
	struct Vector vector2 = {point2.x - point1.x, point2.y - point1.y, point2.z - point1.z};
	struct Vector normal = crossProduct(vector1, vector2);
	return normal;
}

//CHECK
long double whatsTheD(struct Object object){
	long double theD = 0;
	struct Vector point = object.points[0];

	theD -= object.Xc * point.x;
	theD -= object.Yc * point.y;
	theD -= object.Zc * point.z;

	return theD;
}

//PENDING
struct Intersection *polygonIntersection(struct Vector anchor, struct Vector direction, struct Object object){
	struct Vector normal = polygonNormal(object, anchor);

	object.Xc = normal.x;
	object.Yc = normal.y;
	object.Zc = normal.z;
	object.other = whatsTheD(object);

	long double L = getNorm(normal);
	object.Xc /= L;
	object.Yc /= L;
	object.Zc /= L;
	object.other /= L;

	long double numerator = -((anchor.x * object.Xc) + (anchor.y * object.Yc) + (anchor.z * object.Zc));
	long double denominator = (direction.x * object.Xc) + (direction.y * object.Yc) + (direction.z * object.Zc);

	if(denominator = 0){
		return NULL;
	}else{
		long double t = numerator / denominator;
		tempIntersect.distance = t;
		tempIntersect.object = object;
		tempIntersect.Xi = anchor.x + (t * direction.x);
		tempIntersect.Yi = anchor.y + (t * direction.y);
		tempIntersect.Zi = anchor.z + (t * direction.z);

		//PENDING: hacer lo de revisar con intersecciones 2D

		return NULL;
	}
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

//DONE
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
		struct Vector N = normalize(Q.normalVector(Q, intersectVector));
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
				long double distanceToLight = getNorm(light);
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
void createObjectFromData(long double *data, int whichObjectCreate, int quantityData){
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
				printf("Radio de esfera %i : %LF \n",objectIndex, Objects[objectIndex].other);
				objectIndex++;
				
				/*Lo hago como si fuera una pila y ustedes no pueden detenerme.*/
				return;
			}
		case 3: 
			{
				printf("Insertando polígono...");
				struct Object polygon;
				struct Color colorPolygon;
				printf("Color g: %LF \n", data[1]);
				colorPolygon.r = data[0];
				colorPolygon.g = data[1];
				colorPolygon.b = data[2];
				polygon.color =  colorPolygon;
				polygon.Kd = data[3];
				polygon.Ka = data[4];
				polygon.Kn = data[5];
				polygon.Ks = data[6];
				//7 Elementos insertados por el momento
				int numVertexesPolygon = (quantityData-7) / 3;
				polygon.numberVertexes = numVertexesPolygon;
				polygon.points = malloc(sizeof(struct Vector)*numVertexesPolygon);
				int vertexPolygonIndex = 0;
				for (int i =0; i+7 < quantityData;){
					struct Vector vertex;
					vertex.x = data[7+i];
					i++;
					vertex.y = data[7+i];
					i++;
					vertex.z = data[7+i];
					i++;
					polygon.points[vertexPolygonIndex];
					vertexPolygonIndex++;
				}
				Objects[objectIndex] = polygon;
				printf("RGB de poligono 1: %LF %LF %LF \n", Objects[objectIndex].color.r, Objects[objectIndex].color.g, Objects[objectIndex].color.b);
				objectIndex++;

			}
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
				if ((*counterValueSegment)== 0){
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

				}else if((*counterValueSegment)== 1 || //c1
					(*counterValueSegment) == 2 || //c2
					(*counterValueSegment) == 3 || //c3
					(*counterValueSegment) == 4){ //Ip

					values = malloc(sizeof(long double));
					sscanf(lineRead, "%LF", &values[0]);
					if((*counterValueSegment) == 4){//ES Ip, último valor
						(*counterValueSegment) = 0;
					}else{
						(*counterValueSegment)++;
					}
					*numberValuesRead = 1;
					return values;
				}
		case 2: //ESferas
			if ((*counterValueSegment)== 0){ //Nombre o posición de la esfera
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
			}else if((*counterValueSegment)== 1 || //Radio de la esfera
					(*counterValueSegment) == 2 || //Kd de la esfera o coeficiente de reflexion difusa
					(*counterValueSegment) == 3 || //Coeficiente de iluminacion ambiente Ka
					(*counterValueSegment) == 4 || //Factor de atenuación de reflexión especular Kn
					(*counterValueSegment) == 5 ){ //Ks
				values = malloc(sizeof(long double));
				sscanf(lineRead, "%LF", &values[0]);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
			}else if((*counterValueSegment)== 6){ //RGB
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
		case 3: //Poligonos
			if ((*counterValueSegment)== 0){ //Nombre o RGB del poligono
				if (strstr(lineRead,"Poligono")==NULL){ //No es el nombre
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
				}
				*numberValuesRead = 0;
				return NULL;
			}else if((*counterValueSegment)== 1 || //Kd
					(*counterValueSegment) == 2 || //Ka
					(*counterValueSegment) == 3 || //Kn
					(*counterValueSegment) == 4){ //Ks
				values = malloc(sizeof(long double));
				sscanf(lineRead, "%LF", &values[0]);
				(*counterValueSegment)++;
				*numberValuesRead = 1;
				return values;
			}else if((*counterValueSegment) == 5){ //Poligono Vertice
				if (strstr(lineRead,"END_Poligono")!=NULL){
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
		/*case 4: //Cilindros
		case 5: //Conos*/
		
	}
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
				free(valuesRead);
				valuesRead=malloc(sizeof(long double)*37); //Polígonos max size 
				currentTypeObjectReading=3;
				printf("%s",temporalBuffer);
				continue;
			}else if (strstr(temporalBuffer, "Cilindros")!=NULL){
			//Entra en state = 4. Cilindros
				state = 4;
				counterValueSegment = 0;
				indexValuesRead=0;
				free(valuesRead);
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
			//printf("%s \n",temporalBuffer);
			//printf("state %i \n" ,state);
			//printf("counterValueSegment %i \n", counterValueSegment);
			long double *valuesReadTemp = readValueFromLine(state, &counterValueSegment, temporalBuffer, &numberValuesRead);
			
			if (valuesReadTemp == NULL){ //Se devolvió NULL
				//printf("readValueFromLine devolvió NULL \n");
				continue;
			}

			if (state == 0){ //Se leyó iluminación ambiente y no se debe extraer nada.
				free(valuesReadTemp);
				continue;
			}
			int i = 0;
			for (i = 0; i < numberValuesRead; i++){
				
				valuesRead[indexValuesRead+i] = valuesReadTemp[i];
				//printf("Valor devuelto por readValueFromLine: %LF \n", valuesReadTemp[i]);
				indexValuesRead+=1;

			}
			/*if(numberValuesRead>0)
			{
				free(valuesReadTemp);
			}*/
			/*SI counterValueSegment vuelve como un 0, quiere decir que ya
			 se leyeron los datos de dicho objeto/luz, ergo, se crea.*/
			if (counterValueSegment == 0){
				printf("Creando objeto...\n");
				createObjectFromData(valuesRead, currentTypeObjectReading, indexValuesRead);
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


void freeMemoryOfVertexes(){
	/*Libera la memoria allocated en los arreglos dinámicos de vértices perteneciente a los polígonos*/
	for (int i = 0; i < numberObjects; i++){
		if (Objects[i].numberVertexes > 0){
			//Se asegura que sea un polígono
			free (Objects[i].points);
		}
	}
}

//DONE
int main(int argc, char *arcgv[]){
	howManyObjectsLights();
	printf("Luces %i \n", numberLights);
	printf("Objetos %i \n", numberObjects);
	getSceneObjects();
	/*
	printf("")
	printf("Ip de Luz 0 %LF \n", Objects[0].other);
	printf("RAdio de esfera 0 %LF \n", Objects[0].other);
	printf("Color R de esfera 0 %LF \n", Objects[0].color.r);
	int i, j;

	long double L;
	long double Xw, Yw;
	long double Xd, Yd, Zd;

	struct Color color;
	
	struct Vector direction;
	
	long double Xdif = Xmax - Xmin;
	long double Ydif = Ymax - Ymin;
*//*
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
	saveFile();*/
	//freeMemoryOfVertexes();
	free(Objects);
	free(Lights);
}
//////////////END Ray Tracer Stuff