#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "TDA_LDE.h"
#include "structures.h"
#include "vectors.h"
#include "colorLighting.h"
#include "malloc.h"
#include "fileReader.h"


//Esta variable sólo VIVE en este archivo.
static char* escenaFile = "escena1.txt";



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
	long double *values = malloc(sizeof(long double)*20);
	//RECORDAR LIBERARLA TRAS USO.
	switch(state){
		case 0: //Escena
			if (counterValueSegment==0){ //Lee Ilum. AMbiente
 				sscanf(lineRead, "%LF", &Ia);
 				*numberValuesRead = 0;
 				return NULL;
			}
		case 1: //Luces
			switch(*counterValueSegment){
				case 0: //Nombre o posición de la luz
					if (strstr(lineRead,"Luz")==NULL){ //No es el nombre
						long double *positionLight = obtainPointFromString(lineRead);
						values[0] = positionLight[0];
						values[2] = positionLight[1];
						values[1] = positionLight[2];
						memcpy(values, positionLight, 3);
						*numberValuesRead = 3;
						free (positionLight);
						printf("%LF %LF %LF \n", values[0], values[1],values[2]);
						return values;
					}
					*numberValuesRead = 0;
					return NULL;
				case 1: //Coeficiente atenunacion constante c1
					sscanf(lineRead, "%LF", &values[0]);
					*counterValueSegment++;
					*numberValuesRead = 1;
					return values;
				case 2: ////Coeficiente atenunacion lineal c2
					sscanf(lineRead, "%LF", &values[0]);
					*numberValuesRead = 1;
					*counterValueSegment++;
					return values;
				case 3: ////Coeficiente atenunacion cuadrático c3
					sscanf(lineRead, "%LF", &values[0]);
					*numberValuesRead = 1;
					*counterValueSegment++;
					return values;
				case 4: //Intensidad de la luz
					sscanf(lineRead, "%LF", &values[0]);
					*numberValuesRead = 1;
					*counterValueSegment=0;
					return values;
			}
		case 2: //ESferas
			switch(*counterValueSegment){
				case 0: //Nombre o posición de la esfera
					if (strstr(lineRead,"Esfera")==NULL){ //No es el nombre
						long double *positionSphere = obtainPointFromString(lineRead);
						values[0] = positionSphere[0];
						values[2] = positionSphere[1];
						values[1] = positionSphere[2];
						free (positionSphere);
						*numberValuesRead = 3;
						*counterValueSegment++;
						return values;
					}
					*numberValuesRead = 0;
					return NULL;
				case 1: //Radio de la esfera
					sscanf(lineRead, "%LF", &values[0]);
					*counterValueSegment++;
					*numberValuesRead = 1;
					return values;
				case 2: ////Kd de la esfera o coeficiente de reflexion difusa
					sscanf(lineRead, "%LF", &values[0]);
					*counterValueSegment++;
					*numberValuesRead = 1;
					return values;
				case 3: ////Coeficiente de iluminacion ambiente Ka
					sscanf(lineRead, "%LF", &values[0]);
					*counterValueSegment++;
					*numberValuesRead = 1;
					return values;
				case 4: //Factor de atenuación de reflexión especular Kn
					sscanf(lineRead, "%LF", &values[0]);
					*counterValueSegment++;
					*numberValuesRead = 1;
					return values;
				case 5: //Tripleta RGB
					{
					long double *rgbColors = obtainPointFromString(lineRead);
					values[0] = rgbColors[0];
					values[2] = rgbColors[1];
					values[1] = rgbColors[2];
					free (rgbColors);
					*numberValuesRead = 3;
					*counterValueSegment=0;
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


void createObjectFromData(long double *data, int whichObjectCreate){
	/*whichObjectCreate indica qué objeto crear
		SI está en 1, crea luces
		Si está en 2, crea esferas
		Si está en 3, crea polígonos, 
		SI está en 4, crea cilindros
		SI está en 5, crea conos.
	data es el arreglo de valores del objeto a crear. Se asume que estará completo.*/
	switch(whichObjectCreate){
		case 1:
			{
				printf("hola ");
				printf("%LF \n", data[0]);
				struct Light luz;
				luz.Xp=data[0];
				printf("%LF \n", data[1]);
				luz.Yp=data[1];
				luz.Zp=data[2];
				luz.c1=data[3];
				luz.c2=data[4];
				luz.c3=data[5];
				luz.Ip=data[6];

				struct Light *pointerLight = &luz;
				printf("%LF \n", pointerLight->Ip);
				insertPrimeroLuz(&Lights, pointerLight);
				printf("Elemento insertado");
				return;
				/*Lo hago como si fuera una pila y ustedes no pueden detenerme.*/
			}
		case 2:
			{
				struct Object esfera;
				esfera.Xc=data[0];
				esfera.Yc=data[1];
				esfera.Zc=data[2];
				esfera.otherData[0]=data[3]; //Radio
				esfera.Kd=data[4]; 
				esfera.Ka=data[5];
				esfera.Kn=data[6];

				struct Color colorSphere;
				colorSphere.r = data[7];
				colorSphere.g = data[8];
				colorSphere.b = data[9];
				esfera .color=colorSphere;

				struct Object *pointerObject = &esfera;
				insertPrimero(&Objects, pointerObject); 
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
				valuesRead=NULL;
				currentTypeObjectReading=1;
				printf("Luces");
				continue;
			}else if (strstr(temporalBuffer, "Esferas ")!=NULL){
			//Entra en state = 2. Esferas
				state = 2;
				indexValuesRead=0;
				valuesRead=NULL;
				counterValueSegment = 0;
				currentTypeObjectReading=2;
				printf("Esferas");
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
			long double *valuesReadTemp = readValueFromLine(state, &counterValueSegment, temporalBuffer, &numberValuesRead);
			printf("numberValuesRead %i \n",numberValuesRead);
			if (numberValuesRead==0){ //Se devolvió NULL
				continue;
			}
			
			if (state == 0){ //Se leyó iluminación ambiente y no se debe extraer nada.
				free(valuesReadTemp);
				continue;
			}
			int i = 0;
			printf("i %i \n", i);
			//TO FIX- SEG-FAULT AL LEER DE valuesReadTemp.
			printf("Valor devuelto por readValueFromLine: %lf", valuesReadTemp[0]);
			for (i = 0; i < numberValuesRead; i++){
				
				valuesRead[indexValuesRead+i] = valuesReadTemp[i];

			}
			/*SI counterValueSegment vuelve como un 0, quiere decir que ya
			 se leyeron los datos de dicho objeto/luz, ergo, se crea.*/
			if (counterValueSegment == 0){
				createObjectFromData(valuesReadTemp, currentTypeObjectReading);
			}
			free(valuesReadTemp);
		}
	}
	fclose(file);	
	//Pendiente
}