#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structures.h"
#include "TDA_LDE.h"


long double Ia = 0.2;  //Iluminacion ambiente por defecto

//Image resolution
 int Hres = 1366;
 int Vres = 768;

//Eye
 long double Xe = 0;
 long double Ye = 0;
 long double Ze = -200;

//Window
 long double Xmax = 100;
 long double Ymax = 100;
 long double Xmin = -100;
 long double Ymin = -100;

//Others
 long double e = 0.005; //Epsilon 



 struct Color BACKGROUND = {0.3, 0.3, 0.3}; //Gris
 struct Color Framebuffer[768][1366];


 int numberLights = 0;
 int numberObjects = 0;
 struct Color BACKGROUND; //Gris
 struct Color Framebuffer[768][1366];


