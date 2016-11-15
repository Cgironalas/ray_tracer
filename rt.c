#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "malloc.h"
#include <time.h>

// Constants and data types: =====================================
	//Image resolution
	static int Hres;
	static int Vres;

	//Advanced limits
	static int maxAA;
	static long double similar = 0;
	static int maxReflection;
	static int maxTransparency;
	static long double rays = 0;

	//Window
	static long double Xmax;
	static long double Ymax;
	static long double Xmin;
	static long double Ymin;
	static long double Xdif;
	static long double Ydif;

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
	
		long double other;  //Radius of the sphere, cylinder, discs and elipses (k) Also used for the K in quadratics
	  	struct Vector directionVector; //for the cylinder, cone, discs and elipses
		long double extraD; //used in disc and elipses for calculating the d of the normal

	  	long double Xother; //For the elipse
		long double Yother;
	  	long double Zother;

	  	long double Kd;
	  	long double Ka;
	  	long double Kn;
	  	long double Ks;
	  	long double o1;
	  	long double o2; 
	  	long double o3;

	  	long double A; // quadratics
	  	long double B;
	  	long double C;
	  	long double D;
	  	long double E;
	  	long double F;
	  	long double G;
	  	long double H;
	  	long double I;
	  	long double J;
		  
	  	int pointAmount;
		  
	  	long double D1;  
	  	long double D2; 	  //Used to cut the cilinder and the cone
	  	long double K1;
	  	long double K2; 	  //K1 and K2 are onlt for the cone.
	  	long double height;
	  	
	  	struct Color color;
	  	struct Point2D *points2D;
	  	struct Point3D *points3D;
	  	struct Vector (*normalVector)();
	  	struct Intersection (*intersectionFuncion)();
	  	struct Color (*retrieveTextureColor)();
	  	
	  	struct Vector x0y0z0; //Usadas para evaluar el plano de texturas en objetos referentes en un pano
	  	struct Vector x1y1z1;
	  	struct Vector x2y2z2;
	  	struct Vector x3y3z3;
	  	
	  	int numberPlaneCuts;
	  	int numberTextures;
	  	int numberDraftPlanes;
	  	
	  	struct PlaneCut *planeCuts;
	  	struct Texture *textures;
	  	struct DraftPlane *draftPlanes;
	};

	struct PlaneCut{
	  	struct Vector normal;
	  	struct Vector point;
	  	long double d;
	};

	struct Texture{
	  	char* filename;
	  	struct Color **textureMap;
	  	int vRes;
  		int hRes;
	  	struct Vector greenwich; //No asignado en polígonos 
	  	struct Vector north;
	};

	struct DraftPlane{ //Plano de calado
	  	char* Filename;
	  	int ***textureMap;
	  	int vRes;
 		int hRes;
	  	struct Vector greenwich; //No asignado en polígonos 
	  	struct Vector north;
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
	static struct Vector V;
	static struct Vector eye;
	static struct Color **Framebuffer;
	static struct Color background;
	static int numberObjects = 0;
	static int numberLights = 0;
	static int lightIndex = 0;
	static int objectIndex = 0;
	static char* escenaFile = "escena1.txt";

	long double min(long double a, long double b){
	    if(a < b) { return a; }
	    else { return b; }
	}

	long double max(long double a, long double b){
		if(a > b) { return a; }
		else { return b; }
	}
// ===============================================================

// Cut planes ====================================================
	int testPlaneCut(struct PlaneCut plane, long double x, long double y, long double z){
		int val = (plane.normal.x * x) + (plane.normal.y * y) + (plane.normal.z * z) + plane.d;
		if(val > 0){
			return 1;
		}else{
			return 0;
		}
	}

	int testIntersection(long double x, long double y, long double z, struct Object object){
		int k, sign;
		int accept = 1;
		int amount = object.numberPlaneCuts;
		for(k = 0; k < amount; k++){
			sign = testPlaneCut(object.planeCuts[k], x, y, z);
			if(sign == 1){
				accept = 0;
			}
		}

		return accept;
	}
// ===============================================================

// Vectors =======================================================
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
// ===============================================================

// Save ppm ======================================================
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
				int R = (int) 255 * Framebuffer[i][j].r; 
				int G = (int) 255 * Framebuffer[i][j].g;
				int B = (int) 255 * Framebuffer[i][j].b;
				fprintf(file, "%i %i %i   ", R, G, B);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
// ===============================================================

// Lights ========================================================
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
// ===============================================================

// RAY TRACER ====================================================

// Spheres: ======================================================
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
					tempIntersect.distance = t;
					tempIntersect.object = object;
					tempIntersect.Xi = anchor.x + (t * direction.x);
					tempIntersect.Yi = anchor.y + (t * direction.y);
					tempIntersect.Zi = anchor.z + (t * direction.z);

					int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
					if(accept == 0){
						t = max(t1, t2);
						tempIntersect.distance = t;
						tempIntersect.object = object;
						tempIntersect.Xi = anchor.x + (t * direction.x);
						tempIntersect.Yi = anchor.y + (t * direction.y);
						tempIntersect.Zi = anchor.z + (t * direction.z);

						int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
						if(accept == 0){
							tempIntersect.null = 1;
						}
					}
				}else{
					t = t1;
					tempIntersect.distance = t;
					tempIntersect.object = object;
					tempIntersect.Xi = anchor.x + (t * direction.x);
					tempIntersect.Yi = anchor.y + (t * direction.y);
					tempIntersect.Zi = anchor.z + (t * direction.z);

					int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
					if(accept == 0){
						tempIntersect.null = 1;
					}
				}
			}else{
				if(t2 > e){
					t = t2;
					tempIntersect.distance = t;
					tempIntersect.object = object;
					tempIntersect.Xi = anchor.x + (t * direction.x);
					tempIntersect.Yi = anchor.y + (t * direction.y);
					tempIntersect.Zi = anchor.z + (t * direction.z);

					int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
					if(accept == 0){
						tempIntersect.null = 1;
					}
				}else{
					tempIntersect.null = 1;
				}
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

// Polygons: =====================================================
	int getSign(long double v){
		if(v >= 0){ return 1; }
		else{ return 0; }
	}

	struct Intersection polygonIntersection(struct Vector anchor, struct Vector direction, struct Object object){
		
		long double denominator = (direction.x * object.Xc) + (direction.y * object.Yc) + (direction.z * object.Zc);

		struct Intersection tempIntersect;
		tempIntersect.null = 0;

		if(denominator == 0){
			tempIntersect.null = 1;
			return tempIntersect;

		}else{
			long double numerator = -(anchor.x*object.Xc + anchor.y*object.Yc + anchor.z*object.Zc + object.other);
		
			long double t = numerator / denominator;
			tempIntersect.distance = t;
			tempIntersect.object = object;
			tempIntersect.Xi = anchor.x + (t * direction.x);
			tempIntersect.Yi = anchor.y + (t * direction.y);
			tempIntersect.Zi = anchor.z + (t * direction.z);

			long double maxA_B = max(fabs(object.Xc), fabs(object.Yc));
			long double maxA_B_C = max(maxA_B, fabs(object.Zc)); 
			long double u, v;
			
			if(maxA_B_C == fabs(object.Xc)){
				u = tempIntersect.Zi;
				v = tempIntersect.Yi;
			}else if(maxA_B_C == fabs(object.Yc)){
				u = tempIntersect.Xi;
				v = tempIntersect.Zi;
			}else if(maxA_B_C == fabs(object.Zc)){
				u = tempIntersect.Xi;
				v = tempIntersect.Yi;
			}

			int NC = 0;
			int NV = object.pointAmount;

			struct Point2D *points2DArrayTemp = malloc(sizeof(struct Point2D)*NV);
			for(int i = 0; i < NV; i++){ //Traslado a origen

				points2DArrayTemp[i].u =object.points2D[i].u - u;
				points2DArrayTemp[i].v =object.points2D[i].v - v;
			}

			int SH = getSign(points2DArrayTemp[0].v);
			int NSH;
			int a = 0;
			int b = (a+1)%NV;
			for (a = 0; a < NV-1;){
				NSH = getSign(points2DArrayTemp[b].v);
				
				if(SH != NSH){
					if(points2DArrayTemp[a].u > 0 && points2DArrayTemp[b].u > 0){
						NC++;
					}else if(points2DArrayTemp[a].u > 0 || points2DArrayTemp[b].u > 0){
						long double N = (points2DArrayTemp[b].u - points2DArrayTemp[a].u);
						long double D = (points2DArrayTemp[b].v - points2DArrayTemp[a].v);
						if(D != 0){
							if(points2DArrayTemp[a].u - ((points2DArrayTemp[a].v * N)/D) > 0){
								NC++;
							}	
						}
					}
				}

				SH = NSH;
				a++;
				b++;
			}

			if (NC%2 == 0){ tempIntersect.null = 1; }
			else{ tempIntersect.null = 0; }

			int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
			if(accept == 0){
				tempIntersect.null = 1;
			}

			free(points2DArrayTemp);
			return tempIntersect;	
		}
	}

	struct Vector polygonNormal(struct Object object){
		
		struct Point3D point0 = object.points3D[0];
		struct Point3D point1 = object.points3D[1];
		struct Point3D point2 = object.points3D[2];

		struct Vector vector1 = {point1.x - point0.x, point1.y - point0.y, point1.z - point0.z};
		struct Vector vector2 = {point2.x - point1.x, point2.y - point1.y, point2.z - point1.z};

		struct Vector normal = crossProduct(vector1,vector2);
		//normal = normalize(normal);

		return normal;
	}
// ===============================================================

// Cylinders: ====================================================
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

		//Cálculo
		long double xdxq = xd*xq;
		long double ydyq = yd*yq;
		long double zdzq = zd*zq;

		long double xexq = xe*xq;
		long double yeyq = ye*yq;
		long double zezq = ze*zq;

		long double xoxq = xo*xq;
		long double yoyq = yo*yq;
		long double zozq = zo*zq;

		long double coef1 = xdxq*xq + ydyq*xq + zdzq*xq - xd;
		long double coef2 = xdxq*yq + ydyq*yq + zdzq*yq - yd;
		long double coef3 = xdxq*zq + ydyq*zq + zdzq*zq - zd;

		long double coef4 = xo + xexq*xq - xoxq*xq + yeyq*xq - yoyq*xq + zezq*xq - zozq*xq - xe;
		long double coef5 = yo + xexq*yq - xoxq*yq + yeyq*yq - yoyq*yq + zezq*yq - zozq*yq - ye;
		long double coef6 = zo + xexq*zq - xoxq*zq + yeyq*zq - yoyq*zq + zezq*zq - zozq*zq - ze;

		long double A = pow (coef1,2) + pow (coef2,2) + pow (coef3,2);

		long double B =  2*(coef1*coef4 + coef2*coef5 + coef3*coef6);

		long double C = pow (coef4, 2) + 
						pow (coef5 ,2) + 
						pow (coef6 ,2) - 
						pow(radius,2);

		long double discriminant = pow(B, 2) - (4 * A * C);

		long double t, t1, t2;

		if(discriminant >= 0){
			long double root = sqrt(discriminant);

			B *= -1;

			t1 = (B - root)/(2*A);
			t2 = (B + root)/(2*A);

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

							int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
							if(accept == 0){
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

									int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
									if(accept == 0){
										tempIntersect.null = 1;
									}
									return tempIntersect;
								}else{
									tempIntersect.null = 1;
									return tempIntersect;
								}
							}
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

							
						}else{
							tempIntersect.null = 1;
						}
						int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
						if(accept == 0){
							tempIntersect.null = 1;
						}
						return tempIntersect;
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

							int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
							if(accept == 0){
								tempIntersect.null = 1;
							}
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

							int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
							if(accept == 0){
								tempIntersect.null = 1;
							}
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

		//Cálculo
		long double xxoxq = (x-xo)*xq;
		long double yyoyq = (y-yo)*yq;
		long double zzozq = (z-zo)*zq;
		long double parenth = (xxoxq+yyoyq+zzozq);

		long double pxq = 2*(xo + parenth*xq - x);
		long double pyq = 2*(yo + parenth*yq - y);
		long double pzq = 2*(zo + parenth*zq - z);

		normalCilinder.x = pxq*(pow(xq,2) - 1) +
						   pyq*(yq*xq) +
						   pzq*(zq*xq);

		normalCilinder.y = pxq*(xq*yq) +
						   pyq*(pow(yq,2) - 1) +
						   pzq*(zq*yq);

		normalCilinder.z = pxq*(xq*zq) +
						   pyq*(yq*zq) +
						   pzq*(pow(zq,2) - 1);

		normalCilinder = normalize(normalCilinder);

		return normalCilinder;
	}
// ===============================================================

// Cones: ========================================================
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
		//Cálculo
		long double xdxq = xd*xq;
		long double ydyq = yd*yq;
		long double zdzq = zd*zq;

		long double xexq = xe*xq;
		long double yeyq = ye*yq;
		long double zezq = ze*zq;

		long double xoxq = xo*xq;
		long double yoyq = yo*yq;
		long double zozq = zo*zq;

		long double coef1 = xdxq*xq + ydyq*xq + zdzq*xq - xd;
		long double coef2 = xdxq*yq + ydyq*yq + zdzq*yq - yd;
		long double coef3 = xdxq*zq + ydyq*zq + zdzq*zq - zd;

		long double coef4 = xo + xexq*xq - xoxq*xq + yeyq*xq - yoyq*xq + zezq*xq - zozq*xq - xe;
		long double coef5 = yo + xexq*yq - xoxq*yq + yeyq*yq - yoyq*yq + zezq*yq - zozq*yq - ye;
		long double coef6 = zo + xexq*zq - xoxq*zq + yeyq*zq - yoyq*zq + zezq*zq - zozq*zq - ze;

		//En lugar del radio usa una combinación de estos
		long double coefk = pow(k2/k1,2);
		long double coef7 = xdxq + ydyq + zdzq;
		long double coef8 = xexq - xoxq + yeyq - yoyq + zezq - zozq;

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

						int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
						if(accept == 0){
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

								int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
								if(accept == 0){
									tempIntersect.null = 1;
								}
								return tempIntersect;
							}else{
								tempIntersect.null = 1;
								return tempIntersect;
							}
						}
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

							int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
							if(accept == 0){
								tempIntersect.null = 1;
							}
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

							int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
							if(accept == 0){
								tempIntersect.null = 1;
							}
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

							int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
							if(accept == 0){
								tempIntersect.null = 1;
							}
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


		long double xxoxq = (x-xo)*xq;
		long double yyoyq = (y-yo)*yq;
		long double zzozq = (z-zo)*zq;
		long double parenth = (xxoxq+yyoyq+zzozq);

		long double pxq = 2*(xo + parenth*xq - x);
		long double pyq = 2*(yo + parenth*yq - y);
		long double pzq = 2*(zo + parenth*zq - z);

		long double k1sq = pow(k1,2);
		long double k2sq2 = 2*pow(k2,2);

		long double lastFactorDerivedX = (k2sq2*xq*parenth)/k1sq;
		long double lastFactorDerivedY = (k2sq2*yq*parenth)/k1sq;
		long double lastFactorDerivedZ = (k2sq2*zq*parenth)/k1sq;


		normalCone.x = pxq*(pow(xq,2) - 1) +
					   pyq*(yq*xq) +
					   pzq*(zq*xq) - lastFactorDerivedX;

		normalCone.y = pxq*(xq*yq) +
					   pyq*(pow(yq,2) - 1) +
					   pzq*(zq*yq) - lastFactorDerivedY;

		normalCone.z = pxq*(xq*zq) +
					   pyq*(yq*zq) +
					   pzq*(pow(zq,2) - 1) - lastFactorDerivedZ;

		normalCone = normalize(normalCone);

		return normalCone;
	}
// ===============================================================

// Planes: =======================================================
	long double whatsTheD(struct Object object){
		
		struct Point3D point = object.points3D[0];
		long double theD = -((object.Xc * point.x) + (object.Yc * point.y) + (object.Zc * point.z));
		
		return theD;
	}

	long double whatsTheDGeneral(struct Vector normalNotNormalized, struct Vector point){
		
		long double theD = -((normalNotNormalized.x * point.x) + (normalNotNormalized.y * point.y) + (normalNotNormalized.z * point.z));
		
		return theD;
	}
	struct Object getABCD(struct Object object){
		struct Vector normal = polygonNormal(object);

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
// ===============================================================

// Discs: ========================================================
	struct Intersection discIntersection(struct Vector anchor, struct Vector direction, struct Object object){
		//LA normal normalizada se guardará en directionVector
		long double denominator = (direction.x * object.directionVector.x) + (direction.y * object.directionVector.y) + (direction.z * object.directionVector.z);

		struct Intersection tempIntersect;
		tempIntersect.null = 1;

		if(denominator == 0){
			tempIntersect.null = 1;
			return tempIntersect;

		}else{

			long double numerator = -( (anchor.x*object.directionVector.x) + (anchor.y*object.directionVector.y) + (anchor.z*object.directionVector.z) + object.extraD);
			long double t = numerator / denominator;

			tempIntersect.distance = t;
			tempIntersect.object = object;
			tempIntersect.Xi = anchor.x + (t * direction.x);
			tempIntersect.Yi = anchor.y + (t * direction.y);
			tempIntersect.Zi = anchor.z + (t * direction.z);

			long double distanceToCenter = sqrt(pow(tempIntersect.Xi - object.Xc,2) + 
												pow(tempIntersect.Yi - object.Yc,2) + 
												pow(tempIntersect.Zi - object.Zc,2));

			if (distanceToCenter < object.other) {
				tempIntersect.null = 0;
			}
			
			int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
			if(accept == 0){
				tempIntersect.null = 1;
			}
			return tempIntersect;
		}
	}
	struct Vector discNormal (struct Object object){

		return object.directionVector;
	}
// ===============================================================

// Elipses: ======================================================
	struct Intersection elipseIntersection(struct Vector anchor, struct Vector direction, struct Object object){
		//La normal normalizada se guardará en directionVector
		long double denominator = (direction.x * object.directionVector.x) + (direction.y * object.directionVector.y) + (direction.z * object.directionVector.z);

		struct Intersection tempIntersect;
		tempIntersect.null = 1;

		if(denominator == 0){
			tempIntersect.null = 1;
			return tempIntersect;

		}else{

			long double numerator = -( (anchor.x*object.directionVector.x) + (anchor.y*object.directionVector.y) + (anchor.z*object.directionVector.z) + object.extraD);
			long double t = numerator / denominator;

			tempIntersect.distance = t;
			tempIntersect.object = object;
			tempIntersect.Xi = anchor.x + (t * direction.x);
			tempIntersect.Yi = anchor.y + (t * direction.y);
			tempIntersect.Zi = anchor.z + (t * direction.z);

			long double distanceToD1 = sqrt(pow(tempIntersect.Xi - object.Xc,2) + 
											pow(tempIntersect.Yi - object.Yc,2) + 
											pow(tempIntersect.Zi - object.Zc,2));

			long double distanceToD2 = sqrt(pow(tempIntersect.Xi - object.Xother,2) + 
											pow(tempIntersect.Yi - object.Yother,2) + 
											pow(tempIntersect.Zi - object.Zother,2));

			if ( (distanceToD1 + distanceToD2) < object.other) {
				tempIntersect.null = 0;
			}
			else {
				tempIntersect.null = 1;
			}

			int accept = testIntersection(tempIntersect.Xi, tempIntersect.Yi, tempIntersect.Zi, object);
			if(accept == 0){
				tempIntersect.null = 1;
			}
			return tempIntersect;
		}
	}
	struct Vector elipseNormal (struct Object object){

		return object.directionVector;
	}
// ===============================================================

// Textures ======================================================
	long double uRectangle (struct Vector x0y0z0, struct Vector x1y1z1, struct Vector xiyizi) {
		struct Vector U;
		U.x = x1y1z1.x - x0y0z0.x;
		U.y = x1y1z1.y - x0y0z0.y;
		U.z = x1y1z1.z - x0y0z0.z;
		
		struct Vector i0;
		i0.x = xiyizi.x - x0y0z0.x;
		i0.y = xiyizi.y - x0y0z0.y;
		i0.z = xiyizi.z - x0y0z0.z;

		long double H = getNorm(U);
		U = normalize(U);

		return pointProduct(i0,U)/H; 
	}

	long double vRectangle (struct Vector x0y0z0, struct Vector x3y3z3, struct Vector xiyizi) {
		struct Vector V;
		V.x = x3y3z3.x - x0y0z0.x;
		V.y = x3y3z3.y - x0y0z0.y;
		V.z = x3y3z3.z - x0y0z0.z;
		
		struct Vector i0;
		i0.x = xiyizi.x - x0y0z0.x;
		i0.y = xiyizi.y - x0y0z0.y;
		i0.z = xiyizi.z - x0y0z0.z;

		long double L = getNorm(V);
		V = normalize(V);

		return pointProduct(i0,V)/L; 
	}

	struct Color planeTexture (struct Intersection in) {
		struct Object object = in.object;
		struct Vector ipoint;
		
		ipoint.x = in.Xi;
		ipoint.y = in.Yi;
		ipoint.z = in.Zi;

		long double u = uRectangle(object.x0y0z0,object.x1y1z1,ipoint);
		long double v = vRectangle(object.x0y0z0,object.x3y3z3,ipoint);

		struct Color color;
		for (int i = 0; i < object.numberTextures; i++){

			int xs = object.textures[i].hRes*u;
			int ys = object.textures[i].vRes*v;
			return color = object.textures[i].textureMap[xs][ys];
		}

		return color;
	}
// ===============================================================

// Ray Tracer base ===============================================
	struct Intersection getFirstIntersection(struct Vector anchor, struct Vector direction){
		rays +=1;
		int k;
		int objectsAmount = numberObjects;

		long double tmin;
		struct Intersection intersection;
		struct Intersection tempIntersection;
		intersection.null = 1;

		tmin = 10000000;
		tempIntersection.null = 1;
		
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

	struct Color ponderColor(struct Color baseColor, struct Color reflectionColor, long double o1, long double o2){
		struct Color color;

		color.r = baseColor.r * o1 + reflectionColor.r * o2;
		color.g = baseColor.g * o1 + reflectionColor.g * o2;
		color.b = baseColor.b * o1 + reflectionColor.b * o2;

		return color;
	}

	struct Color getColor(struct Vector anchor, struct Vector direction, int rLevel){
		struct Color color;
		struct Intersection intersection;
		struct Intersection *tempIntersection;
		
		intersection = getFirstIntersection(anchor, direction);

		if(intersection.null == 1){
			color = background;
		}else{
			
			int k;
			int lightsAmount = numberLights;
			
			struct Object Q = intersection.object;
			if (Q.numberTextures != 0) {
				color = Q.retrieveTextureColor(intersection);

			}else {
				color = Q.color;
			}
			
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
				L = light;
				L = normalize(light);
				
				long double pp = pointProduct(N, L);
				obstacle = getFirstIntersection(intersectVector, L);

				long double distanceToLight = getNorm(light);
				
				if(obstacle.null == 1 || (obstacle.distance > e && obstacle.distance > distanceToLight)) {
					
					Fatt = getAttenuationFactor(Lights[k], distanceToLight);
					if(pp > 0.0){

						R.x = (2 * N.x * pp) - L.x;
						R.y = (2 * N.y * pp) - L.y;
						R.z = (2 * N.z * pp) - L.z;
						R = normalize(R);
			
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
			color = difusseColor(I, color);

			
			E = min(1.0, E);
			color = specularHighlight(E, color);

			
			if(rLevel > 0){
				long double pNV = pointProduct(N, V);
				R.x = (2 * N.x * pNV) - V.x;
				R.y = (2 * N.y * pNV) - V.y;
				R.z = (2 * N.z * pNV) - V.z;
				R = normalize(R);
				struct Color reflectionColor = getColor(intersectVector, R, rLevel - 1);
				color = ponderColor(color, reflectionColor, Q.o1, Q.o2);

			struct Color transparencyColor = background;
			int levelsAllowed = maxTransparency;
			while (levelsAllowed > 0 && intersection.object.o3 > 0) {
				transparencyColor = getColor(intersectVector, direction, maxReflection);
				levelsAllowed--;
				if (transparencyColor.r == background.r && transparencyColor.g == background.g && transparencyColor.b == background.b){
					break;
				}
			} 
			color = ponderColor(color, transparencyColor, 1, intersection.object.o3);
		}
		
		return (color);
	}
// ===============================================================

// Lectura de archivos ===========================================
	void printPlaneCuts(struct Object objeto){

	    if (objeto.planeCuts == NULL){
	        printf("\n Este objeto no tiene planos de corte asociados. \n");
	        return;
	    }else{
	        int i = 0;
	        struct Vector normal;
	        struct Vector punto;
	        for (i = 0; i < objeto.numberPlaneCuts; i++){
	            printf("Plano %i: \n", i);
	            normal = objeto.planeCuts[i].normal;
	            punto = objeto.planeCuts[i].point;
	            printf("\t Normal unitaria: %LF, %LF, %LF \n", normal.x, normal.y, normal.z);
	            printf("\t Punto base: %LF, %LF, %LF \n\n", punto.x, punto.y, punto.z);
	        }
	    }
	}

	void printTextures(struct Object objeto, int currentTypeObjectReading){

	    if (objeto.textures == NULL){
	        printf("\n Este objeto no tiene texturas asociados. \n");
	        return;
	    }else{
	        int i = 0;
	        struct Vector norte;
	        struct Vector greenwich;
	        
	        for (i = 0; i < objeto.numberTextures; i++){
	            printf("Textura %i: \n", i);
	            printf("Resolución Textura: %ix%i \n", objeto.textures[i].hRes,objeto.textures[i].vRes);
	            printf("Primer texel de la textura: (%LF, %LF, %LF)", objeto.textures[i].textureMap[0][0].r*255,objeto.textures[i].textureMap[0][0].g*255, objeto.textures[i].textureMap[0][0].b*255);
	            int lastX = objeto.textures[i].hRes-1;
				int lastY = objeto.textures[i].vRes-1;
	            printf("Último texel de la textura: (%LF, %LF, %LF)", objeto.textures[i].textureMap[lastX][lastY].r*255,objeto.textures[i].textureMap[lastX][lastY].g*255, objeto.textures[i].textureMap[lastX][lastY].b*255);

	            for (int f = 0; f<objeto.textures[i].hRes; f++ ){
	            	for (int j = 0; j<objeto.textures[i].vRes; j++ ){
	            		printf("Texel [%i][%i]de la textura: (%LF, %LF, %LF) \n", f, j, objeto.textures[i].textureMap[f][j].r*255,objeto.textures[i].textureMap[f][j].g*255, objeto.textures[i].textureMap[f][j].b*255);

	            	}
	            }
	            //printf("\t Localización : %s \n", objeto.textures[i].filename);
	            if (currentTypeObjectReading == 2 ||currentTypeObjectReading == 4 ||currentTypeObjectReading == 5 ||currentTypeObjectReading == 8){
	                greenwich = objeto.textures[i].greenwich;
	                printf("\t Greenwich unitario: %LF, %LF, %LF \n", greenwich.x, greenwich.y, greenwich.z);
	                if (currentTypeObjectReading == 2  || currentTypeObjectReading == 8){
	                    norte = objeto.textures[i].north;
	                printf("\t Norte unitario: %LF, %LF, %LF \n", norte.x, norte.y, norte.z);
	                }

	            }
	        }
	    }
	}

	void createObjectFromData(long double *data, int whichObjectCreate, int quantityData, struct PlaneCut* planeCutsFound, struct Texture* texturesFound, struct DraftPlane *draftPlanesFound, int numberPlaneCuts, int numberTextures, int numberDraftPlanes){
	    /*whichObjectCreate indica qué objeto crear
	        Si está en 0, altera valores de la escena
	        Si está en 1, crea luces
	        Si está en 2, crea esferas
	        Si está en 3, crea polígonos, 
	        Si está en 4, crea cilindros
	        Si está en 5, crea conos.
	        data es el arreglo de valores del objeto a crear. Se asume que estará completo.*/
	    
	    switch(whichObjectCreate){
	        case 0: { //Escena 
	            if (debug == 1) {
	                printf("Insertando datos de escena \n");
	                printf("Reflexiones: %LF, Transparencia: %LF, Anti-aliasing: %LF \n", data[0], data[1], data[2]);
	                printf("Iluminación ambiente: %LF \n", data[3]);
	                printf("Plano de proyección (Xmin, Ymin) (Xmax, Ymax) : (%LF, %LF) (%LF, %LF) \n", data[4], data[5],data[6],data[7]);
	                printf("Resolución:  %LFx%LF \n", data[8], data[9]);
	                printf("Epsilon %LF \n", data[10]);
	                printf("Ojo: (%LF, %LF, %LF) \n", data[11],data[12],data[13]);
	                printf("Color background: (%LF, %LF, %LF) \n", data[14],data[15],data[16]);
	            }

	            maxAA = data[0];
	            maxReflection = data[1];
	            maxTransparency = data[2];
	            Ia = data[3];
	            Xmin = data[4];
	            Ymin = data[5];
	            Xmax = data[6];
	            Ymax = data[7];
	            Hres = data[8];
	            Vres = data[9];
	            e = data[10];

	            eye.x = data[11];
	            eye.y = data[12];
	            eye.z = data[13];

	            background.r = data[14];
	            background.g = data[15];
	            background.b = data[16];

	            Framebuffer[Hres][Vres];
	            Framebuffer = (struct Color **)malloc(Vres * sizeof(struct Color*));
	            for (int i = 0; i<Vres; i++){
	                Framebuffer[i] = (struct Color *)malloc(Hres*sizeof(struct Color));
	            }
	            
	            //printf("Scene data read correctly.\n");
	            return;
	            }
	        case 1: { //Luces
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
	        case 2: { //Esferas
	            if (debug == 1) {
	                printf("Insertando Esfera...");
	                printf("Pos esfera (%LF, %LF, %LF) \n", data[0],data[1],data[2]);
	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[3],data[4],data[5]);
	                printf("Radio esfera: %LF \n", data[6]);
	                printf("Esfera Kd: %LF \n", data[6]);
	                printf("Esfera Ka: %LF \n", data[8]);
	                printf("Esfera Kn: %LF \n", data[9]);
	                printf("Esfera Ks: %LF \n", data[10]);
	                printf("Color esfera (%LF, %LF, %LF) \n", data[11],data[12],data[13]);
	                
	            }
	            
	            struct Object polygon;
	            struct Color colorPolygon;
	            
	            struct Object esfera;
	            esfera.Xc=data[0];
	            esfera.Yc=data[1];
	            esfera.Zc=data[2];
	            esfera.o1 = data[3];
	            esfera.o2 = data[4];
	            esfera.o3 = data[5];
	            esfera.other=data[6]; //Radio
	            esfera.Kd=data[7]; 
	            esfera.Ka=data[8];
	            esfera.Kn=data[9];
	            esfera.Ks=data[10];
	            esfera.normalVector = sphereNormal;
	            esfera.intersectionFuncion = sphereIntersection;
	            struct Color colorSphere;
	            colorSphere.r = data[11];
	            colorSphere.g = data[12];
	            colorSphere.b = data[13];
	            esfera.color=colorSphere;

	            esfera.planeCuts = planeCutsFound;
	            esfera.numberPlaneCuts = numberPlaneCuts;
	            esfera.textures = texturesFound;
	            esfera.numberTextures = numberTextures;

	            Objects[objectIndex]=esfera;
	            if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            
	            objectIndex++;
	            //printf("Sphere processed \n \n");
	            return;
	            }
	        case 3: { //Polígonos

	            int vertexPolygonIndex = 0;
	            int numVertexesPolygon = (quantityData-10-12) / 3 + 1; //+1 para repetir el último vértice
	            
	            int inicioPlano = 10+(numVertexesPolygon-1)*3;
	            if (debug == 1){
	            	printf("quantityData: %i \n", quantityData);
	                printf("numVertexesPolygon: %i \n", numVertexesPolygon);
	                printf("%i inicioPlano \n", inicioPlano);
	                printf("Insertando polígono...");
	                printf("Color polígono (%LF, %LF, %LF) \n", data[0],data[1],data[2]);
	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[3],data[4],data[5]);
	                printf("Poligono Kd: %LF \n", data[6]);
	                printf("Poligono Ka: %LF \n", data[7]);
	                printf("Poligono Kn: %LF \n", data[8]);
	                printf("Poligono Ks: %LF \n", data[9]);
	                printf("Esquina inferior izquierda (%LF, %LF, %LF) \n", data[inicioPlano],data[inicioPlano+1],data[inicioPlano+2]);
	                printf("Esquina inferior derecha (%LF, %LF, %LF) \n", data[inicioPlano+3],data[inicioPlano+4],data[inicioPlano+5]);
	                printf("Esquina superior derecha (%LF, %LF, %LF) \n", data[inicioPlano+6],data[inicioPlano+7],data[inicioPlano+8]);
	                printf("Esquina superior izquierda (%LF, %LF, %LF) \n", data[inicioPlano+9],data[inicioPlano+10],data[inicioPlano+11]);
	            }

	            struct Point3D vertex;
	            struct Point2D squashedVertex;

	            //7 Elementos adicionales a los vertices
	            //Creo un objeto temporal

	            struct Object temp;
	            
	            
	            struct Vector x0y0z0;
	            x0y0z0.x = data[inicioPlano];
	            x0y0z0.y = data[inicioPlano+1];
	            x0y0z0.z = data[inicioPlano+2];
			  	struct Vector x1y1z1;
			  	x1y1z1.x = data[inicioPlano+3];
	            x1y1z1.y = data[inicioPlano+4];
	            x1y1z1.z = data[inicioPlano+5];
			  	struct Vector x2y2z2;
			  	x2y2z2.x = data[inicioPlano+6];
	            x2y2z2.y = data[inicioPlano+7];
	            x2y2z2.z = data[inicioPlano+8];
			  	struct Vector x3y3z3;
			  	x3y3z3.x = data[inicioPlano+9];
	            x3y3z3.y = data[inicioPlano+10];
	            x3y3z3.z = data[inicioPlano+11];

	            temp.points3D = malloc(sizeof(struct Point3D)*3);
	    
	            for (int i =0; i+10 < quantityData-12;){  //Los 12 del plano se ignoran
	                if(vertexPolygonIndex==3){
	                    break;
	                }
	                
	                vertex.x = data[10+i];
	                i++;
	                vertex.y = data[10+i];
	                i++;
	                vertex.z = data[10+i];
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
	            polygon.o1 = data[3];
	            polygon.o2 = data[4];
	            polygon.o3 = data[5];
	            polygon.Kd = data[6];
	            polygon.Ka = data[7];
	            polygon.Kn = data[8];
	            polygon.Ks = data[9];

	            polygon.pointAmount = numVertexesPolygon;
	            polygon.normalVector = polygonNormal;
	            polygon.intersectionFuncion = polygonIntersection;
	            polygon.retrieveTextureColor = planeTexture;
	            
	            long double u;
	            long double v;
	                    
	            long double maxA_B = max(fabs(polygon.Xc), fabs(polygon.Yc)); //maximo entre A y B
	            long double maxA_B_C = max(maxA_B, fabs(polygon.Zc)); //maximo entre los tres
	            
	            int choice = 0;
	            if(maxA_B_C == fabs(polygon.Xc)){ choice = 0; }
	            else if(maxA_B_C == fabs(polygon.Yc) ){ choice = 1; }
	            else if(maxA_B_C == fabs(polygon.Zc)){ choice = 2; } 

	            for (int i =0; i+10 < quantityData-12;){ //Los 12 son los 4 vértices del plano
	                vertex.x = data[10+i];
	                i++;
	                vertex.y = data[10+i];
	                i++;
	                vertex.z = data[10+i];
	                i++;
	                if(debug ==1){
	                	printf("Vertice: (%LF,%LF,%LF) \n", vertex.x, vertex.y, vertex.z);
	                }
	                if(choice == 0){ u = vertex.z; v = vertex.y; }
	                else if(choice  == 1){ u = vertex.x; v = vertex.z; }
	                else if(choice == 2){ u = vertex.x; v = vertex.y; } 
	                    
	                squashedVertex.u = u;
	                squashedVertex.v = v;
	                
	                polygon.points3D[vertexPolygonIndex]=vertex;
	                
	                polygon.points2D[vertexPolygonIndex]=squashedVertex;
	                vertexPolygonIndex++;
	            }

	            vertex.x = data[10];
	            vertex.y = data[11];
	            vertex.z = data[12];
	            
	            if(choice == 0){ u = vertex.z; v = vertex.y; }
	            else if(choice  == 1){ u = vertex.x; v = vertex.z; }
	            else if(choice == 2){ u = vertex.x; v = vertex.y; } 

	            squashedVertex.u = u;
	            squashedVertex.v = v;
	                
	            polygon.points3D[vertexPolygonIndex]=vertex;    
	            polygon.points2D[vertexPolygonIndex]=squashedVertex;    

	            polygon.planeCuts = planeCutsFound;
	            polygon.numberPlaneCuts = numberPlaneCuts;
	            polygon.textures = texturesFound;
	            polygon.numberTextures = numberTextures;

	            polygon.x0y0z0 = x0y0z0;
	            polygon.x1y1z1 = x1y1z1;
	            polygon.x2y2z2 = x2y2z2;
	            polygon.x3y3z3 = x3y3z3;


	            Objects[objectIndex] = polygon;
	            if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            
	            objectIndex++;

	            //printf("Polygon processed. \n \n");
	            return;
	            }
	        case 4: { //Cilindros
	            if (debug == 1) {
	                printf("Insertando cilindro...");

	                printf("Ancla: (%LF, %LF, %LF) \n", data[0], data[1],data[2]);
	                printf("Vector: (%LF, %LF, %LF) \n", data[3], data[4],data[5]);
	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[6],data[7],data[8]);
	                printf("Cilindro Radio: %LF \n", data[9]);
	                printf("Cilindro d1: %LF Cilindro d2: %LF \n", data[10],data[11]);
	                printf("Cilindro Kd: %LF \n", data[12]);
	                printf("Cilindro Ka: %LF \n", data[13]);
	                printf("Cilindro Kn: %LF \n", data[14]);
	                printf("Cilindro Ks: %LF \n", data[15]);
	                printf("RGB Cilindro: (%LF, %LF, %LF) \n", data[16], data[17],data[18]);
	    
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

	            cilinder.o1 = data[6];
	            cilinder.o2 = data[7];
	            cilinder.o3 = data[8];
	            cilinder.other = data[9];
	            cilinder.D1 = data[10];
	            cilinder.D2 = data[11];
	            cilinder.Kd = data[12];
	            cilinder.Ka = data[13];
	            cilinder.Kn = data[14];
	            cilinder.Ks = data[15];
	            cilinder.height = cilinder.D2 - cilinder.D1;
	            cilinder.normalVector = cilinderNormal;
	            cilinder.intersectionFuncion = cilinderIntersection;

	            struct Color cilinderColor;
	            cilinderColor.r = data[16];
	            cilinderColor.g = data[17];
	            cilinderColor.b = data[18];
	            cilinder.color = cilinderColor;
	            
	            cilinder.planeCuts = planeCutsFound;
	            cilinder.numberPlaneCuts = numberPlaneCuts;
	            cilinder.textures = texturesFound;
	            cilinder.numberTextures = numberTextures;

	            Objects[objectIndex] = cilinder;
	           if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            objectIndex++;
	            //printf("Cylinder processed.\n \n");
	            return;
	            }
	        case 5: { //Conos
	            if (debug == 1) {
	                printf("Insertando cono...");
	                printf("Ancla: (%LF, %LF, %LF) \n", data[0], data[1],data[2]);
	                printf("Vector: (%LF, %LF, %LF) \n", data[3], data[4],data[5]);
	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[6],data[7],data[8]);
	                printf("Cono k1: %LF COno k2: %LF \n", data[9],data[10]);
	                printf("Cono d1: %LF COno d2: %LF \n", data[11],data[12]);
	                printf("Cono Kd: %LF \n", data[13]);
	                printf("Cono Ka: %LF \n", data[14]);
	                printf("Cono Kn: %LF \n", data[15]);
	                printf("Cono Ks: %LF \n", data[16]);
	                printf("RGB Cono: (%LF, %LF, %LF) \n", data[17], data[18],data[19]);
	                
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

	            cone.o1 = data[6];
	            cone.o2 = data[7];
	            cone.o3 = data[8];
	            cone.K1 = data[9];
	            cone.K2 = data[10];
	            cone.D1 = data[11];
	            cone.D2 = data[12];
	            cone.height = cone.D2 - cone.D1;
	            cone.Kd = data[13];
	            cone.Ka = data[14];
	            cone.Kn = data[15];
	            cone.Ks = data[16];
	            cone.intersectionFuncion = coneIntersection;
	            cone.normalVector = coneNormal;

	            struct Color coneColor;
	            coneColor.r = data[17];
	            coneColor.g = data[18];
	            coneColor.b = data[19];
	            cone.color = coneColor;

	            cone.planeCuts = planeCutsFound;
	            cone.numberPlaneCuts = numberPlaneCuts;
	            cone.textures = texturesFound;
	            cone.numberTextures = numberTextures;
	            Objects[objectIndex] = cone;
	            if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            objectIndex++;
	            //printf("Cone processed.\n \n");
	            return;
	            }
	        case 6: { //Discos
	        	if (debug == 1) {
	                printf("Insertando disco...");

	                printf("Punto Central: (%LF, %LF, %LF) \n", data[0], data[1],data[2]);
	                printf("Normal: (%LF, %LF, %LF) \n", data[3], data[4],data[5]);
	                printf("Color: (%LF, %LF, %LF) \n", data[6], data[7],data[8]);
	                
	                printf("Disco Radio: %LF \n", data[9]);
	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[10],data[11],data[12]);
	                printf("DIsco Kd: %LF \n", data[13]);
	                printf("DIsco Ka: %LF \n", data[14]);
	                printf("DIsco Kn: %LF \n", data[15]);
	                printf("Disco Ks: %LF \n", data[16]);
	                printf("Esquina inferior izquierda (%LF, %LF, %LF) \n", data[17],data[18],data[19]);
	                printf("Esquina inferior derecha (%LF, %LF, %LF) \n", data[20],data[21],data[22]);
	                printf("Esquina superior derecha (%LF, %LF, %LF) \n", data[23],data[24],data[25]);
	                printf("Esquina superior izquierda (%LF, %LF, %LF) \n",data[26],data[27],data[28]);
	            }


	            struct Object disco;
	            
	            
	            struct Vector x0y0z0;
	            x0y0z0.x = data[17];
	            x0y0z0.y = data[18];
	            x0y0z0.z = data[19];
			  	struct Vector x1y1z1;
			  	x1y1z1.x = data[20];
	            x1y1z1.y = data[21];
	            x1y1z1.z = data[22];
			  	struct Vector x2y2z2;
			  	x2y2z2.x = data[23];
	            x2y2z2.y = data[24];
	            x2y2z2.z = data[25];
			  	struct Vector x3y3z3;
			  	x3y3z3.x = data[26];
	            x3y3z3.y = data[27];
	            x3y3z3.z = data[28];


	            disco.Xc = data[0];
	            disco.Yc = data[1];
	            disco.Zc = data[2];

	            disco.intersectionFuncion = discIntersection;
	            disco.normalVector = discNormal;
	            disco.retrieveTextureColor = planeTexture;

	            struct Vector puntoCentral;
	            puntoCentral.x = data[0];
	            puntoCentral.y = data[1];
	            puntoCentral.z = data[2];

	            struct Vector normalNotNormalized;
	            normalNotNormalized.x = data[3];
	            normalNotNormalized.y = data[4];
	            normalNotNormalized.z = data[5];

	            struct Color colorDisco;
	            colorDisco.r = data[6];
	            colorDisco.g = data[7];
	            colorDisco.b = data[8];
	            disco.color = colorDisco;

	            long double dPlano= whatsTheDGeneral(normalNotNormalized, puntoCentral);
	            dPlano = dPlano / getNorm(normalNotNormalized);
	            disco.extraD = dPlano;
	            normalNotNormalized = normalize(normalNotNormalized);

	            disco.directionVector = normalNotNormalized;
	            disco.other = data[9];
	            disco.o1 = data[10];
	            disco.o2 = data[11];
	            disco.o3 = data[12];
	            disco.Kd = data[13];
	            disco.Ka = data[14];
	            disco.Kn = data[15];
	            disco.Ks = data[16];

	            disco.planeCuts = planeCutsFound;
	            disco.numberPlaneCuts = numberPlaneCuts;
	            disco.textures = texturesFound;
	            disco.numberTextures = numberTextures;
	            Objects[objectIndex] = disco; 
	            if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            objectIndex++;
	            return;
	            }
	        case 7: { //Elipses
	        	if (debug == 1) {
	                printf("Insertando Elipses...");

	                printf("Foco 1: (%LF, %LF, %LF) \n", data[0], data[1],data[2]);
	                printf("Foco 2: (%LF, %LF, %LF) \n", data[3], data[4],data[5]);
	                printf("Normal no normalizada: (%LF, %LF, %LF) \n", data[6], data[7],data[8]);
	                printf("Color: (%LF, %LF, %LF) \n", data[9], data[8],data[9]);
	                
	                printf("K del elipse: %LF \n", data[12]);
	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[13],data[14],data[15]);
	                printf("Elipse Kd: %LF \n", data[16]);
	                printf("Elipse Ka: %LF \n", data[17]);
	                printf("Elipse Kn: %LF \n", data[18]);
	                printf("Elipse Ks: %LF \n", data[19]);
	                printf("Esquina inferior izquierda (%LF, %LF, %LF) \n", data[20],data[21],data[22]);
	                printf("Esquina inferior derecha (%LF, %LF, %LF) \n", data[23],data[24],data[25]);
	                printf("Esquina superior derecha (%LF, %LF, %LF) \n", data[26],data[27],data[28]);
	                printf("Esquina superior izquierda (%LF, %LF, %LF) \n",data[29],data[30],data[31]);
	            }


	            struct Object elipse;
	            
	            
	            struct Vector x0y0z0;
	            x0y0z0.x = data[20];
	            x0y0z0.y = data[21];
	            x0y0z0.z = data[22];
			  	struct Vector x1y1z1;
			  	x1y1z1.x = data[23];
	            x1y1z1.y = data[24];
	            x1y1z1.z = data[25];
			  	struct Vector x2y2z2;
			  	x2y2z2.x = data[26];
	            x2y2z2.y = data[27];
	            x2y2z2.z = data[28];
			  	struct Vector x3y3z3;
			  	x3y3z3.x = data[29];
	            x3y3z3.y = data[30];
	            x3y3z3.z = data[31];



	            elipse.intersectionFuncion = elipseIntersection;
	            elipse.normalVector = elipseNormal;
	            elipse.retrieveTextureColor = planeTexture;

	            //Se crea un punto para obtener el D de la normal
	            struct Vector foco1;
	            foco1.x = data[0];
	            foco1.y = data[1];
	            foco1.z = data[2];
	            //Se ingresa el D1 tal y como se solciita en las funciones de intersección
	            elipse.Xc = data[0];
	            elipse.Yc = data[1];
	            elipse.Zc = data[2];
				//Se ingresa el D2 tal y como se solciita en las funciones de intersección
	            elipse.Xother = data[3];
	            elipse.Yother = data[4];
	            elipse.Zother = data[5];


	            struct Vector normalNotNormalized;
	            normalNotNormalized.x = data[6];
	            normalNotNormalized.y = data[7];
	            normalNotNormalized.z = data[8];


	            struct Color colorElipse;
	            colorElipse.r = data[9];
	            colorElipse.g = data[10];
	            colorElipse.b = data[11];
	            elipse.color = colorElipse;

	            long double dPlano= whatsTheDGeneral(normalNotNormalized, foco1);
	            dPlano = dPlano / getNorm(normalNotNormalized);
	            elipse.extraD = dPlano;
	            normalNotNormalized = normalize(normalNotNormalized);

	            elipse.directionVector = normalNotNormalized;
	            elipse.other = data[12];
	            elipse.o1 = data[13];
	            elipse.o2 = data[14];
	            elipse.o3 = data[15];
	            elipse.Kd = data[16];
	            elipse.Ka = data[17];
	            elipse.Kn = data[18];
	            elipse.Ks = data[19];

	            elipse.planeCuts = planeCutsFound;
	            elipse.numberPlaneCuts = numberPlaneCuts;
	            elipse.textures = texturesFound;
	            elipse.numberTextures = numberTextures;
	            Objects[objectIndex] = elipse; 
	            if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            objectIndex++;
	            return;
	        	}
	        case 8: { //Cuadráticas


	        	if (debug == 1) {
	                printf("Insertando Cuadrática...");

	                
	                printf("Coeficientes: \n\t A: %LF \n\t B: %LF \n\t C: %LF\n\t D: %LF\n\t E: %LF \n\t F: %LF\n\t G: %LF \n\t H: %LF  \n\t I: %LF \n", data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9]);
	                printf("Constante K: %LF \n", data[10]);

	                

	                printf("o1:  %LF, o2: %LF, o3: %LF \n", data[11],data[12],data[13]);
	                printf("Elipse Kd: %LF \n", data[14]);
	                printf("Elipse Ka: %LF \n", data[15]);
	                printf("Elipse Kn: %LF \n", data[16]);
	                printf("Elipse Ks: %LF \n", data[17]);
	                printf("Color: (%LF, %LF, %LF) \n", data[18], data[19],data[20]);

	            }


	        	struct Object cuadratica;



	        	cuadratica.A = data[0];
	        	cuadratica.B = data[1];
	        	cuadratica.C = data[2];
	        	cuadratica.D = data[3];
	        	cuadratica.E = data[4];
	        	cuadratica.F = data[5];
	        	cuadratica.G = data[6];
	        	cuadratica.H = data[7];
	        	cuadratica.I = data[8];
	        	cuadratica.J = data[9];
	        	//El K
	        	cuadratica.other = data[10];


	        	cuadratica.o1 = data[11];
	        	cuadratica.o2 = data[12];
	        	cuadratica.o3 = data[13];

	        	cuadratica.Kd = data[14];
	            cuadratica.Ks = data[15];
	            cuadratica.Kn = data[16];
	            cuadratica.Ks = data[17];

	            struct Color colorCuadratica;
	            colorCuadratica.r = data[18];
	            colorCuadratica.g = data[19];
	            colorCuadratica.b = data[20];
	            cuadratica.color = colorCuadratica;

	            /*
	            cuadratica.intersectionFuncion = quadraticIntersection;
	            cuadratica.normalVector = quadraticNormal;
				*/
	            cuadratica.planeCuts = planeCutsFound;
	            cuadratica.numberPlaneCuts = numberPlaneCuts;
	            cuadratica.textures = texturesFound;
	            cuadratica.numberTextures = numberTextures;
	            Objects[objectIndex] = cuadratica; 
	            if (debug == 1){
	            	printPlaneCuts(Objects[objectIndex]);
	            	printTextures(Objects[objectIndex], whichObjectCreate);
	            }
	            objectIndex++;
	            return;
	        	}
	    	}
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

	//Quita el salto de línea de un string.
	void strip(char *s) {
	    char *p2 = s;
	    while(*s != '\0') {
	        if(*s != '\t' && *s != '\n') {
	            *p2++ = *s++;
	        } else {
	            ++s;
	        }
	    }
	    *p2 = '\0';}

	char *obtainFilenameTexture(char stringLine[]){
	    /*Devuelve el filename de una textura a partir de una línea
	    Filename = "mapas/textura.jpg" */
	    char *token = malloc(sizeof(char)*200);
	    char *search = "=";
	    // Token will point to Filename
	    token = strtok(stringLine, search);
	    // Token will point to "mapas/textura.jpg"
	    token = strtok(NULL, search);
	    strip(token);
	    printf("%s \n", token);
	    return token;}

	struct PlaneCut *readPlaneCuts(long int pos, int *numberPlanes, long int *posAfterReading){
	    //Para este punto ya leyó Planos_Corte: 
	    //ENtonces comenzamos con NumberPlanes
	    //Se lee el número de planos
	    char temporalBuffer[300]; //Aquí se guardará lo leído cada línea
	    struct PlaneCut *planeCutsFound = NULL;
	    long double *datosPlanos; //Siempre serán 6 valores que componen a un plano de corte
	    int indexPlaneCut = -1; //Cual plano estamos leyendo? 
	    FILE *file;
	    if (file = fopen(escenaFile, "r")){
	        fseek(file, pos, SEEK_SET);
	    while (fgets(temporalBuffer, 300, file)!=NULL){ //Mientras el archivo siga teniendo algo
	        //printf("Dentro de readPlaneCuts: %s \n",temporalBuffer);
	            if (temporalBuffer[0] == '\n'){
	                continue;
	            }
	            if (strstr(temporalBuffer, "#")!=NULL){
	                continue;
	            }
	            if (strstr(temporalBuffer, "NumberPlanes")!=NULL){
	                long double numberPlanes= obtainSingleValueFromLine(temporalBuffer);
	                planeCutsFound = malloc(sizeof(struct PlaneCut)*numberPlanes);
	                continue;
	            }else if (strstr(temporalBuffer, "Plano_")!=NULL){
	                indexPlaneCut++;
	                continue;
	            }else if (strstr(temporalBuffer, "END_Planos")!=NULL){
	                *numberPlanes = indexPlaneCut+1;
	                *posAfterReading = ftell(file);
	                return planeCutsFound;
	            }else if (strstr(temporalBuffer, "Punto")!=NULL){
	                datosPlanos = obtainPointFromString(temporalBuffer);

	                struct Vector temp;
	                temp.x = datosPlanos[0];
	                temp.y = datosPlanos[1];
	                temp.z = datosPlanos[2];
	                planeCutsFound[indexPlaneCut].point = temp;
	                free(datosPlanos);
	                continue;
	            }else if (strstr(temporalBuffer, "Normal")!=NULL){
	                datosPlanos = obtainPointFromString(temporalBuffer);
	                struct Vector temp;
	                temp.x = datosPlanos[0];
	                temp.y = datosPlanos[1];
	                temp.z = datosPlanos[2];
	                long double dEquation = whatsTheDGeneral(temp,planeCutsFound[indexPlaneCut].point);
	                dEquation = dEquation/getNorm(temp);
	                temp = normalize(temp); 

	                planeCutsFound[indexPlaneCut].normal = temp;
	                planeCutsFound[indexPlaneCut].d = dEquation;
	                free(datosPlanos);
	                continue;
	            }
	        continue;
	        }
	    }
	}

	struct Color **getTexels(char* pFile, int* hRes, int* vRes){
	    int counter, x, y, i, j;
	    char dump[100];
	    time_t t;
	    struct Color **temp;
	    srand((unsigned) time(&t));
	    FILE* file;
	    if (file = fopen(pFile, "r")) {
	        
	        //Read width and height
	        for (i = 0; i < 11; ++i) { 
	            fscanf(file, "%s", &dump[0]);
	            if (i==8||i==9) {
	                if(i == 8){ //Hres es el octavo string del archivo
	                    sscanf(&dump[0], "%i", hRes);
	                }else{
	                    sscanf(&dump[0], "%i", vRes);
	                }
	                
	                
	            }
	            //printf("%s \n", dump);

	        }
	        //printf("Textura en %s \n HRes: %i \n VRes: %i \n", pFile, *hRes, *vRes);
	        temp =   malloc (sizeof(struct Color*)*(*hRes));
	       	for (int r = 0; r < (*hRes);r++){
	       		temp[r] = malloc(sizeof(struct Color)*(*vRes));
	       	}
	        if (temp==NULL){
	            printf("Devolvió NULL \n");
	        }
	        char temporalBuffer[2000];
	        int i = 0,x = 0, y = 0, counter = 0;
	        while (fgets(temporalBuffer, 200, file)!=NULL){ //Mientras el archivo siga teniendo algo
	        	//strip(temporalBuffer);
	        	if (temporalBuffer[0] == '\n'){
	                continue;
	            }
	            struct Color texel;
	            long double number;
	            sscanf(temporalBuffer, "%LF", &number);
	            
	            //strip(temporalBuffer);
	
	            int xs = *hRes - x - 1;
	            if (i == 0){
	                temp[y][xs].r = (number)/255;
	            }else if(i == 1){
	                temp[y][xs].g = (number)/255;
	            }else if(i == 2){
	                temp[y][xs].b = (number)/255;
	            }
	            i = (i+1)%3;
	            if (i==0){
	                y = (y+1);
	                y = y%(*vRes);
	                if (y==0){
	                    x = (x+1);
	                    x= x%(*hRes);
	                    if (x == 0){
	                        break;
	                    }
	                }
	            }
	            counter=1;     
	        }
	        y = 0;
	        x = 0;
	            while (counter !=5){
	            	//printf("Color leído: (%LF, %LF, %LF) \n",temp[x][y].r*255, temp[x][y].g*255, temp[x][y].b*255 );
	            	counter++;
	            	y++;
	            }
	        
	    }
	    else {
	    	(*vRes) = 128;
	        (*hRes) = 128;
	        temp =   malloc (sizeof(struct Color*)*128);
	       	for (int r = 0; r < (*hRes);r++){
	       		temp[r] = malloc(sizeof(struct Color)*128);
	       	}
	        if (temp==NULL){
	            printf("Devolvió NULL \n");
	        }
	        printf("La textura de %s no pudo abrirse. Se sustituirá por estática\n", pFile);
	        for(x = 0; x < 128; x++){
	            for (y = 0; y < 128; y++){
	            	//int xs = 128 - x - 1;
	                struct Color estatica;
	                estatica.r = ((long double)(rand() % 255))/255;
	                estatica.g = ((long double)(rand() % 255))/255;
	                estatica.b = ((long double)(rand() % 255))/255;
	                temp[x][y] = estatica;
	            }
	        }
	    }
	    return temp;
	}

	struct Texture *readTextures(int currentTypeReading, long int pos, int *numberTextures, long int *posAfterReading){
	    //Para este punto ya leyó Textures: 
	    //ENtonces comenzamos con numberTextures
	    //Se lee el número de planos
	    char temporalBuffer[300]; //Aquí se guardará lo leído cada línea
	    struct Texture *texturesFound = NULL;
	    long double *datosTexture; //Siempre serán 6 valores MÁXIMO que leerá de cada textura
	    int indexTexture = -1; //Cual plano estamos leyendo? 
	    FILE *file;
	    //printf("Dentro de readTextures: Tipo de objeto = %i \n",currentTypeReading);
	    if (file = fopen(escenaFile, "r")){
	        fseek(file, pos, SEEK_SET);
	    while (fgets(temporalBuffer, 300, file)!=NULL){ //Mientras el archivo siga teniendo algo
	        //printf("%s", temporalBuffer);
	            if (temporalBuffer[0] == '\n'){
	                continue;
	            }
	            if (temporalBuffer[0] == '\t'){
	                continue;
	            }
	            if (strstr(temporalBuffer, "#")!=NULL){
	                continue;
	            }
	            if (strstr(temporalBuffer, "NumberTextures")!=NULL || strstr(temporalBuffer, "NumberTexturas")!=NULL){
	                long double numberTextures= obtainSingleValueFromLine(temporalBuffer);
	                texturesFound = malloc(sizeof(struct Texture)*numberTextures);
	                continue;
	            }else if (strstr(temporalBuffer, "Texture_")!=NULL || strstr(temporalBuffer, "Textura_")!=NULL ){
	                indexTexture++;
	                continue;
	            }else if (strstr(temporalBuffer, "END_Textures")!=NULL || strstr(temporalBuffer, "END_Texturas")!=NULL){
	                *numberTextures = indexTexture+1;
	                *posAfterReading = ftell(file);
	                return texturesFound;
	            }else if (strstr(temporalBuffer, "Filename")!=NULL || strstr(temporalBuffer, "filename")!=NULL){
	                char* filename = obtainFilenameTexture(temporalBuffer);
	                int hRes, vRes;
	                texturesFound[indexTexture].filename = filename;
	                struct Color **textureMap=getTexels(texturesFound[indexTexture].filename, &hRes, &vRes);

	                /*for(int x = 0; x < 128; x++){
	            		for ( int y = 0; y < 128; y++){
	                		printf("Color leído fuera:(%LF, %LF, %LF) \n",textureMap[x][y].r, textureMap[x][y].g, textureMap[x][y].b);
	                	}
	                }*/

	                texturesFound[indexTexture].textureMap = textureMap;
	                texturesFound[indexTexture].hRes = hRes;
					texturesFound[indexTexture].vRes = vRes;
	                //printf("FIlename Textura %s  ggg\n", texturesFound[indexTexture].filename);
	                //free(filename);
	                continue;
	            }

	            if (currentTypeReading == 2 || currentTypeReading == 4 || currentTypeReading == 5 || currentTypeReading == 8){
	                //COnos, esferas, cilindros o cuadráticas
	                if (strstr(temporalBuffer, "Greenwich")!=NULL){
	                    datosTexture = obtainPointFromString(temporalBuffer);
	                    struct Vector temp;
	                    temp.x = datosTexture[0];
	                    temp.y = datosTexture[1];
	                    temp.z = datosTexture[2];
	                    temp = normalize(temp);  //NORMALIZO
	                    texturesFound[indexTexture].greenwich = temp;
	                    free(datosTexture);
	                    continue;
	                }
	                if (currentTypeReading == 2 ||currentTypeReading == 8){
	                    //Esferas y cuadráticas
	                    if(strstr(temporalBuffer, "Norte")!=NULL || strstr(temporalBuffer, "North")!=NULL){
	                        datosTexture = obtainPointFromString(temporalBuffer);
	                        struct Vector temp;
	                        temp.x = datosTexture[0];
	                        temp.y = datosTexture[1];
	                        temp.z = datosTexture[2];
	                        temp = normalize(temp); //NORMALIZO
	                        texturesFound[indexTexture].north = temp;
	                        free(datosTexture);
	                        continue;
	                    }
	                }
	            }
	        continue;
	        }
	    }
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
	            if ((*counterValueSegment) >= 0 && (*counterValueSegment) <= 10){ 
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++;
	                *numberValuesRead = 1;
	                return values;
	            }else if ((*counterValueSegment) >= 11 && (*counterValueSegment) <= 12){

	                long double *point = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = point[0];
	                values[1] = point[1];
	                values[2] = point[2];
	                *numberValuesRead = 3;
	                free (point);
	                if((*counterValueSegment) == 12 ){
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
	        case 2: //Esferas
	            if ((*counterValueSegment)== 0 || (*counterValueSegment) ==9){ //Nombre o posición de la esfera
	                long double *positionSphere = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = positionSphere[0];
	                values[1] = positionSphere[1];
	                values[2] = positionSphere[2];

	                free (positionSphere);
	                *numberValuesRead = 3;
	                if((*counterValueSegment)==9){
	                    (*counterValueSegment)=0;
	                }else{
	                    (*counterValueSegment)++;
	                }
	                return values;
	            }else if((*counterValueSegment) >= 1 && (*counterValueSegment )<= 8 ){ 
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
	            }else if((*counterValueSegment) >= 1 && (*counterValueSegment) <= 7){ //Ks
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++;
	                *numberValuesRead = 1;
	                return values;
	            }else if((*counterValueSegment) == 8){ //Poligono Vertice
	                if (strstr(lineRead,"END_Vertices")!=NULL){
	                    //(*counterValueSegment)=0;
	                    (*counterValueSegment)++; //Seguimos leyendo el plano
	                    //printf("Se terminó de leer el plano \n");
	                    return values;
		        	}else{
		                values = malloc(sizeof(long double)*3); 
		                long double *vertexPolygon = obtainPointFromString(lineRead);
		                values[0] = vertexPolygon[0];
		                values[1] = vertexPolygon[1];
		                values[2] = vertexPolygon[2];
		                free (vertexPolygon);
		                //printf("Tripleta: %LF, %LF, %LF \n", values[0], values[1], values[2]);
		                *numberValuesRead = 3;
		                if ((*counterValueSegment) == 12){ //Terminó de leer el plano
		                    //printf("Se terminó de leer el plano \n");
		                    (*counterValueSegment)=0;
		                    return values;
		                }
		                if((*counterValueSegment) != 8){ //Lee el plano 
		                    (*counterValueSegment)++;
		                }
		                //printf("counterValueSegment tras asignar RGB %i \n",(*counterValueSegment));
		                return values;
		            }
		        }else if((*counterValueSegment) >= 9 ){
		        	values = malloc(sizeof(long double)*3); 
		                long double *vertexPolygon = obtainPointFromString(lineRead);
		                values[0] = vertexPolygon[0];
		                values[1] = vertexPolygon[1];
		                values[2] = vertexPolygon[2];
		                free (vertexPolygon);
		                //printf("Tripleta: %LF, %LF, %LF \n", values[0], values[1], values[2]);
		                *numberValuesRead = 3;
		                if ((*counterValueSegment) == 12){ //Terminó de leer el plano
		                    //printf("Se terminó de leer el plano \n");
		                    (*counterValueSegment)=0;
		                    return values;
		                }
		                (*counterValueSegment)++;
		                return values;
		        }
	        case 4: //Cilindros
	            if(*counterValueSegment == 0 || *counterValueSegment == 1|| *counterValueSegment == 12){
	                /*Lee el ancla, color del cilindro o el vector del cilindro. Tripleta*/
	                long double *positionCilinder = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = positionCilinder[0];
	                values[1] = positionCilinder[1];
	                values[2] = positionCilinder[2];
	                //memcpy(values, positionLight, 3);
	                //printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
	                if(*counterValueSegment == 12){
	                    (*counterValueSegment) = 0;
	                }else{
	                    (*counterValueSegment)++;
	                }
	                    
	                *numberValuesRead = 3;
	                free (positionCilinder);
	                //printf("%LF %LF %LF \n", values[0], values[1],values[2]);
	                return values;
	            }else if(*counterValueSegment >= 2 && *counterValueSegment <= 11){
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++;
	                *numberValuesRead = 1;
	                return values;
	                /*Lee radio o d1 o d2 o Kd o Kd o Ka o Kn o Ks*/
	            }
	        case 5: //Conos
	            if(*counterValueSegment == 0 || *counterValueSegment == 1|| *counterValueSegment == 13){
	                /*Lee el ancla, color del cono o el vector del cono. Tripleta*/
	                long double *positionCone = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = positionCone[0];
	                values[1] = positionCone[1];
	                values[2] = positionCone[2];
	                //memcpy(values, positionLight, 3);
	                //printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
	                if(*counterValueSegment == 13){
	                    (*counterValueSegment) = 0;
	                }else{
	                    (*counterValueSegment)++;
	                }
	                *numberValuesRead = 3;
	                free (positionCone);
	                //printf("%LF %LF %LF \n", values[0], values[1],values[2]);
	                return values;
	            }else if(*counterValueSegment >= 2 && *counterValueSegment <= 12){
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++;
	                *numberValuesRead = 1;
	                return values;
	                /*Lee radio o k1 i k2 o d1 o d2 o Kd o Kd o Ka o Kn o Ks*/
	            }
	        case 6:  //Discos
	        	if((*counterValueSegment >= 0 && *counterValueSegment <= 2) || (*counterValueSegment >= 11 && *counterValueSegment <= 14)){
	        		//Lee centro, normal, color y vertices del plano de la textura
	                long double *tripleta = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = tripleta[0];
	                values[1] = tripleta[1];
	                values[2] = tripleta[2];
	                //memcpy(values, positionLight, 3);
	                //printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
	                if(*counterValueSegment == 14){
	                    (*counterValueSegment) = 0;
	                }else{
	                    (*counterValueSegment)++;
	                }
	                *numberValuesRead = 3;
	                free (tripleta);
	                //printf("%LF %LF %LF \n", values[0], values[1],values[2]);
	                return values;
	            }else if((*counterValueSegment >= 3 && *counterValueSegment <= 10)){
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++;
	                *numberValuesRead = 1;
	                return values;
	                /*Lee radio o o1, o2, o3 o o Kd o Kd o Ka o Kn o Ks*/
	            }
	        case 7: //Elipses
	        	if((*counterValueSegment >= 0 && *counterValueSegment <= 3) || (*counterValueSegment >= 12 && *counterValueSegment <= 15)){
	        		//Lee foco1, foco2, normal, color y vertices del plano de la textura
	        		printf("lineRead");
	                long double *tripleta = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = tripleta[0];
	                values[1] = tripleta[1];
	                values[2] = tripleta[2];
	                //memcpy(values, positionLight, 3);
	                //printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
	                if(*counterValueSegment == 15){
	                    (*counterValueSegment) = 0;
	                }else{
	                    (*counterValueSegment)++;
	                }
	                *numberValuesRead = 3;
	                free (tripleta);
	                //printf("%LF %LF %LF \n", values[0], values[1],values[2]);
	                return values;
	            }else if((*counterValueSegment >= 4 && *counterValueSegment <= 11)){
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++;
	                *numberValuesRead = 1;
	                return values;
	                /*Lee radio o1, o2, o3 k1 i k2 o o Kd o Kd o Ka o Kn o Ks*/
	            }

	        case 8: //CUadráticas
	        	if((*counterValueSegment) == 18){
	        		//Lee color
	        		printf("lineRead");
	                long double *tripleta = obtainPointFromString(lineRead);
	                values = malloc(sizeof(long double)*3);
	                values[0] = tripleta[0];
	                values[1] = tripleta[1];
	                values[2] = tripleta[2];
	                //memcpy(values, positionLight, 3);
	                //printf("Pos luz leída (%LF, %LF, %LF) \n", values[0],values[1],values[2]);
	                (*counterValueSegment) = 0; //Termina al leer el color
	                *numberValuesRead = 3;
	                free (tripleta);
	                //printf("%LF %LF %LF \n", values[0], values[1],values[2]);
	                return values;
	            }else if((*counterValueSegment >= 0 && *counterValueSegment <= 17) ){
	                values = malloc(sizeof(long double));
	                values[0] = obtainSingleValueFromLine(lineRead);
	                (*counterValueSegment)++; 
	                *numberValuesRead = 1;
	                return values;
	                /*Lee coeficientes, los k's, los o's y la constante K.*/
	            }

	    }
	}

	int plainCutsFound(long int pos){
	    //Retorna 1 si encuentra que al objeto le siguen planos de corte
	    char temporalBuffer[300]; //Aquí se guardará lo leído cada línea
	    
	    FILE *file;
	    if (file = fopen(escenaFile, "r")){
	        fseek(file, pos, SEEK_SET);
	        while (fgets(temporalBuffer, 300, file)!=NULL){ //Mientras el archivo siga teniendo algo
	            //printf("Dentro de plainCutsFound:  %s \n", temporalBuffer);
	            if (temporalBuffer[0] == '\n'){
	                continue;
	            }
	            if (temporalBuffer[0] == '\t'){
	                continue;
	            }
	            if (strstr(temporalBuffer, "#")!=NULL){
	                continue;
	            }else if (strstr(temporalBuffer, "Planos_Corte:")!=NULL){
	                return 1; //LO ENCONTRÓ. PERFECTO.
	            }else if (strstr(temporalBuffer, "Texturas:")!=NULL){
	            //Entra en state = 2. Esferas
	                return 0;
	            }else if (strstr(temporalBuffer, "Planos_Calado:")!=NULL){
	            //Entra en state = 2. Esferas
	                return 0;
	            }else if (strstr(temporalBuffer, "Sphere_Object")!=NULL){
	            //Entra en state = 2. Esferas
	                return 0;
	            }else if (strstr(temporalBuffer, "Polygon_Object")!=NULL){
	            //Entra en state = 3. Poligonos
	                //printf("Leyó polígonos");
	                return 0;
	            }else if (strstr(temporalBuffer, "Cylinder_Object")!=NULL){
	            //Entra en state = 4. Cilindros
	                return 0;
	            }else if (strstr(temporalBuffer, "Cone_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Disc_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Elipse_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Quadratic_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Scene_Data")!=NULL){
	            //Entra en state = 0. Escena
	                return 0;
	            }else if (strstr(temporalBuffer,"Light_Object")!=NULL){
	                return 0;
	            }
	        }

	    }
	    return 0;}

	int texturesFound(long int pos){
	    //Retorna 1 si encuentra que al objeto le siguen texturas
	    char temporalBuffer[300]; //Aquí se guardará lo leído cada línea
	    
	    FILE *file;
	    if (file = fopen(escenaFile, "r")){
	        fseek(file, pos, SEEK_SET);
	        while (fgets(temporalBuffer, 300, file)!=NULL){ //Mientras el archivo siga teniendo algo
	            //printf("Dentro de texturesFound:  %s \n", temporalBuffer);
	            if (temporalBuffer[0] == '\n'){
	                continue;
	            }
	            if (temporalBuffer[0] == '\t'){
	                continue;
	            }
	            if (strstr(temporalBuffer, "#")!=NULL){
	                continue;
	            }else if (strstr(temporalBuffer, "Texturas:")!=NULL || strstr(temporalBuffer, "Textures:")!=NULL){
	            //Entra en state = 2. Esferas
	                return 1;
	            }else if (strstr(temporalBuffer, "Planos_Calado:")!=NULL){
	            //Entra en state = 2. Esferas
	                return 0;
	            }else if (strstr(temporalBuffer, "Sphere_Object")!=NULL){
	            //Entra en state = 2. Esferas
	                return 0;
	            }else if (strstr(temporalBuffer, "Polygon_Object")!=NULL){
	            //Entra en state = 3. Poligonos
	                //printf("Leyó polígonos");
	                return 0;
	            }else if (strstr(temporalBuffer, "Cylinder_Object")!=NULL){
	            //Entra en state = 4. Cilindros
	                return 0;
	            }else if (strstr(temporalBuffer, "Cone_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Disc_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Elipse_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Quadratic_Object")!=NULL){
	            //Entra en state = 5. Conos
	                return 0;
	            }else if (strstr(temporalBuffer, "Scene_Data")!=NULL){
	            //Entra en state = 0. Escena
	                return 0;
	            }else if (strstr(temporalBuffer,"Light_Object")!=NULL){
	                return 0;
	            }
	        }

	    }
	    return 0;
	}
// ===============================================================

// Leer archivos con la escena ===================================
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
	    
	    char temporalBuffer[300]; //Aquí se guardará lo leído cada línea
	    long double *valuesRead;  //Se guarda los valores del objeto para finalmente crearlo.
	    int indexValuesRead = 0;  //Pos de valuesRead 
	    int currentTypeObjectReading = 1;

	    struct PlaneCut *arrayPlaneCuts = NULL;
	    struct Texture *arrayTextures = NULL;
	    struct DraftPlane *arrayDraftPlanes = NULL;

	    FILE* file; //archivo
	    if (file = fopen(escenaFile, "r")){

	        while (fgets(temporalBuffer, 300, file)!=NULL){ //Mientras el archivo siga teniendo algo
	            //printf("%s \n",temporalBuffer);
	            if (temporalBuffer[0] == '\n'){
	                continue;
	            }
	            if (temporalBuffer[0] == '\t'){
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
	                valuesRead = malloc(sizeof(long double)*22);  
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
	                valuesRead=malloc(sizeof(long double)*22); //Esferas siempre serán 11 valores
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
	                valuesRead=malloc(sizeof(long double)*2000000); //Polígonos max size 
	                currentTypeObjectReading=3;
	                //printf("%s",temporalBuffer);
	                continue;
	            }else if (strstr(temporalBuffer, "Cylinder_Object")!=NULL){
	            //Entra en state = 4. Cilindros
	                state = 4;
	                counterValueSegment = 0;
	                indexValuesRead=0;
	                free(valuesRead);
	                valuesRead = malloc(sizeof(long double)*55);
	                currentTypeObjectReading=4;
	                //printf("%s",temporalBuffer);
	                continue;
	            }else if (strstr(temporalBuffer, "Cone_Object")!=NULL){
	            //Entra en state = 5. Conos
	                state = 5;
	                counterValueSegment = 0;
	                indexValuesRead=0;
	                free(valuesRead);
	                valuesRead = malloc(sizeof(long double)*55);
	                currentTypeObjectReading=5;
	                //printf("%s",temporalBuffer);
	                continue;
	            }else if (strstr(temporalBuffer, "Disc_Object")!=NULL){
	                state = 6;
	                counterValueSegment = 0;
	                indexValuesRead=0;
	                free(valuesRead);
	                valuesRead = malloc(sizeof(long double)*55);
	                currentTypeObjectReading=6;
	                printf("%s",temporalBuffer);
	                continue;
	            }else if (strstr(temporalBuffer, "Elipse_Object")!=NULL){
	                state = 7;
	                counterValueSegment = 0;
	                indexValuesRead=0;
	                free(valuesRead);
	                valuesRead = malloc(sizeof(long double)*55);
	                currentTypeObjectReading=7;
	                //printf("%s",temporalBuffer);
	                continue;
	            }else if (strstr(temporalBuffer, "Quadratic_Object")!=NULL){
	                state = 8;
	                counterValueSegment = 0;
	                indexValuesRead=0;
	                free(valuesRead);
	                valuesRead = malloc(sizeof(long double)*55);
	                currentTypeObjectReading=8;
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
	             se leyeron los datos de dicho objeto/luz, ergo, se procede a revisar
	             si existen planos de corte, texturas o planos de calado definidos*/


	            if (counterValueSegment == 0){
	            	//printf("currentTypeObjectReading %i \n",currentTypeObjectReading );
	                int numberPlaneCuts = 0;
	                int numberTextures = 0;
	                int numberDraftPlanes = 0;
	                long int pos;
	                pos = ftell(file);
	                int areTherePlainCuts = plainCutsFound(pos);
	                //printf("%i \n", areTherePlainCuts);
	                if (areTherePlainCuts==1){
	                    //printf("Planos de corte hallados \n");
	                    long int posAfterReading;
	                    arrayPlaneCuts = readPlaneCuts(pos, &numberPlaneCuts, &posAfterReading);
	                    fseek(file, posAfterReading, SEEK_SET);

	                }
	                pos = ftell(file);
	                int areThereTextures = texturesFound(pos);
	                if (areThereTextures==1){
	                    //printf("Texturas hallados \n");
	                    long int posAfterReading;
	                    arrayTextures = readTextures(currentTypeObjectReading,pos, &numberTextures, &posAfterReading);
	                    fseek(file, posAfterReading, SEEK_SET);

	                }
	                //printf("currentTypeObjectReading: %i \n", currentTypeObjectReading);
	                createObjectFromData(valuesRead, currentTypeObjectReading, indexValuesRead, arrayPlaneCuts, arrayTextures, arrayDraftPlanes, numberPlaneCuts, numberTextures, numberDraftPlanes);
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
	            }else if (strstr(temporalBuffer, "Disc_Object")!=NULL){
	                numberObjects++;
	                continue;
	            }else if (strstr(temporalBuffer, "Elipse_Object")!=NULL){
	                numberObjects++;
	                continue;
	            }else if (strstr(temporalBuffer, "Quadratic_Object")!=NULL){
	                numberObjects++;
	                continue;
	            }

	        }
	    }
	    fclose(file);
	    Objects = malloc(sizeof(struct Object)*numberObjects);
	    Lights= malloc(sizeof(struct Light)*numberLights);
	}
// ===============================================================

// Anti-aliasing =================================================
	struct Vector throwRay(long double x, long double y){
		struct Vector direction;
		long double Xw, Yw;
		
		Xw = (long double) ((x) * Xdif)/Hres + Xmin;
		Yw = (long double) ((y) * Ydif)/Vres + Ymin;

		direction.x = Xw - eye.x;
		direction.y = Yw - eye.y;
		direction.z = -eye.z;

		direction = normalize(direction);

		V.x = -direction.x;
		V.y = -direction.y;
		V.z = -direction.z;

		return direction;
	}

	struct Color ponderAAcolors(struct Color c1, struct Color c2, struct Color c3, struct Color c4){
		struct Color color;

		color.r = (c1.r + c2.r + c3.r + c4.r) / 4;
		color.g = (c1.g + c2.g + c3.g + c4.g) / 4;
		color.b = (c1.b + c2.b + c3.b + c4.b) / 4;

		return color;
	}

	long double getDistanceColors(struct Color c1, struct Color c2){

		return sqrt(pow(c1.r - c2.r,2) + pow(c1.g - c2.g,2) + pow(c1.b - c2.b,2));
	}

	int theyOK(struct Color c1, struct Color c2, struct Color c3, struct Color c4){
		long double dist, dist1, dist2, dist3, dist4, dist5, dist6;
		dist1 = getDistanceColors(c1, c2);
		dist2 = getDistanceColors(c1, c3);
		dist3 = getDistanceColors(c1, c4);
		dist4 = getDistanceColors(c2, c3);
		dist5 = getDistanceColors(c2, c4);
		dist6 = getDistanceColors(c3, c4);
		dist = (dist1 + dist2 + dist3 + dist4 + dist5 + dist6) / 6;
		if(dist <= similar){
			return 1;
		}else{
			return 0;
		}
	}

	struct Color getAAColor(long double x, long double y, int aaLevel){
		//printf("Level: %i - x: %i - y: %i\n", aaLevel, x, y);
		int areOK;
		int level = aaLevel;
		struct Color color, color1, color2, color3, color4;
		
		struct Vector direction;

		long double sum = 1/pow(2, level);
		level++;
		long double nextSum =  1/pow(2, level);
		
		direction = throwRay(x, y);
		color1 = getColor(eye, direction, maxReflection);

		direction = throwRay(x, y + sum);
		color2 = getColor(eye, direction, maxReflection);

		direction = throwRay(x + sum, y);
		color3 = getColor(eye, direction, maxReflection);

		direction = throwRay(x + sum, y + sum);
		color4 = getColor(eye, direction, maxReflection);

		areOK = theyOK(color1, color2, color3, color4);
		
		if(areOK == 0 && level <= maxAA){
			color1 = getAAColor(x          , y          , level);
			color2 = getAAColor(x          , y + nextSum, level);
			color3 = getAAColor(x + nextSum, y          , level);
			color4 = getAAColor(x + nextSum, y + nextSum, level);
			color = ponderAAcolors(color1, color2, color3, color4);
		}else {
			color = ponderAAcolors(color1, color2, color3, color4);
		}
		
		return color;
	}
// ===============================================================

// OP main: ======================================================
	int main(int argc, char *arcgv[]){
		howManyObjectsLights();
		printf("Lights: %i \n", numberLights);
		printf("Objects: %i \n", numberObjects);
		getSceneObjects();
		int i, j;
		
		struct Color color;
		struct Vector direction;
		
		Xdif = Xmax - Xmin;
		Ydif = Ymax - Ymin;
		similar = sqrt(3)/32;

		printf("\nRay Tracing\n...\n...\n");
		for (i = 0; i < Vres; i++){	
			for (j = 0; j < Hres; j++){
				//printf("[%i, %i]\n", j+1, i+1);
				color = getAAColor((long double) j, (long double) i, 0);
				Framebuffer[i][j] = color;
			}
		}

		printf("Rays: %LF\n", rays);
		saveFile();
		free(Objects);
		free(Lights);
		free(Framebuffer);
		printf("\nDONE.\n");
	}
// ===============================================================
