
#ifndef TDA_LDE_H
#define	TDA_LDE_H

#include <stdio.h>
#include <stdlib.h>
#include "structures.h"
#include "malloc.h"
 //LDE simple objetos
/*El tipo de dato se cambiara dependiendo de lo que se necesita*/




typedef  struct Object *Elemento;

 typedef  struct Nodo
{
    Elemento dato;
    struct Nodo *siguiente;
}NODO;

typedef  struct Light *Elemento1;

typedef  struct NodoLuz
{
    Elemento1 dato;
    struct NodoLuz *siguiente;
}NODOLUZ;


extern NODO *Objects;
extern NODOLUZ *Lights;


NODO* crearNodo(Elemento );

void insertPrimero(NODO** , Elemento );

void append(NODO** , Elemento );

Elemento get(NODO**, int);

int getSize(NODO** );

void deleteListObjects(NODO**);
//LDE simple luces


NODOLUZ* crearNodoLuz(Elemento1 );

void insertPrimeroLuz(NODOLUZ**, Elemento1 );

void appendLuz(NODOLUZ**, Elemento1 );

struct Light* getLuzAt(NODOLUZ** , int );

int getSizeLuzCollection(NODOLUZ**);

void deleteListLights(NODOLUZ**);




#endif  /* TDA_LDE_H */