

#include <stdio.h>
#include <stdlib.h>
#include "malloc.h"
#include "structures.h"
#include "TDA_LDE.h"

 


 NODO *Objects = NULL;
 NODOLUZ *Lights= NULL;

NODO* crearNodo(Elemento valor)
{
    NODO *aux;
    aux = (NODO*)malloc(sizeof(NODO));
    aux -> dato = valor;
    aux -> siguiente = NULL;
    return aux;
}

void insertPrimero(NODO** cabeza, Elemento entrada)
{
    NODO *temp;
    temp = crearNodo(entrada);
    temp -> siguiente = *cabeza;
    *cabeza = temp;
}

void append(NODO** cabeza, Elemento entrada)
{
    NODO *ultimo;
    ultimo = *cabeza;
    
    if (ultimo == NULL){
        *cabeza = crearNodo(entrada);
    }else{
        for (; ultimo -> siguiente; )
            
            ultimo = ultimo -> siguiente;
            ultimo -> siguiente = crearNodo(entrada);
    }
}

struct Object* get(NODO** cabeza, int i){
    NODO *ultimo;
    ultimo = *cabeza;
    int counter = 0;
    
    if (i < 0 || i >= getSize(cabeza)){
        printf("IllegalArgument. Índice %i ingresado.", i);
        return NULL;
    }
    
    if (ultimo == NULL){
        return NULL;
    }else{
        //Hay un elemento
        if (i == counter){
            return ultimo->dato; //Devuelve el primer valor
        }
        counter++; 
        for (; ultimo -> siguiente; )
            ultimo = ultimo -> siguiente;
            counter++;
            if (i == counter){
                return ultimo -> dato; //Devuelve el primer valor
            }
    }

}

int getSize(NODO** cabeza){
    int size = 0;
    NODO *ultimo;
    ultimo = *cabeza;
    
    if (ultimo == NULL){
        return size; //No hay elementos
    }else{
        size++; //Hay un elemento
        for (; ultimo -> siguiente; )
            
            size++;
    }
    return size;
}

void deleteListObjects(NODO** cabeza) {
    NODO *temp;
    temp = *cabeza;
    
    if (temp == NULL){
        return;
    }else{
        *cabeza = *cabeza -> siguiente;
        free(temp);
        for (;*cabeza -> siguiente; )
            
            temp = *cabeza;
            *cabeza = *cabeza -> siguiente;
            free(temp);
    }
}


//PARA LUCES


NODOLUZ* crearNodoLuz(Elemento1 valor)
{
    NODOLUZ *aux;
    aux = (NODOLUZ*)malloc(sizeof(NODOLUZ));
    aux -> dato = valor;
    aux -> siguiente = NULL;
    printf("Dentro de crearNodoLuz: %LF \n", valor->Ip);
    return aux;
}

void insertPrimeroLuz(NODOLUZ** cabeza, Elemento1 entrada)
{
    printf("Dentro de insertPrimeroLuz: %LF \n", entrada->Ip);
    NODOLUZ *temp;
    temp = crearNodoLuz(entrada);
    temp -> siguiente = *cabeza;
    *cabeza = temp;
    printf("Dentro de insertPrimeroLuz: %LF \n", temp->dato->Ip);
    return;
}

void appendLuz(NODOLUZ** cabeza, Elemento1 entrada)
{
    NODOLUZ *ultimo;
    ultimo = *cabeza;
    
    if (ultimo == NULL){
        *cabeza = crearNodoLuz(entrada);
    }else{
        for (; ultimo -> siguiente; )
            
            ultimo = ultimo -> siguiente;
            ultimo -> siguiente = crearNodoLuz(entrada);
    }
}

struct Light* getLuzAt(NODOLUZ** cabeza, int i){
    NODOLUZ *ultimo;
    ultimo = *cabeza;
    int counter = 0;
    
    if (i < 0 || i >= getSizeLuzCollection(cabeza)){
        printf("IllegalArgument. Índice %i ingresado.", i);
        return NULL;
    }
    
    if (ultimo == NULL){
        return NULL;
    }else{
        //Hay un elemento
        if (i == counter){
            return ultimo->dato; //Devuelve el primer valor
        }
        counter++; 
        for (; ultimo -> siguiente; )
            ultimo = ultimo -> siguiente;
            counter++;
            if (i == counter){
                return ultimo -> dato; //Devuelve el primer valor
            }
    }

}

int getSizeLuzCollection(NODOLUZ** cabeza){
    int size = 0;
    NODOLUZ *ultimo;
    ultimo = *cabeza;
    
    if (ultimo == NULL){
        return size; //No hay elementos
    }else{
        size++; //Hay un elemento
        for (; ultimo -> siguiente; )
            
            size++;
    }
    return size;
}




void deleteListLights(NODOLUZ** cabeza) {
    NODOLUZ *temp;
    temp = *cabeza;
    
    if (temp == NULL){
        return;
    }else{
        *cabeza = *cabeza -> siguiente;
        free(temp);
        for (;*cabeza -> siguiente; )
            
            temp = *cabeza;
            *cabeza = *cabeza -> siguiente;
            free(temp);
    }
}