Tecnológico de Costa Rica IC8019 - Introducción A Los Gráficos Por Computadora

Profesor: Dr. Francisco Torres

Estudiantes: Carlos Girón Alas Julián J. Méndez O. Daniel Troyo Garro

Programa que grafica objetos básicos matemáticos en una escena tridimensional utilizando la técnica del ray-tracing. El archivo con el código fuente es rt.c , Adicional a estos archivos se encuentran el archivo

makefile
escena1.txt


EL archivo de texto llamado escena1.txt es el que alberga los datos de la escena así como de los objetos presente en ella. EL formato para introducir los datos es el siguiente, permitiendo posicionar los objetos en cualquier orden, aunque el orden de los datos dentro de cada objeto debe ser el mismo y encontrarse de la siguiente manera:

identificador = valor flotante
tripleta = valor valor valor

EL formato es el siguiente, permite introducir comentarios añadiendo un # en la línea y líneas en blanco: 


Scene_Data *Comentario*
    #Datos de la escena
    Ia = IaNum
    Xmin = XminNum
    Ymin = YminNum
    Xmax = XmaxNum
    Ymax = YmaxNum
    Hres = HresNum
    Vres = VresNum
    Epsilon = EpsilonNum
    Eye = Xe Ye Ze
    Background = R, G, B

Light_Object *Comentario*
    Posicion = Xp, Yp, Zp
    c1 = valor
    c2 = valor
    c3 = valor
    Ip = valor


Sphere_Object *Comentario*
    Posicion = Xc, Yc, Zc
    Radio = valor
    Kd = valor
    Ka = valor
    Kn = valor
    Ks = valor
    Color = R, G, B

Polygon_Object *Comentario*
    Color = R, G, B
    Kd = valor
    Ka = valor
    Kn = valor
    Ks = valor
    Vertice = x, y, z
    END_Vertices

Cylinder_Object *Comentario*
    Ancla = Xo, Yo, Zo
    Q = Xq, Yq, Zq
    Radio = valor
    d1 = valor
    d2 = valor
    Kd = valor
    Ka = valor
    Kn = valor
    Ks = valor
    Color = R, G, B

Cone_Object *Comentario*
    Ancla = Xo, Yo, Zo
    Q = Xq, Yq, Zq
    k1 = valor
    k2 = valor
    d1 = valor
    d2 = valor
    Kd = valor
    Ka = valor
    Kn = valor
    Ks = valor
    Color = R, G, B