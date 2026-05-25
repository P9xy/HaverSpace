/*
    Outil de GeoINT pour :
    déterminer les coordonnées latitudes, longitudes d'un endroit inconnu par trilatération
 ou calculer une distance géodésique précise avec les formules d'haversine

    Auteur {P9xy} github -> https://github.com/P9xy/
    Licensed under MIT License 
*/
 
#include "systeme_lineaire.h"
#include <stdio.h>
#include <math.h>

const float PI = 3.141592653589793; // valeur de pi
const float R = 6371000; // valeur moyenne du rayon terrestre en mètres donné par wgs-84

struct Point {
	float longitude;	// longitude d'un point sur terre en degré
	float lattitude;	// lattitude du point en degré
};

struct Params {
	double *c11,*c12,*c13,*c14,*c21,*c22,*c23,*c24,*c31,*c32,*c33,*c34,*r1,*r2,*r3; // pointeurs vers les coefficients et résultats du système de trilatération
	struct Point x; // structure du premier des 3 points connus du problème  
	struct Point y; 
	struct Point z;
	float d1, d2, d3; // distances respectives en mètres de x,y,z au point cherché
};

float haversine(float teta) {
	// Associe au nombre réel teta son nombre d'haversine. X->hav(X) / hav(teta)=sin²(teta/2)
	float hav = ((1-cos(teta))/2);
	return hav;
}

float distanceGeo(struct Point p1, struct Point p2, float r) {
	// r: rayon de la terre considéré en mètre 
	// calcule la distance géodésique en mètres entre deux points sur terre.
	float latDelta = p2.lattitude - p1.lattitude; 	// delta entre les lattitudes des deux points
	float longiDelta = p2.longitude - p1.longitude;	// delta entre les longitudes des deux points
	float distance = 2*r*asin(sqrt(haversine(latDelta) + (cos(p2.lattitude)*cos(p1.lattitude)*haversine(longiDelta))));
	return distance;
}

void menu() {
	for (int i=0; i<100; i++)printf("_");
	printf("\n\n|||  |||     ||||||   \\\\\\           ///\n");
	printf("|||  |||    |||  |||   \\\\\\         ///\n");
	printf("||||||||   |||====|||   \\\\\\       ///\n");
	printf("||||||||   |||====|||	 \\\\\\     ///\n");
	printf("|||  |||  |||      |||	  \\\\\\   ///\n");
	printf("|||  |||  |||      |||     \\\\\\_///\n\n");

	printf("Quelle opération souhaitez-vous effectuer ?\n\t- 0 - Afficher l'aide contextuel.\n\t- 1 - Calcul de distance géodésique entre deux points.\n\t- 2 - 'Trilatération'.\n\t- 3 - Quitter.\n\n");
	for (int i=0; i<100; i++)printf("_");
	printf("\n");
}

void help() {
	printf("\t(1): Permet de calculer la distance (minimale) géodésique entre deux points sur Terre dans le système métrique. (Distance 'à vol d'oiseau.').\n\t(2): Permet de déterminer la position définie par latitude et longitude (ex Paris: 48.85, 2.35), d'un point sachant la position respective de trois autres points et leur distance au point recherché.\n");
}

void modifierCoord(struct Point *p) {
	// demande à l'utilisateur d'entrer les nouvelles coordonnées d'un point donné en paramètre
	float lat=0.0, longi=0.0;
	printf("Veuillez rentrer dans la console les coordonnées du point en degré :\n\tlattitude :>>");
	scanf("%f",&lat);
	printf("\tlongitude :>>");
	scanf("%f",&longi);
	p->lattitude = lat*2*PI/360; // conversion en radian de la lattitude
	p->longitude = longi*2*PI/360; // même chose pour la longitude

	printf("\n");
	for (int i=0; i<100; i++)printf("_");
	printf("\n\t\t\tOpération effectuée !\nNouvelles coordonnées : \n\tLattitude: %.4f rad\n\tLongitude: %.4f rad\n", p->lattitude, p->longitude);
	for (int i=0; i<100; i++)printf("_");
	printf("\n");
	
}

void initCoef(struct Params *p){
	// initialise les coefficients c11<->c34 des inconnues de la matrice et les constantes du résultat
	// { c11 * X + c12*c13 * Y + c12*c14 * Z = r1
	// { c21 * X + c22*c23 * Y + c22*c24 * Z = r2
	// { c31 * X + c32*c33 * Y + c32*c34 * Z = r3
	// R = constante (rayon Terrestre), d1-2-3 = valeur respective de la distance des points 1-2-3 au point recherché
	*p->c11 = sin(p->x.lattitude); 
	*p->c12 = cos(p->x.lattitude); 
	*p->c13 = cos(p->x.longitude); 
	*p->c14 = sin(p->x.longitude);
	*p->c21 = sin(p->y.lattitude); 
	*p->c22 = cos(p->y.lattitude); 
	*p->c23 = cos(p->y.longitude); 
	*p->c24 = sin(p->y.longitude);
	*p->c31 = sin(p->z.lattitude); 
	*p->c32 = cos(p->z.lattitude); 
	*p->c33 = cos(p->z.longitude); 
	*p->c34 = sin(p->z.longitude);
	*p->r1 = cos(p->d1/R); 
	*p->r2 = cos(p->d2/R); 
	*p->r3 = cos(p->d3/R);
	//printf("Système d'équation à résoudre :\n");
	//printf("\t{ %lf * X + %lf*%f * Y + %lf*%f * Z = %f\n",*p->c11,*p->c12,*p->c13,*p->c12,*p->c14,*p->r1);
	//printf("\t{ %lf * X + %lf*%f * Y + %lf*%f * Z = %f\n",*p->c21,*p->c22,*p->c23,*p->c22,*p->c24,*p->r2);
	//printf("\t{ %lf * X + %lf*%f * Y + %lf*%f * Z = %f\n",*p->c31,*p->c32,*p->c33,*p->c32,*p->c34,*p->r3);
}

void initDistance(float *d1, float *d2, float *d3){
	// Procédure pour permettre à l'utilisateur du programme de changer les distances d1,d2,d3
	printf("Veuillez donner les distances respectives en kilomètres des points 1,2 puis 3 au point recherché (Séparées par un espace !):\n\tEx: 130 35.2 81\n");
	printf("\t:>>");
	scanf("%f %f %f",d1,d2,d3);
	*d1 = *d1*1000;
	*d2 = *d2*1000;
	*d3 = *d3*1000;
}

int main(int argc, char *argv[]) {
	char choix = '0';
	do {
		menu();
		printf(":>>");
	        scanf(" %c", &choix);
		while (getchar() != '\n');
		printf("\n");
		switch (choix) {
			case '0':
				help();
				break;
			case '1':
				struct Point p1 = {0.0000, 0.0000}; //initialiser le point à des coordonnées d'origine
				struct Point p2 = {0.0000, 0.0000};
				modifierCoord(&p1); // changer les coordonnées du point
				modifierCoord(&p2);
				printf("Lattitude p1 = %.4f\n", p1.lattitude);
				double d = distanceGeo(p1, p2, R)/1000;
				printf("La distance entre le point 1 (%.4f,%.4f), et le point 2 (%.4f,%.4f) est de %lf km.\n", p1.lattitude, p1.longitude, p2.lattitude, p2.longitude, d);
				break;
			case '2':
				printf("Fonctionnalité en cours de développement.\n");
				//Initialisation des trois points connus
				struct Point x = {0.0000, 0.0000};
				struct Point y = {0.0000, 0.0000};
				struct Point z = {0.0000, 0.0000};
				modifierCoord(&x); 
				modifierCoord(&y);
				modifierCoord(&z);
				double c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,r1,r2,r3;
				float d1=0,d2=0,d3=0;
				initDistance(&d1,&d2,&d3);
				struct Params params = {&c11,&c12,&c13,&c14,&c21,&c22,&c23,&c24,&c31,&c32,&c33,&c34,&r1,&r2,&r3,x,y,z,d1,d2,d3};
				initCoef(&params);
				double system1[3][3] = {{c11,c12*c13,c12*c14},{c21,c22*c23,c22*c24},{c31,c32*c33,c32*c34}};
				double resultat1[3] = {r1,r2,r3};

				int resolution = systeme(3,3,system1, resultat1); // résolution du système non linéaire en posant X=sin(lattitude)
				if (resolution == 0){
					double lattitude = asin(resultat1[0]); // lattitude recherchée en radians
					printf("\n\nCalcul des coordonnées :\n\tLATTITUDE déterminée : %lf°.\n",lattitude*360/(2*PI));
					double cosLat = cos(lattitude);
					double system2[2][2] = {{c12*c13*cosLat,c12*c14*cosLat},{c22*c23*cosLat,c22*c24*cosLat}};
					double s1 = c11*resultat1[0], s2 = c21*resultat1[0];
					double resultat2[2] = {r1-s1,r2-s2};
					int resolu2 = systeme(2,2,system2, resultat2);
					double longitude = asin(resultat2[1]); // longitude recerchée en radians
					printf("\tLONGITUDE déterminée : %lf°.\n",longitude*360/(2*PI));
				}
				else if (resolution == -1){
					for (int i=0; i<100; i++)printf("_");
					printf("\n\nAucunes coordonnées valides correspondantes ...\n");
				}
				else {
					for (int i=0; i<100; i++)printf("_");
					printf("\n\nIl semble y avoir des coordonnées, mais il n'a pas été possible de trouver celles exactes !\n");
				}
				break;
			case '3':
				printf("Au-revoir.\n");
				break;
			default:
				printf("Mauvaise saisie ...\n");
				break;
		}
	} while (choix != '3');
	return 0;
}
