#include "systeme_lineaire.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

void afficheD2(int l, int c, double (*A)[c], double *B) {
	for (int i=0; i<l; i++){
		for (int j=0; j<c; j++) {
			printf("%f ", A[i][j]);
		}
		printf("| %f\n",B[i]);
	}
}

int systeme(int l, int c, double (*A)[c], double *B){
	// résout un système linéaire où les coeff diagonaux sont non nuls
	// A, B tableau des coefficients et résultats du systeme, l/c nbr de lignes/colonnes du système
	short min = (l+c-abs(l-c))/2;
	short hasZeroPivot = 0;
	for (int i=0; i<min; i++) {
		if (A[i][i] == 0){
			hasZeroPivot++;
		}
	}
	// s'il y a autant de pivot nul que d'équation, on échange la première avec la dernière ligne du système
	if (hasZeroPivot == l) {
		double tmpA;
		for (int i=0; i<c; i++) {
			tmpA = A[0][i];
			A[0][i] = A[l-1][i];
			A[l-1][i] = tmpA;
		}
		double tmpB = B[0];
		B[0] = B[l-1];
		B[l-1] = tmpB;
	}
	//afficheD2(l, c, A, B);
	for (int k=0; k<=hasZeroPivot; k++){
		for (int i=0; i<c; i++) {
        		double p = A[i][i]; // "pivot" ou plutôt coefficients de la diagonale principale
			if (p != 0) {
				for (int j=0; j<l; j++) {
					double n = A[j][i]; // numérateur du coefficient qui multiplie la ligne j
					if (i != j) {
						for (int k=0; k<c; k++) {
							A[j][k] = A[j][k] - (A[i][k] * n / p);
						}		
						B[j] = B[j] - (B[i] * n / p);
					}
				}
			}
		
			//printf("\nEtape %d.%d:\n",k,i);
                	//afficheD2(l, c, A, B);
		}
	}
	for (int i=0; i<l; i++) {
		if (A[i][i] != 0) {
                	B[i] = B[i]/A[i][i];
                        A[i][i] = A[i][i]/A[i][i];
		}
	}
	//printf("\n");
	//afficheD2(l, c, A, B);
	bool isZeroCoefficientRow = true; // est vrai si une ligne du système n'a que des coefficients nuls.
	short nbEmptyRow = 0; // nombre de ligne du systeme qui n'a que des zero 
	for (int i=0; i<l; i++) {
		isZeroCoefficientRow = true; // la nouvelle ligne du tableau est considérée remplie de coefficients nuls au départ
		for (int j=0; j<c; j++) {
			if (A[i][j] != 0) isZeroCoefficientRow = false; // si un coefficient de la ligne est different de zero, changer la valeur du booleen à false
		}
		
		// Si une ligne est entièrement nul incrémenter le compteur de ligne nulles
		if (isZeroCoefficientRow && B[i] == 0) {
			nbEmptyRow++;
		}
		
		// Si une ligne n'a que des coeff nuls mais que le résultat n'est pas nul -> système incompatible (retourner -1)
		else if (isZeroCoefficientRow && B[i] != 0) {
			return -1;  
		}
	}	
	// si le système est compatible, retourner le nombre de variable libre
	return c - (l - nbEmptyRow); 
}
