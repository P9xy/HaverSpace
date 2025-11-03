"""
    Outil de GeoINT pour :
    déterminer les coordonnées latitudes, longitudes d'un endroit inconnu par trilatération
 ou calculer une distance géodésique précise avec les formules d'haversine

    Auteur {P9xy} github - https://github.com/P9xy/
    Licensed under MIT License 
"""

from math import *
import re
import traceback
import numpy as np

class BadCoordinates(Exception):
    def __init__(self, *args):
        if args:
        # si args existe
            if isinstance(args[0], dict):
            #verifie que args[0] est un dictionnaire
                dic = args[0]
                self.message = ""
                for key, val in dic.items():
                    if val[1] == False:
                    # Si la validité de la valeur de la clé actuelle key du tableau est mauvaise : ajouter key=valeur au message
                        self.message += f" '{key}'={val[0]}," 
    
                if len(self.message) > 0:
                    self.message = "Attention, les attributs {" + self.message + " } de " + args[1] +" sont incorrects."
            else:
                raise TypeError(f"BadCoordinates: error, '{args[0]}'s type should be dict")
        else:
            self.message = None
    def __str__(self):
        if self.message is not None:
            return f"BadCoordinates: {self.message}"
        return f"BadCoordinates: error bad coordinates"

class Geolocalisation:
    # Permet de déterminer la localisation d'un objet grace à au moins trois lieux connus par trilatération
    # Un objet représente un 'lieu',caractérisé par sa latitude, longitude et distance éventuelle à un autre objet
    def __init__(self, place="", latitude=None, longitude=None, altitude=0, distance=None):
        self.nom = place
        self.lat = latitude
        self.longi = longitude
        self.dist = distance
        self.h = altitude
        if not (self.lat is None or self.longi is None):
            self.vlat = -90 <= self.lat <= 90 # vérifie si la latitude est valide
            self.vlon = -180 <= self.longi <= 180
        else:
            self.vlat, self.vlon = True, True
        if self.dist is not None:
            self.vdist = self.dist > 0
        else:
            self.vdist = True
        self.lst = {'latitude':(self.lat,self.vlat),'longitude':(self.longi,self.vlon),'distance':(self.dist,self.vdist)} # dictionnaire des attributs de Geolocalisation et de leur valeur et validité
        if not (self.vdist & self.vlat & self.vlon):
            raise BadCoordinates(self.lst,self.nom)

    def __str__(self):
        self.verti = 'N' if self.lat>=0 else 'S'
        self.hori =  'E' if self.longi>=0 else 'O'
        return f"{self.nom} est à la position :\nlatitude :'{abs(self.lat)}{self.verti}',\nlongitude :'{abs(self.longi)}{self.hori}'"
    
    def hav(self, teta):
        # retourne le nb d'haversine de teta
        return (sin(teta/2)**2)
    def distance(self,self2, comm=True, km=False, precision=5):
        # calcule la distance minimale entre deux objets par la formule d'haversine
        R = 6378137 # rayon de la terre
        D = 0 # distance entre 1 et 2
        deltaphi = self2.lat-self.lat # difference de latitude entre 1 et 2
        deltalambda = self2.longi-self.longi # difference de longitude entre 1 et 2
        teta1 = (deltaphi * 2 * pi) / 360 # conversion en radians
        teta2 = (deltalambda * 2 * pi) / 360 
        phi1, phi2 = (self.lat * 2 * pi) / 360, (self2.lat * 2 * pi) / 360 
        h = self.hav(teta1) + cos(phi1) * cos(phi2) * self.hav(teta2)
        D = 2 * R  * asin(sqrt(h))
        self.distance = D
        self2.distance = D 
        if comm:
            print(f"\nLa distance entre {self.nom} et {self2.nom} est égale à approximativement {round(D*(10**-3),precision) if km else round(D,precision)} {"km" if km else "mètres"}.\n")
        return D
    
    def sphere2cartesian(self):
        #converti les coordonées spheriques en cartésiennes 
        # pour utiliser dans l'equation de cercle
        R = 6378137 #grand axe en m
        N = R/(sqrt(1-(exp(2)*(sin(self.lat))**2)))
        x = (N+self.h)*cos(self.lat)*cos(self.longi)
        y = (N+self.h)*cos(self.lat)*sin(self.longi)
        z = (N*(1-exp(2))+self.h)*sin(self.lat)
        return (x,y,z)
    
    def trilateration(self,self1,self2,self3,comm=True):
        print(f"Calcul de la position de {self.nom} ...\n")
        # changement du spherique au cartésien pour avoir des equations de cercle
        x1,y1,z1= self1.sphere2cartesian()
        x2,y2,z2= self2.sphere2cartesian()
        x3,y3,z3= self3.sphere2cartesian()
        r1,r2,r3 = self1.dist, self2.dist, self3.dist
       
        # On transforme le système d'équations linéaires en une équation matricielle AX=B
        A = np.array([[2*x2-2*x1,2*y2-2*y1,2*z2-2*z1],
                      [2*x3-2*x2,2*y3-2*y2,2*z3-2*z2],
                      [2*x1-2*x3,2*y1-2*y3,2*z1-2*z3]])
       
        B = np.array([[(r1**2)-(r2**2)+(x2**2)-(x1**2)+(y2**2)-(y1**2)+(z2**2)-(z1**2)],
                     [(r2**2)-(r3**2)+(x3**2)-(x2**2)+(y3**2)-(y2**2)+(z3**2)-(z2**2)],
                     [(r3**2)-(r1**2)+(x1**2)-(x3**2)+(y1**2)-(y3**2)+(z1**2)-(z3**2)]])
        
        resultat = np.linalg.lstsq(A, B, rcond=None)
        X = resultat[0].ravel()
        residuals = resultat[1]
        try:
            self.lat, self.longi=self.cart2sphere(X)
        except Exception as e:
            raise Exception("Erreur dans l'execution de la conversion des coordonnées cartésiennes en sphériques") from e
        if comm:
            print("Résultat (x, y, z) ≈", X.ravel())
            if residuals.size > 0:
                print("Erreur résiduelle :", residuals)
            else:
                print("Aucune erreur résiduelle.")
        
        return X

def main():
    print(r"#################################################")
    print(r"#                                               #")
    print(r"#       ___           ___           ___         #")
    print(r"#      /\__\         /\  \         /\__\        #")
    print(r"#     /:/  /        /::\  \       /:/  /        #")
    print(r"#    /:/__/        /:/\:\  \     /:/  /         #")
    print(r"#   /::\  \ ___   /::\~\:\  \   /:/__/  ___     #")
    print(r"#  /:/\:\  /\__\ /:/\:\ \:\__\  |:|  | /\__\    #")
    print(r"#  \/__\:\/:/  / \/__\:\/:/  /  |:|  |/:/  /    #")
    print(r"#       \::/  /       \::/  /   |:|__/:/  /     #")
    print(r"#       /:/  /        /:/  /     \::::/__/      #")
    print(r"#      /:/  /        /:/  /       ~~~~          #")
    print(r"#      \/__/         \/__/                      #")
    print(r"#                                               #")
    print(r"#################################################")

if __name__ == '__main__':
    optInvalide = True # continue tant que l'utilisateur ne le dit pas explicitement ou que son choix est incorrect 
    main()
    
    while optInvalide:
        print("Choisissez une option:\n\t(1) Calcul de distance minimale géodésique entre deux points.\n\t(2) Calcul de la position d'un point par trilatération.\n\t(Q) Quitter.\n")
        choix = input("Votre choix: ")
        if choix == '1':
            print("Veuillez rentrer le 'nom' la latitude et longitude des 2 endroits. (Ex: Paris, 48.8, 2.35):\n")
            lieuA, latA, lonA = re.split(r'[, ]+', input("Lieu 1 :> "))
            lieuB, latB, lonB = re.split(r'[, ]+', input("Lieu 2 :> "))
            try:
                latA, lonA = float(latA), float(lonA)
                latB, lonB = float(latB), float(lonB)
                P1 = Geolocalisation(lieuA, latA, lonA)
                P2 = Geolocalisation(lieuB, latB, lonB)
                mu = input("Déterminez la précision du résultat à 10^-n, laisser vide par défaut : ")
                if len(mu) != 0 and mu.isdigit():
                        mu = int(mu)
                        P1.distance(P2,km=True,precision=mu)
                else:
                    P1.distance(P2,km=True)
                print("-------------------------------------------------")
            except ValueError as val:
                print(f"'{val}")
        elif choix == '2':
            try:
                print("Veuillez rentrer les informations 'nom', latitude, longitude, altitude et distance en mètre des lieux connus. (Ex: Paris, 48.8, 2.35, 35, 80000):\n")
                lieuA, latA, lonA, h1, dA = re.split(r'[, ]+', input("Lieu 1 :> "))
                lieuB, latB, lonB, h2, dB = re.split(r'[, ]+', input("Lieu 2 :> "))
                lieuC, latC, lonC, h3, dC = re.split(r'[, ]+', input("Lieu 3 :> "))
                latA, lonA, h1, dA = float(latA), float(lonA), float(h1), float(dA)
                latB, lonB, h2, dB = float(latB), float(lonB), float(h2), float(dB)
                latC, lonC, h3, dC = float(latC), float(lonC), float(h3), float(dC)
                P1 = Geolocalisation(lieuA, latA, lonA, h1, dA)
                P2 = Geolocalisation(lieuB, latB, lonB, h2, dB)
                P3 = Geolocalisation(lieuC, latC, lonC, h3, dC)
                x = Geolocalisation("lieu inconnu")
                print(f"Les coordonnées du point {x.nom} sont {x.trilateration(P1,P2,P3)}.\n")
                print("-------------------------------------------------")
            except ValueError as v:
                print(f"{ValueError} Les coordonnées sont des nombres réels")
            except BadCoordinates as bc:
                print(bc)
            except Exception as e:
                print("Une erreur est survenue ",e)
                print("---detail---")
                traceback.print_exc()
        elif choix.upper() == 'Q':
            print("Au-revoir")
            optInvalide = False
        else:
            print("Mauvaise option.\n")
            optInvalide = False