#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 *##### STRUCTURES #####
 */
typedef struct TypePaireExplicative* TypeListePairesExplicatives;

typedef struct TypePaireExplicative {
	int numHaplo1;
	int numHaplo2;
	TypeListePairesExplicatives suiv;
} TypePaireExplicative;

typedef struct TypeIndividu* TypeFichier;

typedef struct TypeIndividu {
    char* nomIndividu;
    int* genotype;
    TypeListePairesExplicatives tetePairesExplicatives;
    TypeFichier suiv;
} TypeIndividu;

/*
 *##### FONCTIONS #####
 */

void recupererParametres(int nbArgument, char* listeArguments[], char** adr_fichier, char** adr_nbIndividu, char** adr_tailleGeno) {

    //DEBUT
    if (nbArgument != 4) {
        printf("Au lancement du programme, vous devez rentrer 3 paramètres : \n");
        printf("- Le nom du fichier contenant les génotypes des individus\n");
        printf("- Le nombre d'individu étudié\n");
        printf("- La longueur du génotype étudié\n");
        printf("Merci\n");
        exit(1);
    } else {
        *adr_fichier = listeArguments[1];
        *adr_nbIndividu = listeArguments[2];
        *adr_tailleGeno = listeArguments[3];
    }
    //FIN
}


void recupererDonneesFichier(TypeFichier tete, int tailleGeno, int nbIndividu, char* fichier) {

    // Variables locales
    int i = 0;
    int j = 0;
    int estNomIndividu = 1;
    FILE* fp = NULL;
    char c[1];
    TypeFichier ptr;

    //DEBUT
    // Ouverture et lecture du fichier
    fp = fopen(fichier, "r");
 
    if (fp == NULL) {
        printf("Le fichier %s n'a pas pu être ouvert\n", fichier);
        exit(1);
    } else {
        printf("Le fichier %s existe\n", fichier);

        ptr = tete;
        ptr -> nomIndividu = (char*)malloc(sizeof(char)*20);
        ptr -> genotype = (int*)malloc(tailleGeno*sizeof(int));
        ptr -> suiv = NULL;

        do {
            c[0] = getc(fp);
            if (c[0] == '\n') {
            	if (i != (nbIndividu - 1)) {
	                ptr -> suiv = (TypeFichier)malloc(sizeof(TypeIndividu));
	                ptr = ptr -> suiv;
	                ptr -> nomIndividu = (char*)malloc(sizeof(char)*20);
	                ptr -> genotype = (int*)malloc(tailleGeno*sizeof(int));
	                ptr -> suiv = NULL;
	                i++;
	                estNomIndividu = 1;
                }
            } else if (c[0] == ' ') {
                estNomIndividu = 0;
                j = 0;
            } else {
                if (estNomIndividu) {   
                    ptr -> nomIndividu = strcat(ptr -> nomIndividu, c);
                } else {
                    ptr -> genotype[j] = atoi(c);
                    j++;
                }
            }

        } while (c[0] != EOF);
    }
    //FIN
}


void afficherListe(TypeFichier tete, int tailleGeno) {

    // Variables locales
    TypeFichier p;
    int i;

    //DEBUT
    printf("\n*************\n");
    printf("Affichage de la liste contenant les génotypes de chaque individu\n");
    p = tete;
    while (p != NULL) {
        printf("nom : %s\ngenotype : ", p->nomIndividu);
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", p-> genotype[i]);
        }
        printf("\n*************\n");
        p=p->suiv;
    }
    //FIN
}


void predictionEspaceRestraintHaplotypes(TypeFichier tete, int tailleGeno, int nbIndividu) {

    // Variables locales
    int i = 0;
    int k = 0;
    int j = 0;
    int t = 0;
    int s = 0;
    int pas = 0;
    int multiplicateur = 0;
    int nbDeux = 0;
    int nbHaploExplicatif = 0;
    TypeFichier p;
    int** tabTemp;

    //DEBUT
    // Allocation mémoire initiale pour le tableau temporaire qui contiendra les paires d'haplotypes possibles permettant d'expliquer le génotype d'un individu
    tabTemp = (int**)malloc(sizeof(int*)*tailleGeno);
        for (i=0 ; i<tailleGeno ; i++) {
            tabTemp[i] = (int*)malloc(sizeof(int));
        }
    // Parcours de la liste pour récupérer les génotypes et prédire l'espace restaint des haplotypes
    p = tete;
    while (p != NULL) { // Pour chaque individu
        // On compte le nombre de 2 dans le génotype
        nbDeux = 0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (p -> genotype[i] == 2) {
                nbDeux++;
            }
        }

        if (nbDeux == 0) {
        	nbDeux = 1;
        }
        nbHaploExplicatif = pow(2, nbDeux);

        // Reallocation mémoire 
        for (i=0 ; i<tailleGeno ; i++) {
            tabTemp[i] = (int*)realloc(tabTemp[i], sizeof(int)*nbHaploExplicatif);
        }
        // Génération des paires d'haplotypes expliquant le génotype
        k = 0; // k : nombre de deux rencontrés jusqu'à maintenant
        for (i=0 ; i<tailleGeno ; i++) { // Parcours du génotype d'un individu
        	if (p -> genotype[i] == 1) {
        		for (j=0 ; j<nbHaploExplicatif ; j++) {
        			tabTemp[i][j] = 1;
        		}
        	} else if (p -> genotype[i] == 0) {
        		for (j=0 ; j<nbHaploExplicatif ; j++) {
        			tabTemp[i][j] = 0;
        		}
        	} else {
        		k++; // On rencontre un 2 dans la séquence de génotype
        		// On calcul le pas et le multiplicateur
        		pas = ((pow(2, nbDeux) / pow(2, k)));
        		multiplicateur = (pow(2, nbDeux) / pow(2, k-1));

        		for (t=0 ; t<(pow(2, k-1)) ; t++) {
        			for (s=(t*multiplicateur) ; s<(t*multiplicateur+pas) ; s++) {
        				tabTemp[i][s] = 1;
        			}
        			for (s=(t*multiplicateur+pas) ; s<(t*multiplicateur+2*pas) ; s++) {
        				tabTemp[i][s] = 0;
        			}
        		}
        	}
        }

        // Affichage des haplotypes générés
        for (i=0 ; i<nbHaploExplicatif ; i++) { // Parcours du génotype d'un individu
        	printf("\n%d : \n", i);
        	for (j=0 ; j<tailleGeno ; j++) {
        			printf("%d", tabTemp[j][i]);
        	}
        	printf("\n");
        }
        printf("############################################\n");

        p=p->suiv;
    }

    // Libération de la mémoire
    for (i=0 ; i<tailleGeno ; i++) {
        free(tabTemp[i]);
    }
    free(tabTemp);


    

    //FIN
}


/*
 *##### MAIN #####
 */

int main(int argc, char* argv[]) {

    //VARIABLE
    char* fichier;
    char* nbIndividu;
    char* tailleGeno;
    TypeFichier ptrFichier;

    //DEBUT

    // Récupération des paramètres
    recupererParametres(argc, argv, &fichier, &nbIndividu, &tailleGeno);

    // Récupération des données présentes dans le fichier
    // On alloue de la mémoire à la première structure de la liste puis on vérifie si l'allocation à fonctionnée
    ptrFichier = (TypeFichier)malloc(sizeof(TypeIndividu));

    if(ptrFichier==NULL){
        printf("L'allocation mémoire a échoué\n");
        exit(1);
    }

    recupererDonneesFichier(ptrFichier, atoi(tailleGeno), atoi(nbIndividu), fichier);
    
    // Afficher les données présentes dans le fichier
    afficherListe(ptrFichier, atoi(tailleGeno));

    // Prédiction de l'espace restraint des haplotypes expliquant les génotypes
    predictionEspaceRestraintHaplotypes(ptrFichier, atoi(tailleGeno), atoi(nbIndividu));

    free(ptrFichier);
    return 0;
    // FIN
}