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
	double proba;
	TypeListePairesExplicatives suiv;
} TypePaireExplicative;

typedef struct TypeIndividu* TypeFichier;

typedef struct TypeIndividu {
    char* nomIndividu;
    int* genotype;
    TypeListePairesExplicatives tetePairesExplicatives;
    TypeFichier suiv;
} TypeIndividu;

typedef struct TypeGenoExplicatif* TypeListeGenosExplicatifs;

typedef struct TypeGenoExplicatif {
    char* nomIndividu;
    int numHaplo2;
    TypeListeGenosExplicatifs suiv;
} TypeGenoExplicatif;

typedef struct TypeHaplotype* TypeListeHaplotypes;

typedef struct TypeHaplotype {
    int numHaplo;
    int* haplotype;
    int cpt_haplo;
    double freq;
    TypeListeGenosExplicatifs teteGenosExplicatifs;
    TypeListeHaplotypes suiv;
} TypeHaplotype;

/*
 *##### FONCTIONS #####
 */

void recupererParametres(int nbArgument, char* listeArguments[], char** adr_fichier, char** adr_nbIndividu, char** adr_tailleGeno) {
		
		printf("############################################\n");
		printf("Début de la fonction recupererParametres\n");
		printf("############################################\n");
		
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

		printf("############################################\n");
		printf("Début de la fonction recupererDonneesFichier\n");
		printf("############################################\n");

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
        ptr -> tetePairesExplicatives = (TypePaireExplicative*)malloc(sizeof(TypePaireExplicative));
        ptr -> suiv = NULL;

        do {
            c[0] = getc(fp);
            if (c[0] == '\n') {
            	if (i != (nbIndividu - 1)) {
	                ptr -> suiv = (TypeFichier)malloc(sizeof(TypeIndividu));
	                ptr = ptr -> suiv;
	                ptr -> nomIndividu = (char*)malloc(sizeof(char)*20);
	                ptr -> genotype = (int*)malloc(tailleGeno*sizeof(int));
	                ptr -> tetePairesExplicatives = (TypePaireExplicative*)malloc(sizeof(TypePaireExplicative));
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

	printf("############################################\n");
	printf("Début de la fonction afficherListe\n");
	printf("############################################\n");

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


int predictionEspaceRestraintHaplotypes(TypeFichier tete, int tailleGeno, int nbIndividu, TypeListeHaplotypes ptrListeHaplo) {

	printf("############################################\n");
	printf("Début de la fonction predictionEspaceRestraintHaplotypes\n");
	printf("############################################\n");

    // Variables locales
    int i = 0;
    int k = 0;
    int j = 0;
    int t = 0;
    int s = 0;
    int l = 0;
    int g = 0;
    int nbHaploTotaux = 0;
    int pas = 0;
    int multiplicateur = 0;
    int nbDeux = 0;
    int nbHaploExplicatif = 0;
    TypeFichier pfic;
    int** tabTemp;
    TypeListeHaplotypes pHaplo;
    TypeListePairesExplicatives pPE;

    //DEBUT
    // Allocation mémoire initiale pour le tableau temporaire qui contiendra les paires d'haplotypes possibles permettant d'expliquer le génotype d'un individu
    tabTemp = (int**)malloc(sizeof(int*)*tailleGeno);
        for (i=0 ; i<tailleGeno ; i++) {
            tabTemp[i] = (int*)malloc(sizeof(int));
        }
    // Parcours de la liste pour récupérer les génotypes et prédire l'espace restaint des haplotypes
    pfic = tete;
    
    // Initialisation de la liste qui contiendra tous les haplotypes explicatifs
        pHaplo = ptrListeHaplo;
        pHaplo -> haplotype = (int*)malloc(sizeof(int)*tailleGeno);
        pHaplo -> freq = 0.0;
        pHaplo -> teteGenosExplicatifs = (TypeGenoExplicatif*)malloc(sizeof(TypeGenoExplicatif));
        pHaplo -> suiv = NULL;
    
    while (pfic != NULL) { // Pour chaque individu
        // On compte le nombre de 2 dans le génotype
        nbDeux = 0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (pfic -> genotype[i] == 2) {
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
        	if (pfic -> genotype[i] == 1) {
        		for (j=0 ; j<nbHaploExplicatif ; j++) {
        			tabTemp[i][j] = 1;
        		}
        	} else if (pfic -> genotype[i] == 0) {
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
        
        // On garde les haplotypes dans une liste chainée TypeListeHaplotype et les paires explicatives dans TypeListePairesExplicatives
        l = 0;
        g = nbHaploTotaux;
        l = (nbHaploTotaux + nbHaploExplicatif - 1);
        for (i=0 ; i<(nbHaploExplicatif) ; i++) {
        	pPE = pfic -> tetePairesExplicatives;
    		pPE -> suiv = NULL;
        	if (i < (nbHaploExplicatif/2)) {
        	    pPE -> numHaplo1 = i + g;
      			printf("pfic -> nomIndividu : %s\n", pfic -> nomIndividu);
        	    printf("pPE -> numHaplo1 : %d\n", pPE -> numHaplo1);
        	    pPE -> numHaplo2 = l;
        	    printf("pPE -> numHaplo2 : %d\n", pPE -> numHaplo2);
        	    l--;
        	    if (i != (nbHaploExplicatif/2 - 1)) {
        	        pPE -> suiv = (TypePaireExplicative*)malloc(sizeof(TypePaireExplicative));
        	        pPE = pPE -> suiv;
        	        pPE -> suiv = NULL;
        	    }
        	}
        	
        	pHaplo -> numHaplo = nbHaploTotaux;
        	pHaplo -> cpt_haplo = 1;
        	
        	for (j=0 ; j<tailleGeno ; j++) {
        		pHaplo -> haplotype[j] = tabTemp[j][i];
        	}
        	
        	if (pfic -> suiv != NULL || i != (nbHaploExplicatif - 1)) {
        		pHaplo -> suiv = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
        		pHaplo = pHaplo -> suiv;
        		pHaplo -> haplotype = (int*)malloc(sizeof(int)*tailleGeno);
        		pHaplo -> freq = 0.0;
        		pHaplo -> teteGenosExplicatifs = (TypeGenoExplicatif*)malloc(sizeof(TypeGenoExplicatif));
        		pHaplo -> suiv = NULL;
        	}
        	
        	nbHaploTotaux++;
        }
        
        // On garde 
        /*for (i=0 ; i<(nbHaploExplicatif/2) ; i++) {
        	l = nbHaploExplicatif - 1;
        	for (j=0 ; j<tailleGeno ; j++) {
        		
        			printf("tabTemp[j][i] : %d\n", tabTemp[j][i]);
        			printf("tabTemp[j][l] : %d\n", tabTemp[j][l]);
        	}
        	l--;
        	printf("\n");
        }
        printf("############################################\n");*/

        pfic = pfic -> suiv;
    }

    // Libération de la mémoire
    for (i=0 ; i<tailleGeno ; i++) {
        free(tabTemp[i]);
    }
    free(tabTemp); 
	
	return nbHaploTotaux;
    //FIN
}


void initialisation_freq_haplotype (TypeListeHaplotypes tete, int nb_haplo, int tailleGeno){
	
    // Variables locales
	double freqCalc;
	int i;
	TypeListeHaplotypes p;

    //DEBUT
    p = tete;
	freqCalc=(1.0/nb_haplo);
	//parcour de la liste des haplotype 
	while(p != NULL ){
		//on rempli la valeur de la frequence
		p -> freq = freqCalc;
		p = p -> suiv;
	}	
	// Affichage de la liste des haplotypes explicatifs
	printf("\n*************\n");
	printf("Affichage de la liste des haplotypes\n");
	printf("\n*************\n");
	p = tete;
	while (p != NULL) {
		printf("numHaplo : %d\n", p->numHaplo);
		for (i=0; i<tailleGeno ; i++) {
			printf("%d", p -> haplotype[i]);
		}
		printf("\ncpt_haplo : %d\n", p->cpt_haplo);
		printf("freq : %le\n", p->freq);
		printf("\n*************\n");
		p=p->suiv;
	}
	
    //FIN
}

/*
 *##### MAIN #####
 */

int main(int argc, char* argv[]) {

	printf("############################################\n");
	printf("Début du programme Inference_Haplotype\n");
	printf("############################################\n");

    //VARIABLE
    int nbHaploTotaux;
    char* fichier;
    char* nbIndividu;
    char* tailleGeno;
    TypeFichier ptrFichier;
    TypeListeHaplotypes ptrListeHaplo;

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
    ptrListeHaplo = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
    if(ptrListeHaplo==NULL){
        printf("L'allocation mémoire a échoué\n");
        exit(1);
    }
    
    nbHaploTotaux = predictionEspaceRestraintHaplotypes(ptrFichier, atoi(tailleGeno), atoi(nbIndividu), ptrListeHaplo);
    
    // On ajoute dans la liste des haplotypes la fréquence initiale
    initialisation_freq_haplotype(ptrListeHaplo, nbHaploTotaux, atoi(tailleGeno));

    free(ptrFichier);
    free(ptrListeHaplo);
    return 0;
    // FIN
}
