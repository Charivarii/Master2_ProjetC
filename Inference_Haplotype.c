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
    double proba;
    int cpt;
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
    int cpt;
    double freq;
    TypeListeGenosExplicatifs teteGenosExplicatifs;
    TypeListeHaplotypes suiv;
} TypeHaplotype;

/*
 * ##### SIGNATURES FONCTIONS #####
 */

void recupererParametres(int nbArgument, char* listeArguments[], char** adr_fichier, char** adr_nbIndividu, char** adr_tailleGeno);

void recupererDonneesFichier(TypeFichier tete, int tailleGeno, int nbIndividu, char* fichier);

int recherche_redondance_geno(TypeFichier tete, int tailleGeno, int* tab);

int predictionEspaceRestraintHaplotypes(TypeFichier tete, int tailleGeno, int nbIndividu, TypeListeHaplotypes ptrListeHaplo);

void insertionTete(int numHaplo, TypeListeGenosExplicatifs* adr_teteGE, char* nomIndividu, TypeListePairesExplicatives tetePE);

int recherche_redondance_haplo(TypeListeHaplotypes tete, int tailleGeno, int* tab);

void initialisation_freq_haplotype (TypeListeHaplotypes tete, int nb_haplo, int tailleGeno);

void afficherListes(TypeFichier teteGeno, TypeListeHaplotypes teteHaplo, int tailleGeno);


/*
 *##### MAIN #####
 */

int main(int argc, char* argv[]) {

    printf("############################################\n");
    printf("Début du programme Inference_Haplotype\n");
    printf("############################################\n\n");

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
        printf("L'allocation mémoire a échoué\n\n");
        exit(1);
    }

    recupererDonneesFichier(ptrFichier, atoi(tailleGeno), atoi(nbIndividu), fichier);

    // Prédiction de l'espace restraint des haplotypes expliquant les génotypes
    ptrListeHaplo = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
    if(ptrListeHaplo==NULL){
        printf("L'allocation mémoire a échoué\n\n");
        exit(1);
    }
    
    nbHaploTotaux = predictionEspaceRestraintHaplotypes(ptrFichier, atoi(tailleGeno), atoi(nbIndividu), ptrListeHaplo);
    
    // On ajoute dans la liste des haplotypes la fréquence initiale
    initialisation_freq_haplotype(ptrListeHaplo, nbHaploTotaux, atoi(tailleGeno));

    // Afficher les données des structures
    afficherListes(ptrFichier, ptrListeHaplo, atoi(tailleGeno));

    free(ptrFichier);
    free(ptrListeHaplo);
    return 0;
    // FIN
}

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
        printf("Merci\n\n");
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
    int estRedondant = 0;
    int estNomIndividu = 0;
    FILE* fp = NULL;
    char c[1];
    TypeFichier ptr;

    //DEBUT
    // Ouverture et lecture du fichier
    fp = fopen(fichier, "r");
 
    if (fp == NULL) {
        printf("Le fichier %s n'a pas pu être ouvert\n\n", fichier);
        exit(1);
    } else {
        printf("Le fichier %s existe\n\n", fichier);

        ptr = tete;
        ptr -> nomIndividu = (char*)malloc(sizeof(char)*20);
        ptr -> genotype = (int*)malloc(tailleGeno*sizeof(int));
        ptr -> proba = 0.0;
        ptr -> cpt = 1;
        ptr -> tetePairesExplicatives = NULL;
        ptr -> suiv = NULL;

        do {
            c[0] = getc(fp);

            //On a atteint la fin de la ligne
            if (c[0] == '\n') {
                //On a écrit un génotype et on regarde s'il existe déjà dans notre liste TypeFichier, si oui on augmente le compteur
                if (tete -> suiv != NULL) {
                    estRedondant = recherche_redondance_geno(tete, tailleGeno, ptr -> genotype);
                    if (estRedondant) {
                        ptr -> cpt += 1;
                    }
                }
                //On passe à la ligne suivante donc on créer une nouvelle structure à notre liste
                if (i != (nbIndividu - 1)) {
                    ptr -> suiv = (TypeFichier)malloc(sizeof(TypeIndividu));
                    ptr = ptr -> suiv;
                    ptr -> nomIndividu = (char*)malloc(sizeof(char)*20);
                    ptr -> genotype = (int*)malloc(tailleGeno*sizeof(int));
                    ptr -> proba = 0.0;
                    ptr -> cpt = 1;
                    ptr -> tetePairesExplicatives = NULL;
                    ptr -> suiv = NULL;
                    i++;
                    estNomIndividu = 1;
                }

            //On a atteint l'espace qui sépare le nom de l'individu de son génotype
            } else if (c[0] == ' ') {
                estNomIndividu = 0;
                j = 0;

            //On a rencontré un caractère à garder donc on l'écrit soit dans le nom soit dans le génotype en fonction de la variable estNomindividu
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


int recherche_redondance_geno(TypeFichier tete, int tailleGeno, int* tab) {

    //VARIABLES LOCALES
    TypeFichier p = tete;
    int identique, i;
    int estRedondant = 0;

    //DEBUT
    while (p != NULL && p -> suiv != NULL) { //parcours la liste
        identique=0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (tab[i] == p -> genotype[i]) {
                identique++;
            }
        }
        if (identique == tailleGeno) {
            p -> cpt++;
            estRedondant = 1;
        }
        p = p -> suiv;
    }
    return estRedondant;
    //FIN
}


int predictionEspaceRestraintHaplotypes(TypeFichier tete, int tailleGeno, int nbIndividu, TypeListeHaplotypes ptrListeHaplo) {

    printf("\nDébut de la fonction predictionEspaceRestraintHaplotypes\n");
    printf("############################################\n\n");

    //VARIABLES LOCALES
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
    int estRedondant = 0;
    TypeFichier pfic;
    int* vecTemp;
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

        // On remplit la structure des paires explicatives pour chaque génotype
        pfic -> tetePairesExplicatives = (TypePaireExplicative*)malloc(sizeof(TypePaireExplicative));
        pPE = pfic -> tetePairesExplicatives;
        pPE -> suiv = NULL;
        for (i=0 ; i<(nbHaploExplicatif) ; i++) {
            /*if (i < (nbHaploExplicatif/2)) {
                pPE -> numHaplo1 = i + g;
                pPE -> numHaplo2 = l;
                l--;
                if (i != (nbHaploExplicatif/2 - 1)) {
                    pPE -> suiv = (TypePaireExplicative*)malloc(sizeof(TypePaireExplicative));
                    pPE = pPE -> suiv;
                    pPE -> suiv = NULL;
                }
            }*/
            vecTemp = (int*)malloc(sizeof(int)*(tailleGeno));
            for (j=0 ; j<tailleGeno ; j++) {
                vecTemp[j] = tabTemp[j][i];
            }

            //####### A REVOIR #####
            // Appelle de la fonction de redondance
            if (nbHaploTotaux != 0) {
                estRedondant = recherche_redondance_haplo(ptrListeHaplo, tailleGeno, vecTemp);
                printf("estRedondant : %d\n", estRedondant);
            }
            if (!estRedondant) {
            	if (nbHaploTotaux > 0) {
            		pHaplo = pHaplo -> suiv;
            		printf("pseudo1\n");
            	}
            	printf("1\n");
            	pHaplo -> haplotype = (int*)malloc(sizeof(int)*tailleGeno);
            	printf("2\n");
        		pHaplo -> freq = 0.0;
        		pHaplo -> cpt = 1;
        		pHaplo -> teteGenosExplicatifs = (TypeGenoExplicatif*)malloc(sizeof(TypeGenoExplicatif));
        		pHaplo -> suiv = NULL;

                pHaplo -> numHaplo = nbHaploTotaux;
                for (j=0 ; j<tailleGeno ; j++) {
                    pHaplo -> haplotype[j] = vecTemp[j];
                }

                //####### A REVOIR #####

                // On remplit la structure genoExplicatif, insertion en tete
                //insertionTete(pHaplo -> numHaplo, &pHaplo -> teteGenosExplicatifs, pfic -> nomIndividu, pfic -> tetePairesExplicatives);

                nbHaploTotaux++;

            } else {

    

            }

            free(vecTemp);
            
        }

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



int recherche_redondance_haplo(TypeListeHaplotypes tete, int tailleGeno, int* tab) {

    //VARIABLES LOCALES
    TypeListeHaplotypes p = tete;
    int identique, i;
    int estRedondant = 0;

    //DEBUT
    while (p != NULL && p -> suiv != NULL) { //parcours la liste
        identique=0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (tab[i] == p -> haplotype[i]) {
                identique++;
            }
        }
        if (identique == tailleGeno) {
            p -> cpt++;

            estRedondant = 1;
        }
        p = p -> suiv;
    }
    return estRedondant;
    //FIN
}



void insertionTete(int numHaplo, TypeListeGenosExplicatifs* adr_teteGE, char* nomIndividu, TypeListePairesExplicatives tetePE) {

    //VARIABLES LOCALES
    TypeListeGenosExplicatifs p_new = (TypeGenoExplicatif*)malloc(sizeof(TypeGenoExplicatif));
    TypeListePairesExplicatives ptr = tetePE;
    //DEBUT
    p_new -> nomIndividu = nomIndividu;
    while (ptr != NULL) {
        if (numHaplo == ptr -> numHaplo1) {
            p_new -> numHaplo2 = ptr -> numHaplo2;
            p_new -> suiv = *adr_teteGE;
            *adr_teteGE = p_new;
            break;
        }
        if (numHaplo == ptr -> numHaplo2) {
            p_new -> numHaplo2 = ptr -> numHaplo1;
            p_new -> suiv = *adr_teteGE;
            *adr_teteGE = p_new;
            break;
        }
        ptr = ptr -> suiv;
    }
}


void initialisation_freq_haplotype(TypeListeHaplotypes tete, int nb_haplo, int tailleGeno) {
    
    //VARIABLES LOCALES
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
    //FIN
}



void afficherListes(TypeFichier teteGeno, TypeListeHaplotypes teteHaplo, int tailleGeno) {

    //VARIABLES LOCALES
    TypeFichier pGeno;
    TypeListeHaplotypes pHaplo;
    TypeListePairesExplicatives pPE;
    TypeListeGenosExplicatifs pGE;
    int i;

    //DEBUT
    printf("Affichage de la liste contenant les génotypes de chaque individu\n");
    printf("********************************************************************\n");
    pGeno = teteGeno;
    while (pGeno != NULL) {
        printf("nom : %s\ngenotype : ", pGeno->nomIndividu);
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", pGeno-> genotype[i]);
        }
        printf("\nproba : %le\n", pGeno->proba);
        printf("cpt : %d\n", pGeno->cpt);
        printf("\nAffichage des paires explicatives\n");
        pPE = pGeno -> tetePairesExplicatives;
        while (pPE != NULL) {
            printf("numHaplo1 : %d\n", pPE -> numHaplo1);
            printf("numHaplo2 : %d\n", pPE -> numHaplo2);
            pPE = pPE -> suiv;
        }
        printf("-----------------------------------------\n");
        pGeno=pGeno->suiv;
    }

    printf("\nAffichage de la liste contenant les haplotypes de l'espace restraint\n");
    printf("********************************************************************\n");
    pHaplo = teteHaplo;
    while (pHaplo != NULL) {
        printf("numHaplo : %d\nhaplotype : ", pHaplo->numHaplo);
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", pHaplo-> haplotype[i]);
        }
        printf("\nfreq : %le\n", pHaplo->freq);
        printf("cpt : %d\n", pHaplo->cpt);
        printf("\nAffichage des genotypes expliquant l'haplotype\n");
        pGE = pHaplo -> teteGenosExplicatifs;
        while (pGE != NULL) {
            printf("nomIndividu : %s\n", pGE -> nomIndividu);
            printf("numHaplo2 : %d\n", pGE -> numHaplo2);
            pGE = pGE -> suiv;
        }
        printf("-----------------------------------------\n");
        pHaplo=pHaplo->suiv;
    }

    //FIN
}