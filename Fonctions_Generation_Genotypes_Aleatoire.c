/*
 * Fonctions_Generation_Genotypes_Aleatoire.c
 ********************************************
 * Auteurs : Lea Boulongne et Estelle Geffard
 * Université de Nantes
 * Master Bioinformatique
 * UE Algorithmique et programmation avancée pour les biologistes
 * UE Méthodes et algorithmes pour la bioinformatique
 * Encadrement C. Sinoquet
 * 2016 - 2017
 */

/* LIBRAIRIES ======================================================================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Signatures_Generation_Genotypes_Aleatoire.h"

/* FONCTIONS ======================================================================================================================== */
void recupererArguments(int nbArgument, char* listeArguments[], char** adr_nbIndividu, char** adr_tailleGeno , char** adr_pc_ambigue) {

    /* DEBUT */

    // L'utilisateur doit rentrer 3 arguments
    if (nbArgument != 4) {
        printf("Au lancement du programme, vous devez rentrer 3 arguments : \n");
        printf("- Le nombre d'individu souhaité\n");
        printf("- La longueur du génotype voulue\n");   
        printf("- Le pourcentage d'ambigue maximum : (50 max) \n");
        printf("Merci\n");
        exit(1);
    }

    // nbArguments == 4
    *adr_nbIndividu = listeArguments[1];
    *adr_tailleGeno = listeArguments[2];
    *adr_pc_ambigue = listeArguments[3];

    /* FIN */
}

/* ******************************************************************************************************************************** */
void generation_haplotype(int nbInd, int tailleGeno, int nbAmbigue, int* haplo1, int* haplo2) {

    /* VARIABLES LOCALES */

    int i, j, k, l, nbDeux;
    int* tabTemp;

    /* DEBUT */

    // Generation aleatoire du premier haplotype avec des 0 et des 1 et copie dans le deuxieme haplotype
    for (j = 0 ; j < tailleGeno ; j++) {
        haplo1[j] = rand()%2;
        haplo2[j] = haplo1[j];
    }

    // Tirage aleatoire du nombre de loci ambigues par rapport au nombre maximum
    nbDeux = rand()%(nbAmbigue + 1);

    // Si le tirage aleatoire sort un nombre de loci ambigues different de deux, on modifie le deuxieme haplotype
    if (nbDeux != 0) {

        // Allocation de memoire
        tabTemp = (int*)malloc(nbDeux*sizeof(int));

        // Verification des allocations de memoire dynamique
        if(tabTemp == NULL) {
            printf("Problème d'allocation de mémoire dynamique");
            exit(1);
        }

        // Tirage aleatoire des positions a changer dans le deuxieme haplotype
        for (k = 0 ; k < nbDeux ; k++) {
            tabTemp[k] = rand()%(tailleGeno);
        }

        // Modification de l'haplotype a la position tiree aleatoirement
        for (k = 0 ; k < nbDeux ; k++) {
            for (l = 0 ; l < k ; l++) {
                if (tabTemp[k] == tabTemp[l]) {
                    break; // Si la position est la meme, on quitte la boucle
                }
            }
            if (haplo2[tabTemp[k]] == 0) {
                haplo2[tabTemp[k]] = 1;
            } else {
                haplo2[tabTemp[k]] = 0;
            }
        }

        // Liberation de la memoire
        free(tabTemp);
    }

    /* FIN */   
}

/* ******************************************************************************************************************************** */
void generation_genotype(int tailleGeno, int* haplo1, int* haplo2, int* geno) {

    /* VARIABLES LOCALES */

    int i, j;

    /* DEBUT */

    for (i = 0 ; i < tailleGeno ; i++) {
        if (haplo2[i] != haplo1[i]) {
            geno[i] = 2;
        } else if (haplo2[i] == 1) {
            geno[i] = 1;
        } else {
            geno[i] = 0;
        }
    }

    /* FIN */
}

/* ******************************************************************************************************************************** */
void ecriture_fichier_geno_haplo(int i, int tailleGeno, int* geno, int* haplo1, int* haplo2) {

    /* VARIABLES LOCALES */

    FILE* fichier = NULL;
    int j;

    /* DEBUT */

    fichier = fopen("genotypesHaplotypesAleatoires.txt", "a");

    if (fichier != NULL) {
        fprintf(fichier, "/ind%d geno ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", geno[j]);
        }
        fprintf(fichier, "\n");
        fprintf(fichier, "/ind%d real haplo1 ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", haplo1[j]);
        }
        fprintf(fichier, "\n");
        fprintf(fichier, "/ind%d real haplo2 ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", haplo2[j]);
        }
        fprintf(fichier, "\n");
        fclose(fichier);
    }

    /* FIN */
}

/* ******************************************************************************************************************************** */
void ecriture_fichier_geno(int i, int tailleGeno, int* geno) {

    /* VARIABLES LOCALES */

    FILE* fichier = NULL;
    int j;

    /* DEBUT */

    fichier = fopen("genotypesAleatoires.txt", "a");

    if (fichier != NULL) {
        fprintf(fichier, "/ind%d ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", geno[j]);
        }
        fprintf(fichier, "\n");
        fclose(fichier);
    }

    /* FIN */
}
