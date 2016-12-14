/*
 * Inference_Haplotype.c
 ***********************
 * Auteurs : Lea Boulongne et Estelle Geffard
 * Universit?de Nantes
 * Master Bioinformatique
 * UE Algorithmique et programmation avancee pour les biologistes
 * UE M?hodes et algorithmes pour la bioinformatique
 * Encadrement C. Sinoquet
 * 2016 - 2017
 */

/* LIBRAIRIES ======================================================================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Structures_Inference_Haplotype.h"
#include "Signatures_Inference_Haplotype.h"

/* MAIN ============================================================================================================================= */
int main(int argc, char* argv[]) {

    /* VARIABLES GLOBALES */

    char* fichier;
    int i;
    int nbIndividu;
    int tailleGeno;
    int nbHaploTotaux;
    int nbEtapeMax = 20;
    double* tabFreq;
	TPtrHaplotypes* tabHaplo;
    TPtrGenotypes teteGenotypes;
    TPtrHaplotypes teteHaplotypes;
    TPtrIndividus teteIndividus;

    // Allocation et verification de la memoire dynamique
    teteGenotypes = (TPtrGenotypes)malloc(sizeof(TypeGenotype));

    if(teteGenotypes==NULL){
        printf("L'allocation memoire a echouee\n");
        exit(1);
    }

    teteHaplotypes = (TPtrHaplotypes)malloc(sizeof(TypeHaplotype));

    if(teteHaplotypes==NULL){
        printf("L'allocation memoire a echouee\n");
        exit(1);
    }

    teteIndividus = (TPtrIndividus)malloc(sizeof(TypeIndividu));

    if(teteIndividus==NULL){
        printf("L'allocation memoire a echouee\n");
        exit(1);
    }

    /* DEBUT */

    printf("\nInference d'haplotypes\n");
    printf("**********************\n\n");

    // Recuperation des arguments
    recupererArguments(argc, argv, &fichier, &nbIndividu, &tailleGeno);

    printf("Arguments utilises : \n");
    printf("- Lecture du fichier : %s\n", fichier);
    printf("- Nombre d'individu etudie : %d\n", nbIndividu);
    printf("- Taille du genotype utilisee : %d\n", tailleGeno);

    // Recuperation des donnees du fichier
    recupererDonneesFichier(teteIndividus, teteGenotypes, tailleGeno, nbIndividu, fichier);

    // Predication de l'espace restreint des haplotypes expliquant les genotypes
    nbHaploTotaux = predictionEspaceRestraintHaplotypes(teteGenotypes, teteHaplotypes, tailleGeno, nbIndividu);

    // Initialisation des fr?uences des haplotypes
    initialisationFreqHaplotype(teteHaplotypes, nbHaploTotaux, tailleGeno);

    // Calcul de la probabilite du genotype
    calculProbaGenotype(teteGenotypes);

    // Affichage des donnees
    printf("\n\n############################################LISTES AVANT############################################\n\n");
    afficher(teteIndividus, teteGenotypes, teteHaplotypes, tailleGeno, nbIndividu);

    // D?ut de l'agorithme it?atif EM
    algoEM(teteGenotypes, teteHaplotypes, nbEtapeMax, nbIndividu, tailleGeno);

    tabFreq = (double*)malloc(sizeof(double)*nbHaploTotaux);

    if(tabFreq==NULL){
        printf("L'allocation m?oire a ?hou?n\n");
        exit(1);
    }

    tabHaplo = (TPtrHaplotypes*)malloc(sizeof(TPtrHaplotypes*)*nbHaploTotaux);

    if(tabHaplo==NULL){
        printf("L'allocation m?oire a ?hou?n\n");
        exit(1);
    }

    // Trie des haplotypes par rapport leur frequence
    trieFreqHaplo(teteHaplotypes, tabFreq, tabHaplo, tailleGeno, nbHaploTotaux);

    //Ecriture des haplotypes tries par frequences dans un fichier
    ecritureListeHaploTrieFreqDecr(tabFreq , tabHaplo, nbHaploTotaux, tailleGeno);

    rechercheProbaMaxPairHaplo(teteIndividus, teteGenotypes, tailleGeno);

    // Affichage des donnees
    printf("\n\n############################################LISTES APRES############################################\n\n");
    afficher(teteIndividus, teteGenotypes, teteHaplotypes, tailleGeno, nbIndividu);

    // Liberation de la memoire
    free(tabFreq);
	free(tabHaplo);
    free(teteGenotypes);
    free(teteHaplotypes);
    free(teteIndividus);

    /* FIN */
}