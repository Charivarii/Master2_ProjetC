/*
 * Generation_Genotypes_Aleatoire.c
 ***********************************
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

/* MAIN ============================================================================================================================= */
int main(int argc, char* argv[]) {

    /* VARIABLES GLOBALES */

    int i, j, nbAmbigue, nbInd, tailleGeno, pcAmbigue;
    int* haplo1;
    int* haplo2;
    int* geno;
    char*  nbIndC;
    char* tailleGenoC;
    char* pcAmbigueC;

    /* DEBUT */

    printf("\nGeneration de genotypes aleatoires\n");
    printf("**********************************\n\n");

    // Recuperation des arguments
    recupererArguments(argc, argv, &nbIndC, &tailleGenoC, &pcAmbigueC);

    // Passage de char en int
    nbInd = atoi(nbIndC);
    tailleGeno = atoi(tailleGenoC);
    pcAmbigue = atoi(pcAmbigueC);
	
    printf("Arguments : \n");
    printf("- Nombre d'individus : %d\n", nbInd);
    printf("- Taille du genotype : %d\n", tailleGeno);
    printf("- Pourcentage d'ambiguite : %d\n", pcAmbigue);

    // Allocation dynamique de memoire
    haplo1 = (int*)malloc(tailleGeno*sizeof(int));
    haplo2 = (int*)malloc(tailleGeno*sizeof(int));
    geno = (int*)malloc(tailleGeno*sizeof(int));

    // Verification des allocations
    if(haplo1 == NULL || haplo2 == NULL  || geno == NULL ){
        printf("Probleme d'allocation de memoire dynamique\n");
        exit(1);
    }

    // Supression des fichiers de sauvegarde s'ils existent deja
    remove ("genotypesAleatoires.txt");
    remove ("genotypesHaplotypesAleatoires.txt");

    if (pcAmbigue > 50) { // On verifie si le nombre d'ambigue n'est pas trop grand
        printf("Pourcentage d'ambigue choisit trop grand\n");
        exit(1);
    } else {
        nbAmbigue = (pcAmbigue * tailleGeno)/100;
        printf("- Nombre d'ambigues dans le genotype : %d\n", nbAmbigue);
    }

    // Pour chaque individu, on appelle les fonctions pour generer les haplotypes et le genotype
    for (i=0 ; i<nbInd ; i++) {
        generation_haplotype(nbInd, tailleGeno, nbAmbigue, haplo1, haplo2);
        generation_genotype(tailleGeno, haplo1, haplo2, geno);
        ecriture_fichier_geno_haplo(i, tailleGeno, geno, haplo1, haplo2);
        ecriture_fichier_geno(i, tailleGeno, geno);
    }

    printf("\nGeneration des haplotypes et genotypes effectues\n\n");

    // Liberation de la memoire
    free(haplo1);
    free(haplo2);
    free(geno);

    return 0;

    /* FIN */
}