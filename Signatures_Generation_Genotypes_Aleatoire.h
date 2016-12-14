/*
 * Signatures_Generation_Genotypes_Aleatoire.h
 *********************************************
 * Auteurs : Lea Boulongne et Estelle Geffard
 * Université de Nantes
 * Master Bioinformatique
 * UE Algorithmique et programmation avancée pour les biologistes
 * UE Méthodes et algorithmes pour la bioinformatique
 * Encadrement C. Sinoquet
 * 2016 - 2017
 */

/* SIGNATURES ======================================================================================================================= */
void recupererArguments(int nbArgument, char* listeArguments[], char** adr_nbIndividu, char** adr_tailleGeno , char** adr_pc_ambigue);

void generation_haplotype(int nbInd, int tailleGeno, int nbAmbigue, int* haplo1, int* haplo2);

void generation_genotype(int tailleGeno, int* haplo1, int* haplo2, int* geno);

void ecriture_fichier_geno_haplo(int i, int tailleGeno, int* geno, int* haplo1, int* haplo2);

void ecriture_fichier_geno(int i, int tailleGeno, int* geno);