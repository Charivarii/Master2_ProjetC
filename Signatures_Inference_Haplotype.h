/*
 * Signatures_Inference_Haplotype.h
 **********************************
 * Auteurs : Lea Boulongne et Estelle Geffard
 * Universit?de Nantes
 * Master Bioinformatique
 * UE Algorithmique et programmation avanc? pour les biologistes
 * UE M?hodes et algorithmes pour la bioinformatique
 * Encadrement C. Sinoquet
 * 2016 - 2017
 */

/* LIBRAIRIES ======================================================================================================================= */
#include <stdio.h>
#include <stdlib.h>

/* SIGNATURES ======================================================================================================================= */
void recupererArguments(int nbArgument, char* listeArguments[], char** adr_fichier, int* adr_nbIndividu, int* adr_tailleGeno);

void recupererDonneesFichier(TPtrIndividus teteIndividus, TPtrGenotypes teteGenotypes, int tailleGeno, int nbIndividu, char* fichier);

int rechercheRedondanceGeno(TPtrGenotypes teteGenotypes, TPtrIndividus pIndividus, int tailleGeno, int* tempGenotype);

int predictionEspaceRestraintHaplotypes(TPtrGenotypes teteGenotypes, TPtrHaplotypes teteHaplotypes, int tailleGeno, int nbIndividu);

int rechercheRedondanceHaplo(TPtrHaplotypes teteHaplotypes, int tailleGeno, int* tab);

void creationListeGenoExplicatif(TPtrHaplotypes teteHaplotypes, TPtrGenotypes teteGenotypes);

void insertionTete(TPtrGenosExplicatifs* adr_teteGE, TPtrHaplotypes pHaplotypes, TPtrGenotypes pGenotypes);

void initialisationFreqHaplotype(TPtrHaplotypes teteHaplotypes, int nb_haplo, int tailleGeno);

void calculProbaGenotype(TPtrGenotypes teteGenotype);

double calculProba (TPtrHaplosExplicatifs pHE);

void algoEM(TPtrGenotypes teteGenotype, TPtrHaplotypes pHaplo, int nbEtapeMax, int nbIndividu, int tailleGeno);

void maximisation (int nbIndividu, TPtrHaplotypes teteHaplo);

double estimationEsperance(TPtrGenotypes teteGenotype);

void trieFreqHaplo(TPtrHaplotypes teteHaplo, double* tabFreq, TPtrHaplotypes* tabHaplo, int tailleGeno, int nbHaploTotaux);

void quickSort(double* tabFreq, TPtrHaplotypes* tabHaplo, int g, int d, int tailleGeno);

int partitionner(double* tabFreq, TPtrHaplotypes* tabHaplo, int g, int d, int tailleGeno);

void echange(double* tabFreq, TPtrHaplotypes* tabHaplo, int i, int j, int tailleGeno);

void ecritureListeHaploTrieFreqDecr(double* vect_freq , TPtrHaplotypes* mat_haplo, int nbHaploTotaux, int tailleGeno);

void rechercheProbaMaxPairHaplo(TPtrIndividus teteIndividu, TPtrGenotypes teteGenotypes, int tailleGeno);

void ecritureListeIndGenoPaireHaplo(TPtrIndividus pInd, int tailleGeno);

void afficher(TPtrIndividus teteIndividus, TPtrGenotypes teteGenotypes, TPtrHaplotypes teteHaplotypes, int tailleGeno, int nbIndividu);