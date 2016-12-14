/*
 * Structures_Inference_Haplotype.h
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

/* STRUCTURES ======================================================================================================================= */
typedef struct TypeGenoExplicatif* TPtrGenosExplicatifs;

typedef struct TypeHaplotype* TPtrHaplotypes;

typedef struct TypeHaplosExplicatifs* TPtrHaplosExplicatifs;

typedef struct TypeGenotype* TPtrGenotypes;

typedef struct TypeIndividu* TPtrIndividus;


typedef struct TypeGenoExplicatif {
    TPtrGenotypes genotype;
    TPtrHaplotypes haplotype;
    TPtrGenosExplicatifs suiv;
} TypeGenoExplicatif;


typedef struct TypeHaplotype {
    int* haplotype;
    int cpt;
    double freq;
    TPtrGenosExplicatifs teteGenosExplicatifs;
    TPtrHaplotypes suiv;
} TypeHaplotype;


typedef struct TypeHaplosExplicatifs {
	TPtrHaplotypes haplo1;
	TPtrHaplotypes haplo2;
	double p_part;
	TPtrHaplosExplicatifs suiv;
} TypeHaplosExplicatifs;


typedef struct TypeGenotype {
    int* genotype;
    double proba;
    int cpt;
    TPtrHaplosExplicatifs pPaireProbaMax;
    TPtrHaplosExplicatifs teteHaplosExplicatifs;
    TPtrGenotypes suiv;
} TypeGenotype;


typedef struct TypeIndividu {
	char* nomIndividu;
	TPtrGenotypes genotype;
	TPtrIndividus suiv;
} TypeIndividu;