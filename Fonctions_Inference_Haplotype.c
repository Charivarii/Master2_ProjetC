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

/* FONCTIONS ======================================================================================================================== */
void recupererArguments(int nbArgument, char* listeArguments[], char** adr_fichier, int* adr_nbIndividu, int* adr_tailleGeno) {

    /* DEBUT */

    // L'utilisateur doit rentrer 3 arguments
    if (nbArgument != 4) {
        printf("Au lancement du programme, vous devez rentrer 3 arguments : \n");
        printf("- Le nom du fichier contenant les genotypes des individus\n");
        printf("- Le nombre d'individu eudie");
        printf("- La longueur du genotype etudie\n");
        printf("Merci\n\n");
        exit(1);
    }

    // nbArguments == 4
    *adr_fichier = listeArguments[1];
    *adr_nbIndividu = atoi(listeArguments[2]);
    *adr_tailleGeno = atoi(listeArguments[3]);

    /* FIN */
}

/* ******************************************************************************************************************************** */
void recupererDonneesFichier(TPtrIndividus teteIndividus, TPtrGenotypes teteGenotypes, int tailleGeno, int nbIndividu, char* fichier) {

    /* VARIABLES LOCALES */

    int i = 0;
    int j = 0;
    int nbChar = 1;
    int estRedondant = 0;
    int estNomIndividu = 1;
    int estPremierElement = 1;
    int* tempGenotype;
    char c[1];
    FILE* fp = NULL;
    TPtrIndividus pIndividus;
    TPtrGenotypes pGenotypes;
    TPtrGenotypes pGenoRedondant;

    /* DEBUT */

    // Ouverture et lecture du fichier
    fp = fopen(fichier, "r");

    if (fp == NULL) {
        printf("Le fichier %s n'a pas pu etre ouvert\n\n", fichier);
        exit(1);
    }

    pIndividus = teteIndividus;
    pIndividus -> nomIndividu = (char*)calloc(nbChar, sizeof(char));
    pIndividus -> genotype = (TPtrGenotypes)malloc(sizeof(TypeGenotype));
    pIndividus -> suiv = NULL;

    pGenotypes = teteGenotypes;
    pGenotypes -> genotype = (int*)malloc(tailleGeno*sizeof(int));
    pGenotypes -> proba = 0.0;
    pGenotypes -> cpt = 1;
    pGenotypes -> pPaireProbaMax = NULL;
    pGenotypes -> teteHaplosExplicatifs = NULL;
    pGenotypes -> suiv = NULL;

    // On recupere le premier caractere du fichier
    c[0] = getc(fp);

    while (c[0] != EOF) { // On parcours jusqu'a la fin du fichier
        if (c[0] == '\n') { // Fin de la ligne
            if (estPremierElement) { // Premier element dans la liste, donc on ecrit la structure
            	pIndividus -> genotype = pGenotypes;
                for (j=0 ; j<tailleGeno ; j++) {
                	pGenotypes -> genotype[j] = tempGenotype[j];
                }
                estPremierElement = 0;
            } else { // Deja des elements dans la liste, on verifie la redondance
				estRedondant = rechercheRedondanceGeno(teteGenotypes, pIndividus, tailleGeno, tempGenotype);
                if (!estRedondant) {
                	pGenotypes -> suiv = (TPtrGenotypes)malloc(sizeof(TypeGenotype));
                	pGenotypes = pGenotypes -> suiv;
                	pGenotypes -> genotype = (int*)malloc(tailleGeno*sizeof(int));
    				pGenotypes -> proba = 0.0;
    				pGenotypes -> cpt = 1;
    				pGenotypes -> pPaireProbaMax = NULL;
    				pGenotypes -> teteHaplosExplicatifs = NULL;
    				pGenotypes -> suiv = NULL;
                	pIndividus -> genotype = pGenotypes;
                	for (j=0 ; j<tailleGeno ; j++) {
                		pGenotypes -> genotype[j] = tempGenotype[j];
                	}
            	}
            }

            // Ligne suivante donc creation d'une nouvelle structure a notre liste
            if (i != (nbIndividu - 1)) {
                pIndividus -> suiv = (TPtrIndividus)malloc(sizeof(TypeIndividu));
                pIndividus = pIndividus -> suiv;
                pIndividus -> nomIndividu = (char*)calloc(nbChar, sizeof(char));
    			pIndividus -> genotype = (TPtrGenotypes)malloc(sizeof(TypeGenotype));
    			pIndividus -> suiv = NULL;
                i++;
                estNomIndividu = 1;
            }

            free(tempGenotype);

        // Le caractere espace separe le nom de l'individu de son genotype
        } else if (c[0] == ' ') {
        	tempGenotype = (int*)malloc(sizeof(int)*tailleGeno);
        	j = 0;
            estNomIndividu = 0;

        //Caractere a conserver donc ecriture soit dans le nom soit dans le genotype en fonction de la variable estNomindividu
        } else {
            if (estNomIndividu) {
            	nbChar++;
                pIndividus -> nomIndividu = (char*)realloc(pIndividus -> nomIndividu, nbChar);
                pIndividus -> nomIndividu = strcat(pIndividus -> nomIndividu, c);
            } else {
                tempGenotype[j] = atoi(c);
                j++;
            }
        }
        c[0] = getc(fp);
    }

    /* FIN */

}

/* ******************************************************************************************************************************** */
int rechercheRedondanceGeno(TPtrGenotypes teteGenotypes, TPtrIndividus pIndividus, int tailleGeno, int* tempGenotype) {

    /* VARIABLES LOCALES */

    TPtrGenotypes pGeno = teteGenotypes;
    int identique, i;
    int estRedondant = 0;

    /* DEBUT */

    while (pGeno != NULL && pGeno -> suiv != NULL) { // Parcours la liste
        identique=0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (tempGenotype[i] == pGeno -> genotype[i]) {
                identique++;
            }
        }
        if (identique == tailleGeno) {
            pGeno -> cpt += 1;
            pIndividus -> genotype = pGeno;
            estRedondant = 1;
        }
        pGeno = pGeno -> suiv;
    }
    return estRedondant;

    /* FIN */

}

/* ******************************************************************************************************************************** */
int predictionEspaceRestraintHaplotypes(TPtrGenotypes teteGenotypes, TPtrHaplotypes teteHaplotypes, int tailleGeno, int nbIndividu) {

	/* VARIABLES LOCALES */

    int i = 0;
    int k = 0;
    int j = 0;
    int t = 0;
    int s = 0;
    int l = 0;
    int identique1, identique2;
    int nbHaploTotaux = 0;
    int pas = 0;
    int multiplicateur = 0;
    int nbDeux = 0;
    int nbHaploExplicatif = 0;
    int estRedondant = 0;
    int* vecTemp;
    int* vecTempHaplo2;
    int** tabTemp;
    TPtrGenotypes pGenotypes;
    TPtrHaplotypes pHaplotypes, pHaplopHE;
    TPtrHaplosExplicatifs pHE;

    /* DEBUT */

    // Allocation memoire du tableau temporaire des paires d'haplotypes possibles expliquant le genotype
    tabTemp = (int**)malloc(sizeof(int*)*tailleGeno);
    for (i=0 ; i<tailleGeno ; i++) {
        tabTemp[i] = (int*)malloc(sizeof(int));
    }

    // Parcours de la liste pour recuperer les genotypes et predire l'espace restaint des haplotypes
    pGenotypes = teteGenotypes;

    // Initialisation de la liste qui contiendra les haplotypes explicatifs
    pHaplotypes = teteHaplotypes;

    while (pGenotypes != NULL) { // Pour chaque genotype
        // On compte le nombre de 2 dans le genotype
        nbDeux = 0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (pGenotypes -> genotype[i] == 2) {
                nbDeux++;
            }
        }

        if (nbDeux == 0) {
            nbDeux = 1;
        }
        nbHaploExplicatif = pow(2, nbDeux);

        // Reallocation memoire
        for (i=0 ; i<tailleGeno ; i++) {
            tabTemp[i] = (int*)realloc(tabTemp[i], sizeof(int)*nbHaploExplicatif);
        }

        // Generation des paires d'haplotypes expliquant le genotype
        k = 0; // k : nombre de deux rencontres jusqu'a maintenant
        for (i=0 ; i<tailleGeno ; i++) { // Parcours du g?otype
            if (pGenotypes -> genotype[i] == 1) {
                for (j=0 ; j<nbHaploExplicatif ; j++) {
                    tabTemp[i][j] = 1;
                }
            } else if (pGenotypes -> genotype[i] == 0) {
                for (j=0 ; j<nbHaploExplicatif ; j++) {
                    tabTemp[i][j] = 0;
                }
            } else {
                k++; // On rencontre un 2 dans la sequence de genotype
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

        // On garde les haplotypes dans une liste chainee TPtrHaplotypes et les paires explicatives dans TPtrPairesExplicatives

        // On remplit la liste chainee regroupant tous les haplotypes de l'espace restreint des solutions
        for (i=0 ; i<(nbHaploExplicatif) ; i++) {

            vecTemp = (int*)malloc(sizeof(int)*(tailleGeno));

            for (j=0 ; j<tailleGeno ; j++) {
                vecTemp[j] = tabTemp[j][i];
            }

            // Appelle de la fonction de redondance
            if (nbHaploTotaux != 0) {
                estRedondant = rechercheRedondanceHaplo(teteHaplotypes, tailleGeno, vecTemp);
            }
            if (!estRedondant) {
            	if (nbHaploTotaux > 0) {
            		pHaplotypes -> suiv = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
            		pHaplotypes = pHaplotypes -> suiv;
            	}
            	pHaplotypes -> haplotype = (int*)malloc(sizeof(int)*(tailleGeno));
        		pHaplotypes -> freq = 0.0;
        		pHaplotypes -> cpt = 1;
        		pHaplotypes -> teteGenosExplicatifs = NULL;
        		pHaplotypes -> suiv = NULL;

                for (j=0 ; j<tailleGeno ; j++) {
                    pHaplotypes -> haplotype[j] = vecTemp[j];
                }

                nbHaploTotaux++;

            }

            free(vecTemp);

        }

        // On initialise la liste chainee qui contiendra les structures avec les paires d'haplotypes explicatifs
        pGenotypes -> teteHaplosExplicatifs = (TypeHaplosExplicatifs*)malloc(sizeof(TypeHaplosExplicatifs));
        pHE = pGenotypes -> teteHaplosExplicatifs;
        pHE -> haplo1 = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
        pHE -> haplo2 = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
        pHE -> p_part = 0.0;
        pHE -> suiv = NULL;

        l = 0;
        l = (nbHaploExplicatif - 1);

        // On remplit la structure des paires explicatives pour chaque genotype
        for (i=0 ; i<(nbHaploExplicatif/2) ; i++) {

        	vecTemp = (int*)malloc(sizeof(int)*(tailleGeno));
            vecTempHaplo2 = (int*)malloc(sizeof(int)*(tailleGeno));
           	for (j=0 ; j<tailleGeno ; j++) {
           		vecTemp[j] = tabTemp[j][i];
               	vecTempHaplo2[j] = tabTemp[j][l];
            }

            if (i != 0) {
            	pHE -> suiv = (TypeHaplosExplicatifs*)malloc(sizeof(TypeHaplosExplicatifs));
                pHE = pHE -> suiv;
                pHE -> haplo1 = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
        		pHE -> haplo2 = (TypeHaplotype*)malloc(sizeof(TypeHaplotype));
                pHE -> p_part = 0.0;
                pHE -> suiv = NULL;
            }
            pHaplopHE = teteHaplotypes;
            while (pHaplopHE != NULL) {
            	identique1 = 0;
            	identique2 = 0;
        		for (j=0 ; j<tailleGeno ; j++) {
            		if (vecTemp[j] == pHaplopHE -> haplotype[j]) {
                		identique1++;
            		}
            		if (vecTempHaplo2[j] == pHaplopHE -> haplotype[j]) {
                		identique2++;
            		}
			    }
			    if (identique1 == tailleGeno) {
			        pHE -> haplo1 = pHaplopHE;
			    }
			    if (identique2 == tailleGeno) {
			        pHE -> haplo2 = pHaplopHE;
			    }
            	pHaplopHE = pHaplopHE -> suiv;
            }
            free(vecTemp);
           	free(vecTempHaplo2);
           	l--;
        }

        pGenotypes = pGenotypes -> suiv;
    }

    // Liberation de la memoire
    for (i=0 ; i<tailleGeno ; i++) {
        free(tabTemp[i]);
    }
    free(tabTemp);

    // On remplit la structure genoExplicatif, insertion en tete
    creationListeGenoExplicatif(teteHaplotypes, teteGenotypes);

    return nbHaploTotaux;

    /* FIN */
}

/* ******************************************************************************************************************************** */
int rechercheRedondanceHaplo(TPtrHaplotypes teteHaplotypes, int tailleGeno, int* tab) {

    /* VARIABLES LOCALES */

    TPtrHaplotypes pHaplotypes = teteHaplotypes;
    int identique, i;
    int estRedondant = 0;

    /* DEBUT */

    while (pHaplotypes != NULL) { // Parcours la liste

        identique=0;
        for (i=0 ; i<tailleGeno ; i++) {
            if (tab[i] == pHaplotypes -> haplotype[i]) {
                identique++;
            }
        }
        if (identique == tailleGeno) {
            pHaplotypes -> cpt++;
            estRedondant = 1;
            break;
        }
        pHaplotypes = pHaplotypes -> suiv;
    }
    return estRedondant;

    /* FIN */

}

/* ******************************************************************************************************************************** */
void creationListeGenoExplicatif(TPtrHaplotypes teteHaplotypes, TPtrGenotypes teteGenotypes) {

	/* VARIABLES LOCALES */

	TPtrHaplotypes pHaplotypes;
	TPtrGenotypes pGenotypes;
	TPtrHaplosExplicatifs pHE;
	pGenotypes = teteGenotypes;
	pHaplotypes = teteHaplotypes;


	/* DEBUT */

	while (pGenotypes != NULL) {
		pHE = pGenotypes -> teteHaplosExplicatifs;
		while (pHE != NULL) {
			insertionTete(&pHE -> haplo1 -> teteGenosExplicatifs, pHE -> haplo2, pGenotypes);
			insertionTete(&pHE -> haplo2 -> teteGenosExplicatifs, pHE -> haplo1, pGenotypes);
			pHE = pHE -> suiv;
		}
		pGenotypes = pGenotypes -> suiv;
	}

	/* FIN */

}

/* ******************************************************************************************************************************** */
void insertionTete(TPtrGenosExplicatifs* adr_teteGE, TPtrHaplotypes pHaplotypes, TPtrGenotypes pGenotypes) {

    /* VARIABLES LOCALES */

    TPtrGenosExplicatifs p_new = (TypeGenoExplicatif*)malloc(sizeof(TypeGenoExplicatif));

    /* DEBUT */

    p_new -> genotype = pGenotypes;
    p_new -> haplotype = pHaplotypes;
    p_new -> suiv = *adr_teteGE;
	*adr_teteGE = p_new;

	/* FIN */

}


/* ******************************************************************************************************************************** */
void initialisationFreqHaplotype(TPtrHaplotypes teteHaplotypes, int nb_haplo, int tailleGeno) {

    /* VARIABLES LOCALES */

    double freqCalc;
    int i;
    TPtrHaplotypes pHaplotypes;

    /* DEBUT */

    pHaplotypes = teteHaplotypes;
    freqCalc = (1.0/nb_haplo);

    while (pHaplotypes != NULL ){ // Parcours de la liste des haplotype
        // Calcul de la frequence et initialisation dans la liste d'haplotypes
        pHaplotypes -> freq = freqCalc;
        pHaplotypes = pHaplotypes -> suiv;
    }

    /* FIN */
}

/* ******************************************************************************************************************************** */
void calculProbaGenotype(TPtrGenotypes teteGenotype) {

	/* VARIABLES LOCALES */

	double proba = 0.0;
	TPtrGenotypes pGenotype = teteGenotype;

	/* DEBUT */

	while (pGenotype != NULL) { // Parcours la liste de genotype
		proba = calculProba(pGenotype -> teteHaplosExplicatifs);
		pGenotype -> proba = proba ;
		pGenotype = pGenotype -> suiv;
	}

	/* FIN */

}

/* ******************************************************************************************************************************** */
double calculProba (TPtrHaplosExplicatifs pHE) {

	/* VARIABLES LOCALES */

	double proba = 0.0;
	double p_part, freqH1, freqH2;
	TPtrHaplosExplicatifs pHaploExplicatifs = pHE;

	/* DEBUT */

	while (pHaploExplicatifs!= NULL) { // Parcours la liste de paire d'haplotypes explicatifs pour un g?otype
		if (pHaploExplicatifs -> haplo1 == pHaploExplicatifs -> haplo2) {
			p_part = pow(pHaploExplicatifs -> haplo1 -> freq, 2);
		} else {
	 		p_part = 2*(pHaploExplicatifs -> haplo1 -> freq) * (pHaploExplicatifs -> haplo2 -> freq);
		}
		pHaploExplicatifs -> p_part = p_part;
		proba += p_part;
		pHaploExplicatifs = pHaploExplicatifs -> suiv;
	}
	return(proba);

	/* FIN */

}

/* ******************************************************************************************************************************** */
void algoEM(TPtrGenotypes teteGenotype, TPtrHaplotypes pHaplo, int nbEtapeMax, int nbIndividu, int tailleGeno) {

	/* VARIABLES LOCALES */

	double log_likelihood, log_likelihood_prec;
	double seuil = 0.0001;
	int convergence = 0;
	int nbEtape = 0;

	/* DEBUT */

	log_likelihood_prec = 0;
	log_likelihood = 0;

	while (convergence != 1 && nbEtape < nbEtapeMax) {
		maximisation(nbIndividu, pHaplo);
		log_likelihood = estimationEsperance(teteGenotype);
		if ((fabs(log_likelihood_prec - log_likelihood)) < seuil) { //Si convergence
			convergence = 1;
		} else {
			log_likelihood_prec = log_likelihood;
		}
		nbEtape++;
	}
	printf("nbEtape : %d\n", nbEtape);

	/* FIN */

}

/* ******************************************************************************************************************************** */
void maximisation (int nbIndividu, TPtrHaplotypes teteHaplo) {

	/* VARIABLES LOCALES */

	double Ng, N, freq, freq_prec_H1, freq_prec_H2, prob_prec, contribution;
	TPtrGenotypes pGenotype;
	TPtrGenosExplicatifs pGenotypeExplicatif;
	TPtrHaplotypes pHaplo1 = teteHaplo;
	TPtrHaplotypes pHaplo2;

	N = (double)nbIndividu;

	/* DEBUT */

	while (pHaplo1 != NULL) { // Parcours la liste d'haplotype
		freq_prec_H1 = pHaplo1 -> freq;
		freq = 0.0;
		pGenotypeExplicatif = pHaplo1 -> teteGenosExplicatifs;
		while (pGenotypeExplicatif != NULL) { // Parcours liste de genotype correspondant ?un haplotype pour recuperer H2
			pHaplo2 = pGenotypeExplicatif -> haplotype;
			pGenotype = pGenotypeExplicatif -> genotype;
			Ng = (double)(pGenotype -> cpt);
			prob_prec = pGenotype -> proba;

			if (pHaplo1 == pHaplo2){ // Calcul de la contribution en fonction de l'homozygotie ou l'heterozygotie
				contribution = (2* ( (pow(freq_prec_H1, 2)) / (prob_prec) ) * (Ng/N) );
			} else {
				freq_prec_H2 = pHaplo2 -> freq;
				contribution = ( ( (2*(freq_prec_H1*freq_prec_H2)) / (prob_prec) ) * (Ng/N) );
			}
			freq += contribution;
			pGenotypeExplicatif =  pGenotypeExplicatif -> suiv;
		}

		freq = freq/2;
		pHaplo1 -> freq = freq;
		pHaplo1 = pHaplo1 -> suiv;
	}

	/* FIN */
}

/* ******************************************************************************************************************************** */
double estimationEsperance(TPtrGenotypes teteGenotype) {

	/* VARIABLES LOCALES */

	double log_likelihood = 0.0;
	double proba = 0.0;
	TPtrGenotypes pGeno = teteGenotype;

	/* DEBUT */

	while (pGeno != NULL) { // Parcours la liste de genotype
		proba = calculProba(pGeno -> teteHaplosExplicatifs);
		pGeno -> proba = proba;
		log_likelihood += (pGeno -> cpt) * log(proba);
		pGeno = pGeno -> suiv;
	}
	return log_likelihood;

	/* FIN */

}

/* ******************************************************************************************************************************** */
void trieFreqHaplo(TPtrHaplotypes teteHaplo, double* tabFreq, TPtrHaplotypes* tabHaplo, int tailleGeno, int nbHaploTotaux) {

	/* VARIABLES LOCALES */

	int j = 0;

	TPtrHaplotypes pHaplo = teteHaplo;

	/* DEBUT */

	j = 0;
	while (pHaplo != NULL) { // Parcours la liste d'haplotype
		tabFreq[j] = pHaplo -> freq; // On garde les frequences dans un tableau
		tabHaplo[j] = pHaplo; // On garde les haplotypes
		j++;
		pHaplo = pHaplo -> suiv;
	}
	quickSort(tabFreq, tabHaplo, 0, nbHaploTotaux-1, tailleGeno); // Appel de la fonction de trie

	/* FIN */
}

/* ******************************************************************************************************************************** */
//procedure de trie rapide pour trier les frequences
void quickSort(double* tabFreq, TPtrHaplotypes* tabHaplo, int g, int d, int tailleGeno) {

	/* VARIABLES LOCALES */

	double	indPivot;

	/* DEBUT */

	if(g < d){
		indPivot = partitionner(tabFreq, tabHaplo, g, d, tailleGeno);
		quickSort(tabFreq, tabHaplo, g, indPivot-1, tailleGeno);
		quickSort(tabFreq, tabHaplo, indPivot+1, d, tailleGeno);
	}

	/* FIN */
}

/* ******************************************************************************************************************************** */
int partitionner(double* tabFreq, TPtrHaplotypes* tabHaplo, int g, int d, int tailleGeno) {

	/* VARIABLES LOCALES */

	double valPivot;
	int i, j;

	/* DEBUT */

	valPivot = tabFreq[g];
	i = g+1;
	j = d;
	while(i <= j) {
		if(tabFreq[i] >= valPivot) { // Trie decroissant
			i++;
		} else { // v[i] attend son transfert ?droite
			while(tabFreq[j] < valPivot) { // v[j] bien placé
				j--;
			}
			if(i < j) {
				echange(tabFreq, tabHaplo, i, j, tailleGeno);
				i++;
				j--;
			}
		}
	}
	echange(tabFreq, tabHaplo, g, j, tailleGeno);
	return j;

	/* FIN */

}

/* ******************************************************************************************************************************** */
void echange(double* tabFreq, TPtrHaplotypes* tabHaplo, int i, int j, int tailleGeno) {

	/* VARIABLES LOCALES */

	TPtrHaplotypes tempHi, tempHj;
	double tempFi, tempFj;

	/* DEBUT */

	tempFi = tabFreq[i];
	tempFj = tabFreq[j];
	tabFreq[i] = tempFj;
	tabFreq[j] = tempFi;

    tempHi=tabHaplo[i];
    tempHj=tabHaplo[j];
    tabHaplo[i]=tempHj;
    tabHaplo[j]=tempHi;


	/* FIN */

}

/* ******************************************************************************************************************************** */
void ecritureListeHaploTrieFreqDecr(double* vect_freq , TPtrHaplotypes* tabHaplo, int nbHaploTotaux, int tailleGeno) {

	/* VARIABLES LOCALES */

	int i, j;
    FILE* fichier = NULL;

    /* DEBUT */

    // Supression du fichier de sauvegarde s'il existe deja
    remove ("haplotypesTriesFrequence.txt");

	// Ouverture du fichier dans lequel on va ecrire les resultats
    fichier = fopen("haplotypesTriesFrequence.txt", "a");

    if (fichier != NULL) {
		for(i=0;i<nbHaploTotaux;i++){
			for(j=0;j<tailleGeno;j++){
				fprintf(fichier, "%d", tabHaplo[i] -> haplotype[j]);
			}
			fprintf(fichier, "\t");
			fprintf(fichier, "%f \n", (vect_freq[i]*100));
		}
	}
    fclose(fichier);

    /* FIN */

}

/* ******************************************************************************************************************************** */
//fonction qui parcours la liste d'individus pour rechercher la paire ayant la plus grande proba
void rechercheProbaMaxPairHaplo(TPtrIndividus teteIndividu, TPtrGenotypes teteGenotypes, int tailleGeno) {

	/* VARIABLES LOCALES */

    TPtrIndividus pInd = teteIndividu;
	TPtrGenotypes pGenotype = teteGenotypes;
	TPtrHaplosExplicatifs pHaploE;

    /* DEBUT */

    // Supression du fichier de sauvegarde s'il existe deja
 	remove ("listeIndGenoPaireHaplo.txt");

    while (pGenotype != NULL) { // Parcours la liste des genotypes
	// Pour chaque genotype on va parcourir la liste de paireexplicatives pour trouver celle ayant le p_part le plus grand
	pHaploE = pGenotype -> teteHaplosExplicatifs;
	pGenotype -> pPaireProbaMax = pHaploE; // Initialisation du pointeur vers la paire ayant le plus grand_part avec la premiere paire
    	while (pHaploE != NULL) { // Parcours la liste des paires explicatives
			if (pHaploE -> p_part > pGenotype -> pPaireProbaMax -> p_part) { // Test si p_part plus grand alors on stock le pointeur vers cette paire
			 	pGenotype -> pPaireProbaMax = pHaploE;
			}
       		pHaploE = pHaploE -> suiv; // On avance a la paire suivante
		}
        pGenotype = pGenotype -> suiv; // On avance au genotype suivant
    }

    while (pInd != NULL) { // Parcours la liste des individus
        ecritureListeIndGenoPaireHaplo(pInd, tailleGeno); // On ecrit dans le fichier les informations de l'individu : nom, haplo1 haplo2 geno
        pInd = pInd -> suiv; // On avance a l'individu suivant
    }

    /* FIN */

}

/* ******************************************************************************************************************************** */
void ecritureListeIndGenoPaireHaplo(TPtrIndividus pInd, int tailleGeno) {

	/* VARIABLES LOCALES */

	int j = 0;
    FILE* fichier = NULL;

    /* DEBUT */

	// Ouverture du fichier dans lequel on va ecrire les r?ultats
    fichier = fopen("listeIndGenoPaireHaplo.txt", "a");

    if (fichier != NULL) {
			fprintf(fichier, "%s geno ", pInd -> nomIndividu);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", pInd -> genotype -> genotype[j]);
		}
			fprintf(fichier, "\n%s haplo1 ", pInd -> nomIndividu);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", pInd -> genotype -> pPaireProbaMax -> haplo1 -> haplotype[j]);
		}
			fprintf(fichier, "\n%s haplo2 ", pInd -> nomIndividu);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", pInd -> genotype -> pPaireProbaMax -> haplo2 -> haplotype[j]);
		}
    }
    fprintf(fichier, "\n");
    fclose(fichier);
}

/* ******************************************************************************************************************************** */
void afficher(TPtrIndividus teteIndividus, TPtrGenotypes teteGenotypes, TPtrHaplotypes teteHaplotypes, int tailleGeno, int nbIndividu) {

	/* VARIABLES LOCALES */

	int i;
    TPtrIndividus pIndividus = teteIndividus;
    TPtrGenotypes pGenotypes = teteGenotypes;
    TPtrHaplotypes pHaplotypes = teteHaplotypes;
    TPtrHaplosExplicatifs pHE;
    TPtrGenosExplicatifs pGE;

    /* DEBUT */

	printf("\n***********************************\n");
    printf("\nAffichage de la liste des individus\n");
    printf("\n-----------------------------------------\n");
    while (pIndividus != NULL) {
        printf("- nom : %s\n- genotype : ", pIndividus -> nomIndividu);
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", pIndividus -> genotype -> genotype[i]);
        }
        printf("\n-----------------------------------------\n");
        pIndividus = pIndividus -> suiv;
    }

    printf("\n***********************************\n");
    printf("\nAffichage de la liste des genotypes\n");
    printf("\n-----------------------------------------\n");
    while (pGenotypes != NULL) {
        printf("- genotype : ");
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", pGenotypes -> genotype[i]);
        }
        printf("\n- proba : %le\n", pGenotypes -> proba);
        printf("- cpt : %d\n", pGenotypes -> cpt);
        printf("- haplotypes explicatifs : \n");
        pHE = pGenotypes -> teteHaplosExplicatifs;
        while (pHE != NULL) {
        	printf("    * haplo1 : ");
        	for (i=0; i<tailleGeno ; i++) {
            	printf("%d", pHE -> haplo1 -> haplotype[i]);
        	}
        	printf("\n");
        	printf("    * haplo2 : ");
        	for (i=0; i<tailleGeno ; i++) {
            	printf("%d", pHE -> haplo2 -> haplotype[i]);
        	}
        	printf("\n");
        	pHE = pHE -> suiv;
        }
        printf("-----------------------------------------\n");
        pGenotypes = pGenotypes -> suiv;
    }

    printf("\n***********************************\n");
    printf("\nAffichage de la liste des haplotypes\n");
    printf("\n-----------------------------------------\n");
    while (pHaplotypes != NULL) {
        printf("- haplotype : ");
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", pHaplotypes -> haplotype[i]);
        }
        printf("\n- freq : %le\n", pHaplotypes -> freq);
        printf("- cpt : %d\n", pHaplotypes -> cpt);
        printf("- genotypes explicatifs : \n");
        pGE = pHaplotypes -> teteGenosExplicatifs;
        while (pGE != NULL) {
        	printf("    * geno : ");
        	for (i=0; i<tailleGeno ; i++) {
            	printf("%d", pGE -> genotype -> genotype[i]);
        	}
        	printf("\n");
        	printf("    * haplo : ");
        	for (i=0; i<tailleGeno ; i++) {
            	printf("%d", pGE -> haplotype -> haplotype[i]);
        	}
        	printf("\n");
        	pGE = pGE -> suiv;
        }
        printf("-----------------------------------------\n");
        pHaplotypes = pHaplotypes -> suiv;
    }

    /* FIN */

}