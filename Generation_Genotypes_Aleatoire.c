#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void recupererParametres(int nbArgument, char* listeArguments[], char** adr_nbIndividu, char** adr_tailleGeno , char** adr_pc_ambigue) {

    //DEBUT
    if (nbArgument != 4) {
        printf("Au lancement du programme, vous devez rentrer 3 paramètres : \n");
        printf("- Le nombre d'individu souhaité\n");
        printf("- La longueur du génotype voulu\n");   
	printf("- Le pourcentage d'ambigue max : (0.5 max) \n");
        printf("Merci\n");
        exit(1);
    } else {
        *adr_nbIndividu = listeArguments[1];
        *adr_tailleGeno = listeArguments[2];
        *adr_pc_ambigue = listeArguments[3];
    }
    //FIN
}

void generation_haplotype(int nbInd, int tailleGeno, int nbAmbigue, int* haplo1, int* haplo2) {
    int i, j, k, l, nbDeux;
    int* tabTemp;
    for (j = 0 ; j < tailleGeno ; j++) {
        haplo1[j] = rand()%2;
        haplo2[j] = haplo1[j];
    }

    nbDeux = rand()%(nbAmbigue + 1);//nombre de loci ambigue

    tabTemp = (int*)malloc(nbDeux*sizeof(int));

    //Vérification des allocations de mémoire dynamique
    if(tabTemp == NULL) {
        printf("Problème d'allocation de mémoire dynamique");
        exit(1);
    }

    for (k = 0 ; k < nbDeux ; k++) {//tirage aléatoire de nbChangeLoci position dans haplo2
        tabTemp[k] = rand()%(tailleGeno + 1);
    }

    for (k = 0 ; k < nbDeux ; k++) {//on va modifier les positions  variables
        for (l = 0 ; l < k ; l++) {
            if (tabTemp[k] == tabTemp[l]) {
                break; // on quitte cette boucle
            }
        }
        if (haplo2[tabTemp[k]] == 0) {
            haplo2[tabTemp[k]] = 1;
        } else {
            haplo2[tabTemp[k]] = 0;
        }
    }

    // On n'a plus besoin de la mémoire, on la libère
    free(tabTemp);
}


void generation_genotype(int tailleGeno, int* haplo1, int* haplo2, int* geno) {
    int i, j;
    for (i = 0 ; i < tailleGeno ; i++) {
        if (haplo2[i] != haplo1[i]) {
            geno[i] = 2;
        } else if (haplo2[i] == 1) {
            geno[i] = 1;
        } else {
            geno[i] = 0;
        }
    }
}


void ecriture_fichier_geno_haplo(int i, int tailleGeno, int* geno, int* haplo1, int* haplo2) {
    FILE* fichier = NULL;
    int j;

    fichier = fopen("geno_haplo.txt", "a");

    if (fichier != NULL) {
        fprintf(fichier, "/ind%d geno ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", geno[j]);
        }
        fprintf(fichier, "\n");
        fprintf(fichier, "/ind%d haplo1 ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", haplo1[j]);
        }
        fprintf(fichier, "\n");
        fprintf(fichier, "/ind%d haplo2 ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", haplo2[j]);
        }
        fprintf(fichier, "\n");
        fclose(fichier);
    }
}


void ecriture_fichier_geno(int i, int tailleGeno, int* geno) {
    FILE* fichier = NULL;
    int j;

    fichier = fopen("geno.txt", "a");

    if (fichier != NULL) {
        fprintf(fichier, "/ind%d ", i);
        for (j = 0 ; j < tailleGeno ; j++) {
            fprintf(fichier, "%d", geno[j]);
        }
        fprintf(fichier, "\n");
        fclose(fichier);
    }
}


int main(int argc, char* argv[]) {
    //déclaration des variables
    int i, j, nbAmbigue;
    char*  nbIndC;
    char* tailleGenoC;
    char* pcAmbigueC;
    int  nbInd, tailleGeno,pcAmbigue;
    int* haplo1;
    int* haplo2;
    int* geno;

    // Récupération des paramètres
    recupererParametres(argc, argv, &nbIndC, &tailleGenoC, &pcAmbigueC);

    //passage de char en int
    nbInd=atoi(nbIndC);
    tailleGeno=atoi(tailleGenoC);
    pcAmbigue=atoi(pcAmbigueC);
	
printf("%d nb ind \n",nbInd);
printf("%d taille geno \n",tailleGeno);
printf("%d ambigue \n",pcAmbigue);

    //allocation dynamique de mémoire
    haplo1 =(int*)malloc(tailleGeno*sizeof(int));
    haplo2 =(int*)malloc(tailleGeno*sizeof(int));
    geno =(int*)malloc(tailleGeno*sizeof(int));

    //Vérification des allocations de mémoire dynamique
    if(haplo1 == NULL || haplo2 == NULL  || geno == NULL ){
        printf("Problème d'allocation de mémoire dynamique");
        exit(1);
    }

    //supression des fichiers de sauvegarde s'ils existent deja
    remove ("geno.txt");
    remove ("geno_haplo.txt");

    if (pcAmbigue > 50) {//verifier ambigue pas trop grand
        printf("Pourcentage d'ambigue choisit trop grand");
        exit(1);
    } else {
        nbAmbigue = (pcAmbigue * tailleGeno)/100;
printf("%d nbambigue \n",nbAmbigue);
    }
    for (i = 0 ; i < nbInd ; i++) {
        generation_haplotype(nbInd, tailleGeno, nbAmbigue, haplo1, haplo2);
        generation_genotype(tailleGeno, haplo1, haplo2, geno);
        ecriture_fichier_geno_haplo(i, tailleGeno, geno, haplo1, haplo2);
        ecriture_fichier_geno(i, tailleGeno, geno);
    }

    printf("Generation des haplotypes et genotypes effectues\n");

    // On n'a plus besoin de la mémoire, on la libère
    free(haplo1);
    free(haplo2);
    free(geno);

    return 0;
}
