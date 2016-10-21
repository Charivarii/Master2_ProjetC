#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 *##### STRUCTURES #####
 */
typedef struct TypeIndividu* TypeFichier;

typedef struct TypeIndividu {
    char* nomIndividu;
    int* genotype;
    TypeFichier suiv;
} TypeIndividu;

/*
 *##### FONCTIONS #####
 */
void recupererParametres(int nbArgument, char* listeArguments[], char** adr_fichier, char** adr_nbIndividu, char** adr_tailleGeno) {
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
}

void recupererDonneesFichier(TypeFichier tete, int tailleGeno, int nbIndividu, char* fichier) {
    // Variables globales
    int i = 0;
    int j = 0;
    int estNomIndividu = 1;
    FILE* fp = NULL;
    char c[1];
    TypeFichier ptr;

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
        ptr -> suiv = NULL;

        do {
            c[0] = getc(fp);
            if (c[0] == '\n' && i != (nbIndividu - 1)) {
                ptr -> suiv = (TypeFichier)malloc(sizeof(TypeIndividu));
                ptr = ptr -> suiv;
                ptr -> nomIndividu = (char*)malloc(sizeof(char)*20);
                ptr -> genotype = (int*)malloc(tailleGeno*sizeof(int));
                ptr -> suiv = NULL;
                i++;
                estNomIndividu = 1;
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
}

void afficherListe(TypeFichier tete, int tailleGeno) {
    TypeFichier p;
    int i;
    p = tete;
    while (p != NULL) {
        printf("nom : %s\ngenotype : ", p->nomIndividu);
        for (i=0; i<tailleGeno ; i++) {
            printf("%d", p-> genotype[i]);
        }
        printf("\n*************\n");
        p=p->suiv;
    }
}

/*
 *##### MAIN #####
 */
int main(int argc, char* argv[]) {
    //VARIABLE
    char* fichier;
    char* nbIndividu;
    char* tailleGeno;
    TypeFichier ptrFichier;

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

    return 0;
}