#include <stdio.h>
#include <stdlib.h>

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

//***************************************************************
int main(int argc, char* argv[]) {
    //VARIABLE
    FILE* fp = NULL;
    char c;
    int i = 0;
    int j = 0;
    char tab[20][20];
    char* fichier;
    char* nbIndividu;
    char* tailleGeno;

    //DEBUT

    // Récupération des paramètres
    recupererParametres(argc, argv, &fichier, &nbIndividu, &tailleGeno);


    fp = fopen(fichier, "r");
 
    if (fp == NULL) {
        printf("Le fichier %s n'a pas pu être ouvert\n", fichier);
        exit(1);
    } else {
        printf("Le fichier %s existe\n", fichier);
        do {
            c = getc(fp);
            if (c == '\n') {
                j = -1;
                i++;
            }
            
            tab[i][j] = c;
            j++;
        } while (c != EOF);
    }

    for (i=0 ; i<20 ; i++) {
        for (j=0 ; j<20 ; j++) {
            printf("i : %d - j : %d - c : %c\n", i, j, tab[i][j]);
        }
    }

    return 0;
}