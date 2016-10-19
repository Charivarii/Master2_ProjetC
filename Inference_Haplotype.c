#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    //Declaration
    FILE* fp = NULL;
    char c;
    int i = 0;
    int j = 0;
    char tab[20][20];

    //Programme
    fp = fopen("geno.txt", "r");
 
    if (fp == NULL) {
        printf("Le fichier geno.txt n'a pas pu être ouvert\n");
        return EXIT_FAILURE;
    } else {
        printf("Le fichier geno.txt existe\n");
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

