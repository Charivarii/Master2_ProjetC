#!/bin/bash

nbIndividu=10
tailleGeno=10
pcAmbiguite=30

echo "Compilation du fichier de fonctions de Generation Aleatoire de genotypes en fichier .o"
/usr/bin/gcc -c Fonctions_Generation_Genotypes_Aleatoire.c
echo ""

echo "Creation de l'executable de Generation_Genotypes_Aleatoire.c"
/usr/bin/gcc -o Generation_Genotypes_Aleatoire Generation_Genotypes_Aleatoire.c Fonctions_Generation_Genotypes_Aleatoire.o
echo ""

echo "Execution de Generation_Genotypes_Aleatoire $nbIndividu $tailleGeno $pcAmbiguite"
./Generation_Genotypes_Aleatoire $nbIndividu $tailleGeno $pcAmbiguite
echo ""

echo "Compilation du fichier de fonctions de Inference Haplotype en fichier .o"
/usr/bin/gcc -c Fonctions_Inference_Haplotype.c
echo ""

echo "Creation de l'executable de Inference_Haplotype.c"
/usr/bin/gcc -o Inference_Haplotype Inference_Haplotype.c Fonctions_Inference_Haplotype.o -lm -g
echo ""

echo "Inference_Haplotype genotypeAleatoires.txt $nbIndividu $tailleGeno"
./Inference_Haplotype genotypesAleatoires.txt $nbIndividu $tailleGeno
echo ""


