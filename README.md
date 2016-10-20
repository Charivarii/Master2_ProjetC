########################################################
# Master2_ProjetC : Inférence d'haplotypes à partir de génotypes (Algorithme EM itératif)
# Master 2 Bioinformatique - Université de Nantes - 2016-2017
# UEs "Algorithmique et Programmation avancée pour les biologistes" et "Méthodes et algorithmes pour la bioinformatique"
# Encadrement C. Sinoquet
# Réalisé par Estelle Geffard et Léa Boulongne
########################################################

# Script Generation_Genotypes_Aleatoire.c
------------------------------------------

# Script Inference_Haplotype.c
##############################

- Compilation et lancement :
----------------------------

gcc -o Inference_Haplotype Inference_Haplotype.c
./Inference_Haplotype <fichier> <nombreIndividus> <tailleGenotype>

<fichier> - nom du fichier contenant les individus et leur génotype sous le format ci-dessous : 
/ind0 011202
/ind1 112012
...
/ind10 001211

<nombreIndividu> - nombre d'individu contenu dans le fichier (dans l'exemple ci-dessus, 11)
<tailleGenotype> - taille des génotypes des individus (dans l'exemple ci-dessus, 6)
