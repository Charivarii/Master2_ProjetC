########################################################
# Master2_ProjetC : Inférence d'haplotypes à partir de génotypes (Algorithme EM itératif)
# Master 2 Bioinformatique - Université de Nantes - 2016-2017
# UEs "Algorithmique et Programmation avancée pour les biologistes" et "Méthodes et algorithmes pour la bioinformatique"
# Encadrement C. Sinoquet
# Réalisé par Estelle Geffard et Léa Boulongne

# Execution rapide
###################

- Execution du script bash :
----------------------------

Changer les valeurs des arguments dans le script <nbIndividu> <tailleGeno> <pcAmbiguite>
gedit exe.sh

Executer le script (si besoin changer les permissions accorder au script avec chmod 755 exe.sh)
./exe.sh


# Script Generation_Genotypes_Aleatoire.c
#########################################

- Compilation et lancement :
----------------------------

gcc -o Generation_Genotypes_Aleatoire Generation_Genotypes_Aleatoire.c
./Generation_Genotypes_Aleatoire <nombreIndividus> <tailleGenotype> <pourcentageDe2>

<nombreIndividu> - nombre d'individu contenu dans le fichier
<tailleGenotype> - taille des génotypes des individus
<pourcentageDe2> - détermine le nombre de 2 dans la séquence (nombre d'hétérozygote)
									 doit être inférieur à 0.5

# Script Inference_Haplotype.c
##############################

- Compilation et lancement :
----------------------------

gcc -o Inference_Haplotype Inference_Haplotype.c -lm
./Inference_Haplotype <fichier> <nombreIndividus> <tailleGenotype>

<fichier> - nom du fichier contenant les individus et leur génotype sous le format ci-dessous : 
/ind0 011202
/ind1 112012
...
/ind10 001211

<nombreIndividu> - nombre d'individu contenu dans le fichier (dans l'exemple ci-dessus, 11)
<tailleGenotype> - taille des génotypes des individus (dans l'exemple ci-dessus, 6)
