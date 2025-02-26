# Multiplication distribuée de matrices avec MPI

## Description
`smult` est un programme en C permettant d'effectuer une multiplication de matrices en parallèle en utilisant MPI (Message Passing Interface). Ce projet est conçu pour lire une matrice au format `.mtx` et effectuer des calculs de multiplication en exploitant les architectures distribuées.

## Prérequis
Avant de compiler et d'exécuter ce programme, assurez-vous d'avoir installé les bibliothèques suivantes :
- Un compilateur compatible avec MPI (ex: `mpicc` de `OpenMPI` ou `MPICH`).

## Compilation
Pour compiler le programme, utilisez la commande :
```sh
make
```
Cela génère un exécutable nommé `smult`.

## Environnement Python (Optionnel)
Un environnement Python peut être crée pour la génération des graphes, avec :
```sh
python -m venv env
source env/bin/activate
```

## Exécution
Le programme s'exécute en spécifiant un nombre de processus et un fichier de matrice :
```sh
mpirun -n 4 ./smult mtx/bcsstk03.mtx 100
```
ou :
- `4` est le nombre de processus MPI à utiliser.
- `mtx/bcsstk03.mtx` est le fichier de la matrice au format `.mtx`.
- `100` est un paramètre supplémentaire (à adapter selon l'implémentation).

## Nettoyage
Pour supprimer l'exécutable généré, utilisez :
```sh
make clean
```

## Structure du projet
- `src/main.c` : Code source principal du programme.
- `src/csr.c` : Code source de la structure de données CSR.
- `src/kernel.c` : Code source des noyaux de calcul utilisés dans le programme.
- `Makefile` : Fichier pour automatiser la compilation et l'exécution.
- `mtx/` : Dossier contenant les matrices d'entrée au format `.mtx`.
- `src/` : Dossier contenant les codes sources du programme.
- `include/` : Dossier contenant les headers du programme.
- `obj/` : Dossier contenant les fichiers objets du programme.
- `assets/` : Dossier contenant divers fichiers tels que des données et plot.
