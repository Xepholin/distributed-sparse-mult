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

## Exécution
Le programme s'exécute en spécifiant un nombre de processus et un fichier de matrice :
```sh
mpirun -n 4 ./smult data/bcsstk03.mtx 100
```
ou :
- `4` est le nombre de processus MPI à utiliser.
- `data/bcsstk03.mtx` est le fichier de la matrice au format `.mtx`.
- `100` est un paramètre supplémentaire (à adapter selon l'implémentation).

## Nettoyage
Pour supprimer l'exécutable généré, utilisez :
```sh
make clean
```

## Structure du projet
- `main.c` : Code source principal du programme.
- `Makefile` : Fichier pour automatiser la compilation et l'exécution.
- `data/` : Dossier contenant les matrices d'entrée au format `.mtx`.
