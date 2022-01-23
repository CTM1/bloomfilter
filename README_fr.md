[English please!](https://github.com/CTM1/bloomfilter)
# Bloomfilter
Test de candidature pour un stage à l'Institut Pasteur. Implémentation d'un filtre de Bloom.

Cette implémentation crée un [filtre de Bloom](https://fr.wikipedia.org/wiki/Filtre_de_Bloom) 
à partir d'un fichier de format [FASTA](https://fr.wikipedia.org/wiki/FASTA_(format_de_fichier)) 
qui ne contient qu'un seul identifiant, 
et exécute ensuite un nombre donné de requêtes, chacune vérifiant 
si un kmer aléatoire est présent dans le filtre de Bloom.

Le logiciel peut prendre deux ensembles d'arguments, vérifier [Exécution](#exécution) pour plus d'informations.

1. [Compilation](#compilation)
2. [Exécution](#exécution)
3. [Documentation](#documentation)
4. [Notes](#notes)

# Compilation
Vous pouvez compiler à l'aide de `g++` en tapant `make`. D'autres options sont disponibles:

> `make clang` :
    Mêmes options que `g++`, mais pour `clang++`.

> `make lto_off` :
    Désactive la "[Link Time Optimization](https://gcc.gnu.org/wiki/LinkTimeOptimization)" de `g++`.
    Petit impact négatif sur la performance, exécutable légèrement plus petit.

> `make unoptimized` :
    Désactive les "[options d'optimisation](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)" de `g++`, `LTO` inclus.
    Grand impact sur la performance, exécutable le plus grand.

> `make clang_lto_off` :
    Même comportement pour `clang`
    Impact sur la performance négligeable, exécutable légèrement plus petit.

> `make clang_unoptimized` :
    Même comportement pour `clang`
    Grand impact sur la performance, exécutable le plus grand.

# Exécution

Le programme prendra toujours en argument un **fichier**, un **k**, et un **nombre de requêtes** (`r`). 

Par défaut, vous devez spécifier:

- La **Taille du vecteur de bits** en bits (max 2^34) (`size`).
- Un **nombre** de fonctions de **hachage** (max 64) (`nh`).

Tels que:
```bash
  #                  file        k   size  nh  r
  ./bloomfilter data/ecoli.fasta 31 456637 3 1000000
```

L'option `-p`\* prend:
- **La probabilité voulue de faux-positifs** en float (`fp`).
- **Une estimation du nombre d'items à stocker** dans le filtre de Bloom (`n`).

```bash
  #                     file        k   fp      n       r 
  ./bloomfilter -p data/ecoli.fasta 31 0.0001 4967535 1000000
```
Notes: 
- **Cette option ne respecte pas ce [README](https://github.com/yoann-dufresne/bloomtest/)!**,

Spécifiquement la taille de ce filtre de Bloom peut être au-dessus de $2^{34}$ et il peut y avoir
plus de 64 fonctions de hachage pour respecter le taux de faux positifs.

- **Comportement non-défini pour une taille de filtre supérieure à $2^{63} - 1$**. 

C'est la [taille maximum d'un vecteur en C++](https://en.cppreference.com/w/c/types/size_t) et
aussi d'un uint64_t.

J'ai estimé que $10^{13}$ éléments seraient necéssaires avec un taux de faux-positifs de $10^{-17}$, 
pour dépasser cette limite. Pour cette raison je ne l'ai pas pris en compte.


\* Le parsing des arguments n'est pas très robuste, n'importe quelle autre option aura le même comportement.

# Documentation

Vous pouvez générer la documentation avec `doxygen`, pour le faire, installez doxygen et tapez:

```bash
doxygen Doxyfile
```
Des dossiers `html` et `latex` seront créés.

Ouvrez `html/index.html` sur votre navigateur. sur GNOME:

```bash
xdg-open html/index.html
```

`make clean` supprimera les deux dossiers.

# Notes

### Petite anomalie
Le [README](https://github.com/yoann-dufresne/bloomtest/) spécifiait que le xorshift
hachera toujours 0 à 0. 

Pour tout autre chiffe, [il ne retourne pas 0](https://stackoverflow.com/questions/44753463/can-xorshift-return-zero).

Comme chaque `uint64_t kmer` est unique, changer sa valeur créerait une collision involontaire.

Dans cette implémentation, tout `k-mer` ne contenant que le nucléotide `'A'` hachera toujours à zéro,
En effet dans le filtre de Bloom, le seul bit répresentant la présence de ce kmer sera le premier.

Je n'ai pas trouvé d'astuce pour l'éviter. Néanmoins,
il n'y a pas de possibilité de collision.

---

### J'aurais aimé ajouter
- Une CLI plus robuste (le parsing n'existe presque pas).
- Parser un fichier FASTA avec plusieurs identifications, des options pour parser des identifications spécifiques avec leurs filtres de Bloom respectifs.
- Une façon de stocker des kmers avec un $k > 31$ avec le même [code Huffman](https://fr.wikipedia.org/wiki/Codage_de_Huffman).

---
### Merci !
C'était plutôt amusant.