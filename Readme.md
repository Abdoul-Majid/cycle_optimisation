Optimisation de Cycles dans des Complexes Simpliciaux
=====================================================

Ce projet implémente une suite d'algorithmes visant à optimiser les cycles dans un complexe simplicial. L'objectif principal est de réduire la longueur totale des cycles tout en préservant leurs propriétés topologiques (homologie). Pour ce faire, le projet utilise des techniques d'optimisation locale (greedy, optimisation par chemin le plus court, recuit simulé) et de raffinage par déplacement vers le centre de gravité.

Table des Matières
------------------

*   [Présentation](#présentation)
    
*   [Fonctionnalités](#fonctionnalités)
    
*   [Structure du Projet](#structure-du-projet)
    
*   [Prérequis](#prérequis)
    
*   [Installation](#installation)
    
*   [Utilisation](#utilisation)
    
*   [Détails Algorithmiques](#détails-algorithmiques)
    
    *   [Construction du Complexe et Représentation des Simplexes](#construction-du-complexe-et-représentation-des-simplexes)
        
    *   [Optimisation des Chaînes et des Cycles](#optimisation-des-chaînes-et-des-cycles)
        
    *   [Méthodes du Plus Court Chemin et de la Bordure](#méthodes-du-plus-court-chemin-et-de-la-bordure)
        
*   [Résultats et Perspectives](#résultats-et-perspectives)
    
*   [Améliorations Futures](#améliorations-futures)
    
*   [Licence](#licence)
    

Présentation
------------

Ce projet vise à optimiser des cycles présents dans un complexe simplicial. Un cycle peut, par exemple, représenter le contour d'un trou dans un maillage 3D. L'algorithme tente de réduire la distance totale du cycle tout en respectant la classe d'homologie, afin d'éviter de "sauter" des trous ou d'altérer la topologie de la structure.

Fonctionnalités
---------------

*   **Construction de Complexes Simpliciaux :**
    
    *   Lecture de fichiers de tétraèdres (\*.ele) et de sommets (\*.node) pour créer un complexe.
        
    *   Marquage des simplexes "inside" ou "outside" en fonction d'un critère (par exemple, le quatrième sommet de la tétraèdre).
        
*   **Représentation des Simplexes :**
    
    *   Chaque simplex (noeud, arête, face, tétraèdre) est représenté par ses sommets, ses faces et ses cofaces.
        
*   **Optimisation de Cycles :**
    
    *   Algorithmes pour raccourcir et raffiner les cycles, tels que :
        
        *   **minPathOptmisationStep** : Remplace des segments du cycle par le chemin le plus court entre deux sommets.
            
        *   **approachCycleToCenter** : Ajuste le cycle pour le rapprocher de son centre de gravité en remplaçant certains segments par deux segments passant par un sommet intermédiaire.
            
        *   **simulated\_annealing** : Recuit simulé pour optimiser les cycles en échappant aux minima locaux.
            
*   **Calcul de la Bordure et Cobordure :**
    
    *   Méthodes permettant de calculer la frontière (boundary) et la cobordure d'une chaîne de simplexes.
        
*   **Visualisation :**
    
    *   Affichage interactif des complexes et des cycles optimisés à l'aide de la librairie PyVista.
        

Structure du Projet
-------------------

Plain textANTLR4BashCC#CSSCoffeeScriptCMakeDartDjangoDockerEJSErlangGitGoGraphQLGroovyHTMLJavaJavaScriptJSONJSXKotlinLaTeXLessLuaMakefileMarkdownMATLABMarkupObjective-CPerlPHPPowerShell.propertiesProtocol BuffersPythonRRubySass (Sass)Sass (Scss)SchemeSQLShellSwiftSVGTSXTypeScriptWebAssemblyYAMLXML`   pythonCopier.  ├── data  │   ├── socket.1.ele       # Fichier de tétraèdres  │   └── socket.1.node      # Fichier de sommets  ├── complex.py             # Contient les classes Simplex et Complex  ├── chain_optimization.py  # Contient la classe ChainOptimization et ses méthodes  ├── main.py                # Point d'entrée : construction du complexe, optimisation du cycle et affichage  └── README.md              # Ce fichier   `

Prérequis
---------

*   **Python 3.7+**
    
*   **Bibliothèques requises :**
    
    *   numpy
        
    *   networkx
        
    *   pyvista
        
    *   math
        
    *   itertools
        
    *   collections
        
    *   copy
        
    *   json
        
    *   random
        

Installez les bibliothèques nécessaires via pip :

Plain textANTLR4BashCC#CSSCoffeeScriptCMakeDartDjangoDockerEJSErlangGitGoGraphQLGroovyHTMLJavaJavaScriptJSONJSXKotlinLaTeXLessLuaMakefileMarkdownMATLABMarkupObjective-CPerlPHPPowerShell.propertiesProtocol BuffersPythonRRubySass (Sass)Sass (Scss)SchemeSQLShellSwiftSVGTSXTypeScriptWebAssemblyYAMLXML`   pip install numpy networkx pyvista   `

Installation
------------

1.  git clone cd
    
2.  **Vérifier que Python et les bibliothèques requises sont installés.**
    
3.  **Placer les fichiers de données (socket.1.ele et socket.1.node) dans le dossier data.**
    

Utilisation
-----------

Exécutez le script principal pour construire le complexe, optimiser un cycle et afficher le résultat :

Plain textANTLR4BashCC#CSSCoffeeScriptCMakeDartDjangoDockerEJSErlangGitGoGraphQLGroovyHTMLJavaJavaScriptJSONJSXKotlinLaTeXLessLuaMakefileMarkdownMATLABMarkupObjective-CPerlPHPPowerShell.propertiesProtocol BuffersPythonRRubySass (Sass)Sass (Scss)SchemeSQLShellSwiftSVGTSXTypeScriptWebAssemblyYAMLXML`   bashCopierpython main.py   `

Le script effectue les opérations suivantes :

*   Lit les fichiers de tétraèdres et de sommets.
    
*   Construit le complexe simplicial.
    
*   Extrait un cycle initial (liste d'indices d'arêtes) qui encercle un trou.
    
*   Applique plusieurs étapes d'optimisation pour réduire la distance totale du cycle sans modifier sa topologie.
    
*   Affiche le complexe et le cycle optimisé avec PyVista.
    

Détails Algorithmiques
----------------------

### Construction du Complexe et Représentation des Simplexes

*   **Classe Simplex :**Chaque simplex possède :
    
    *   Un index unique.
        
    *   Une liste de sommets (indices).
        
    *   Un booléen inside indiquant s'il appartient au complexe KKK.
        
    *   Des listes pour les faces (boundary) et cofaces (coboundary).
        
*   **Classe Complex :**
    
    *   Lit les fichiers de données pour créer les simplexes et les sommets.
        
    *   Calcule les relations de boundary et coboundary.
        
    *   Fournit des méthodes pour extraire la surface, afficher le complexe et calculer le centre de gravité d'un cycle.
        

### Optimisation des Chaînes et des Cycles

*   **ChainOptimization :**Regroupe plusieurs méthodes pour optimiser les cycles :
    
    *   **minPathOptmisationStep :**Découpe le cycle en segments (défini par un paramètre jump) et remplace chaque segment par le chemin le plus court entre ses extrémités (utilisation de Dijkstra).
        
    *   **approachCycleToCenter :**Réoriente le cycle en cherchant,
        
    *   **simulated\_annealing :**Méthode probabiliste d'optimisation qui permet d'éviter les minima locaux en acceptant occasionnellement des solutions moins optimales.
        

### Méthodes du Plus Court Chemin et de la Bordure

*   **find\_shortest\_path :**Construit un graphe à partir des arêtes "inside" du complexe et utilise l'algorithme de Dijkstra pour trouver le chemin de coût minimal entre deux sommets (coût = distance euclidienne).
    
*   **remove\_lonely\_edges :**Élimine les arêtes isolées, c’est-à-dire celles pour lesquelles l’un des sommets n’est pas connecté à un autre sommet par une autre arête. Cette méthode permet de garantir la cohérence du cycle.
    
*   **is\_boundary\_edge :**Détermine si une arête est considérée comme faisant partie du bord (par exemple, si elle a exactement un coface "inside").
    

Résultats et Perspectives
-------------------------

*   **Preuve de Concept sur Structures Complexes :**Les tests sur des maillages complexes (ex. : le fichier socket.1.ele) démontrent que l'algorithme parvient à réduire la longueur totale des cycles tout en préservant leur topologie.
    
*   **Possibilités d'Amélioration et d'Adaptation :**
    
    *   **En imagerie :** Optimisation des contours et segmentation d'images.
        
    *   **En réseaux :** Optimisation des chemins dans des réseaux de communication ou de transport.
        
    *   **En modélisation 3D :** Amélioration de la qualité des maillages, simplification de surfaces, etc.
        

Ces perspectives ouvrent la voie à des applications dans divers domaines, renforçant l'intérêt de combiner optimisation géométrique et topologique.

Améliorations Futures
---------------------

*   Ajustement des paramètres (ex. : pénalités dans la recherche du chemin le plus court) pour mieux respecter les contraintes topologiques.
    
*   Développement d'une interface interactive pour visualiser et manipuler les cycles optimisés.
    
*   Extension de la méthode pour traiter des complexes de dimensions supérieures ou des applications en temps réel.
    
*   Intégration d'approches d'apprentissage automatique pour guider l'optimisation.