# Studies Of The Foulis Randall Product
This repository contains numerous examples related to studies of the Foulis-Randall product.

See each of the accompanying Jupyter notebooks for implementations of the examples:
* __BellCHSHContextualityExperiment.ipynb__ : Here we conduct a standard Bell-CHSH experiment involving probabilistic outcomes detailed in our paper (see Figure 2). The results present the Foulis-Randall product of the system, and carry out the specification detailed in Section 4.1 to determine that the system is in fact non-classical.

* __CalculatingCardinalitiesOfFSets.ipynb__ : Here we implement the algorithm 'CombinationAssist' and part of the algorithm 'FoulisRandallProducts', to examine how (with the increase of the number of parties involved in the experiment), the  𝔉  set and  𝔽  set grow in size (see Section 5 of our paper).

* __FoulisRandallProductForThreeParties.ipynb)__ : Here we consider an experiment involving three parties, Alice, Bob, and Charlie (Subsection 4.2.2 of our paper). We calculate all possible Foulis-Randall products on the scenario. Our results reveal that all three possible products on the scenario are exactly equivalent.

* __NonAssociativityOfFoulisRandallProduct.ipynb__ : Here we demonstrate non-associativity of the example presented in Subsection 4.2.3 of our paper. Highlighted in red in the results is an edge that is only found in one of the three possible Foulis-Randall products on the scenario.
