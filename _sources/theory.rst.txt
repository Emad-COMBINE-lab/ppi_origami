Theory
======

Preparing datasets for the purpose of training models that predict protein-protein interactions is a deceptively fraught
process.

Several studies have outlined that insufficiently controlling for protein identity between cross-validation splits can
lead to serious over-fitting of PPI prediction methods [:ref:`1-3 <References>`].

PPI Origami uses the notation from Park and Marcotte, which defines three types of PPI cross-validation datasets [:ref:`1 <References>`]:

- **C3** - Proteins that constitute interactions in one split (*i.e.*: training, validation, or test) are not to be found in any other split.
- **C2** - No more than one protein in a given interaction may be found in another split.
- **C1** - No restriction on protein split membership. Interactions are randomly assigned to a split.

In addition to this, PPI Origami ensures that **C3 datasets** meet the following two criteria oultined in the INTREPPPID
manuscript.

First, let's begin by defining :math:`\{P_{\textsf{Tr}}, P_{\textsf{Te}}, P_{\textsf{V}}\}`, which are the set of proteins
present in the interactions found in the Training, Testing, and Validation split, respectively.

Further, let's define :math:`\mathcal{P}` as the collection of protein sets :math:`\{P_{\textsf{Tr}}, P_{\textsf{Te}}, P_{\textsf{V}}\}`

**Criterion 1 - Distinct Protein Identity** The protein sets :math:`\{P_{\textsf{Tr}}, P_{\textsf{Te}}, P_{\textsf{V}}\}` must be mutually disjoint:

.. math::

   \forall Q, R \in \mathcal{P}, Q \cap R = \varnothing \textsf{ if } Q \neq R.

**Criterion 2 - Distinct Sequence Identity** 

.. math::

   \forall Q, R \in \mathcal{P}, \;\;\;\; \forall q \in Q, \;\;\;\; \forall r \in R,\;\;\;\; f(q,r) \leq 90\% \;\;\;\; \textsf{ if } \;\;\;\; Q \neq R,

where :math:`f` is some sequence similarity metric. We use UniRef cluster membership for sequence similarity.

References
----------

1. Park, Yungki and Edward M. Marcotte. “`A flaw in the typical evaluation scheme for pair-input computational predictions <https://doi.org/10.1038/nmeth.2259>`_.” *Nature methods* 9 (2012): 1134 - 1136.
2. Hamp, Tobias and Burkhard Rost. “`More challenges for machine-learning protein interactions <https://doi.org/10.1093/bioinformatics%2Fbtu857>`_.” *Bioinformatics* 31 10 (2015): 1521-5 .
3. Bernett, Judith, David B. Blumenthal and Markus List. “`Cracking the black box of deep sequence-based protein-protein interaction prediction <https://doi.org/10.1101/2023.01.18.524543>`_.” *bioRxiv* (2023): n. pag.