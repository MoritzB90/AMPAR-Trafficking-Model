Introduction
============

The maintenance of synaptic changes resulting from long-term potentiation (LTP) is essential for brain function such as memory and learning. Different LTP phases have been associated with diverse molecular processes and pathways, and the molecular underpinnings of LTP on the short, as well as long time scales, are well established. However, the principles on the intermediate time scale of 1-6 hours that mediate the early phase of LTP (E-LTP) remain elusive. We hypothesize that the interplay between specific features of postsynaptic AMPA receptor (AMPAR) trafficking is responsible for sustaining synaptic changes during this LTP phase. We test this hypothesis by formalizing a biophysical model that integrates several experimentally-motivated mechanisms. The model captures a wide range of experimental findings and predicts that synaptic changes are preserved for hours when the receptor dynamics are shaped by the interplay of structural changes of the spine in conjunction with increased trafficking from recycling endosomes and the cooperative binding of receptors. In the following sections, the basics of the model are introduced. For a more detailed description we refer to our paper (PLOS CB: `https://doi.org/10.1371/journal.pcbi.1008813 <https://doi.org/10.1371/journal.pcbi.1008813>`_).

Rate Model
----------

The rate model integrates the core processes of AMPAR trafficking in the dendritic spine by a set of differential equations with rate parameters: lateral diffusion to and from the dendritic membrane, exo- and endocytosis within the spine, and binding to and disassociation from scaffolding proteins in the postsynaptic density (PSD). We formalize the coaction of these processes by considering the number of AMPARs that freely move on the membrane of the spine denoted by :math:`U` and of AMPARs that are bound in the PSD labeled by :math:`B`. The surface area of the spine is denoted by :math:`A_{spine}`. On the synaptic membrane, :math:`U` is increased by exocytosis at a rate :math:`k_{exo}\, S_{exo}` where :math:`S_{exo}` denotes the size of an exocytic event and :math:`k_{exo}` the rate at which exocytosis events occur. :math:`U` is further increased by lateral diffusion from the dendritic to the spine membrane at a rate :math:`k_{in}`, and by the release of AMPARs from the PSD described by :math:`k_{BU}\, B`. :math:`U` is decreased by endocytosis at a rate :math:`k_{endo} \,\frac{U}{A_{spine}}` and lateral diffusion to the dendrite at a rate :math:`k_{out}\, \frac{U}{A_{spine}}`, as well as by AMPARs binding to the PSD scaffold described by :math:`k_{UB}(P-B)\,\frac{U}{A_{spine}}`, where :math:`P` denotes the number of slots formed by scaffolding proteins. Thus, the change in the number of unbound AMPARs can be formalized by

.. math::
   :label: dUdt

   \frac{dU}{dt}=k_{exo}\, S_{exo}+k_{in}+k_{BU}\, B-\left(k_{endo}+k_{out}+k_{UB}(P-B) \right)\,\frac{U}{A_{spine}}.

We do not consider direct endo- and exocytosis of bound AMPARs such that the number is only regulated by the interchange with freely moving AMPARs. 

.. math::
   :label: dBdt

   \frac{dB}{dt}=k_{UB}(P-B)\, \frac{U}{A_{spine}}-k_{BU} B.

Also, we assume that the size of an exocytosis event :math:`S_{exo}` depends on the overal spine volume :math:`V_{spine}`.

.. math::
   :label: dSexodt

   \frac{dS_{exo}}{dt}=k_{in}^{RE}-k_{out}^{RE}\, \frac{S_{exo}}{V_{spine}(t)},

Equations :eq:`dUdt`, :eq:`dBdt` and :eq:`dSexodt` describe the AMPAR-dynamics on the spine's membrane under basal condition. We derived the corresponding set of parameter values of an average spine by estimating the values from several experimental studies that are mainly conducted in hippocampal cultures and slices. For details see publication.
In addition, the model can account for cooperative receptor binding. In this case, the binding and unbinding rate :math:`k_{UB}` and :math:`k_{BU}` are exchanged by functions of :math:`P` and :math:`B`, which we here denote by :math:`k_{UB}^{coop}(P,B)` and :math:`k_{BU}^{coop}(P,B)`.

.. math::
   :label: kUBcoop

   k_{UB}^{coop}=k_{UB}\,(m(P) B^{0.8}+1),

.. math::
   :label: kBUcoop

   k_{BU}^{coop}(P,B)=k_{BU}\, (\frac{\lambda(P)}{\beta(P)+B}-0.5),

with functions :math:`m(P)`, :math:`\lambda(P)`, :math:`\beta(P)`. For details see publication. We implemented this model of AMPAR trafficking in python. The corresponding module "rate_model" containing all necessary functions can be found in the "ampar_trafficking" package. Doumentation of the code can be found under "Code documentation".

Stochastic Model
----------------

To model cooperative binding of AMPARs in the PSD we utilize a recent theoretical model capturing the essence of cooperativity given a few general assumptions [1]. In our adaptation of this model the scaffolding proteins or slots are spatially organized on a grid where each mobile receptor :math:`U` binds to a free slot with a certain rate, dependent on the number of occupied nearest-neighbor slots (Equation :eq:`kUBcoop_stochastic`). Similarly, the rate of unbinding decreases with more bound receptors :math:`B` being in the surrounding (Equation :eq:`kBUcoop_stochastic`). Note that in this model the grid size (:math:`N \times N`) is analog to :math:`P` in the rate model (Equations :eq:`dBdt` and :eq:`dUdt`).

.. math::
   :label: kUBcoop_stochastic

   \hat{k}_{UB}^{coop}=k_{UB}(\alpha \chi +1),

.. math::
   :label: kBUcoop_stochastic

   \hat{k}_{BU}^{coop}=k_{BU}(1- \chi).

Here, :math:`\chi` is the fraction of occupied nearest neighbours on the grid (i.e. :math:`n/8` with :math:`n=0,1,...,8`) and :math:`\alpha` is a measure for the cooperativity, where a big value coincides with a stronger influence of cooperativity. For details see publication. We implemented this stochastic model of AMPAR trafficking on a grid in python. The corresponding module "stochastic_model" containing all necessary functions can be found in the "ampar_trafficking" package. Doumentation of the code can be found under "Code documentation".


References
----------

[1] A. Shomar, L. Geyrhofer, N. E. Ziv, and N. Brenner. Cooperative stochastic binding and unbinding explain synaptic size dynamics and statistics. PLOS Computational Biology, 13 (7): e1005668, jul 2017. doi: `10.1371/journal.pcbi.1005668 <https://doi.org/10.1371/journal.pcbi.1005668>`_.
