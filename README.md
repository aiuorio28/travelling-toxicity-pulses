This is the code used to produce Fig. 19-22 in the paper:
"Travelling pulses on three spatial scales in a Klausmeier-type vegetation-autotoxicity model" by Paul Carter, Arjen Doelman, Annalisa Iuorio, and Frits Veerman.

The main code to run in order to produce such figures is "solve_k_pde_1D_fin(tend,K)". This solves the "Klausmeier autoxicity" PDE Eq. (1.1) using
the right-hand side provided in the file k_pde_1D_rhs.m and the Jacobian in the file Dk_pde_1D_rhs.m.

Figure (1) reproduces a 3D space-time plot of the pulse dynamics (in Fig. 21-22 in the paper we focus on its 2D projection).
Figure (2) shows a plot of the travelling pulse in (V,V_x,S)-phase space for case (i) (resp. (V,U,S)-phase space for case (ii)) corresponding to Fig. 19-20 in the paper.
