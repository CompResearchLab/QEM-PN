# QEM-PN: Quantum Element Method for Periodic Nanostructures
The Quantum Element Method (QEM), first conceived by Ming-Cheng Cheng in [1], combines  reduced order modeling techniques — such as Proper Orthogonal Decomposition and Galerkin Projection — with domain decomposition and generic elements to form an effective simulation methodology for quantum nanostructures. This code, which accompanies [3], demonstrates the effectiveness of QEM for periodic structures by computing the wavefunctions and energy levels of nanostructures like those shown in Fig. 1. For details on the methodology of the QEM check out the papers below.
<div align="center">
<figure>
  <p align="center">
  <img src="Images/20x20Nanostructure.png" alt="Diagram" width="400">
     </p>
  <figcaption>Fig. 1: Contour of a nanostructure potential with 400 elements.</figcaption>
</figure>
</div>



## Instructions
First, add the src directory to your matlab path; this directory contains all scripts used by the  the runners in the main directory. Next read the Input.txt file, which explains the  significance and importance of each input parameter. After modifying the inputs to suit your needs, you can run the training and verification stages all together using the script RunAll.m. If everyting is set up correctly your QEM least square error (LSE) plot should resemble Fig. 2. 

<figure>
 <p align="center">
  <img src="Images/LSE.jpg" alt="Diagram" width="400">
  </p>
  <figcaption>Fig. 2: LSE of QSs 1-15, as a function of the maximum number of modes in each of the 400 elements. (b) Total number of modes (or DoF) used for the first 15 states needed in post processing.</figcaption> 
  
</figure>

<br>


<p>Run Verification.m if you want to generate a new verification nanostructure inaddition to generating the QEM results. This is usefull if RANDt=true (want to generate a random nanostructure).</p>  

To verify the model using the same nanostructure you have previously used for verification run QEMrunner.m. This allows you to change QEM parameters to optimize the LSE.

## References

[1] Cheng, M.-C. (2020). Quantum element method for quantum eigenvalue problems derived from projection-based model order reduction. AIP Advances, 10(11). https://doi.org/10.1063/5.0018698

[2] Veresko, M., & Cheng, M.-C. (2023). Quantum element method for multi-dimensional nanostructures enabled by a projection-based learning algorithm. Solid-State Electronics, 202, 108610. https://doi.org/10.1016/j.sse.2023.108610

[3] Veresko, M., Liu, Y., Hou, D., & Cheng, M.-C. (2025). Physics-Aware POD-Based Learning for Ab initio QEM-Galerkin Simulations of Periodic Nanostructures (Version 1). arXiv. https://doi.org/10.48550/ARXIV.2501.09089
