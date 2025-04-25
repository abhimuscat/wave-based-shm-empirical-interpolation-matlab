# wave-based-shm-empirical-interpolation-matlab
A MATLAB toolbox for real-time digital-twin structural health monitoring: implements discrete and continuous empirical interpolation (DEIM/EIM) accelerated by Gaussian-process surrogates for wave-based damage detection.

==
                 Wave-Based SHM – Empirical Interpolation ==
==

PAPER
-----
Title: Surrogate-Accelerated Empirical Interpolation for Digital Twin Wave-Based SHM  
Authors: Abhilash Sreekumar, Linjun Zhong, U. Arasan, Dimitrios Chronopoulos  
Published: xxxxxxxxxxxxx
DOI: xxxxxxxxxxx
Preprint/GitHub: xxxxxxxxxxxxxxxxxx

DESCRIPTION
-----------
This MATLAB code accompanies the above paper. It implements two digital-twin
workflows for guided-wave structural health monitoring (SHM):
  • DEIM + Gaussian-process surrogate
  • Continuous EIM + Gaussian-process surrogate

There are two example cases (Example 1 and Example 2) for each method.

REQUIREMENTS
------------
• MATLAB R2018b or later  
• ooDACE-1.4 toolbox (included)  
• No other toolboxes required  

FOLDER CONTENTS
---------------
ooDACE-1.4/             — Kriging toolbox (dual-licensed, see LICENSE-OO)  
wavePropagation1d.m     — 1D wave solver  
f.m                     — helper function  
main_Example1_DEIM.m    — Example 1 entry point (DEIM)  
main_Example1_EIM.m     — Example 1 entry point (EIM)  
main_Example2_DEIM.m    — Example 2 entry point (DEIM)  
main_Example2_EIM.m     — Example 2 entry point (EIM)  

HOW TO RUN
----------
1. In MATLAB, cd into this folder.  
2. Type one of:
     >> main_Example1_DEIM
     >> main_Example1_EIM
     >> main_Example2_DEIM
     >> main_Example2_EIM  
   Each script will:
     • Run training simulations  
     • Build DEIM or EIM basis  
     • Fit kriging surrogates  
     • Run test cases and report errors  
     • Produce plots (convergence, signal fits, coefficients)  

TECTNICAL REPORT / PAPER LINK
-----------------------------
PDF of full paper: xxxxxxxxxxxxxxx
Source & code on GitHub: xxxxxxxxxxxxxxx

CITATION
--------
If you use this code in published work, please cite:
  Sreekumar, A., Zhong, L., Arasan, U., Chronopoulos, D. 
  “Surrogate-Accelerated Empirical Interpolation for Digital Twin Wave-Based SHM,” 
  Journal Not Specified (2025). https://doi.org/10.3390/1010000  

ooDACE TOOLBOX LICENSE (abridged)
---------------------------------
Revision: 4818; Last changed: 2009-03-23 by dgorissen  
Dual-licensed under GNU Affero GPL v3 / Commercial license  

This program is free software: you can redistribute it and/or modify  
it under the terms of the GNU Affero General Public License v3  
<http://www.gnu.org/licenses/>.  

This program is distributed WITHOUT ANY WARRANTY;  
without even the implied warranty of MERCHANTABILITY or FITNESS  
FOR A PARTICULAR PURPOSE.  

Appropriate Legal Notices must retain the “ooDACE Toolbox” text  
and homepage. When mentioning the toolbox in publications,  
reference its original publication:

  Gorissen, D., Couckuyt, I., Demeester, P., Dhaene, T., Verstichel, S.  
  “A Surrogate Modeling and Adaptive Sampling Toolbox for Computer  
   Experiments in MATLAB,” _Journal of Statistical Software_, 2010.

For full terms, see _ooDACE-1.4/LICENSE-AGPL.txt_ or  
http://www.sumowiki.intec.ugent.be/index.php/ooDACE:License_terms  

CONTACT
-------
Questions?  
  Abhilash Sreekumar — abhilash.sreekumar@kuleuven.be  
  GitHub Issues: https://github.com/xxxxxxxxxx/WaveSHM-EIM-DEIM/issues  
