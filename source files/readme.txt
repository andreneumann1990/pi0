title:
Estimating the proportion of true null hypotheses under dependency: 
a marginal bootstrap approach

authors: 
André Neumann(*), Taras Bodnar, and Thorsten Dickhaus
(*)Author responsible for the code (email: neumann@uni-bremen.de).

notation:
SS - Schweder-Spjotvoll
BSS - Bootstrap-SS

Section 1: Motivation
 Example 1.1:
	independence:
		The source code is contained in the folder "2018-10-30-pi0-histogram2-Independence". Only the SS results are used in the paper.
	rho = 0.5:
		The source code is contained in the folder "2018-10-24-pi0-histogram2-rho-0.5". The SS results are used in this example.

Section 2: Algorithm
 Example 2.2:
	case B = 1: 
		The source code is contained in the folder "2018-10-29-pi0-histogram2-B-1". Only the BSS results are used in the paper.
	case B = 1000:
		The source code is contained in the folder "2018-10-24-pi0-histogram2-rho-0.5". The BSS results are used in this example.
 Remark 2.3:
	The implementation using directly the distribution of the bootstrap p-values is contained in the folder "2018-11-19-pi0-histogram2-BSS-Shrinkage". The setting is the same as in Example 1.1 and 2.2, i.e., rho=0.5 and B=1000.

Section 3: Simulation Study
	The source code and the numerical results are contained in the folder "2018-10-08-pi0-simulation" and the figures are contained in the folder "2018-10-09-pi0-plots".

Section 5: Conclusions
	The comparison results of the parametric and non-parametric bootstrap are contained in the folder "2018-11-02-pi0-simulation-Comparison-P-NP-Bootstrap". In addition, some code examples for alternative pi0 estimation methods are contained in the file "conclusions-code-examples.R". These methods can be used to replace the SS estimator in the implementation. They might lead to better performance in specific situations.