---
title: 'Appendix: Mendelian randomization with multiple exposures: The importance
  of thinking about time'
csl: ije.csl
header-includes:
- \usepackage{tikz}
- \usepackage{amsmath}
- \usepackage{float}
- \usepackage{cite}
- \usepackage{caption}
- \usepackage{subcaption}
- \usepackage{fixltx2e}
- \usepackage{longtable}
- \usetikzlibrary{shapes, decorations,calc}
- \usepackage{pdflscape}
- \newcommand{\blandscape}{\begin{landscape}}
- \newcommand{\elandscape}{\end{landscape}}
- \usepackage{multirow}
- \usepackage{booktabs}
- \newcolumntype{L}{<{\centering\arraybackslash}m{9cm}}
- \usepackage{float}
- \floatstyle{plaintop}
- \restylefloat{table}
- \usepackage{longtable}
output:
  html_document:
    df_print: paged
    keep_md: TRUE
    toc: TRUE
    toc_float: true
  pdf_document:
    keep_tex: yes
bibliography: V:/HomeDir/495055(J. Labrecque)/Analyses/library.bib
---
<!-- /Users/jeremylabrecque/Google Drive/epidemiology -->
<!-- V:/HomeDir/495055(J. Labrecque)/Analyses/library.bib -->



# Simulation of longitudinal datasets under different causal structures and analyses (network MR, multivariable MR and factorial MR)

The code in this repository accompanies an article about how network MR, multivariable MR and factorial perform when the longitudinal nature of data generation for Mendelian randomization studies is considered. Find more information on these methods in the following references:


Burgess S, Daniel RM, Butterworth AS, Thomspon SG, EPIC-InterAct Consortium. [Network Mendelian randomization: Using genetic variants as instrumental variables to investigate mediation in causal pathways.](https://www.ncbi.nlm.nih.gov/pubmed/25150977) Int J of Epidemiol. 2015;44(2):484-495.

Sanderson E, Davey Smith G, Windmeijer F, Bowden J. [An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings.](https://www.ncbi.nlm.nih.gov/pubmed/30535378) Int J Epidemiol. 2018;1-25. 

Burgess S, Thompson SG. [Multivariable Mendelian randomization: The use of pleiotropic genetic variants to estimate causal effects.](https://www.ncbi.nlm.nih.gov/pubmed/25632051) Am J of Epidemiol. 2015;181(4):251-260. 

Rees JMB, Foley C, Burgess S. [Factorial Mendelian randomization: using genetic variants to assess interactions.](https://www.biorxiv.org/content/10.1101/531228v1) bioarXiv. 2019.

# Appendix A - Simulations

All the code for the simulations and analyses in the text can be found at: [https://github.com/jalabrecque/MR_multiple_exposures](https://github.com/jalabrecque/MR_multiple_exposures)

## Network MR simulation

The simulation was modeled on the simulations in Burgess et al 2015 [@Burgess2015c].

\begin{align*}
a_0 &= \alpha_{G_A}*g_A + u_1 + u_2 + \epsilon_{A_0} \\
b_0 &= \beta_{G_B}*g_B + u_1 + u_3 + \epsilon_{B_0} \\
a_1 &= \alpha_{A_0}*a_0 + \alpha_{B_0}*b_0 + u_1 + u_2 + \epsilon_{A_1} \\
b_1 &= \beta_{B_0}*b_0 + \beta_{A_0}*a_0 + \beta_{A_1}*a_1 + u_1 + u_3 + \epsilon_{B_1} \\
y_1 &= \gamma_{A_0}*a_0 + \gamma_{B_0}*b_0 + \gamma_{A_1}*a_1 + \gamma_{B_1}*b_1 + u_2 + u_3 + \epsilon_Y \\
g_A, g_B &\sim \textrm{Binomial}(2, 0.3) \textrm{ independently} \\
\alpha_{G_A}, \alpha_{G_B} 
u_1, u_2, u_3 &\sim \mathcal{N}(0,1) \textrm{ independently} \\
\epsilon_{A_0}, \epsilon_{A_1}, \epsilon_{B_0}, \epsilon_{B_1}, \epsilon_{Y_1} &\sim \mathcal{N}(0,1) \textrm{ independently}
\end{align*}

where $A_t$ represents variable $A$ at time $t$ and $B_t$ represents variable $M$ at time $t$. $G_A$ and $G_B$ are modelled as biallelic genetic variants with a minor allele frequency of 0.3. $U_1$ is a confounder of $A$ and $B$, $U_2$ is a confounder of $A$ and $Y$ and $U_3$ is a confounder of $B$ and $Y$. Both $\alpha_{A_0}$ and $\beta_{B_0}$ were set to 1 so the effect of the genetic variants did not vary with time. $\alpha_{G_A}$ was set to 0.6 and $\alpha_{G_B}$ was set to 1.2. These are higher values than in the original paper. This was done because the addition of a second time point introduced additional random variation which in turn reduced the F-statistics and explained variation. The higher $\alpha_{G_A}$ and $\alpha_{G_B}$ values restored the F-statistics to around 65 and the variance explained to 1.3% in each model as in the original paper. The four gray arrows in Figure \ref{fig:timevarying_multivariable} which correspond to $\beta_{A_0}$, $\alpha_{B_0}$, $\gamma_{A_0}$ and $\gamma_{B_0}$ in the simulation were varied to determine whether they biased the network MR analysis. The effect of $A_1$ was set to 1, The code for the simulations are available at [https://github.com/jalabrecque/MR_multiple_exposures](https://github.com/jalabrecque/MR_multiple_exposures) and the user has the ability to change all parameters in the model. 


\newpage

## Multivariable MR and factorial MR simulations.

The simulation was modeled on the simulations in Sanderson et al 2019 [@Sanderson2018b]. Although only the mediation scenario was presented in the text, the confounding, collider and pleiotropy scenarios were also simulated and the results included in Appendix A. For the multivariable MR simulations, all the interactions are set to 0. 

\textbf{For all scenarios:}

\begin{align*}
g_1, g_2, ..., g_{30} &\sim \textrm{Binomial}(2, 0.3) \textrm{independently} \\
\alpha_{G_{A_{1-10}}} &= (0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50) \\
\alpha_{G_{A_{11-20}}} &= (0.50,0.45,0.40,0.35,0.30,0.25,0.20,0.15,0.10,0.05) \\
\alpha_{G_{A_{21-30}}} &= 0 \\
\alpha_{G_{B_{1-10}}} &= 0 \\
\alpha_{G_{B_{11-20}}} &= (0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50) \\
\alpha_{G_{B_{21-30}}} &= (0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50) \\
u &\sim \mathcal{N}(0,1) \textrm{ independently}
\end{align*}


\textbf{Confounding scenario}

\begin{align*}
b_0 &= \sum_{i=1}^{30}\beta_{G_{B_i}}*g_i + u\\
a_0 &= \sum_{i=1}^{30}\alpha_{G_{A_i}}*g_i + \alpha_{B_0}*b_0 + u \\
b_1 &= \beta_{b_0}*b_0 + \beta_{A_0}*a_0 + u \\
a_1 &= \alpha_{A_0}*a_0 + \alpha_{B_0}*b_0 + \alpha_{B_1}*b_1 + u \\
y_1 &= \gamma_{A_0}*a_0 + \gamma_{b_0}*b_0 + \gamma_{a_0*b_0}*a_0*b_0 + \gamma_{A_1}*a_1 + \gamma_{b_1}*b_1 + \gamma_{a_1*b_1}*a_1*b_1 + u 
\end{align*}


\textbf{Collider scenario}

\begin{align*}
a_0 &= \sum_{i=1}^{30}\alpha_{G_{A_i}}*g_i + u \\
b_0 &= \sum_{i=1}^{30}\beta_{G_{B_i}}*g_i + \beta_{A_0}*a_0 u\\
a_1 &= \alpha_{A_0}*a_0 + \alpha_{B_0}*b_0 + u \\
y_1 &= \gamma_{A_0}*a_0 + \gamma_{B_0}*b_0 + \gamma_{a_0*b_0}*a_0*b_0 + \gamma_{A_1}*a_1 + u \\
b_1 &= \beta_{B_0}*b_0 + \beta_{A_0}*a_0 + \beta_{A_1}*a_1 + \beta_Y*y + u
\end{align*}


\textbf{Pleiotropy scenario}

\begin{align*}
a_0 &= \sum_{i=1}^{30}\alpha_{G_{A_i}}*g_i + u \\
b_0 &= \sum_{i=1}^{30}\beta_{G_{B_i}}*g_i + u\\
a_1 &= \sum_{i=1}^{30}\alpha_{G_{A_i}}*g_i + \alpha_{B_0}*b_0 + u \\
b_1 &= \beta_{B_0}*b_0 + \beta_{A_0}*a_0 + u \\
y_1 &= \gamma_{A_0}*a_0 + \gamma_{B_0}*b_0 + \gamma_{a_0*b_0}*a_0*b_0 + \gamma_{A_1}*a_1 + \gamma_{B_1}*b_1 + \gamma_{a_1*b_1}*a_1*b_1 + u 
\end{align*}



\textbf{Mediation scenario}

\begin{align*}
a_0 &= \sum_{i=1}^{30}\alpha_{G_{A_i}}*g_i + u \\
b_0 &= \sum_{i=1}^{30}\beta_{G_{B_i}}*g_i + \beta_{A_0}*a_0 + u\\
a_1 &= \alpha_{A_0}*a_0 + \alpha_{B_0}*b_0 + u \\
b_1 &= \beta_{B_0}*b_0 + \beta_{A_0}*a_0 + \beta_{A_1}*a_1 + u \\
y_1 &= \gamma_{A_0}*a_0 + \gamma_{B_0}*b_0 + \gamma_{a_0*b_0}*a_0*b_0 + \gamma_{A_1}*a_1 + \gamma_{B_1}*b_1 + \gamma_{a_1*b_1}*a_1*b_1 + u 
\end{align*}












\newpage

# Appendix B - Full results

## Network MR

All results were included in the text.

\newpage

## Multivariable MR

<!--Table: MVMR SIMULATION -->

\begin{longtable}{lccccccclcc}
\caption{\label{tab:mvmr_longitudinal_simulation_full}\label{tab:mvmr_results_appendix} Estimates and bias using multivariate MR to estimate longitudinal parameters from simulated data with exposure and mediator measured at two time points in a scenario where variable B confounds variable A (confounding), where $B_1$ is caused by Y (collider), where variables A and B represent separate pleiotropic paths (Pleiotropy) and where variable B mediates the effect of variable A (Mediation).}\\
\toprule
\multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{3}{c}{Variable A} & \multicolumn{3}{c}{Variable B} \\
\cmidrule(l{3pt}r{3pt}){6-8} \cmidrule(l{3pt}r{3pt}){9-11}
Data structure & $A_0 \rightarrow Y$ & $B_0 \rightarrow Y$ & $A_0 \rightarrow B_1$ & $B_0 \rightarrow A_1$ & True & Estimate & Bias & True & Estimate & Bias\\
\midrule
Confounding & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Collider & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0\\
Pleiotropy & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Mediation & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
\addlinespace
Confounding & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Collider & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0\\
Pleiotropy & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Mediation & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
\addlinespace
Confounding & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Collider & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0\\
Pleiotropy & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Mediation & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
\addlinespace
Confounding & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Collider & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0\\
Pleiotropy & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
Mediation & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0\\
\addlinespace
Confounding & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 1.0 & -0.5\\
Collider & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 0.7 & -0.8 & 0.5 & 0.3 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 1.5 & 0.0\\
Mediation & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.0 & -0.5 & 1.5 & 1.5 & 0.0\\
\addlinespace
Confounding & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.2 & -0.2 & 1.5 & 1.2 & -0.2\\
Collider & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 0.5 & -1.0 & 0.5 & 0.3 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.2 & -0.2 & 1.5 & 1.5 & 0.0\\
Mediation & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 0.8 & -0.8 & 1.5 & 1.5 & 0.0\\
\addlinespace
Confounding & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 1.3 & -0.2 & 1.5 & 1.2 & -0.3\\
Collider & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 1.0 & -0.5 & 0.5 & 0.2 & -0.3\\
Pleiotropy & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 1.5 & 0.0 & 1.5 & 1.2 & -0.2\\
Mediation & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 1.2 & -0.2 & 1.5 & 1.2 & -0.2\\
\addlinespace
Confounding & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.2 & -0.3 & 1.5 & 1.3 & -0.2\\
Collider & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 0.8 & -0.8 & 0.5 & 0.2 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.3 & -0.2 & 1.5 & 1.3 & -0.2\\
Mediation & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.0 & -0.5 & 1.5 & 1.3 & -0.2\\
\bottomrule
\end{longtable}

\newpage

## Factorial MR

\newpage


\begin{landscape}\begin{table}[t]

\caption{\label{tab:FMR_no_intx}Full results for factorial MR when the interaction at time 0 = 0 and the interaction at time 1 = 0.5.}
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{lcccccccccclcc}
\toprule
\multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{3}{c}{Variable A} & \multicolumn{3}{c}{Variable B} & \multicolumn{3}{c}{Interaction} \\
\cmidrule(l{3pt}r{3pt}){6-8} \cmidrule(l{3pt}r{3pt}){9-11} \cmidrule(l{3pt}r{3pt}){12-14}
Data structure & $A_0 \rightarrow Y$ & $B_0 \rightarrow Y$ & $A_0 \rightarrow B_1$ & $B_0 \rightarrow A_1$ & True & Estimate & Bias & True & Estimate & Bias & True & Estimate & Bias\\
\midrule
Confounding & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Collider & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Mediation & 0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Confounding & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Collider & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Mediation & 0.0 & 0.0 & 0.5 & 0.0 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Confounding & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Collider & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Mediation & 0.0 & 0.0 & 0.0 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Confounding & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Collider & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Mediation & 0.0 & 0.0 & 0.5 & 0.5 & 1.0 & 1.0 & 0.0 & 1.0 & 1.0 & 0.0 & 0.5 & 0.5 & 0\\
Confounding & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 1.0 & -0.5 & 0.5 & 0.5 & 0\\
Collider & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 0.7 & -0.8 & 0.5 & 0.3 & -0.2 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 1.5 & 0.0 & 0.5 & 0.5 & 0\\
Mediation & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.0 & -0.5 & 1.5 & 1.5 & 0.0 & 0.5 & 0.5 & 0\\
Confounding & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.3 & -0.2 & 1.5 & 1.2 & -0.3 & 0.5 & 0.5 & 0\\
Collider & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.0 & -0.5 & 0.5 & 0.2 & -0.3 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 1.3 & -0.2 & 0.5 & 0.5 & 0\\
Mediation & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.3 & -0.2 & 1.5 & 1.3 & -0.2 & 0.5 & 0.5 & 0\\
Confounding & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 1.3 & -0.2 & 1.5 & 1.3 & -0.2 & 0.5 & 0.5 & 0\\
Collider & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 0.5 & -1.0 & 0.5 & 0.3 & -0.2 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 1.3 & -0.2 & 1.5 & 1.5 & 0.0 & 0.5 & 0.5 & 0\\
Mediation & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 0.8 & -0.7 & 1.5 & 1.5 & 0.0 & 0.5 & 0.5 & 0\\
Confounding & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.2 & -0.3 & 1.5 & 1.3 & -0.2 & 0.5 & 0.5 & 0\\
Collider & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 0.7 & -0.8 & 0.5 & 0.3 & -0.2 & 0.0 & 0.0 & 0\\
Pleiotropy & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.3 & -0.2 & 1.5 & 1.3 & -0.2 & 0.5 & 0.5 & 0\\
Mediation & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.0 & -0.5 & 1.5 & 1.3 & -0.2 & 0.5 & 0.5 & 0\\
\bottomrule
\end{tabular}}
\end{table}
\end{landscape}

\newpage


\begin{landscape}\begin{table}[t]

\caption{\label{tab:FMR_w_intx}Full results for factorial MR when the interaction at time 0 = 0.2 and the interaction at time 1 = 0.5. Only the simulations with effect of A and B at time zero are shown because interactions are not expected (but are possible) when the variables have no effects (i.e. $A_0 \rightarrow Y=B_0 \rightarrow Y=0$).}
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{lcccccccccclcc}
\toprule
\multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{3}{c}{Variable A} & \multicolumn{3}{c}{Variable B} & \multicolumn{3}{c}{Interaction} \\
\cmidrule(l{3pt}r{3pt}){6-8} \cmidrule(l{3pt}r{3pt}){9-11} \cmidrule(l{3pt}r{3pt}){12-14}
Data structure & $A_0 \rightarrow Y$ & $B_0 \rightarrow Y$ & $A_0 \rightarrow B_1$ & $B_0 \rightarrow A_1$ & True & Estimate & Bias & True & Estimate & Bias & True & Estimate & Bias\\
\midrule
Confounding & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.8 & 0.3 & 1.5 & 0.5 & -1.0 & 0.7 & 0.6 & -0.1\\
Collider & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & -0.2 & -1.7 & 0.5 & 0.4 & -0.1 & 0.2 & 0.0 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 1.5 & 0.0 & 0.7 & 0.7 & 0.0\\
Mediation & 0.5 & 0.5 & 0.0 & 0.0 & 1.5 & 0.5 & -1.0 & 1.5 & 1.8 & 0.3 & 0.7 & 0.6 & -0.1\\
Confounding & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.5 & 0.0 & 1.5 & 0.8 & -0.7 & 0.7 & 0.6 & -0.1\\
Collider & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.5 & 0.0 & 0.5 & 0.1 & -0.4 & 0.2 & 0.0 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.7 & 0.2 & 1.5 & 0.9 & -0.6 & 0.7 & 0.6 & -0.1\\
Mediation & 0.5 & 0.5 & 0.5 & 0.0 & 1.5 & 1.9 & 0.4 & 1.5 & 0.9 & -0.6 & 0.7 & 0.6 & -0.1\\
Confounding & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 0.9 & -0.6 & 1.5 & 1.9 & 0.4 & 0.7 & 0.6 & -0.1\\
Collider & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & -0.5 & -2.0 & 0.5 & 0.4 & -0.1 & 0.2 & 0.0 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & 0.9 & -0.6 & 1.5 & 1.7 & 0.2 & 0.7 & 0.6 & -0.1\\
Mediation & 0.5 & 0.5 & 0.0 & 0.5 & 1.5 & -0.1 & -1.6 & 1.5 & 1.8 & 0.3 & 0.7 & 0.6 & -0.1\\
Confounding & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 0.9 & -0.6 & 1.5 & 2.0 & 0.5 & 0.7 & 0.5 & -0.2\\
Collider & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.2 & -0.3 & 0.5 & 0.1 & -0.4 & 0.2 & 0.0 & -0.2\\
Pleiotropy & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.3 & -0.2 & 1.5 & 1.3 & -0.2 & 0.7 & 0.6 & -0.1\\
Mediation & 0.5 & 0.5 & 0.5 & 0.5 & 1.5 & 1.7 & 0.2 & 1.5 & 1.0 & -0.5 & 0.7 & 0.5 & -0.2\\
\bottomrule
\end{tabular}}
\end{table}
\end{landscape}

\newpage

# Appendix C: Sanderson et al 2019 example [@Sanderson2018b]

We assume the causal strucure below:

\begin{figure}[H]
\centering
\begin{tikzpicture}
\node[text centered] at (-4,0) (ga) {$G_{EA}$};
\node[text centered] at (-4,-1) (gb) {$G_{CA}$};
\node[text centered] at (0,0) (a) {$EA$};
\node[text centered] at (-1.5,-1) (b0) {$CA_{\textrm{pre-EA}}$};
\node[text centered] at (1.5,-1) (b1) {$CA_{\textrm{post-EA}}$};
\node[text centered] at (4,-0.5) (y) {$BMI$};
\draw[->, line width= 1] (ga) -- (a);
\draw[->, line width= 1] (gb) -- (a);
\draw[->, line width= 1] (ga) -- (b0);
\draw[->, line width= 1] (gb) -- (b0);
\draw[->, line width= 1] (b0) -- (a);
\draw[->, line width= 1] (a) -- (b1);
\draw[->, line width= 1] (b0) -- (b1);
\draw[->, line width= 1] (a) -- (y);
\draw[->, line width= 1] (b1) -- (y);
\draw[->, line width= 1] (b0) to [out=315,in=270, looseness=0.75] (y);
\end{tikzpicture}
\caption{caption here}
\label{fig:sanderson_proof}
\end{figure}

where $EA$ is educational attainment, $CA_{\textrm{pre-EA}}$ is cognitive ability before educational attainment is reached, $CA_{\textrm{post-EA}}$ is cognitive ability after educational attainment is reached, $G_{EA}$ are the genetic variants related to educational attainment, $G_{CA}$ are the genetic variants related to cognitive ability and $BMI$ is body mass index.The goal was not to numerically reproduce the example from Sanderson et al, but to determine whether this causal structure would produce biased estimates of the total effects of educational attainment and/or cognitive ability.

We found that the total effect of cognitive ability was estimated without bias but that the total effect of educational attainment is biased to the degree that $CA_{\textrm{pre-EA}}$ has a direct effect on $BMI$. Therefore, if $CA_{\textrm{pre-EA}}$ has no direct effect on $BMI$, both parameters are estimated without bias. 

The code to run this example can be found at: [https://github.com/jalabrecque/MR_multiple_exposures](https://github.com/jalabrecque/MR_multiple_exposures)

\newpage

# Appendix references

<div id="refs"></div>
