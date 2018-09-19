# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 16:27:51 2018

@author: s4395102
"""

\documentclass{article}

\usepackage{amssymb}
\usepackage{amsmath}

\usepackage{graphicx}

\usepackage{hyperref}

\setlength{\parindent}{0pt}
\usepackage[margin=3cm]{geometry}

\usepackage{listings}
\usepackage{xcolor}
\definecolor{mygreen}{rgb}{0,0.6,0}
\definecolor{mygray}{rgb}{0.5,0.5,0.5}
\definecolor{mymauve}{rgb}{0.58,0,0.82}
\lstset{ %
  backgroundcolor=\color{white},   % choose the background color
  basicstyle=\footnotesize,        % size of fonts used for the code
  breaklines=true,                 % automatic line breaking only at whitespace
  captionpos=b,                    % sets the caption-position to bottom
  commentstyle=\color{mygreen},    % comment style
  escapeinside={\%*}{*)},          % if you want to add LaTeX within your code
  keywordstyle=\color{blue},       % keyword style
  stringstyle=\color{mymauve},     % string literal style
  numbers=left,
  stepnumber=1,
}

\begin{document}
\title{Assignment 1 Example Template}
\date{Due Date: 5pm, 30-Aug-2018}
\author{\textbf{Name:} Travis Mitchell \\
        \textbf{Student Number:} 43214321}

\maketitle
\tableofcontents
\begin{abstract}
\noindent
    This assignment looks to address finding numerical solutions of multivariate system using Newton's method and Fourier series. It was found that Newton's method can be used to solve complex systems and Fourier transforms can be used in finding a way to numerically express a system using sine waves and cosine waves. Question 1 overlooks application of Newton's method in mutivariate system whilst Question 2 and 3 assesses Numerical implementation of Fourier series.
\end{abstract}
\newpage

\section{Question 1:}
The aim of this question is to develop a code to solve the complex, multivariate problem of modeling the equilibrium hanging position of a cable. This will be done with the use of Newton's method to solve a system of equations that describe the total energy of the hanging system.\\
\subsection{Newton's Method and Forward Difference}
Newton's method used in this assessment is as below.\\
From 
\begin{align}
    \mathbf{F}(t,u,v,w) = \begin{bmatrix} f_{1} \\
                                          f_{2} \\
                                          f_{3} \\
                                          f_{4}
                          \end{bmatrix}
\end{align}
A 4x4 Jacobian matrix consisting of partial derivatives of each f is calculated.\\
\begin{center}
{J}=\\
  $\left[\begin{array}{cccc}
  \dfrac{df_{1}}{dt} & \dfrac{df_{1}}{du} & \dfrac{df_{1}}{dv} & \dfrac{df_{1}}{dw}	\\
  \\
  \dfrac{df_{2}}{dt} & \dfrac{df_{2}}{du} & \dfrac{df_{2}}{dv} & \dfrac{df_{2}}{dw}		\\
  \\
  \dfrac{df_{3}}{dt} & \dfrac{df_{3}}{du} & \dfrac{df_{3}}{dv} & \dfrac{df_{3}}{dw}		\\
  \\
  \dfrac{df_{4}}{dt} & \dfrac{df_{4}}{du} & \dfrac{df_{4}}{dv} & \dfrac{df_{4}}{dw}		\\
  \\
  \end{array}\right]$
  \end{center}
Then, using below equation, next iteration of Newton's method is derived.
\begin{equation}
\begin{aligned}
\underline{x\textsuperscript{(n+1)}} = \underline{x\textsuperscript{(n)}}-{[J(\underline{x\textsuperscript{(n)}})]\textsuperscript{-1}}\underline{f}(\underline{x\textsuperscript{(n)}})
\end{aligned}
\end{equation}
\\
Forward difference is calculated below.
\begin{equation}
\begin{aligned}
%\underline{x\textsuperscript{(n+1)}} = \underline{x\textsuperscript{(n)}}-{[J(\underline{x\textsuperscript{(n)}})]\textsuperscript{-1}}\underline{f}(\underline{x\textsuperscript{(n)}})
f\textsubscript{x+h} = \dfrac{f(x\textsubscript{1}+h,x\textsubscript{2},x\textsubscript{3},x\textsubscript{4})-f(x\textsubscript{1},x\textsubscript{2},x\textsubscript{3}+h,x\textsubscript{4})}{h}
\end{aligned}
\end{equation}
\subsection{Step 1}
The system of equations we are looking to solve is as below.
%\begin{align}
%    \mathbf{F}(t,u,v,w) = \begin{bmatrix} \dots\\
%                                          \dots \\
%                                         \dots \\
%                                         \dots
%                         \end{bmatrix}
%end{align}
\begin{center}

{F}=
  $\left[\begin{array}{c}
  t^4+u^4-1\\
  t^2-u^2+1\\
  v^4+w^4-1\\
  v^2-w^2+1\\
  \end{array}\right]$\\
Initially, each entry in F is 1. Jacobian matrix is as below.\\
{J}=
  $\left[\begin{array}{cccc}
  4t^3 & 4u^3 & 0 & 0	\\
  \\
  2t^3 & -2u & 0 & 0	\\
  \\
  0 & 0 & 4v^3 & 4w^3	\\
  \\
  0 & 0 & 2v & -2w	\\
  \\
  \end{array}\right]$\\
Numerical value is as below.\\
{J}=
  $\left[\begin{array}{cccc}
  4 & 4 & 0 & 0	\\
  \\
  2 & -2 & 0 & 0	\\
  \\
  0 & 0 & 4 & 4	\\
  \\
  0 & 0 & 2 & -2	\\
  \\
  \end{array}\right]$\\
Applying this to Newton's method, result of iteration is as below. This process can be repeated to get accurate solution.\\
{x1}=
  $\left[\begin{array}{c}
  5/8\\
  9/8\\
  5/8\\
  9/8\\
  \end{array}\right]$\\
  \end{center}



\subsection{Step 2}
Using Q1v3.py, result has been obtained.The steps are initially using part der function to get partial derivative of each f and construct a Jacobian matrix, then using sub num df function to substitute numbers within each element in Jacobian matrix, then repeatedly using newtf functions (applies newton's method using forward difference) to converge the solution. After 100 iterations, following result was gained.\\
't': 0.000222351171308398, 'u': 1.00000000000517, 'v': 0.000222351171308398, 'w': 1.00000000000517\\
Its difference with answer from step 1 is due to numerous iterations, deriving more accurate result with each iteration unlike step 1 which was 1 iteration of newton's method.
\subsection{Step 3}
For step 3, following result was gained.By activating line 28~31 instead of 23~26, step 3 was conducted on the same code.\\
't': -0.275900847665785, 'u': 1.11680539026281, 'v': -0.985572641693287, 'w': -0.242765956626871\\
It is preferred to used numerical method rather than exact calculation because whilst computational time is advantageous in terms of the numerical method, it improves accuracy with more iterations, meaning that with sufficienct iteration numbers, accurate result can be gained.


\subsection{Step 4}
Total energy equation and graphical representation of the system is presented below.
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{springsystem.PNG}
    \caption{Graphical representation of the system}
    \label{fig:my_label}
\end{figure}
%    \[\dfrac{k_{01}\left(\sqrt{\left(x_1-x_0\right)^2+\left(y_1-y_0\right)^2}-l_{01}\right)^2}{2}+\dfrac{k_{12}\left(\sqrt{\left(x_2-x_1\right)^2+\left(y_2-y_1\right)^2}-l_{12}\right)^2}{2}+\dfrac{k_{23}\left(\sqrt{\left(y_3-y_2\right)^2+\left(x_3-x_2\right)^2}-l_{23}\right)^2}{2}+gm_2y_2+gm_1y_1 \]
\begin{equation}
\begin{aligned}
E ={} & \ gm_2y_2+gm_1y_1 \\
      & \ +\dfrac{k_{01}\left(\sqrt{\left(x_1-x_0\right)^2+\left(y_1-y_0\right)^2}-l_{01}\right)^2}{2} \\
      & \ +\dfrac{k_{12}\left(\sqrt{\left(x_2-x_1\right)^2+\left(y_2-y_1\right)^2}-l_{12}\right)^2}{2} \\
      & \ +\dfrac{k_{23}\left(\sqrt{\left(y_3-y_2\right)^2+\left(x_3-x_2\right)^2}-l_{23}\right)^2}{2}
\end{aligned}
\end{equation}
\\
From Figure 1 and total energy equation, we can see that at equilibrium, variables $x_{1},y_{1},x_{2},y_{2}$ will be minimised. 
Its partial derivatives are as following:
\begin{equation}
\begin{aligned}
\dfrac{dE}{d x_{1}} ={} & \dfrac{k_{01}\left(x_1-x_0\right)\left(\sqrt{\left(x_1-x_0\right)^2+\left(y_1-y_0\right)^2}-l_{01}\right)}{\sqrt{\left(x_1-x_0\right)^2+\left(y_1-y_0\right)^2}}-\dfrac{k_{12}\left(\sqrt{\left(x_2-x_1\right)^2+\left(y_2-y_1\right)^2}-l_{12}\right)\left(x_2-x_1\right)}{\sqrt{\left(x_2-x_1\right)^2+\left(y_2-y_1\right)^2}} \\
\dfrac{dE}{d y_{1}} ={} & \dfrac{k_{01}\left(y_1-y_0\right)\left(\sqrt{\left(y_1-y_0\right)^2+\left(x_1-x_0\right)^2}-l_{01}\right)}{\sqrt{\left(y_1-y_0\right)^2+\left(x_1-x_0\right)^2}}-\dfrac{k_{12}\left(\sqrt{\left(y_2-y_1\right)^2+\left(x_2-x_1\right)^2}-l_{12}\right)\left(y_2-y_1\right)}{\sqrt{\left(y_2-y_1\right)^2+\left(x_2-x_1\right)^2}}+gm_1 \\
\dfrac{dE}{d x_{2}} ={} & \dfrac{k_{12}\left(x_2-x_1\right)\left(\sqrt{\left(x_2-x_1\right)^2+\left(y_2-y_1\right)^2}-l_{12}\right)}{\sqrt{\left(x_2-x_1\right)^2+\left(y_2-y_1\right)^2}}-\dfrac{k_{23}\left(\sqrt{\left(x_3-x_2\right)^2+\left(y_3-y_2\right)^2}-l_{23}\right)\left(x_3-x_2\right)}{\sqrt{\left(x_3-x_2\right)^2+\left(y_3-y_2\right)^2}}\\      
\dfrac{dE}{d y_{2}} ={} & \dfrac{k_{12}\left(y_2-y_1\right)\left(\sqrt{\left(y_2-y_1\right)^2+\left(x_2-x_1\right)^2}-l_{12}\right)}{\sqrt{\left(y_2-y_1\right)^2+\left(x_2-x_1\right)^2}}-\dfrac{k_{23}\left(\sqrt{\left(y_3-y_2\right)^2+\left(x_3-x_2\right)^2}-l_{23}\right)\left(y_3-y_2\right)}{\sqrt{\left(y_3-y_2\right)^2+\left(x_3-x_2\right)^2}}+gm_2 \\
\end{aligned}
\end{equation}
\\
\subsection{Step 5}
For step 5, newton's method was used. By calculating Jacobian matrix based on partial derivative equations described above, approximate equilibrium positions for each case are described.
\subsubsection{case 1}
Jacobian matrix was calculated as below.

\begin{align}
    \mathbf{J}= \begin{bmatrix} 9.99265& -2.13718 & -0.1 & -1.69003e-6 \\
                                 -2.13718&100.01886 & -1.69003e-6& -0.06954\\
                                          -0.1&-1.69003e-6 &9.96122 &1.88976\\
                                         -1.69003e-6 &-0.06954 &1.88975 &100.02993
                          \end{bmatrix}
\end{align}
Based on the jacobian matrix and initial guess of $x_{1}$:1.1,$y_{1}$::0.9,$x_{2}$::1.7,$y_{2}$:0.8, gained positions are as below.\\
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c}
        \hline
         \textbf{variables} & \textit{values } \\
        \hline
         $x_{1}$&   0.02096 \\
         $y_{1}$&   0.00219 \\
         $x_{2}$&   2.97904\\
         $y_{2}$&   0.00219\\
        \hline
    \end{tabular}
    \caption{positions of m1 and m2 for case 1}
    \label{tab:my_label}
\end{table}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{s5case1.png}
    \caption{position of masses on case 2}
    \label{fig:my_label}
\end{figure}
\subsubsection{case 2}
Jacobian matrix was calculated as below.

\begin{align}
    \mathbf{J}= \begin{bmatrix} 199.99996	& -0.024499 & -100.00000 & 0 \\
                                 -0.024499	&160		& 0			& -80\\
                                 -100		&0 			&200 		&	0.024499\\
                                    0 &-80 &0.024499 &160
                          \end{bmatrix}
\end{align}
Based on the jacobian matrix and initial guess of $x_{1}$:1.1,$y_{1}$::0.9,$x_{2}$::1.7,$y_{2}$:0.8, gained positions are as below.\\
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c}
        \hline
         \textbf{variables} & \textit{values } \\
        \hline
         $x_{1}$&   1 \\
         $y_{1}$&   1 \\
         $x_{2}$&   2\\
         $y_{2}$&   1\\
        \hline
    \end{tabular}
    \caption{positions of m1 and m2 for case 2}
    \label{tab:my_label}
\end{table}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{s5case2.png}
    \caption{position of masses on case 2}
    \label{fig:my_label}
\end{figure}

\subsection{Step 6}
With given test values of 37.2, 4.98, -67.6, 25.1 as f, at 
%$x_{1}$:1.1,$y_{1}$:0.9,$x_{2}$::1.7,$y_{2}$:0.8 as positions, solution was checked using code named Q1_p2_v2.py, applying symbolic calculation (sympy) module, established forawrd difference method. \\
$x_{1}$:1.1,$y_{1}$::0.9,$x_{2}$::1.7,$y_{2}$:0.8, solution was checked using code Q1p2v2.py, applying symbolic calculation (sympy) module, established forawrd difference method. The result is as below.\\
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c|c|c|c}
        \hline
         \textbf{Harmonics} &  \textit{Test f}&\textit{Method 1} & \textit{Method 2} & \textit{Method 3} \\ 
        \hline
         $x_{1}$&   37.2&   37.25345 &37.24416 & 37.25345 \\
         $y_{1}$&   4.98&   4.98094 &4.98147 & 4.98094 \\
         $x_{2}$&   -67.6&   -67.55392 &-67.56269 & -67.55392\\
         $y_{2}$&   25.1&   25.06764 &25.06712 & 25.06764\\
        \hline
    \end{tabular}
    \caption{Comparison of results with test f}
    \label{tab:my_label}
\end{table}
Appendix 3 has shown the result lines on the program for step 6.
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{s6.png}
    \caption{position of masses on case 2}
    \label{fig:my_label}
\end{figure}
\newpage
\subsection{Step 7}
\subsubsection{case 1}
x values and y values gained using code is provided below.
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c|c|c}
        \hline
         \textit{$x_{0}$,$y_{0}$} &  \textit{$x_{1}$,$y_{1}$}&\textit{$x_{2}$,$y_{2}$} & \textit{$x_{3}$,$y_{3}$}\\ 
        \hline
         0.0& 0.0209603& 2.979045& 3.0 \\
         1.0& 0.002195& 0.0021984& 1.0 \\
        \hline
    \end{tabular}
    \caption{Reuslt gained using code for case 1}
    \label{tab:my_label}
\end{table}
\subsubsection{case 2}
x values and y values gained using code is provided below.
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c|c|c}
        \hline
         \textit{$x_{0}$,$y_{0}$} &  \textit{$x_{1}$,$y_{1}$}&\textit{$x_{2}$,$y_{2}$} & \textit{$x_{3}$,$y_{3}$}\\ 
        \hline
         0.0& 0.99999994& 2.0000000& 3.0 \\
         1.0& 0.998775& 0.998775& 1.0 \\
        \hline
    \end{tabular}
    \caption{Reuslt gained using code for case 2}
    \label{tab:my_label}
\end{table}
Comparing these results to analytical solution in step 5 which used symbollic calculations, it can be seen that the code is accurately working.
\subsection{Step 8}
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c|c|c}
        \hline
         \textit{$x_{0}$,$y_{0}$} &  \textit{$x_{1}$,$y_{1}$}&\textit{$x_{2}$,$y_{2}$} & \textit{$x_{3}$,$y_{3}$}\\ 
        \hline
         0.0& 1.0901590& 2.14688108& 3.0 \\
         1.0& 0.408660& 0.2202267& 1.0 \\
        \hline
    \end{tabular}
    \caption{Reuslt gained using code for values in Step 6}
    \label{tab:my_label}
\end{table}
Plotting this, following graph is obtained.
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{s8.png}
    \caption{position of the masses on Step 6}
    \label{fig:my_label}
\end{figure}
\subsection{Step 9}
Generally speaking about the code, I can say the code is accurate, especially from step 5 which showed accurate result using 3 methods. However, the problem with codes I used for Question 1 is that they are optimised to calculate 4 variables. To calculate results using 20 variables (20 masses), possible ways include setting additional functions with f numbers, improving function responsible for executing newton's method based on method 2 $(sympy)$ I have implemented. Since $part_der$ function always takes 1 function anyways, by partially deriving each function (f1~f20) and applying them to newton's method function, the code can be improved.
\newpage
\section{Question 2:}
\begin{itemize}
    \item What is the aim of this question?
    \item What are the numerical methods and a basic introduction?
\end{itemize}
Aim of the question is to apply Fourier series on a given set of data for effective data analysis.
\newpage


\subsection{Part (a)}
The real components can be seen in Fig.~\ref{fig:2_a}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{Q2parl.png}
    \caption{Real solution of series}
    \label{fig:2_a}
\end{figure}\\
The imaginary components can be seen in Fig.~\ref{fig:2_b}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{Q2paimg.png}
    \caption{Imaginary solution of series}
    \label{fig:2_b}
\end{figure}\\
The real components can be seen in Fig.~\ref{fig:2_c}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{Q2parlvsimg.png}
    \caption{Real vs Imaginary solution of series}
    \label{fig:2_c}
\end{figure}\\
Polar representation is also available in Fig.~\ref{fig:2_d}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{Q2papolar.png}
    \caption{Polar representation}
    \label{fig:2_d}
\end{figure}\\
\\
\newpage
\subsection{Part (b)}
Below is the graph showing frequency and its corresponding magnitude.
Note that frequency is n/40 (n is each entry index) since the signal is 40 seconds long.\\
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{secb.png}
    \caption{Frequency vs magnitude}
    \label{fig:2_d}
\end{figure}
\begin{table}[h!]
    \centering
    \begin{tabular}{c|c}
        \hline
         \textbf{Frequency} & \textit{Magnitude} \\
        \hline
         & 485.57568  \\
         &\\
         &369.28656\\
         &\\
        \hline
    \end{tabular}
    \caption{Frequencies and magnitudes of outliers}
    \label{tab:my_label}
\end{table}
\newpage
\subsection{Part (c)}
Compressed signal from part b is shown below. Comparing this with original signal, we can see that through neglecting proportion of data, it has smoothened the series, enabling easier visual assessment.
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{secc.png}
    \caption{Frequency vs magnitude of compressed signal}
    \label{fig:2_d}
\end{figure}
\newpage
\section{Question 3:}



\newpage
\begin{thebibliography}{00}
\bibitem{Answers} Mitchell, T., \textit{Development of MECH3750 solutions keys}, myDocuments, \textit{In Press}.


\newpage
\section{Appendix}
\subsection{Question 1 Codes}
You should describe what each of your functions do! \\

Function: \textbf{F(xvec)} \\
This function defines the system of equations to be solved....
\subsection{Appendix 3}
\begin{figure}[h!]
    \centering
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio]{s6_result.PNG}
    \label{fig:my_label}
\end{figure}
\begin{lstlisting}[language=Python]
def F(xvec):
    This is how we can include code in LaTeX

def includeCodeFromFile(file):
    Or if we want to link this to the file we have written we can use
    \lstinputlisting[language=Python]{./Code/AllMyAssignmentSolutions.py}
\end{lstlisting}

\end{thebibliography}
\end{document}
