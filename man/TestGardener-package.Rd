\docType{package}
\name{TestGardener}
\alias{TestGardener}
\title{Analyses of Tests and Rating Scales using Information or Surprisal}
\description{
TestGardener is designed to permit the analysis of choice data from multiple choice tests and rating scales using information as an alternative to the usual models based on probability of choice. 

Probabability and information are related by the simple transformation "information = -log probability".  Another term for information is "surprisal."

The advantage of information methodology, often used in the engineering and physical sciences, is that measurabe, and therefore is on what is called a "ratio scale" in the social sciences.  That is, information or rurprisal has a lower limit of zero, is unbounded above, and can be added, subtracted and rescaled with a positive multiplier. 

The disadvantage of probability as a basis for representing choice is that differences near its two boundaries are on very different scales than those near 0.5, and our visual and other sensory systems, which are adapted to mangitudes, have many problems in assessing the nonlinear probability continuum.

TestGardener uses highly adaptable and computationally efficient spline basis functions to represent item characteristic curves for both probability and surprisal.  Splines bases permit as much flexibility as the task requires, and also can control the smoothness and the order of differentiation.  

The higher variability revealed by information or surprisal curves reveals many more insights into choice behavior than the usual simple curve employed in standard probability-based item response theory.

The use of information as a measure also implies a measure of inter-item covariation called mutual entropy.  Entropy a function whose value at any point is the average across surprisal curves produced by summing over curves for a given item of the product of probability and surprisal.

Graphical display is a large part of the TestGardener capacity, with extensive use of the ggplot2 and plotly packages.
}