\name{biogeom}
\encoding{UTF-8}
\alias{biogeom}
\docType{package}
\title{
Biological Geometries
}
\description{
Is used to simulate and fit biological geometries. 'biogeom' incorporates several novel universal parametric equations 
that can generate the profiles of bird eggs, flowers, linear and lanceolate leaves, seeds, starfish, 
and tree-rings (Gielis, 2003; Shi et al., 2020), three growth-rate curves representing 
the ontogenetic growth trajectories of animals and plants against time, 
and the axially symmetrical and integral forms of all these functions (Shi et al., 2017, 2021). 
The optimization method proposed by Nelder and Mead (1965) was used to estimate model parameters. 
'biogeom' includes several real data sets of the boundary coordinates of natural shapes, 
including avian eggs, fruit, lanceolate and ovate leaves, 
tree rings, seeds, and sea stars,and can be potentially applied to other natural shapes. 
'biogeom' can quantify the conspecific or interspecific similarity of natural outlines, and provides information 
with important ecological and evolutionary implications for the growth and form of living organisms. Please see 
Shi et al. (2022) for details.  
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{biogeom}
\packageIndices{biogeom}
}
\note{
We are deeply thankful to Cang Hui, Yang Li, Valeriy G. Narushin, Ülo Niinemets, Karl J. Nikas, Honghua Ruan, 
David A. Ratkowsky, Julian Schrader, Rolf Turner, Lin Wang, and Victoria Wimmer for 
their valuable help during creating this package. This work was supported by 
the National Key Research and Development Program of China (Grant No. 2021YFD02200403) and 
Simon Stevin Institute for Geometry (Antwerpen, Belguim).
}
\author{
Peijian Shi [aut, cre], Johan Gielis [aut], Brady K. Quinn [aut]

Maintainer: Peijian Shi <pjshi@njfu.edu.cn>
}
\references{
Gielis, J. (2003) A generic geometric transformation that unifies a wide range of natural 
and abstract shapes. \emph{American Journal of Botany} 90, 333\eqn{-}338. \doi{10.3732/ajb.90.3.333}

Nelder, J.A., Mead, R. (1965). A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Ratkowsky, D.A., Gielis, J. (2020) The generalized Gielis geometric 
equation and its application. \emph{Symmetry} 12, 645. \doi{10.3390/sym12040645}

Shi, P., Yu, K., Niklas, K.J., Schrader, J., Song, Y., Zhu, R., Li, Y., Wei, H., Ratkowsky, D.A. (2021)  
A general model for describing the ovate leaf shape. \emph{Symmetry}, 13, 1524. \doi{10.3390/sym13081524}
}

\keyword{ package }


