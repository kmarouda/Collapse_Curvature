## Curvature and dynamical spacetimes

arxiv version---> https://arxiv.org/abs/2209.12589

# Producing the pdf of the paper
pdflatex ./TeX/main.tex 

# Mathematica Notebooks
SphGR.nb: uses mathematica's basic functions to produce the Einstein equations in the gauge choice we make including a minimally coupled scalar field stress energy tensor and its equation of motion. The discretization tricks used in the code are worked out in the mathematica file as well.

CurvatureInvariants.nb: uses xAct xCoba to compute the Kretschmann invariant.

# Codes written in Julia

Crevolution.jl: contains all the main functions needed for the evolution.
CrConvergence.jl: does a convergence test running 3 different resolutions
CrPlotsConvergence.jl: produces the L-2 norm convergence and pointwise convergence plots
curvatureinv.jl: inputs the form of the Kretschmann invariant into the code, as computed in the relevant mathematica file
Kruns.jl: Runs different initial data and computes the Kretschmann invariant throughout the evolution. 
KrContour.jl: Creates a contour plot of the Kretschmann invariant once the Kruns.jl file has produced the data
