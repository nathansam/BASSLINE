# SMLN
Objective Bayesian survival analysis using shape mixtures of log-normal distributions.

A project to convert code written by Vallejos & Steel[1] into a library using C++ to speed-up functions.  

### Build states


![Windows build state](https://github.com/nathansam/SMLN/workflows/Win-build/badge.svg)
![macOS build state](https://github.com/nathansam/SMLN/workflows/macOS-build/badge.svg)
![Linux Build state](https://github.com/nathansam/SMLN/workflows/Linux-build/badge.svg)

### Installation

The recommended way to install SMLN is via pak which will handle intalling all dependencies: 
```{R}
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
pak::pkg_install("nathansam/SMLN")
```

Alternatively, SMLN can be install via the devtools package

```{R}
devtools::install_github("nathansam/SMLN")
```

### References 
- [1] <a href="http://dx.doi.org/10.1080/01621459.2014.923316">Vallejos & Steel (2015). Journal of the American Statistical Association. </a>
