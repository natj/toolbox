# toolbox.jl
### Small tools and snippets that I happen to use with julia



## Installation

Installation is easy, just type:
```julia
Pkg.clone("https://github.com/natj/toolbox.jl.git")
```
and to load it put on top of your source code
```julia
using toolbox
```

## Info
### Mathematical & numerical stuff

Integrate vector using weighted mean of two parabolas on each side
* `integ(x, y)`

Cumulative integral using `integ`. Also possible to start from zero using different extrapolations (default is `:lin`)
* `cuminteg(x, y[; extrapolate_zero=:lin, :quad, :plaw, :none])`

First derivative using weighted parabolas
* `deriv(x, y)`

Cubic or linear interpolation and linear extrapolation
* `interp(x, y, val[; method=:cubic, :lin])`

Binary search for arrays
* `locate(x, val)`

Smooth vector using Gaussian kernel `N` times. Also possible to define offset so that only `x[offs+1:end-offs-1]` is smoothed ensuring proper boundary conditions.
* `smooth(x[, N=1; offs=3])` and `smooth!(x[, N=1; offs=3])`

Smooth vector with B-splines (DeBoor's algorithm)
* `şmooth_spline(x, [weights], smoothfactor)`

Smooth vector that has a powerlaw behaviour by not smoothing the vector values directly but by smoothing the ratio between it and some reference vector `ref`.
* `şmooth_plaw(x, ref, smoothfactor/N[; offs=3, method=:kernel/:spline])`

Nodes and weights for `N` point Gaussian quadrature. Returns tuple of `(nodes, weights)`
* `gauss_laguerre_nw(N)`
* `gauss_legendre_nw(N)`

Exponential integral E_N(x) = int_1^infty e^(-x t) dt/ t^N for positive arguments
* `expi(N, x)`

Indexes of an array fulfilling given criteria `expr`
* `@where expr` (TODO)

### IO & System

Read parameters from config file. Searching for line where `param=val` and returns `val`.
* `ReadConf(file, params...)`

Throw an error if x has `NaN`s in it.
* `catch_NaN(x)`


* Modified `read/writeddlm` (TODO)


