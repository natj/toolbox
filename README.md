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

* `integ`
* `interp`
* `locate`
* `deriv`
* `smooth` and `smooth!`
* `gauss_laguerre_nw`
* `gauss_legendre_nw`
* `expi`
* `@where` (TODO)

### IO & System

* Modified `read/writeddlm` (TODO)
* `catch_NaN`
* `ReadConf`
