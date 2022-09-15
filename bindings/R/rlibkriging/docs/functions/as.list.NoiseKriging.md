# `as.list.NoiseKriging`

Coerce a `NoiseKriging` Object into a List


## Description

Coerce a `NoiseKriging` Object into a List


## Usage

```r
list(list("as.list"), list("NoiseKriging"))(x, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     An object with class `"NoiseKriging"` .
`...`     |     Ignored


## Value

A list with its elements copying the content of the
  `NoiseKriging` object fields: `kernel` , `optim` ,
  `objective` , `theta` (vector of ranges),
  `sigma2` (variance), `X` , `centerX` ,
  `scaleX` , `y` , `centerY` , `scaleY` ,
  `regmodel` , `F` , `T` , `M` , `z` ,
  `beta` .


## Author

Yann Richet yann.richet@irsn.fr


## Examples

```r
f <- function(x) 1 - 1 / 2 * (sin(12 * x) / (1 + x ) + 2 * cos(7 * x) * x^5 + 0.7)
set.seed(123)
X <- as.matrix(runif(5))
y <- f(X) + 0.01*rnorm(nrow(X))
r <- NoiseKriging(y, rep(0.01^2,nrow(X)), X, kernel = "gauss")
l <- as.list(r)
cat(paste0(names(l), " =" , l, collapse = "\n"))
```

