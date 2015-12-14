Overview
--------

Tool box for experimenting with TV types.

Requirements
------------

1. multimap
2. numbers
3. TODO

Maxima
------

TODO

```
(n k) = 2^n / sqrt(1/2*n*pi) * exp(-(k-(n/2))^2 / n/2)
```

```c++
////////////////////
// Using binomial //
////////////////////

// P(x+X successes in n+k trials | x successes in n trials)
```

```maxima
P(n, x, k, X) := ((n+1)*binomial(k, X)*binomial(n, x))
              / ((k+n+1)*binomial(k+n, X+x));
```
```c++
// Using p=x/n and s=(x+X)/(n+k)
//
// p = x/n
// x = p*n
// s = (x+X) / (n+k) = (p*n+X) / (n+k)
// s*(n+k) = p*n+X
// X = s*(n+k) - p*n
// X+x = s*(n+k)
// s(X) = (p*n+X) / (n+k)
// s(X+1) = (p*n+X+1) / (n+k)
// s(X) - s(X+1) = 1/(n+k)
```

```maxima
P_prob(n, p, k, s) := ((n+1)*binomial(k, s*(n+k) - p*n)*binomial(n, p*n))
                        / ((k+n+1)*binomial(k+n, s*(n+k)));

pdf_k(n, p, k, s) := P_prob(n, p, k, s) * (n+k);
pdf(n, p, s) = lim k->inf pdf_k(n, p, k, s)

limit(pdf(n, p, k, s), k, inf);
```

```c++
//////////////////////////////////
// Using binomial approximation //
//////////////////////////////////

// Shitty: when n is high and k is low, the approximation is really bad
```
```maxima
binomial_approx(n, k) := 2^n / sqrt(1/2*n*%pi) * %e^(-(k-(n/2))^2/(n/2));
P_approx(n, x, k, X) := ((n+1)*binomial_approx(k, X)*binomial_approx(n, x))
                  / ((k+n+1)*binomial_approx(k+n, X+x));
```

```c++
/////////////////////////
// Using beta function //
/////////////////////////
```

```maxima
binomial_beta(n, k) := 1 / ((n+1) * beta(n-k+1, k+1));
P_beta(n, x, k, X) := ((n+1)*binomial_beta(k, X)*binomial_beta(n, x))
                  / ((k+n+1)*binomial_beta(k+n, X+x));

P_beta_prob(n, p, k, s) := ((n+1)*binomial_beta(k, s*(n+k) - p*n)*binomial_beta(n, p*n))
                        / ((k+n+1)*binomial_beta(k+n, s*(n+k)));

                   beta(- X - x + n + k + 1, X + x + 1)
P_beta = ---------------------------------------------------------
         (k + 1) beta(- x + n + 1, x + 1) beta(- X + k + 1, X + 1)

                                  beta(- (n + k) s + n + k + 1, (n + k) s + 1)
P_beta_prob = -----------------------------------------------------------------------------------------
              (k + 1) beta(- n p + n + 1, n p + 1) beta(- (n + k) s + n p + k + 1, (n + k) s - n p + 1)

pdf_beta_k(n, p, k, s) := P_beta_prob(n, p, k, s) * (n+k);

pdf_beta(n, p, s) = lim k->inf pdf_beta(n, p, k, s)

// Plotting
plot3d(pdf_beta_k(500, 0.3, k, s), [s, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(100, 0.3, k, s), [s, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(50, 0.3, k, s), [s, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(10, 0.3, k, s), [s, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(5, 0.3, k, s), [s, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(n, 0.3, 500, s), [s, 0, 1], [n, 1, 500], [grid, 100, 100], [z, 0, 100]);

// Typical b
b = 0.9
```

Alpha and beta parameter of a Beta distribution using mean and mode
-------------------------------------------------------------------

```math
mean = alpha / (alpha + beta)
mean*alpha + mean*beta = alpha
beta = (alpha - mean*alpha) / mean
beta = (alpha * (1 - mean)) / mean
beta = alpha * (1/mean - 1)

mode = (alpha - 1) / (alpha + beta - 2)
mode * alpha + mode * beta - 2*mode = alpha - 1
beta = (alpha - 1 - mode*alpha + 2*mode) / mode

(alpha - mean*alpha) / mean = (alpha - 1 - mode*alpha + 2*mode) / mode
alpha/mean - alpha = alpha/mode - 1/mode - alpha + 2
alpha/mean - alpha - alpha/mode + alpha = - 1/mode + 2
alpha (1/mean - 1/mode) = -1/mode + 2
alpha = (-1/mode + 2) * (1/mean - 1/mode)
alpha = (-1/mode + 2*mode/mode) * (mode/(mean*mode) - mean/(mean*mode))
alpha = ( (2*mode-1) / mode ) * ( (mode - mean) / (mean*mode) )
alpha = ((2*mode-1) * (mode - mean)) / (mode^2*mean)

beta = ((2*mode-1) * (mode - mean)) / (mode^2*mean) * (1/mean - 1)
beta = - ((mean - 1) (mode - mean) (2*mode - 1)) / (mean^2*mode^2)
```
