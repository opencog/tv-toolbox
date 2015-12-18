Overview
--------

Tool box for experimenting with TV types.

Requirements
------------

* ghc version 7.8.3 or above
* Libraries
  1. multimap
  2. numbers
  3. memoize
  4. gamma
  5. exact-combinatorics
  6. text-format-simple
  7. gnuplot
TODO

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
// Using s=x/n and p=(x+X)/(n+k)
//
// s = x/n
// x = s*n
// p = (x+X) / (n+k) = (s*n+X) / (n+k)
// p*(n+k) = s*n+X
// X = p*(n+k) - s*n
// X+x = p*(n+k)
// p(X) = (s*n+X) / (n+k)
// p(X+1) = (s*n+X+1) / (n+k)
// p(X) - p(X+1) = 1/(n+k)
```

```maxima
P_prob(n, s, k, p) := ((n+1)*binomial(k, p*(n+k) - s*n)*binomial(n, s*n))
                        / ((k+n+1)*binomial(k+n, p*(n+k)));

pdf_k(n, s, k, p) := P_prob(n, s, k, p) * (n+k);
pdf(n, s, p) = lim k->inf pdf_k(n, s, k, p)

limit(pdf(n, s, k, p), k, inf);
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

P_beta_prob(n, s, k, p) := ((n+1)*binomial_beta(k, p*(n+k) - s*n)*binomial_beta(n, s*n))
                        / ((k+n+1)*binomial_beta(k+n, p*(n+k)));

// Discard nonsensical p
P_beta_prob_condition(n, s, k, p) := if (0 <= (p*(n + k) - s*n) and (p*(n + k) - s*n <= k)) then (((n + 1)*binomial_beta(k, p*(n + k) - s*n)*binomial_beta(n, s*n))/((k + n + 1)*binomial_beta(k + n, p*(n + k)))) else 0;

// Simplified by Maxima
(%i44) P_beta_prob(n, s, k, p);
(%o44) beta(- (n + k) p + n + k + 1, (n + k) p + 1)
/((k + 1) beta(- n s + n + 1, n s + 1) beta(- n s + (n + k) p + 1, 
n s - (n + k) p + k + 1))

pdf_beta_k(n, s, k, p) := P_beta_prob(n, s, k, p) * (n+k);
pdf_beta_condition_k(n, s, k, p) := P_beta_prob_condition(n, s, k, p) * (n+k);

pdf_beta(n, s, p) = lim k->inf pdf_beta(n, s, k, p)

// Plotting
plot3d(pdf_beta_k(500, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(100, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(50, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(10, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(5, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(n, 0.3, 500, p), [p, 0, 1], [n, 1, 500], [grid, 100, 100], [z, 0, 100]);

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
