# Overview

Tool box for experimenting with TV types.

It mostly contains a Haskell module TVToolBox with a collection of
functions to experiment with TVs, in particular TV conversion (Simple
<-> Indefinite <-> Distributional), and PLN formula application.

It also contains a bit of Maxima code pasted directly in the document
(see Section Maxima further below).

# Motivation

The reason for developing this TV toolbox, and writing it in Haskell,
is twofold (letting aside the fact that Haskell is so awesome to work
with!):

  1. I'm not entirely sure where I'm going. In spite of the PLN book
     dwelling pretty deep, there are still knowledge holes to fill, so
     I need a convenient platform to experiment on, which Haskell
     provides, due to its multi-sided scripting, native, interpreting
     nature, as well as it's gnuplot binding.

  2. Haskell has pretty much unlimited integer and floating precision,
     which is very convenient to tell apart numerical errors from
     phenomenal discrepencies.

# Requirements

* ghc version 7.8.3 or above
* Libraries
  1. multimap
  2. numbers
  3. memoize
  4. gamma
  5. exact-combinatorics
  6. text-format-simple
  7. gnuplot

# Usage

To use the module see directly TVToolBox.hs for the functions it
implements. Besides that there are several files using it which
demonstrate the main functions. Each file can be directly executed as
a script (no need to compile them). Here's a short description of each
of them. For a indepth description see Section Description.

* Plot Distributional TVs given STVs

```
dtv-exp.hs
```

* Plot Indefinite TVs given STVs

```
itv-exp.sh
```

* Plot relationship between the mean of a distribution and it mode

```
mean-mode-exp.hs
```

* Plot deduction experiment using Distributional TVs

```
deduction-exp.hs
```

# Experiment Report

This section reports of a series of experiments on TV conversion and
deduction formula calculation.

## Disclaimer

I cannot make any sense of the start of Chapter 6 of the PLN book when
2 levels of sampling using Beta distributions are used to generate a
distributional TV estimate (EDIT: I think I'm partially understanding
it, see Proof p.60). Chapter 4, introducing the foundation of
distributional TV, makes perfect sense, and so this work is based on
this chapter.

Regardless of whether I eventually get it, I do intend to compare side
by side the methodology introduced in Chapter 6 with the methodology
used here. Maybe this will help me to see the light, and if not I'll
still gather knowledge to make informed implementational decisions
regarding TV conversion and formula calculation in the C++ code.

## TV conversion

### STV to DTV

Conversion from STV to DTV is done using the formula introduced in
chapter 4, Section 4.5.1 of the PLN book.

```
P(x+X successes in n+k trials | x successes in n trials)
=
(n+1)*(choose k X)*(choose n x)
--------------------------------
 (k+n+1)*(choose (k+n) (X+x))
```

where k is the lookahead. X ranges from 0 to k (included), so the
distribution is limited to k+1 data points. For small k, it may be
convenient to generate more data points, like in order to get smoother
distributions on premises and conclusions, or calculate more precise
indefinite intervals, or even just compare with the methodology of
Chapter 6. For that we use a continuous version of the binomial
coefficient using the beta function (maybe it is a bad idea, although
according to my experiments it seems coherent).

```haskell
choose_beta n k = 1.0 / ((n+1.0) * beta (n-k+1) (k+1))
```

which gives as probability

```
P((p*100)% success in n+k trials | (s*100)% success in n trials)
=
               beta(-(n+k)*p+n+k+1, (n+k)*p+1)
------------------------------------------------------------------
(k+1)*beta(-n*s+n+1, n*s+1)*beta(- n*s+(n+k)*p+1, n*s-(n+k)*p+k+1)
```

(EDIT: It might be equivalent to the formula in Chapter 6. TODO:
compare them numerically).

When n+k it large, this continuous version gets very imprecise (this
is an implementation limit, not a limit on Haskell's floating number
representation), so we use Haskell's binomial function `choose` that
works well on arbitrary large integers. In other words, when n+k is
small we use `choose_beta` that provides some interpolation, when n+k
is large we use `choose` that remains very accurate. That way we can
experiment with both small and large n+k.

See figures (obtained with `dtv-exp.hs`)

![](https://github.com/opencog/tv-toolbox/blob/master/plots/pdf_s_0.5_n_10_k_1_100.png)
![](https://github.com/opencog/tv-toolbox/blob/master/plots/pdf_s_0.7_n_1000_k_1_1000.png)
![](https://github.com/opencog/tv-toolbox/blob/master/plots/pdf_s_0.2_n_100000_k_1_1000.png)

for 3d plots of distributions obtained from simple TVs, varying the
strength, the count and the lookahead. As you may see the
corresponding distributions are a narrower when k is low. But I don't
think this alone a good reason to set k as a low value while doing
inferences. The real question is whether the width of the infered TVs
grows at a higher rate as measure as inferences progress when k is
high versus low. The PLN book seems to say that is the case, I
personally would like to measure that first hand.

### DTV to ITV

In order to convert a Distributional TV to an Indefinite TV we need to
find an indefinite interval [L, U] so that U-L, the width, is as small
as possible, yet takes up about the portion set by the indefinite
confidence level, b.

To do that we perform a optimization on L so that U-L is minimized and
the cumulative probability between L and U is equal to b or
greater. The optimization algo is a mere bisective search (with some
parameters to overcome the presence of noise).

See figures (obtained with `itv-exp.hs`)

TODO

for the results of this optimization.

### DTV to STV

In order to convert a DTV back to an STV we need to find the strength
and the count that fits the most a given distributional TV.

It may be quite important to get an accurate conversion from DTV to
STV when applying formula using DTV internally but STV externally. The
problem of course is that the DTV of a conclusion may not fit well the
DTV corresponding to some STV. According to preliminary experiments it
seems it does fit well enough though, indicating that STV on steroid
(see Section STV on Steroid) might work well enough.

For the strength we take the mode of the distributional TV, rather
than its mean. For an explanation as to why we take the mode see the
following section.

For the count we find it by optimization. The optimization here again
is a mere bisective search, and no effort is put on run-time
efficiency at that point. All fitness functions rely on generating an
STV with strength equal to the mode of the DTV, and the count being
evaluated for fitness. Three fitnesses have been tried:

1. Jensen-Shannon divergence between the DTV and the generated STV.
2. Standard deviation distance between the DTV and the generated STV.
3. Width distance between the corresponding ITV and the generated STV.

The Jensen-Shannon divergence based fitness function is supposed to be
optimal. However it is noisy (EDIT: I think this is due to some shit I
need to fix) and I suppose would only be truly optimal if the STV
strength was adjusted as well (perhaps micro-adjusted, but still).

The indefinite TV width based fitness function is the one suggested in
the PLN book, it works but it is somewhat noisy.

The standard deviation based fitness function is not noisy at all. It
may also be computationally cheaper than the other two. For that
reason it is the fitness of choice for count search.

Although our optimization algo can cope with noise, it requires a bit
of tweaking on the user's part, while optimizing a non noisy function
just works.

See figures (obtained with `deduction-exp.hs`)

![](https://github.com/opencog/tv-toolbox/blob/master/plots/JSDSqrtProfile-k_5000.png)
![](https://github.com/opencog/tv-toolbox/blob/master/plots/StdDevDiffProfile-k_5000.png)
![](https://github.com/opencog/tv-toolbox/blob/master/plots/widthDiffStdProfile-k_5000.png)

for profile comparisons of the 3 fitness functions for a deduction
conclusion AC. See figures

![](https://github.com/opencog/tv-toolbox/blob/master/plots/ACWithCounts-Zoom-K_5000.png)

for distributional TV and then estimated nearest STVs according to
these 3 fitness functions. As you may see the stdDev based count
optimization yield the best result in this case, clearly just due to
the fact that it doesn't have to deal with the noise.

### Mode vs mean

The mean of a distribution associated to an STV differs from it's
strength. This difference is very significant when n is low, as shown
in the figure below (obtained with `mean-mode-exp.hs`).

![](https://github.com/opencog/tv-toolbox/blob/master/plots/modeToMean-n_0_30_k100.png)

TODO: add fig of distribution with low n with low or high strength. 

However, the strength is equal to its mode whenever k is even or k
tends to infinity. See

```
P(x+X successes in n+k trials | x successes in n trials)
=
(n+1)*(choose k X)*(choose n x)
--------------------------------
 (k+n+1)*(choose (k+n) (X+x))
```

Removing constant factors we obtain

```
max_X (choose k X) / (choose (k+n) (X+x))
```

as choose n k is maximized when k=n/2 (EDIT: I didn't manage to proove
it yet, but I'm pretty sure it's true).

It should be noted therefore that for low n, the naive PLN formulae
are probably wrong on the strength calculations as well, not just the
count. The premises means should probably used instead of the
strengths. The conclusion's mean should then be translated back into a
strength (i.e. a mode).

## Deduction using DTV

TODO.

## Further work

TODO (in brief):
- [ ] Look at whether the rate of distribution widening increases when
  k is greater (rate of n, and U-L).
- [ ] Attempt to infer some STV formulae using DTV internally (STV
  formulae on steroid, although I've realized this might not be
  feasable for formulae involving more than a couple of premises, the
  deduction formula has already 5 premises, though, but anyway)
- [ ] Re-implement C++ double sampling in Haskell and compare with
  current single sampling results
- [ ] Implement Monte-Carlos instead of full product

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
