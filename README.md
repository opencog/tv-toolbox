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
     nature, as well as its gnuplot binding.

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
of them. For a indepth description see Section Experiment Report.

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

## TV Conversion

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

It might be equivalent to the formula in Chapter 6. Unfortunately I
cannot make any sense of that one (s is multipled by k, while p is
multipled by n, contrary to here and previous notations in the PLN
book, and it's not just a swap, p should be multipled by (n+k). I
don't know how to fix/use it). TODO: attempt to make sense of it
anyway, or reimplement here the double sampling in the obsolete C++
PLN code (see the obsolete-C++-PLN tag on the opencog repo).

When n+k it large, this continuous version gets very imprecise (this
is an implementation limit, not a limit on Haskell's floating number
representation), so we use Haskell's binomial function `choose` that
works well on arbitrary large integers. In other words, when n+k is
small we use `choose_beta` that provides some interpolation, when n+k
is large we use `choose` that remains very accurate. That way we can
experiment with both small and large n+k.

See figures (obtained with `dtv-exp.hs`)

![](plots/pdf-s_0.5_n_10_k_1_100.png)
![](plots/pdf-s_0.7_n_1000_k_1_1000.png)
![](plots/pdf-s_0.2_n_100000_k_1_1000.png)

for 3d plots of distributions obtained from simple TVs, varying the
strength, the count and the lookahead. As you may see the
corresponding distributions are a narrower when k is low. But I don't
think this alone is a good reason to set k as a low value while doing
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

See figures (obtained with `itv-exp.hs`) showing the indefinite intervals for different strengths, counts and lookaheads.

![](plots/itv-s_0.3_n_100_k_5000.png)
![](plots/itv-s_0.5_n_10_k_5000.png)
![](plots/itv-s_0.7_n_1000_k_5000.png)
![](plots/itv-s_0.3_n_100_k_50.png)
![](plots/itv-s_0.5_n_10_k_50.png)
![](plots/itv-s_0.7_n_1000_k_50.png)

### DTV to STV

In order to convert a DTV back to an STV we need to find the strength
and the count that best fits a given distributional TV.

It may be quite important to get an accurate conversion from DTV to
STV when applying formula using DTV internally but STV externally. The
problem of course is that the DTV of a conclusion may not fit well the
an STV. However, according to preliminary experiments it seems it does
fit well, indicating that STV inference on steroid (see Section STV
Inference on Steroid) might work well enough.

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

![](plots/JSDSqrtProfile-k_5000.png)
![](plots/StdDevDiffProfile-k_5000.png)
![](plots/widthDiffStdProfile-k_5000.png)

for profile comparisons of the 3 fitness functions for a deduction
conclusion AC. See figures

![](plots/ACWithCounts-Zoom-K_5000.png)

for distributional TV and then estimated nearest STVs according to
these 3 fitness functions. As you may see the stdDev based count
optimization yield the best result in this case, clearly just due to
the fact that it doesn't have to deal with noise.

### Mode vs mean

The mean of a distribution associated to an STV differs from it's
strength. This difference is very significant when n is low, as shown
in the figure below (obtained with `mean-mode-exp.hs`).

![](plots/modeToMean-n_0_30_k100.png)
![](plots/dtv-s_0.1_n_10_k_100.png)
![](plots/dtv-s_0.1_n_100_k_100.png)

It's been numerically observed that strength of an STV is equal to the
mode of its DTV (whenever k*x is divisible by n, or k tends to
infinity).

Here's a half-finished proof of it.

Given

```
    x
s = -
    n

    x+X
p = ---
    n+k
```

we want to proove that P (from Section 4.5.1 of the PLN book) is
maximized when `p = s`, that is

```
    k*x
X = ---
     n
```

Let's look at P:

```
P(x+X successes in n+k trials | x successes in n trials)
=
(n+1)*(choose k X)*(choose n x)
-------------------------------
 (k+n+1)*(choose (k+n) (X+x))
```

Removing constant factors we obtain

```
    choose k X
------------------
choose (k+n) (X+x)
```

We want to show that following

```
    k*x
X = ---
     n
```

maximizes it. That is for all X in [0, k]

```
    choose k X            choose k (k*x/n)
------------------ <= ------------------------
choose (k+n) (X+x)    choose (k+n) ((k*x/n)+x)
```

which is equivalent to

```
     k!/(X!(k-X)!)                   k!/((k*x/n)!(k-k*x/n)!)
------------------------- <= -------------------------------------
(k+n)!/((X+x)!(k+n-X-x)!)    (k+n)!/(((k*x/n)+x)!(k+n-(k*x/n)-x)!)
```

equivalent to

```
(X+x)!(k+n-X-x)!    ((k*x/n)+x)!(k+n-(k*x/n)-x)!
---------------- <= ----------------------------
    X!(k-X)!             (k*x/n)!(k-k*x/n)!
```

TODO: complete the proof (not very important though).

## Deduction Using DTV

We focus solely on the deduction rule for now. It is probably the most
important one in PLN, but also one of the hardest to "put on steroid"
because it has many premises, 5 in total, AB, BC, A, B, C.

To obtain the DTV of a conclusion we

1. Generate the DTV corresponding to each premises (see Subsection STV
   to DTV).
2. Sample the DTVs of all premises.
3. Check the consistency of the sampled probabilities.
4. If they are consistent, calculate the resulting strength using the
   deduction formula. Use this value as sample of the DTV's
   conclusion.
5. Convert the resulting DTV back to an STV (see Subsection DTV to
   STV).

Note that for now sampling the DTV merely is done in brute force, just
enumerating all discretized strengths weighted by their second order
probabilities, thus the resulting strength is also weighted by the
product of their second order probabilities.

TODO: possibly implement Monte Carlos sampling.

![](plots/deduction-premises.png)
![](plots/deduction-conclusion.png)

## STV Inference on Steroid

The idea of STV inference on steroid is to generate simple, efficient
and hopefully accurate inference formulae over STVs based on DTVs
internal computations. Think of a digital mixing consoles performing
internal signal processing in 48-bit while inputting and outputting
signal in 24-bit coupled with super-compilation.

Given a certain inference rule accompanied with its formula `f` it
would

1. Randomly pick input STVs, varying strengths and confidences and
   compute the resulting STV
   1. Turn them into DTVs
   2. Compute the resulting DTV using `f`
   3. Converte the resulting DTV into an STV
2. Record in a dataset all input STVs (strength and confidence) and
   their corresponding output STV. Try to gather a large number of
   individuals, ideally around `r*exp(2*m)`, with `r=10` and `m` equal
   to the number of premises.
3. Given the dataset use multivariate polynomial regression to find a
   model that accurately compute the resulting STV given the STVs of
   the premises. The model would be composed of 2 multivariate
   polynomials, one for the strength, one for the confidence.

## Further Work

Here's a description of the work that remains to be done, in more or
less the right order to go about them.

### Improve mode estimation

Naively computing the mode of a sampled distribution (picking up the
most frequent one) is lot less robust than the mean (averaging). I
noticed that this leads to more jitters and bader fits than using the
mean when the distribution is narrow. Yet the mode clearly is the
correct statistics for the strength of the corresponding STV.

One way to improve that would probably to try to fit a beta
distribution and use it's estimated mode, in fact there might just be
a formula relating the mode, the mean and standard deviation of a beta
distribution that could be used, see Subsection Relating parameters of
a Beta distribution with its mean and mode. Note that one wouldn't
need to compute the beta function at all, bypassing the beta function
limited precision issue.

The other simpler but less accurate option could be to use the mean as
an estimate of the mode when the standard deviation is sufficiently
low. But that's a bit of a hack, better go with the suggestion above.

### Make discretization have a varying step-size

So far the discretization has a constant step-size `1/resolution`. We
could make it more efficient, with larger step-sizes for low
probability regions, and shorter step-size for high probability
regions. This would increase the resolution where it matters and speed
up formula computation.

### Optimize STV to DTV for very high counts

When the count of an STV is very high, that is corresponding to a
confidence that approaches 1.0, calculating the binomial and thus the
probability distribution takes a very long time (about 10 minutes for
n=4M and k=10K). It would be good to optimize that. Clearly when a
distribution is so narrow you don't need as many datapoints as
k. Maybe this could be sampled. Alternatively the binomial itself
could be approximated and optimized. In the same vein it would be
useful to have a very efficient and accurate gamma or beta functions,
possibly replacing the need for an exact binomial.

### Implement Monte-Carlos simulation instead of full convolution

So far the DTV of a conclusion is calculated by combining virtually
all consistent discretized strengths, which is very costly.

Instead on may perform Monte-Carlos simulations, i.e. using second
order probability of a strength as the probability to pick it up.

### Try to make sense of Chapter 7 on distributional TV inference

In Section 7.2. of the PLN book it is explained how to use matrix
product and matrix inversion to perform DTV deduction and inversion. I
can't honnestly make sense of it. I supposed it is not contradictory
with the approach implemented here so far, but it feels it could
simplify and and possibly optimize computation a lot. So understanding
that part and implementing it here would be a must.

### Re-implement C++ double sampling in Haskell and compare

The obsolete C++ PLN uses a double sampling of beta distributions to
generate the DTVs. The code has been removed a long while ago but can
be found in the opencog git repository under the tag
`obsolete-C++-PLN`, under the directory `opencog/reasoning/pln`. It
would interesting to reimplement the idea here and compare it with
current single sampling approach.

### Estimate rate of DTV widening w.r.t. the lookahead

The PLN book seems to imply that the indefinite intervals tend to
widen at a more rapid rate if k, the lookahead, is large. It would be
interesting to verify that firt hand. For that we would just chain n-1
inference steps using the conclusions of the previous steps as
premises of the following steps and calculate w1/wn, where w1
corresponds to the standard deviation of the premise (or it's average
if there are more than one premise) of step 1, and wn would correspond
to the standard deviation of the conclusion of the last step.

We would obtain an average of multiple inferences of w1/wn varying the
STV of the initial premises, then look at how mean(w1/wn) varies
w.r.t. k.

### Implement STV inference on steroid

See Section STV Inference on Steriod for the details.

# Maxima

Here's a record of some work done with Maxima, symbolic manipulations
and plotting. Mostly trying to calculate and approach

```math
P(x+X successes in n+k trials | x successes in n trials)
```

and see how far ahead we can go.

## Using binomial coefficients

Here's the maxima code to estimate the pdf of DTV using the binomial
coefficients.

```maxima
P(n, x, k, X) := ((n+1)*binomial(k, X)*binomial(n, x))
              / ((k+n+1)*binomial(k+n, X+x));
```
```
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

## Using the Beta function

Same as above but the Beta function is used instead of the binomial
coefficients. Maxima is able to perform some symbolic simplifications
along the way.

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
```

## Approximating the binomial

Given maxima limited large integer support we tried the following
approximation
https://en.wikipedia.org/wiki/Binomial_coefficient#Approximations

```maxima
binomial_approx(n, k) := 2^n / sqrt(1/2*n*%pi) * %e^(-(k-(n/2))^2/(n/2));
P_approx(n, x, k, X) := ((n+1)*binomial_approx(k, X)*binomial_approx(n, x))
                  / ((k+n+1)*binomial_approx(k+n, X+x));
```

turns out it is really bad, especially when n is high and k is
low. Not worth redoing.

## Some plotting commands

```maxima
plot3d(pdf_beta_k(500, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(100, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(50, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(10, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(5, 0.3, k, p), [p, 0, 1], [k, 1, 500], [grid, 100, 100], [z, 0, 100]);
plot3d(pdf_beta_k(n, 0.3, 500, p), [p, 0, 1], [n, 1, 500], [grid, 100, 100], [z, 0, 100]);
```

# Others

## Relating parameters of a Beta distribution with its mean and mode

This may turn out to be handy.

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

TODO: Do the same for the standard deviation. The mean and standard
deviation are rather robust statistics, so this could be then handy
for estimating the mode given the beta function parameters.

## NOTES

X = a*k

(n+k) * (n + 1) * (k! / (X! (k-X)!) * (n! / x! (n-x)!)
------------------------------------------------------
      (k+n+1) ((k+n)! / ((X+x)! (k+n-X-x)!)

(n+k) * (n + 1) * (k! / (ak! (k-ak)!) * (n! / x! (n-x)!)
------------------------------------------------------
      (k+n+1) ((k+n)! / ((ak+x)! (k+n-ak-x)!)

                    (n+k) * (n + 1) * (k! / (ak! (k-ak)!) * (n! / x! (n-x)!)
pdf?(p=a) = k->inf ---------------------------------------------------------
                            (k+n+1) ((k+n)! / ((ak+x)! (k+n-ak-x)!)

f(p=a) = (n x)p^x*(1-p)^(n-x)

pdf(p=a) == f(p=a)
