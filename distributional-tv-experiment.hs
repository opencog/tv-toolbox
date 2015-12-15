#!/usr/bin/env runhaskell

import Control.Concurrent (threadDelay)
import Math.Gamma (gamma)
import Data.Maybe (fromJust)
import Data.Ratio ((%))
import Data.Number.BigFloat (BigFloat, Prec10, Prec50)
import Data.Function.Memoize (memoize)
import Math.Combinatorics.Exact.Binomial (choose)
import System.Environment (getArgs)
import Text.Format (format)
import Data.MultiMap (MultiMap, fromList, toMap, fromMap, mapKeys)
import Data.Map (Map, map, toList, foldr, size, lookup, unionWith,
                 foldWithKey, empty, insertWith, filter)
-- import Data.List (length)
import Debug.Trace (trace)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange),
                                PlotStyle, defaultStyle,
                                lineSpec, LineSpec(CustomStyle),
                                LineAttr(LineTitle))

-----------
-- Types --
-----------

-- type MyFloat = Float
type MyFloat = Double
-- type MyFloat = BigFloat Prec10
-- type MyFloat = BigFloat Prec50

-- Distribution type, maps a first order probability to its second
-- order probability (or strength to probability).
type MultiDist = MultiMap MyFloat MyFloat

-- Like MultiDist but each strength is unique
type Dist = Map MyFloat MyFloat

---------------
-- Constants --
---------------

defaultK :: Integer
defaultK = 5000

defaultResolution :: Integer
defaultResolution = 100

---------------
-- Functions --
---------------

showDist :: Dist -> String
showDist h = format "size = {0}, total = {1}, data = {2}"
             [show (size h), show (distSum h), show (Data.Map.toList h)]

defaultTitle :: Attribute
defaultTitle = Title (format "Simple TV distribution (k={0})" [show defaultK])

-- Plot distributions, specifying whether to enable zoom or not.
plotDists :: Bool -> [(String, Dist)] -> IO ()
plotDists zoom nhs = plotPathsStyle attributes (Prelude.map fmt nhs)
    where attributes = [defaultTitle, XLabel "Strength", YLabel "Probability"]
                       ++ if zoom then [] else [XRange (0.0, 1.0)]
          fmt (n, h) = (defaultStyle {lineSpec = CustomStyle [LineTitle n]}, toPath h)

-- Like plotDists but plot only one distribution
plotDist :: Bool -> String -> Dist -> IO ()
plotDist zoom name h = plotDists zoom [(name, h)]

-- Turn a distribution into a plotable path
toPath :: Dist -> [(Double, Double)]
toPath h = [(realToFrac s, realToFrac p) | (s, p) <- (Data.Map.toList h)]

-- Using the fact that c = n / (n+k) we infer that n = c*k / (1-c)
confidenceToCount :: MyFloat -> Integer -> Integer
confidenceToCount c k = round (c*(fromInteger k) / (1 - c))

-- x = s*n
strengthToCount :: MyFloat -> Integer -> Integer
strengthToCount s n = round (s * (fromInteger n))

-- Return the sum of the probabilities of the distribution. It should
-- normally be equal to 1.0
distSum :: Dist -> MyFloat
distSum d = Data.Map.foldr (+) 0 d

-- Given a simple TV <s, c> and a lookahead k, generate the
-- corresponding (multi)-distribution.
genMultiDist :: MyFloat -> Integer -> Integer -> MultiDist
genMultiDist s n k =
  fromList [(fromRational ((x+cx) % (n+k)), prob n x k cx) | cx <- [0..k]]
  where x = strengthToCount s n

-- Like genDist but output a distribution directly
genDist :: MyFloat -> Integer -> Integer -> Dist
genDist s n k = toDist (genMultiDist s n k)

-- Multiply the probabilities of a distribution by a given value
scale :: MyFloat -> Dist -> Dist
scale c h = Data.Map.map ((*) c) h

-- Normalize a distribution so that it sums up to 1
normalize :: Dist -> Dist
normalize h = scale (1.0 / (distSum h)) h

-- Add 2 distributions
add :: Dist -> Dist -> Dist
add = unionWith (+)

-- Compute the average of 2 distributions, (hP + hQ) / 2
average :: Dist -> Dist -> Dist
average hP hQ = scale 0.5 (add hP hQ)

-- Return the nearest (lower) bin corresponding to a strength
bin :: Integer -> MyFloat -> MyFloat
bin n s = fromRational ((floor (s * (fromInteger n))) % n)

-- Turn a multi-distribution into a distribution (sum up the
-- duplicated probabilities).
toDist :: MultiDist -> Dist
toDist d = Data.Map.map sum (toMap d)

-- Turn a distribution into a multi-distribution
toMultiDist :: Dist -> MultiDist
toMultiDist = Data.MultiMap.fromList . Data.Map.toList

-- Discretize a distribution in n bins
discretize :: Integer -> Dist -> Dist
discretize n h = toDist (mapKeys_fix (bin n) (toMultiDist h))

-- Discard probabilities under a given value
trim :: Double -> Dist -> Dist
trim e h = Data.Map.filter ((<=) e) h

-- P(x+X successes in n+k trials | x successes in n trials)
prob :: Integer -> Integer -> Integer -> Integer -> MyFloat
prob n x k cx = fromRational (num % den)
  where num = ((n+1)*(choose k cx)*(choose n x))
        den = ((k+n+1)*(choose (k+n) (cx+x)))

-- Deduction consistency, minimum value P(B|A) considering the
-- smallest possible intersection between A and B. That is
-- max((sA+sB-1)/sA, 0)
deductionSmallestIntersection :: MyFloat -> MyFloat -> MyFloat
deductionSmallestIntersection 0 sB = 10.0   -- To deal with division
                                            -- by zero, this will make
                                            -- sure that the condition
                                            -- is not met. (Just in
                                            -- case lazyness doesn't
                                            -- take care of that.)
deductionSmallestIntersection sA sB = max ((sA + sB -1) / sA) 0.0

-- Deduction consistency, maximum value P(B|A) considering the largest
-- possible intersection between A and B. That is min(1, sB/sA)
deductionLargestIntersection :: MyFloat -> MyFloat -> MyFloat
deductionLargestIntersection 0.0 sB = -10.0 -- To deal with division
                                            -- by zero, this will make
                                            -- syre that the condition
                                            -- is not met. (Just in
                                            -- case lazyness doesn't
                                            -- take care of that.)
deductionLargestIntersection sA sB = min (sB / sA) 1.0

-- Deduction consistency max((sA+sB-1)/sA, 0) <= sAB <= min(sB/sA, 1)
deductionConsistency :: MyFloat -> MyFloat -> MyFloat -> Bool
deductionConsistency sA sB sAB =
  and [0 < sA,
       deductionSmallestIntersection sA sB <= sAB,
       sAB <= deductionLargestIntersection sA sB]

-- Deduction formula using single probabilities
-- sAC = [sAB.sBC + ((1-sAB)(sC - sBsBC))/(1-sB)]
deductionFormula :: MyFloat -> MyFloat -> MyFloat -> MyFloat -> MyFloat -> MyFloat
deductionFormula sA sB sC sAB sBC =
  let result = 
        if 0.9999 < sB
        then sC                         -- avoid division by zero
        else sAB * sBC + ((1 - sAB) * (sC - sB * sBC) / (1 - sB))
  in
    -- trace (format "Trace deductionFormula(sA={0},sB={1},sC={2},sAB={3},sBC={4})={5}"
    --        (Prelude.map show [sA, sB, sC, sAB, sBC, result]))
    result

fullDeduction :: Dist-> Dist -> Dist -> Dist -> Dist -> Dist
fullDeduction hA hB hC hAB hBC = toDist dAC
    where dAC = fromList [ (deductionFormula sA sB sC sAB sBC, pAC) |
                           (sA, pA) <- (Data.Map.toList hA),
                           (sB, pB) <- (Data.Map.toList hB),
                           (sC, pC) <- (Data.Map.toList hC),
                           (sAB, pAB) <- (Data.Map.toList hAB),
                           (sBC, pBC) <- (Data.Map.toList hBC),
                           let pAC = pA * pB * pC * pAB * pBC, 1e-10 < pAC,
                           deductionConsistency sA sB sAB,
                           deductionConsistency sB sC sBC ]

-- Perform Deduction by sampling m times the distributions
monteCarloDeduction :: Integer -> Dist-> Dist -> Dist -> Dist -> Dist -> Dist
monteCarloDeduction m hA hB hC hAB hBC = undefined -- TODO



-- MultiMap mapKeys fix
mapKeys_fix f m = fromMap (foldWithKey f' empty m')
    where m' = toMap m
          f' k a b = insertWith (++) (f k) a b

-- Base 2 log
log2 = logBase 2

-- Compute the Kullback-Leibler divergence from hP to hQ. It assumes
-- that hQ has the same support as hP or greater, and that their
-- probabilities are non-null.
kld :: Dist -> Dist -> MyFloat
kld hP hQ = sum [p * (log2 (p / (hQp s))) | (s, p) <- toList hP]
    where hQp s = fromJust (Data.Map.lookup s hQ)

-- Compute the Jensen-Shannon divergence between hP and hQ. hM must be
-- the average distribution between hP and hQ.
jsd :: Dist -> Dist -> Dist -> MyFloat
jsd hP hQ hM = ((kld hP hM) + (kld hQ hM)) / 2.0

-- Compute the square root of the Jensen-Shannon divergence between hP
-- and hQ (to be a true metric). Generate hM on the fly as well.
sqrtJsd :: Dist -> Dist -> MyFloat
sqrtJsd hP hQ = sqrt (jsd hP hQ (average hP hQ))

-- This function takes a function: strength x probability -> value,
-- and distribution, and accumulate the result of this function over a
-- distribution.
accumulateWith :: (MyFloat -> MyFloat -> MyFloat) -> Dist -> MyFloat
accumulateWith f = foldWithKey (\s p r -> r + (f s p)) 0.0

-- Compute the mean of a distribution.
mean :: Dist -> MyFloat
mean = accumulateWith (*)

-- Compute the mode of a distribution (the strength with the highest
-- probability). -1.0 if the distribution is empty.
mode :: Dist -> MyFloat
mode h = fst (foldWithKey max_p (-1.0, -1.0) h)
    where max_p s p (s_max_p, p_max_p) = if p > p_max_p
                                       then (s, p)
                                       else (s_max_p, p_max_p)

-- Compute the variance of a distribution.
variance :: Dist -> MyFloat
variance h = accumulateWith (\s p -> p*(s - m)**2.0) h
    where m = mean h

-- Compute the standard deviation of a distribution.
stdDev :: Dist -> MyFloat
stdDev = sqrt . variance

-- Compute the interval [L, U] of a distribution as to minimize U-L
-- and such that (b*100)% of it is in this interval.
indefiniteInterval :: MyFloat -> Dist -> (MyFloat, MyFloat)
indefiniteInterval b h = (low, up)
    where min_low = 0
          max_low = mode h
          guess = max (max_low - b*stdDev h) min_low
          fun l = (toUp l) - l
          jump = floatMiddle
          step = 2.0 / (fromInteger defaultResolution)
          toUp l = fst (foldWithKey f (0.0, 0.0) h)
              where f s p (u, a) | s < l || b <= a = (u, a)
                                 | otherwise = (s, a + p)
          low = optimize fun jump step min_low max_low guess
          up = toUp low

-- Return the width of
indefiniteIntervalWidth :: MyFloat -> Dist -> MyFloat
indefiniteIntervalWidth b h = up - low
    where (low, up) = indefiniteInterval b h

-- Find the middle between 2 integers using the division function d
-- (I'm sure I can do better by defining a division that works for
-- both Integral and Fractional types).
middle :: Num a => (a -> a -> a) -> a -> a -> a
middle d l u = d (l + u) 2

integerMiddle = middle div
floatMiddle = middle (/)

-- Find the value x that minimzes fun x, assuming x is in [low, up],
-- and given an initial guess. It is strongly adviced to memoize fun
-- beforehand.
--
-- jump is a function that take the new interval and returns a new
-- guess (it would typically be the middle of the new interval).
--
-- step is the smaller meaningful change in x. If the function is
-- noisy setting a larger step can be useful.
optimize :: (Num a, Ord a, Show a) =>
            (a -> MyFloat) -> (a -> a -> a) -> a -> a -> a -> a -> a
optimize fun jump step low up guess
    | width < step = guess
    | up < rs = if fun guess <= fun ls then guess else rec_optimize low ls lj
    | ls < low = if fun guess <= fun rs then guess else rec_optimize rs up rj
    | otherwise = if fun guess <= fun ls && fun guess <= fun rs then guess
                  else if (fun rs) - (fun ls) >= 0
                       then rec_optimize low ls lj
                       else rec_optimize rs up rj
    where width = up - low
          ls = guess - step  -- step to the left
          rs = guess + step  -- step to the right
          lj = jump low ls  -- jump to the left
          rj = jump rs up   -- jump to the right
          rec_optimize = optimize' fun jump step -- Simplified
                                                -- recursive call of
                                                -- optimize

optimize' fun jump step low up guess =
    trace (format "optimize fun {0} {1} {2} = {3}"
           (Prelude.map show [low, up, guess, result]))
    result where result = optimize fun jump step low up guess

----------
-- Main --
----------

main :: IO ()
main = do
  -- Parse arguments
  -- [sA_str, cA_str, sB_str, cB_str, sC_str, cC_str,
  --  sAB_str, cAB_str, sBC_str, cBC_str, k_str] <- getArgs
  let
    -- sA = read sA_str :: MyFloat
    -- cA = read cA_str :: MyFloat
    -- sB = read sB_str :: MyFloat
    -- cB = read cB_str :: MyFloat
    -- sC = read sC_str :: MyFloat
    -- cC = read cC_str :: MyFloat
    -- sAB = read sAB_str :: MyFloat
    -- cAB = read cAB_str :: MyFloat
    -- sBC = read sBC_str :: MyFloat
    -- cBC = read cBC_str :: MyFloat
    -- k = read k_str :: Integer
    -- Directly enter them
    sA = 0.5
    nA = 500
    sB = 0.3
    nB = 200
    sC = 0.2
    nC = 10
    sAB = 0.8
    nAB = 400
    sBC = 0.9
    nBC = 500
    k = defaultK
    resolution = defaultResolution       -- number of bins in the distribution

  -- Calculate corresponding counts
  let
    xA = strengthToCount sA nA
    xB = strengthToCount sB nB
    xC = strengthToCount sC nC
    xAB = strengthToCount sAB nAB
    xBC = strengthToCount sBC nBC
  putStrLn (format "xA = {0}, xB = {1}, xC = {2}, xAB = {3}, xBC = {4}"
            (Prelude.map show [xA, xB, xC, xAB, xBC]))

  -- Generate corresponding distributions
  let
    trimDis = (trim 1e-10) . (discretize resolution)
    genTrim s n = trimDis (genDist s n k)
    hA = genTrim sA nA
    hB = genTrim sB nB
    hC = genTrim sC nC
    hAB = genTrim sAB nAB
    hBC = genTrim sBC nBC

  putStrLn ("hA: " ++ (showDist hA))
  putStrLn ("hB: " ++ (showDist hB))
  putStrLn ("hC: " ++ (showDist hC))
  putStrLn ("hAB: " ++ (showDist hAB))
  putStrLn ("hBC: " ++ (showDist hBC))

  let lineTitle name strength count =
          format "{0}.tv(s={1}, n={2})" [name, show strength, show count]
  plotDists False [(lineTitle "A" sA nA, hA),
                   (lineTitle "B" sB nB, hB),
                   (lineTitle "C" sC nC, hC),
                   (lineTitle "AB" sAB nAB, hAB),
                   (lineTitle "BC" sBC nBC, hBC)]

  -- Compute the result of deduction
  let
    hAC = trimDis (fullDeduction hA hB hC hAB hBC)
  putStrLn ("hAC: " ++ (showDist hAC))

  -- Normalize the distribution
  let
    hACnorm = normalize hAC
  putStrLn ("hACnorm: " ++ (showDist hACnorm))

  -- Find the resulting distribution count
  let
    sAC = mean hACnorm
    nAClow = 1
    nACup = nA + nB + nC + nAB + nBC
    nACguess = min nAB nBC

    -- Using sqrt JSD as metric
    nToSqrtJsd n = sqrtJsd (genTrim sAC n) hACnorm
    memNToSqrtJsd = memoize nToSqrtJsd
    nAC = optimize memNToSqrtJsd integerMiddle 10 nAClow nACup nACguess
    hACstv = genTrim sAC nAC

    -- Using std dev distance as metric
    stdDevAC = stdDev hACnorm
    nToStdDevDiff n = abs ((stdDev (genTrim sAC n)) - stdDevAC)
    memNToStdDevDiff = memoize nToStdDevDiff
    nACStdDev = optimize memNToStdDevDiff integerMiddle 10 nAClow nACup nACguess
    hACStdDevStv = genTrim sAC nACStdDev

    -- Using U-L distance as metric
    b = 0.9
    widthAC = indefiniteIntervalWidth b hACnorm
    nToWidthDiff n =  abs ((indefiniteIntervalWidth b (genTrim sAC n)) - widthAC)
    memNToWidthDiff = memoize nToWidthDiff
    nACWidth = optimize memNToWidthDiff integerMiddle 10 nAClow nACup nACguess
    hACWidthStv = genTrim sAC nACWidth

  -- Plot the distributions using the found counts of all metrics
  plotDists False [("AC", hACnorm),
                   (lineTitle "AC" sAC nAC, hACstv),
                   (lineTitle "ACStdDev" sAC nACStdDev, hACStdDevStv),
                   (lineTitle "ACWidth" sAC nACWidth, hACWidthStv)]
  plotDists True [("(zoom) AC", hACnorm),
                  (lineTitle "(zoom) AC" sAC nAC, hACstv),
                  (lineTitle "(zoom) ACStdDev" sAC nACStdDev, hACStdDevStv),
                  (lineTitle "(zoom) ACWidth" sAC nACWidth, hACWidthStv)]

  -- Plot the profile of each metric
  let
    n2funProfile fun = [(fromIntegral n, fun n) | n <- [nAClow,5..nACup]]
    sqrtJsdProfile = n2funProfile memNToSqrtJsd
    stdDevDiffProfile = n2funProfile memNToStdDevDiff
    widthDiffProfile = n2funProfile memNToWidthDiff

  -- Sqrt JSD
  putStrLn (format "sqrtJsdProfile: size = {0}, data = {1}"
            [show (length sqrtJsdProfile), show sqrtJsdProfile])
  plotPathStyle [Title "JSD w.r.t. nAC", XLabel "nAC", YLabel "JSD sqrt"]
                (defaultStyle {lineSpec = CustomStyle [LineTitle "JSD sqrt"]})
                sqrtJsdProfile

  -- Std dev distance
  putStrLn (format "stdStdDevDiffProfile: size = {0}, data = {1}"
            [show (length stdDevDiffProfile), show stdDevDiffProfile])
  plotPathStyle [Title "Std dev diff w.r.t. nAC",
                 XLabel "nAC", YLabel "Std dev diff"]
                (defaultStyle {lineSpec = CustomStyle [LineTitle "Std dev diff"]})
                stdDevDiffProfile

  -- Width distance
  putStrLn (format "stdWidthDiffProfile: size = {0}, data = {1}"
            [show (length widthDiffProfile), show widthDiffProfile])
  plotPathStyle [Title "Width diff w.r.t. nAC",
                 XLabel "nAC", YLabel "Width diff"]
                (defaultStyle {lineSpec = CustomStyle [LineTitle "Width diff"]})
                widthDiffProfile

  threadDelay 100000000000
