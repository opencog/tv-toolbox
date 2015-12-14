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
type Dist = MultiMap MyFloat MyFloat

-- Like Dist but each strength is unique
type Hist = Map MyFloat MyFloat

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

showHist :: Hist -> String
showHist h = format "size = {0}, total = {1}, data = {2}"
             [show (size h), show (histSum h), show (Data.Map.toList h)]

defaultTitle :: Attribute
defaultTitle = Title (format "Simple TV distribution (k={0})" [show defaultK])

-- Plot distributions, specifying whether to enable zoom or not.
plotHists :: Bool -> [(String, Hist)] -> IO ()
plotHists zoom nhs = plotPathsStyle attributes (Prelude.map fmt nhs)
    where attributes = [defaultTitle, XLabel "Strength", YLabel "Probability"]
                       ++ if zoom then [] else [XRange (0.0, 1.0)]
          fmt (n, h) = (defaultStyle {lineSpec = CustomStyle [LineTitle n]}, toPath h)

-- Like plotHists but plot only one distribution
plotHist :: Bool -> String -> Hist -> IO ()
plotHist zoom name h = plotHists zoom [(name, h)]

-- Turn a histogram into a plotable path
toPath :: Hist -> [(Double, Double)]
toPath h = [(realToFrac s, realToFrac p) | (s, p) <- (Data.Map.toList h)]

-- Using the fact that c = n / (n+k) we infer that n = c*k / (1-c)
confidenceToCount :: MyFloat -> Integer -> Integer
confidenceToCount c k = round (c*(fromInteger k) / (1 - c))

-- x = s*n
strengthToCount :: MyFloat -> Integer -> Integer
strengthToCount s n = round (s * (fromInteger n))

-- Return the sum of the probabilities of the distribution. It should
-- normally be equal to 1.0
histSum :: Hist -> MyFloat
histSum d = Data.Map.foldr (+) 0 d

-- Given a simple TV <s, c> and a lookahead k, generate the
-- corresponding distribution.
genDist :: MyFloat -> Integer -> Integer -> Dist
genDist s n k =
  fromList [(fromRational ((x+cx) % (n+k)), prob n x k cx) | cx <- [0..k]]
  where x = strengthToCount s n

-- Like genDist but output a histogram directly
genHist :: MyFloat -> Integer -> Integer -> Hist
genHist s n k = toHist (genDist s n k)

-- Multiply the probabilities of a distribution by a given value
scale :: MyFloat -> Hist -> Hist
scale c h = Data.Map.map ((*) c) h

-- Normalize a histogram so that it sums up to 1
normalize :: Hist -> Hist
normalize h = scale (1.0 / (histSum h)) h

-- Add 2 distributions
add :: Hist -> Hist -> Hist
add = unionWith (+)

-- Compute the average of 2 distributions, (hP + hQ) / 2
average :: Hist -> Hist -> Hist
average hP hQ = scale 0.5 (add hP hQ)

-- Return the nearest (lower) bin corresponding to a strength
bin :: Integer -> MyFloat -> MyFloat
bin n s = fromRational ((floor (s * (fromInteger n))) % n)

-- Turn a distribution into a histogram (sum up the duplicated
-- probabilities).
toHist :: Dist -> Hist
toHist d = Data.Map.map sum (toMap d)

-- Turn a histogram into a distribution
toDist :: Hist -> Dist
toDist = Data.MultiMap.fromList . Data.Map.toList

-- Turn a distribution into a histogram of n bins
toBinHist :: Integer -> Dist -> Hist
toBinHist n d = toHist (mapKeys_fix (bin n) d)

-- Discretize a distribution in n bins
discretize :: Integer -> Hist -> Hist
discretize n h = toHist (mapKeys_fix (bin n) (toDist h))

-- Discard probabilities under a given value
trim :: Double -> Hist -> Hist
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

fullDeduction :: Hist-> Hist -> Hist -> Hist -> Hist -> Hist
fullDeduction hA hB hC hAB hBC = toHist dAC
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
monteCarloDeduction :: Integer -> Hist-> Hist -> Hist -> Hist -> Hist -> Hist
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
kld :: Hist -> Hist -> MyFloat
kld hP hQ = sum [p * (log2 (p / (hQp s))) | (s, p) <- toList hP]
    where hQp s = fromJust (Data.Map.lookup s hQ)

-- Compute the Jensen-Shannon divergence between hP and hQ. hM must be
-- the average distribution between hP and hQ.
jsd :: Hist -> Hist -> Hist -> MyFloat
jsd hP hQ hM = ((kld hP hM) + (kld hQ hM)) / 2.0

-- Compute the square root of the Jensen-Shannon divergence between hP
-- and hQ (to be a true metric). Generate hM on the fly as well.
sqrtJsd :: Hist -> Hist -> MyFloat
sqrtJsd hP hQ = sqrt (jsd hP hQ (average hP hQ))

-- This function takes a function: strength x probability -> value,
-- and distribution, and accumulate the result of this function over a
-- distribution.
accumulateWith :: (MyFloat -> MyFloat -> MyFloat) -> Hist -> MyFloat
accumulateWith f = foldWithKey (\s p r -> r + (f s p)) 0.0

-- Compute the mean of a distribution.
mean :: Hist -> MyFloat
mean = accumulateWith (*)

-- Compute the mode of a distribution (the strength with the highest
-- probability). -1.0 if the distribution is empty.
mode :: Hist -> MyFloat
mode h = fst (foldWithKey max_p (-1.0, -1.0) h)
    where max_p s p (s_max_p, max_p) = if p > max_p
                                       then (s, p)
                                       else (s_max_p, max_p)

-- Compute the standard deviation of a distribution.
stdev :: Hist -> MyFloat
stdev h = sqrt (accumulateWith (\s p -> p*(s - m)**2.0) h)
    where m = mean h

-- Compute the interval [L, U] of a distribution such that (b*100)% of
-- it is within this interval as to minimize U-L.
indefiniteInterval :: MyFloat -> Hist -> (MyFloat, MyFloat)
indefiniteInterval b h = optimizeFloat f 0 m lest
    where m = mode h
          lest = undefined -- Estimate of L assuming a gaussian, of course it
                 -- would be much better to assume a Beta
                 -- distribution, but then finding L and U isn't that
                 -- easier.
          f = undefined
          optimizeFloat = undefined

-- Find the middle between 2 integers
middle :: Integer -> Integer -> Integer
middle l u = div (l + u) 2

-- Find the Integer value x so that f(x) is minimized, assuming x is
-- in [l, u], and given an initial guess. It is strongly adviced to
-- memoize f beforehand. Note that [l, u] has nothing to do with the
-- indefinite interval [L, U].
optimize :: (Integer -> MyFloat) -> Integer -> Integer -> Integer -> Integer
optimize f l u m | m == l && m == u = m
                 | m == u = if f (m-1) >= f m then m
                            else optimize' f l (m-1) (middle l (m-1))
                 | m == l = if f m <= f (m+1) then m
                            else optimize' f (m+1) u (middle (m+1) u)
                 | otherwise = if f (m-1) >= f m && f m <= f (m+1) then m
                               else if (f (m+1)) - (f (m-1)) < 0
                                    then optimize' f (m+1) u (middle (m+1) u)
                                    else optimize' f l (m-1) (middle l (m-1))
optimize' f l u m = trace (format "optimize f {0} {1} {2} = {3}"
                                  (Prelude.map show [l, u, m, result]))
                    result where result = optimize f l u m

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
    resolution = defaultResolution       -- number of bins in the histogram

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
    hA = trimDis (genHist sA nA k)
    hB = trimDis (genHist sB nB k)
    hC = trimDis (genHist sC nC k)
    hAB = trimDis (genHist sAB nAB k)
    hBC = trimDis (genHist sBC nBC k)

  putStrLn ("hA: " ++ (showHist hA))
  putStrLn ("hB: " ++ (showHist hB))
  putStrLn ("hC: " ++ (showHist hC))
  putStrLn ("hAB: " ++ (showHist hAB))
  putStrLn ("hBC: " ++ (showHist hBC))

  let lineTitle name strength count =
          format "{0}.tv(s={1}, n={2})" [name, show strength, show count]
  plotHists False [(lineTitle "A" sA nA, hA),
                   (lineTitle "B" sB nB, hB),
                   (lineTitle "C" sC nC, hC),
                   (lineTitle "AB" sAB nAB, hAB),
                   (lineTitle "BC" sBC nBC, hBC)]

  -- Compute the result of deduction
  let
    hAC = trimDis (fullDeduction hA hB hC hAB hBC)
  putStrLn ("hAC: " ++ (showHist hAC))

  -- Normalize the distribution
  let
    hACnorm = normalize hAC
  putStrLn ("hACnorm: " ++ (showHist hACnorm))

  -- Find the resulting distribution count
  let
    sAC = mean hACnorm
    countToSqrtJsd n = sqrtJsd (trimDis (genHist sAC n k)) hACnorm
    memCountToSqrtJsd = memoize countToSqrtJsd
    nACl = 0
    nACu = nA + nB + nC + nAB + nBC
    nACm = min nAB nBC
    dsts = [(fromIntegral n, memCountToSqrtJsd n) | n <- [nACl..nACu]]
    nAC = optimize' memCountToSqrtJsd nACl nACu nACm
    hACstv = trimDis (genHist sAC nAC k)

  -- Plot the distribution
  putStrLn ("dsts = " ++ (show dsts))
  plotPathStyle [Title "JSD w.r.t. nAC", XLabel "nAC", YLabel "JSD sqrt"]
                (defaultStyle {lineSpec = CustomStyle [LineTitle "JSD sqrt"]}) dsts
  plotHists False [("AC", hACnorm), (lineTitle "AC" sAC nAC, hACstv)]
  plotHists True [("(zoom) AC", hACnorm), (lineTitle "(zoom) AC" sAC nAC, hACstv)]

  threadDelay 1000000000
