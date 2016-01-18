#!/usr/bin/env runhaskell

{-# LANGUAGE BangPatterns #-}

import Control.Concurrent (threadDelay)
import Data.Function.Memoize (memoize)
import System.Environment (getArgs)
import Text.Format (format)
import Data.MultiMap (MultiMap, fromList, toMap, fromMap, mapKeys)
import Data.Map (Map, map, toList, foldr, size, lookup, unionWith,
                 foldWithKey, empty, insertWith, filter)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange, PNG),
                                PlotStyle, defaultStyle,
                                lineSpec, LineSpec(CustomStyle),
                                LineAttr(LineTitle))

import TVToolBox

---------------
-- Functions --
---------------

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
    sA = 0.05
    cA = 0.8
    nA = confidenceToCount cA 800
    sB = 0.2
    cB = 0.9999
    nB = confidenceToCount cB 800
    sC = 1.0
    cC = 0.8
    nC = confidenceToCount cC 800
    sAB = 1.0
    cAB = 0.9999
    nAB = confidenceToCount cAB 800
    sBC = 0.3
    cBC = 0.8
    nBC = confidenceToCount cBC 800
    k = 10000 -- defaultK
    resolution = 400 -- defaultResolution       -- number of bins in the distribution

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
    trimDis = (trim 1e-6) . (discretize resolution)
    genTrim s n = trimDis (genDist s n k)
    !hA = genTrim sA nA
    !hB = genTrim sB nB
    !hC = genTrim sC nC
    !hAB = genTrim sAB nAB
    !hBC = genTrim sBC nBC

  putStrLn ("hA: " ++ (showDist hA))
  putStrLn ("hB: " ++ (showDist hB))
  putStrLn ("hC: " ++ (showDist hC))
  putStrLn ("hAB: " ++ (showDist hAB))
  putStrLn ("hBC: " ++ (showDist hBC))

  let lineTitle name strength count =
          format "{0}.tv(s={1}, n={2})"
          [name, show strength, show count]
  plotDists
    [(lineTitle "A" sA nA, hA),
     (lineTitle "B" sB nB, hB),
     (lineTitle "C" sC nC, hC),
     (lineTitle "AB" sAB nAB, hAB),
     (lineTitle "BC" sBC nBC, hBC)]
    "plots/deduction-premises" False True

  -- Compute the result of deduction
  let
    !hAC = fullDeduction hA hB hC hAB hBC
  putStrLn ("hAC: " ++ (showDist hAC))

  -- Normalize the distribution
  let
    !hACnorm = (trimDis . normalize) hAC
  putStrLn ("hACnorm: " ++ (showDist hACnorm))

  -- Find the resulting distribution count
  let
    -- We use the mode rather than the mean to specify the strength of
    -- the conclusion (see README.md for the explanation)
    sAC = mode hACnorm
    -- sAC = mean hACnorm

    meanForCountSearch = sAC
    nAClow = 0
    nACup = maximum [nA, nB, nC, nAB, nBC]
    nACguess = min nAB nBC
    optimizeCount f = optimizeDbg f integerMiddle 20 nAClow nACup nACguess

    -- -- Using sqrt JSD as metric
    -- nToSqrtJsd n = sqrtJsd (genTrim meanForCountSearch n) hACnorm
    -- memNToSqrtJsd = memoize nToSqrtJsd
    -- !nAC = optimizeCount memNToSqrtJsd
    -- !hACstv = genTrim sAC nAC

    -- Using std dev distance as metric
    stdDevAC = stdDev hACnorm
    nToStdDevDiff n = abs ((stdDev (genTrim meanForCountSearch n)) - stdDevAC)
    memNToStdDevDiff = memoize nToStdDevDiff
    !nACStdDev = optimizeCount memNToStdDevDiff
    !hACStdDevStv = genTrim sAC nACStdDev

    -- -- Using U-L distance as metric
    -- b = 0.9
    -- widthAC = indefiniteIntervalWidth b hACnorm
    -- nToWidthDiff n =  abs ((indefiniteIntervalWidth b (genTrim meanForCountSearch n)) - widthAC)
    -- memNToWidthDiff = memoize nToWidthDiff
    -- !nACWidth = optimizeCount memNToWidthDiff
    -- !hACWidthStv = genTrim sAC nACWidth

  -- Plot the distributions using only StdDev metric
  plotDists
    [("AC", hACnorm),
     (lineTitle "AC-STV-approximation" meanForCountSearch nACStdDev, hACStdDevStv)]
    "plots/deduction-conclusion" False True

  -- -- Plot the distributions using the found counts of all metrics
  -- plotDists
  --   [("AC", hACnorm),
  --    (lineTitle "AC" sAC nAC, hACstv),
  --    (lineTitle "ACStdDev" meanForCountSearch nACStdDev, hACStdDevStv),
  --    (lineTitle "ACWidth" meanForCountSearch nACWidth, hACWidthStv)]
  --   ("plots/ACWithCounts-Raw-k_" ++ show k) False True
  -- plotDists
  --   [("(zoom) AC", hACnorm),
  --    (lineTitle "(zoom) AC" sAC nAC, hACstv),
  --    (lineTitle "(zoom) ACStdDev" meanForCountSearch nACStdDev, hACStdDevStv),
  --    (lineTitle "(zoom) ACWidth" meanForCountSearch nACWidth, hACWidthStv)]
  --   ("plots/ACWithCounts-Raw-Zoom-k_" ++ show k) True True

  -- -- Plot the profile of each metric
  -- let
  --   n2funProfile fun = [(fromIntegral n, fun n) | n <- [nAClow,1..nACup]]
  --   !sqrtJsdProfile = n2funProfile memNToSqrtJsd
  --   !stdDevDiffProfile = n2funProfile memNToStdDevDiff
  --   !widthDiffProfile = n2funProfile memNToWidthDiff

  -- -- Sqrt JSD
  -- putStrLn (format "sqrtJsdProfile: size = {0}, data = {1}"
  --           [show (length sqrtJsdProfile), show sqrtJsdProfile])
  -- plotPathStyle [Title "JSD w.r.t. nAC", XLabel "nAC", YLabel "JSD sqrt",
  --                PNG (format "plots/JSDSqrtProfile-k_{0}.png" [show k])]
  --               (defaultStyle {lineSpec = CustomStyle [LineTitle "JSD sqrt"]})
  --               sqrtJsdProfile

  -- -- Std dev distance
  -- putStrLn (format "stdStdDevDiffProfile: size = {0}, data = {1}"
  --           [show (length stdDevDiffProfile), show stdDevDiffProfile])
  -- plotPathStyle [Title "Std dev diff w.r.t. nAC",
  --                XLabel "nAC", YLabel "Std dev diff",
  --                PNG (format "plots/StdDevDiffProfile-k_{0}.png" [show k])]
  --               (defaultStyle {lineSpec = CustomStyle [LineTitle "Std dev diff"]})
  --               stdDevDiffProfile

  -- -- Width distance
  -- putStrLn (format "WidthDiffProfile: size = {0}, data = {1}"
  --           [show (length widthDiffProfile), show widthDiffProfile])
  -- plotPathStyle [Title "Width diff w.r.t. nAC",
  --                XLabel "nAC", YLabel "Width diff",
  --                PNG (format "plots/widthDiffStdProfile-k_{0}.png" [show k])]
  --               (defaultStyle {lineSpec = CustomStyle [LineTitle "Width diff"]})
  --               widthDiffProfile

  -- threadDelay 100000000000
