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
                                Attribute(Title, XLabel, YLabel, XRange),
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
          format "{0}.tv(s={1}, n={2})" [name, show strength, show count]
  plotDists False [(lineTitle "A" sA nA, hA),
                   (lineTitle "B" sB nB, hB),
                   (lineTitle "C" sC nC, hC),
                   (lineTitle "AB" sAB nAB, hAB),
                   (lineTitle "BC" sBC nBC, hBC)]

  -- Compute the result of deduction
  let
    !hAC = trimDis (fullDeduction hA hB hC hAB hBC)
  putStrLn ("hAC: " ++ (showDist hAC))

  -- Normalize the distribution
  let
    !hACnorm = normalize hAC
  putStrLn ("hACnorm: " ++ (showDist hACnorm))

  -- Find the resulting distribution count
  let
    sAC = mean hACnorm
    nAClow = 1
    nACup = nA + nB + nC + nAB + nBC
    nACguess = min nAB nBC
    optimizeCount f = optimizeDbg f integerMiddle 20 nAClow nACup nACguess

    -- Using sqrt JSD as metric
    nToSqrtJsd n = sqrtJsd (genTrim sAC n) hACnorm
    memNToSqrtJsd = memoize nToSqrtJsd
    !nAC = optimizeCount memNToSqrtJsd
    !hACstv = genTrim sAC nAC

    -- Using std dev distance as metric
    stdDevAC = stdDev hACnorm
    nToStdDevDiff n = abs ((stdDev (genTrim sAC n)) - stdDevAC)
    memNToStdDevDiff = memoize nToStdDevDiff
    !nACStdDev = optimizeCount memNToStdDevDiff
    !hACStdDevStv = genTrim sAC nACStdDev

    -- -- Using U-L distance as metric
    b = 0.9
    widthAC = indefiniteIntervalWidth b hACnorm
    nToWidthDiff n =  abs ((indefiniteIntervalWidth b (genTrim sAC n)) - widthAC)
    memNToWidthDiff = memoize nToWidthDiff
    !nACWidth = optimizeCount memNToWidthDiff
    !hACWidthStv = genTrim sAC nACWidth

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
    n2funProfile fun = [(fromIntegral n, fun n) | n <- [nAClow,20..nACup]]
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
