#!/usr/bin/env runhaskell

import Control.Concurrent (threadDelay)
import Math.Gamma (gamma)
import Data.Ratio ((%))
import Data.Number.BigFloat (BigFloat, Prec10, Prec50)
import Math.Combinatorics.Exact.Binomial (choose)
import System.Environment (getArgs)
import Text.Format (format)
import Data.MultiMap (MultiMap, fromList, toMap, fromMap, mapKeys)
import Data.Map (Map, map, toList, foldr, size,
                 foldWithKey, empty, insertWith, filter)
import Debug.Trace (trace)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange),
                                PlotStyle, defaultStyle,
                                lineSpec, LineSpec(CustomStyle),
                                LineAttr(LineTitle))

-- type MyFloat = Float
type MyFloat = Double
-- type MyFloat = BigFloat Prec10
-- type MyFloat = BigFloat Prec50

-- Distribution type, maps a first order probability to its second
-- order probability (or strength to probability).
type Dist = MultiMap MyFloat MyFloat

-- Like Dist but each strength is unique
type Hist = Map MyFloat MyFloat

defaultK :: Integer
defaultK = 10000

defaultResolution :: Integer
defaultResolution = 150

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
  plotHists [(lineTitle "A" sA nA, hA),
             (lineTitle "B" sB nB, hB),
             (lineTitle "C" sC nC, hC),
             (lineTitle "AB" sAB nAB, hAB),
             (lineTitle "BC" sBC nBC, hBC)]

  -- Compute the result of deduction
  let
    hAC = trimDis (deduction hA hB hC hAB hBC)
  putStrLn ("hAC: " ++ (showHist hAC))

  -- Normalize the distribution
  let
    hACnorm = normalize hAC
  putStrLn ("hACnorm: " ++ (showHist hACnorm))

  -- Plot the distribution
  plotHist "AC.tv" hACnorm
  plotHistZoom "AC - zoom" hACnorm

  threadDelay 100000000

showHist :: Hist -> String
showHist h = format "size = {0}, total = {1}, data = {2}"
             [show (size h), show (histSum h), show (Data.Map.toList h)]

defaultTitle :: Attribute
defaultTitle = Title (format "Simple TV distribution (k={0})" [show defaultK])

-- Plot a distribution, provided its name, xrange and histogram
plotHistRange :: String -> Double -> Double -> Hist -> IO ()
plotHistRange name l u h =
    plotPathStyle [defaultTitle, XLabel "Strength", YLabel "Probability", XRange (l, u)]
                  (defaultStyle {lineSpec = CustomStyle [LineTitle name]})
                  (toPath h)

-- Same as plotHistsRange but plot several histograms
plotHistsRange :: [(String, Hist)] -> Double -> Double -> IO ()
plotHistsRange nhs l u =
    plotPathsStyle [defaultTitle, XLabel "Strength", YLabel "Probability", XRange (l, u)]
                   (Prelude.map fmt nhs)
        where fmt (n, h) = (defaultStyle {lineSpec = CustomStyle [LineTitle n]}, toPath h)

-- Same as plotHistRange but the xrange is (0.0, 1.0)
plotHist :: String -> Hist -> IO ()
plotHist name h = plotHistRange name 0.0 1.0 h

-- Same as plotHists for plots several histograms
plotHists :: [(String, Hist)] -> IO ()
plotHists nhs = plotHistsRange nhs 0.0 1.0

-- Same as plotHistRange but without specifying the xrange
plotHistZoom :: String -> Hist -> IO ()
plotHistZoom name h =
    plotPathStyle [defaultTitle, XLabel "Strength", YLabel "Probabilitiy"]
                  (defaultStyle {lineSpec = CustomStyle [LineTitle name]})
                  (toPath h)

-- Turn a histogram into a plotable path
toPath :: Hist -> [(Double, Double)]
toPath h = [(realToFrac s, realToFrac p) | (s, p) <- (Data.Map.toList h)]

-- Using the fact that c = n / (n+k) we infer that n = c*k / (1-c)
confidenceToCount :: MyFloat -> Integer -> Integer
confidenceToCount c k = floor (c*(fromInteger k) / (1 - c))

-- x = s*n
strengthToCount :: MyFloat -> Integer -> Integer
strengthToCount s n = floor (s * (fromInteger n))

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

-- Normalize a histogram so that it sums up to 1
normalize :: Hist -> Hist
normalize h = Data.Map.map (\x -> x / hs) h
  where hs = histSum h

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

deduction :: Hist-> Hist -> Hist -> Hist -> Hist -> Hist
deduction hA hB hC hAB hBC = toHist dAC
    where dAC = fromList [ (deductionFormula sA sB sC sAB sBC, pAC) |
                           (sA, pA) <- (Data.Map.toList hA),
                           (sB, pB) <- (Data.Map.toList hB),
                           (sC, pC) <- (Data.Map.toList hC),
                           (sAB, pAB) <- (Data.Map.toList hAB),
                           (sBC, pBC) <- (Data.Map.toList hBC),
                           let pAC = pA * pB * pC * pAB * pBC, 1e-10 < pAC,
                           deductionConsistency sA sB sAB,
                           deductionConsistency sB sC sBC ]

-- MultiMap mapKeys fix
mapKeys_fix f m = fromMap (Data.Map.foldWithKey f' empty m')
    where m' = toMap m
          f' k a b = insertWith (++) (f k) a b
