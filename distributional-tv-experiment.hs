#!/usr/bin/env runhaskell

import Control.Concurrent (threadDelay)
import Math.Gamma (gamma)
import Data.Ratio ((%))
import Data.Number.CReal (CReal)
import Data.Number.BigFloat (BigFloat, Prec10, Prec50)
import Math.Combinatorics.Exact.Binomial (choose)
import System.Environment (getArgs)
import Text.Format (format)
import Data.MultiMap (MultiMap, fromList, toMap, mapKeys)
import Data.Map (Map, map, toList, foldr, size)
import Debug.Trace (trace)
import Graphics.Gnuplot.Simple (plotPath, Attribute(XRange))

-- type MyFloat = Float
type MyFloat = Double
-- type MyFloat = BigFloat Prec10
-- type MyFloat = BigFloat Prec50

-- Distribution type, maps a first order probability to its second
-- order probability (or strength to probability).
type Dist = MultiMap MyFloat MyFloat

-- Like Dist but each strength is unique
type Hist = Map MyFloat MyFloat

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
    nA = 50000
    sB = 0.3
    nB = 200
    sC = 0.2
    nC = 300
    sAB = 0.3
    nAB = 400
    sBC = 0.4
    nBC = 500
    k = 10
    hres = 100000                  -- number of bins in the histogram

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
    hA = genHist sA nA k
    hB = genHist sB nB k
    hC = genHist sC nC k
    hAB = genHist sAB nAB k
    hBC = genHist sBC nBC k

  putStrLn (format "hA = {0}" [(show . Data.Map.toList) hA])
  plotPath [XRange (0.0, 1.0)] (toPath hA)

  putStrLn (format "hB = {0}" [(show . Data.Map.toList) hB])
  plotPath [XRange (0.0, 1.0)] (toPath hB)

  putStrLn (format "hC = {0}" [(show . Data.Map.toList) hC])
  plotPath [XRange (0.0, 1.0)] (toPath hC)

  putStrLn (format "hAB = {0}" [(show . Data.Map.toList) hAB])
  plotPath [XRange (0.0, 1.0)] (toPath hAB)

  putStrLn (format "hBC = {0}" [(show . Data.Map.toList) hBC])
  plotPath [XRange (0.0, 1.0)] (toPath hBC)

  threadDelay 10000000

  -- Compute the result of deduction
  let
    hAC = deduction hA hB hC hAB hBC
  putStrLn (format "hAC = {0}" [(show . Data.Map.toList) hAC])
  putStrLn (format "hAC total = {0}, size = {1}"
            [show (histSum hAC), (show . Data.Map.size) hAC])

  -- Normalize the distribution
  let
    hACnorm = normalize hAC
  putStrLn (format "hACnorm = {0}" [(show . Data.Map.toList) hACnorm])
  putStrLn (format "hACnorm total = {0}, size = {1}"
            [show (histSum hACnorm), (show . Data.Map.size) hACnorm])

-- Turn a histogram into a plotable path
toPath h = [(realToFrac s :: Double, realToFrac p :: Double)
            | (s, p) <- (Data.Map.toList h)]

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

-- Turn a distribution into a histogram of n bins
toBinHist :: Dist -> Integer -> Hist
toBinHist d n = toHist (mapKeys (bin n) d)

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
                           deductionConsistency sA sB sAB,
                           deductionConsistency sB sC sBC,
                           let pAC = pA * pB * pC * pAB * pBC ]
