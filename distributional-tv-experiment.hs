#!/usr/bin/env runhaskell

import Data.Number.CReal (CReal)
import Math.Combinatorics.Exact.Binomial (choose)
import System.Environment (getArgs)
import Text.Format (format)
import Data.MultiMap (MultiMap, fromList, toList, foldr, map, size)
import Debug.Trace (trace)

type MyFloat = CReal

-- Distribution type, maps a first order probability to its second
-- order probability (the probability of this probability).
type Dist = MultiMap MyFloat MyFloat

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
    cA = 0.2
    sB = 0.3
    cB = 0.3
    sC = 0.2
    cC = 0.2
    sAB = 0.3
    cAB = 0.1
    sBC = 0.4
    cBC = 0.2
    k = 10

  -- Calculate corresponding counts
  let
    nA = confidenceToCount cA k
    xA = strengthToCount sA nA
    nB = confidenceToCount cB k
    xB = strengthToCount sB nA
    nC = confidenceToCount cC k
    xC = strengthToCount sC nA
    nAB = confidenceToCount cAB k
    xAB = strengthToCount sAB nAB
    nBC = confidenceToCount cBC k
    xBC = strengthToCount sBC nBC
  putStrLn (format "nA = {0}, xA = {1}, nB = {2}, xB = {3}, nC = {4}, xC = {5}, nAB = {6}, xAB = {7}, nBC = {8}, xBC = {9}" (Prelude.map show [nA, xA, nB, xB, nC, xC, nAB, xAB, nBC, xBC]))

  -- Generate corresponding distributions
  let
    dA = genDist sA cA k
    dB = genDist sB cB k
    dC = genDist sC cC k
    dAB = genDist sAB cAB k
    dBC = genDist sBC cBC k
  putStrLn (format "dA = {0}\ndB = {1}\ndC = {2}\ndAB = {3}\ndBC = {4}" (Prelude.map (show . toList) [dA, dB, dC, dAB, dBC]))

  -- Compute the result of deduction
  let
    dAC = deduction dA dB dC dAB dBC
  putStrLn (format "dAC = {0}" [show (toList dAC)])
  putStrLn (format "dAC total = {0}, size = {1}" [show (distSum dAC),
                                                  show (Data.MultiMap.size dAC)])

  -- Normalize the distribution
  let
    dACnorm = normDist dAC
  putStrLn (format "dACnorm = {0}" [show (toList dACnorm)])
  putStrLn (format "dACnorm total = {0}, size = {1}"
            [show (distSum dACnorm), show (Data.MultiMap.size dACnorm)])

-- Using the fact that c = n / (n+k) we infer that n = c*k / (1-c)
confidenceToCount :: MyFloat -> Integer -> Integer
confidenceToCount c k = round (c*(fromInteger k) / (1 - c))

-- x = s*n
strengthToCount :: MyFloat -> Integer -> Integer
strengthToCount s n = round (s * (fromInteger n))

-- Return the sum of the probabilities of the distribution. It should
-- normally be equal to 1.0
distSum :: Dist -> MyFloat
distSum d = Data.MultiMap.foldr (+) 0 d

-- Given a simple TV <s, c> and a lookahead k, generate the
-- corresponding distribution.
genDist :: MyFloat -> MyFloat -> Integer -> Dist
genDist s c k =
  fromList [(((xreal + (cxreal)) / (nreal + kreal)), p n x k cx)
           | cx <- [0..k], let cxreal = fromIntegral cx]
  where
    n = confidenceToCount c k
    nreal = fromInteger n
    x = strengthToCount s n
    xreal = fromInteger x
    kreal = fromInteger k

-- Normalize a distribution so that it sums up to 1
normDist :: Dist -> Dist
normDist d = Data.MultiMap.map (\x -> x / ds) d 
  where ds = distSum d

-- P(x+X successes in n+k trials | x successes in n trials)
p :: Integer -> Integer -> Integer -> Integer -> MyFloat
p n x k cx = (fromIntegral num) / (fromIntegral den)
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

deduction :: Dist-> Dist -> Dist -> Dist -> Dist -> Dist
deduction dA dB dC dAB dBC =
  fromList [ (deductionFormula sA sB sC sAB sBC, pA * pB * pC * pAB * pBC)  |
             (sA, pA) <- (toList dA),
             (sB, pB) <- (toList dB),
             (sC, pC) <- (toList dC),
             (sAB, pAB) <- (toList dAB),
             (sBC, pBC) <- (toList dBC),
             deductionConsistency sA sB sAB,
             deductionConsistency sB sC sBC]
