#!/usr/bin/env runhaskell

import Data.Number.CReal (CReal)
import Math.Combinatorics.Exact.Binomial (choose)
import System.Environment (getArgs)
import Text.Format (format)
import Data.MultiMap (MultiMap, fromList, toList, foldr)

type Dist = MultiMap CReal CReal

main :: IO ()
main = do
  [n_str, x_str, k_str] <- getArgs
  let n = read n_str :: Integer
      x = read x_str :: Integer
      k = read k_str :: Integer
      dist = fromList [((fromIntegral cx) / (fromIntegral k), p n x k cx)
                      | cx <- [0..k]]
  putStrLn (format "p({0}, {1}, {2}) = {3}"
            [n_str, x_str, k_str, show (toList dist)])
  putStrLn (format "Total = {0} (should be 1.0)"
            [show (distSum dist)])

distSum :: Dist -> CReal
distSum d = Data.MultiMap.foldr (+) 0 d

-- P(x+X successes in n+k trials | x successes in n trials)
p :: Integer -> Integer -> Integer -> Integer -> CReal
p n x k cx = (fromIntegral num) / (fromIntegral den)
  where num = ((n+1)*(choose k cx)*(choose n x))
        den = ((k+n+1)*(choose (k+n) (cx+x)))

-- Deduction consistency max
-- max((sA+sB-1)/sA, 0)
deductionConsistencyMax sA sB = max ((sA + sB -1) / sA) 0

-- Deduction consistency min
-- min(1, sB/sA)
deductionConsistencyMin sA sB = min (sB / sA) 1

-- Deduction consistency 1
-- sAB-max((sA+sB-1)/sA, 0)
deductionConsistency1 sA sB sAB = sAB - (deductionConsistencyMax sA sB)

-- Deduction consistency 2
-- min(1, sB/sA)-max((sA+sB-1)/sA, 0)
deductionConsistency2 sA sB sAB = (deductionConsistencyMin sA sB)
                                  - (deductionConsistencyMax sA sB)

-- Deduction formula using single probabilities
-- sAC = [sAB.sBC + ((1-sAB)(sC - sBsBC))/(1-sB)]
deductionFormula :: CReal -> CReal -> CReal -> CReal -> CReal -> CReal
deductionFormula sA sB sC sAB sBC =
  if 0.9999 < sB
  then sC                         -- avoid division by zero
  else sAB * sBC + ((1 - sAB) * (sC - sB * sBC) / (1 - sB))

deduction :: Dist-> Dist -> Dist -> Dist -> Dist -> Dist
deduction dA dB dC dAB dBC = undefined -- TODO
