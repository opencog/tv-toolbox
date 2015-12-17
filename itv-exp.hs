#!/usr/bin/env runhaskell

-- Experiment with <[L, U], b> search

{-# LANGUAGE BangPatterns #-}

import Control.Concurrent (threadDelay)
import Data.Function.Memoize (memoize)
import Data.Map (fromList, lookupLE, keys)
import System.Environment (getArgs)
import Text.Format (format)
import Data.Maybe (fromJust)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange),
                                PlotStyle, defaultStyle,
                                lineSpec, LineSpec(CustomStyle),
                                LineAttr(LineTitle),
                                plotFunc3d)

import TVToolBox

----------
-- Main --
----------

main :: IO ()
main = do
  let
    sA = 0.3
    nA = 100
    k = defaultK
    resolution = defaultResolution       -- number of bins in the distribution

  -- Calculate corresponding count
  let
    xA = strengthToCount sA nA
  putStrLn ("xA = " ++ show xA)

  -- Generate corresponding distribution
  let
    trimDis = (trim 1e-10) . (discretize resolution)
    genTrim s n = trimDis (genDist s n k)
    !hA = genTrim sA nA

  putStrLn ("hA: " ++ (showDist hA))

  -- Search [L, U] interval
  let
    b = 0.9
    (lA, uA) = indefiniteInterval b hA

  -- -- Plot profile of width, to see if the optimization is effective
  -- let
  --   !widthProfile = fromList [(s, w * 0.1) | s <- keys hA,
  --                             let w = (toUp hA b s) - s]
  
  -- Plot the distribution with the interval
  let lineTitle name strength count =
        format "{0}.tv(s={1}, n={2})" [name, show strength, show count]
      ulTitle l u =
        format "<[{0}, {1}], {2}>" (Prelude.map show [l, u, b])
      ulDist l u h = Data.Map.fromList [(l, uly), (u, uly)]
        where
          getp (Just (s, p)) = p
          uly = (getp (lookupLE l h) + getp (lookupLE u h)) / 2

  putStrLn ("hA: " ++ (showDist hA))
  plotDists False [(lineTitle "A" sA nA, hA),
                   (ulTitle lA uA, ulDist lA uA hA)]
                   -- (ulTitle lA uA, ulDist lA uA hA),
                   -- ("Width profile (*0.1)", widthProfile)]

  threadDelay 100000000000

