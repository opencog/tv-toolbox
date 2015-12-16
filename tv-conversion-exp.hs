#!/usr/bin/env runhaskell

import Control.Concurrent (threadDelay)
import Data.Function.Memoize (memoize)
import System.Environment (getArgs)
import Text.Format (format)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange),
                                PlotStyle, defaultStyle,
                                lineSpec, LineSpec(CustomStyle),
                                LineAttr(LineTitle))

import TVToolBox

----------
-- Main --
----------

main :: IO ()
main = do
  let
    sA = 0.5
    nA = 500
    k = defaultK
    resolution = defaultResolution       -- number of bins in the distribution

  -- Calculate corresponding counts
  let
    xA = strengthToCount sA nA
  putStrLn ("xA = " ++ show xA)

  -- Generate corresponding distributions
  let
    trimDis = (trim 1e-10) . (discretize resolution)
    genTrim s n = trimDis (genDist s n k)
    dA = genTrim sA nA

  putStrLn ("dA: " ++ (showDist dA))

  let lineTitle name strength count =
          format "{0}.tv(s={1}, n={2})" [name, show strength, show count]
  plotDists False [(lineTitle "A" sA nA, dA)]

  threadDelay 100000000000
