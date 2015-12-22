#!/usr/bin/env runhaskell

import Control.Concurrent (threadDelay)
import Data.Function.Memoize (memoize)
import Data.Map (fromList)
import System.Environment (getArgs)
import Text.Format (format)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange, PNG),
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
    k = 100
    resolution = defaultResolution       -- number of bins in the distribution
    trimDis = (trim 1e-10) -- . (discretize resolution)
    genTrim s n = trimDis (genDist_beta s n k)
    modeToMean s n = mean (genTrim s n) -- s, the strength is in fact the mode

  -- Plot in 3d the mapping from mode to mean varying n
  let
    min_n = 0
    max_n = 100
    step_n = 1
  plotFunc3d [Title (format "Mean w.r.t. mode (k={0})" [show k]),
              PNG (format "plots/modeToMean-n_{0}_{1}_k{2}.png"
                   [show min_n, show max_n, show k])]
    [] [0.0,0.1..1.0] [min_n,min_n+step_n..max_n] modeToMean

  -- threadDelay 100000000000
