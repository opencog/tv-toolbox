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
    sA = 0.6
    nA = 10
    k = 100000
    resolution = 1000 -- defaultResolution       -- number of bins in the distribution

  -- Calculate corresponding count
  let
    xA = strengthToCount sA nA
  putStrLn ("xA = " ++ show xA)

  -- Generate corresponding distribution
  let
    pdf_A = fromList [(p, pdf_beta nA sA k p) | p <- [0.0,0.01..1.0]]
    pdf_chapter_4_A = fromList [(p, pdf_chapter_4 nA xA p) | p <- [0.0,0.01..1.0]]


  putStrLn ("pdf_A: " ++ (showDist pdf_A))
  putStrLn ("pdf_chapter_4_A: " ++ (showDist pdf_chapter_4_A))

  plotDists [("pdf_A", pdf_A), ("pdf_chapter_4_A", pdf_chapter_4_A)] "" False False
  threadDelay 100000000000

  -- let lineTitle name strength count mean =
  --         format "{0}.tv(s={1}, n={2}, k={3}, mean={4})"
  --         [name, show strength, show count, show k, show mean]
  -- plotDists
  --   -- [(lineTitle "A" sA nA meanA, dA),
  --    [(lineTitle "Abeta" sA nA meanAbeta, dAbeta)]
  --   (format "plots/dtv-s_{0}_n_{1}_k_{2}" [show sA, show nA, show k])
  --   False True

  -- -- Plot in 3d the pdf varying k
  -- let
  --   mink = 1
  --   maxk = 1000
  --   strs = [show sA, show nA, show mink, show maxk]
  -- plotFunc3d [Title (format "pdf(<s={0},n={1}>) k=[{2}..{3}]" strs),
  --             PNG (format "plots/pdf-s_{0}_n_{1}_k_{2}_{3}.png" strs)]
  --   [] [0.0,0.001..1.0] [mink..maxk] (\p k -> pdf_beta nA sA k p)
  -- -- plotFunc3d [] [] [0.0,0.001..1.0] [1..500] (\p k -> pdf_only_beta nA sA k p)

  -- threadDelay 100000000000
