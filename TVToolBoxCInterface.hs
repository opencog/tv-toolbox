module TVToolBoxCIntegererface where

import TVToolBox

import Foreign.C.Types

--Compute Jensen-Shannon divergence given bin,TV,Count and lookahead values.
sqrtJsdC :: Integer -> MyFloat -> Integer -> Integer -> Integer -> MyFloat -> Integer -> Integer -> MyFloat
sqrtJsdC b1 s1 c1 k1 b2 s2 c2 k2 = sqrtJsd d1 d2
	where d1 = discretize b1 (genDist s1 c1 k1)
              d2 = discretize b2 (genDist s2 c2 k2)

sqrtJsdC_hs :: CInt -> CDouble -> CInt -> CInt -> CInt -> CDouble -> CInt -> CInt ->CDouble 
sqrtJsdC_hs  b1 s1 c1 k1 b2 s2 c2 k2 = realToFrac $ sqrtJsdC (fromIntegral b1) (realToFrac s1)
                                       (fromIntegral c1) (fromIntegral k1) (fromIntegral b2)
                                       (realToFrac s2) (fromIntegral c2) (fromIntegral k2)

foreign export ccall sqrtJsdC_hs :: CInt -> CDouble -> CInt -> CInt -> CInt -> CDouble -> CInt -> CInt ->CDouble 

