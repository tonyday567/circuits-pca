-- | Linear array ops for PCA: center, Gram (dagger-compose), BLAS mult.
--
-- Data layout is always /samples × features/ (rows = observations).
--
-- The dagger promise used here is the array transpose: covariance on the
-- feature side is @X† X@ via transpose-compose (up to scaling). Spectral
-- work is delegated to hmatrix (same bridge pattern as @Harpie.Hmatrix@,
-- which is flag-gated in harpie).
module Circuit.PCA.Lin
  ( -- * Layout
    nSamples,
    nFeatures,

    -- * Center
    columnMeans,
    centerColumns,

    -- * Dagger-compose
    gramFeatures,

    -- * BLAS
    multM,
    transpose2,

    -- * hmatrix bridge
    toMatrix,
    fromMatrix,
  )
where

import Data.Vector.Storable qualified as VS
import Data.Vector.Unboxed qualified as VU
import Harpie.Array.Storable (Array)
import Harpie.Array.Storable qualified as A
import Numeric.LinearAlgebra (Matrix, cols, flatten, reshape, rows, (<>))
import Numeric.LinearAlgebra qualified as LA
import Prelude hiding ((<>))

-- | Number of samples (rows).
nSamples :: Array Double -> Int
nSamples a = case VU.toList (A.shape a) of
  (n : _) -> n
  _ -> 0

-- | Number of features (columns).
nFeatures :: Array Double -> Int
nFeatures a = case VU.toList (A.shape a) of
  [_, p] -> p
  [p] -> p
  _ -> 0

-- | Rank-2 storable array → hmatrix (Nothing if not rank 2).
toMatrix :: Array Double -> Maybe (Matrix Double)
toMatrix a =
  case VU.toList (A.shape a) of
    [_, c] -> Just (reshape c (A.asVector a))
    _ -> Nothing

-- | hmatrix → rank-2 storable array.
fromMatrix :: Matrix Double -> Array Double
fromMatrix m = A.array [rows m, cols m] (flatten m)

-- | BLAS matrix multiply for rank-2 arrays.
multM :: Array Double -> Array Double -> Array Double
multM a b =
  case (toMatrix a, toMatrix b) of
    (Just ma, Just mb) -> fromMatrix (ma <> mb)
    _ -> error "Circuit.PCA.Lin.multM: expected rank-2 arrays"

-- | Column means as a length-@p@ vector (shape @[p]@).
columnMeans :: Array Double -> Array Double
columnMeans x =
  case toMatrix x of
    Nothing -> error "Circuit.PCA.Lin.columnMeans: expected rank-2 samples×features"
    Just m ->
      let n = fromIntegral (rows m) :: Double
       in columnMeansM m n

columnMeansM :: Matrix Double -> Double -> Array Double
columnMeansM m n =
  let p = cols m
      means =
        VS.generate p $ \j ->
          LA.sumElements (m LA.¿ [j]) / n
   in A.array [p] means

-- | Center columns; returns @(centered, means)@.
centerColumns :: Array Double -> (Array Double, Array Double)
centerColumns x =
  case toMatrix x of
    Nothing -> error "Circuit.PCA.Lin.centerColumns: expected rank-2 samples×features"
    Just m ->
      let n = fromIntegral (rows m) :: Double
          mu = columnMeansM m n
          muFlat = A.asVector mu
          xc = m - LA.asRow muFlat
       in (fromMatrix xc, mu)

-- | Feature Gram matrix @X† X@ (features × features), unnormalised.
--
-- Dagger-compose: reverse wire is transpose, compose is BLAS mult.
-- Scale by @1/(n-1)@ for sample covariance if desired.
gramFeatures :: Array Double -> Array Double
gramFeatures x =
  case toMatrix x of
    Nothing -> error "Circuit.PCA.Lin.gramFeatures: expected rank-2"
    Just m -> fromMatrix (LA.tr m <> m)

-- | Rank-2 transpose.
transpose2 :: Array Double -> Array Double
transpose2 a =
  case toMatrix a of
    Nothing -> error "Circuit.PCA.Lin.transpose2: expected rank-2"
    Just m -> fromMatrix (LA.tr m)
