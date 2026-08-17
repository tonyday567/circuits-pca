{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE PatternSynonyms #-}

-- | Differentiable PCA operations.
--
-- The fitted 'PCAModel' is treated as constant: this module differentiates
-- through the /use/ of a PCA (centering, projection, reconstruction, loss),
-- not yet through the spectral fit itself.
--
-- Each operation is packaged as a 'Circuit.Diff.Diff' arrow, so it composes
-- with the rest of the circuits-ad stack.
module Circuit.PCA.Diff
  ( -- * Differentiable linear maps
    scoresD,
    reconstructD,
    projectRowsD,

    -- * Differentiable loss
    reconstructionLossD,
  )
where

import Circuit.Diff.Circuit (Diff, pattern Diff)
import Circuit.PCA (PCAModel (..), projectRows)
import Circuit.PCA.Lin (fromMatrix, multM, toMatrix, transpose2)
import Control.Category
import Harpie.Array.Storable (Array)
import Numeric.LinearAlgebra qualified as LA
import Prelude hiding (id, (.))

-- | Subtract a constant mean row.  The mean is a parameter, not an input.
subtractMeanD :: Array Double -> Diff (Array Double) (Array Double)
subtractMeanD mean =
  Diff $ \x ->
    case (toMatrix x, toMatrix mean) of
      (Just m, Just muM) ->
        let muFlat = LA.flatten muM
            xc = fromMatrix (m - LA.asRow muFlat)
         in (xc, id)
      _ -> error "Circuit.PCA.Diff.subtractMeanD: rank-2 input and rank-1 mean expected"

-- | Add a constant mean row.
addMeanD :: Array Double -> Diff (Array Double) (Array Double)
addMeanD mean =
  Diff $ \x ->
    case (toMatrix x, toMatrix mean) of
      (Just m, Just muM) ->
        let muFlat = LA.flatten muM
            xbar = fromMatrix (m + LA.asRow muFlat)
         in (xbar, id)
      _ -> error "Circuit.PCA.Diff.addMeanD: rank-2 input and rank-1 mean expected"

-- | Multiply on the right by a fixed matrix.
multMD :: Array Double -> Diff (Array Double) (Array Double)
multMD w =
  Diff $ \x ->
    let y = multM x w
     in (y, \dy -> multM dy (transpose2 w))

-- | Differentiable scores: center, then project onto principal axes.
scoresD :: PCAModel -> Diff (Array Double) (Array Double)
scoresD PCAModel {pcaMean, pcaComponents} =
  multMD pcaComponents . subtractMeanD pcaMean

-- | Differentiable reconstruction from scores.
reconstructD :: PCAModel -> Diff (Array Double) (Array Double)
reconstructD PCAModel {pcaMean, pcaComponents} =
  addMeanD pcaMean . multMD (transpose2 pcaComponents)

-- | Full differentiable project-through-model: @reconstruct . scores@.
projectRowsD :: PCAModel -> Diff (Array Double) (Array Double)
projectRowsD model = reconstructD model . scoresD model

-- | Mean-squared reconstruction error, differentiable with respect to the
-- input array.
--
-- Forward: @loss = (1/n) * ||x - projectRows model x||²@
-- Pullback: @dloss/dx = (2/n) * (x - projectRows model x)@
reconstructionLossD :: PCAModel -> Diff (Array Double) Double
reconstructionLossD model =
  Diff $ \x ->
    case toMatrix x of
      Nothing -> error "Circuit.PCA.Diff.reconstructionLossD: rank-2 input expected"
      Just m ->
        let xhM = case toMatrix (projectRows model x) of Just m' -> m'; Nothing -> error "projectRows not rank-2"
            rM = m - xhM
            n = fromIntegral (LA.rows m * LA.cols m) :: Double
            loss = LA.sumElements (LA.cmap (\e -> e * e) rM) / n
         in ( loss,
              \dl -> fromMatrix (LA.cmap (\e -> 2 * dl / n * e) rM)
            )
