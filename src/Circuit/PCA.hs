{-# LANGUAGE NamedFieldPuns #-}

-- | Principal component analysis as optic residual ownership over harpie arrays.
--
-- Pipeline (see circuits @examples\/pca.md@):
--
-- @
-- center → dagger-compose (Gram / SVD of centered data) → spectral cut
--       → Optic (,) with major = focus, minor = owned residual
-- @
--
-- 'Dagger' from circuits promises reverse wires; here reverse is matrix
-- transpose and the spectral cut is hmatrix. The keep\/discard protocol is
-- 'Circuit.PCA.Optic.Lens' — hinge table cell A (own + finite residual).
--
-- 'circuits-ad' is a sibling dependency for reverse-mode work on surrounding
-- nets; this module does not yet differentiate through the spectral cut.
module Circuit.PCA
  ( -- * Model
    PCAModel (..),
    fit,
    fitSVD,

    -- * Project / reconstruct (optic face)
    scores,
    reconstruct,
    projectRows,

    -- * Dagger face
    gramDagger,
    gramViaDagger,

    -- * Single-row optic
    rowLens,
    viewMajor,
    setMajor,

    -- * Re-exports
    module Circuit.PCA.Optic,
    module Circuit.PCA.Lin,
  )
where

-- circuits-ad: reverse-mode sibling; kept linked for diff-through-PCA work.
import Circuit.AD ()
import Circuit.Dagger (Dagger (..), transpose)
import Circuit.PCA.Lin
import Circuit.PCA.Optic
import Data.Vector.Storable qualified as VS
import Harpie.Array.Storable (Array)
import Harpie.Array.Storable qualified as A
-- harpie-numhask: numhask orphans for fixed-shape harpie (streaming/fixed path later).
import Harpie.NumHask ()
import Numeric.LinearAlgebra qualified as LA
import Prelude

-- | Fitted PCA: column means, loadings (features × k), singular values.
data PCAModel = PCAModel
  { -- | Column means, shape @[p]@.
    pcaMean :: !(Array Double),
    -- | Principal axes as columns, shape @[p, k]@.
    pcaComponents :: !(Array Double),
    -- | Top-k singular values of the centered data matrix (length k).
    pcaSingularValues :: !(VS.Vector Double),
    -- | Requested component count.
    pcaK :: !Int
  }
  deriving stock (Show)

-- | Fit PCA with @k@ components via thin SVD of centered data.
--
-- Data shape: samples × features. Uses right singular vectors as loadings
-- (same subspace as eigendecomposition of the feature Gram).
fit :: Int -> Array Double -> PCAModel
fit = fitSVD

-- | Explicit SVD path (BLAS via hmatrix).
fitSVD :: Int -> Array Double -> PCAModel
fitSVD k x
  | k < 1 = error "Circuit.PCA.fitSVD: k must be >= 1"
  | otherwise =
      case toMatrix x of
        Nothing -> error "Circuit.PCA.fitSVD: expected rank-2 samples×features"
        Just m ->
          let n = LA.rows m
              p = LA.cols m
              k' = min k (min n p)
              (xc, mu) = centerColumns x
           in case toMatrix xc of
                Nothing -> error "Circuit.PCA.fitSVD: centered data not rank-2"
                Just mc ->
                  -- thin SVD: mc == u <> diag s <> tr v; columns of v = right SVs
                  let (_u, s, v) = LA.thinSVD mc
                      vk = v LA.?? (LA.All, LA.Take k')
                      sTop = VS.fromList (take k' (LA.toList s))
                   in PCAModel
                        { pcaMean = mu,
                          pcaComponents = fromMatrix vk,
                          pcaSingularValues = sTop,
                          pcaK = k'
                        }

meanRow :: Array Double -> LA.Vector Double
meanRow v = A.asVector v

-- | Project centered rows to k-dimensional scores. Shape: samples × k.
scores :: PCAModel -> Array Double -> Array Double
scores PCAModel {pcaMean, pcaComponents} x =
  case (toMatrix x, toMatrix pcaComponents) of
    (Just m, Just vk) ->
      let mu = meanRow pcaMean
          xc = m - LA.asRow mu
       in fromMatrix (xc LA.<> vk)
    _ -> error "Circuit.PCA.scores: rank-2 expected"

-- | Reconstruct data from scores: @scores × components† + mean@.
reconstruct :: PCAModel -> Array Double -> Array Double
reconstruct PCAModel {pcaMean, pcaComponents} sc =
  case (toMatrix sc, toMatrix pcaComponents) of
    (Just s, Just vk) ->
      let mu = meanRow pcaMean
          xh = s LA.<> LA.tr vk
       in fromMatrix (xh + LA.asRow mu)
    _ -> error "Circuit.PCA.reconstruct: rank-2 expected"

-- | Full project-through-model: @reconstruct . scores@.
projectRows :: PCAModel -> Array Double -> Array Double
projectRows model x = reconstruct model (scores model x)

-- | Dagger-shaped Gram as an explicit reverse wire on rank-2 arrays.
--
-- @fwd@ multiplies on the left by @X†@ (features←samples); @bwd@ is the
-- dual map. Illustrates why 'Dagger' is in circuits: covariance is
-- reverse-then-compose, not a one-way arrow.
gramDagger :: Array Double -> Dagger (->) (Array Double) (Array Double)
gramDagger x =
  case toMatrix x of
    Nothing -> error "Circuit.PCA.gramDagger: expected rank-2"
    Just m ->
      let xt = LA.tr m
          fr a = case toMatrix a of
            Just ma -> fromMatrix (xt LA.<> ma)
            Nothing -> error "gramDagger.front: rank-2"
          ba b = case toMatrix b of
            Just mb -> fromMatrix (m LA.<> mb)
            Nothing -> error "gramDagger.back: rank-2"
       in Dagger fr ba

-- | Feature Gram via dagger-compose: @front d x@.
gramViaDagger :: Array Double -> Array Double
gramViaDagger x =
  let d = gramDagger x
   in front d x

-- | Product optic on one ambient row.
--
-- Focus = major scores (@k@); residual = reconstruction error in ambient
-- space (owned). @put@ replaces scores and re-adds the previous residual.
rowLens ::
  PCAModel ->
  Lens (VS.Vector Double) (VS.Vector Double) (VS.Vector Double) (VS.Vector Double)
rowLens PCAModel {pcaComponents, pcaMean} =
  fromClassical $ \row ->
    let meanV = A.asVector pcaMean
        xc = VS.zipWith (-) row meanV
     in case toMatrix pcaComponents of
          Nothing -> error "Circuit.PCA.rowLens: components not rank-2"
          Just vk ->
            let sc = LA.flatten (LA.asRow xc LA.<> vk)
                recon = VS.zipWith (+) meanV (LA.flatten (LA.asRow sc LA.<> LA.tr vk))
                err = VS.zipWith (-) row recon
             in ( sc,
                  \sc' ->
                    let recon' = VS.zipWith (+) meanV (LA.flatten (LA.asRow sc' LA.<> LA.tr vk))
                     in VS.zipWith (+) recon' err
                )

-- | View major scores of one row.
viewMajor :: PCAModel -> VS.Vector Double -> VS.Vector Double
viewMajor model row = view (rowLens model) row

-- | Set major scores of one row (residual error held from @row@).
setMajor :: PCAModel -> VS.Vector Double -> VS.Vector Double -> VS.Vector Double
setMajor model sc row = set (rowLens model) sc row
