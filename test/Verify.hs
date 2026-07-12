-- | Smoke tests for circuits-pca.
module Main where

import Circuit.PCA
import Data.Vector.Storable qualified as VS
import Harpie.Array.Storable qualified as A
import Numeric.LinearAlgebra qualified as LA
import System.Exit (exitFailure, exitSuccess)
import Prelude

assert :: String -> Bool -> IO ()
assert msg True = putStrLn ("ok  — " ++ msg)
assert msg False = do
  putStrLn ("FAIL — " ++ msg)
  exitFailure

-- | Simple 2D cloud elongated on x = y.
--
-- Four points roughly along the diagonal.
sample2d :: A.Array Double
sample2d =
  A.array
    [4, 2]
    ( [ 1,
        1,
        2,
        2,
        3,
        3,
        4,
        4
      ] ::
        [Double]
    )

almostEq :: Double -> Double -> Bool
almostEq a b = abs (a - b) < 1e-6

almostEqV :: VS.Vector Double -> VS.Vector Double -> Bool
almostEqV a b =
  VS.length a == VS.length b
    && VS.and (VS.zipWith almostEq a b)

main :: IO ()
main = do
  putStrLn "circuits-pca verify"

  let (xc, mu) = centerColumns sample2d
  assert "mean shape [2]" (A.shape mu == A.shape (A.array [2] ([0, 0] :: [Double])))
  -- mean should be (2.5, 2.5)
  let meanV = A.asVector mu
  assert "column mean ~ 2.5" (almostEq (meanV VS.! 0) 2.5 && almostEq (meanV VS.! 1) 2.5)

  let g = gramFeatures xc
  assert "gram is 2×2" (case toMatrix g of Just m -> LA.rows m == 2 && LA.cols m == 2; _ -> False)

  let model = fit 1 sample2d
  assert "k=1" (pcaK model == 1)
  assert "one singular value" (VS.length (pcaSingularValues model) == 1)

  let sc = scores model sample2d
  assert "scores shape 4×1" $
    case toMatrix sc of
      Just m -> LA.rows m == 4 && LA.cols m == 1
      _ -> False

  let xh = projectRows model sample2d
  assert "project shape 4×2" $
    case toMatrix xh of
      Just m -> LA.rows m == 4 && LA.cols m == 2
      _ -> False

  -- First PC should align with the diagonal (≈ equal components)
  case toMatrix (pcaComponents model) of
    Just vk -> do
      let a = vk LA.! 0 LA.! 0
          b = vk LA.! 1 LA.! 0
      assert "PC1 ~ equal |loadings|" (almostEq (abs a) (abs b))
    Nothing -> assert "components matrix" False

  -- row lens: view major, set major preserves residual structure
  let row0 = VS.fromList [1, 1]
      maj = viewMajor model row0
      row0' = setMajor model maj row0
  assert "setMajor . viewMajor ~ id on diagonal row" (almostEqV row0 row0')

  -- optic laws smoke: set then view
  let maj2 = VS.map (* 2) maj
      row2 = setMajor model maj2 row0
      maj2' = viewMajor model row2
  assert "view . set ~ id on major" (almostEqV maj2 maj2')

  putStrLn "all ok"
  exitSuccess
