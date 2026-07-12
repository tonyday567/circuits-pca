# circuits-pca

Principal component analysis as a **residual-ownership** protocol on the
circuits stack, with **harpie** arrays and **hmatrix** BLAS/SVD underneath.

Design notes live in the parent library:

- [`circuits/examples/pca.md`](../circuits/examples/pca.md) — hinge table, why `Dagger` is here  
- [`circuits/examples/tambara.md`](../circuits/examples/tambara.md) — where `Optic (,)` comes from  

## hinge (cell A)

```text
                 residual ownership
              own ──────────────── seal
    finite    Optic / PCA            Arr / open wire
    ∞         coinductive optic      Knot / Hyper
```

PCA **owns** a finite residual (minor axes / reconstruction error). It does
not `trace`-seal feedback. Spectral work is an interpreter; the packing is
the free-Tambara product optic.

## pipeline

```text
data (samples × features)
  │ centerColumns
  ▼
centered X
  │ thin SVD  (same axes as eigendecomp of X†X)
  ▼
loadings V_k , singular values
  │ Optic (,)
  ▼
focus = scores (major)     residual = minor / error (owned)
```

Dagger-shaped step available explicitly as `gramFeatures` = feature Gram
`X†X` (transpose compose). Fit uses SVD of centered data (more stable,
same subspace).

## install / build

Local monorepo style (`cabal.project` pins siblings):

```bash
cd ~/haskell/circuits-pca
cabal build all
cabal test
cabal repl circuits-pca
```

Dependencies (local packages in `cabal.project`):

| package | role |
|---------|------|
| `circuits` | `Dagger` / traced vocabulary (protocol home) |
| `circuits-ad` | reverse-mode sibling (diff-through-PCA later) |
| `numhask` 0.14 | numeric hierarchy |
| `harpie` | array API |
| `harpie-numhask` | numhask orphans / fixed-shape bridge |
| `hmatrix` | BLAS mult + thin SVD (via `Harpie.Hmatrix`) |

## usage

```haskell
import Circuit.PCA
import Harpie.Array.Storable qualified as A

x :: A.Array Double
x = A.array [4, 2] ([1, 1, 2, 2, 3, 3, 4, 4] :: [Double])

model = fit 1 x
sc    = scores model x          -- 4×1
xh    = projectRows model x     -- 4×2 reconstruction
```

Single-row optic (major scores focus, residual error owned):

```haskell
import Data.Vector.Storable qualified as VS

row = VS.fromList [2, 2]
maj = viewMajor model row
row' = setMajor model maj row   -- ≈ row when row is in the span
```

## performance

Hot path is hmatrix/BLAS (`multM`, `thinSVD`). Protocol types (`Optic`,
`PCAModel`) are thin. Meter later with `circuits-meter` / `harpie-perf`
patterns (`~/haskell/harpie-perf`).

## status

0.1 — batch PCA, product optic, smoke tests. Not yet: streaming PCA,
diff through SVD, fixed-shape harpie-numhask path, randomised SVD.
