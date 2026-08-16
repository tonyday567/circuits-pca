{-# LANGUAGE GADTs #-}
{-# LANGUAGE RankNTypes #-}

-- | Free-Tambara product optic packing.
--
-- The promoted canonical shape puts the base arrow first:
--
-- @
-- data Optic arr mon s u a b where
--   Optic :: arr s (mon m a) -> arr (mon m b) u -> Optic arr mon s u a b
-- @
--
-- * @mon@ is the monoidal action on the residual × interface.
-- * @s, u@ are the outer state endpoints (@s -> u@).
-- * @a, b@ are play / coplay (the interface).
--
-- This module is the function-arrow special case (@arr = (->)@):
--
-- @
-- Optic mon s u a b  ≅  ∃m. (s -> mon m a) × (mon m b -> u)
-- @
--
-- Provenance: circuits @examples/tambara.md@ / Milewski /Tambara Equipment/.
--
-- @
-- FreeTamb t j  →  j = Rep a b  →  Optic t  →  t = (,)  →  Lens
-- @
--
-- Residual @m@ is /owned data/ (hinge table cell A). PCA keeps the principal
-- summand as focus and the minor complement as residual — never seals it
-- with 'Circuit.Trace.trace'.
module Circuit.PCA.Optic
  ( Optic (..),
    Lens,
    OneShot,
    fromClassical,
    toClassical,
    view,
    set,
    over,

    -- * Adapter to polynomial morphisms
    morphismAsLens,
    lensAsMorphism,
  )
where

import Circuit.Poly (Mono, Morphism, applyLens, lens)

-- $setup
-- >>> import Circuit.Poly
-- >>> import Circuit.PCA.Optic

-- | Existential residual optic for monoidal action @mon@ over @(->)@.
--
-- The promoted/canonical shape is @Optic arr mon s u a b@; here the base
-- arrow is fixed to @(->)@.
--
-- @
-- Optic mon s u a b  ≅  ∃m. (s -> mon m a) × (mon m b -> u)
-- @
data Optic mon s u a b where
  Optic :: (s -> mon m a) -> (mon m b -> u) -> Optic mon s u a b

-- | Product-action optic = classical lens packing.
--
-- In classical notation the outer ends are named @s, t@ rather than @s, u@.
type Lens s t a b = Optic (,) s t a b

-- | One-shot product-residual optic.
--
-- This is the shape that a single state-changing morphism takes. It is /not/
-- the same as 'Circuit.Poly.System': a 'System' has a fixed carrier and is
-- iterable, whereas a 'OneShot' morphism may change its outer state @s -> u@.
type OneShot s u a b = Optic (,) s u a b

-- | Yoneda form @s -> (a, b -> t)@ into existential residual form.
fromClassical :: (s -> (a, b -> t)) -> Lens s t a b
fromClassical f =
  Optic
    (\s -> let (a, k) = f s in (k, a))
    (\(k, b) -> k b)

-- | Existential residual form into Yoneda form.
toClassical :: Lens s t a b -> s -> (a, b -> t)
toClassical (Optic get put) s =
  let (m, a) = get s
   in (a, \b -> put (m, b))

-- | Read the focus, discarding residual.
view :: Lens s s a a -> s -> a
view l s = fst (toClassical l s)

-- | Replace the focus, keeping residual from @s@.
set :: Lens s t a b -> b -> s -> t
set l b s = snd (toClassical l s) b

-- | Map the focus.
over :: Lens s t a b -> (a -> b) -> s -> t
over l f s =
  let (a, k) = toClassical l s
   in k (f a)

-- ---------------------------------------------------------------------------
-- Adapter: polynomial monomial morphism <-> classical state-preserving lens
-- ---------------------------------------------------------------------------

-- | A polynomial morphism @Mono s s -> Mono i o@ is exactly a state-preserving
-- classical lens @Lens s s o i@.
--
-- Both pack the same data: @s -> (o, i -> s)@.
--
-- >>> let m = lens (\s -> s + 1) (\s i -> s + i) :: Morphism (Mono Int Int) (Mono Int Int)
-- >>> view (morphismAsLens m) 5
-- 6
morphismAsLens :: Morphism (Mono s s) (Mono i o) -> Lens s s o i
morphismAsLens m = fromClassical $ \s ->
  let (o, put) = applyLens m s
   in (o, put)

-- | Inverse of 'morphismAsLens'.
--
-- >>> let l = fromClassical (\s -> (s * 2, \i -> s + i)) :: Lens Int Int Int Int
-- >>> let m = lensAsMorphism l :: Morphism (Mono Int Int) (Mono Int Int)
-- >>> let (o, put) = applyLens m 5 in (o, put 3)
-- (10,8)
lensAsMorphism :: Lens s s o i -> Morphism (Mono s s) (Mono i o)
lensAsMorphism l = lens get put
  where
    get s = fst (toClassical l s)
    put s i = snd (toClassical l s) i
