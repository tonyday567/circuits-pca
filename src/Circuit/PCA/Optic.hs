{-# LANGUAGE GADTs #-}
{-# LANGUAGE RankNTypes #-}

-- | Free-Tambara product optic packing (card-local type promoted to API).
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
    fromClassical,
    toClassical,
    view,
    set,
    over,
  )
where

-- | Existential residual optic for monoidal action @t@.
--
-- @
-- Optic t s u a b  ≅  ∃m. (s → t m a) × (t m b → u)
-- @
data Optic t s u a b where
  Optic :: (s -> t m a) -> (t m b -> u) -> Optic t s u a b

-- | Product-action optic = classical lens packing.
type Lens s t a b = Optic (,) s t a b

-- | Yoneda form @s → (a, b → t)@ into existential residual form.
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
