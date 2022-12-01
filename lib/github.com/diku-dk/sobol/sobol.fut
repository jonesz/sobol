-- | Quasi-random Sobol number generation.

module type sobol_dir = {
  val n: i64
  val k: i64
  val a: [n]u32
  val s: [n]i32
  val m: [n][k]u32
}

-- | Output type to be converted to.
module type sobol_output = {
  type t

  val i64: i64 -> t
  val u32: u32 -> t
  val f32: f32 -> t
  val **: t -> t -> t
  val /: t -> t -> t
}

module type sobol = {
  -- | Scalar type.
  type t
  -- | Dimensionality of sequence.
  val D : i64 
  -- | The value `2**32`.
  val norm : t
  -- | `independent i` returns the `i`'th sobol vector (in u32)
  val independent : i32 -> [D]u32
  -- | `recurrent i v` returns the `i`'th sobol vector given `v` is the
  -- `i-1`'th sobol vector
  val recurrent : i32 -> [D]u32 -> [D]u32
  -- | `chunk i n` returns the array `[v(i),...,v(i+n-1)]` of sobol
  -- vectors where `v(j)` is the `j`'th D-dimensional sobol vector
  val chunk : i32 -> (n:i64) -> [n][D]t
  --| `sobol n` generates `n` `D`-dimensional pointwise normalised
  -- sobol vectors.
  val sobol : (n:i64) -> [n][D]t
}

module Sobol (DM: sobol_dir) (X: { val D : i64 }) (T: sobol_output) : sobol with t = T.t = {
  let D = X.D
  type t = T.t

  -- Compute direction vectors. In general, some work can be saved if
  -- we know the number of sobol numbers (N) up front. Here, however,
  -- we calculate sufficiently sized direction vectors to work with
  -- upto N = 2^L, where L=32 (i.e., the maximum number of bits
  -- needed).

  let L = 32i64

  -- direction vector for dimension j
  let dirvec (j:i64) : [L]u32 =
    if j == 0 then
       map (\i -> 1u32 << (u32.i64 L-u32.i64 (i+1))) (iota L)
    else
       let s = i64.i32 DM.s[j-1]
       let a = DM.a[j-1]
       let V = map (\i -> if i >= s then 0u32
                          else DM.m[j-1,i] << (u32.i64 L-u32.i64(i+1))
                   ) (iota L)
       in loop (V) for i' < L-s do
            let i = i'+s
            let v = V[i-s]
            let vi0 = v ^ (v >> (u32.i64 s))
            let vi =
              loop vi = vi0 for k' < s-1 do
                let k = k'+1
                in vi ^ (((a >> u32.i64(s-1-k)) & 1u32) * V[i-k])
            in V with [i] = vi

  let index_of_least_significant_0 (x:i32) : i32 =
    loop i = 0 while i < 32 && ((x>>i)&1) != 0 do i + 1

  let norm = (T.f32 2.0) T.** T.i64 L

  let grayCode (x: i32): i32 = (x >> 1) ^ x

  let testBit (n: i32) (ind:i32) : bool =
    let t = (1 << ind) in (n & t) == t

  let dirvecs : [D][L]u32 =
    map dirvec (iota D)

  let recSob (i:i32) (dirvec:[L]u32) (x:u32) : u32 =
    (x ^ dirvec[index_of_least_significant_0 i])

  let recurrent (i:i32) (xs:[D]u32) : [D]u32 =
    map2 (recSob (i-1)) dirvecs xs

  let indSob (n:i32) (dirvec:[L]u32) : u32 =
    let reldv_vals = map2 (\dv i -> if testBit (grayCode n) (i32.i64 i) then dv
                                    else 0u32)
                         dirvec (iota L)
    in reduce (^) 0u32 reldv_vals

  let independent (i:i32) : [D]u32 =
    map (indSob i) dirvecs

  -- utils
  let recM (i:i32) : [D]u32 =
    let bit = index_of_least_significant_0 i
    in map (\row -> row[bit]) dirvecs

  -- computes sobol numbers: offs,..,offs+n-1
  let chunk (offs:i32) (n:i64) : [n][D]t =
    let sob_beg = independent offs
    let contrbs = map (\k: [D]u32 ->
                       if k==0 then sob_beg
                       else recM (i32.i64 k+offs-1))
                      (iota n)
    let vct_ints = scan (map2 (^)) (replicate D 0u32) contrbs
    in map (\xs -> map (\x -> (T.u32 x) T./ norm) xs)
           vct_ints

  let sobol (n:i64) : [n][D]t =
    tabulate n (i32.i64 >-> independent >-> map T.u32 >-> map (T./ norm))
}
