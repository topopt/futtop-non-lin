-- n-dimensional map
def map_2d 'a 'x [n][m] (f: a -> x) (as: [n][m]a): [n][m]x = map (map f) as
def map_3d 'a 'x [n][m][l] (f: a -> x) (as: [n][m][l]a): [n][m][l]x = map (map_2d f) as
def map_4d 'a 'x [n][m][l][k] (f: a -> x) (as: [n][m][l][k]a): [n][m][l][k]x = map (map_3d f) as

-- n-dimensional map2
def map2_2d 'a 'b 'x [n][m] (f: a -> b -> x) (as: [n][m]a) (bs: [n][m]b): [n][m]x =
  map2 (map2 f) as bs
def map2_3d 'a 'b 'x [n][m][l] (f: a -> b -> x) (as: [n][m][l]a) (bs: [n][m][l]b): [n][m][l]x =
  map2 (map2_2d f) as bs
def map2_4d 'a 'b 'x [n][m][l][k] (f: a -> b -> x) (as: [n][m][l][k]a) (bs: [n][m][l][k]b): [n][m][l][k]x =
  map2 (map2_3d f) as bs

-- n-dimensional map3
def map3_3d 'a 'b 'c 'x [n][m][l] (f: a -> b -> c -> x) (as: [n][m][l]a) (bs: [n][m][l]b) (cs: [n][m][l]c) : [n][m][l]x =
  map3 (map3 (map3 f)) as bs cs
def map3_4d 'a 'b 'c 'x [n][m][l][k] (f: a -> b -> c -> x) (as: [n][m][l][k]a) (bs: [n][m][l][k]b) (cs: [n][m][l][k]c) : [n][m][l][k]x =
  map3 (map3_3d f) as bs cs

-- 4-dimensional map4
def map4_2d 'a 'b 'c 'd 'x [n][m] (f: a -> b -> c -> d -> x) (as: [n][m]a) (bs: [n][m]b) (cs: [n][m]c) (ds: [n][m]d): [n][m]x =
  map4 (map4 f) as bs cs ds
def map4_3d 'a 'b 'c 'd 'x [n][m][l] (f: a -> b -> c -> d -> x) (as: [n][m][l]a) (bs: [n][m][l]b) (cs: [n][m][l]c) (ds: [n][m][l]d): [n][m][l]x =
  map4 (map4_2d f) as bs cs ds
def map4_4d 'a 'b 'c 'd 'x [n][m][l][k] (f: a -> b -> c -> d -> x) (as: [n][m][l][k]a) (bs: [n][m][l][k]b) (cs: [n][m][l][k]c) (ds: [n][m][l][k]d): [n][m][l][k]x =
  map4 (map4_3d f) as bs cs ds


def replicate_2d a b v = replicate a (replicate b v)
def replicate_3d a b c v = replicate a (replicate_2d b c v)
def replicate_4d a b c d v = replicate a (replicate_3d b c d v)

def matmul [n][m][p] 'a (add: a -> a -> a) (mul: a -> a -> a) (zero: a) (A: [n][m]a) (B: [m][p]a) : [n][p]a =
  map (\A_row ->
         map (\B_col ->
                reduce add zero (map2 mul A_row B_col))
             (transpose B))
      A

def vecmul [n][m] 'a (add: a -> a -> a) (mul: a -> a -> a) (zero: a) (A: [n][m]a) (x: [m]a) : [n]a =
  map (\A_row -> reduce add zero (map2 mul A_row x)) A

def matmul_f64 = matmul (+) (*) 0f64
def vecmul_f64 = vecmul (+) (*) 0f64
def matmul_f32 = matmul (+) (*) 0f32
def vecmul_f32 = vecmul (+) (*) 0f32

def matsum_f64 [n][m] (A: [n][m]f64) (B: [n][m]f64) : [n][m]f64 =
  map2 (\A_row B_row -> map2 (+) A_row B_row) A B

def matsum_f32 [n][m] (A: [n][m]f32) (B: [n][m]f32) : [n][m]f32 =
  map2 (\A_row B_row -> map2 (+) A_row B_row) A B

def outer_product [n][m] 'a (mul: a -> a -> a) (A: [n]a) (B: [m]a) : [n][m]a =
  map (\A_elem -> map (\B_elem -> mul A_elem B_elem) B) A

def det33_f64 (m :[3][3]f64) =
  let a = #[unsafe] m[0,0]
  let b = #[unsafe] m[0,1]
  let c = #[unsafe] m[0,2]
  let d = #[unsafe] m[1,0]
  let e = #[unsafe] m[1,1]
  let f = #[unsafe] m[1,2]
  let g = #[unsafe] m[2,0]
  let h = #[unsafe] m[2,1]
  let i = #[unsafe] m[2,2]
  in a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)

def inv33_f64 (m :[3][3]f64) :[3][3]f64 =
  let a = #[unsafe] m[0,0]
  let b = #[unsafe] m[0,1]
  let c = #[unsafe] m[0,2]
  let d = #[unsafe] m[1,0]
  let e = #[unsafe] m[1,1]
  let f = #[unsafe] m[1,2]
  let g = #[unsafe] m[2,0]
  let h = #[unsafe] m[2,1]
  let i = #[unsafe] m[2,2]
  let det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
  let invdet = 1/det
  in [[e*i-f*h, c*h-b*i, b*f-c*e],
      [f*g-d*i, a*i-c*g, c*d-a*f],
      [d*h-e*g, b*g-a*h, a*e-b*d]]
      |> map_2d (*invdet)
