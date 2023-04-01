import "src/nonlinConjugateGradient"
import "src/utility"

-- entry main [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
--   let uold = map_4d f64.f32 uold
--   let f = map_4d f64.f32 f
--   in nonlinConjugateGradient x f uold
--
-- entry newton [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
--   let uold = map_4d f64.f32 uold
--   let f = map_4d f64.f32 f
--   in nonlinNewton x f uold

entry jacobi [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
  let uold = map_4d f64.f32 uold
  let f = map_4d f64.f32 f
  in nonlinConjugateGradientJacobi 1e-5 10000 x f uold

entry cg [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
  let uold = map_4d f64.f32 uold
  let f = map_4d f64.f32 f
  in nonlinConjugateGradientNoPre 1e-5 10000 x f uold

entry mix [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
  let u = map_4d f64.f32 u
  let f = map_4d f64.f32 f
  let (u, _, cgits) = nonlinConjugateGradientNoPre 1e-3 1000 x f u
  let (u, nres, nits) = nonlinNewton x f u
  in (u, nres, cgits, nits)

entry generate_test_input (nelx :i64) (nely :i64) (nelz :i64) =
  let nx = nelx+1
  let ny = nely+1
  let nz = nelz+1
  let x = replicate nelx (replicate nely (replicate nelz 1f32))
  let u = replicate nx (replicate ny (replicate nz [0f32,0,0]))
  let fval = -0.1
  let fvals = replicate ny [0,0,fval]
  let fidx = zip3 (replicate ny nelx) (iota ny) (replicate ny 0)
  let f = scatter_3d (copy u) fidx fvals
  let f = f with [nelx,0,0,2] = fval/2
  let f = f with [nelx,nely,0,2] = fval/2
  in (x, u, f)

entry generate_test_input_val (nelx :i64) (nely :i64) (nelz :i64) (fval :f32) =
  let nx = nelx+1
  let ny = nely+1
  let nz = nelz+1
  let x = replicate nelx (replicate nely (replicate nelz 1f32))
  let u = replicate nx (replicate ny (replicate nz [0f32,0,0]))
  let fvals = replicate ny [0,0,fval]
  let fidx = zip3 (replicate ny nelx) (iota ny) (replicate ny 0)
  let f = scatter_3d (copy u) fidx fvals
  let f = f with [nelx,0,0,2] = fval/2
  let f = f with [nelx,nely,0,2] = fval/2
  in (x, u, f)

-- ==
-- entry: jacobi
-- compiled script input { generate_test_input 8i64 8i64 8i64 }
-- compiled script input { generate_test_input 32i64 16i64 16i64 }
-- compiled script input { generate_test_input 64i64 32i64 32i64 }

-- ==
-- entry: cg
-- compiled script input { generate_test_input 8i64 8i64 8i64 }
-- compiled script input { generate_test_input 32i64 16i64 16i64 }
-- compiled script input { generate_test_input 64i64 32i64 32i64 }

-- ==
-- entry: mix
-- compiled script input { generate_test_input 8i64 8i64 8i64 }
-- compiled script input { generate_test_input 32i64 16i64 16i64 }
-- compiled script input { generate_test_input 64i64 32i64 32i64 }
