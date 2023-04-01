import "src/explicitTimestepping"
import "src/utility"

def main [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
  let f = map_4d f64.f32 f
  in explicitQuasiStaticSolve x f

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

-- ==
-- entry: main
-- compiled script input { generate_test_input 8i64 8i64 8i64 }
-- compiled script input { generate_test_input 64i64 64i64 64i64 }
-- compiled script input { generate_test_input 128i64 128i64 128i64 }
