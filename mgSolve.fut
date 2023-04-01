import "libmultigrid"
import "src/explicitTimestepping"
import "src/utility"

def main [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
  let data = assembleMultigridData x
  in solveSystem x data uold f

-- entry diff [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
--   let (sol,_,_) = main x uold f
--   let (ex_sol,_) = explicitQuasiStaticSolve x f
--   let diff = map2_4d (\a b -> a-b) sol ex_sol
--   let diff_norm = map_4d (\v -> v*v) diff
--     |> flatten_3d
--     |> flatten
--     |> reduce f32.max 0
--   let true_norm = map_4d (\v -> v*v) sol
--     |> flatten_3d
--     |> flatten
--     |> reduce f32.max 0
--   in diff_norm / true_norm

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
