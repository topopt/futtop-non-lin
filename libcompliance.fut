import "src/utility"
import "src/solvers"
import "src/nonlinConjugateGradient"
import "libmultigrid"

entry computeCompliance f u =
     map2_4d (*) u f |> map_3d f64.sum |> map_2d f64.sum |> map f64.sum |> f64.sum

entry solveLinear [nelx][nely][nelz][nx][ny][nz] (xPhys :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) = 
    let f   = map_4d f32.f64 f
    let uold = map_4d f32.f64 uold
    let mgData = assembleMultigridData xPhys
    let (u,res,its) = cgSolveMG 1e-8 1000 xPhys mgData f uold
    let u = map_4d f64.f32 u
    let res = f64.f32 res
    in (u,res,1i32,its)

entry solveNonlinearLoadStep [nelx][nely][nelz][nx][ny][nz] (xPhys :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (loadscale :f64) = 
    let f = map_4d (*loadscale) f
    let (u,res,its,lits) = nonlinNewton 1e-8 xPhys f uold
    in (u,res,its,lits)