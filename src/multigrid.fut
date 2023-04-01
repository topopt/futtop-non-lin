import "smoothers"
import "applyStiffnessMatrix"
import "assembly"
import "assembleNonlinear"
import "projection"
import "utility"
import "sor"


-- type~ mgL2Data = ([][][][3]f64, [][][][3][81]f64)
-- type~ mgL1Data = [][][][3][3]f64
-- type~ mgL0Data = [][][][3][3]f64
-- type~ multigridData = (mgL0Data, mgL1Data, mgL2Data)

type~ mgL3Data = ([][][][3]f64, [][][][3][81]f64)
type~ mgL2Data = [][][][3][81]f64
type~ mgL1Data = [][][][3][81]f64
type~ mgL0Data = [][][][3][81]f64
type~ multigridData = (mgL0Data, mgL1Data, mgL2Data, mgL3Data)

def generateMultigridData [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :multigridData =
  -- let d0 = [[[replicate 3 (replicate 81 0)]]] -- dummy
  -- let d1 = [[[replicate 3 (replicate 81 0)]]] -- dummy
  -- let m2 = assembleStiffnessMatrix 2 x
  -- let d2 = extractInverseDiagonal m2
  -- in (d0, d0, (d2,m2))
  let m0 = assembleStiffnessMatrix 0 x
  let m1 = assembleStiffnessMatrix 1 x
  let m2 = assembleStiffnessMatrix 2 x
  let m3 = assembleStiffnessMatrix 3 x
  let d3 = extractInverseDiagonal m3
  in (m0, m1, m2, (d3,m3))

def generateMultigridDataNonlinear [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :multigridData =
  -- let d0 = getNonlinearBlockDiagonal x u
  -- let d1 = assembleBlockDiagonalNonlinear 1 x u
  -- let m2 = assembleStiffnessMatrixNonlinear 2 x u
  -- let d2 = extractInverseDiagonal m2
  -- in (d0, d1, (d2,m2))
  let m0 = assembleStiffnessMatrixNonlinear 0 x u
  let m1 = assembleStiffnessMatrixNonlinear 1 x u
  let m2 = assembleStiffnessMatrixNonlinear 2 x u
  let m3 = assembleStiffnessMatrixNonlinear 3 x u
  let d3 = extractInverseDiagonal m3
  in (m0, m1, m2, (d3,m3))

def innerProduct [n][m][l][k] (a :[n][m][l][k]f64) (b :[n][m][l][k]f64) :f64 =
  map2_4d (*) a b
  |> map_3d f64.sum
  |> map_2d f64.sum
  |> map    f64.sum
  |>        f64.sum

def norm [n][m][l][k] (a :[n][m][l][k]f64) :f64 =
  map_4d (\x -> x*x) a
  |> map_3d f64.sum
  |> map_2d f64.sum
  |> map    f64.sum
  |>        f64.sum
  |> f64.sqrt

-- def cgSolveJacSubspace [nx][ny][nz] (data :multigridData) (b :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
--   let blockIt = 50
--   let maxIt   = 800
--   let tol     = 1e-10
--   let invD    = data.2.0 :> [nx][ny][nz][3]f64
--   let matrix  = data.2.1 :> [nx][ny][nz][3][81]f64
--   let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))
--   let bnorm = norm b

--   let inner_iteration (uold, rold, pold, rhoold) = 
--     let z      = map2_4d (*) invD rold -- jacobi preconditioner
--     let rho    = innerProduct rold z
--     let beta   = rho / rhoold
--     let p      = map2_4d (\pp zz -> beta * pp + zz) pold z
--     let q      = applyAssembledStiffnessMatrix matrix p
--     let alpha  = rho / (innerProduct p q)
--     let u      = map2_4d (\uu pp -> uu + alpha * pp) uold p
--     let r      = map2_4d (\rr qq -> rr - alpha * qq) rold q
--     in (u, r, p, rho)

--   let  (u, _, _, _, res, _) = 
--     loop (u, r, p, rho, res, it) = (zero_4d, b, zero_4d, 1f64, 1, 0) while res > tol && it < maxIt do
--       let (u, r, p, rho) = (iterate blockIt inner_iteration) (u, r, p, rho)
--       let res = (norm r) / bnorm
--       in (u, r, p, rho, res, it+blockIt)

--   let res = #[trace(cres)] res/tol
--   let u = map_4d (*res) u
--   let u = map_4d (/res) u
--   in u


def bicgstabSubspace [nx][ny][nz] (data :multigridData) (b :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let blockIt = 25
  let maxIt   = 4000
  let tol     = 1e-10
  let invD    = data.3.0 :> [nx][ny][nz][3]f64
  let matrix  = data.3.1 :> [nx][ny][nz][3][81]f64
  let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))
  let bnorm = norm b
  let dnorm = norm (map_4d (\v->1/v) invD)

  let scale =  dnorm/bnorm
  let inverse_scale = 1/scale
  let b = map_4d (*scale) b

  let r0hat = b
  let applyPrec = map2_4d (*) invD 
  let applyMat  = applyAssembledStiffnessMatrix matrix

  let inner_iteration (xold, rold, pold, vold, rhoold, omegaold, alphaold) = 
    let rho    = innerProduct r0hat rold
    let beta   = (rho/rhoold) * (alphaold / omegaold)
    let p      = map3_4d (\ro po vo -> ro + beta*(po-omegaold*vo)) rold pold vold
    let y      = applyPrec p -- jacobi preconditioner
    let v      = applyMat y
    let alpha  = rho / (innerProduct r0hat v)
    let s      = map2_4d (\rr vv -> rr - alpha * vv) rold v
    let z      = applyPrec p -- jacobi preconditioner
    let t      = applyMat z
    let omega  = (innerProduct t s) / (innerProduct t t)
    let x      = map3_4d (\xo zz yy -> xo + omega * zz + alpha * yy) xold z y
    let r      = map2_4d (\ss tt -> ss - omega * tt) s t
    in (x, r, p, v, rho, omega, alpha)

  let  (u,_,_,_,_,_,_,res,_) = 
    loop (u, r, p, v, rho, omega, alpha, res, it) = (zero_4d, b, zero_4d, zero_4d, 1f64, 1f64, 1f64, 1, 0) while res > tol && it < maxIt do
      let (u, r, p, v, rho, omega, alpha) = (iterate blockIt inner_iteration) (u, r, p, v, rho, omega, alpha)
      let res = (norm r) * (inverse_scale / bnorm)
      in (u, r, p, v, rho, omega, alpha, res, it+blockIt)

  -- let res = #[trace(cres)] res/tol
  -- let u = map_4d (*res) u
  -- let u = map_4d (/res) u
  in map_4d (*inverse_scale) u

def vcycle_l2 [nx][ny][nz] (data :multigridData) (f :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let matrix = data.2 :> [nx][ny][nz][3][81]f64
  let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))
  let z = sorAssembled matrix f zero_4d
  let v = z
    |> applyAssembledStiffnessMatrix matrix
    |> map2_4d (\ff dd -> ff - dd) f
    |> projectToCoarser
    |> bicgstabSubspace data
    |> projectToFiner
  let z = map2_4d (+) z (v :> [nx][ny][nz][3]f64)
  in sorAssembled matrix f z

def vcycle_l1 [nx][ny][nz] (data :multigridData) (f :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let matrix = data.1 :> [nx][ny][nz][3][81]f64
  let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))
  let z = sorAssembled matrix f zero_4d
  let v = z
    |> applyAssembledStiffnessMatrix matrix
    |> map2_4d (\ff dd -> ff - dd) f
    |> projectToCoarser
    |> vcycle_l2 data
    |> projectToFiner
  let z = map2_4d (+) z (v :> [nx][ny][nz][3]f64)
  in sorAssembled matrix f z

def vcycle_l0 [nelx][nely][nelz][nx][ny][nz] (data :multigridData) (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) :[nx][ny][nz][3]f32 =
  let f = map_4d f64.f32 f
  let zero_4d = replicate nx (replicate ny (replicate nz [0,0,0f64]))
  let matrix = data.0 :> [nx][ny][nz][3][81]f64
  let z = sorAssembled matrix f zero_4d
  let v = z
    |> applyAssembledStiffnessMatrix matrix
    |> map2_4d (\ff dd -> ff - dd) f -- compute r
    |> projectToCoarser
    |> vcycle_l1 data
    |> projectToFiner
  let z = map2_4d (+) z (v :> [nx][ny][nz][3]f64)
  in sorAssembled matrix f z |> map_4d f32.f64


def vcycle_l0_nonlinear [nelx][nely][nelz][nx][ny][nz] (data :multigridData) (x :[nelx][nely][nelz]f32) (u0 :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let zero_4d = replicate nx (replicate ny (replicate nz [0,0,0f64]))
  let matrix = data.0 :> [nx][ny][nz][3][81]f64
  let z = sorAssembled matrix f zero_4d
  let v = z
    |> applyAssembledStiffnessMatrix matrix
    |> map2_4d (\ff dd -> ff - dd) f -- compute r
    |> projectToCoarser
    |> vcycle_l1 data
    |> projectToFiner
  let z = map2_4d (+) z (v :> [nx][ny][nz][3]f64)
  in sorAssembled matrix f z
