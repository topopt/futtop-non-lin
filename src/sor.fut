import "utility"
import "indexUtilities"
import "assembly"
import "boundaryConditions"
import "smoothers"
import "material"
import "applyStiffnessMatrix"
import "keConstants"
import "nonlinearElement"
import "assemblyUtilities"
import "assembleNonlinear"

def omega         :f64 = 0.6

def sorLocalMatrix (f :[3]f64) (S :[3]f64) (M :[3][3]f64) (u :[3]f64) =
  let rx     = #[unsafe] M[0,1]*u[1] + M[0,2]*u[2]
  let ux_new = #[unsafe] (1/M[0,0]) * (f[0]-S[0]-rx)
  let ry     = #[unsafe] M[1,0]*ux_new + M[1,2]*u[2]
  let uy_new = #[unsafe] (1/M[1,1]) * (f[1]-S[1]-ry)
  let rz     = #[unsafe] M[2,0]*ux_new + M[2,1]*uy_new
  let uz_new = #[unsafe] (1/M[2,2]) * (f[2]-S[2]-rz)
  let unew = [ux_new, uy_new, uz_new]
  in #[sequential] map2 (\un uo -> omega*un + (1-omega)*uo) unew u

def sorLocalMatrixBack (f :[3]f64) (S :[3]f64) (M :[3][3]f64) (u :[3]f64) =
  let rz     = #[unsafe] M[2,0]*u[0] + M[2,1]*u[1]
  let uz_new = #[unsafe] (1/M[2,2]) * (f[2]-S[2]-rz)
  let ry     = #[unsafe] M[1,0]*u[0] + M[1,2]*uz_new
  let uy_new = #[unsafe] (1/M[1,1]) * (f[1]-S[1]-ry)
  let rx     = #[unsafe] M[0,1]*uy_new + M[0,2]*uz_new
  let ux_new = #[unsafe] (1/M[0,0]) * (f[0]-S[0]-rx)
  let unew = [ux_new, uy_new, uz_new]
  in #[sequential] map2 (\un uo -> omega*un + (1-omega)*uo) unew u

def sorForwardAssembled (mat :[3][81]f64) (uStencil :*[81]f64) (f :[3]f64) :*[3]f64 =
  let uold = #[unsafe] [uStencil[39], uStencil[40], uStencil[41]]
  let uStencil = uStencil with [39] = 0
  let uStencil = uStencil with [40] = 0
  let uStencil = uStencil with [41] = 0
  let S = map (\row -> (f64.sum (map2 (*) uStencil row))) mat
  let M = mat[:3,39:42] :> [3][3]f64
  in sorLocalMatrix f S M uold

def sorBackwardAssembled (mat :[3][81]f64) (uStencil :*[81]f64) (f :[3]f64) :*[3]f64 =
  let uold = #[unsafe] [uStencil[39], uStencil[40], uStencil[41]]
  let uStencil = uStencil with [39] = 0
  let uStencil = uStencil with [40] = 0
  let uStencil = uStencil with [41] = 0
  let S = map (\row -> (f64.sum (map2 (*) uStencil row))) mat
  let M = mat[:3,39:42] :> [3][3]f64
  in sorLocalMatrixBack f S M uold

def ssorSweepAssembled [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let uhalf = tabulate_3d nx ny nz (\i j k ->
    let uloc = getLocalU u {x=i,y=j,z=k}
    in #[unsafe] sorForwardAssembled mat[i,j,k] uloc f[i,j,k])
  let unew = tabulate_3d nx ny nz (\i j k ->
    let uloc = getLocalU uhalf {x=i,y=j,z=k}
    in #[unsafe] sorBackwardAssembled mat[i,j,k] uloc f[i,j,k])
  in unew

entry sorAssembled [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let number_sweeps :i32 = 4
  in (iterate number_sweeps (ssorSweepAssembled mat f)) u

-- let generateNodeOffsetWithoutCenter (eo :index) =
--   let no = getNodeIndices eo
--   let no = filter (\i -> i.x != 0 || i.y != 0 || i.z != 0) no
--   let no = no :> [7]index
--   in zip (replicate 7 eo) no
--
-- let allOffsetPairs =
--   [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
--    {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
--    {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
--   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}] |> map generateNodeOffsetWithoutCenter |> flatten_to (7*8)

def allOffsetPairs = [
({x = 0i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 1i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 1i64, z = 1i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 0i64, z = 1i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 0i64, y = 0i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64},{x = -1i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = 0i64},{x = -1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = -1i64, y = 1i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64},{x = 0i64, y = 0i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = -1i64, y = 0i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = -1i64, z = 0i64}, {x = 1i64, y = -1i64, z = 0i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 0i64, y = -1i64, z = 0i64}),
({x = 0i64, y = -1i64, z = 0i64}, {x = 0i64, y = 0i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 1i64, y = 0i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64}, {x = 1i64, y = -1i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 0i64, y = -1i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = 0i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = 0i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = 0i64, y = 0i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = 0i64, y = -1i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = -1i64, z = 1i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64},{x = 1i64, y = 1i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 1i64, y = 0i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64},{x = 0i64, y = 0i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = -1i64},{x = 1i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = 0i64, z = -1i64},{x = -1i64, y = 1i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = 0i64, y = 0i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = -1i64, y = 0i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = -1i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = -1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 0i64, y = 0i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 1i64, y = 0i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 1i64, y = -1i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = -1i64, z = -1i64},{x = 1i64, y = -1i64, z = 0i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = 0i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = 0i64, y = 0i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = -1i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = -1i64, z = 0i64})]

def getSContribution [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (nodeIndex :index) (elementOffset :index, nodeOffset :index) :[3]f64 =
  let elementIndex = addIndices nodeIndex elementOffset
  let elementScale = getElementYoungsModule x elementIndex
  let localMatrix  = getLocalMatrix(elementOffset,nodeOffset)
  let ulocal       = getInputVectorSingle(nodeIndex,nodeOffset,u)
  let ulocal       = #[sequential] map f64.f32 ulocal
  in multiplyScaledLocalMatrix(localMatrix, ulocal, elementScale)

def build_S [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (nodeIndex :index) :[3]f64 =
  copy allOffsetPairs
  |> (#[sequential] map (getSContribution x u nodeIndex))
  |> transpose
  |> (#[sequential] map f64.sum)

def getOwnLocalMatrix (elementOffset :index) =
  let recievingNodeOffset :index = {x=(-elementOffset.x), y=(-elementOffset.y), z=(-elementOffset.z)}
  let li = getLocalNodeIndex(recievingNodeOffset)
  in getke_l0 (li,li)

def getMContribution [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (nodeIndex :index) (elementOffset :index) =
  let elementIndex = addIndices nodeIndex elementOffset
  let elementScale = getElementYoungsModule x elementIndex
  let localMatrix  = getOwnLocalMatrix elementOffset
  in scaleLocalMatrix localMatrix elementScale

def build_M [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) esize (nodeIndex :index) =
  let m = [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
   {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
   {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}]
   |> (#[sequential] map (getMContribution x nodeIndex))
   |> (#[sequential] reduce addLocalMatrix {xx=0,xy=0,xz=0,yx=0,yy=0,yz=0,zx=0,zy=0,zz=0})
  in scaleLocalMatrix m esize

def sorNodeFree  [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) esize (u :[nx][ny][nz][3]f32) (f :[3]f32) (nodeIndex :index) :[3]f32 =
 -- extract value of own node, and zero
 let ni = nodeIndex
 let ux_old = f64.f32 u[ni.x,ni.y,ni.z,0]
 let uy_old = f64.f32 u[ni.x,ni.y,ni.z,1]
 let uz_old = f64.f32 u[ni.x,ni.y,ni.z,2]

 let f = #[sequential] map f64.f32 f
 let S = build_S x u ni
 let M = build_M x esize ni

 let rx = M.xy*uy_old + M.xz*uz_old
 let ux_new = (1/M.xx) * (f[0]-S[0]-rx)
 let ry = M.yx*ux_new + M.yz*uz_old
 let uy_new = (1/M.yy) * (f[1]-S[1]-ry)
 let rz = M.zx*ux_new + M.zy*uy_new
 let uz_new = (1/M.zz) * (f[2]-S[2]-rz)

 let uold = [ux_old, uy_old, uz_old]
 let unew = [ux_new, uy_new, uz_new]
 let usmoothed = #[sequential] map2 (\un uo -> omega*un + (1-omega)*uo) unew uold
 in #[sequential] map f32.f64 usmoothed

def sorSweepFree [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
  let esize = computeElementSize nely
  in tabulate_3d nx ny nz (\i j k ->
    sorNodeFree x esize u f[i,j,k] {x=i,y=j,z=k})

entry sorMatrixFree [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
  let number_sweeps :i32 = 2
  let usmoothed = (iterate number_sweeps (sorSweepFree x f)) u
  in setBCtoInput u usmoothed

-- entry jacobiReference [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
--  let invD = assembleInverseDiagonal 0 x |> map_4d f32.f64 :> [nx][ny][nz][3]f32
--  in jacobiSmootherSingle x b invD u

-- ==
-- entry: sorMatrixFree
-- nobench input @../testData/sor1.txt auto output
-- nobench input @../testData/sor1.txt output @../testData/sor1Out.txt
-- nobench input @../testData/sor2.txt output @../testData/sor2Out.txt
-- compiled random input { [64][64][64]f32 [65][65][65][3]f32 [65][65][65][3]f32 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f32 [129][129][129][3]f32 }
-- compiled random input { [256][256][256]f32 [257][257][257][3]f32 [257][257][257][3]f32 }

def sorNonlinearForward [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) = 
  let Ku = applyLinearisedStiffnessMatrix x u0 u
  let Kd = map2_3d vecmul_f64 bdiag u
  let S  = map2_4d (-) Ku Kd
  let unew = map4_3d sorLocalMatrix f S bdiag u
  in setBCtoZero 0 unew

def sorNonlinearBackward [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) = 
  let Ku = applyLinearisedStiffnessMatrix x u0 u
  let Kd = map2_3d vecmul_f64 bdiag u
  let S  = map2_4d (-) Ku Kd
  let unew = map4_3d sorLocalMatrixBack f S bdiag u
  in setBCtoZero 0 unew

def sorNonlinearSweep [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) = 
  let uhalf = sorNonlinearForward x bdiag u0 f u
  in sorNonlinearBackward x bdiag u0 f uhalf

entry sorNonlinear [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let number_sweeps :i32 = 2
  in (iterate number_sweeps (sorNonlinearSweep x bdiag u0 f)) u

-- entry sorNLtest [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) = 
--   let f = map_4d f64.f32 f
--   let u = map_4d f64.f32 u
--   let u0 = replicate nx (replicate ny (replicate nz [0,0,0f64]))
--   let diag = #[trace(M)] getNonlinearBlockDiagonal x u0
--   in sorNonlinear x diag u0 f u |> map_4d f32.f64

-- -- ==
-- -- entry: sorNLtest
-- -- nobench input @../testData/sor1.txt auto output
-- -- nobench input @../testData/sor1.txt output @../testData/sor1Out.txt
-- -- nobench input @../testData/sor2.txt output @../testData/sor2Out.txt

def sorNonlinearForwardLevel1 [nelx][nely][nelz][nx][ny][nz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (bdiag :[nxc][nyc][nzc][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nxc][nyc][nzc][3]f64) (u :[nxc][nyc][nzc][3]f64) = 
  let Ku = applyLinearisedStiffnessMatrixLevel1 x u0 u
  let Kd = map2_3d vecmul_f64 bdiag u
  let S  = map2_4d (-) Ku Kd
  let unew = map4_3d sorLocalMatrix f S bdiag u
  in setBCtoZero 0 unew

def sorNonlinearBackwardLevel1 [nelx][nely][nelz][nx][ny][nz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (bdiag :[nxc][nyc][nzc][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nxc][nyc][nzc][3]f64) (u :[nxc][nyc][nzc][3]f64) = 
  let Ku = applyLinearisedStiffnessMatrixLevel1 x u0 u
  let Kd = map2_3d vecmul_f64 bdiag u
  let S  = map2_4d (-) Ku Kd
  let unew = map4_3d sorLocalMatrixBack f S bdiag u
  in setBCtoZero 0 unew

def sorNonlinearSweepLevel1 [nelx][nely][nelz][nx][ny][nz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (bdiag :[nxc][nyc][nzc][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nxc][nyc][nzc][3]f64) (u :[nxc][nyc][nzc][3]f64) = 
  let uhalf = sorNonlinearForwardLevel1 x bdiag u0 f u
  in sorNonlinearBackwardLevel1 x bdiag u0 f uhalf

def sorNonlinearLevel1 [nelx][nely][nelz][nx][ny][nz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (bdiag :[nxc][nyc][nzc][3][3]f64) (u0 :[nx][ny][nz][3]f64) (f :[nxc][nyc][nzc][3]f64) (u :[nxc][nyc][nzc][3]f64) =
  let number_sweeps :i32 = 2
  in (iterate number_sweeps (sorNonlinearSweepLevel1 x bdiag u0 f)) u