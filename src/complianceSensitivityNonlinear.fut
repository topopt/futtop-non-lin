import "indexUtilities"
import "nonlinearElement"
import "utility"
import "material"
import "nonlinConjugateGradient"


entry getElementSensitivity (x :f64) (esize :f64) (u :*[24]f64) (lambda :*[24]f64) =
  let young  = getYoungsModule x
  let Dyoung = getYoungsModuleDerivative x
  let frac   = getNonlinearFraction x
  let Dfrac  = getNonlinearFractionDerivative x
  let u_nl   = map (*frac) u
  let r_nl   = getNonlinearResidual esize u_nl
  let r_l    = keprod 1 (copy u)
  let r      = map2 (\f_nl f_l -> frac*f_nl + (1-frac*frac)*f_l) r_nl r_l
  let dr1 = map (*Dyoung) r
  let dr2 = map2 (\f_nl f_l -> young*Dfrac*(f_nl - 2*frac*f_l)) r_nl r_l
  let dr = map2 (+) dr1 dr2
  in map2 (*) lambda dr |> f64.sum |> (* -1) |> f32.f64

def getIndexComplianceAndSensitivity [nelx][nely][nelz][nx][ny][nz] (esize :f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (lambda :[nx][ny][nz][3]f64) (elementIndex :index) =
  let uloc   = getLocalState u elementIndex
  let xloc   = getDensityUnsafe x elementIndex
  let lambdaloc   = getLocalState lambda elementIndex
  in getElementSensitivity xloc esize uloc lambdaloc

def solveAdjointProblem = flexiblecgLinearised 1e-5 50 

entry getComplianceSensitivityNonlinear [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (lambda_old :*[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) =
  let esize = computeElementSize nely
  let (lambda,_) = solveAdjointProblem x f u lambda_old
  let c = map2_4d (*) f u |> map_3d f64.sum |> map_2d f64.sum |> map f64.sum |> f64.sum |> f32.f64
  let dc = #[incremental_flattening(no_outer)] tabulate_3d nelx nely nelz (\i j k -> getIndexComplianceAndSensitivity esize x u lambda {x=i,y=j,z=k})
  let dc = map_3d (f32.min 0) dc
  in (c, dc, lambda)

-- ==
-- entry: getComplianceSensitivity
-- nobench input @../testData/getComplianceSens.txt output @../testData/getComplianceSensOut.txt
-- compiled random input { [64][64][64]f32 [65][65][65][3]f32 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f32 }


-- entry getFDsens (u :[24]f64) (x :f64) (dx :f64) = 
--   let (c0,_) = getElementComplianceAndSensitivity (copy u) (x-dx)
--   let (c1,_) = getElementComplianceAndSensitivity (copy u) (x+dx)
--   let c0 = f64.f32 c0
--   let c1 = f64.f32 c1
--   in (c0-c1)/(2*dx)

-- entry FDcheck (dx :f64) = 
--   let u = [0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64,
--             0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 
--             0.0f64, 0.0f64, 0.0f64, 0.0f64, 1.0f64, 0.0f64]
--   let x = 0.5f64
--   let (_,dc) = getElementComplianceAndSensitivity (copy u) x
--   let dc_fd = getFDsens u x dx
--   let dc = f64.f32 dc
--   let diff = (dc - dc_fd) / dc_fd
--   in (dc_fd, dc, diff)