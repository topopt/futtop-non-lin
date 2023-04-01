import "src/assembleMass"
import "src/explicitTimestepping"
import "src/utility"
import "src/assembleResidual"
import "src/keUtilities"

entry assembleMass [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) =
  assembleLumpedMassMatrix x |> map_4d (f32.f64)

entry TimeStep [nx][ny][nz][nelx][nely][nelz] (totalSteps :i64) (currentStep :i64) (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (M :[nx][ny][nz][3]f32) (Minv :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) (uold :[nx][ny][nz][3]f32) =
  explicitTimeStep totalSteps currentStep x f M Minv u uold

entry TimeStepNL [nx][ny][nz][nelx][nely][nelz] (totalSteps :i64) (currentStep :i64) (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (M :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) (uold :[nx][ny][nz][3]f32) =
  let loadscale = f32.f64 (getCurrentLoadScale totalSteps currentStep)
  in (loadscale, explicitTimeStepNonlin loadscale x f M u uold)

entry quasiStaticSolve [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) =
  explicitQuasiStaticSolve x f

entry stoppingCriteria [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32)  (M :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) (uold :[nx][ny][nz][3]f32) :f32 =
  stoppingCriteria x f M u uold

def elementStress (u :[24]f32) =
  let iso = {xi=0,eta=0,zeta=0}
  let B0 = getB0 iso
  let BL = getBL iso u
  let eps_0 = vecmul_f32 B0 u
  let eps_L = vecmul_f32 BL u
  let eps = map2 (\e0 eL -> e0 + 0.5*eL) eps_0 eps_L
  let sigma = getStress eps
  in sigma

def elementVMstress [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (elementIndex :index) :f32 =
  let uloc = getLocalState_f32 u elementIndex
  let es = getElementScale_f32 x elementIndex
  let sigma = elementStress uloc
    |> map (*es)    -- scale for youngs module
  let s11 = sigma[0]
  let s22 = sigma[1]
  let s33 = sigma[2]
  let s12 = sigma[3]
  let s13 = sigma[4]
  let s23 = sigma[5]
  in 0.5* ((s11-s22)*(s11-s22) + (s22-s33)*(s22-s33) + (s11-s33)*(s11-s33)) + 3*(s12*s12+s13*s13+s23*s23)
    |> f32.sqrt

entry computeVMstress [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) =
  #[incremental_flattening(no_outer)] tabulate_3d nelx nely nelz (\i j k -> elementVMstress x u {x=i,y=j,z=k})
