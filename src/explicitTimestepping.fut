import "applyStiffnessMatrix"
import "assembleMass"
import "utility"
import "assembleNonlinear"


-- https://en.wikipedia.org/wiki/Smoothstep#7th-order_equation
def smootherstep (x :f64) :f64 =
  let lim :f64 = 1.0
  in if (x > lim)
    then 1.0
  else
    let xs = x / lim
    in xs*xs*xs*(10-15*xs+6*xs*xs)

def getCurrentLoadScale (totalSteps :i64) (currentStep :i64) :f64 =
  let fictuousTime = (f64.i64 currentStep) / (f64.i64 totalSteps)
  in smootherstep fictuousTime

def alpha = 0.001f64
def beta  = 0.001f64
def dt    = 0.65f64 --0.65 @ M=1

entry explicitTimeStepNonlin [nx][ny][nz][nelx][nely][nelz] (loadScale :f64) (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f64) (M :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) (uold :[nx][ny][nz][3]f64) =
  let VelocityNhalf = map2_4d (\d dold -> (1/dt)*(d-dold)) u uold
  let Rext = map_4d (*loadScale) f
  let residualrhs = map2_4d (\uu vv -> uu + beta*vv) u VelocityNhalf
  let Rint =  assembleNonlinearResidual x residualrhs
  let massDamping = map2_4d (*) M VelocityNhalf |> map_4d (*alpha)
  let accelerationRhs = map2_4d (\uu vv -> uu+dt*vv) u VelocityNhalf
  let accelerationScale = 1/(dt*dt)
  let accelerationTerm = map2_4d (*) M accelerationRhs
    |> map_4d (*accelerationScale)
  let rhs = map4_4d (\re ri mt ct -> re-ri+mt-ct) Rext Rint accelerationTerm massDamping
  in map2_4d (\m r -> (dt*dt*r)/m) M rhs

def stoppingCriteria [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f64)  (M :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) (uold :[nx][ny][nz][3]f64) :f64 =
  let r = assembleNonlinearResidual x u
  let VelocityNhalf = map2_4d (\d dold -> (1/dt)*(d-dold)) u uold
  let internalEnergy = map2_4d (*) f r
    |> flatten_3d
    |> flatten
    |> f64.sum
  let kineticEnergy = map2_4d (\v m -> v*m*v) VelocityNhalf M
    |> flatten_3d
    |> flatten
    |> f64.sum
  in (kineticEnergy / internalEnergy)

def norm [n][m][l][k] (a :[n][m][l][k]f64) :f64 =
  map_4d (\x -> x*x) a
  |> flatten_3d
  |> flatten_to (n*m*l*k)
  |> f64.sum
  |> f64.sqrt

def stoppingCriteriaClassic [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :f64 =
  let r = assembleNonlinearResidual x u
    |> map2_4d (\bb rr -> rr - bb) f
  in (norm r) / (norm f)

entry explicitQuasiStaticSolve [nx][ny][nz][nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f64) =
  let M = assembleLumpedMassMatrix x :> [nx][ny][nz][3]f64
  let initial = replicate nx (replicate ny (replicate nz [0,0,0f64]))
  let inital_steps = 2000
  let max_steps = 10000
  let (u, uold) =
  loop (u,uold) = (initial,initial) for i < inital_steps do
    let ls = getCurrentLoadScale inital_steps i
    let unew = explicitTimeStepNonlin ls x f M u uold
    in (unew, u)

  let inner_steps = 100
  let tol = 1e-5

  let innerLoop (u,uold) = 
    let unew = explicitTimeStepNonlin 1 x f M u uold
    in (unew, u)

  let (u, _, it, res) =
  loop (u,uold,it,res) = (u,uold,inital_steps,1) while (res > tol && it < max_steps) do
    let (u, uold) = (iterate (i32.i64 inner_steps) innerLoop) (u,uold)
    let res = #[trace(res)] stoppingCriteriaClassic x f u
    in (u, uold, it+inner_steps, res)
  in (u,res,i32.i64 it)

-- ==
-- entry: explicitTimeStepNonlin
-- compiled random input { 0i64 100i64 [64][64][64]f32 [65][65][65][3]f64 [65][65][65][3]f64 [65][65][65][3]f64 [65][65][65][3]f64 [65][65][65][3]f64 }
-- compiled random input { 0i64 100i64 [128][128][128]f32 [129][129][129][3]f64 [129][129][129][3]f64 [129][129][129][3]f64 [129][129][129][3]f64 [129][129][129][3]f64 }

-- ==
-- entry: explicitQuasiStaticSolve
-- compiled random input { [64][8][8]f32 [65][9][9][3]f32 }
-- compiled random input { [64][64][64]f32 [65][65][65][3]f32 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f32 }
