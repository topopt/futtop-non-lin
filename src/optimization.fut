import "solvers"
import "multigrid"
import "densityFilter"
import "complianceSensitivity"
import "utility"

def getVolume [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :f32 =
  x |> flatten_to (nelx*nely)
    |> flatten_to (nelx*nely*nelz)
    |> reduce (+) 0

def getNewDensity [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (dc :[nelx][nely][nelz]f32) (dv :[nelx][nely][nelz]f32) (lambda :f32) =
  let move = 0.1f32
  let setMoveLimit v vn = f32.max 0 (f32.max (v-move) (f32.min 1 (f32.min (v+move) vn)))
  in map3_3d (\xx dc dv -> setMoveLimit xx (xx*(((-dc)/(dv*lambda)) ** 0.2))) x dc dv

def bisect (f :f32 -> f32) l1 l2 = 
  let (l1,l2) = loop (l1,l2) while (((l2 - l1) / (l1 + l2)) > 1e-6) do
    let lmid = 0.5 * (l2 + l1)
    let v = f lmid
    in if (v > 0) then
      (lmid,l2)
    else
      (l1,lmid)
  in 0.5 * (l2 + l1)

def updateDensignOC [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (dc :[nelx][nely][nelz]f32) (dv :[nelx][nely][nelz]f32) (g :f32) =
  let computeConstraint lambda = 
    let xtest = getNewDensity x dc dv lambda
    in map3_3d (\xn xo dv -> (xn-xo)*dv) xtest x dv |> flatten_3d |> reduce (+) 0 |> (+g)

  let lmid = bisect computeConstraint 0f32 1e9f32
  let xnew = getNewDensity x dc dv lmid
  let change = map2_3d (\xn xo -> f32.abs(xn-xo)) xnew x |> flatten_3d |> reduce (f32.max) 0
  in (change,xnew)

def innerProduct [n][m][l] (a :[n][m][l]f32) (b :[n][m][l]f32) :f32 =
  map2_3d (*) a b
  |> map_2d f32.sum
  |> map    f32.sum
  |>        f32.sum

def updateDesignMMA [nelx][nely][nelz] (iter :i32) 
  (x :[nelx][nely][nelz]f32) (xold :[nelx][nely][nelz]f32) (xold2 :[nelx][nely][nelz]f32)
  (dc :[nelx][nely][nelz]f32) (dv :[nelx][nely][nelz]f32) 
  (asL :[nelx][nely][nelz]f32) (asU :[nelx][nely][nelz]f32) (g :f32) =
  let moveLim = 0.1
  let asyIncr = 1.2 -- increase rate for asymptotes when things are going well
  let asyDecr = 0.7 -- decrease rate for asymptotes when things are not going well
  let beta    = 0.5 -- asymptote initialization, default 0.5, lower for tighter initial asymptotes

  let xL = map_3d (\v -> f32.max 0f32 (v-moveLim)) x -- lower bound of x
  let xU = map_3d (\v -> f32.min 1f32 (v+moveLim)) x -- upper bound of x
  let (asL,asU) = 
    if (iter < 3) then
      let update = map2_3d (\u l -> (u-l)/(beta+1)) xU xL
      let L = map2_3d (\xx up-> xx-0.5*up) x update
      let U = map2_3d (\xx up-> xx+0.5*up) x update
      in (L,U)
    else
      let tmp = map3_3d (\xx xo xo2 -> (xx-xo)*(xo-xo2)) x xold xold2
      let gm  = map_3d (\v -> if v > 0 then asyIncr else if v < 0 then asyDecr else 1) tmp 
      let L   = map4_3d (\xx xo gg low -> xx + gg * (low-xo)) x xold gm asL
      let U   = map4_3d (\xx xo gg upp -> xx + gg * (upp-xo)) x xold gm asU
      in (L,U)

  -- adaptive bounds
  let xL = map3_3d (\xx xl low -> f32.max xl (0.9*low+0.1*xx)) x xL asL
  let xU = map3_3d (\xx xu upp -> f32.min xu (0.9*upp+0.1*xx)) x xU asU

  -- split fields to positive and negative parts
  let p0_0 = map_3d (\v -> if v > 0 then v else 0) dc
  let q0_0 = map_3d (\v -> if v < 0 then v else 0) dc
  let p1_0 = map_3d (\v -> if v > 0 then v else 0) dv
  let q1_0 = map_3d (\v -> if v < 0 then v else 0) dv

  let p0 = map3_3d (\v xx upp -> v*(upp-xx)*(upp-xx)) p0_0 x asU
  let p1 = map3_3d (\v xx upp -> v*(upp-xx)*(upp-xx)) p1_0 x asU
  let q0 = map3_3d (\v xx low -> -v*(xx-low)*(xx-low)) q0_0 x asL
  let q1 = map3_3d (\v xx low -> -v*(xx-low)*(xx-low)) q1_0 x asL

  -- define methods for computing primal and dual map
  let primalProj lm = 
    let pval = map2_3d (\p0 p1 -> f32.sqrt (p0+lm*p1)) p0 p1
    let qval = map2_3d (\q0 q1 -> f32.sqrt (q0+lm*q1)) q0 q1
    let xupdate = map4_3d (\pv qv low upp -> (pv*low + qv*upp)/(pv+qv)) pval qval asL asU
    in map3_3d (\xx xl xu -> f32.min xu (f32.max xl xx)) xupdate xL xU

  let dualProj lm =
    let primal = primalProj lm
    let p0con = innerProduct (map2_3d (-) asU x) p1_0
    let q0con = innerProduct (map2_3d (-) asL x) q1_0
    let p1con = map3_3d (\p1 upp prim -> p1 / (f32.max 1e-12 (upp-prim))) p1 asU primal |> map_2d f32.sum |> map f32.sum |> f32.sum
    let q1con = map3_3d (\q1 low prim -> q1 / (f32.max 1e-12 (prim-low))) q1 asL primal |> map_2d f32.sum |> map f32.sum |> f32.sum
    in g - p0con - q0con + p1con + q1con

  let lmUpp = 1e6
  let dual0 = dualProj 0
  let dualUpp = dualProj lmUpp

  let lambda = 
    if dual0 * dualUpp < 0
    then bisect dualProj 0 lmUpp
    else if dual0 < 0
    then 0
    else lmUpp -- if dualProj lmUpp > 0

  let xnew = (primalProj lambda)
  in (xnew, asL, asU)




def performDesignIteration [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (xPhys :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (volfrac :f32) (rmin :f32) =
  let (compliance,dc) = getComplianceSensitivity xPhys u
  let dc = backwardDensityFilter rmin dc
  let dv = replicate nelx (replicate nely (replicate nelz 1f32))
        |> backwardDensityFilter rmin
  let vol = getVolume xPhys
  let volScaled = vol / (f32.i64 (nelx*nely*nelz))
  let g = vol - (volfrac * (f32.i64 (nelx*nely*nelz)))
  let (change,x) = updateDensignOC x dc dv g
  in (compliance,change,volScaled,x)


