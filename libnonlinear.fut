import "src/nonlinConjugateGradient"
import "src/complianceSensitivityNonlinear"
import "src/optimization"
import "src/densityFilter"
import "src/utility"

-- Update design variables x given displacement field and volume fraction.
entry designIteration [nelx][nely][nelz][nx][ny][nz] (x :*[nelx][nely][nelz]f32) (xPhys :*[nelx][nely][nelz]f32) (lambda :*[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (volfrac :f32) (rmin :f32) =
  let (compliance,dc,lambda) = getComplianceSensitivityNonlinear xPhys lambda u f
  let dc = backwardDensityFilter rmin dc
  let dv = replicate nelx (replicate nely (replicate nelz 1f32))
        |> backwardDensityFilter rmin
  let vol = getVolume xPhys
  let volScaled = vol / (f32.i64 (nelx*nely*nelz))
  let g = vol - (volfrac * (f32.i64 (nelx*nely*nelz)))
  let (change,x) = updateDensignOC x dc dv g
  in (compliance,change,volScaled,x,lambda)

entry designIterationMMA [nelx][nely][nelz][nx][ny][nz] (iter :i32)  (x :[nelx][nely][nelz]f32) (xold :[nelx][nely][nelz]f32) (xold2 :[nelx][nely][nelz]f32) (xPhys :*[nelx][nely][nelz]f32) 
(u :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (lambda :*[nx][ny][nz][3]f64) (asL :*[nelx][nely][nelz]f32) (asU :*[nelx][nely][nelz]f32)
(volfrac :f32) (rmin :f32) =
  let (compliance,dc,lambda) = getComplianceSensitivityNonlinear xPhys lambda u f
  let dc = backwardDensityFilter rmin dc
  let dv = replicate nelx (replicate nely (replicate nelz 1f32))
        |> backwardDensityFilter rmin
  let vol = getVolume xPhys
  let volScaled = vol / (f32.i64 (nelx*nely*nelz))
  let g = vol - (volfrac * (f32.i64 (nelx*nely*nelz)))
  let (xnew,asL,asU) = updateDesignMMA iter x xold xold2 dc dv asL asU g
  let change = map2_3d (\xn xo -> f32.abs(xn-xo)) xnew x |> flatten_3d |> reduce (f32.max) 0
  in (compliance,change,volScaled,xnew,x,xold,asL,asU,lambda)

-- Applies the density filter with radius rmin to the structured grid variables given in x.
entry densityFilter [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  forwardDensityFilter rmin x

-- Solve the system of eqations for u K(xPhys) u = f, using uold as initial guess.
entry solveSystem [nelx][nely][nelz][nx][ny][nz] (xPhys :[nelx][nely][nelz]f32) (uold :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) =
  nonlinNewton 1e-5 xPhys f uold

entry scaleVector (s :f64) = map_4d (*s)


entry generate_solve_input (nelx :i64) (nely :i64) (nelz :i64) =
  let nx = nelx+1
  let ny = nely+1
  let nz = nelz+1
  let x = replicate_3d nelx nely nelz 1f32
  let u = replicate_4d nx ny nz 3 0f64
  let fval = -0.1
  let fvals = replicate ny [0,0,fval]
  let fidx = zip3 (replicate ny nelx) (iota ny) (replicate ny 0)
  let f = scatter_3d (copy u) fidx fvals
  let f = f with [nelx,0,0,2] = fval/2
  let f = f with [nelx,nely,0,2] = fval/2
  in (x, u, f)

-- ==
-- entry: solveSystem
-- compiled script input { generate_solve_input 8i64 8i64 8i64 }
-- compiled script input { generate_solve_input 32i64 32i64 32i64 }
-- compiled script input { generate_solve_input 64i64 64i64 64i64 }
