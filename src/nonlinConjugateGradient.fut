import "assembleNonlinear"
import "multigrid"
import "utility"
import "assembly"

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

def linesearch [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u0 :[nx][ny][nz][3]f64) (g0 :[nx][ny][nz][3]f64) (d :[nx][ny][nz][3]f64) =
  let a0 = 0
  let f0 = innerProduct d g0
  let a1 = 1
  let u1 = map2_4d (+) u0 d
  let g1 = assembleNonlinearResidual x u1
    |> map2_4d (\bb rr -> rr - bb) b
  let f1 = innerProduct d g1
  let a2 = a1 - f1*(a1-a0)/(f1-f0)
  in a2

def linesearchFine [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u0 :[nx][ny][nz][3]f64) (g0 :[nx][ny][nz][3]f64) (d :[nx][ny][nz][3]f64) =
  let a0 = 0
  let f0 = innerProduct d g0
  let a1 = 1
  let u1 = map2_4d (+) u0 d
  let g1 = assembleNonlinearResidual x u1
    |> map2_4d (\bb rr -> rr - bb) b
  let f1 = innerProduct d g1
  let a2 = a1 - f1*(a1-a0)/(f1-f0)
  let u2 = map2_4d (\uu dd -> uu+a2*dd) u0 d
  let g2 = assembleNonlinearResidual x u2
    |> map2_4d (\bb rr -> rr - bb) b
  let f2 = innerProduct d g2
  let a3 = a2 - f2*(a2-a1)/(f2-f1)
  let u3 = map2_4d (\uu dd -> uu+a3*dd) u0 d
  let g3 = assembleNonlinearResidual x u3
    |> map2_4d (\bb rr -> rr - bb) b
  let f3 = innerProduct d g3
  let a4 = a3 - f3*(a3-a2)/(f3-f2)
  in a4

def backtrackingLineSearch [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (dfx :[nx][ny][nz][3]f64) (d :[nx][ny][nz][3]f64) =
  let alpha = 1
  let c     = 1e-2
  let rho   = 0.9
  let getNextU a = map2_4d (\uu dd -> uu + a*dd) u d
  let getfx    a = (getNonlinearEnergy x (getNextU a)) - f64.sum (flatten_4d (map2_4d (*) f (getNextU a)))
  let fx    = getfx 0
  let slope = c * (innerProduct dfx d)
  let res_l = getfx alpha
  let res_r = fx+alpha*slope
  let it = 0
  let (alpha,_,_,_) = loop (alpha, it, res_l, res_r) while (res_l > res_r && it < 25) do
    let alpha = rho*alpha
    let res_l = getfx alpha
    let res_r = fx+alpha*slope
    in (alpha, it+1, res_l, res_r)
  in alpha

def linesearchZoom (getphi: f64 -> f64) (getdphi: f64 -> f64) phi0 dphi0 c1 c2 alphaLo alphaHi =
  let (_,alpha,_,_) = loop (alphaLo, alphaHi, it, done) = (alphaLo, alphaHi, 1, false) while (it < 5 && !done) do
    let alphaTrial = 0.5*alphaLo + 0.5*alphaHi
    let phi = getphi alphaTrial
    in if (phi >= (getphi alphaLo))
      then (alphaLo, alphaTrial, it+1, false)
    else
      let dphi = getdphi alphaTrial
      in if (f64.abs(dphi) <= -c2*dphi0)
        then (alphaTrial, alphaTrial, it+1, true)
      else if (dphi*(alphaHi-alphaLo) >= 0)
        then (alphaTrial, alphaLo, it+1, false)
      else (alphaTrial, alphaHi, it+1, false)
  in alpha

def wolfeLineSearch [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (f :[nx][ny][nz][3]f64) (g0 :[nx][ny][nz][3]f64) (d :[nx][ny][nz][3]f64) =
  let getU     a = map2_4d (\uu dd -> uu + a*dd) u d
  let getR     a = assembleNonlinearResidual x (getU a) |> map2_4d (\bb rr -> rr - bb) f
  let getphi   a = (getNonlinearEnergy x (getU a)) - f64.sum (flatten_4d (map2_4d (*) f (getU a)))
  let getdphi  a = innerProduct d (getR a)
  let phi0       = getphi 0
  let dphi0      = innerProduct g0 d
  let alpha_max  = 1.5
  let c1         = 1e-4
  let c2         = 0.9
  let zoom       = linesearchZoom getphi getdphi phi0 dphi0 c1 c2

  let (alpha,_,_,_) = loop (alpha, alphaOld, it, done) = (1, 0, 1, false) while (it < 5 && !done) do
    let phi  = getphi alpha
    in 
      if (phi > phi0 + c1*alpha*dphi0)
        then 
          let alpha = #[trace(zoom1)] zoom alphaOld alpha
          in (alpha, 0, it+1, true)
      else 
        let dphi = getdphi alpha
      in if (f64.abs(dphi) <= -c2*dphi0)
        then (alpha, 0, it+1, true)
      else if (dphi >= 0)
        then
          let alpha = #[trace(zoom2)] zoom alpha alphaOld
          in (alpha, 0, it+1, true)
      else
        let alphaNew = #[trace(nextTrial)] 0.5*alpha + 0.5*alpha_max
        in (alphaNew, alpha, it+1, false)
  in alpha

def cgLinearised [nelx][nely][nelz][nx][ny][nz] (tol :f64) (maxIt :i32) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u0 :[nx][ny][nz][3]f64) =
  let bnorm   = norm b
  let zero_4d = replicate_4d nx ny nz 3 0f64

  let mgData = generateMultigridDataNonlinear x u0
  let r = applyLinearisedStiffnessMatrix x u0 u0
    |> map2_4d (\bb rr -> bb - rr) b

  let  (u, _, _, _, its, _) =
    loop (uold, rold, pold, res, its, rhoold) = (u0, r, zero_4d, 1f64, 0i32, 1f64) while (res > tol && its < maxIt) do
      let z      = vcycle_l0_nonlinear mgData x u0 rold
      let rho    = innerProduct rold z
      let beta   = rho / rhoold
      let p      = map2_4d (\pp zz -> beta * pp + zz) pold z
      let q      = applyLinearisedStiffnessMatrix x u0 p
      let alpha  = rho / (innerProduct p q)
      let u      = map2_4d (\uu pp -> uu + alpha * pp) uold p
      let r      = map2_4d (\rr qq -> rr - alpha * qq) rold q
      let relres = (norm r) / bnorm
    in (u, r, p, relres, its+1, rho)
  in (u,its)

def flexiblecgLinearised [nelx][nely][nelz][nx][ny][nz] (reltol :f64) (maxIt :i32) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u_lin :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let bnorm   = norm b
  let abstol  = 1e-50*reltol
  let zero_4d = replicate_4d nx ny nz 3 0f64

  let mgData = generateMultigridDataNonlinear x u_lin

  let K = mgData.0 :> [nx][ny][nz][3][81]f64
  let applyStiffnessMatrix = applyAssembledStiffnessMatrix K

  let r0 = applyStiffnessMatrix u
    |> map2_4d (\bb rr -> bb - rr) b

  let  (u,_,_,_,_,_,_,its) =
    loop (u, r, pold, zold, rhoold, relres, absres, its) = (u, r0, zero_4d, zero_4d, 1, 1, 1, 0i32) while (relres > reltol && absres > abstol && its < maxIt) do
      let z      = vcycle_l0_nonlinear mgData x u_lin r
      let beta   = 
        if its==0 
          then 0
          else (innerProduct r (map2_4d (-) z zold)) / rhoold
      let p      = map2_4d (\pp zz -> beta * pp + zz) pold z
      let q      = applyStiffnessMatrix p
      let rho    = innerProduct r z
      let alpha  = rho / (innerProduct p q)
      let u      = map2_4d (\pp uu -> uu + alpha * pp) p u
      let r      = map2_4d (\qq rr -> rr - alpha * qq) q r
      let relres = (norm r) / bnorm
      let absres = (norm r) 
      in (u, r, p, z, rho, relres, absres, its+1)
  in (u,its)

def bicgstabLinearised [nelx][nely][nelz][nx][ny][nz] (reltol :f64) (maxIt :i32) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u_lin :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let zero_4d = replicate_4d nx ny nz 3 0f64
  let bnorm = norm b

  let mgData = generateMultigridDataNonlinear x u_lin
  let r0 = applyLinearisedStiffnessMatrix x u_lin u
    |> map2_4d (\bb rr -> bb - rr) b
  let r0hat = r0

  let applyPrec = vcycle_l0_nonlinear mgData x u_lin
  let applyMat  = applyLinearisedStiffnessMatrix x u_lin

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

  let  (u,_,_,_,_,_,_,_,its) = 
    loop (u, r, p, v, rho, omega, alpha, res, it) = (zero_4d, b, zero_4d, zero_4d, 1f64, 1f64, 1f64, 1, 0) while res > reltol && it < maxIt do
      let (u, r, p, v, rho, omega, alpha) = inner_iteration (u, r, p, v, rho, omega, alpha)
      let res = (norm r) / bnorm
      in (u, r, p, v, rho, omega, alpha, res, it+2)

  in (u,its)


def nonlinConjugateGradient [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let tol     = 1e-5
  let maxIt   = 50
  let approximateSolver = cgLinearised 1e-2 50
  let bnorm   = norm b

  let g = assembleNonlinearResidual x u
    |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual

  let  (u, _, _, _, _, relres, it, lits) =
    loop (uold, gold, gold2, yold, dold, res, its, lits) = (u, g, g, g, g, 1, 0i32, 0) while (res > tol && its < maxIt) do

      let (y,linear_it) = approximateSolver x gold uold -- preconditioned gradient

      -- update search direction by conjugacy
      let g_update = map2_4d (-) gold gold2

      -- Fletcher-Reeves, simple, fast, flaky
      -- let beta = (innerProduct g y) / (innerProduct gold yold)

      -- Polak-Ribiere-Polyak, robust due to self-correction
      let beta = if (its == 0) 
        then 0
        else (innerProduct g_update y) / (innerProduct gold2 yold)

      -- Hestensen-Steifel, robust independently on accuracy of line-search
      -- equivalent to HS if line-search is exact
      -- let beta = (innerProduct g_update y) / (innerProduct g_update dold)

      let d = map2_4d (\yy dd -> -1*yy+beta*dd) y dold

      -- perform step
      let alpha = linesearchFine x b uold gold d
      let u = map2_4d (\uu dd -> uu+alpha*dd) uold d

      -- find new search direction
      let g = assembleNonlinearResidual x u
        |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual

      -- compute stopping criteria
      let relres = (norm g) / bnorm

    in (u, g, gold, y, d, relres, its+1, lits+linear_it)
  in (u, relres, it, lits)

def nonlinConjugateGradientJacobi [nelx][nely][nelz][nx][ny][nz] (tol :f64) (maxIt :i32) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let bnorm   = norm b

  let g = assembleNonlinearResidual x u
    |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual

  let y = map2_4d (\v d -> v/d) g (getNonlinearStiffnessDiagonal x u)
  let d = map_4d (* -1) y -- search direction

  let  (u, _, _, _, relres, it) =
    loop (uold, gold, yold, dold, res, its) = (u, g, y, d, 1, 0i32) while (res > tol && its < maxIt) do

      -- perform step
      -- let alpha :f32 = 1 -- step length
      let alpha = linesearch x b uold gold dold
      let u = map2_4d (\uu dd -> uu+alpha*dd) uold dold

      -- find new search direction
      let g = assembleNonlinearResidual x u
        |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual
      let y = map2_4d (\v d -> v/d) g (getNonlinearStiffnessDiagonal x u)

      -- update search direction by conjugacy
      let g_update = map2_4d (-) g gold

      -- Fletcher-Reeves, simple, fast, flaky
      -- let beta = (innerProduct g y) / (innerProduct gold yold)

      -- Polak-Ribiere-Polyak, robust due to self-correction
      let beta = (innerProduct g_update y) / (innerProduct gold yold)

      -- Hestensen-Steifel, robust independently on accuracy of line-search
      -- equivalent to HS if line-search is exact
      -- let beta = (innerProduct g_update y) / (innerProduct g_update dold)

      let d = map2_4d (\yy dd -> -1*yy+beta*dd) y dold

      -- compute stopping criteria
      let its = (its+1)
      let relres = ((norm g) / bnorm)

    in (u, g, y, d, relres, its)
  in (u, relres, it, 0i32)

def nonlinConjugateGradientNoPre [nelx][nely][nelz][nx][ny][nz] (tol :f64) (maxIt :i32) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let bnorm   = norm b
  let g = assembleNonlinearResidual x u
    |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual
  let d = map_4d (* -1) g -- search direction
  let  (u, _, _, relres, it) =
    loop (uold, gold, dold, res, its) = (u, g, d, 1, 0i32) while (res > tol && its < maxIt) do
      let alpha = linesearch x b uold gold dold
      let u = map2_4d (\uu dd -> uu+alpha*dd) uold dold
      let g = assembleNonlinearResidual x u
        |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual
      let g_update = map2_4d (-) g gold
      let beta = (innerProduct g_update g) / (innerProduct g_update dold)
      let d = map2_4d (\yy dd -> -1*yy+beta*dd) g dold
      let its = (its+1)
      let relres = (norm g) / bnorm
    in (u, g, d, relres, its)
  in (u, relres, it, 0i32)

def nonlinNewton [nelx][nely][nelz][nx][ny][nz] (tol :f64) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let maxIt   = 50
  let bnorm   = norm b
  -- let tol_alpha = 2
  -- let tol_gamma = 0.9
  let max_tol   = 1e-1

  let g = assembleNonlinearResidual x u
        |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual

  let relres = (norm g) / bnorm

  let  (u, _, _, _, relres, it, lits) =
    loop (uold, gold, goldnorm, lintol, res, its, lits) = (u, g,  relres*bnorm, max_tol, relres, 0i32, 0) while (res > tol && its < maxIt) do

      -- find new search direction
      let (d, linear_it) = flexiblecgLinearised lintol 150 x gold uold uold -- search direction
      let d = map_4d (* -1) d

      -- perform step
      let alpha = backtrackingLineSearch x uold b gold d
      -- let alpha = #[trace(alpha)] wolfeLineSearch x uold b gold d
      let step = map_4d (*alpha) d
      let u = map2_4d (\uu ss -> uu+ss) uold step

      -- update residual
      let g = assembleNonlinearResidual x u
        |> map2_4d (\bb rr -> rr - bb) b -- gradient / residual

      -- compute stopping criteria
      let gnorm = norm g
      let relres = (gnorm / bnorm)

      -- compute linear tolerance for next iteration
      let linear_direction = applyLinearisedStiffnessMatrix x uold step
      let decrease = map3_4d (\rn ro rl -> rn - ro - rl) g gold linear_direction
      let eta = (norm decrease) / goldnorm
      let safeguard_eta = lintol ** ((1+(f64.sqrt 5))/2)

      -- let eta = tol_gamma*((gnorm/goldnorm)**tol_alpha)
      -- let safeguard_eta = tol_gamma*(lintol**tol_alpha)

      let lintol_new = 
        if safeguard_eta > 0.1
          then f64.max eta safeguard_eta
          else eta

      let lintol_new = f64.min lintol_new max_tol

    in (u, g, gnorm, lintol_new, relres, its+1, lits+linear_it)
  in (u, relres, it, lits)
