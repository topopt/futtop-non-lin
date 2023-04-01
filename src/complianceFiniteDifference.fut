import "complianceSensitivityNonlinear"
import "nonlinConjugateGradient"
import "utility"

def finiteDifferenceSingleElement [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) pert (ix :i64,iy :i64,iz :i64) =
    let zeros = replicate_4d nx ny nz 3 0f64
    let xpert1 = (copy x) with [ix,iy,iz] = x[ix,iy,iz] + pert
    let xpert2 = (copy x) with [ix,iy,iz] = x[ix,iy,iz] - pert
    let (u1,_,_,_) = nonlinNewton 1e-12 xpert1 b zeros 
    let (u2,_,_,_) = nonlinNewton 1e-12 xpert2 b zeros 
    let (c1,_) = getComplianceSensitivityNonlinear x u1 b
    let (c2,_) = getComplianceSensitivityNonlinear x u2 b
    let c1 = f64.f32 c1
    let c2 = f64.f32 c2
    let grad = (c1-c2) / (2* (f64.f32 pert) )
    in grad

def fdTest_single [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) idx pert = 
    let zeros = replicate_4d nx ny nz 3 0f64
    let (u,_,_,_) = nonlinNewton 1e-12 x b zeros
    let (_,dc) = getComplianceSensitivityNonlinear x u b
    let dc_fd = finiteDifferenceSingleElement x b pert idx
    let dc_e = f64.f32 dc[idx.0,idx.1,idx.2]
    in (dc_e, dc_fd, (dc_e - dc_fd) / dc_fd)

def fdTest_all [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) pert = 
    let zeros = replicate_4d nx ny nz 3 0f64
    let (u,_,_,_) = nonlinNewton 1e-12 x b zeros
    let (_,dc) = getComplianceSensitivityNonlinear x u b
    in loop maxerr = 0 for i < nelx do 
       loop maxerr for j < nely do 
       loop maxerr for k < nelz do 
        let idx = (i,j,k)
        let dc_fd = finiteDifferenceSingleElement x b pert idx
        let dc_e = f64.f32 dc[idx.0,idx.1,idx.2]
        let err = (dc_e - dc_fd) / dc_fd
        in f64.max maxerr (f64.abs(err))

entry generate_input (nelx :i64) (nely :i64) (nelz :i64) =
  let nx = nelx+1
  let ny = nely+1
  let nz = nelz+1
  let x = replicate nelx (replicate nely (replicate nelz 0.5f32))
  let f = replicate nx (replicate ny (replicate nz [0f64,0,0]))
  let fval = -0.01
  let fvals = replicate ny [0,0,fval]
  let fidx = zip3 (replicate ny nelx) (iota ny) (replicate ny 0)
  let f = scatter_3d f fidx fvals
  let f = f with [nelx,0,0,2] = fval/2
  let f = f with [nelx,nely,0,2] = fval/2
  in (x, f)

entry main (nelx :i64) (nely :i64) (nelz :i64) (pert :f32) = 
    let (x,f) = generate_input nelx nely nelz
    in fdTest_all x f pert

-- entry main (nelx :i64) (nely :i64) (nelz :i64) (tx :i64) (ty :i64) (tz :i64) (pert :f32) = 
--     let (x,f) = generate_input nelx nely nelz
--     in fdTest_single x f (tx,ty,tz) pert
