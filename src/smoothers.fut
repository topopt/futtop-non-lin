import "applyStiffnessMatrix"
import "utility"
import "assembleNonlinear"

def jacobiSmoother [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (invD :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let nsweeps = 2
  let omega   = 0.6
  let smooth u = applyCoarseStiffnessMatrix l x u
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) u b invD
  in (iterate nsweeps smooth) u
      
def jacobiSmootherSingle [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f32) (invD :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) :[nx][ny][nz][3]f32 =
  let nsweeps = 2
  let omega   = 0.6
  let smooth u = applyStiffnessMatrixSingle x u
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) u b invD
  in (iterate nsweeps smooth) u

def jacobiSmootherNonlinear [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (invD :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let nsweeps = 2
  let omega   = 0.6
  let smooth u = assembleNonlinearResidual x u
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) u b invD
  in (iterate nsweeps smooth) u

def jacobiSmootherLinearised [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (invD :[nx][ny][nz][3]f64) (u0 :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let nsweeps = 1
  let omega   = 0.6
  let smooth u = applyLinearisedStiffnessMatrix x u0 u
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) u b invD
  in (iterate nsweeps smooth) u

def jacobiSmootherLinearisedLevel1 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u0 :[][][][3]f64) (b :[nxc][nyc][nzc][3]f64) (invD :[nxc][nyc][nzc][3]f64) (u :[nxc][nyc][nzc][3]f64) :[nxc][nyc][nzc][3]f64 =
  let nsweeps = 1
  let omega   = 0.6
  let smooth v = applyLinearisedStiffnessMatrixLevel1 x u0 v
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) v b invD
  in (iterate nsweeps smooth) u

def jacobiPreconditioner [nx][ny][nz] (invD :[nx][ny][nz][3]f64) (r :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  map2_4d (*) invD r


entry jacobiSmoother_test [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f64) (invD :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64  =
  jacobiSmoother l x b invD u

-- ==
-- entry: jacobiSmoother_test
-- compiled random input { 0u8 [64][64][64]f32 [65][65][65][3]f64 [65][65][65][3]f64 [65][65][65][3]f64 }
-- compiled random input { 0u8 [128][128][128]f32 [129][129][129][3]f64 [129][129][129][3]f64 [129][129][129][3]f64 }
-- compiled random input { 0u8 [256][256][256]f32 [257][257][257][3]f64 [257][257][257][3]f64 [257][257][257][3]f64 }
