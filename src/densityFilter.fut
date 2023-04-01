type index = {x: i64, y: i64, z: i64}


def forwardDensityFilter [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  let gridRadius = i64.max 0 (i64.f32 (f32.ceil (rmin-1)))
  let boxSize    = 2*gridRadius+1

  let isInsideDomain(idx :index) :bool =
    idx.x >= 0 && idx.y >= 0 && idx.z >= 0 &&
    idx.x < nelx && idx.y < nely && idx.z < nelz

  let getFilterWeight (ownIdx :index, neighIdx :index, rmin :f32) :f32 =
    f32.max 0 (rmin - f32.sqrt( f32.i64 (
      (ownIdx.x-neighIdx.x)*(ownIdx.x-neighIdx.x) +
      (ownIdx.y-neighIdx.y)*(ownIdx.y-neighIdx.y) +
      (ownIdx.z-neighIdx.z)*(ownIdx.z-neighIdx.z))))

  in tabulate_3d nelx nely nelz (\i j k->
    let ownIdx :index = {x=i,y=j,z=k}

    let weightSum =
      reduce (+) 0 (
        tabulate boxSize (\ii ->
          reduce (+) 0 (
            tabulate boxSize (\jj ->
              reduce (+) 0 (
                tabulate boxSize (\kk ->
                  let neighIdx :index = {x=i+ii-gridRadius,y=j+jj-gridRadius,z=k+kk-gridRadius} in
                  if isInsideDomain(neighIdx) then getFilterWeight(ownIdx, neighIdx, rmin) else 0 ))))))

    let scaledDensity =
      reduce (+) 0 (
        tabulate boxSize (\ii ->
          reduce (+) 0 (
            tabulate boxSize (\jj ->
              reduce (+) 0 (
                tabulate boxSize (\kk ->
                  let neighIdx :index = {x=i+ii-gridRadius,y=j+jj-gridRadius,z=k+kk-gridRadius} in
                  if isInsideDomain(neighIdx) then
                  let wgt      = getFilterWeight(ownIdx, neighIdx, rmin)
                  in #[unsafe] (wgt * x[neighIdx.x,neighIdx.y,neighIdx.z])
                  else 0 ))))))

    in scaledDensity / weightSum)




def backwardDensityFilter [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  let gridRadius = i64.max 0 (i64.f32 (f32.ceil (rmin-1)))
  let boxSize    = 2*gridRadius+1

  let isInsideDomain(idx :index) :bool =
    idx.x >= 0 && idx.y >= 0 && idx.z >= 0 &&
    idx.x < nelx && idx.y < nely && idx.z < nelz

    let getFilterWeight (ownIdx :index, neighIdx :index, rmin :f32) :f32 =
      f32.max 0 (rmin - f32.sqrt( f32.i64 (
        (ownIdx.x-neighIdx.x)*(ownIdx.x-neighIdx.x) +
        (ownIdx.y-neighIdx.y)*(ownIdx.y-neighIdx.y) +
        (ownIdx.z-neighIdx.z)*(ownIdx.z-neighIdx.z))))

  let weigthSums =
  tabulate_3d nelx nely nelz (\i j k->
    let ownIdx :index = {x=i,y=j,z=k} in
    reduce (+) 0 (
      tabulate boxSize (\ii ->
        reduce (+) 0 (
          tabulate boxSize (\jj ->
            reduce (+) 0 (
              tabulate boxSize (\kk ->
                let neighIdx :index = {x=i+ii-gridRadius,y=j+jj-gridRadius,z=k+kk-gridRadius} in
                if isInsideDomain(neighIdx)
                then getFilterWeight(ownIdx, neighIdx, rmin)
                else 0 )))))))

  in tabulate_3d nelx nely nelz (\i j k->
    let ownIdx :index = {x=i,y=j,z=k}

    in reduce (+) 0 (
      tabulate boxSize (\ii ->
        reduce (+) 0 (
          tabulate boxSize (\jj ->
            reduce (+) 0 (
              tabulate boxSize (\kk ->
                let neighIdx :index = {x=i+ii-gridRadius,y=j+jj-gridRadius,z=k+kk-gridRadius} in
                if isInsideDomain(neighIdx) then
                let wgt = getFilterWeight(ownIdx, neighIdx, rmin)
                in #[unsafe](wgt * (x[neighIdx.x,neighIdx.y,neighIdx.z] / weigthSums[neighIdx.x,neighIdx.y,neighIdx.z]))
                else 0)))))))


entry forwardDensityFilter_test [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  forwardDensityFilter rmin x

entry backwardDensityFilter_test [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  backwardDensityFilter rmin x

-- ==
-- entry: forwardDensityFilter_test
-- nobench input @../testData/densityFilter.txt output @../testData/densityFilterForwardOut.txt
-- compiled random input { 1.5f32 [64][64][64]f32 }
-- compiled random input { 2.5f32 [64][64][64]f32 }
-- compiled random input { 5.5f32 [64][64][64]f32 }
-- compiled random input { 1.5f32 [128][128][128]f32 }

-- ==
-- entry: backwardDensityFilter_test
-- nobench input @../testData/densityFilter.txt output @../testData/densityFilterBackwardOut.txt
-- compiled random input { 1.5f32 [64][64][64]f32}
-- compiled random input { 2.5f32 [64][64][64]f32}
-- compiled random input { 5.5f32 [64][64][64]f32}
-- compiled random input { 1.5f32 [128][128][128]f32}
