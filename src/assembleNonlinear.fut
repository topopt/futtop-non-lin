import "utility"
import "indexUtilities"
import "boundaryConditions"
import "nonlinearElement"
import "assemblyUtilities"
import "material"
import "projection"

-- #[noinline]
entry assembleNonlinearResidual [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let esize = computeElementSize nely
  in elementwiseAssembly x u (getRegularizedResidual esize) (setBCtoInput u)

-- ==
-- entry: assembleNonlinearResidual
-- compiled random input { [8][8][8]f32 [9][9][9][3]f64 }
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f64 }

-- #[noinline]
entry applyLinearisedStiffnessMatrix [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u_exp :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let esize = computeElementSize nely
  in elementwiseAssembly x u (getLinearisedStiffnessProduct esize u_exp) (setBCtoZero 0)

entry applyLinearisedStiffnessMatrixLevel1 [nelx][nely][nelz][nx][ny][nz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u_exp :[nx][ny][nz][3]f64) (u :[nxc][nyc][nzc][3]f64) :[nxc][nyc][nzc][3]f64 =
  let ufine = projectToFiner u
  let rcoarse = applyLinearisedStiffnessMatrix x u_exp (ufine :> [nx][ny][nz][3]f64)
  in projectToCoarser rcoarse :> [nxc][nyc][nzc][3]f64

-- #[noinline]
entry getNonlinearStiffnessDiagonal [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let esize = computeElementSize nely
  in elementwiseAssembly x u (getElementStiffnessDiagonal esize) (setBCtoZero 1)

-- #[noinline]
entry getNonlinearBlockDiagonal [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  let esize = computeElementSize nely
  in elementwiseAssemblyNodal x u (getElementStiffnessBlockDiagonal esize) (replicate 3 (replicate 3 0)) matsum_f64 id

-- #[noinline]
entry getNonlinearEnergy [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) = 
  let esize = computeElementSize nely
  in tabulate_3d nelx nely nelz (\i j k -> (getRegularizedEnergy esize) x u {x=i,y=j,z=k})
  |> flatten_3d
  |> f64.sum

-- ==
-- entry: applyLinearisedStiffnessMatrix
-- compiled random input { [32][32][32]f32 [33][33][33][3]f64 [33][33][33][3]f64 }
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 [65][65][65][3]f64 }

-- ==
-- entry: getNonlinearStiffnessDiagonal
-- compiled random input { [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 }

def getFineValueNonlinear [nx][ny][nz][nelx][nely][nelz] esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (dofNumber :i64) (cellIndex :index) (w :nodalWeights) :[24]f64 =
  let w_onBoundary = applyBoundaryConditionsToWeightsInverse cellIndex w ((nelx+1),(nely+1),(nelz+1)) dofNumber
  let l_onBoundary = generateLoad dofNumber w_onBoundary
  in if indexIsInside (nelx, nely, nelz) cellIndex
    then
      let w_inDomain   = applyBoundaryConditionsToWeights cellIndex w ((nelx+1),(nely+1),(nelz+1)) dofNumber
      let l_inDomain   = generateLoad dofNumber w_inDomain
      let ke = getRegularizedTangent esize x u cellIndex
      let kevalues = vecmul_f64 ke l_inDomain
      in map2 (+) kevalues l_onBoundary
    else
      l_onBoundary

def getCellContributionCellBased [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (cellIndex :index) (nodeOffset :index) :[3][24]f64 =
  let weights = getInitialWeights nodeOffset
  let valsX = getCoarseCellContribution (getFineValueNonlinear esize x u 0) l cellIndex weights
  let valsY = getCoarseCellContribution (getFineValueNonlinear esize x u 1) l cellIndex weights
  let valsZ = getCoarseCellContribution (getFineValueNonlinear esize x u 2) l cellIndex weights
  in [valsX,valsY,valsZ]

def getAllCellContributions [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (cellIndex :index) :[8][3][24]f64 =
  getCoarseCellContributionAllValues (getRegularizedTangent esize x u) (nelx,nely,nelz) l cellIndex

def getCellContribution [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (nodeIndex :index) (elementOffset :index) :[3][24]f64 =
  let cellIndex = addIndices nodeIndex elementOffset
  let nodeOffset = {x=(-elementOffset.x),y=(-elementOffset.y),z=(-elementOffset.z)}
  in getCellContributionCellBased l esize x u cellIndex nodeOffset

def getNodeAssembledRow [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (nodeIndex :index) :[3][81]f64 =
  let cellValues  = map (getCellContribution l esize x u nodeIndex) elementOffsets
    |> transpose
    |> map (flatten_to 192)
  in map (\v -> #[sequential] reduce_by_index (replicate 81 0f64) (+) 0 elementAssembledOffsets v) cellValues

-- #[noinline]
def assembleStiffnessMatrixNonlinear [nx][ny][nz][nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[][][][3][81]f64 =
  let ncell = 2**(i64.u8 l)
  let nx = (nelx/ncell)+1
  let ny = (nely/ncell)+1
  let nz = (nelz/ncell)+1
  let esize = computeElementSize nely
  in #[incremental_flattening(only_inner)]
    tabulate_3d nx ny nz (\i j k -> getNodeAssembledRow l esize x u {x=i,y=j,z=k})

entry assembleStiffnessMatrixNonlinear_new [nx][ny][nz][nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[][][][3][81]f64 =
  let ncell = 2**(i64.u8 l)
  let nelxc = (nelx/ncell)
  let nelyc = (nely/ncell)
  let nelzc = (nelz/ncell)
  let nxc = (nelx/ncell)+1
  let nyc = (nely/ncell)+1
  let nzc = (nelz/ncell)+1
  let esize = computeElementSize nely
  let cellVals = #[incremental_flattening(only_inner)] 
    tabulate_3d (nelxc+2) (nelyc+2) (nelzc+2) (\i j k -> getAllCellContributions l esize x u {x=i-1,y=j-1,z=k-1})
  in #[incremental_flattening(only_inner)] 
    tabulate_3d nxc nyc nzc (\i j k -> gatherNodeVals cellVals {x=i+1,y=j+1,z=k+1})

entry testAlt [nx][ny][nz][nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  let ncell = 2**(i64.u8 l)
  let nxc = (nelx/ncell)+1
  let nyc = (nely/ncell)+1
  let nzc = (nelz/ncell)+1
  let size = nxc*nyc*nzc*3*81
  let a = assembleStiffnessMatrixNonlinear_new l x u |> flatten_4d |> flatten_to size
   |> map (\v -> if f64.isnan v then -10000 else v)
  let b = assembleStiffnessMatrixNonlinear l x u |> flatten_4d |> flatten_to size
   |> map (\v -> if f64.isnan v then -10000 else v)
  in map2 (-) a b |> map f64.abs |> f64.sum

-- ==
-- entry: testAlt
-- nobench compiled random input { 0u8 [8][8][8]f32 [9][9][9][3]f64 } output { 0f64 }
-- nobench compiled random input { 1u8 [8][8][8]f32 [9][9][9][3]f64 } output { 0f64 }
-- nobench compiled random input { 2u8 [8][8][8]f32 [9][9][9][3]f64 } output { 0f64 }
-- nobench compiled random input { 3u8 [8][8][8]f32 [9][9][9][3]f64 } output { 0f64 }

-- ==
-- entry: assembleStiffnessMatrixNonlinear
-- compiled random input { 0u8 [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { 1u8 [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { 2u8 [32][32][32]f32 [33][33][33][3]f64 }
  -- compiled random input { 0u8 [64][64][64]f32 [65][65][65][3]f64 }
  -- compiled random input { 1u8 [64][64][64]f32 [65][65][65][3]f64 }
  -- compiled random input { 2u8 [64][64][64]f32 [65][65][65][3]f64 }

-- ==
-- entry: assembleStiffnessMatrixNonlinearAlt
-- compiled random input { 0u8 [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { 1u8 [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { 2u8 [32][32][32]f32 [33][33][33][3]f64 }
  -- compiled random input { 0u8 [64][64][64]f32 [65][65][65][3]f64 }
  -- compiled random input { 1u8 [64][64][64]f32 [65][65][65][3]f64 }
  -- compiled random input { 2u8 [64][64][64]f32 [65][65][65][3]f64 }

def getDiagonalCellContribution [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (nodeIndex :index) (elementOffset :index) :[3]f64 =
  let vals = getCellContribution l esize x u nodeIndex elementOffset
  let nodeOffset = {x=(-elementOffset.x),y=(-elementOffset.y),z=(-elementOffset.z)}
  let li = i64.i32 (getLocalNodeIndex nodeOffset)
  in [vals[0,3*li+0],vals[1,3*li+1],vals[2,3*li+2]]

def getNodeDiagonalValues [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (nodeIndex :index) :[3]f64 =
  map (getDiagonalCellContribution l esize x u nodeIndex) elementOffsets
   |> transpose
   |> map (reduce (+) 0)

def assembleDiagonalNonlinear [nx][ny][nz][nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[][][][3]f64 =
  let ncell = 2**(i64.u8 l)
  let nx = (nelx/ncell)+1
  let ny = (nely/ncell)+1
  let nz = (nelz/ncell)+1
  let esize = computeElementSize nely
  in #[incremental_flattening(only_inner)]
    tabulate_3d nx ny nz (\i j k -> getNodeDiagonalValues l esize x u {x=i,y=j,z=k})

def assembleInverseDiagonalNonlinear [nx][ny][nz][nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[][][][3]f64 =
  assembleDiagonalNonlinear l x u |> map_4d (\x -> 1/x)

-- ==
-- entry: assembleInverseDiagonalNonlinear
-- compiled random input { 0u8 [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { 1u8 [32][32][32]f32 [33][33][33][3]f64 }
-- compiled random input { 2u8 [32][32][32]f32 [33][33][33][3]f64 }
  -- compiled random input { 0u8 [64][64][64]f32 [65][65][65][3]f64 }
  -- compiled random input { 1u8 [64][64][64]f32 [65][65][65][3]f64 }
  -- compiled random input { 2u8 [64][64][64]f32 [65][65][65][3]f64 }

def getBlockDiagonalCellContribution [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (nodeIndex :index) (elementOffset :index) :[3][3]f64 =
  let vals = getCellContribution l esize x u nodeIndex elementOffset
  let nodeOffset = {x=(-elementOffset.x),y=(-elementOffset.y),z=(-elementOffset.z)}
  let li = i64.i32 (getLocalNodeIndex nodeOffset)
  in [[vals[0,3*li+0],vals[0,3*li+1],vals[0,3*li+2]],
      [vals[1,3*li+0],vals[1,3*li+1],vals[1,3*li+2]],
      [vals[2,3*li+0],vals[2,3*li+1],vals[2,3*li+2]]]

def getNodeBlockDiagonalValues [nx][ny][nz][nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (nodeIndex :index) :[3][3]f64 =
  map (getBlockDiagonalCellContribution l esize x u nodeIndex) elementOffsets
   |> transpose
   |> map transpose
   |> map_2d f64.sum

def assembleBlockDiagonalNonlinear [nx][ny][nz][nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[][][][3][3]f64 =
  let ncell = 2**(i64.u8 l)
  let nx = (nelx/ncell)+1
  let ny = (nely/ncell)+1
  let nz = (nelz/ncell)+1
  let esize = computeElementSize nely
  in #[incremental_flattening(only_inner)]
    tabulate_3d nx ny nz (\i j k -> getNodeBlockDiagonalValues l esize x u {x=i,y=j,z=k})
