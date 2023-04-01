import "indexUtilities"
import "keConstants"
import "boundaryConditions"
import "utility"
import "assemblyUtilities"
import "material"
import "applyStiffnessMatrix"

-- takes a elementindex-weght values pair, and generates the resulting row
-- contribution of the matrix for a given dofNumber.
def getFineValue [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) esize (dofNumber :i64) (cellIndex :index) (w :nodalWeights) :[24]f64 =
  let elementScale = (getElementYoungsModule x cellIndex) * esize
  let w_onBoundary = applyBoundaryConditionsToWeightsInverse cellIndex w ((nelx+1),(nely+1),(nelz+1)) dofNumber
  let w_inDomain   = applyBoundaryConditionsToWeights        cellIndex w ((nelx+1),(nely+1),(nelz+1)) dofNumber
  let l_onBoundary = generateLoad dofNumber w_onBoundary
  let l_inDomain   = generateLoad dofNumber w_inDomain
  let kevalues = keprod elementScale l_inDomain
  in #[sequential] map2 (+) kevalues l_onBoundary

def getDiagonalCellContribution [nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (nodeIndex :index) (elementOffset :index) :[3]f64 =
  let cellIndex = addIndices nodeIndex elementOffset
  let nodeOffset = {x=(-elementOffset.x),y=(-elementOffset.y),z=(-elementOffset.z)}
  let li = i64.i32 (getLocalNodeIndex nodeOffset)
  let weights = getInitialWeights nodeOffset
  let valsX = getCoarseCellContribution (getFineValue x esize 0) l cellIndex weights
  let valsY = getCoarseCellContribution (getFineValue x esize 1) l cellIndex weights
  let valsZ = getCoarseCellContribution (getFineValue x esize 2) l cellIndex weights
  in [valsX[3*li+0],valsY[3*li+1],valsZ[3*li+2]]

def getNodeDiagonalValues [nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (nodeIndex :index) :[3]f64 =
   map (getDiagonalCellContribution l esize x nodeIndex) elementOffsets
   |> transpose
   |> map (reduce (+) 0)

def assembleDiagonal [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3]f64 =
  let ncell = 2**(i64.u8 l)
  let nx = (nelx/ncell)+1
  let ny = (nely/ncell)+1
  let nz = (nelz/ncell)+1
  let esize = computeElementSize nely
  in #[incremental_flattening(only_inner)]
    tabulate_3d nx ny nz (\i j k -> getNodeDiagonalValues l esize x {x=i,y=j,z=k})

def assembleInverseDiagonal [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3]f64 =
  assembleDiagonal l x |> map_4d (\x -> 1/x)

def extractDiagonal [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) :[nx][ny][nz][3]f64 =
  map_3d (\nodevals ->
    let x = nodevals[0,39]
    let y = nodevals[1,40]
    let z = nodevals[2,41]
    in [x,y,z]
    ) mat

def extractInverseDiagonal [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) :[nx][ny][nz][3]f64 =
  extractDiagonal mat |> map_4d (\x -> 1/x)

def getFineMatrix [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (cellIndex :index) :[24][24]f64 =
  let young = getElementYoungsModule x cellIndex
  in map_2d (*young) (copy keconst)

def getAllCellContributions [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (cellIndex :index) :[8][3][24]f64 =
  getCoarseCellContributionAllValues (getFineMatrix x) (nelx,nely,nelz) l cellIndex

def getCellContributionCellBased [nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (cellIndex :index) (nodeOffset :index) :[3][24]f64 =
  let weights = getInitialWeights nodeOffset
  let valsX = getCoarseCellContribution (getFineValue x esize 0) l cellIndex weights
  let valsY = getCoarseCellContribution (getFineValue x esize 1) l cellIndex weights
  let valsZ = getCoarseCellContribution (getFineValue x esize 2) l cellIndex weights
  in [valsX,valsY,valsZ]

def getCellContribution [nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (nodeIndex :index) (elementOffset :index) :[3][24]f64 =
  let cellIndex = addIndices nodeIndex elementOffset
  let nodeOffset = {x=(-elementOffset.x),y=(-elementOffset.y),z=(-elementOffset.z)}
  in getCellContributionCellBased l esize x cellIndex nodeOffset

def getNodeAssembledRow [nelx][nely][nelz] (l :u8) esize (x :[nelx][nely][nelz]f32) (nodeIndex :index) :[3][81]f64 =
  let cellValues  = map (getCellContribution l esize x nodeIndex) elementOffsets
    |> transpose
    |> map (flatten_to 192)
  in map (\v -> #[sequential] reduce_by_index (replicate 81 0f64) (+) 0 elementAssembledOffsets v) cellValues

entry assembleStiffnessMatrix [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3][81]f64 =
  let ncell = 2**(i64.u8 l)
  let nx = (nelx/ncell)+1
  let ny = (nely/ncell)+1
  let nz = (nelz/ncell)+1
  let esize = computeElementSize nely
  in #[incremental_flattening(only_inner)]
    tabulate_3d nx ny nz (\i j k -> getNodeAssembledRow l esize x {x=i,y=j,z=k})


entry assembleStiffnessMatrixAlt [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3][81]f64 =
  let ncell = 2**(i64.u8 l)
  let nelxc = (nelx/ncell)
  let nelyc = (nely/ncell)
  let nelzc = (nelz/ncell)
  let nxc = (nelx/ncell)+1
  let nyc = (nely/ncell)+1
  let nzc = (nelz/ncell)+1
  
  let cellVals = #[incremental_flattening(only_inner)] 
    tabulate_3d (nelxc+2) (nelyc+2) (nelzc+2) (\i j k -> getAllCellContributions l x {x=i-1,y=j-1,z=k-1})
  in #[incremental_flattening(only_inner)] 
    tabulate_3d nxc nyc nzc (\i j k -> gatherNodeVals cellVals {x=i+1,y=j+1,z=k+1})



    

-- let uoffsets: ([27]i64,[27]i64,[27]i64) =
--   let nodeOffsetsX = [replicate 9 (-1i64), replicate 9 0i64, replicate 9 1i64]
--     |> flatten_to 27
--   let nodeOffsetsY = replicate 3 ([replicate 3 (-1i64), replicate 3 0i64, replicate 3 1i64])
--     |> flatten
--     |> flatten_to 27
--   let nodeOffsetsZ = replicate 3 (replicate 3 ([-1,0,1]))
--     |> flatten
--     |> flatten_to 27
--   in (nodeOffsetsX,nodeOffsetsY,nodeOffsetsZ)

def uoffsets: ([27]i64,[27]i64,[27]i64) =
  ([-1i64, -1i64, -1i64, -1i64, -1i64, -1i64, -1i64, -1i64, -1i64, 0i64, 0i64,
  0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 1i64, 1i64, 1i64, 1i64, 1i64, 1i64,
  1i64, 1i64, 1i64], [-1i64, -1i64, -1i64, 0i64, 0i64, 0i64, 1i64, 1i64, 1i64,
  -1i64, -1i64, -1i64, 0i64, 0i64, 0i64, 1i64, 1i64, 1i64, -1i64, -1i64, -1i64,
  0i64, 0i64, 0i64, 1i64, 1i64, 1i64], [-1i64, 0i64, 1i64, -1i64, 0i64, 1i64,
  -1i64, 0i64, 1i64, -1i64, 0i64, 1i64, -1i64, 0i64, 1i64, -1i64, 0i64, 1i64,
  -1i64, 0i64, 1i64, -1i64, 0i64, 1i64, -1i64, 0i64, 1i64])

def getLocalU [nx][ny][nz] (u :[nx][ny][nz][3]f64) (nodeIndex :index) :*[81]f64 =
  let (nodeOffsetsX,nodeOffsetsY,nodeOffsetsZ) = uoffsets
  let ni = nodeIndex
  in (map3 (\i j k ->
    if (indexIsInside (nx,ny,nz) {x=ni.x+i,y=ni.y+j,z=ni.z+k}) then
      #[unsafe] (u[ni.x+i,ni.y+j,ni.z+k])
    else
      [0,0,0]) nodeOffsetsX nodeOffsetsY nodeOffsetsZ)
    |> flatten_to 81

def applyAssembledStiffnessMatrix [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) (u :[nx][ny][nz][3]f64) =
  tabulate_3d nx ny nz (\i j k ->
    let uloc = getLocalU u {x=i,y=j,z=k}
    in map (\row -> (reduce (+) 0 (map2 (\a b -> a*b) uloc row))) (mat[i,j,k])
    )




entry assembleInverseDiagonal_test [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3]f64 =
  assembleInverseDiagonal l x

entry applyAssembledStiffnessMatrix_test [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) (u :[nx][ny][nz][3]f64) =
  applyAssembledStiffnessMatrix mat u

entry assembleStiffnessMatrix_test [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3][81]f64 =
  assembleStiffnessMatrix l x

entry assembleStiffnessMatrix_fixed1 [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :[][][][3][81]f64 =
  assembleStiffnessMatrix 1 x

entry assembleStiffnessMatrix_fixed2 [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :[][][][3][81]f64 =
  assembleStiffnessMatrix 2 x

def testAssembledStiffnessMatrix [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  let mat = (assembleStiffnessMatrix l x) :> [nx][ny][nz][3][81]f64
  in applyAssembledStiffnessMatrix mat u

def testAssembledStiffnessMatrixAlt [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  let mat = (assembleStiffnessMatrixAlt l x) :> [nx][ny][nz][3][81]f64
  in applyAssembledStiffnessMatrix mat u

entry testAssembledStiffnessMatrix_fixed1 [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  testAssembledStiffnessMatrix 1 x u

entry testAssembledStiffnessMatrix_fixed2 [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  testAssembledStiffnessMatrix 2 x u

entry testAssembledStiffnessMatrixAlt_fixed1 [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  testAssembledStiffnessMatrixAlt 1 x u

entry testAssembledStiffnessMatrixAlt_fixed2 [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) =
  testAssembledStiffnessMatrixAlt 2 x u

  -- ==
  -- entry: assembleInverseDiagonal_test
  -- nobench input @../testData/assembleInverseDiagonal1.txt output @../testData/assembleInverseDiagonal1Out.txt
  -- nobench input @../testData/assembleInverseDiagonal2.txt output @../testData/assembleInverseDiagonal2Out.txt
  -- nobench input @../testData/assembleInverseDiagonal3.txt output @../testData/assembleInverseDiagonal3Out.txt
  -- compiled random input { 0u8 [128][128][128]f32 }
  -- compiled random input { 1u8 [128][128][128]f32 }
  -- compiled random input { 2u8 [128][128][128]f32 }

  -- ==
  -- entry: applyAssembledStiffnessMatrix_test
  -- compiled random input { [129][129][129][3][81]f64 [129][129][129][3]f64 }
  -- compiled random input { [65][65][65][3][81]f64 [65][65][65][3]f64 }
  -- compiled random input { [33][33][33][3][81]f64 [33][33][33][3]f64 }
  -- compiled random input { [17][17][17][3][81]f64 [17][17][17][3]f64 }

-- ==
-- entry: assembleStiffnessMatrix
-- compiled random input { 1u8 [128][128][128]f32 }
-- compiled random input { 2u8 [128][128][128]f32 }

-- ==
-- entry: assembleStiffnessMatrixAlt
-- compiled random input { 1u8 [128][128][128]f32 }
-- compiled random input { 2u8 [128][128][128]f32 }

-- ==
-- entry: assembleStiffnessMatrix_fixed1
-- compiled random input { [128][128][128]f32 }

  -- ==
  -- entry: assembleStiffnessMatrix_fixed2
  -- compiled random input { [128][128][128]f32 }

-- ==
-- entry: testAssembledStiffnessMatrix_fixed1
-- nobench input @../testData/applyStateOperatorFold3.txt output @../testData/applyStateOperatorFold3Out.txt
 -- nobench input @../testData/applyStateOperatorFold2.txt output @../testData/applyStateOperatorFold2Out.txt

-- ==
-- entry: testAssembledStiffnessMatrix_fixed2
-- nobench input @../testData/applyStateOperatorFold4.txt output @../testData/applyStateOperatorFold4Out.txt

-- ==
-- entry: testAssembledStiffnessMatrixAlt_fixed1
-- nobench input @../testData/applyStateOperatorFold3.txt output @../testData/applyStateOperatorFold3Out.txt
 -- nobench input @../testData/applyStateOperatorFold2.txt output @../testData/applyStateOperatorFold2Out.txt

-- ==
-- entry: testAssembledStiffnessMatrixAlt_fixed2
-- nobench input @../testData/applyStateOperatorFold4.txt output @../testData/applyStateOperatorFold4Out.txt
