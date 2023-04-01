import "assemblyWeights"
import "indexUtilities"
import "boundaryConditions"
import "utility"

type localState   = [24]f64
type nodalWeights = [8]f64

def elementOffsets :[8]index =
  [{x=(-1),y=( 0),z=(-1)}, {x=( 0),y=( 0),z=(-1)},
   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)},
   {x=(-1),y=( 0),z=( 0)}, {x=( 0),y=( 0),z=( 0)},
   {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)}]

def nodeOffsets :[8]index =
  [{x=1,y=0,z=1},{x=0,y=0,z=1},
   {x=0,y=1,z=1},{x=1,y=1,z=1},
   {x=1,y=0,z=0},{x=0,y=0,z=0},
   {x=0,y=1,z=0},{x=1,y=1,z=0}]

-- generates initial weights with zeros and 1 on the local index
def getInitialWeights (nodeOffset :index) :nodalWeights =
 let li = getLocalNodeIndex nodeOffset
 in (replicate 8 0) with [li] = 1

--   map getInitialWeights nodeOffsets
def initialWeights = 
  [[0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 1f64, 0f64],
   [0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 1f64],
   [0f64, 0f64, 0f64, 0f64, 1f64, 0f64, 0f64, 0f64],
   [0f64, 0f64, 0f64, 0f64, 0f64, 1f64, 0f64, 0f64],
   [0f64, 0f64, 1f64, 0f64, 0f64, 0f64, 0f64, 0f64],
   [0f64, 0f64, 0f64, 1f64, 0f64, 0f64, 0f64, 0f64],
   [1f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64],
   [0f64, 1f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64]]

-- zeroes out indices on boundary
def applyBoundaryConditionsToWeights (elementIndex :index) (w :nodalWeights) (nx :i64, ny :i64, nz :i64) (dofNumber :i64) :nodalWeights =
 getNodeIndices elementIndex
   |> (#[sequential] map (\x -> isOnBoundary (nx,ny,nz) x dofNumber))
   |> (#[sequential] map2 (\v b -> if b then 0 else v) w)

-- zeroes out indices not on boundary
def applyBoundaryConditionsToWeightsInverse (elementIndex :index) (w :nodalWeights) (nx :i64, ny :i64, nz :i64) (dofNumber :i64) :nodalWeights =
 getNodeIndices elementIndex
   |> (#[sequential] map (\x -> isOnBoundary (nx,ny,nz) x dofNumber && (indexIsInside (nx,ny,nz) x)))
   |> (#[sequential] map2 (\v b -> if b then (v/8) else 0) w)

def generateLoad (o :i64) (w :nodalWeights) :*localState =
 #[sequential] scatter (replicate 24 0) [0+o,3+o,6+o,9+o,12+o,15+o,18+o,21+o] w

def prolongateCellIndices (cellIndex :index) =
  getNodeIndices {x=2*cellIndex.x, y=2*cellIndex.y, z=2*cellIndex.z}

-- prolongs a set of eight nodal weights and the cell index
def prolongateCellValues (w :nodalWeights) :[8]nodalWeights =
  #[sequential] map (\pMat -> vecmul_f64 pMat w) prolongationWeights

-- combines eight values to one. Is the transpose operation of prolongateCellValues
def restrictCell (vals :[8][24]f64) :[24]f64 =
  let indX   = [0,3,6, 9,12,15,18,21]
  let indY   = [1,4,7,10,13,16,19,22]
  let indZ   = [2,5,8,11,14,17,20,23]
  let indAll = [0,3,6, 9,12,15,18,21,1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23]

  let valX = #[sequential] map (\v -> map (\i -> #[unsafe](v[i])) indX) vals
    |> map2 (\pMat a -> vecmul_f64 (transpose pMat) a) prolongationWeights
    |> transpose
    |> map (reduce (+) 0)

  let valY = #[sequential] map (\v -> map (\i -> #[unsafe](v[i])) indY) vals
    |> map2 (\pMat a -> vecmul_f64 (transpose pMat) a) prolongationWeights
    |> transpose
    |> map (reduce (+) 0)

  let valZ = #[sequential] map (\v -> map (\i -> #[unsafe](v[i])) indZ) vals
    |> map2 (\pMat a -> vecmul_f64 (transpose pMat) a) prolongationWeights
    |> transpose
    |> map (reduce (+) 0)

  in scatter (replicate 24 0) indAll ([valX,valY,valZ] |> flatten_to 24)

-- computes the contribution of a coarse grid cell to a
def getCoarseCellContribution (getFineValue: index -> nodalWeights -> localState) (l :u8) (cellIndex :index) (w :nodalWeights) =
  -- this is implemented as a case-match to allow the compiler to reason about
  -- array sizes, and realize that we will always restrict to an array of 24
  -- values.
  match l
    case 0 ->
      getFineValue cellIndex w
    case 1 ->
      let fineCellWeights = prolongateCellValues w
      let fineCellIndices = prolongateCellIndices cellIndex
      let fineValues = map2 getFineValue fineCellIndices fineCellWeights
      in restrictCell fineValues
    case 2 ->
      let fineCellWeights = prolongateCellValues w
        |> map prolongateCellValues
      let fineCellIndices = prolongateCellIndices cellIndex
        |> map prolongateCellIndices
      let fineValues = map2_2d getFineValue fineCellIndices fineCellWeights
      in fineValues
        |> map restrictCell
        |>     restrictCell
    case 3 ->
      let fineCellWeights = prolongateCellValues w
        |> map    prolongateCellValues
        |> map_2d prolongateCellValues
      let fineCellIndices = prolongateCellIndices cellIndex
        |> map    prolongateCellIndices
        |> map_2d prolongateCellIndices
      let fineValues = map2_3d getFineValue fineCellIndices fineCellWeights
      in fineValues
        |> map_2d restrictCell
        |> map    restrictCell
        |>        restrictCell
    case _ ->
      replicate 24 0 -- not implemented

def applyFineMatrix (nelx,nely,nelz) (cellIndex :index) (ke :[24][24]f64) (w :nodalWeights) (dofNumber :i64) :[24]f64 =
  let w_onBoundary = applyBoundaryConditionsToWeightsInverse cellIndex w ((nelx+1),(nely+1),(nelz+1)) dofNumber
  let l_onBoundary = generateLoad dofNumber w_onBoundary
  in if indexIsInside (nelx, nely, nelz) cellIndex
    then
      let w_inDomain   = applyBoundaryConditionsToWeights cellIndex w ((nelx+1),(nely+1),(nelz+1)) dofNumber
      let l_inDomain   = generateLoad dofNumber w_inDomain
      let kevalues = vecmul_f64 ke l_inDomain
      in map2 (+) kevalues l_onBoundary
    else
      l_onBoundary

def getFineCellContribution (getFineMatrix: index -> [24][24]f64) (nelx,nely,nelz) (cellIndex :index) (weights :[8]nodalWeights) :[8][3][24]f64 =
  let ke = getFineMatrix cellIndex
  in map (\w -> map (applyFineMatrix (nelx,nely,nelz) cellIndex ke w) (iota 3)) weights

-- computes the contribution of a coarse grid cell to a
def getCoarseCellContributionAllValues (getFineMatrix: index -> [24][24]f64) (nelx,nely,nelz) (l :u8) (cellIndex :index) =
  -- this is implemented as a case-match to allow the compiler to reason about
  -- array sizes, and realize that we will always restrict to an array of 24
  -- values.
  let w = (copy initialWeights)
  let getFineMethod = getFineCellContribution getFineMatrix (nelx,nely,nelz)
  in match l
     case 0 ->
       getFineMethod cellIndex w
     case 1 ->
       let fineCellWeights = map prolongateCellValues w 
          |> transpose
       let fineCellIndices = prolongateCellIndices cellIndex
       let fineValues = map2 getFineMethod fineCellIndices fineCellWeights |> transpose |> map transpose
       in map_2d restrictCell fineValues
     case 2 ->
       let fineCellWeights = map (prolongateCellValues >-> map prolongateCellValues) w 
          |> transpose |> map transpose
       let fineCellIndices = (prolongateCellIndices >-> map prolongateCellIndices) cellIndex
       let fineValues = map2_2d getFineMethod fineCellIndices fineCellWeights 
          |> map transpose |> transpose |> map_2d transpose |> map transpose
       in map_2d (map restrictCell >-> restrictCell) fineValues
    case 3 ->
      let fineCellWeights = map (prolongateCellValues >-> map prolongateCellValues >-> map_2d prolongateCellValues) w
        |> transpose |> map transpose |> map_2d transpose
      let fineCellIndices = (prolongateCellIndices >-> map prolongateCellIndices >-> map_2d prolongateCellIndices) cellIndex
      let fineValues = map2_3d getFineMethod fineCellIndices fineCellWeights 
        |> map_2d transpose |> map transpose |> transpose |> map_3d transpose |> map_2d transpose |> map transpose
      in map_2d (map_2d restrictCell >-> map restrictCell >-> restrictCell) fineValues
    case _ ->
      replicate 8 (replicate 3 (replicate 24 0)) -- not implemented


-- takes an offset with values [-1..1], and maps each compination to the three local indices
-- let offsetToAssembledOffset (idx :index) =
--   let ind = (idx.x+1)*9 + (idx.y+1)*3 + (idx.z+1)
--   in [3*ind+0,3*ind+1,3*ind+2]
--
-- let elementAssembledOffsets :[192]i64 = elementOffsets
--                                         |> map getNodeIndices
--                                         |> flatten_to 64
--                                         |> map offsetToAssembledOffset
--                                         |> flatten_to 192

def elementAssembledOffsets :[192]i64 =
  [18i64, 19i64, 20i64, 45i64, 46i64, 47i64, 36i64, 37i64, 38i64, 9i64, 10i64,
  11i64, 21i64, 22i64, 23i64, 48i64, 49i64, 50i64, 39i64, 40i64, 41i64, 12i64,
  13i64, 14i64, 45i64, 46i64, 47i64, 72i64, 73i64, 74i64, 63i64, 64i64, 65i64,
  36i64, 37i64, 38i64, 48i64, 49i64, 50i64, 75i64, 76i64, 77i64, 66i64, 67i64,
  68i64, 39i64, 40i64, 41i64, 36i64, 37i64, 38i64, 63i64, 64i64, 65i64, 54i64,
  55i64, 56i64, 27i64, 28i64, 29i64, 39i64, 40i64, 41i64, 66i64, 67i64, 68i64,
  57i64, 58i64, 59i64, 30i64, 31i64, 32i64, 9i64, 10i64, 11i64, 36i64, 37i64,
  38i64, 27i64, 28i64, 29i64, 0i64, 1i64, 2i64, 12i64, 13i64, 14i64, 39i64,
  40i64, 41i64, 30i64, 31i64, 32i64, 3i64, 4i64, 5i64, 21i64, 22i64, 23i64,
  48i64, 49i64, 50i64, 39i64, 40i64, 41i64, 12i64, 13i64, 14i64, 24i64, 25i64,
  26i64, 51i64, 52i64, 53i64, 42i64, 43i64, 44i64, 15i64, 16i64, 17i64, 48i64,
  49i64, 50i64, 75i64, 76i64, 77i64, 66i64, 67i64, 68i64, 39i64, 40i64, 41i64,
  51i64, 52i64, 53i64, 78i64, 79i64, 80i64, 69i64, 70i64, 71i64, 42i64, 43i64,
  44i64, 39i64, 40i64, 41i64, 66i64, 67i64, 68i64, 57i64, 58i64, 59i64, 30i64,
  31i64, 32i64, 42i64, 43i64, 44i64, 69i64, 70i64, 71i64, 60i64, 61i64, 62i64,
  33i64, 34i64, 35i64, 12i64, 13i64, 14i64, 39i64, 40i64, 41i64, 30i64, 31i64,
  32i64, 3i64, 4i64, 5i64, 15i64, 16i64, 17i64, 42i64, 43i64, 44i64, 33i64,
  34i64, 35i64, 6i64, 7i64, 8i64]

def gatherNodeVals (vals :[][][][8][3][24]f64) (nodeIndex :index) :[3][81]f64 =
  let cellIdx = map (addIndices nodeIndex) elementOffsets
  let cellValues = tabulate 8 (\i -> 
    let idx = #[unsafe] cellIdx[i]
    in #[unsafe] vals[idx.x,idx.y,idx.z,i] ) 
    |> transpose
    |> map (flatten_to 192)
  in map (\v -> reduce_by_index (replicate 81 0f64) (+) 0 elementAssembledOffsets v) cellValues

def getFlatElementIndices (nodeSize :index) (elementIndex :index) :[24]i64 =
  let ny = nodeSize.y
  let nz = nodeSize.z
  in getNodeIndices elementIndex
    |> map (\n -> 3*nz*ny*n.x+3*nz*n.y+3*n.z)
    |> map (\ni -> map (\d -> ni+d) (iota 3))
    |> flatten_to 24

def elementwiseAssembly [nelx][nely][nelz][nx][ny][nz] 
  (x :[nelx][nely][nelz]f32) 
  (u :[nx][ny][nz][3]f64) 
  (elementOperation : [nelx][nely][nelz]f32 -> [nx][ny][nz][3]f64 -> index -> [24]f64) 
  (postOperation : [nx][ny][nz][3]f64 -> [nx][ny][nz][3]f64) = 
  let fvals = tabulate_3d nelx nely nelz (\i j k -> elementOperation x u {x=i,y=j,z=k})
    |> flatten_3d
    |> flatten_to (24*nelx*nely*nelz)
  let domainSize = {x=nx,y=ny,z=nz}
  let ivals = tabulate_3d nelx nely nelz (\i j k -> getFlatElementIndices domainSize {x=i,y=j,z=k})
    |> flatten_3d
    |> flatten_to (24*nelx*nely*nelz)
  let f = reduce_by_index (replicate (3*nx*ny*nz) 0) (+) 0 ivals fvals
    |> unflatten (nx*ny*nz) 3
    |> unflatten_3d nx ny nz
  in postOperation f

def getNodeElementIndices (nodeSize :index) (elementIndex :index) :[8]i64 =
  let ny = nodeSize.y
  let nz = nodeSize.z
  in getNodeIndices elementIndex
    |> map (\n -> nz*ny*n.x+nz*n.y+n.z)

def elementwiseAssemblyNodal [nelx][nely][nelz][nx][ny][nz] 'a
  (x :[nelx][nely][nelz]f32) 
  (u :[nx][ny][nz][3]f64) 
  (elementOperation : [nelx][nely][nelz]f32 -> [nx][ny][nz][3]f64 -> index -> [8]a)
  (zeroElement : a) 
  (sumOperation : a -> a -> a) 
  (postOperation : [nx][ny][nz]a -> [nx][ny][nz]a) = 
  let fvals = tabulate_3d nelx nely nelz (\i j k -> elementOperation x u {x=i,y=j,z=k})
    |> flatten_3d
    |> flatten_to (8*nelx*nely*nelz)
  let domainSize = {x=nx,y=ny,z=nz}
  let ivals = tabulate_3d nelx nely nelz (\i j k -> getNodeElementIndices domainSize {x=i,y=j,z=k})
    |> flatten_3d
    |> flatten_to (8*nelx*nely*nelz)
  let f = reduce_by_index (replicate (nx*ny*nz) zeroElement) sumOperation zeroElement ivals fvals
    |> unflatten_3d nx ny nz
  in postOperation f