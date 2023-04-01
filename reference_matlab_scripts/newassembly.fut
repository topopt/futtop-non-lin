#import



def coarseToFineIndex (l :u8) (nodeIndex :index) :index =
  let l = i64.u8 l
  let ncell = 2**l
  in {x=ncell*nodeIndex.x,y=ncell*nodeIndex.y,z=ncell*nodeIndex.z}  

def cellOffsetToSubcells (l :u8) (cellOffset :index) :[][][]index =
  let l = i64.u8 l
  let ncell = 2**l
  let multIndex i v = {x=v*i.x,y=v*i.y,z=v*i.z}
  let inds = replicate_3d ncell ncell ncell (multIndex cellOffset ncell)
  let offsets = map (\x-> map (\y-> map (\z-> {x=x,y=y,z=z}) (iota ncell)) (iota ncell)) (iota ncell)
  in map2_3d addIndices inds offsets

-- alternate approach using preintegration.
def getElementDiagonalContribution [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (nodeIndex :index) (cellOffset :index) (keslices :[8][3][24]f64) (elementOffset :index) :[3]f64 =
  let elementIndex = addIndices (coarseToFineIndex l nodeIndex) elementOffset
  let E = getElementYoungsModule x elementIndex
  let recievingNodeOffset :index = {x=(-cellOffset.x), y=(-cellOffset.y), z=(-cellOffset.z)}
  let li = i64.i32 (getLocalNodeIndex recievingNodeOffset)
  let diag = map2 (\ke_row offset -> ke_row[3*li + offset]) keslices[li] (iota 3)
  in map (*E) diag

def getCellDiagonalContribution [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (nodeIndex :index) (cellOffset :index) :[3]f64 =
  match l
    case 0 -> getElementDiagonalContribution l x nodeIndex cellOffset (copy keslices) cellOffset
    case 1 -> 
      let offsets = (cellOffsetToSubcells 1 cellOffset) :> [2][2][2]index
      let vals = map2_3d (getElementDiagonalContribution l x nodeIndex cellOffset) (copy keslices_l1) offsets
      in vals |> map_2d transpose |> map transpose |> transpose |> map_3d f64.sum |> map_2d f64.sum |> map f64.sum
    case _ -> replicate 3 0

def getNodeDiagonal [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) (nodeIndex :index) :[3]f64 =
   [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
    {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
    {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
    {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}]
    |> map (getCellDiagonalContribution l x nodeIndex)
    |> transpose
    |> map f64.sum

def assembleDiagonalNew [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3]f64 =
  let ncell = 2**(i64.u8 l)
  let nx = (nelx/ncell)+1
  let ny = (nely/ncell)+1
  let nz = (nelz/ncell)+1
  in #[incremental_flattening(only_inner)]
    tabulate_3d nx ny nz (\i j k -> getNodeDiagonal l x {x=i,y=j,z=k})

entry assembleInverseDiagonalNew [nelx][nely][nelz] (l :u8) (x :[nelx][nely][nelz]f32) :[][][][3]f64 =
  assembleDiagonalNew l x |> map_4d (\x -> 1/x)

-- ==
-- entry: assembleInverseDiagonalNew
-- compiled random input { 0u8 [128][128][128]f32 }
-- compiled random input { 1u8 [128][128][128]f32 }
-- compiled random input { 2u8 [128][128][128]f32 }