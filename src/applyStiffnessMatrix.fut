import "keConstants"
import "indexUtilities"
import "boundaryConditions"
import "projection"
import "utility"
import "material"

def getSendingNode (elementOffset :index, nodeOffset :index) :(i32, i32) =
  let recievingNodeOffset :index = {x=(-elementOffset.x), y=(-elementOffset.y), z=(-elementOffset.z)}
  let sendingNodeOffset :index = {x=nodeOffset.x-elementOffset.x, y=nodeOffset.y-elementOffset.y, z=nodeOffset.z-elementOffset.z}
  in (getLocalNodeIndex(recievingNodeOffset), getLocalNodeIndex(sendingNodeOffset))

def getInputVector [nx][ny][nz] (nodeIndex :index, nodeOffset :index, u :[nx][ny][nz][3]f64) :[3]f64 =
  let loadIndex :index = {x=nodeIndex.x+nodeOffset.x,y=nodeIndex.y+nodeOffset.y,z=nodeIndex.z+nodeOffset.z} in
  tabulate 3 (\d ->
    if ((indexIsInside (nx,ny,nz) loadIndex) && !(isOnBoundary (nx,ny,nz) loadIndex d)) then
      #[unsafe] (u[loadIndex.x,loadIndex.y,loadIndex.z,d])
    else
      0)

def getInputVectorSingle [nx][ny][nz] (nodeIndex :index, nodeOffset :index, u :[nx][ny][nz][3]f32) :[3]f32 =
  let loadIndex :index = {x=nodeIndex.x+nodeOffset.x,y=nodeIndex.y+nodeOffset.y,z=nodeIndex.z+nodeOffset.z} in
  tabulate 3 (\d ->
    if ((indexIsInside (nx,ny,nz) loadIndex) && !(isOnBoundary (nx,ny,nz) loadIndex d)) then
      #[unsafe] (u[loadIndex.x,loadIndex.y,loadIndex.z,d])
    else
      0)

def multiplyScaledLocalMatrix(m :localMatrix, a :[3]f64, s :f64) :[3]f64 =
  [(s*(m.xx*a[0]+m.xy*a[1]+m.xz*a[2])),
   (s*(m.yx*a[0]+m.yy*a[1]+m.yz*a[2])),
   (s*(m.zx*a[0]+m.zy*a[1]+m.zz*a[2]))]

def scaleLocalMatrix(m :localMatrix) (s :f64) :localMatrix =
 {xx=s*m.xx,xy=s*m.xy,xz=s*m.xz,
  yx=s*m.yx,yy=s*m.yy,yz=s*m.yz,
  zx=s*m.zx,zy=s*m.zy,zz=s*m.zz}

def addLocalMatrix(a :localMatrix) (b :localMatrix) :localMatrix =
 {xx=a.xx+b.xx,xy=a.xy+b.xy,xz=a.xz+b.xz,
  yx=a.yx+b.yx,yy=a.yy+b.yy,yz=a.yz+b.yz,
  zx=a.zx+b.zx,zy=a.zy+b.zy,zz=a.zz+b.zz}

def getLocalMatrix (elementOffset :index, nodeOffset :index) :localMatrix =
  let (recieve,send) = getSendingNode(elementOffset, nodeOffset)
  in getke_l0 (recieve,send)

def getStencilSpoke [nx][ny][nz] (nodeIndex :index) (u :[nx][ny][nz][3]f64) (elementOffset :index) (elementScale :f64) (nodeOffset :index) =
  let localMatrix = getLocalMatrix(elementOffset,nodeOffset)
  let inputVector = getInputVector(nodeIndex,nodeOffset,u)
  in multiplyScaledLocalMatrix(localMatrix, inputVector, elementScale)

def getElementContribution [nx][ny][nz][nelx][nely][nelz] (nodeIndex :index, u :[nx][ny][nz][3]f64, x :[nelx][nely][nelz]f32) (elementOffset :index) :[3]f64 =
  let eo = elementOffset
  let elementIndex = addIndices nodeIndex elementOffset
  let E = getElementYoungsModule x elementIndex
  in [{x=(eo.x+0),y=(eo.y+0),z=(eo.z+0)}, {x=(eo.x+1),y=(eo.y+0),z=(eo.z+0)},
     {x=(eo.x+0),y=(eo.y+1),z=(eo.z+0)}, {x=(eo.x+1),y=(eo.y+1),z=(eo.z+0)},
     {x=(eo.x+0),y=(eo.y+0),z=(eo.z+1)}, {x=(eo.x+1),y=(eo.y+0),z=(eo.z+1)},
     {x=(eo.x+0),y=(eo.y+1),z=(eo.z+1)}, {x=(eo.x+1),y=(eo.y+1),z=(eo.z+1)}]
     |> map (getStencilSpoke nodeIndex u eo E)
     |> transpose
     |> map (f64.sum)

def getElementContributionAlt [nx][ny][nz][nelx][nely][nelz] (nodeIndex :index, u :[nx][ny][nz][3]f64, x :[nelx][nely][nelz]f32) (elementOffset :index) :[3]f64 =
  let elementIndex = addIndices nodeIndex elementOffset
  let ulocal = getLocalState u elementIndex
  let E = getElementYoungsModule x elementIndex
  let recievingNodeOffset :index = {x=(-elementOffset.x), y=(-elementOffset.y), z=(-elementOffset.z)}
  let li = i64.i32 (getLocalNodeIndex recievingNodeOffset)
  in #[unsafe] map (\ke_row -> (f64.sum (map2(\ke u -> ke*E*u) ke_row ulocal))) keslices[li]

def applyStencilOnNode [nelx][nely][nelz][nx][ny][nz] (nodeIndex :index, u :[nx][ny][nz][3]f64, x :[nelx][nely][nelz]f32) :[3]f64 =
   [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
    {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
    {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
    {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}]
    |> map (\eo -> getElementContributionAlt (nodeIndex, u, x) eo)
    |> transpose
    |> map (\x -> f64.sum x)

def computeElementSize (nely :i64) =
  1f64 / (f64.i64 nely) -- 1000mm total, divided into nely elements

-- returns the matrix-vector product K(x)*u
#[noinline]
def applyStiffnessMatrix [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let es = computeElementSize nely
  let f = tabulate_3d nx ny nz (\i j k -> applyStencilOnNode({x=i,y=j,z=k}, u, x))
  let f = map_4d (*es) f
  in setBCtoInput u f

-- returns the matrix-vector product K(x)*u for f32
def applyStiffnessMatrixSingle [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) :[nx][ny][nz][3]f32 =
  map_4d f64.f32 u
  |> applyStiffnessMatrix x
  |> map_4d f32.f64

-- returns the matrix-vector product (P^T*K(x)*P)*u for any subspace
def applyCoarseStiffnessMatrix [nelx][nely][nelz][nxc][nyc][nzc] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) :[nxc][nyc][nzc][3]f64 =
  let ufine =
    #[unroll] loop u for i < (i16.u8 l) do
      projectToFiner u

  let rfine = applyStiffnessMatrix x ufine

  let rcoarse =
    #[unroll] loop r = rfine for i < (i16.u8 l) do
      projectToCoarser r

  in (rcoarse :> [nxc][nyc][nzc][3]f64)

entry applyStiffnessMatrix_test [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  applyStiffnessMatrix x u

entry applyCoarseStiffnessMatrix_test [nelx][nely][nelz][nxc][nyc][nzc] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix l x u

entry applyCoarseStiffnessMatrix_test0 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 0 x u

entry applyCoarseStiffnessMatrix_test1 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 1 x u

entry applyCoarseStiffnessMatrix_test2 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 2 x u

entry applyCoarseStiffnessMatrix_test3 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 3 x u

-- ==
-- entry: applyStiffnessMatrix_test
-- nobench input @../testData/applyStateOperator1.txt output @../testData/applyStateOperator1Output.txt
-- nobench input @../testData/applyStateOperator2.txt output @../testData/applyStateOperator2Output.txt
-- nobench input @../testData/applyStateOperator3.txt output @../testData/applyStateOperator3Output.txt
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f64 }
-- compiled random input { [256][256][256]f32 [257][257][257][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test
-- nobench input @../testData/applyStateOperatorFold1.txt output @../testData/applyStateOperatorFold1Out.txt
-- nobench input @../testData/applyStateOperatorFold2.txt output @../testData/applyStateOperatorFold2Out.txt
-- nobench input @../testData/applyStateOperatorFold3.txt output @../testData/applyStateOperatorFold3Out.txt
-- nobench input @../testData/applyStateOperatorFold4.txt output @../testData/applyStateOperatorFold4Out.txt

-- ==
-- entry: applyCoarseStiffnessMatrix_test0
-- compiled random input { [128][128][128]f32 [129][129][129][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test1
-- compiled random input { [128][128][128]f32 [65][65][65][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test2
-- compiled random input { [128][128][128]f32 [33][33][33][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test3
-- compiled random input { [128][128][128]f32 [17][17][17][3]f64 }
