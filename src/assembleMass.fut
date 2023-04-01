import "indexUtilities"
import "utility"

def elementOffsets :[8]index =
  [{x=(-1),y=( 0),z=(-1)}, {x=( 0),y=( 0),z=(-1)},
   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)},
   {x=(-1),y=( 0),z=( 0)}, {x=( 0),y=( 0),z=( 0)},
   {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)}]

def M    :f64 = 1
def Mmin :f64 = 1e-3

def getElementMass [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (elementIndex :index) :f64 =
  if (indexIsInside (nelx,nely,nelz) elementIndex) then
    let rho :f64 = #[unsafe](f64.f32 x[elementIndex.x,elementIndex.y,elementIndex.z])
    in (Mmin + (M-Mmin)*rho)
  else 0

def getElementLumpedMassContribution [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (nodeIndex :index) (elementOffset :index) :[3]f64 =
  let elementIndex = addIndices nodeIndex elementOffset
  let mass = getElementMass x elementIndex
  in replicate 3 (0.125*mass)

def getNodeLumpedMass [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (nodeIndex :index) :[3]f64 =
  copy elementOffsets
   |> map (getElementLumpedMassContribution x nodeIndex)
   |> transpose
   |> map (reduce (+) 0)

entry assembleLumpedMassMatrix [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :[][][][3]f64 =
  tabulate_3d (nelx+1) (nely+1) (nelz+1) (\i j k -> getNodeLumpedMass x {x=i,y=j,z=k})
