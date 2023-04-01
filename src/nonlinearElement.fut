import "utility"
import "indexUtilities"
import "material"
import "keConstants"

def pt :f64 = 0.5773502691896257645091 -- 1f64/f64.sqrt(3f64)
def quadpoints :[8]isoCoordinates =
  [{xi=pt,eta=pt,zeta=pt},{xi=(-pt),eta=pt,zeta=pt},
  {xi=pt,eta=(-pt),zeta=pt},{xi=(-pt),eta=(-pt),zeta=pt},
  {xi=pt,eta=pt,zeta=(-pt)},{xi=(-pt),eta=pt,zeta=(-pt)},
  {xi=pt,eta=(-pt),zeta=(-pt)},{xi=(-pt),eta=(-pt),zeta=(-pt)}]

def computeElementSize (nely :i64) =
  1f64 / (f64.i64 nely) -- 1000mm total, divided into nely elements

def computeJacobiDeterminant (esize :f64) = 0.125*esize*esize*esize
def computeJacobiInverse (esize :f64) = 2/esize

def getNodalShapeDiff (iso :isoCoordinates) =
  let t1 = 1 + iso.eta
  let t2 = 1 - iso.zeta
  let t3 = 0.1e1 / 0.8e1 - iso.xi / 8
  let t4 = t3 * t2
  let t5 = t1 / 8
  let t6 = t5 * t2
  let t1 = t3 * t1
  let t7 = 1 + iso.xi
  let t8 = t5 * t7
  let t7 = t7 / 8
  let t9 = t7 * t2
  let t10 = 1 - iso.eta
  let t11 = t7 * t10
  let t12 = t10 / 8
  let t2 = t12 * t2
  let t10 = t3 * t10
  let t13 = 1 + iso.zeta
  let t3 = t3 * t13
  let t5 = t5 * t13
  let t7 = t7 * t13
  let t12 = t12 * t13
  in [[-t6,t4,-t1],[t6,t9,-t8],[t2,-t9,-t11],[-t2,-t4,-t10],[-t5,t3,t1],[t5,t7,t8],[t12,-t7,t11],[-t12,-t3,t10]]

def getDmatrix (Jinv :f64) (u :[24]f64) (nodeShapeDiff :[8][3]f64) =
  let u_nodal = unflatten 8 3 u
  in matmul_f64 (transpose nodeShapeDiff) u_nodal
    |> map_2d (*Jinv) -- multiply with J^-1

def getPhi (Jinv :f64) (u :[24]f64) (nodeShapeDiff :[8][3]f64) = 
  getDmatrix Jinv u nodeShapeDiff
    |> flatten_to 9

def getDeformationGradient (Jinv :f64) (u :[24]f64) (nodeShapeDiff :[8][3]f64) :[3][3]f64 =
  getDmatrix Jinv u nodeShapeDiff
    |> transpose
    |> map2_2d (+) (copy I33) -- +I

def rho0_nl :f64 = 0.1
def beta_nl :f64 = 500

def getNonlinearFraction (x :f64) :f64 =
  let reg = f64.tanh(beta_nl*rho0_nl)
  in (reg + f64.tanh(beta_nl*(x-rho0_nl))) / (reg + f64.tanh(beta_nl*(1-rho0_nl)))

def sech x = 1 / f64.cosh(x)

def getNonlinearFractionDerivative (x :f64) :f64 =
  let s = sech(beta_nl*(x-rho0_nl))
  let a = beta_nl*s*s
  let b = f64.tanh(beta_nl*rho0_nl) + f64.tanh(beta_nl*(1-rho0_nl))
  in a / b

def getB0 (esize :f64) (iso :isoCoordinates) :[6][24]f64 =
  let t1 = 1 - iso.zeta
  let t2 = 1 + iso.zeta
  let t3 = 0.1e1 / esize / 4
  let t4 = t3 * (1 + iso.eta)
  let t5 = t4 * t1
  let t6 = t4 * t2
  let t7 = t3 * (1 - iso.eta)
  let t8 = t7 * t1
  let t9 = t7 * t2
  let t10 = 1 - iso.xi
  let t11 = 1 + iso.xi
  let t12 = t3 * t10
  let t13 = t12 * t1
  let t12 = t12 * t2
  let t3 = t3 * t11
  let t1 = t3 * t1
  let t2 = t3 * t2
  let t3 = t4 * t10
  let t10 = t7 * t10
  let t4 = t4 * t11
  let t7 = t7 * t11
  in [[-t5,0,0,t5,0,0,t8,0,0,-t8,0,0,-t6,0,0,t6,0,0,t9,0,0,-t9,0,0],
      [0,t13,0,0,t1,0,0,-t1,0,0,-t13,0,0,t12,0,0,t2,0,0,-t2,0,0,-t12,0],
      [0,0,-t3,0,0,-t4,0,0,-t7,0,0,-t10,0,0,t3,0,0,t4,0,0,t7,0,0,t10],
      [0,-t3,t13,0,-t4,t1,0,-t7,-t1,0,-t10,-t13,0,t3,t12,0,t4,t2,0,t7,-t2,0,t10,-t12],
      [-t3,0,-t5,-t4,0,t5,-t7,0,t8,-t10,0,-t8,t3,0,-t6,t4,0,t6,t7,0,t9,t10,0,-t9],
      [t13,-t5,0,t1,t5,0,-t1,t8,0,-t13,-t8,0,t12,-t6,0,t2,t6,0,-t2,t9,0,-t12,-t9,0]]

def getG (esize :f64) (iso :isoCoordinates) :[9][24]f64 =
  let t1 = 1 - iso.zeta
  let t2 = 1 + iso.zeta
  let t3 = 0.1e1 / esize / 4
  let t4 = t3 * (1 + iso.eta)
  let t5 = t4 * t1
  let t6 = t4 * t2
  let t7 = t3 * (1 - iso.eta)
  let t8 = t7 * t1
  let t9 = t7 * t2
  let t10 = 1 - iso.xi
  let t11 = 1 + iso.xi
  let t12 = t3 * t10
  let t13 = t12 * t1
  let t12 = t12 * t2
  let t3 = t3 * t11
  let t1 = t3 * t1
  let t2 = t3 * t2
  let t3 = t4 * t10
  let t10 = t7 * t10
  let t4 = t4 * t11
  let t7 = t7 * t11
  in [[-t5,0,0,t5,0,0,t8,0,0,-t8,0,0,-t6,0,0,t6,0,0,t9,0,0,-t9,0,0],
      [0,-t5,0,0,t5,0,0,t8,0,0,-t8,0,0,-t6,0,0,t6,0,0,t9,0,0,-t9,0],
      [0,0,-t5,0,0,t5,0,0,t8,0,0,-t8,0,0,-t6,0,0,t6,0,0,t9,0,0,-t9],
      [t13,0,0,t1,0,0,-t1,0,0,-t13,0,0,t12,0,0,t2,0,0,-t2,0,0,-t12,0,0],
      [0,t13,0,0,t1,0,0,-t1,0,0,-t13,0,0,t12,0,0,t2,0,0,-t2,0,0,-t12,0],
      [0,0,t13,0,0,t1,0,0,-t1,0,0,-t13,0,0,t12,0,0,t2,0,0,-t2,0,0,-t12],
      [-t3,0,0,-t4,0,0,-t7,0,0,-t10,0,0,t3,0,0,t4,0,0,t7,0,0,t10,0,0],
      [0,-t3,0,0,-t4,0,0,-t7,0,0,-t10,0,0,t3,0,0,t4,0,0,t7,0,0,t10,0],
      [0,0,-t3,0,0,-t4,0,0,-t7,0,0,-t10,0,0,t3,0,0,t4,0,0,t7,0,0,t10]]

def getA (Phi :[9]f64) :[6][9]f64 =
  #[unsafe] ([[Phi[0],Phi[1],Phi[2],0,0,0,0,0,0],
              [0,0,0,Phi[3],Phi[4],Phi[5],0,0,0],
              [0,0,0,0,0,0,Phi[6],Phi[7],Phi[8]],
              [0,0,0,Phi[6],Phi[7],Phi[8],Phi[3],Phi[4],Phi[5]],
              [Phi[6],Phi[7],Phi[8],0,0,0,Phi[0],Phi[1],Phi[2]],
              [Phi[3],Phi[4],Phi[5],Phi[0],Phi[1],Phi[2],0,0,0]])

def getMmatrix (s :[6]f64) =
  #[unsafe] [[s[0],0,0,s[5],0,0,s[4],0,0],
             [0,s[0],0,0,s[5],0,0,s[4],0],
             [0,0,s[0],0,0,s[5],0,0,s[4]],
             [s[5],0,0,s[1],0,0,s[3],0,0],
             [0,s[5],0,0,s[1],0,0,s[3],0],
             [0,0,s[5],0,0,s[1],0,0,s[3]],
             [s[4],0,0,s[3],0,0,s[2],0,0],
             [0,s[4],0,0,s[3],0,0,s[2],0],
             [0,0,s[4],0,0,s[3],0,0,s[2]]]

def getQuadraturePointEnergy (esize :f64) (u :[24]f64) iso =
  let nodeShapeDiff = getNodalShapeDiff iso
  let Jinv = computeJacobiInverse esize
  let F = getDeformationGradient Jinv u nodeShapeDiff
  let Fdet = det33_f64 F
  let green = matmul_f64 (transpose F) F
  in neoHookeanEnergy green Fdet

-- st. vernaint
-- def getQuadraturePointNonlinearResidual (u :[24]f64) (iso :isoCoordinates) :[24]f64 =
--   let B0 = getB0 iso
--   let BL = getBL iso u
--   let eps_0 = vecmul_f64 B0 u
--   let eps_L = vecmul_f64 BL u
--   let eps = map2 (\e0 eL -> e0 + 0.5*eL) eps_0 eps_L
--   let sigma = vecmul_f64 C_linear eps
--   let Bbar = matsum_f64 B0 BL
--   in vecmul_f64 (transpose Bbar) sigma

-- neo-hookean
def getQuadraturePointNonlinearResidual (esize :f64) (u :[24]f64) iso :[24]f64 =
  let G = getG esize iso
  let B0 = getB0 esize iso
  let nodeShapeDiff = getNodalShapeDiff iso
  let Jinv = computeJacobiInverse esize
  let phi = getPhi Jinv u nodeShapeDiff
  let A = getA phi
  let BL = matmul_f64 A G
  let Bbar = matsum_f64 B0 BL
  let F = getDeformationGradient Jinv u nodeShapeDiff
  let C = matmul_f64 (transpose F) F
  let Cinv = inv33_f64 C
  let Fdet = det33_f64 F
  let sigma = getNeoHookeanStress Cinv Fdet
    |> toVoigtVector
  in vecmul_f64 (transpose Bbar) sigma

def getQuadraturePointStiffnessMatrix (esize :f64) (u :[24]f64) iso :[24][24]f64 =
  let G = getG esize iso
  let B0 = getB0 esize iso
  let nodeShapeDiff = getNodalShapeDiff iso
  let Jinv = computeJacobiInverse esize
  let F = getDeformationGradient Jinv u nodeShapeDiff
  let phi = getPhi Jinv u nodeShapeDiff
  let A = getA phi
  let BL = matmul_f64 A G
  let green = matmul_f64 (transpose F) F
  let greenInv = inv33_f64 green
  let Fdet = det33_f64 F
  let stress = getNeoHookeanStress greenInv Fdet 
    |> toVoigtVector
  let M = getMmatrix stress
  let inter = matmul_f64 (transpose G) M
  let K_stress = matmul_f64 inter G
  let Bbar = matsum_f64 B0 BL
  let C = getNeoHookeanComplianceTensor greenInv Fdet
    |> toVoigtMatrix
  let inter_geom = matmul_f64 (transpose Bbar) C
  let K_geom = matmul_f64 inter_geom Bbar
  in matsum_f64 K_stress K_geom


entry getElementEnergy (esize :f64) (u :[24]f64) =
  let JDet = computeJacobiDeterminant esize
  in map (getQuadraturePointEnergy esize u) quadpoints
    |> f64.sum
    |> (*JDet)

entry getNonlinearResidual (esize :f64) (u :[24]f64) :[24]f64 =
  let JDet = computeJacobiDeterminant esize
  in #[sequential_outer] map (getQuadraturePointNonlinearResidual esize u) quadpoints
    |> transpose
    |> map f64.sum
    |> map (*JDet)

entry getElementNonlinearTangentMatrix (esize :f64) (u :[24]f64) :[24][24]f64 =
  let JDet = computeJacobiDeterminant esize
  in #[sequential_outer] map (getQuadraturePointStiffnessMatrix esize u) quadpoints
    |> transpose
    |> map transpose
    |> map_2d f64.sum
    |> map_2d (*JDet)

-- ==
-- entry: getElementNonlinearStiffnessMatrix
-- nobench input @../testData/kelinear_input.txt output @../testData/kelinear.txt

def getElementRegularizedEnergy (esize :f64) (u :[24]f64) (x :f64) =
  let young = getYoungsModule x
  let frac  = getNonlinearFraction x
  let u_nl  = map (*frac) u
  let e_nl  = getElementEnergy esize u_nl
  let e_l   =  (*0.5*esize) (f64.sum (map2 (*) (vecmul_f64 keconst_unscaled u) u))
  let e_l2  =  (*0.5*esize) (f64.sum (map2 (*) (vecmul_f64 keconst_unscaled u_nl) u_nl))
  in young*(e_nl+e_l-e_l2)

def getElementRegularizedNonlinearResidual (esize :f64) (u :[24]f64) (x :f64) =
  let young = getYoungsModule x
  let frac  = getNonlinearFraction x
  let u_nl  = map (*frac) u
  let r_nl  = getNonlinearResidual esize u_nl 
    |> map  (*frac)
  let r_l   = vecmul_f64 keconst_unscaled u
    |> map  (*esize)
    |> map  (*(1-frac*frac))
  let r     = map2 (+) r_nl r_l
  in map (*young) r   -- scale for youngs module

def getElementRegularizedNonlinearStiffness (esize :f64) (u :[24]f64) (x :f64) =
  let young     = getYoungsModule x
  let frac      = getNonlinearFraction x
  let u_nl      = map (*frac) u
  let ke_nonlin = getElementNonlinearTangentMatrix esize u_nl
    |> map_2d (*(frac*frac))
  let ke_lin    = (copy keconst_unscaled)
    |> map_2d (*esize)
    |> map_2d (*(1-frac*frac))
  let ke = map2_2d (+) ke_nonlin ke_lin
  in map_2d (*young) ke

def getRegularizedEnergy [nx][ny][nz][nelx][nely][nelz] (esize :f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (elementIndex :index) =
  let uloc = getLocalState u elementIndex
  let xloc = getDensityUnsafe x elementIndex
  in getElementRegularizedEnergy esize uloc xloc

def getRegularizedResidual [nx][ny][nz][nelx][nely][nelz] (esize :f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (elementIndex :index) :[24]f64 =
  let uloc = readLocalState u elementIndex
  let xloc = getDensityUnsafe x elementIndex
  in getElementRegularizedNonlinearResidual esize uloc xloc

def getRegularizedTangent [nx][ny][nz][nelx][nely][nelz] (esize :f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (elementIndex :index) :[24][24]f64 =
  let uloc = readLocalState u elementIndex
  let xloc = getDensityUnsafe x elementIndex
  in getElementRegularizedNonlinearStiffness esize uloc xloc


def getLinearisedStiffnessProduct [nx][ny][nz][nelx][nely][nelz] (esize :f64) (u_expansion :[nx][ny][nz][3]f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (elementIndex :index) :[24]f64 =
  let uloc = readLocalState u elementIndex
  let kt   = getRegularizedTangent esize x u_expansion elementIndex
  in vecmul_f64 kt uloc

def getElementStiffnessDiagonal [nx][ny][nz][nelx][nely][nelz] (esize :f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (elementIndex :index) :[24]f64 =
  let kt = getRegularizedTangent esize x u elementIndex
  in map (\i -> #[unsafe] kt[i,i]) (iota 24)

def getElementStiffnessBlockDiagonal [nx][ny][nz][nelx][nely][nelz] (esize :f64) (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) (elementIndex :index) :[8][3][3]f64 =
  let kt = getRegularizedTangent esize x u elementIndex
  in tabulate 8 (\n -> 
    tabulate_2d 3 3 (\i j -> 
      let ii = 3*n + i
      let jj = 3*n + j
      in #[unsafe] kt[ii,jj]
    ))