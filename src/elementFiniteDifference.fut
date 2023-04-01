import "utility"
import "material"
import "nonlinearElement"
import "indexUtilities"

-- finite difference utilities
def pert :f64 = 1e-6

def elementSize :f64 = 1.24

def finiteDiffPoint [n] f (v :[n]f64) (i :i64) = 
  let x1 = (copy v) with [i] = v[i] + pert
  let x2 = (copy v) with [i] = v[i] - pert
  let f1 = f x1
  let f2 = f x2
  in (f1-f2) / (2*pert)

def finiteDiff [n] f (v :[n]f64) = 
  tabulate n (finiteDiffPoint f v)

def finiteDiffPoint_2d [n][m] f (v :[n][m]f64) (i :i64) (j :i64) = 
  let x1 = (copy v) with [i,j] = v[i,j] + pert
  let x2 = (copy v) with [i,j] = v[i,j] - pert
  let f1 = f x1
  let f2 = f x2
  in (f1-f2) / (2*pert)

def finiteDiff_2d [n][m] f (v :[n][m]f64) = 
  tabulate_2d n m (finiteDiffPoint_2d f v)

def finiteDiffPoint_22 [n][m] f (v :[n][m]f64) (i :i64) (j :i64) = 
  let x1 = (copy v) with [i,j] = v[i,j] + pert
  let x2 = (copy v) with [i,j] = v[i,j] - pert
  in map2_2d (\c1 c2 -> (c1-c2) / (2*pert)) (f x1) (f x2) 

def finiteDiff_22 [n][m] f (v :[n][m]f64) = 
  tabulate_2d n m (finiteDiffPoint_22 f v)

def finiteDiffPoint_11 [n] f (v :[n]f64) (i :i64) = 
  let x1 = (copy v) with [i] = v[i] + pert
  let x2 = (copy v) with [i] = v[i] - pert
  in map2 (\c1 c2 -> (c1-c2) / (2*pert)) (f x1) (f x2) 

def finiteDiff_11 [n] f (v :[n]f64) = 
  tabulate n (finiteDiffPoint_11 f v)

def greenFromStrain =
  map2_2d (\i e->2*e+i) I33

-- finite difference on material level
def energyFromStrain (E :[3][3]f64) =
  let C = greenFromStrain E
  let Fdet = f64.sqrt(det33_f64 C)
  in neoHookeanEnergy C Fdet

def stressFromStrain (E :[3][3]f64) = 
  let C = greenFromStrain E
  let Fdet = f64.sqrt(det33_f64 C)
  let Cinv = inv33_f64 C
  in getNeoHookeanStress Cinv Fdet

def stressFiniteDifference (E :[3][3]f64) =
    finiteDiff_2d energyFromStrain E

entry stressFiniteDifferenceDiff (E :[3][3]f64) = 
    let stress    = stressFromStrain E
    let stress_fd = stressFiniteDifference E
    in map2_2d (-) stress stress_fd
        |> flatten
        |> map f64.abs
        |> f64.sum

-- ==
-- entry: stressFiniteDifferenceDiff
-- nobench input @../testData/strain1.txt output { 0f64 }
-- nobench input @../testData/strain2.txt output { 0f64 }

def complianceFromStrain E =
    let C = greenFromStrain E
    let Fdet = f64.sqrt(det33_f64 C)
    let Cinv = inv33_f64 C
    in getNeoHookeanComplianceTensor Cinv Fdet

def complianceFiniteDifference = 
    finiteDiff_22 stressFromStrain

entry complianceFiniteDifferenceDiff (E :[3][3]f64) =
    let C = complianceFromStrain E
    let C_fd = complianceFiniteDifference E
    in map2_4d (-) C C_fd
        |> flatten_4d
        |> map f64.abs
        |> f64.sum

    -- ==
    -- entry: complianceFiniteDifferenceDiff
    -- nobench input @../testData/strain1.txt output { 0f64 }
    -- nobench input @../testData/strain2.txt output { 0f64 }

entry testSym1 (E :[3][3]f64) = 
  let C = complianceFromStrain E
  let sym c = map2_2d (-) c (transpose c)
  in map_2d sym C 
    |> flatten_4d 
    |> map f64.abs
    |> f64.sum

-- ==
-- entry: testSym1
-- nobench input @../testData/strain1.txt output { 0f64 }
-- nobench input @../testData/strain2.txt output { 0f64 }

entry testSym2 (E :[3][3]f64) = 
  let C = complianceFromStrain E
  in map2_4d (-) C (transpose C)
    |> flatten_4d 
    |> map f64.abs
    |> f64.sum

-- ==
-- entry: testSym2
-- nobench input @../testData/strain1.txt output { 0f64 }
-- nobench input @../testData/strain2.txt output { 0f64 }


def applyCtensor (E :[3][3]f64) = 
  let C = complianceFromStrain E 
  let s = map_2d (\cij -> f64.sum (flatten (map2_2d (*) cij E))) C 
  in toVoigtVector s
    
def applyCmatrix (E :[3][3]f64) = 
  let Cv = complianceFromStrain E |> toVoigtMatrix
  let Ev = toVoigtVector E
  let Evv = [Ev[0],Ev[1],Ev[2],2*Ev[3],2*Ev[4],2*Ev[5]]
  in vecmul_f64 Cv Evv

entry testVoigt (E :[3][3]f64) =
  let a = applyCtensor E
  let b = applyCmatrix E
  in map2 (-) a b
    |> map f64.abs
    |> f64.sum
  
-- ==
-- entry: testVoigt
-- nobench input @../testData/strain1.txt output { 0f64 }
-- nobench input @../testData/strain2.txt output { 0f64 }


-- finite difference on element level
def residualFiniteDifference =
    finiteDiff (getElementEnergy elementSize)

entry residualFiniteDifferenceDiff (u :[24]f64) =
    let r = getNonlinearResidual elementSize u
    let r_fd = residualFiniteDifference u
    in map2 (-) r r_fd
      |> map f64.abs
      |> f64.sum

-- ==
-- entry: residualFiniteDifferenceDiff
-- nobench random input { [24]f64 } output { 0f64 }

-- finite difference on element level
def regularizedResidualFiniteDifference (x :f64) =
    finiteDiff (\v -> getElementRegularizedEnergy elementSize v x)
    
entry regularizedResidualFiniteDifferenceDiff (u :[24]f64) (x :f64) =
    let r = getElementRegularizedNonlinearResidual elementSize u x
    let r_fd = regularizedResidualFiniteDifference x u
    in map2 (-) r r_fd
      |> map f64.abs
      |> f64.sum

-- ==
-- entry: regularizedResidualFiniteDifferenceDiff
-- nobench random input { [24]f64 0.0f64 } output { 0f64 }
-- nobench random input { [24]f64 0.09f64 } output { 0f64 }
-- nobench random input { [24]f64 0.1f64 } output { 0f64 }
-- nobench random input { [24]f64 0.11f64 } output { 0f64 }
-- nobench random input { [24]f64 0.5f64 } output { 0f64 }
-- nobench random input { [24]f64 0.75f64 } output { 0f64 }
-- nobench random input { [24]f64 1.0f64 } output { 0f64 }

def tangentFiniteDifference = 
  finiteDiff_11 (getNonlinearResidual elementSize)

entry tangentFiniteDifferenceDiff (u :[24]f64) =
    let kt = getElementNonlinearTangentMatrix elementSize u
    let kt_fd = tangentFiniteDifference u
    in map2_2d (-) kt kt_fd
         |> flatten
         |> map f64.abs
         |> f64.sum

-- ==
-- entry: tangentFiniteDifferenceDiff
-- nobench random input { [24]f64 } output { 0f64 }


def regularizedTangentFiniteDifference (x :f64) = 
  finiteDiff_11 (\v -> getElementRegularizedNonlinearResidual elementSize v x)

entry regularizedTangentFiniteDifferenceDiff (u :[24]f64) (x :f64) =
    let kt = getElementRegularizedNonlinearStiffness elementSize u x
    let kt_fd = regularizedTangentFiniteDifference x u
    in map2_2d (-) kt kt_fd
         |> flatten
         |> map f64.abs
         |> f64.sum

-- ==
-- entry: regularizedTangentFiniteDifferenceDiff
-- nobench random input { [24]f64 0.0f64 } output { 0f64 }
-- nobench random input { [24]f64 0.09f64 } output { 0f64 }
-- nobench random input { [24]f64 0.1f64 } output { 0f64 }
-- nobench random input { [24]f64 0.11f64 } output { 0f64 }
-- nobench random input { [24]f64 0.5f64 } output { 0f64 }
-- nobench random input { [24]f64 0.75f64 } output { 0f64 }
-- nobench random input { [24]f64 1.0f64 } output { 0f64 }