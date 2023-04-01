import "indexUtilities"
import "keConstants"
import "utility"

type isoCoordinates = {xi: f64, eta: f64, zeta: f64}

def E    :f64 = 1 -- 3GPa, Nylon 
def Emin :f64 = 1e-6*E

def E_f32    :f32 = f32.f64 E
def Emin_f32 :f32 = f32.f64 Emin

def getYoungsModule (x :f64) = Emin + (E-Emin)*x*x*x
def getYoungsModuleDerivative (x :f64) = 3*(E-Emin)*x*x

def nu     :f64 = 0.3 -- Nylon
def lambda :f64 = nu / ((1+nu)*(1-2*nu))
def mu     :f64 = 1 / (2*(1+nu))

def getElementYoungsModule [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (elementIndex :index) :f64 =
  if (indexIsInside (nelx,nely,nelz) elementIndex) then
    let xloc :f64 = getDensityUnsafe x elementIndex
    in getYoungsModule xloc
  else 0

def keprod (s :f64) (f :*[24]f64) :*[24]f64 =
  map (\keRow -> reduce (+) 0 (map2(\ff kr -> ff*kr*s) f keRow)) keconst

def I33 :[3][3]f64 = [[1,0,0],[0,1,0],[0,0,1]]

def neoHookeanEnergy (C :[3][3]f64) (Fdet :f64) = 
  let lnJ = f64.log(Fdet)
  let trC = C[0,0]+C[1,1]+C[2,2]
  in 0.5*lambda*lnJ*lnJ + 0.5*mu*(trC-3) - mu*lnJ

-- Second P-K stress
def getNeoHookeanStress (Cinv :[3][3]f64) (Fdet :f64) = 
  map2_2d (\i cinv -> lambda*f64.log(Fdet)*cinv + mu*(i-cinv)) (copy I33) Cinv

def getNeoHookeanComplianceTensor (Cinv :[3][3]f64) (Fdet :f64) =
  let a = map_2d (\Cij -> 
            map_2d (\Ckl -> 
              Cij*Ckl*lambda
            ) Cinv
          ) Cinv
  let bfact = mu-lambda*f64.log(Fdet)
  let b = map (\Cix -> 
            map (\Cjx ->
              map2 (\Cik Cjk -> 
                map2 (\Cil Cjl -> 
                  (Cil*Cjk + Cik*Cjl)*bfact
                ) Cix Cjx
              ) Cix Cjx
            ) Cinv 
          ) Cinv
  in map2_4d (+) a b