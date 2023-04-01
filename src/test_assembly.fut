import "assembleNonlinear"


-- #[noinline]
entry nonlinearResidual = assembleNonlinearResidual
entry linearisedTangent = applyLinearisedStiffnessMatrix

-- ==
-- entry: nonlinearResidual
-- compiled random input { [32][32][32]f32 [33][33][33][3]f64 } 
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 }

-- ==
-- entry: linearisedTangent
-- compiled random input { [32][32][32]f32 [33][33][33][3]f64 [33][33][33][3]f64 } 
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 [65][65][65][3]f64 } 



