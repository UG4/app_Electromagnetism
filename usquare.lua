--------------------------------------------------------------------------------
--[[!
-- \file apps/electromagnetism/usquare.lua
-- \author Dmitry Logashenko
-- \brief Lua-Script for the test simulation with the electromagnetism module
]]--
--------------------------------------------------------------------------------

ug_load_script ("ug_util.lua")

-- constants
dim        = 2; -- the problem is formulated in 2d (for debugging)
numPreRefs = util.GetParamNumber ("-numPreRefs", 0, "number of refinements before parallel distribution")
numRefs    = util.GetParamNumber ("-numRefs",    1, "number of refinements")

gridName = "grids/usquare.ugx"

print (" Choosen Parameter:")
print ("    numRefs    = " .. numRefs)
print ("    numPreRefs = " .. numPreRefs)

-- Init UG for dimension and algebra
InitUG (dim, AlgebraType("CPU", 2)); -- note: the block size should be 2

--------------------------------------------------------------------------------
-- Domain Setup
--------------------------------------------------------------------------------

-- Create, load, refine and distribute the domain
neededSubsets = {"Inner", "Bottom", "Top", "Left", "Right"}
dom = util.CreateAndDistributeDomain (gridName, numRefs, numPreRefs, neededSubsets)

-- Set the electromagnetic parameters to the subdomains
em = EMaterial (dom)
em:add ("Inner", 1.0, 0.0)
em:close ()

-- Create the edge-centered approximation space for E
print ("--> Edge-centered DoF distribution")
edgeApproxSpace = ApproximationSpace (dom)
edgeApproxSpace:add_fct ("r", "Nedelec", 1) -- Re
edgeApproxSpace:add_fct ("i", "Nedelec", 1) -- Im

edgeApproxSpace:init_levels ()
edgeApproxSpace:init_top_surface ()
edgeApproxSpace:print_statistic ()

-- Create the vertex-centered approximation space for the potential
print ("--> Vertex-centered DoF distribution")
vertApproxSpace = ApproximationSpace (dom, AlgebraType("CPU", 1))
vertApproxSpace:add_fct ("u", "Lagrange", 1)

vertApproxSpace:init_levels ()
vertApproxSpace:init_top_surface ()
vertApproxSpace:print_statistic ()

--------------------------------------------------------------------------------
-- User Data Setup
--------------------------------------------------------------------------------

function orig_E (x, y, t)
	return y, -x
end

--------------------------------------------------------------------------------
-- FE Disc setup
--------------------------------------------------------------------------------

-- Local discretization
elemDisc = EddyCurrent_E_Nedelec ("r,i", em, 1) -- last argument = frequency omega

-- Dirichlet BC
dirichletBC = NedelecDirichletBC ("r,i")
dirichletBC:add ("orig_E", "r", "Bottom")
dirichletBC:add ("orig_E", "r", "Top")
dirichletBC:add ("orig_E", "r", "Left")
dirichletBC:add ("orig_E", "r", "Right")

-- Global discretization
domainDisc = DomainDiscretization (edgeApproxSpace)
domainDisc:add (elemDisc)
domainDisc:add (dirichletBC)

--------------------------------------------------------------------------------
--  Linear solver for the Nedelec-element-based discretization
--------------------------------------------------------------------------------

-- matrix and vectors
A = MatrixOperator ()
u = GridFunction (edgeApproxSpace)
b = GridFunction (edgeApproxSpace)

-- edge-centered smoother in the coarse grid solver
edgeBaseSmoother = GaussSeidel ()

-- vertex-centered smoother in the coarse grid solver (should be based on a scalar algebra)
vertBaseSmoother = GaussSeidelCPU1 ()

-- hybrid smoother in the coarse grid solver
baseHybridSmoother = HiptmairHybridSmoother (vertApproxSpace, edgeBaseSmoother, vertBaseSmoother)
baseHybridSmoother:set_Dirichlet (dirichletBC)

-- convergence check for the coarse solver
baseConvCheck = ConvCheck ()
baseConvCheck:set_maximum_steps (1024)
baseConvCheck:set_minimum_defect (1e-10)
baseConvCheck:set_reduction (1e-10)
baseConvCheck:set_verbose (false)

-- coarse grid solver
baseSolver = BiCGStab ()
baseSolver:set_preconditioner (baseHybridSmoother)
baseSolver:set_convergence_check (baseConvCheck)

-- edge-centered smoother
edgeSmoother = GaussSeidel ()

-- vertex-centered smoother (should be based on a scalar algebra)
vertSmoother = GaussSeidelCPU1 ()

-- hybrid smoother
hybridSmoother = HiptmairHybridSmoother (vertApproxSpace, edgeSmoother, vertSmoother)
hybridSmoother:set_Dirichlet (dirichletBC)

-- transfer operators
transferOp = NedelecTransfer (edgeApproxSpace)

-- geometric multigrid method
gMG = GeometricMultiGrid (edgeApproxSpace)
gMG:set_discretization (domainDisc)
gMG:set_base_level (0)
gMG:set_base_solver (baseSolver)
gMG:set_smoother (hybridSmoother)
gMG:set_transfer (transferOp)
gMG:set_cycle_type (1)
gMG:set_num_presmooth (1)
gMG:set_num_postsmooth (1)

-- convergence check
ConvCheck = ConvCheck ()
ConvCheck:set_maximum_steps (1024)
ConvCheck:set_minimum_defect (1e-10)
ConvCheck:set_reduction (1e-10)
ConvCheck:set_verbose (true)

-- solver for the discretization
linSolver = LinearSolver ()
--linSolver = BiCGStab ()
linSolver:set_preconditioner (gMG)
--linSolver:set_preconditioner (hybridSmoother)
linSolver:set_convergence_check (ConvCheck)

--------------------------------------------------------------------------------
--  Projection
--------------------------------------------------------------------------------

-- vertex-centered smoother for the coarse grid solver of the projection
projBaseSmoother = ILUCPU1 ()

-- convergence check for the coarse grid solver of the projection
projBaseConvCheck = ConvCheckCPU1 ()
projBaseConvCheck:set_maximum_steps (1024)
projBaseConvCheck:set_minimum_defect (1e-10)
projBaseConvCheck:set_reduction (1e-10)
projBaseConvCheck:set_verbose (false)

-- coarse grid solver for the projection
projBaseSolver = CGCPU1 ()
projBaseSolver:set_preconditioner (projBaseSmoother)
projBaseSolver:set_convergence_check (projBaseConvCheck)

-- vertex-centered smoother for the projection
projSmoother = ILUCPU1 ()

-- geometric multigrid method for the projection
projGMG = _G ["GeometricMultiGrid"..dim.."dCPU1"] (vertApproxSpace)
projGMG:set_base_level (1)
projGMG:set_base_solver (projBaseSolver)
projGMG:set_smoother (projSmoother)
projGMG:set_cycle_type (1)
projGMG:set_num_presmooth (1)
projGMG:set_num_postsmooth (1)

-- convergence check for the projection
projConvCheck = ConvCheckCPU1 ()
projConvCheck:set_maximum_steps (1024)
projConvCheck:set_minimum_defect (1e-10)
projConvCheck:set_reduction (1e-10)
projConvCheck:set_verbose (true)

-- linear solver for the projection
projSolver = CGCPU1 ()
projSolver:set_preconditioner (projGMG)
projSolver:set_convergence_check (projConvCheck)

-- projection
projection = NedelecProject (em, vertApproxSpace, projSolver)
projection:set_Dirichlet (dirichletBC)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

print ("--> Solving the discretized problem")

-- 1. assemble matrix and rhs
domainDisc:assemble_linear (A, b)

-- 2. set dirichlet values and start iterate
u:set (0.0)
--ComputeNedelecDoFs ("orig_E", u, "r")
domainDisc:adjust_solution (u)

--SaveVectorForConnectionViewer (u, "InitSol2d.vec")
--SaveVectorForConnectionViewer (b, "RHS2d.vec")
--SaveMatrixForConnectionViewer (u, A, "Matrix2d.mat")

-- 3. init solver for linear Operator
linSolver:init (A, u)

-- 4. apply solver
linSolver:apply_return_defect (u, b)

--SaveVectorForConnectionViewer (u, "Sol2d.vec")

--------------------------------------------------------------------------------
--  Apply projection
--------------------------------------------------------------------------------

print ("--> Projection of the solution")

projection:apply (u, "r,i")

--SaveVectorForConnectionViewer (u, "ProjSol2d.vec")

--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------

print ("--> Output")

out = VTKOutput ()
out:clear_selection ()

-- electric field
ReEData = NedelecGridFunctionData (u, "r")
ImEData = NedelecGridFunctionData (u, "i")
out:select_element (ReEData, "ReE");
out:select_element (ImEData, "ImE");

-- curl of the electric field
ReCurlE = NedelecCurlData (u, "r")
ImCurlE = NedelecCurlData (u, "i")
out:select_element (ReCurlE, "ReCurlE");
out:select_element (ImCurlE, "ImCurlE");

out:print ("Solution2d", u)

-- End of File
