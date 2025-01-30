include("Resource_Model/DA_birmmals_with_pi.jl")
include("Resource_Model/Functions/generate_competition_matrix.jl")
include("Resource_Model/Functions/species_dict.jl")

include("Resource_Model/ecosystem_dynamics!.jl")
include("Resource_Model/FI_functions.jl")
include("Resource_Model/extract_H0_DA.jl")
include("Resource_Model/Functions/attempt_setup_community.jl")
include("Resource_Model/Functions/Callbacks_function.jl")
include("Resource_Model/npp_DA_relative_to_1000.jl")
include("Resource_Model/Functions/Computing_metrics.jl")
include("HerbivoresVsPredators/Exploring HerbPred metaweb composition.jl")
include("Resource_Model/Functions/attempt_feasibility.jl")
global EXTINCTION_THRESHOLD = 1e-6
global T_ext               = 250.0
include("Resource_Model/SingleModelRun.jl")
include("Resource_Model/species_distributions.jl")
