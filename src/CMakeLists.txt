add_library(GEN_MC SHARED
            algorithm.cpp
            analyzer.cpp
            apprx_algorithms.cpp
            generator.cpp
            permutation.cpp
            )
target_compile_options(GEN_MC PUBLIC ${comp_flags})