
add_executable(experiment experiment.cpp)
target_link_libraries(experiment PUBLIC GEN_MC)

add_executable(max_pen max_pen.cpp)
target_link_libraries(max_pen PUBLIC GEN_MC)

add_executable(single_instance single_instance.cpp)
target_link_libraries(single_instance PUBLIC GEN_MC)