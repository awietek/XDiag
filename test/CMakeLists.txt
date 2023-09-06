add_executable(tests
  tests.cpp
  bitops/test_bitops.cpp
  bitops/test_bits_static.cpp

  io/test_args.cpp
  io/test_file_toml.cpp
  io/test_file_h5.cpp

  combinatorics/test_binomial.cpp
  combinatorics/test_subsets.cpp
  combinatorics/test_bit_patterns.cpp
  combinatorics/test_combinations.cpp
  combinatorics/test_combinations_index.cpp

  indexing/test_lin_table.cpp
  indexing/test_fermi_table.cpp
  indexing/spinhalf/test_spinhalf_indexing_sublattice.cpp
  indexing/spinhalf/test_spinhalf_indexing.cpp

  symmetries/group_action/test_group_action.cpp
  symmetries/group_action/test_group_action_lookup.cpp
  symmetries/group_action/test_group_action_sublattice.cpp
  symmetries/operations/test_symmetry_operations.cpp
  # symmetries/operations/test_representative_list_omp.cpp
  symmetries/test_fermi_sign.cpp
  symmetries/test_permutation.cpp
  symmetries/test_permutation_group.cpp
  symmetries/test_generated_group.cpp
  symmetries/test_qn.cpp

  operators/test_symmetrized_operator.cpp
  operators/test_non_branching_bonds.cpp

  blocks/spinhalf/testcases_spinhalf.cpp
  blocks/spinhalf/test_spinhalf_matrix.cpp
  blocks/spinhalf/test_spinhalf_apply.cpp

  blocks/spinhalf_symmetric/test_spinhalf_symmetric.cpp
  blocks/spinhalf_symmetric/test_spinhalf_symmetric_matrix.cpp
  blocks/spinhalf_symmetric/test_spinhalf_symmetric_apply.cpp
  blocks/spinhalf_symmetric/test_kitaev_gamma.cpp

  blocks/tj/testcases_tj.cpp
  blocks/tj/test_tj_utils.cpp
  blocks/tj/test_tj_matrix.cpp
  blocks/tj/test_tj_apply.cpp
  blocks/tj/test_tj_raiselower.cpp

  blocks/tj_symmetric/test_tj_symmetric.cpp
  blocks/tj_symmetric/test_tj_symmetric_matrix.cpp
  blocks/tj_symmetric/test_tj_symmetric_apply.cpp

  blocks/electron/testcases_electron.cpp
  blocks/electron/test_electron_matrix.cpp
  blocks/electron/test_electron_apply.cpp
  blocks/electron/test_electron_raiselower.cpp

  blocks/electron_symmetric/test_electron_symmetric.cpp
  blocks/electron_symmetric/test_electron_symmetric_matrix.cpp
  blocks/electron_symmetric/test_electron_symmetric_apply.cpp

  algorithms/lanczos/test_lanczos_eigenvalues.cpp
  algorithms/lanczos/test_lanczos_pro.cpp
  algorithms/lanczos/test_lanczos_eigenvector.cpp
  algorithms/arnoldi/test_arnoldi.cpp
  algorithms/test_exp_sym_v.cpp
  algorithms/test_norm_estimate.cpp
  algorithms/time_evolution/test_time_evolution.cpp
  algorithms/time_evolution/test_pade.cpp

  states/test_random_state.cpp
  states/test_product_state.cpp
)
target_link_libraries(tests PRIVATE hydra)
add_test(NAME HydraTest COMMAND tests)
