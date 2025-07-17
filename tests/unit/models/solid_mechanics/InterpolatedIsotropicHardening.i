[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/T state/dummy_lo state/ep1 state/ep2 state/ep3 state/dummy_hi'
    input_Scalar_values = '150 60 60 70 80 80'
    output_Scalar_names = 'state/internal/k'
    output_Scalar_values = '60'
  []
[]

[Models]
  [model]
    type = InterpolatedIsotropicHardening
    temp_argument = 'state/T'
    hf_temperatures = '0 100 200 300 400'
    stresses = 'state/dummy_lo state/ep1 state/ep2 state/ep3 state/dummy_hi'
  []
[]