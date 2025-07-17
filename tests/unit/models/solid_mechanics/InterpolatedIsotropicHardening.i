[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/T state/ep1 state/ep2'
    input_Scalar_values = '140 60 70'
    output_Scalar_names = 'state/internal/k'
    output_Scalar_values = '64'
  []
[]

[Models]
  [model]
    type = InterpolatedIsotropicHardening
    temp_argument = 'state/T'
    hf_temperatures = '100 200'
    stresses = 'state/ep1 state/ep2'
  []
[]