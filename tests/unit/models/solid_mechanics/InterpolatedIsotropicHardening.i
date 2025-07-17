[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/T state/internal/ep1 state/internal/ep2 state/internal/ep3'
    input_Scalar_values = '130 0 1000 2000'
    output_Scalar_names = 'state/internal/k'
    output_Scalar_values = '300'
  []
[]

[Models]
  [model0]
    type = InterpolatedIsotropicHardening
    temp_argument = 'state/T'
    hf_temperatures = '100 200 300'
    stresses = 'state/internal/ep1 state/internal/ep2 state/internal/ep3'
  []
  [model]
    type = ComposedModel
    models = 'model0'
  []
[]