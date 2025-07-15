[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/T state/internal/ep1 state/internal/ep2'
    input_Scalar_values = '250 600 1000'
    output_Scalar_names = 'state/internal/k'
    output_Scalar_values = '700'
  []
[]

[Models]
  [model0]
    type = InterpolatedIsotropicHardening
    temp_argument = 'state/T'
    hf_temperatures = '200 400'
    stresses = 'state/internal/ep1 state/internal/ep2'
  []
  [model]
    type = ComposedModel
    models = 'model0'
  []
[]