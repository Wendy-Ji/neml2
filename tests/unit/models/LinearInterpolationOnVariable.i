[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/X state/A state/B state/C'
    input_Scalar_values = '250 60 70 80'
    output_Scalar_names = 'state/outsub/C'
    output_Scalar_values = '75'
  []
[]

[Models]
  [model]
    type = LinearInterpolationOnVariable
    argument = 'state/X'
    abscissa = 'X'
    ordinate_vector = 'state/A state/B state/C'
    output = 'state/outsub/C'
  []
[]

[Tensors]
  [X]
    type = Scalar
    values = '100 200 300'
    batch_shape = '(3)'
  []
[]