[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/T state/A state/B state/C state/D state/E'
    input_Scalar_values = '150 60 60 70 80 80'
    output_Scalar_names = 'state/outsub/C'
    output_Scalar_values = '65'
  []
[]

[Tensors]
  [X]
    type = Scalar
    values = '0 100 200 300 400'
  []
[]

[Models]
  [model]
    type = LinearInterpolationOnVariable
    argument = 'state/T'
    abscissa = 'X'
    ordinate_vector = 'state/A state/B state/C state/D state/E'
    output = 'state/outsub/C'
  []
[]