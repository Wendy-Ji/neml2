// Copyright 2024, UChicago Argonne, LLC
// All Rights Reserved
// Software Name: NEML2 -- the New Engineering material Model Library, version 2
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE (MIT)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "neml2/models/LinearInterpolationOnVariable.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/indexing.h"

namespace neml2
{

register_NEML2_object(LinearInterpolationOnVariable);

OptionSet
LinearInterpolationOnVariable::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Voce isotropic hardening model, \\f$ h = R \\left[ 1 - \\exp(-d \\varepsilon_p) "
                  "\\right] \\f$, where \\f$ R \\f$ is the isotropic hardening upon saturation, "
                  "and \\f$ d \\f$ is the hardening rate.";

  options.set<bool>("define_second_derivatives") = true;

  options.set_input("argument");
  options.set("argument").doc() = "Argument used to query the interpolant";

  options.set_parameter<std::vector<TensorName<Scalar>>>("abscissa_vector");
  options.set("abscissa_vector").doc() =
      "Vector of scalars defining the abscissa values of the interpolant";

  options.set<std::vector<VariableName>>("ordinate_vector");
  options.set("ordinate_vector").doc() =
      "Vector of variables defining the ordinate values of the interpolant";

  options.set_output("output");
  options.set("output").doc() = "output";

  return options;
}

LinearInterpolationOnVariable::LinearInterpolationOnVariable(const OptionSet & options)
  : Model(options),
    _x(declare_input_variable<Scalar>("argument")),
    _output(declare_output_variable<Scalar>("output"))
{
  for (const auto & y : options.get<std::vector<VariableName>>("ordinate_vector"))
    _Y.push_back(&declare_input_variable<Scalar>(y));

  const auto X_refs = options.get<std::vector<TensorName<Scalar>>>("abscissa_vector");
  _X.resize(_Y.size());
  for (std::size_t i = 0; i < _Y.size(); i++)
  {
    _X[i] = &declare_buffer<Scalar>("X_" + std::to_string(i), X_refs[i]);
  }
}

void
LinearInterpolationOnVariable::set_value(bool out, bool dout_din, bool d2out_din2)
{
  auto slope = Scalar::zeros(_Y.size() - 1);
  auto X0 = Scalar::zeros(_Y.size() - 1);
  auto X1 = Scalar::zeros(_Y.size() - 1);
  auto Y0 = Scalar::zeros(_Y.size() - 1);

  // partial derivatives for ordinate variables
  auto Y_dev_lo = Scalar::zeros(_Y.size() - 1);
  auto Y_dev_hi = Scalar::zeros(_Y.size() - 1);
  auto Y_dev = Scalar::zeros(_Y.size());

  for (std::size_t i = 0; i < _Y.size() - 1; i++)
  {
    slope[i] = (*_Y[i + 1] - *_Y[i]) / (*_X[i + 1] - *_X[i]);
    X0[i] = *_X[i];
    X1[i] = *_X[i + 1];
    Y0[i] = *_Y[i];

    Y_dev_lo[i] = (*_X[i + 1] - _x) / (*_X[i + 1] - *_X[i]);
    Y_dev_hi[i] = 1 - Y_dev_lo[i];
  }

  const auto loc = Scalar(at::logical_and(at::ge(_x, X0), at::lt(_x, X1)));
  const auto si = Scalar(slope.index({loc}));

  Y_dev.batch_index({indexing::Slice(indexing::None, -1)})
      .index_put_({loc}, Scalar(Y_dev_lo.index({loc})));
  Y_dev.batch_index({indexing::Slice(1, indexing::None)})
      .index_put_({loc}, Scalar(Y_dev_hi.index({loc})));

  if (out)
  {
    const auto X0i = Scalar(X0.index({loc}));
    const auto Y0i = Scalar(Y0.index({loc}));

    _output = Y0i + si * (_x - X0i);
  }

  if (dout_din)
  {
    if (_x.is_dependent())
      _output.d(_x) = si;

    for (std::size_t i = 0; i < _Y.size(); i++)
    {
      if (_Y[i]->is_dependent())
      {
        _output.d(*_Y[i]) = Scalar(Y_dev[i]);
      }
    }
  }

  if (d2out_din2)
  {
  }
}
} // namespace neml2