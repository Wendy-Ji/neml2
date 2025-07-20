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
#include "neml2/models/LinearInterpolation.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/indexing.h"
#include "neml2/tensors/functions/stack.h"
#include "neml2/tensors/functions/diff.h"

namespace neml2
{

register_NEML2_object(LinearInterpolationOnVariable);

OptionSet
LinearInterpolationOnVariable::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "This object performs a linear interpolation, with scalars as the abscissa and "
                  "variables as the ordinates.";

  options.set<bool>("define_second_derivatives") = true;

  options.set_input("argument");
  options.set("argument").doc() = "Argument used to query the interpolant";

  options.set<TensorName<Scalar>>("abscissa");
  options.set("abscissa").doc() = "Scalar defining the abscissa values of the interpolant";

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
    _X(declare_buffer<Scalar>("X", "abscissa")),
    _output(declare_output_variable<Scalar>("output"))
{
  for (const auto & y : options.get<std::vector<VariableName>>("ordinate_vector"))
    _Y.push_back(&declare_input_variable<Scalar>(y));
}

void
LinearInterpolationOnVariable::set_value(bool out, bool dout_din, bool d2out_din2)
{
  std::vector<Scalar> ytens(_Y.size());
  for (std::size_t i = 0; i < _Y.size(); i++)
    ytens[i] = _Y[i]->value();
  const auto Y = batch_stack(ytens, -1);

  const auto slope = diff(Y) / diff(_X);
  const auto X0 = _X.batch_index({indexing::Ellipsis, indexing::Slice(indexing::None, -1)});
  const auto X1 = _X.batch_index({indexing::Ellipsis, indexing::Slice(1, indexing::None)});
  const auto Y0 = Y.batch_index({indexing::Ellipsis, indexing::Slice(indexing::None, -1)});

  const auto x = Scalar(_x);
  const auto loc =
      Scalar(at::logical_and(at::gt(x.batch_unsqueeze(-1), X0), at::le(x.batch_unsqueeze(-1), X1)));
  const auto si = LinearInterpolation<Scalar>::mask(slope, loc);

  auto Y_dev = Scalar::zeros(_Y.size());
  const auto Ydev0 = (X1 - x.batch_unsqueeze(-1)) / diff(_X);
  std::cout << "x.batch_unsqueeze(-1): " << x.batch_unsqueeze(-1).batch_sizes() << " "
            << x.batch_unsqueeze(-1).base_sizes() << std::endl;
  std::cout << "X0: " << X0.batch_sizes() << " " << X0.base_sizes() << std::endl;
  std::cout << "X1: " << X1.batch_sizes() << " " << X1.base_sizes() << std::endl;
  std::cout << "Ydev: " << Y_dev.batch_sizes() << " " << Y_dev.base_sizes() << std::endl;
  std::cout << "Ydev0: " << Ydev0.batch_sizes() << " " << Ydev0.base_sizes() << std::endl;
  std::cout << "si: " << si.batch_sizes() << " " << si.base_sizes() << std::endl;
  Y_dev.batch_unsqueeze(-1)
      .batch_index({indexing::Slice(indexing::None, -1)})
      .index_put_({loc}, LinearInterpolation<Scalar>::mask(Ydev0, loc));
  std::cout << "wow" << std::endl;
  Y_dev.batch_unsqueeze(-1)
      .batch_index({indexing::Slice(1, indexing::None)})
      .index_put_({loc}, 1 - LinearInterpolation<Scalar>::mask(Ydev0, loc));
  std::cout << "wow2" << std::endl;

  if (out)
  {
    const auto X0i = LinearInterpolation<Scalar>::mask(X0, loc);
    const auto Y0i = LinearInterpolation<Scalar>::mask(Y0, loc);

    _output = Y0i + si * (x - X0i);
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