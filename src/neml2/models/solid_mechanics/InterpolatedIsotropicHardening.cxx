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

#include "neml2/models/solid_mechanics/InterpolatedIsotropicHardening.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/indexing.h"
// #include "neml2/tensors/functions/cat.h"

namespace neml2
{
register_NEML2_object(InterpolatedIsotropicHardening);

OptionSet
InterpolatedIsotropicHardening::expected_options()
{
  OptionSet options = IsotropicHardening::expected_options();
  options.doc() = "Voce isotropic hardening model, \\f$ h = R \\left[ 1 - \\exp(-d \\varepsilon_p) "
                  "\\right] \\f$, where \\f$ R \\f$ is the isotropic hardening upon saturation, "
                  "and \\f$ d \\f$ is the hardening rate.";

  options.set<bool>("define_second_derivatives") = true;

  options.set_input("temp_argument");
  options.set("temp_argument").doc() = "Temperature argument used to query the interpolant";

  options.set<std::vector<VariableName>>("stresses");
  options.set("stresses").doc() = "Stress list";

  options.set_parameter<std::vector<TensorName<Scalar>>>("hf_temperatures");
  options.set("hf_temperatures").doc() = "List of temperatures for each stress";

  options.set_output("isotropic_hardening") = VariableName(STATE, "internal", "k");
  options.set("isotropic_hardening").doc() = "Isotropic hardening";

  return options;
}

InterpolatedIsotropicHardening::InterpolatedIsotropicHardening(const OptionSet & options)
  : Model(options),
    _temp_var(declare_input_variable<Scalar>("temp_argument")),
    _h(declare_output_variable<Scalar>("isotropic_hardening"))
{
  for (const auto & fv : options.get<std::vector<VariableName>>("stresses"))
    _stress.push_back(&declare_input_variable<Scalar>(fv));

  const auto hf_temp_refs = options.get<std::vector<TensorName<Scalar>>>("hf_temperatures");
  _hf_temps.resize(_stress.size());
  for (std::size_t i = 0; i < _stress.size(); i++)
  {
    _hf_temps[i] = &declare_buffer<Scalar>("hf_temp_" + std::to_string(i), hf_temp_refs[i]);
  }
}

void
InterpolatedIsotropicHardening::set_value(bool out, bool dout_din, bool d2out_din2)
{
  auto slope = Scalar::zeros(_stress.size() - 1);
  auto X0 = Scalar::zeros(_stress.size() - 1);
  auto X1 = Scalar::zeros(_stress.size() - 1);
  auto Y0 = Scalar::zeros(_stress.size() - 1);

  auto stress_dev_lo = Scalar::zeros(_stress.size() - 1);
  auto stress_dev_hi = Scalar::zeros(_stress.size() - 1);
  auto stress_dev = Scalar::zeros(_stress.size());

  for (std::size_t i = 0; i < _stress.size() - 1; i++)
  {
    slope[i] = (*_stress[i + 1] - *_stress[i]) / (*_hf_temps[i + 1] - *_hf_temps[i]);
    X0[i] = *_hf_temps[i];
    X1[i] = *_hf_temps[i + 1];
    Y0[i] = *_stress[i];

    stress_dev_lo[i] = (*_hf_temps[i + 1] - _temp_var) / (*_hf_temps[i + 1] - *_hf_temps[i]);
    stress_dev_hi[i] = (_temp_var - *_hf_temps[i]) / (*_hf_temps[i + 1] - *_hf_temps[i]);
  }

  // cannot use gt,le as the finite differencing uses the higher point's gradient at abrupt changes
  // in gradient, e.g. changing from 1 to 0, it will use the 0
  const auto loc = Scalar(at::logical_and(at::ge(_temp_var, X0), at::lt(_temp_var, X1)));
  const auto si = Scalar(slope.index({loc}));

  // partial derivatives for stresses
  const auto stress_dev_lo_i = Scalar(stress_dev_lo.index({loc}));
  const auto stress_dev_hi_i = Scalar(stress_dev_hi.index({loc}));
  stress_dev.batch_index({indexing::Slice(indexing::None, -1)}).index_put_({loc}, stress_dev_lo_i);
  stress_dev.batch_index({indexing::Slice(1, indexing::None)}).index_put_({loc}, stress_dev_hi_i);

  if (out)
  {
    const auto X0i = Scalar(X0.index({loc}));
    const auto Y0i = Scalar(Y0.index({loc}));

    _h = Y0i + si * (_temp_var - X0i);
  }

  if (dout_din)
  {
    if (_temp_var.is_dependent())
      _h.d(_temp_var) = si;

    for (std::size_t i = 0; i < _stress.size(); i++)
    {
      if (_stress[i]->is_dependent())
      {
        _h.d(*_stress[i]) = Scalar(stress_dev[i]);
      }
    }
  }

  if (d2out_din2)
  {
  }
}
} // namespace neml2