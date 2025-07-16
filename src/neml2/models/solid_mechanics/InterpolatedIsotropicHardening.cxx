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

  // auto stress_dev_lo = Scalar::zeros(_stress.size() - 1);
  // auto stress_dev = Scalar::zeros(_stress.size());

  for (std::size_t i = 0; i < _stress.size() - 1; i++)
  {
    slope[i] = (*_stress[i + 1] - *_stress[i]) / (*_hf_temps[i + 1] - *_hf_temps[i]);
    X0[i] = *_hf_temps[i];
    X1[i] = *_hf_temps[i + 1];
    Y0[i] = *_stress[i];

    // stress_dev_lo = (*_hf_temps[i + 1] - _temp_var) / (*_hf_temps[i + 1] - *_hf_temps[i]);
  }

  const auto loc = Scalar(at::logical_and(at::gt(_temp_var, X0), at::le(_temp_var, X1)));

  const auto si = Scalar(slope.index({loc}));

  // auto ind_lo = Scalar::linspace();

  // stress_dev.index({loc}) = stress_dev_lo.index({loc});

  // works for two values
  // const auto ind_lo = 0;
  // const auto ind_hi = 1;
  // const auto frac = (_temp_var - *_hf_temps[ind_lo]) / (*_hf_temps[ind_hi] - *_hf_temps[ind_lo]);

  if (out)
  {
    const auto X0i = Scalar(X0.index({loc}));
    const auto Y0i = Scalar(Y0.index({loc}));
    _h = Y0i + si * (_temp_var - X0i);
    // works for two values
    // _h = (1 - frac) * (*_stress[ind_lo]) + frac * (*_stress[ind_hi]);
  }

  if (dout_din)
  {
    if (_temp_var.is_dependent())
      _h.d(_temp_var) = si;

    // for (std::size_t i = 0; i < _stress.size(); i++)
    // {
    //   if (_stress[i]->is_dependent())
    //   {
    //     _h.d(*_stress[i]) = Scalar(stress_dev[i]);
    //   }
    // }

    // the derivatives
    if (_stress[0]->is_dependent())
      _h.d(*_stress[0]) = (*_hf_temps[1] - _temp_var) / (*_hf_temps[1] - *_hf_temps[0]);

    if (_stress[1]->is_dependent())
      _h.d(*_stress[1]) = (_temp_var - *_hf_temps[0]) / (*_hf_temps[1] - *_hf_temps[0]);

    // works for two values
    // for (std::size_t i = 0; i < _stress.size(); i++)
    // {
    //   if (_stress[i]->is_dependent())
    //   {
    //     if (i == ind_lo)
    //       _h.d(*_stress[i]) = 1 - frac;
    //     else if (i == ind_hi)
    //       _h.d(*_stress[i]) = frac;
    //     else
    //       _h.d(*_stress[i]) = neml2::Scalar::full(0.0);
    //   }
    // }

    // if (_temp_var.is_dependent())
    //   _h.d(_temp_var) =
    //       (*_stress[ind_hi] - *_stress[ind_lo]) / (*_hf_temps[ind_hi] - *_hf_temps[ind_lo]);
  }

  if (d2out_din2)
  {
  }
}
} // namespace neml2