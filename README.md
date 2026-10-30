# Native proFit

A collection of native Fortran study primitives and model building blocks.

This independent package shares the study goals of
[Python proFit](https://github.com/redmod-team/profit). It has no Python runtime
dependency.

The initial API scores Gaussian predictions, summarizes independent study
replicates, writes study tables with consumer-supplied provenance, and evaluates
Gaussian expected improvement. The study API and model building blocks do not
themselves provide application-specific training objectives or a full study
runner, or certify a global optimum. Gaussian prediction and inference belong in
model and probability libraries; application equations belong with their study.

Build with a Fortran 2008 compiler and [fpm](https://fpm.fortran-lang.org/):

```sh
fpm test
fpm run --example score_study
```

All numerical arguments use `real(real64)` from `iso_fortran_env`. A `status`
of `profit_success` (zero) means success. Other named values in `profit_status`
distinguish shape, input, numerical range and file errors. Numerical outputs
are zeroed on failure. Large independent loops use OpenMP; the number of threads
is controlled by `OMP_NUM_THREADS`.

```fortran
use profit_scores, only: gaussian_scores
call gaussian_scores(y, mean, stddev, z, rmse, nll, coverage, status)
```

`y`, `mean` and `stddev` are equally sized vectors in the same physical units;
`stddev` must be strictly positive. `rmse` is root mean squared error in those
units. `nll` is the mean negative log marginal Gaussian density, using natural logarithms
and the supplied units. It does not evaluate a correlated joint likelihood. Rescaling all three vectors by a positive factor `a`
multiplies `rmse` by `a` and adds `log(a)` to `nll`. `coverage` is the fraction
with `abs(y-mean) <= z*stddev`. The caller selects `z`: for example,
`1.959963984540054` denotes a central nominal 95% Gaussian interval. Empirical
coverage alone does not establish calibration, especially on correlated data.

```fortran
use profit_replicates, only: replicate_summary
call replicate_summary(values, critical, mean, se, lower, upper, status)
```

`values` contains at least two **independent replicate estimates of the same
quantity**, such as one held-out RMSE per independent simulation. `se` is the
sample standard deviation divided by the square root of the replicate count.
`lower` and `upper` are `mean +/- critical*se`, describing uncertainty in the
replicate mean. They are distinct from predictive intervals for an observation.
The caller supplies a positive critical value: a Student-t quantile with one
fewer degree of freedom than replicates gives the usual interval for independent
normally distributed estimates; a normal quantile requires an appropriate
large-sample justification. The routine cannot establish independence or those
distributional assumptions from a vector.

```fortran
use profit_design, only: expected_improvement
call expected_improvement(mean, stddev, incumbent, minimize, ei, status, &
                          dmean, dstddev)
```

The input vectors hold marginal Gaussian means and nonnegative standard
deviations at candidate locations. `minimize=.true.` measures improvement below
`incumbent`; `.false.` measures improvement above it. `ei` has the objective's
units. Optional output vectors `dmean` and `dstddev` give its local derivatives
for positive standard deviations. At zero standard deviation the value is the
positive part of deterministic improvement. At a zero-variance tie the mean
derivative is undefined: the API returns zero by convention and returns the
right derivative with respect to standard deviation. Beyond 40 standardized
deviations the limiting deterministic value and derivatives are used. The
routine evaluates an acquisition function; candidate optimization and noisy,
batch or constrained acquisitions are separate concerns.

```fortran
use profit_output, only: write_study_csv
call write_study_csv(path, keys, provenance, columns, values, status)
```

`values(column,row)` is a finite numerical table. The output is rectangular CSV
(comma-separated values), with three leading columns `record_type`,
`metadata_key`, `metadata_value`. Each `metadata` record contains a key and
value, followed by empty numerical fields; each `data` record has empty metadata
fields followed by numerical values. String fields are quoted and embedded
quotes are doubled. The consumer supplies units in column names and provenance
such as data origin, revision, seed, interval definition and partition. No
provenance is inferred. A successful call replaces the destination; a file
error may leave a partial file.

Pointwise numerical leaves are committed under `src/generated/`. Ordinary
builds need no symbolic runtime. Regeneration uses the separate pinned
[FortSym](https://github.com/lazy-fortran/fortsym) development package:

```sh
cd gen
fpm run -- ../src/generated
git diff -- ../src/generated
```

The generator checks the two expected-improvement derivative identities
symbolically before emitting their stable `erfc` form. Tests independently
check those derivatives by central differences and the acquisition value by
Gaussian quadrature integration.

## Tree and smooth models

`profit_models_tree` fits a shallow weighted regression tree using histogram
splits. Inputs have shape `(nfeature, nrow)`. Fit bin thresholds once using an
explicit fitting mask, then reuse them when fitting successive residuals:

```fortran
use profit_models_tree, only: tree_options_t, tree_bins_t, regression_tree_t
type(tree_options_t) :: options
type(tree_bins_t) :: bins
type(regression_tree_t) :: tree

call bins%fit(features, fit_mask, options, status)
call bins%transform(features, binned, status)
call tree%fit(bins, binned, response, weights, fit_mask, options, status)
call tree%predict(features, prediction, status)
```

Check `status` after every call. `binned` is a caller-allocated integer array
with the same shape as `features`; response, weights and masks have one entry
per row. The defaults are `max_depth=3`, `min_leaf=20`, `max_bins=32` and
`leaf_l2=0`. These set maximum depth, minimum fitting rows per leaf, bin count
and a nonnegative ridge penalty that shrinks leaf values toward zero.
`fit_mask=.true.` selects fitting rows. Weights specify weighted squared-error
fitting: fitting weights must be finite and nonnegative, with at least one
positive weight. Zero-weight fitting rows still count toward `min_leaf`. Thresholds depend only
on fitting features, and fitting ignores non-fitting responses and weights.
An optional `valid` mask on `transform` excludes unavailable input rows.
Prediction requires finite inputs. Ties use a fixed feature and threshold
order. The tree exposes leaf indices and editable leaf values for applications
that differentiate a loss with respect to fixed-tree parameters. Hard splits
are discontinuous in their inputs, so these are not continuous input
derivatives. Versioned `save` and `load` methods support native-format
round trips; this is not a portable interchange format.

`profit_models_network` evaluates an affine model or a single hidden layer
with smooth hyperbolic-tangent activations and a scalar output:

```fortran
use profit_models_network, only: network_parameter_count, network_predict

nparameter = network_parameter_count(nfeature, nhidden)
call network_predict(nfeature, nhidden, theta, features, prediction, status)
```

Allocate `theta(nparameter)`, `features(nfeature,nrow)` and
`prediction(nrow)`. For `nhidden > 0`, `theta` contains the column-major input-weight matrix
`(nfeature, nhidden)`, hidden biases, output weights and output bias, in that
order. With `nhidden = 0`, it contains the affine input weights followed by
the output bias. `network_jvp` applies parameter and input directional
derivatives; `network_vjp` applies their transposes to output weights. These
products use the following interfaces; the transpose product sums parameter
contributions across rows:

```fortran
call network_jvp(nfeature, nhidden, theta, features, theta_direction, &
                 feature_direction, dvalues, status)
call network_vjp(nfeature, nhidden, theta, features, output_adjoint, &
                 theta_adjoint, feature_adjoint, status)
```

Derivative arrays match their corresponding input or output shapes. These
products permit local covariance propagation when the caller supplies a
justified joint input or parameter covariance. They do not supply a posterior
or establish interval calibration. `network_positive_map` supplies stable
softplus values and derivatives, including the floating-point zero limit.

The model modules provide fitting or evaluation building blocks. Applications
own their training objective, temporal dynamics, model selection and predictive
uncertainty. OpenMP parallelizes independent feature histograms and batched
evaluation; fitting-bin construction is separate preprocessing. The native
tests check tree splits against exhaustive small-data enumeration and neural
derivative products against finite differences. The examples report measured
costs for their stated workloads.
