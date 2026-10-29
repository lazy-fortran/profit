# Native proFit

An initial collection of native Fortran study primitives.

This independent package shares the study goals of
[Python proFit](https://github.com/redmod-team/profit). It has no Python runtime
dependency.

The initial API scores Gaussian predictions, summarizes independent study
replicates, writes study tables with consumer-supplied provenance, and evaluates
Gaussian expected improvement. It does not fit a surrogate or certify a global
optimum. Gaussian prediction and inference belong in model and probability
libraries; application equations belong with their study.

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
