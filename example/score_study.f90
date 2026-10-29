program score_study
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use profit_scores, only: gaussian_scores
    use profit_replicates, only: replicate_summary
    use profit_output, only: write_study_csv
    implicit none
    character(24), parameter :: keys(3) = &
        [character(24) :: "source", "interval", "replicate_interval"]
    character(80), parameter :: metadata(3) = [character(80) :: &
        "synthetic known-answer fixture", &
        "central nominal 95% Gaussian predictive interval", &
        "Student-t 95%, df=2; illustrative independent normal estimates"]
    character(24), parameter :: columns(3) = &
        [character(24) :: "RMSE [K]", "NLL [log K]", "coverage [1]"]
    real(dp) :: rmse, nll, coverage, average, se, lower, upper, table(3, 1)
    integer :: status
    call gaussian_scores([0.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], &
        [1.0_dp, 1.0_dp], 1.959963984540054_dp, rmse, nll, coverage, status)
    if (status /= 0) error stop "prediction scoring failed"
    call replicate_summary([1.0_dp, 2.0_dp, 3.0_dp], 4.302652729696142_dp, &
        average, se, lower, upper, status)
    if (status /= 0) error stop "replicate summary failed"
    table(1, 1) = rmse
    table(2, 1) = nll
    table(3, 1) = coverage
    call write_study_csv("synthetic-study.csv", keys, metadata, columns, table, status)
    if (status /= 0) error stop "study output failed"
    print '(a,f8.5)', "Independent replicate mean [K]: ", average
    print '(a,f8.5)', "Standard error of that mean [K]: ", se
    print '(a,2f10.5)', "Illustrative confidence interval [K]: ", lower, upper
end program score_study
