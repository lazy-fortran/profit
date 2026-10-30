module profit_models_tree
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use profit_status, only: profit_success, profit_invalid_shape, &
        profit_invalid_input, profit_range_error, profit_io_error
    use profit_generated_tree_leaf, only: tree_leaf
    implicit none
    private

    type, public :: tree_options_t
        integer :: max_depth = 3
        integer :: min_leaf = 20
        integer :: max_bins = 32
        real(dp) :: leaf_l2 = 0
    end type tree_options_t

    type, public :: tree_bins_t
        private
        integer :: n_features = 0
        integer :: max_bins = 0
        integer, allocatable :: n_cuts(:)
        real(dp), allocatable :: cuts(:, :)
    contains
        procedure :: fit => bins_fit
        procedure :: transform => bins_transform
    end type tree_bins_t

    type, public :: regression_tree_t
        private
        integer :: n_features = 0
        integer :: node_count = 0
        integer :: leaf_count = 0
        integer, allocatable :: feature(:), left(:), right(:), leaf(:)
        real(dp), allocatable :: threshold(:), leaf_value(:)
    contains
        procedure :: fit => tree_fit
        procedure :: predict => tree_predict
        procedure :: leaf_indices => tree_leaf_indices
        procedure :: get_leaf_values => tree_get_leaf_values
        procedure :: set_leaf_values => tree_set_leaf_values
        procedure :: n_leaves => tree_n_leaves
        procedure :: save => tree_save
        procedure :: load => tree_load
    end type regression_tree_t

contains

    subroutine bins_fit(self, features, fit_mask, options, status)
        class(tree_bins_t), intent(out) :: self
        real(dp), intent(in) :: features(:, :)
        logical, intent(in) :: fit_mask(:)
        type(tree_options_t), intent(in) :: options
        integer, intent(out) :: status
        real(dp), allocatable :: values(:), scratch(:)
        integer :: p, n, nf, f, i, j, k, cut_index

        status = profit_invalid_shape
        p = size(features, 1)
        n = size(features, 2)
        if (p == 0 .or. n == 0) return
        if (size(fit_mask) /= n) return
        status = profit_invalid_input
        if (.not. valid_options(options)) return
        nf = count(fit_mask)
        if (nf == 0) return
        do i = 1, n
            if (.not. fit_mask(i)) cycle
            if (.not. all(ieee_is_finite(features(:, i)))) return
        end do
        self%n_features = p
        self%max_bins = options%max_bins
        allocate(self%n_cuts(p), self%cuts(options%max_bins - 1, p))
        allocate(values(nf), scratch(nf))
        self%n_cuts = 0
        self%cuts = 0
        do f = 1, p
            k = 0
            do i = 1, n
                if (.not. fit_mask(i)) cycle
                k = k + 1
                values(k) = features(f, i)
            end do
            call merge_sort(values, scratch)
            do j = 1, options%max_bins - 1
                cut_index = int(real(j, dp)*real(nf, dp)/real(options%max_bins, dp))
                cut_index = max(1, min(nf - 1, cut_index))
                if (nf == 1) exit
                if (values(cut_index) >= values(nf)) cycle
                k = self%n_cuts(f)
                if (k > 0) then
                    if (values(cut_index) <= self%cuts(k, f)) cycle
                end if
                self%n_cuts(f) = k + 1
                self%cuts(k + 1, f) = values(cut_index)
            end do
        end do
        status = profit_success
    end subroutine bins_fit

    subroutine bins_transform(self, features, binned, status, valid)
        class(tree_bins_t), intent(in) :: self
        real(dp), intent(in) :: features(:, :)
        integer, intent(out) :: binned(:, :)
        integer, intent(out) :: status
        logical, optional, intent(in) :: valid(:)
        integer :: i, f, lo, hi, mid
        logical :: finite_inputs

        binned = 0
        status = profit_invalid_shape
        if (self%n_features == 0) return
        if (size(features, 1) /= self%n_features) return
        if (any(shape(binned) /= shape(features))) return
        if (present(valid)) then
            if (size(valid) /= size(features, 2)) return
        end if
        status = profit_invalid_input
        finite_inputs = .true.
        !$omp parallel do default(shared) private(i) reduction(.and.:finite_inputs)
        do i = 1, size(features, 2)
            if (present(valid)) then
                if (.not. valid(i)) cycle
            end if
            if (.not. all(ieee_is_finite(features(:, i)))) finite_inputs = .false.
        end do
        !$omp end parallel do
        if (.not. finite_inputs) return
        !$omp parallel do default(shared) private(i, f, lo, hi, mid)
        do i = 1, size(features, 2)
            if (present(valid)) then
                if (.not. valid(i)) cycle
            end if
            do f = 1, self%n_features
                lo = 1
                hi = self%n_cuts(f) + 1
                do while (lo < hi)
                    mid = lo + (hi - lo)/2
                    if (features(f, i) <= self%cuts(mid, f)) then
                        hi = mid
                    else
                        lo = mid + 1
                    end if
                end do
                binned(f, i) = lo
            end do
        end do
        !$omp end parallel do
        status = profit_success
    end subroutine bins_transform

    subroutine tree_fit(self, bins, binned, response, weights, fit_mask, options, status)
        class(regression_tree_t), intent(out) :: self
        type(tree_bins_t), intent(in) :: bins
        integer, intent(in) :: binned(:, :)
        real(dp), intent(in) :: response(:), weights(:)
        logical, intent(in) :: fit_mask(:)
        type(tree_options_t), intent(in) :: options
        integer, intent(out) :: status
        integer, allocatable :: rows(:), work(:), hist_count(:, :)
        real(dp), allocatable :: hist_weight(:, :), hist_sum(:, :)
        integer :: p, n, nf, i, j, max_nodes, root
        logical :: has_weight, valid_inputs

        status = profit_invalid_shape
        p = size(binned, 1)
        n = size(binned, 2)
        if (p == 0 .or. n == 0) return
        if (bins%n_features /= p) return
        if (size(response) /= n .or. size(weights) /= n) return
        if (size(fit_mask) /= n) return
        status = profit_invalid_input
        if (.not. valid_options(options)) return
        nf = count(fit_mask)
        if (nf == 0) return
        has_weight = .false.
        valid_inputs = .true.
        !$omp parallel do default(shared) private(i) &
        !$omp reduction(.and.:valid_inputs) reduction(.or.:has_weight)
        do i = 1, n
            if (.not. fit_mask(i)) cycle
            if (.not. ieee_is_finite(response(i))) valid_inputs = .false.
            if (.not. ieee_is_finite(weights(i))) valid_inputs = .false.
            if (weights(i) < 0) valid_inputs = .false.
            if (weights(i) > 0) has_weight = .true.
            if (any(binned(:, i) < 1)) valid_inputs = .false.
            if (any(binned(:, i) > bins%n_cuts + 1)) valid_inputs = .false.
        end do
        !$omp end parallel do
        if (.not. valid_inputs) return
        if (.not. has_weight) return
        self%n_features = p
        max_nodes = 2**(options%max_depth + 1) - 1
        allocate(self%feature(max_nodes), self%left(max_nodes), &
            self%right(max_nodes), self%leaf(max_nodes), &
            self%threshold(max_nodes), self%leaf_value(max_nodes))
        self%feature = 0
        self%left = 0
        self%right = 0
        self%leaf = 0
        self%threshold = 0
        self%leaf_value = 0
        allocate(rows(nf), work(nf))
        allocate(hist_count(bins%max_bins, p), hist_weight(bins%max_bins, p), &
            hist_sum(bins%max_bins, p))
        j = 0
        do i = 1, n
            if (.not. fit_mask(i)) cycle
            j = j + 1
            rows(j) = i
        end do
        status = profit_success
        call build(1, nf, 0, root)
        if (status /= profit_success) then
            self%node_count = 0
            self%leaf_count = 0
        end if
    contains
        recursive subroutine build(first, last, depth, node)
            integer, intent(in) :: first, last, depth
            integer, intent(out) :: node
            integer :: f, j, row, b, best_feature, best_cut, left_count
            integer :: nl, nr, split, lc, rc, out_left, out_right
            real(dp) :: total_weight, total_sum, parent_score, unused
            real(dp) :: wl, wr, yl, yr, sl, sr, gain, best_gain

            self%node_count = self%node_count + 1
            node = self%node_count
            total_weight = 0
            total_sum = 0
            do j = first, last
                row = rows(j)
                total_weight = total_weight + weights(row)
                total_sum = total_sum + weights(row)*response(row)
            end do
            status = profit_range_error
            if (.not. ieee_is_finite(total_weight)) return
            if (.not. ieee_is_finite(total_sum)) return
            status = profit_success
            best_feature = 0
            best_cut = 0
            if (total_weight > 0) then
                call tree_leaf(total_sum, total_weight, options%leaf_l2, unused, parent_score)
                if (.not. ieee_is_finite(parent_score)) then
                    status = profit_range_error
                    return
                end if
                if (depth < options%max_depth .and. &
                    options%min_leaf <= (last - first + 1)/2) then
                    hist_count = 0
                    hist_weight = 0
                    hist_sum = 0
                    !$omp parallel do default(shared) private(f, j, row, b)
                    do f = 1, p
                        do j = first, last
                            row = rows(j)
                            b = binned(f, row)
                            hist_count(b, f) = hist_count(b, f) + 1
                            hist_weight(b, f) = hist_weight(b, f) + weights(row)
                            hist_sum(b, f) = hist_sum(b, f) + weights(row)*response(row)
                        end do
                    end do
                    !$omp end parallel do
                    best_gain = 0
                    do f = 1, p
                        wl = 0
                        yl = 0
                        nl = 0
                        do b = 1, bins%n_cuts(f)
                            wl = wl + hist_weight(b, f)
                            yl = yl + hist_sum(b, f)
                            nl = nl + hist_count(b, f)
                            nr = last - first + 1 - nl
                            if (nl < options%min_leaf .or. nr < options%min_leaf) cycle
                            wr = total_weight - wl
                            yr = total_sum - yl
                            if (wl <= 0 .or. wr <= 0) cycle
                            call tree_leaf(yl, wl, options%leaf_l2, unused, sl)
                            call tree_leaf(yr, wr, options%leaf_l2, unused, sr)
                            gain = sl + sr - parent_score
                            if (.not. ieee_is_finite(gain)) then
                                status = profit_range_error
                                return
                            end if
                            if (gain > best_gain) then
                                best_gain = gain
                                best_feature = f
                                best_cut = b
                            end if
                        end do
                    end do
                end if
            end if
            if (best_feature == 0) then
                self%leaf_count = self%leaf_count + 1
                self%leaf(node) = self%leaf_count
                if (total_weight > 0) then
                    call tree_leaf(total_sum, total_weight, options%leaf_l2, &
                        self%leaf_value(self%leaf_count), unused)
                end if
                return
            end if
            self%feature(node) = best_feature
            self%threshold(node) = bins%cuts(best_cut, best_feature)
            left_count = 0
            do j = first, last
                if (binned(best_feature, rows(j)) <= best_cut) left_count = left_count + 1
            end do
            split = first + left_count - 1
            lc = first
            rc = split + 1
            do j = first, last
                row = rows(j)
                if (binned(best_feature, row) <= best_cut) then
                    work(lc) = row
                    lc = lc + 1
                else
                    work(rc) = row
                    rc = rc + 1
                end if
            end do
            rows(first:last) = work(first:last)
            call build(first, split, depth + 1, out_left)
            if (status /= profit_success) return
            call build(split + 1, last, depth + 1, out_right)
            if (status /= profit_success) return
            self%left(node) = out_left
            self%right(node) = out_right
        end subroutine build
    end subroutine tree_fit

    subroutine tree_leaf_indices(self, features, indices, status)
        class(regression_tree_t), intent(in) :: self
        real(dp), intent(in) :: features(:, :)
        integer, intent(out) :: indices(:)
        integer, intent(out) :: status
        integer :: i, node, f
        logical :: finite_inputs

        indices = 0
        status = profit_invalid_shape
        if (size(features, 1) /= self%n_features) return
        if (size(features, 2) /= size(indices)) return
        status = profit_invalid_input
        if (self%node_count == 0) return
        finite_inputs = .true.
        !$omp parallel do default(shared) private(i) reduction(.and.:finite_inputs)
        do i = 1, size(features, 2)
            if (.not. all(ieee_is_finite(features(:, i)))) finite_inputs = .false.
        end do
        !$omp end parallel do
        if (.not. finite_inputs) return
        !$omp parallel do default(shared) private(i, node, f)
        do i = 1, size(features, 2)
            node = 1
            do
                f = self%feature(node)
                if (f == 0) exit
                if (features(f, i) <= self%threshold(node)) then
                    node = self%left(node)
                else
                    node = self%right(node)
                end if
            end do
            indices(i) = self%leaf(node)
        end do
        !$omp end parallel do
        status = profit_success
    end subroutine tree_leaf_indices

    subroutine tree_predict(self, features, values, status)
        class(regression_tree_t), intent(in) :: self
        real(dp), intent(in) :: features(:, :)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status
        integer, allocatable :: indices(:)
        integer :: i

        values = 0
        status = profit_invalid_shape
        if (size(values) /= size(features, 2)) return
        allocate(indices(size(values)))
        call self%leaf_indices(features, indices, status)
        if (status /= profit_success) return
        !$omp parallel do default(shared) private(i)
        do i = 1, size(values)
            values(i) = self%leaf_value(indices(i))
        end do
        !$omp end parallel do
    end subroutine tree_predict

    subroutine tree_get_leaf_values(self, values, status)
        class(regression_tree_t), intent(in) :: self
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status

        values = 0
        status = profit_invalid_shape
        if (size(values) /= self%leaf_count) return
        status = profit_invalid_input
        if (self%node_count == 0) return
        values = self%leaf_value(:self%leaf_count)
        status = profit_success
    end subroutine tree_get_leaf_values

    subroutine tree_set_leaf_values(self, values, status)
        class(regression_tree_t), intent(inout) :: self
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: status

        status = profit_invalid_shape
        if (size(values) /= self%leaf_count) return
        status = profit_invalid_input
        if (self%node_count == 0) return
        if (.not. all(ieee_is_finite(values))) return
        self%leaf_value(:self%leaf_count) = values
        status = profit_success
    end subroutine tree_set_leaf_values

    integer function tree_n_leaves(self) result(n)
        class(regression_tree_t), intent(in) :: self
        n = self%leaf_count
    end function tree_n_leaves

    subroutine tree_save(self, path, status)
        class(regression_tree_t), intent(in) :: self
        character(*), intent(in) :: path
        integer, intent(out) :: status
        integer :: unit, ios, close_ios, n

        status = profit_invalid_input
        n = self%node_count
        if (n == 0) return
        status = profit_io_error
        open(newunit=unit, file=path, status="replace", form="unformatted", &
            access="stream", iostat=ios)
        if (ios /= 0) return
        write(unit, iostat=ios) "PROFIT_TREE_1", self%n_features, n, self%leaf_count
        if (ios == 0) write(unit, iostat=ios) self%feature(:n), self%left(:n), &
            self%right(:n), self%leaf(:n), self%threshold(:n), &
            self%leaf_value(:self%leaf_count)
        close(unit, iostat=close_ios)
        if (ios /= 0 .or. close_ios /= 0) return
        status = profit_success
    end subroutine tree_save

    subroutine tree_load(self, path, status)
        class(regression_tree_t), intent(out) :: self
        character(*), intent(in) :: path
        integer, intent(out) :: status
        character(13) :: magic
        integer :: unit, ios, p, n, nl, i
        integer, allocatable :: seen(:), parents(:)

        status = profit_io_error
        open(newunit=unit, file=path, status="old", form="unformatted", &
            access="stream", iostat=ios)
        if (ios /= 0) return
        read(unit, iostat=ios) magic, p, n, nl
        if (ios /= 0) then
            close(unit)
            return
        end if
        if (magic /= "PROFIT_TREE_1" .or. p < 1 .or. p > 100000 .or. &
            n < 1 .or. n > 2047 .or. nl < 1 .or. nl > n) then
            close(unit)
            return
        end if
        allocate(self%feature(n), self%left(n), self%right(n), self%leaf(n), &
            self%threshold(n), self%leaf_value(nl), seen(nl), parents(n))
        read(unit, iostat=ios) self%feature, self%left, self%right, self%leaf, &
            self%threshold, self%leaf_value
        close(unit)
        if (ios /= 0) return
        if (.not. all(ieee_is_finite(self%threshold))) return
        if (.not. all(ieee_is_finite(self%leaf_value))) return
        if (any(self%feature < 0) .or. any(self%feature > p)) return
        seen = 0
        parents = 0
        do i = 1, n
            if (self%feature(i) == 0) then
                if (self%left(i) /= 0 .or. self%right(i) /= 0) return
                if (self%leaf(i) < 1 .or. self%leaf(i) > nl) return
                seen(self%leaf(i)) = seen(self%leaf(i)) + 1
            else
                if (self%leaf(i) /= 0) return
                if (self%left(i) <= i) return
                if (self%left(i) > n) return
                if (self%right(i) <= i) return
                if (self%right(i) > n) return
                parents(self%left(i)) = parents(self%left(i)) + 1
                parents(self%right(i)) = parents(self%right(i)) + 1
            end if
        end do
        if (any(seen /= 1)) return
        if (parents(1) /= 0) return
        if (any(parents(2:) /= 1)) return
        self%n_features = p
        self%node_count = n
        self%leaf_count = nl
        status = profit_success
    end subroutine tree_load

    logical function valid_options(options) result(ok)
        type(tree_options_t), intent(in) :: options

        ok = .false.
        if (options%max_depth < 0 .or. options%max_depth > 10) return
        if (options%min_leaf < 1) return
        if (options%max_bins < 2 .or. options%max_bins > 256) return
        if (.not. ieee_is_finite(options%leaf_l2)) return
        if (options%leaf_l2 < 0) return
        ok = .true.
    end function valid_options

    subroutine merge_sort(values, scratch)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(out) :: scratch(:)
        integer :: width, first, middle, last, i, j, k, n

        n = size(values)
        width = 1
        do while (width < n)
            first = 1
            do while (first <= n)
                middle = min(first + width - 1, n)
                last = min(first + 2*width - 1, n)
                i = first
                j = middle + 1
                do k = first, last
                    if (i > middle) then
                        scratch(k) = values(j)
                        j = j + 1
                    else if (j > last) then
                        scratch(k) = values(i)
                        i = i + 1
                    else if (values(i) <= values(j)) then
                        scratch(k) = values(i)
                        i = i + 1
                    else
                        scratch(k) = values(j)
                        j = j + 1
                    end if
                end do
                first = last + 1
            end do
            values = scratch
            width = 2*width
        end do
    end subroutine merge_sort
end module profit_models_tree
