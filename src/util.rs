/// Check if a number is a power of two.
///
/// This is useful in several places, since we require that various objects
/// have a nice size, in order to make operations easier.
pub fn is_power_of_two(x: usize) -> bool {
    // A power of two will have exactly one set bit.
    // This will actually generate very good assembly.
    // https://rust.godbolt.org/z/WaE9rvbrG
    x.count_ones() == 1
}
