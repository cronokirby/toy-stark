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

/// Calculate the binary logarithm of some number.
///
/// This is returned as u32, which is convenient to use for different numbers.
pub fn lg(x: usize) -> u32 {
    x.trailing_zeros()
}

/// Reverse the first few bits of a number.
///
/// We take an extra parameter denoting the number of bits we want to reverse.
/// So, for a 4 bit number, we can just touch only the first 4 bits.
pub fn reverse(x: usize, bits: u32) -> usize {
    let shift = usize::BITS
        .checked_sub(bits)
        .expect("bit count is greater than usize");
    x.reverse_bits() >> shift
}
