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
