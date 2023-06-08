#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;


use std::alloc;

static BYTE_LUT: [u8; 128] = {
    let mut lut = [0u8; 128];
    lut[b'a' as usize] = 0b00;
    lut[b't' as usize] = 0b10;
    lut[b'u' as usize] = 0b10;
    lut[b'c' as usize] = 0b01;
    lut[b'g' as usize] = 0b11;
    lut[b'A' as usize] = 0b00;
    lut[b'T' as usize] = 0b10;
    lut[b'U' as usize] = 0b10;
    lut[b'C' as usize] = 0b01;
    lut[b'G' as usize] = 0b11;
    lut
};

static BITS_LUT: [u8; 4] = {
    let mut lut = [0u8; 4];
    lut[0b00] = b'A';
    lut[0b10] = b'T';
    lut[0b01] = b'C';
    lut[0b11] = b'G';
    lut
};

/// Encode `{A, T/U, C, G}` from the byte string into pairs of bits (`{00, 10, 01, 11}`) packed into 64-bit integers,
/// by using a naive scalar method.
pub fn n_to_bits_lut(n: &[u8]) -> Vec<u64> {
    let mut res = vec![0u64; (n.len() >> 5) + if n.len() & 31 == 0 { 0 } else { 1 }];

    unsafe {
        for i in 0..n.len() {
            let offset = i >> 5;
            let shift = (i & 31) << 1;
            *res.get_unchecked_mut(offset) = *res.get_unchecked(offset)
                | ((*BYTE_LUT.get_unchecked(*n.get_unchecked(i) as usize) as u64) << shift);
        }
    }

    res
}

/// Decode pairs of bits from packed 64-bit integers to get a byte string of `{A, T/U, C, G}`, by using a naive scalar
/// method.
pub fn bits_to_n_lut(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() << 5) {
        panic!("The length is greater than the number of nucleotides!");
    }

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(len, 1);
        let res_ptr = alloc::alloc(layout);

        for i in 0..len {
            let offset = i >> 5;
            let shift = (i & 31) << 1;
            let curr = *bits.get_unchecked(offset);
            *res_ptr.offset(i as isize) = *BITS_LUT.get_unchecked(((curr >> shift) & 0b11) as usize);
        }

        Vec::from_raw_parts(res_ptr, len, len)
    }
}

/// Encode `{A, T/U, C, G}` from the byte string into pairs of bits (`{00, 10, 01, 11}`) packed into 64-bit integers,
/// by using a vectorized method with the `permute4x64`, `unpack`, and `movemask` instructions.
///
/// Requires AVX2 support.
pub fn n_to_bits_movemask(n: &[u8]) -> Vec<u64> {
    let ptr = n.as_ptr() as *const __m256i;
    let end_idx = n.len() >> 5;
    let len = end_idx + if n.len() & 31 == 0 {0} else {1};

    unsafe {
        let layout = alloc::Layout::from_size_align_unchecked(len << 3, 8);
        let res_ptr = alloc::alloc(layout) as *mut u64;
        for i in 0..end_idx as isize {
            let v = _mm256_loadu_si256(ptr.offset(i));

            // permute because unpacks works on the low/high 64 bits in each lane
            let v = _mm256_permute4x64_epi64(v, 0b11011000);

            // shift each group of two bits for each nucleotide to the end of each byte
            let lo = _mm256_slli_epi64(v, 6);
            let hi = _mm256_slli_epi64(v, 5);

            // interleave bytes then extract the bit at the end of each byte
            let a = _mm256_unpackhi_epi8(lo, hi);
            let b = _mm256_unpacklo_epi8(lo, hi);

            // zero extend after movemask
            let a = (_mm256_movemask_epi8(a) as u32) as u64;
            let b = (_mm256_movemask_epi8(b) as u32) as u64;

            *res_ptr.offset(i) = (a << 32) | b;
        }

        if n.len() & 31 > 0 {
            *res_ptr.offset(end_idx as isize) = *n_to_bits_lut(&n[(end_idx << 5)..]).get_unchecked(0);
        }

        Vec::from_raw_parts(res_ptr, len, len)
    }
}