#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::alloc;
use std::ffi::OsStr;
use std::path::Path;
use thiserror::Error;

pub fn encode(nuc: &[u8]) -> Vec<u64> {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { encode_movemask_avx(nuc) };
        } else if is_x86_feature_detected!("sse2") {
            return unsafe { encode_movemask_sse(nuc) };
        }
    }

    encode_lut(nuc)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn encode_movemask_avx(nuc: &[u8]) -> Vec<u64> {
    let ptr = nuc.as_ptr() as *const __m256i;
    let end_idx = nuc.len() / 32;
    let len = end_idx + if nuc.len() % 32 == 0 { 0 } else { 1 };

    let layout = alloc::Layout::from_size_align_unchecked(len * 8, 32);
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

    if nuc.len() % 32 > 0 {
        *res_ptr.offset(end_idx as isize) = *encode_lut(&nuc[(end_idx * 32)..]).get_unchecked(0);
    }

    Vec::from_raw_parts(res_ptr, len, len)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn encode_movemask_sse(nuc: &[u8]) -> Vec<u64> {
    let ptr = nuc.as_ptr() as *const __m128i;
    let end_idx = nuc.len() / 16;
    let len = nuc.len() / 32 + if nuc.len() % 32 == 0 { 0 } else { 1 };

    let layout = alloc::Layout::from_size_align_unchecked(len * 8, 16);
    let res_ptr = alloc::alloc(layout) as *mut u32;

    for i in 0..end_idx as isize {
        let v = _mm_loadu_si128(ptr.offset(i));

        // shift each group of two bits for each nucleotide to the end of each byte
        let lo = _mm_slli_epi64(v, 6);
        let hi = _mm_slli_epi64(v, 5);

        // interleave bytes then extract the bit at the end of each byte
        let a = _mm_unpackhi_epi8(lo, hi);
        let b = _mm_unpacklo_epi8(lo, hi);
        let a = _mm_movemask_epi8(a);
        let b = _mm_movemask_epi8(b);

        *res_ptr.offset(i) = ((a << 16) | b) as u32;
    }

    if nuc.len() % 16 > 0 {
        *res_ptr.offset(end_idx as isize) =
            *encode_lut(&nuc[(end_idx * 16)..]).get_unchecked(0) as u32;
    }

    Vec::from_raw_parts(res_ptr as *mut u64, len, len)
}

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

fn encode_lut(nuc: &[u8]) -> Vec<u64> {
    let mut res = vec![0u64; (nuc.len() / 32) + if nuc.len() % 32 == 0 { 0 } else { 1 }];

    for i in 0..nuc.len() {
        let offset = i / 32;
        let shift = (i % 32) << 1;

        unsafe {
            *res.get_unchecked_mut(offset) = *res.get_unchecked(offset)
                | ((*BYTE_LUT.get_unchecked(*nuc.get_unchecked(i) as usize) as u64) << shift);
        }
    }

    res
}

pub fn decode(bits: &[u64], len: usize) -> Vec<u8> {
    if len > (bits.len() * 32) {
        panic!(
            "The length {} is greater than the number of nucleotides!",
            len
        );
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { decode_shuffle_avx(bits, len) };
        } else if is_x86_feature_detected!("sse4.1") {
            return unsafe { decode_shuffle_sse(bits, len) };
        }
    }

    decode_lut(bits, len)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn decode_shuffle_avx(bits: &[u64], len: usize) -> Vec<u8> {
    let layout = alloc::Layout::from_size_align_unchecked(bits.len() * 32, 32);
    let ptr = alloc::alloc(layout) as *mut __m256i;

    let shuffle_mask = _mm256_set_epi32(
        0x07070707, 0x06060606, 0x05050505, 0x04040404, 0x03030303, 0x02020202, 0x01010101,
        0x00000000,
    );
    let lo_mask = _mm256_set1_epi16(0b0000110000000011);
    let lut_i32 =
        (b'A' as i32) | ((b'C' as i32) << 8) | ((b'T' as i32) << 16) | ((b'G' as i32) << 24);
    let lut = _mm256_set_epi32(
        b'G' as i32,
        b'T' as i32,
        b'C' as i32,
        lut_i32,
        b'G' as i32,
        b'T' as i32,
        b'C' as i32,
        lut_i32,
    );

    for i in 0..bits.len() {
        let curr = *bits.get_unchecked(i) as i64;
        let v = _mm256_set1_epi64x(curr);

        // duplicate each byte four times
        let v1 = _mm256_shuffle_epi8(v, shuffle_mask);

        // separately right shift each 16-bit chunk by 0 or 4 bits
        let v2 = _mm256_srli_epi16(v1, 4);

        // merge together shifted chunks
        let v = _mm256_blend_epi16(v1, v2, 0b10101010i32);

        // only keep two bits in each byte
        // either 0b0011 or 0b1100
        let v = _mm256_and_si256(v, lo_mask);

        // use lookup table to convert nucleotide bits to bytes
        let v = _mm256_shuffle_epi8(lut, v);
        _mm256_store_si256(ptr.offset(i as isize), v);
    }

    Vec::from_raw_parts(ptr as *mut u8, len, bits.len() * 32)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn decode_shuffle_sse(bits: &[u64], len: usize) -> Vec<u8> {
    let layout = alloc::Layout::from_size_align_unchecked(bits.len() * 32, 16);
    let ptr = alloc::alloc(layout) as *mut __m128i;

    let bits_ptr = bits.as_ptr() as *const i32;

    let shuffle_mask = _mm_set_epi32(0x03030303, 0x02020202, 0x01010101, 0x00000000);
    let lo_mask = _mm_set1_epi16(0b0000110000000011);
    let lut_i32 =
        (b'A' as i32) | ((b'C' as i32) << 8) | ((b'T' as i32) << 16) | ((b'G' as i32) << 24);
    let lut = _mm_set_epi32(b'G' as i32, b'T' as i32, b'C' as i32, lut_i32);

    for i in 0..(bits.len() * 2) as isize {
        let curr = *bits_ptr.offset(i);
        let v = _mm_set1_epi32(curr);

        // duplicate each byte four times
        let v1 = _mm_shuffle_epi8(v, shuffle_mask);

        // separately right shift each 16-bit chunk by 0 or 4 bits
        let v2 = _mm_srli_epi16(v1, 4);

        // merge together shifted chunks
        let v = _mm_blend_epi16(v1, v2, 0b10101010i32);

        // only keep two bits in each byte
        // either 0b0011 or 0b1100
        let v = _mm_and_si128(v, lo_mask);

        // use lookup table to convert nucleotide bits to bytes
        let v = _mm_shuffle_epi8(lut, v);
        _mm_store_si128(ptr.offset(i), v);
    }

    Vec::from_raw_parts(ptr as *mut u8, len, bits.len() * 32)
}

static BITS_LUT: [u8; 4] = {
    let mut lut = [0u8; 4];
    lut[0b00] = b'A';
    lut[0b10] = b'T';
    lut[0b01] = b'C';
    lut[0b11] = b'G';
    lut
};

fn decode_lut(bits: &[u64], len: usize) -> Vec<u8> {
    let layout = unsafe { alloc::Layout::from_size_align_unchecked(len, 1) };
    let res_ptr = unsafe { alloc::alloc(layout) };

    for i in 0..len {
        let offset = i >> 5;
        let shift = (i & 31) << 1;
        let curr = unsafe { *bits.get_unchecked(offset) };

        unsafe {
            *res_ptr.offset(i as isize) =
                *BITS_LUT.get_unchecked(((curr >> shift) & 0b11) as usize);
        }
    }

    unsafe { Vec::from_raw_parts(res_ptr, len, len) }
}

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum CliError {
    /// Indicates that a string cannot be parsed into a [`CompressionFormat`](#compressionformat).
    #[error("{0} is not a valid output format")]
    InvalidCompression(String),
}

pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

impl CompressionExt for niffler::compression::Format {
    /// Attempts to infer the compression type from the file extension. If the extension is not
    /// known, then Uncompressed is returned.
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("zst")) => Self::Zstd,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

pub fn parse_compression_format(s: &str) -> Result<niffler::compression::Format, CliError> {
    match s {
        "b" | "B" => Ok(niffler::Format::Bzip),
        "g" | "G" => Ok(niffler::Format::Gzip),
        "l" | "L" => Ok(niffler::Format::Lzma),
        "u" | "U" => Ok(niffler::Format::No),
        "z" | "Z" => Ok(niffler::Format::Zstd),
        _ => Err(CliError::InvalidCompression(s.to_string())),
    }
}

/// A utility function to validate compression level is in allowed range
#[allow(clippy::redundant_clone)]
pub fn parse_level(s: &str) -> Result<niffler::Level, String> {
    let lvl = match s.parse::<u8>() {
        Ok(1) => niffler::Level::One,
        Ok(2) => niffler::Level::Two,
        Ok(3) => niffler::Level::Three,
        Ok(4) => niffler::Level::Four,
        Ok(5) => niffler::Level::Five,
        Ok(6) => niffler::Level::Six,
        Ok(7) => niffler::Level::Seven,
        Ok(8) => niffler::Level::Eight,
        Ok(9) => niffler::Level::Nine,
        Ok(10) => niffler::Level::Ten,
        Ok(11) => niffler::Level::Eleven,
        Ok(12) => niffler::Level::Twelve,
        Ok(13) => niffler::Level::Thirteen,
        Ok(14) => niffler::Level::Fourteen,
        Ok(15) => niffler::Level::Fifteen,
        Ok(16) => niffler::Level::Sixteen,
        Ok(17) => niffler::Level::Seventeen,
        Ok(18) => niffler::Level::Eighteen,
        Ok(19) => niffler::Level::Nineteen,
        Ok(20) => niffler::Level::Twenty,
        Ok(21) => niffler::Level::TwentyOne,
        _ => return Err(format!("Compression level {} not in the range 1-21", s)),
    };
    Ok(lvl)
}
