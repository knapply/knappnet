//! The BLAKE2b hash function (Aumasson, Neves, Wilcox-O'Hearn, and
//! Winnerlein 2013, `aumasson2013blake2`), as specified in RFC 7693
//! (Saarinen and Aumasson 2015, `saarinen2015blake2`), unkeyed, with any
//! digest length from 1 to 64 bytes.
//!
//! The digest length is part of the parameter block, so a 16-byte digest
//! is not a prefix of the 64-byte one. This is the function behind
//! Python's `hashlib.blake2b(data, digest_size = k)`, which NetworkX's
//! Weisfeiler-Lehman graph hashes use; [`hex`] returns the same lowercase
//! hexadecimal string as its `hexdigest()`.

const IV: [u64; 8] = [
    0x6a09_e667_f3bc_c908,
    0xbb67_ae85_84ca_a73b,
    0x3c6e_f372_fe94_f82b,
    0xa54f_f53a_5f1d_36f1,
    0x510e_527f_ade6_82d1,
    0x9b05_688c_2b3e_6c1f,
    0x1f83_d9ab_fb41_bd6b,
    0x5be0_cd19_137e_2179,
];

const SIGMA: [[usize; 16]; 12] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
    [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
    [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
    [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
    [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
    [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
    [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
    [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
];

#[inline(always)]
fn mix(v: &mut [u64; 16], a: usize, b: usize, c: usize, d: usize, x: u64, y: u64) {
    v[a] = v[a].wrapping_add(v[b]).wrapping_add(x);
    v[d] = (v[d] ^ v[a]).rotate_right(32);
    v[c] = v[c].wrapping_add(v[d]);
    v[b] = (v[b] ^ v[c]).rotate_right(24);
    v[a] = v[a].wrapping_add(v[b]).wrapping_add(y);
    v[d] = (v[d] ^ v[a]).rotate_right(16);
    v[c] = v[c].wrapping_add(v[d]);
    v[b] = (v[b] ^ v[c]).rotate_right(63);
}

fn compress(h: &mut [u64; 8], block: &[u8; 128], t: u128, last: bool) {
    let mut m = [0u64; 16];
    for (i, w) in m.iter_mut().enumerate() {
        let mut b = [0u8; 8];
        b.copy_from_slice(&block[8 * i..8 * i + 8]);
        *w = u64::from_le_bytes(b);
    }
    let mut v = [0u64; 16];
    v[..8].copy_from_slice(h);
    v[8..].copy_from_slice(&IV);
    v[12] ^= t as u64;
    v[13] ^= (t >> 64) as u64;
    if last {
        v[14] = !v[14];
    }
    for s in &SIGMA {
        mix(&mut v, 0, 4, 8, 12, m[s[0]], m[s[1]]);
        mix(&mut v, 1, 5, 9, 13, m[s[2]], m[s[3]]);
        mix(&mut v, 2, 6, 10, 14, m[s[4]], m[s[5]]);
        mix(&mut v, 3, 7, 11, 15, m[s[6]], m[s[7]]);
        mix(&mut v, 0, 5, 10, 15, m[s[8]], m[s[9]]);
        mix(&mut v, 1, 6, 11, 12, m[s[10]], m[s[11]]);
        mix(&mut v, 2, 7, 8, 13, m[s[12]], m[s[13]]);
        mix(&mut v, 3, 4, 9, 14, m[s[14]], m[s[15]]);
    }
    for i in 0..8 {
        h[i] ^= v[i] ^ v[i + 8];
    }
}

/// The unkeyed BLAKE2b digest of `data`, `out.len()` bytes long (1 to 64).
///
/// # Panics
/// If `out` is empty or longer than 64 bytes.
pub fn digest_into(data: &[u8], out: &mut [u8]) {
    let nn = out.len();
    assert!((1..=64).contains(&nn), "BLAKE2b digests are 1 to 64 bytes");
    let mut h = IV;
    // Parameter block: digest length, key length 0, fanout 1, depth 1.
    h[0] ^= 0x0101_0000 ^ nn as u64;
    let mut block = [0u8; 128];
    let mut t: u128 = 0;
    let mut rest = data;
    // Every block but the last is compressed as it is seen; the last
    // (possibly partial, possibly empty) block carries the final flag.
    while rest.len() > 128 {
        block.copy_from_slice(&rest[..128]);
        t += 128;
        compress(&mut h, &block, t, false);
        rest = &rest[128..];
    }
    block = [0u8; 128];
    block[..rest.len()].copy_from_slice(rest);
    t += rest.len() as u128;
    compress(&mut h, &block, t, true);
    for (i, b) in out.iter_mut().enumerate() {
        *b = (h[i / 8] >> (8 * (i % 8))) as u8;
    }
}

/// The digest of `data` with `size` bytes as a lowercase hexadecimal
/// string of 2 `size` characters (Python's `hexdigest()`).
pub fn hex(data: &[u8], size: usize) -> String {
    let mut out = vec![0u8; size];
    digest_into(data, &mut out);
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut s = String::with_capacity(2 * size);
    for b in out {
        s.push(DIGITS[(b >> 4) as usize] as char);
        s.push(DIGITS[(b & 15) as usize] as char);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rfc_7693_test_vector() {
        // Appendix A of RFC 7693: BLAKE2b-512("abc").
        assert_eq!(
            hex(b"abc", 64),
            "ba80a53f981c4d0d6a2797b69f12f6e94c212f14685ac4b74b12bb6fdbffa2d1\
             7d87c5392aab792dc252d5de4533cc9518d38aa8dbf1925ab92386edd4009923"
        );
    }

    #[test]
    fn known_digests_of_python_hashlib() {
        // hashlib.blake2b(b"", digest_size = 64), a 16-byte digest, and
        // 200 and 128 bytes "a" (two blocks and exactly one).
        assert_eq!(
            hex(b"", 64),
            "786a02f742015903c6c6fd852552d272912f4740e15847618a86e217f71f5419\
             d25e1031afee585313896444934eb04b903a685b1448b755d56f701afe9be2ce"
        );
        assert_eq!(hex(b"abc", 16), "cf4ab791c62b8d2b2109c90275287816");
        let a = [b'a'; 200];
        assert_eq!(hex(&a, 16), "d63a46941ee0449fd1a7292005e08096");
        // Exactly one block: the final flag goes on the first block.
        assert_eq!(hex(&a[..128], 16), "1271d28e731afb57bdfba6bdaf3501be");
    }
}
