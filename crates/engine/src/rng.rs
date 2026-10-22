//! Counter-based random numbers.
//!
//! Philox4x32-10 (Salmon et al. 2011, `salmon2011parallel`) turns a key
//! and a counter into random bits with no shared state. Every random
//! stream is named by (seed, stream index), so a permutation test gives
//! identical results whether it runs on one thread or many, and on any
//! platform. The same generator is standard on GPUs, which keeps results
//! reproducible if accelerated backends are added later.

const M0: u32 = 0xD251_1F53;
const M1: u32 = 0xCD9E_8D57;
const W0: u32 = 0x9E37_79B9;
const W1: u32 = 0xBB67_AE85;

#[inline]
fn mulhilo(a: u32, b: u32) -> (u32, u32) {
    let p = a as u64 * b as u64;
    ((p >> 32) as u32, p as u32)
}

/// One Philox4x32-10 block: 128 random bits for a counter and key.
pub fn philox4x32_10(counter: [u32; 4], key: [u32; 2]) -> [u32; 4] {
    let (mut c, mut k) = (counter, key);
    for round in 0..10 {
        if round > 0 {
            k[0] = k[0].wrapping_add(W0);
            k[1] = k[1].wrapping_add(W1);
        }
        let (hi0, lo0) = mulhilo(M0, c[0]);
        let (hi1, lo1) = mulhilo(M1, c[2]);
        c = [hi1 ^ c[1] ^ k[0], lo1, hi0 ^ c[3] ^ k[1], lo0];
    }
    c
}

/// An independent random stream identified by a seed and a stream index.
#[derive(Debug, Clone)]
pub struct Stream {
    key: [u32; 2],
    stream: u64,
    block: u64,
    buffer: [u32; 4],
    used: usize,
}

impl Stream {
    /// Stream `index` of the generator seeded with `seed`.
    pub fn new(seed: u64, index: u64) -> Stream {
        Stream {
            key: [seed as u32, (seed >> 32) as u32],
            stream: index,
            block: 0,
            buffer: [0; 4],
            used: 4,
        }
    }

    /// The next 32 random bits.
    pub fn next_u32(&mut self) -> u32 {
        if self.used == 4 {
            let counter = [
                self.block as u32,
                (self.block >> 32) as u32,
                self.stream as u32,
                (self.stream >> 32) as u32,
            ];
            self.buffer = philox4x32_10(counter, self.key);
            self.block += 1;
            self.used = 0;
        }
        let x = self.buffer[self.used];
        self.used += 1;
        x
    }

    /// The next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        (self.next_u32() as u64) << 32 | self.next_u32() as u64
    }

    /// A uniform integer in `0..n` (n > 0), without modulo bias.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "below(0)");
        // Rejection from the largest multiple of n that fits in 64 bits.
        let zone = u64::MAX - (u64::MAX % n + 1) % n;
        loop {
            let x = self.next_u64();
            if x <= zone {
                return x % n;
            }
        }
    }

    /// A uniform double in [0, 1) with 53 random bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Fills `perm` with a uniformly random permutation of `0..perm.len()`
    /// (Fisher and Yates; Durstenfeld 1964, `durstenfeld1964algorithm`).
    pub fn permutation(&mut self, perm: &mut [usize]) {
        for (i, p) in perm.iter_mut().enumerate() {
            *p = i;
        }
        for i in (1..perm.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            perm.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_answers() {
        // Known-answer vectors published with Random123 (Salmon et al. 2011).
        assert_eq!(
            philox4x32_10([0; 4], [0; 2]),
            [0x6627_e8d5, 0xe169_c58d, 0xbc57_ac4c, 0x9b00_dbd8]
        );
        assert_eq!(
            philox4x32_10([u32::MAX; 4], [u32::MAX; 2]),
            [0x408f_276d, 0x41c8_3b0e, 0xa20b_c7c6, 0x6d54_51fd]
        );
        assert_eq!(
            philox4x32_10(
                [0x243f_6a88, 0x85a3_08d3, 0x1319_8a2e, 0x0370_7344],
                [0xa409_3822, 0x299f_31d0]
            ),
            [0xd16c_fe09, 0x94fd_cceb, 0x5001_e420, 0x2412_6ea1]
        );
    }

    #[test]
    fn streams_are_reproducible_and_distinct() {
        let a: Vec<u64> = (0..5).map(|_| Stream::new(7, 3).next_u64()).collect();
        assert!(a.windows(2).all(|w| w[0] == w[1]));
        assert_ne!(Stream::new(7, 3).next_u64(), Stream::new(7, 4).next_u64());
        assert_ne!(Stream::new(7, 3).next_u64(), Stream::new(8, 3).next_u64());
    }

    #[test]
    fn permutations_are_valid_and_roughly_uniform() {
        let mut counts = [[0u32; 4]; 4];
        let mut perm = [0usize; 4];
        for i in 0..40_000 {
            Stream::new(1, i).permutation(&mut perm);
            let mut seen = [false; 4];
            for (pos, &v) in perm.iter().enumerate() {
                assert!(!seen[v]);
                seen[v] = true;
                counts[pos][v] += 1;
            }
        }
        // Each value lands in each position about 10,000 times.
        for row in counts {
            for c in row {
                assert!((9_400..10_600).contains(&c), "{c}");
            }
        }
    }
}
