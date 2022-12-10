use ark_bls12_381::Bls12_381;
use ark_ec::AffineCurve;
use ark_serialize::CanonicalDeserialize;
use ark_std::{UniformRand, Zero};
use prompt::{puzzle, welcome};

pub mod algorithms;
pub mod data_structures;
use data_structures::*;

pub mod attack;

fn main() {
    welcome();
    puzzle(PUZZLE_DESCRIPTION);
    // Supports committing to vectors of length up to 512.
    let ck = data_structures::CommitmentKey::<Bls12_381>::deserialize_unchecked(SRS).unwrap();
    let attack = attack(&ck, SUPPORTED_DIM);
    attack.assert_attack_works(&ck, SUPPORTED_DIM);
}

pub fn attack<E: ark_ec::PairingEngine>(ck: &CommitmentKey<E>, dim: usize) -> attack::Attack<E> {
    // The first thing we note is the length of the first part of the SRS
    // It should only only have dim+1 elements for the beta^(0..dim) terms
    assert_eq!(dim+2, ck.powers_of_beta_g_first.len());

    // We can check that this last element is indeed g^(beta^(dim+1)) using pairings
    // (We assume that the other values are constructed correctly)
    assert_eq!(
        E::pairing(
            ck.powers_of_beta_g_first[dim],
            ck.powers_of_beta_h[1]
        ),
        E::pairing(
            ck.powers_of_beta_g_first[dim+1],
            ck.powers_of_beta_h[0]
        )
    );

    // Having access to the beta^(dim+1) term means that a prover can fake a witness
    // to any inner product m_fake (for any committed vector a and public vector b), like so:
    // 
    // - Compute g^(AB(beta)). This can be done directly, or by computing the witness 
    //   to the actual inner product and adding the dim+1 term back in. 
    //      - This latter case is interesting since it means that upon receipt of any witness
    //        that verifies, "real" or not, anyone else can also fake a proof for commit(a) and 
    //        b without having to know the original vector that was committed to.
    // - Subtract out (g^(beta^(dim+1)))^m_fake and set the result as your witness W
    // - Return m_fake and W to the verifier

    // If we commit to the zero vector, then A is the zero polynomial and so is AB.
    // In this special case, we can fake a witness without even knowing b beforehand, since AB is constant.
    let a = (0..dim).map(|_| E::Fr::zero()).collect::<Vec<_>>();
    let cm = algorithms::ILV::commit(&ck, &a);
    assert_eq!(cm.0, E::G1Affine::zero()); // A bit conspicuous o_O

    let mut rng = ark_std::test_rng();
    let claimed_inner_product = E::Fr::rand(&mut rng); // Anything we want
    let proof = Proof(
        ck.powers_of_beta_g_first[dim+1].mul(-claimed_inner_product).into()
    );

    attack::Attack {
        a,
        commitment: cm,
        claimed_inner_product,
        proof,
    }
}

const SRS: &'static [u8] = include_bytes!("../ck.srs");
const SUPPORTED_DIM: usize = 512;

const PUZZLE_DESCRIPTION: &str = r"
Bob was catching up on the latest in zkSNARK research, and came across the
Vampire paper [1]. In that paper, he found a reference to an inner-product
commitment scheme [2], which allows committing to a vector and later proving
that its inner-product with another (public) vector is equal to a claimed value.
Bob was intrigued by this scheme, and decided to implement it in Rust.

Bob was delighted with the performance of the resulting implementation, and so
decided to deploy it. The scheme requires a universal Powers-of-Tau-type trusted setup, 
and so Bob generated a SRS using an MPC ceremony.

Things were going smoothly for a while, but then Bob received an anonymous email that 
contained a full break of the scheme! Unfortunately for Bob, the email didn't contain
any details about the break. Can you help Bob figure out the issue, and fix his scheme?

[1]: https://ia.cr/2022/406
[2]: http://www.lsv.fr/Publis/PAPERS/PDF/ILV-imacc11-long.pdf
";
