// const modulus = 2**256 - 2**32 * 351 + 1;
// const non_residue = 7;
// const spot_check_security_factor = 80;
// const extension_factor = 8;

pub fn mimc(inp: i32, steps: i32, round_constants: &[i32]) -> i32 {
  let modulus = 2**256 - 2**32 * 351 + 1;
  // start_time = time.time();
  for i in 0..(steps - 1) {
    inp = (inp.pow(3) + round_constants[i % round_constants.len()]) % modulus;
  }
  // print("MIMC computed in %.4f sec" % (time.time() - start_time))
  inp
}
